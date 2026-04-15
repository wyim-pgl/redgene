#!/usr/bin/env python3
"""Step 5 — Targeted insert assembly via soft-clip junction detection,
strand-aware k-mer extension, and Pilon iterative gap filling.

Pipeline:
  Phase 1: Scan host BAM for bidirectional soft-clip clusters → insertion sites
  Phase 2: Extract candidate reads from junction regions + unmapped pairs
  Phase 3: Iterative k-mer extension + minimap2 soft-clip extension + Pilon gap fill (max 15 rounds)
  Phase 4: Annotate via local element_db BLAST + remote NCBI nt BLAST

No external assembler (SPAdes, SSAKE, TASR) is used.
External tools: minimap2, samtools, Pilon, blastn.

Inputs:
  - Host BAM from step 7 (or junctions.tsv from step 6 as fallback)
  - Host reference FASTA
  - Element database FASTA for annotation

Outputs:
  - insert_only.fasta   — assembled insert sequence(s)
  - element_annotation.tsv — BLAST hits along insert
  - border_hits.tsv      — T-DNA border motif locations
  - insert_report.txt    — human-readable linear map
  - s05_stats.txt        — convergence and assembly statistics
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pysam

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

STEP = "s05_insert_assembly"

# ---------------------------------------------------------------------------
# Host-endogenous exclusion thresholds (used in classify_site_tiers)
# ---------------------------------------------------------------------------
# DB-level: BLAST transgene_db entries vs host genome to find host-derived
# elements (e.g., P-Act1-rice in a rice genome). Two tiers handle cultivar
# drift — P-Act1 hits Nipponbare at ~78%/39% due to Xiushui vs Nipponbare.
HOST_ENDO_T1_PIDENT = 90.0   # Tier 1 (clean host): min % identity
HOST_ENDO_T1_QCOVS = 50.0    # Tier 1: min query coverage %
HOST_ENDO_T2_PIDENT = 75.0   # Tier 2 (divergent host): min % identity
HOST_ENDO_T2_QCOVS = 30.0    # Tier 2: min query coverage %

# Per-clip host verification: stricter than DB-level because individual clips
# are short (~20-50 bp), so even moderate-identity host hits are significant.
# DB-level Tier 2 uses 75% because full-length elements (500-3000 bp) may
# diverge across cultivars; per-clip uses 95% because a 30 bp clip hitting
# host at <95% is too noisy to trust as evidence of host origin.
CLIP_HOST_PIDENT = 95.0       # Per-clip host check: min % identity
CLIP_HOST_MIN_LEN = 30        # Per-clip host check: min alignment length

# Post-assembly host-fraction filter: BLAST assembled insert vs host genome.
# A chimeric assembly artifact has high host fraction AND only a tiny non-host
# gap.  True insertions also have high host fraction (junction contigs span
# host→T-DNA→host), but the non-host gap is large (≥1 kb of T-DNA).
# Both conditions must hold to call FALSE_POSITIVE:
#   host_fraction ≥ INSERT_HOST_FRACTION  AND  largest_gap < INSERT_MIN_FOREIGN_GAP
# Example: insertion_22966 = 96% host, 461bp gap → FP.
#          insertion_32461 = 87% host, 2800bp gap → true CANDIDATE.
INSERT_HOST_FRACTION = 0.80     # host coverage threshold
INSERT_MIN_FOREIGN_GAP = 500    # non-host gap must be ≥ this to be real T-DNA
INSERT_HOST_MIN_PIDENT = 90.0   # min identity for host alignment to count

# Construct-flanking filter: if a construct reference entry contains host
# genomic DNA at its ends (common when constructs are cloned with flanking),
# sites at those host coordinates are false detections.
CONSTRUCT_FLANK_PIDENT = 95.0   # min identity for construct→host flanking hit
CONSTRUCT_FLANK_MIN_LEN = 50    # min alignment length
CONSTRUCT_FLANK_SLOP = 500      # bp slop when checking site overlap

# Multi-locus chimeric filter: if an assembled insert's host-aligned portions
# map to ≥2 different chromosomes (besides the site's own), the assembly
# merged reads from unrelated loci sharing element homology.
# Uses strict identity (≥98%) to distinguish actual chimeric DNA pieces from
# low-level element homologies (e.g., 35S promoter paralogs at 80-90%).
CHIMERIC_MIN_PIDENT = 98.0      # strict identity for chimeric detection
CHIMERIC_MIN_OFFTARGET_BP = 150 # min bp on off-target chromosome to count

# Filter D: Construct-host coverage — if the assembled insert is fully
# explained by construct + host sequences (high combined coverage) AND the
# construct portion is large (≥30%), the insert is host genomic DNA with
# construct-element homology, not a real T-DNA insertion.
# Real T-DNA inserts have low construct coverage (~10%) because only border
# regions match; FPs have high construct coverage (~50%) because the insert
# IS the construct fragment.
CONSTRUCT_HOST_MIN_COMBINED = 0.85  # min combined coverage (construct + host)
CONSTRUCT_MIN_FRACTION = 0.25       # min construct coverage to suspect FP
CONSTRUCT_HOST_MIN_PIDENT = 80.0    # identity threshold for construct BLAST

# UNKNOWN → FALSE_POSITIVE auto-reclassification: if insert has no element
# annotations but is mostly host DNA with negligible construct match,
# it is host genomic DNA, not a real transgene insertion.
UNKNOWN_HOST_MIN_FRACTION = 0.85   # min host fraction to classify as host-only
UNKNOWN_MAX_CONSTRUCT_FRAC = 0.05  # max construct fraction for host-only


def log(msg: str) -> None:
    print(f"[{STEP}] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Reverse complement
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class JunctionCluster:
    """A cluster of soft-clipped reads at a genomic position."""
    host_chr: str
    position: int           # median position of cluster
    clip_direction: str     # 'right' (fwd read 3' clip) or 'left' (rev read 5' clip)
    clipped_seqs: list[str] = field(default_factory=list)
    consensus_clip: str = ""
    n_reads: int = 0
    maps_to_host: bool = False
    element_hits: list[str] = field(default_factory=list)


@dataclass
class InsertionSite:
    """An insertion site with paired 5'/3' junction clusters."""
    site_id: str
    host_chr: str
    junction_5p: JunctionCluster | None = None  # forward reads clipped on right
    junction_3p: JunctionCluster | None = None  # reverse reads clipped on left
    confidence: str = "low"
    seed_5p: str = ""       # consensus clip from 5' junction (insert start)
    seed_3p: str = ""       # consensus clip from 3' junction (insert end)
    is_validated: bool = False

    # Validation details
    clips_are_different: bool = False
    clips_not_in_host: bool = False
    has_element_hits: bool = False

    # Positions for read extraction
    pos_5p: int = 0
    pos_3p: int = 0


@dataclass
class LegacyJunction:
    """Junction from step 6 junctions.tsv (fallback mode)."""
    contig_name: str
    contig_len: int
    host_chr: str
    host_start: int
    host_end: int
    host_strand: str
    construct_element: str
    construct_start: int
    construct_end: int
    junction_pos_host: int
    junction_type: str
    confidence: str
    host_mapq: int


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = ""
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name:
        seqs[name] = "".join(parts)
    return seqs


def write_fasta(path: Path, name: str, seq: str, wrap: int = 80) -> None:
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), wrap):
            fh.write(seq[i:i + wrap] + "\n")


def _read_fq_seqs(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs: list[str] = []
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip().upper())
    return seqs


# ---------------------------------------------------------------------------
# Phase 1: Soft-clip Junction Detection
# ---------------------------------------------------------------------------

def _build_consensus(seqs: list[str], direction: str) -> str:
    """Build majority-vote consensus from a list of clipped sequences.

    For 'right' clips: sequences are aligned from their left (start of insert).
    For 'left' clips: sequences are aligned from their right (end of insert).
    """
    if not seqs:
        return ""

    if direction == "right":
        # Align from left
        votes: dict[int, Counter] = defaultdict(Counter)
        for seq in seqs:
            for i, base in enumerate(seq):
                votes[i][base.upper()] += 1
        consensus = []
        for p in range(max(votes.keys()) + 1):
            if p not in votes:
                break
            v = votes[p]
            total = sum(v.values())
            best, cnt = v.most_common(1)[0]
            if total >= 2 and cnt / total >= 0.51:
                consensus.append(best)
            else:
                break
        return "".join(consensus)
    else:
        # Align from right
        votes = defaultdict(Counter)
        max_len = max(len(s) for s in seqs)
        for seq in seqs:
            offset = max_len - len(seq)
            for i, base in enumerate(seq):
                votes[offset + i][base.upper()] += 1
        consensus = []
        for p in range(max_len):
            if p not in votes:
                consensus.append("N")
                continue
            v = votes[p]
            total = sum(v.values())
            best, cnt = v.most_common(1)[0]
            if total >= 2 and cnt / total >= 0.51:
                consensus.append(best)
            else:
                consensus.append("N")
        # Trim leading/trailing N
        result = "".join(consensus).strip("N")
        return result


def _batch_check_maps_to_host(seqs: dict[str, str], host_ref: Path, workdir: Path,
                              min_identity: float = 0.90,
                              min_coverage: float = 0.80) -> set[str]:
    """Batch check which sequences map to host genome using one minimap2 call.

    Args:
        seqs: dict of {name: sequence} to check
        Returns: set of names that map to host
    """
    if not seqs:
        return set()

    query_fa = workdir / "_clip_check_batch.fa"
    with open(query_fa, "w") as fh:
        for name, seq in seqs.items():
            if len(seq) >= 20:
                fh.write(f">{name}\n{seq}\n")

    result = subprocess.run(
        ["minimap2", "-c", "--secondary=no", "-t", "4",
         str(host_ref), str(query_fa)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    query_fa.unlink(missing_ok=True)

    maps_to_host: set[str] = set()
    for line in result.stdout.splitlines():
        cols = line.strip().split("\t")
        if len(cols) < 12:
            continue
        qname = cols[0]
        q_len = int(cols[1])
        match_bp = int(cols[9])
        block_len = int(cols[10])
        if block_len == 0:
            continue
        identity = match_bp / block_len
        coverage = (int(cols[3]) - int(cols[2])) / q_len if q_len > 0 else 0
        if identity >= min_identity and coverage >= min_coverage:
            maps_to_host.add(qname)
    return maps_to_host


def _batch_check_element_hits(
    seqs: dict[str, str],
    element_db: Path,
    workdir: Path,
    extra_dbs: list[Path] | None = None,
) -> dict[str, list[str]]:
    """Batch check which sequences hit element DB(s) using blastn.

    When ``extra_dbs`` is provided (e.g., the always-on ``common_payload.fa``
    and/or per-sample SPAdes contigs from s04b), each query is BLASTed
    against every DB in turn and hits are merged. This lets shared
    transgene payloads and sample-specific assemblies contribute to
    element annotation without rebuilding the shared DB.
    """
    if not seqs:
        return {}

    query_fa = workdir / "_clip_element_batch.fa"
    with open(query_fa, "w") as fh:
        for name, seq in seqs.items():
            if len(seq) >= 20:
                fh.write(f">{name}\n{seq}\n")

    hits: dict[str, list[str]] = defaultdict(list)
    dbs = [element_db]
    for edb in (extra_dbs or []):
        if edb is not None and edb.exists() and edb.stat().st_size > 0:
            dbs.append(edb)

    for i, db in enumerate(dbs):
        blast_out = workdir / f"_clip_blast_{i}_{db.stem}.tsv"
        subprocess.run(
            ["blastn", "-query", str(query_fa), "-subject", str(db),
             "-outfmt", "6 qseqid sseqid pident length",
             "-evalue", "1e-3", "-max_target_seqs", "5",
             "-out", str(blast_out)],
            stderr=subprocess.DEVNULL,
        )
        if blast_out.exists():
            with open(blast_out) as fh:
                for line in fh:
                    cols = line.strip().split("\t")
                    if len(cols) >= 4 and int(cols[3]) >= 20:
                        hits[cols[0]].append(cols[1])
            blast_out.unlink(missing_ok=True)

    query_fa.unlink(missing_ok=True)
    return dict(hits)


def find_softclip_junctions(
    host_bam: Path,
    host_ref: Path,
    element_db: Path | None,
    workdir: Path,
    min_clip: int = 20,
    cluster_window: int = 50,
    extra_dbs: list[Path] | None = None,
) -> list[InsertionSite]:
    """Scan host BAM for insertion sites using bidirectional soft-clip analysis.

    Three conditions for a validated insertion site:
    1. Bidirectional soft-clips at the same genomic position
    2. The two clip sequences are DIFFERENT (insert, not SV)
    3. Neither clip maps to host genome (foreign sequence)
    """
    workdir.mkdir(parents=True, exist_ok=True)
    log("Phase 1: Scanning host BAM for soft-clip junctions...")

    # Step 1: Collect all soft-clipped reads
    right_clips: dict[str, list[tuple[int, str]]] = defaultdict(list)  # chr -> [(pos, seq)]
    left_clips: dict[str, list[tuple[int, str]]] = defaultdict(list)

    bam = pysam.AlignmentFile(str(host_bam), "rb")
    n_right = 0
    n_left = 0

    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        cigar = read.cigartuples
        if cigar is None:
            continue
        seq = read.query_sequence
        if seq is None:
            continue

        chrom = read.reference_name

        # Right clip (3' end clipped — forward read approaching junction)
        if cigar[-1][0] == 4 and cigar[-1][1] >= min_clip:
            clip_len = cigar[-1][1]
            clip_pos = read.reference_end  # position where clip starts
            clip_seq = seq[-clip_len:]
            right_clips[chrom].append((clip_pos, clip_seq))
            n_right += 1

        # Left clip (5' end clipped — reverse read approaching junction)
        if cigar[0][0] == 4 and cigar[0][1] >= min_clip:
            clip_len = cigar[0][1]
            clip_pos = read.reference_start  # position where clip ends
            clip_seq = seq[:clip_len]
            left_clips[chrom].append((clip_pos, clip_seq))
            n_left += 1

    bam.close()
    log(f"  Collected {n_right:,} right-clips, {n_left:,} left-clips")

    # Step 2: Cluster by position (min_depth=3 to filter noise)
    MIN_CLUSTER_DEPTH = 3

    def _cluster(clips: list[tuple[int, str]], window: int) -> list[tuple[int, list[str]]]:
        if not clips:
            return []
        clips.sort(key=lambda x: x[0])
        clusters = []
        current_seqs = [clips[0][1]]
        current_positions = [clips[0][0]]

        for pos, seq in clips[1:]:
            if pos - current_positions[0] <= window:
                current_seqs.append(seq)
                current_positions.append(pos)
            else:
                if len(current_seqs) >= MIN_CLUSTER_DEPTH:
                    median_pos = sorted(current_positions)[len(current_positions) // 2]
                    clusters.append((median_pos, current_seqs))
                current_seqs = [seq]
                current_positions = [pos]
        if len(current_seqs) >= MIN_CLUSTER_DEPTH:
            median_pos = sorted(current_positions)[len(current_positions) // 2]
            clusters.append((median_pos, current_seqs))

        return clusters

    # Step 3: Pair forward/reverse clusters at same position
    sites: list[InsertionSite] = []
    used_ids: set[str] = set()

    for chrom in set(list(right_clips.keys()) + list(left_clips.keys())):
        r_clusters = _cluster(right_clips.get(chrom, []), cluster_window)
        l_clusters = _cluster(left_clips.get(chrom, []), cluster_window)

        # Try to pair right and left clusters within window
        used_l = set()
        for r_pos, r_seqs in r_clusters:
            best_l = None
            best_dist = cluster_window + 1
            for li, (l_pos, l_seqs) in enumerate(l_clusters):
                if li in used_l:
                    continue
                dist = abs(r_pos - l_pos)
                if dist <= cluster_window and dist < best_dist:
                    best_l = li
                    best_dist = dist

            if best_l is not None:
                l_pos, l_seqs = l_clusters[best_l]
                used_l.add(best_l)

                r_consensus = _build_consensus(r_seqs, "right")
                l_consensus = _build_consensus(l_seqs, "left")

                if len(r_consensus) < min_clip or len(l_consensus) < min_clip:
                    continue

                anchor_pos = min(r_pos, l_pos)
                sid = f"insertion_{chrom}_{anchor_pos}"
                if sid in used_ids:
                    suffix = 2
                    while f"{sid}_{suffix}" in used_ids:
                        suffix += 1
                    sid = f"{sid}_{suffix}"
                used_ids.add(sid)
                jc_5p = JunctionCluster(
                    host_chr=chrom, position=r_pos,
                    clip_direction="right",
                    clipped_seqs=r_seqs, consensus_clip=r_consensus,
                    n_reads=len(r_seqs),
                )
                jc_3p = JunctionCluster(
                    host_chr=chrom, position=l_pos,
                    clip_direction="left",
                    clipped_seqs=l_seqs, consensus_clip=l_consensus,
                    n_reads=len(l_seqs),
                )

                site = InsertionSite(
                    site_id=sid,
                    host_chr=chrom,
                    junction_5p=jc_5p,
                    junction_3p=jc_3p,
                    seed_5p=r_consensus,
                    seed_3p=l_consensus,
                    pos_5p=min(r_pos, l_pos),
                    pos_3p=max(r_pos, l_pos),
                )
                sites.append(site)
            else:
                # Single-direction cluster: only keep high-depth (≥5 reads)
                if len(r_seqs) >= 5:
                    r_consensus = _build_consensus(r_seqs, "right")
                    if len(r_consensus) >= min_clip:
                        sid = f"insertion_{chrom}_{r_pos}"
                        if sid in used_ids:
                            suffix = 2
                            while f"{sid}_{suffix}" in used_ids:
                                suffix += 1
                            sid = f"{sid}_{suffix}"
                        used_ids.add(sid)
                        jc_5p = JunctionCluster(
                            host_chr=chrom, position=r_pos,
                            clip_direction="right",
                            clipped_seqs=r_seqs, consensus_clip=r_consensus,
                            n_reads=len(r_seqs),
                        )
                        site = InsertionSite(
                            site_id=sid,
                            host_chr=chrom,
                            junction_5p=jc_5p,
                            seed_5p=r_consensus,
                            pos_5p=r_pos,
                            confidence="low",
                        )
                        sites.append(site)

        # Unpaired left clusters (only high-depth)
        for li, (l_pos, l_seqs) in enumerate(l_clusters):
            if li in used_l:
                continue
            if len(l_seqs) < 5:
                continue
            l_consensus = _build_consensus(l_seqs, "left")
            if len(l_consensus) >= min_clip:
                sid = f"insertion_{chrom}_{l_pos}"
                if sid in used_ids:
                    suffix = 2
                    while f"{sid}_{suffix}" in used_ids:
                        suffix += 1
                    sid = f"{sid}_{suffix}"
                used_ids.add(sid)
                jc_3p = JunctionCluster(
                    host_chr=chrom, position=l_pos,
                    clip_direction="left",
                    clipped_seqs=l_seqs, consensus_clip=l_consensus,
                    n_reads=len(l_seqs),
                )
                site = InsertionSite(
                    site_id=sid,
                    host_chr=chrom,
                    junction_3p=jc_3p,
                    seed_3p=l_consensus,
                    pos_5p=l_pos,
                    confidence="low",
                )
                sites.append(site)

    log(f"  Found {len(sites)} candidate insertion sites")

    # Step 4: Validate — batched external tool calls for speed
    # Separate paired (bidirectional) vs single-direction sites
    paired_sites = [s for s in sites
                    if s.junction_5p is not None and s.junction_3p is not None]
    single_sites = [s for s in sites
                    if s.junction_5p is None or s.junction_3p is None]

    log(f"  Paired (bidirectional): {len(paired_sites)}, "
        f"Single-direction: {len(single_sites)}")

    # Condition 2 filter on paired sites: clips must be different
    cond2_passed: list[InsertionSite] = []
    for site in paired_sites:
        cmp_len = min(20, len(site.seed_5p), len(site.seed_3p))
        site.clips_are_different = site.seed_5p[:cmp_len] != site.seed_3p[:cmp_len]
        if not site.clips_are_different:
            continue
        cond2_passed.append(site)

    log(f"  After clip-difference filter: {len(cond2_passed)} paired sites")

    # Condition 3: batch host-mapping check (one minimap2 call)
    host_check_seqs: dict[str, str] = {}
    for site in cond2_passed:
        host_check_seqs[f"{site.site_id}_5p"] = site.seed_5p
        host_check_seqs[f"{site.site_id}_3p"] = site.seed_3p

    maps_to_host_set = _batch_check_maps_to_host(host_check_seqs, host_ref, workdir)
    log(f"  Host-mapping check: {len(maps_to_host_set)} clips map to host")

    validated_sites: list[InsertionSite] = []
    for site in cond2_passed:
        maps_5p = f"{site.site_id}_5p" in maps_to_host_set
        maps_3p = f"{site.site_id}_3p" in maps_to_host_set
        if site.junction_5p:
            site.junction_5p.maps_to_host = maps_5p
        if site.junction_3p:
            site.junction_3p.maps_to_host = maps_3p
        site.clips_not_in_host = not maps_5p and not maps_3p
        if not site.clips_not_in_host:
            continue
        site.is_validated = True
        site.confidence = "high"
        validated_sites.append(site)

    # Batch element-hit check for all validated + single-direction sites
    element_check_seqs: dict[str, str] = {}
    for site in validated_sites:
        element_check_seqs[f"{site.site_id}_5p"] = site.seed_5p
        element_check_seqs[f"{site.site_id}_3p"] = site.seed_3p
    for site in single_sites:
        clip_seq = site.seed_5p or site.seed_3p
        if clip_seq:
            element_check_seqs[f"{site.site_id}_clip"] = clip_seq

    if element_db and element_check_seqs:
        all_element_hits = _batch_check_element_hits(
            element_check_seqs, element_db, workdir, extra_dbs=extra_dbs)

        for site in validated_sites:
            hits_5p = all_element_hits.get(f"{site.site_id}_5p", [])
            hits_3p = all_element_hits.get(f"{site.site_id}_3p", [])
            if site.junction_5p:
                site.junction_5p.element_hits = hits_5p
            if site.junction_3p:
                site.junction_3p.element_hits = hits_3p
            site.has_element_hits = bool(hits_5p or hits_3p)

        for site in single_sites:
            clip_hits = all_element_hits.get(f"{site.site_id}_clip", [])
            if clip_hits:
                if site.junction_5p:
                    site.junction_5p.element_hits = clip_hits
                if site.junction_3p:
                    site.junction_3p.element_hits = clip_hits
                site.has_element_hits = True
                site.confidence = "medium"
                validated_sites.append(site)
                log(f"  {site.site_id}: single-direction but element hit "
                    f"({clip_hits[0]}) → medium confidence")

    for site in validated_sites:
        log(f"  {site.site_id}: {site.host_chr}:{site.pos_5p}-{site.pos_3p} "
            f"VALIDATED [{site.confidence}] "
            f"(5p={len(site.seed_5p)}bp, 3p={len(site.seed_3p)}bp, "
            f"element={'yes' if site.has_element_hits else 'no'})")

    # Sort: element-hit sites FIRST (primary), then paired > single, then read support
    # Element DB hits are the strongest signal for T-DNA insertion vs SV
    validated_sites.sort(key=lambda s: (
        0 if s.has_element_hits else 1,
        0 if s.confidence == "high" else (1 if s.confidence == "medium" else 2),
        -(s.junction_5p.n_reads if s.junction_5p else 0)
        - (s.junction_3p.n_reads if s.junction_3p else 0),
    ))

    log(f"  Validated: {len(validated_sites)} insertion sites "
        f"({sum(1 for s in validated_sites if s.has_element_hits)} with element hits)")

    return validated_sites


# ---------------------------------------------------------------------------
# Transgene-positive identification
# ---------------------------------------------------------------------------
# Strategy: Instead of asking "is this clip in the host?" (fails with wrong
# cultivar reference), ask "is this clip from a transgene?"
# Uses blastn-short against transgene_db (element DB + UniVec vectors).
# Any transgene hit = assembly target. No hit = skip.

@dataclass
class TierResult:
    """Classification result for one insertion site."""
    site_id: str
    chrom: str = ""
    pos: int = 0
    transgene_positive: bool = False
    clip_5p_len: int = 0
    clip_3p_len: int = 0
    hit_5p: str = ""            # best transgene_db hit for 5p clip
    hit_5p_identity: float = 0
    hit_5p_aln_len: int = 0
    hit_5p_source: str = ""     # "element_db" or "univec"
    hit_3p: str = ""
    hit_3p_identity: float = 0
    hit_3p_aln_len: int = 0
    hit_3p_source: str = ""


def _filter_host_endogenous(
    transgene_db: Path,
    host_ref: Path,
    tier_dir: Path,
    workdir: Path,
    threads: int,
) -> tuple[Path, set[str]]:
    """BLAST transgene_db vs host genome to remove host-derived entries.

    Returns (blast_db, exclude_ids): blast_db is the filtered DB path
    (or original if nothing excluded), exclude_ids are the removed entries.
    """
    blast_db = transgene_db
    exclude_ids: set[str] = set()
    exclude_details: dict[str, str] = {}

    # Ensure host BLAST DB exists
    host_blast_db_exists = (
        host_ref.with_suffix(".fa.ndb").exists()
        or Path(str(host_ref) + ".ndb").exists()
    )
    if not host_blast_db_exists:
        log("  Host BLAST DB not found — creating with makeblastdb...")
        result = subprocess.run(
            ["makeblastdb", "-in", str(host_ref), "-dbtype", "nucl"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        if result.returncode == 0:
            host_blast_db_exists = True
        else:
            log("  WARNING: makeblastdb failed, skipping host-endogenous exclusion")

    if not host_blast_db_exists:
        log("  WARNING: Host BLAST DB unavailable, skipping host-endogenous exclusion")
        return blast_db, exclude_ids

    host_hits_out = tier_dir / "host_endogenous_hits.tsv"
    log("  BLAST transgene_db vs host genome to find host-endogenous entries...")
    result = subprocess.run(
        ["blastn", "-task", "blastn",
         "-query", str(transgene_db), "-db", str(host_ref),
         "-outfmt", "6 qseqid qlen sseqid pident length qcovs",
         "-evalue", "1e-10", "-max_target_seqs", "1",
         "-num_threads", str(threads),
         "-out", str(host_hits_out)],
        stderr=subprocess.DEVNULL,
    )
    if result.returncode != 0:
        log(f"  WARNING: host-endogenous BLAST failed (rc={result.returncode})")
        return blast_db, exclude_ids

    if host_hits_out.exists():
        with open(host_hits_out) as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 6:
                    continue
                qseqid = cols[0]
                pident = float(cols[3])
                qcovs = float(cols[5])
                tier1 = pident >= HOST_ENDO_T1_PIDENT and qcovs >= HOST_ENDO_T1_QCOVS
                tier2 = pident >= HOST_ENDO_T2_PIDENT and qcovs >= HOST_ENDO_T2_QCOVS
                if tier1 or tier2:
                    tier_label = "T1" if tier1 else "T2"
                    exclude_ids.add(qseqid)
                    exclude_details[qseqid] = (
                        f"[{tier_label}] pident={pident:.1f}%, qcovs={qcovs:.0f}%"
                    )

    if exclude_ids:
        filtered_db = tier_dir / "transgene_db_clean.fa"
        n_total = 0
        n_excluded = 0
        with open(transgene_db) as fin, open(filtered_db, "w") as fout:
            write = True
            for line in fin:
                if line.startswith(">"):
                    n_total += 1
                    seq_id = line[1:].split()[0]
                    if seq_id in exclude_ids:
                        n_excluded += 1
                        write = False
                    else:
                        write = True
                        fout.write(line)
                elif write:
                    fout.write(line)
        log(f"  Host-endogenous exclusion: removed {n_excluded}/{n_total} entries")
        for eid in sorted(exclude_ids):
            log(f"    excluded: {eid} ({exclude_details[eid]})")
        subprocess.run(
            ["makeblastdb", "-in", str(filtered_db), "-dbtype", "nucl"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        blast_db = filtered_db
    else:
        log("  Host-endogenous exclusion: no entries matched host genome")

    # Write for debugging/auditability (Phase 4 receives exclude_ids directly)
    persist_path = workdir / "host_endogenous_ids.txt"
    with open(persist_path, "w") as fh:
        for eid in sorted(exclude_ids):
            fh.write(f"{eid}\t{exclude_details.get(eid, '')}\n")

    return blast_db, exclude_ids


def classify_site_tiers(
    sites: list[InsertionSite],
    element_db: Path,
    host_ref: Path,
    workdir: Path,
    threads: int = 4,
    min_identity: float = 80.0,
    min_aln_len: int = 20,
    extra_transgene_dbs: list[Path] | None = None,
) -> tuple[list[InsertionSite], list[InsertionSite], list[TierResult], set[str]]:
    """Classify sites by transgene-positive identification.

    BLASTs all clip sequences against transgene_db (element DB + UniVec)
    using blastn-short (optimized for 20-80bp queries).

    TRANSGENE-POSITIVE: at least one clip hits transgene_db → assemble
    TRANSGENE-NEGATIVE: no clip hits transgene_db → skip

    Returns (assembly_sites, skip_sites, all_tier_results, host_endo_ids).
    """
    if not sites:
        return [], [], [], set()

    tier_dir = workdir / "_tier_classification"
    tier_dir.mkdir(parents=True, exist_ok=True)

    # Locate transgene_db (element_db + UniVec combined)
    transgene_db = element_db.parent / "transgene_db.fa"
    if not transgene_db.exists():
        # Build transgene_db by combining element_db + univec
        univec_db = Path(__file__).resolve().parent.parent / "db" / "univec_vectors.fa"
        log("  Building transgene_db (element DB + UniVec)...")
        with open(transgene_db, "w") as out:
            if element_db.exists():
                out.write(element_db.read_text())
            if univec_db.exists():
                out.write(univec_db.read_text())
        # Build BLAST index
        subprocess.run(
            ["makeblastdb", "-in", str(transgene_db), "-dbtype", "nucl"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )

    # Ensure BLAST DB index exists
    if not transgene_db.with_suffix(".fa.ndb").exists() and \
       not Path(str(transgene_db) + ".ndb").exists():
        subprocess.run(
            ["makeblastdb", "-in", str(transgene_db), "-dbtype", "nucl"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )

    blast_db, exclude_ids = _filter_host_endogenous(
        transgene_db, host_ref, tier_dir, workdir, threads
    )

    # Step 1: Collect all clip sequences
    clip_seqs: dict[str, str] = {}
    for site in sites:
        if site.seed_5p and len(site.seed_5p) >= 15:
            clip_seqs[f"{site.site_id}_5p"] = site.seed_5p
        if site.seed_3p and len(site.seed_3p) >= 15:
            clip_seqs[f"{site.site_id}_3p"] = site.seed_3p

    log(f"  BLASTing {len(clip_seqs)} clip sequences against transgene_db "
        f"(blastn-short, {blast_db.name})...")

    # Step 2: Write clips and BLAST against transgene_db
    clip_fa = tier_dir / "all_clips.fa"
    with open(clip_fa, "w") as fh:
        for name, seq in clip_seqs.items():
            fh.write(f">{name}\n{seq}\n")

    blast_out = tier_dir / "transgene_blast.tsv"
    result = subprocess.run(
        ["blastn", "-task", "blastn-short",
         "-query", str(clip_fa), "-db", str(blast_db),
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
         "-evalue", "1e-5", "-max_target_seqs", "3",
         "-num_threads", str(threads),
         "-out", str(blast_out)],
        stderr=subprocess.DEVNULL,
    )
    if result.returncode != 0:
        log(f"  WARNING: transgene BLAST failed (rc={result.returncode})")

    # Step 3: Parse hits — best hit per query with identity/length filter
    hits: dict[str, dict] = {}
    if blast_out.exists():
        with open(blast_out) as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 6:
                    continue
                qname = cols[0]
                pident = float(cols[2])
                aln_len = int(cols[3])
                bitscore = float(cols[5])

                if pident < min_identity or aln_len < min_aln_len:
                    continue

                if qname not in hits or bitscore > hits[qname]["bitscore"]:
                    sseqid = cols[1]
                    source = "univec" if sseqid.startswith("univec|") else "element_db"
                    hits[qname] = {
                        "element": sseqid,
                        "identity": pident,
                        "aln_length": aln_len,
                        "bitscore": bitscore,
                        "source": source,
                    }

    log(f"  {len(hits)} clips hit transgene_db (identity>={min_identity}%, "
        f"aln>={min_aln_len}bp)")

    # Step 3b: Optional BLAST against extra transgene DBs. These include
    # the always-on shared payload DB (common_payload.fa) and per-sample
    # SPAdes contigs from s04b. Merged into `hits` so both shared and
    # sample-specific transgene references contribute to transgene-positive
    # classification.
    extra_dbs_list = [
        edb for edb in (extra_transgene_dbs or [])
        if edb is not None and edb.exists() and edb.stat().st_size > 0
    ]
    for i, edb in enumerate(extra_dbs_list):
        extra_blast_out = tier_dir / f"transgene_blast_extra_{i}_{edb.stem}.tsv"
        log(f"  BLASTing clips against extra transgene DB "
            f"({edb.name}, blastn-short)...")
        result_extra = subprocess.run(
            ["blastn", "-task", "blastn-short",
             "-query", str(clip_fa), "-subject", str(edb),
             "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
             "-evalue", "1e-5", "-max_target_seqs", "3",
             "-out", str(extra_blast_out)],
            stderr=subprocess.DEVNULL,
        )
        if result_extra.returncode != 0:
            log(f"  WARNING: extra transgene BLAST failed "
                f"(rc={result_extra.returncode})")
        n_extra = 0
        if extra_blast_out.exists():
            with open(extra_blast_out) as fh:
                for line in fh:
                    cols = line.strip().split("\t")
                    if len(cols) < 6:
                        continue
                    qname = cols[0]
                    pident = float(cols[2])
                    aln_len = int(cols[3])
                    bitscore = float(cols[5])
                    if pident < min_identity or aln_len < min_aln_len:
                        continue
                    if qname not in hits or bitscore > hits[qname]["bitscore"]:
                        hits[qname] = {
                            "element": cols[1],
                            "identity": pident,
                            "aln_length": aln_len,
                            "bitscore": bitscore,
                            "source": "element_db",
                        }
                        n_extra += 1
        log(f"  Extra transgene DB {edb.name} contributed/updated "
            f"{n_extra} clip hits (total now {len(hits)})")

    # Step 3.5: Per-clip host verification
    # Some divergent host elements (e.g., I-actin rice intron) escape DB-level
    # exclusion because their full-entry qcovs is low, but the SHORT clip we
    # match against them happens to be a 100% conserved window. To catch this,
    # BLAST each transgene-hit clip back against the host genome with
    # blastn-short. If the clip ALSO hits the host at >=95% identity over
    # >=30bp, it's an endogenous conserved sequence — drop it.
    host_blast_db_exists = (
        host_ref.with_suffix(".fa.ndb").exists()
        or Path(str(host_ref) + ".ndb").exists()
    )
    if hits and host_blast_db_exists:
        hit_clip_fa = tier_dir / "hit_clips.fa"
        with open(hit_clip_fa, "w") as fh:
            for qname in hits:
                fh.write(f">{qname}\n{clip_seqs[qname]}\n")

        host_clip_blast = tier_dir / "hit_clips_vs_host.tsv"
        log(f"  Per-clip host verification: BLAST {len(hits)} hit clips "
            f"vs host genome (blastn-short)...")
        result = subprocess.run(
            ["blastn", "-task", "blastn-short",
             "-query", str(hit_clip_fa), "-db", str(host_ref),
             "-outfmt", "6 qseqid sseqid pident length evalue",
             "-evalue", "1e-5", "-max_target_seqs", "1",
             "-num_threads", str(threads),
             "-out", str(host_clip_blast)],
            stderr=subprocess.DEVNULL,
        )
        if result.returncode != 0:
            log(f"  WARNING: per-clip host BLAST failed (rc={result.returncode})")

        host_endogenous_clips: dict[str, tuple[float, int]] = {}
        if host_clip_blast.exists():
            with open(host_clip_blast) as fh:
                for line in fh:
                    cols = line.strip().split("\t")
                    if len(cols) < 5:
                        continue
                    qname = cols[0]
                    pident = float(cols[2])
                    length = int(cols[3])
                    if pident >= CLIP_HOST_PIDENT and length >= CLIP_HOST_MIN_LEN:
                        # Keep best (longest) hit per clip
                        prev = host_endogenous_clips.get(qname)
                        if prev is None or length > prev[1]:
                            host_endogenous_clips[qname] = (pident, length)

        if host_endogenous_clips:
            log(f"  Per-clip host filter: dropping {len(host_endogenous_clips)} "
                f"clips that match host (pident>={CLIP_HOST_PIDENT}%, length>={CLIP_HOST_MIN_LEN}bp)")
            for qname, (pid, ln) in sorted(host_endogenous_clips.items()):
                hit = hits[qname]
                log(f"    dropped: {qname} → {hit['element']} "
                    f"(host hit: {pid:.0f}%/{ln}bp)")
                del hits[qname]
            log(f"  After host verification: {len(hits)} transgene-only clips remain")
        else:
            log("  Per-clip host filter: no clips matched host")

    # Step 4: Classify each site
    tier_results: list[TierResult] = []
    assembly_sites: list[InsertionSite] = []
    skip_sites: list[InsertionSite] = []

    for site in sites:
        key_5p = f"{site.site_id}_5p"
        key_3p = f"{site.site_id}_3p"

        hit_5p = hits.get(key_5p, {})
        hit_3p = hits.get(key_3p, {})

        # Require at least one element_db hit (not univec-only).
        # UniVec-only matches at this step are typically short (20-31bp)
        # alignments that match native plant DNA by chance. Real T-DNA
        # insertions always have at least one characteristic element
        # (promoter, selection marker, terminator) matching element_db.
        has_element_hit = (hit_5p.get("source") == "element_db") or \
                           (hit_3p.get("source") == "element_db")
        is_positive = has_element_hit

        tr = TierResult(
            site_id=site.site_id,
            chrom=site.host_chr,
            pos=site.pos_5p or site.pos_3p or 0,
            transgene_positive=is_positive,
            clip_5p_len=len(site.seed_5p) if site.seed_5p else 0,
            clip_3p_len=len(site.seed_3p) if site.seed_3p else 0,
            hit_5p=hit_5p.get("element", ""),
            hit_5p_identity=hit_5p.get("identity", 0),
            hit_5p_aln_len=hit_5p.get("aln_length", 0),
            hit_5p_source=hit_5p.get("source", ""),
            hit_3p=hit_3p.get("element", ""),
            hit_3p_identity=hit_3p.get("identity", 0),
            hit_3p_aln_len=hit_3p.get("aln_length", 0),
            hit_3p_source=hit_3p.get("source", ""),
        )
        tier_results.append(tr)

        if is_positive:
            assembly_sites.append(site)
        else:
            skip_sites.append(site)

    shutil.rmtree(tier_dir, ignore_errors=True)

    n_pos = len(assembly_sites)
    n_neg = len(skip_sites)
    log(f"  Transgene-positive (assemble): {n_pos}")
    log(f"  Transgene-negative (skip):     {n_neg}")

    # Log positive site details
    for tr in tier_results:
        if tr.transgene_positive:
            parts = []
            if tr.hit_5p:
                parts.append(f"5p={tr.hit_5p} ({tr.hit_5p_identity:.0f}%/{tr.hit_5p_aln_len}bp)")
            if tr.hit_3p:
                parts.append(f"3p={tr.hit_3p} ({tr.hit_3p_identity:.0f}%/{tr.hit_3p_aln_len}bp)")
            log(f"    {tr.site_id} {tr.chrom}:{tr.pos}: {', '.join(parts)}")

    return assembly_sites, skip_sites, tier_results, exclude_ids


def write_tier_classification(
    tier_results: list[TierResult],
    output_path: Path,
) -> None:
    """Write site_tier_classification.tsv."""
    with open(output_path, "w") as fh:
        fh.write("site_id\tchrom\tpos\ttransgene_positive\t"
                 "clip_5p_len\tclip_3p_len\t"
                 "hit_5p\thit_5p_identity\thit_5p_aln_len\thit_5p_source\t"
                 "hit_3p\thit_3p_identity\thit_3p_aln_len\thit_3p_source\n")
        for tr in tier_results:
            fh.write(f"{tr.site_id}\t{tr.chrom}\t{tr.pos}\t{tr.transgene_positive}\t"
                     f"{tr.clip_5p_len}\t{tr.clip_3p_len}\t"
                     f"{tr.hit_5p}\t{tr.hit_5p_identity}\t{tr.hit_5p_aln_len}\t{tr.hit_5p_source}\t"
                     f"{tr.hit_3p}\t{tr.hit_3p_identity}\t{tr.hit_3p_aln_len}\t{tr.hit_3p_source}\n")
    log(f"  Classification written: {output_path}")


# ---------------------------------------------------------------------------
# Fallback: Parse junctions.tsv from step 6
# ---------------------------------------------------------------------------

def parse_legacy_junctions(tsv_path: Path) -> list[LegacyJunction]:
    """Parse junctions.tsv from step 6 as fallback."""
    raw: list[LegacyJunction] = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            raw.append(LegacyJunction(
                contig_name=row["contig_name"],
                contig_len=int(row["contig_len"]),
                host_chr=row["host_chr"],
                host_start=int(row["host_start"]),
                host_end=int(row["host_end"]),
                host_strand=row["host_strand"],
                construct_element=row["construct_element"],
                construct_start=int(row["construct_start"]),
                construct_end=int(row["construct_end"]),
                junction_pos_host=int(row["junction_pos_host"]),
                junction_type=row["junction_type"],
                confidence=row["confidence"],
                host_mapq=int(row["host_mapq"]),
            ))

    # Deduplicate
    best: dict[tuple[str, str, int], LegacyJunction] = {}
    for j in raw:
        key = (j.contig_name, j.host_chr, j.junction_pos_host)
        if key not in best or j.host_mapq > best[key].host_mapq:
            best[key] = j
    return list(best.values())


def legacy_junctions_to_sites(
    junctions: list[LegacyJunction],
    host_bam: Path,
    min_clip: int = 15,
    window: int = 10,
) -> list[InsertionSite]:
    """Convert step 6 junctions to InsertionSites with soft-clip seeds."""
    if not junctions:
        return []

    # Group by chromosome
    by_chr: dict[str, list[LegacyJunction]] = {}
    for j in junctions:
        by_chr.setdefault(j.host_chr, []).append(j)

    sites: list[InsertionSite] = []
    used_ids: set[str] = set()

    for chrom, juncs in by_chr.items():
        juncs.sort(key=lambda j: j.junction_pos_host)

        # Pair within 50kb
        used = set()
        for i, j1 in enumerate(juncs):
            if i in used:
                continue
            paired_idx = None
            for k, j2 in enumerate(juncs):
                if k in used or k == i:
                    continue
                if abs(j2.junction_pos_host - j1.junction_pos_host) <= 50000:
                    paired_idx = k
                    break

            anchor_pos = j1.junction_pos_host
            sid = f"insertion_{chrom}_{anchor_pos}"
            if sid in used_ids:
                suffix = 2
                while f"{sid}_{suffix}" in used_ids:
                    suffix += 1
                sid = f"{sid}_{suffix}"
            used_ids.add(sid)
            site = InsertionSite(
                site_id=sid,
                host_chr=chrom,
            )

            positions = [j1.junction_pos_host]
            if paired_idx is not None:
                used.add(i)
                used.add(paired_idx)
                positions.append(juncs[paired_idx].junction_pos_host)
            else:
                used.add(i)

            site.pos_5p = min(positions)
            site.pos_3p = max(positions) if len(positions) > 1 else 0

            # Extract seeds from soft-clips at junction positions
            _extract_seeds_at_positions(site, host_bam, positions, chrom,
                                        min_clip, window)
            sites.append(site)

    return sites


def _extract_seeds_at_positions(
    site: InsertionSite,
    host_bam: Path,
    positions: list[int],
    chrom: str,
    min_clip: int,
    window: int,
) -> None:
    """Extract soft-clip seeds at known junction positions."""
    bam = pysam.AlignmentFile(str(host_bam), "rb")

    for jpos in positions:
        right_clips: list[str] = []
        left_clips: list[str] = []

        for read in bam.fetch(chrom, max(0, jpos - 200), jpos + 200):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            cigar = read.cigartuples
            if cigar is None:
                continue
            seq = read.query_sequence
            if seq is None:
                continue

            if cigar[-1][0] == 4 and cigar[-1][1] >= min_clip:
                if abs(read.reference_end - jpos) <= window:
                    right_clips.append(seq[-cigar[-1][1]:])

            if cigar[0][0] == 4 and cigar[0][1] >= min_clip:
                if abs(read.reference_start - jpos) <= window:
                    left_clips.append(seq[:cigar[0][1]])

        if right_clips and len(right_clips) >= 2:
            consensus = _build_consensus(right_clips, "right")
            if len(consensus) >= min_clip:
                site.junction_5p = JunctionCluster(
                    host_chr=chrom, position=jpos, clip_direction="right",
                    clipped_seqs=right_clips, consensus_clip=consensus,
                    n_reads=len(right_clips),
                )
                if not site.seed_5p or len(consensus) > len(site.seed_5p):
                    site.seed_5p = consensus

        if left_clips and len(left_clips) >= 2:
            consensus = _build_consensus(left_clips, "left")
            if len(consensus) >= min_clip:
                site.junction_3p = JunctionCluster(
                    host_chr=chrom, position=jpos, clip_direction="left",
                    clipped_seqs=left_clips, consensus_clip=consensus,
                    n_reads=len(left_clips),
                )
                if not site.seed_3p or len(consensus) > len(site.seed_3p):
                    site.seed_3p = consensus

    bam.close()


# ---------------------------------------------------------------------------
# Phase 2: Candidate Read Extraction
# ---------------------------------------------------------------------------

def extract_candidate_reads(
    host_bam: Path,
    site: InsertionSite,
    out_r1: Path,
    out_r2: Path,
    flank: int = 5000,
    threads: int = 4,
    min_clip: int = 20,
    junction_window: int = 20,
) -> int:
    """Extract reads for assembly from junction regions, phased by allele.

    For hemizygous insertions, mixing WT-allele and insertion-allele reads
    causes assembly ambiguity. We phase reads using soft-clip diagnostics:

    1. Scan junction position(s) ±junction_window
    2. INSERTION-allele: reads with soft-clip at junction position
       (these span the host→T-DNA boundary)
    3. WT-allele: reads that span the junction position WITHOUT a clip
       at that boundary (they read straight through → no insertion)
    4. AMBIGUOUS: distal reads (>junction_window from junction).
       These are kept ONLY IF their mate is not classified as WT.
       A distal read paired with a WT-classified mate is excluded
       (the WT evidence from the mate overrides the distal's ambiguity).
    5. Mate propagation: if any read in a pair is insertion-allele,
       the whole pair is insertion-allele (insertion wins ties).
    6. Output: insertion-allele + ambiguous-without-WT-mate pairs,
       excluding pairs where any mate has WT evidence and no mate
       has insertion evidence (exclude = wt_reads - insertion_reads).
    """
    tmp_dir = out_r1.parent

    # Build regions
    regions = []
    positions = []
    if site.pos_5p > 0:
        positions.append(site.pos_5p)
    if site.pos_3p > 0:
        positions.append(site.pos_3p)
    if not positions:
        return 0

    for pos in positions:
        start = max(1, pos - flank)
        end = pos + flank
        regions.append(f"{site.host_chr}:{start}-{end}")

    # Extract reads from junction regions
    region_bam = tmp_dir / f"_{site.site_id}_regions.bam"
    with open(region_bam, "wb") as bam_fh:
        subprocess.run(
            ["samtools", "view", "-b", "-@", str(threads),
             str(host_bam)] + regions,
            stdout=bam_fh,
            stderr=subprocess.DEVNULL, check=True,
        )

    # Phase reads by allele using pysam
    insertion_reads: set[str] = set()
    wt_reads: set[str] = set()
    distal_reads: set[str] = set()

    with pysam.AlignmentFile(str(region_bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            if qname is None:
                continue
            cigar = read.cigartuples
            if not cigar:
                distal_reads.add(qname)
                continue

            ref_start = read.reference_start  # 0-based
            ref_end = read.reference_end       # exclusive
            if ref_start is None or ref_end is None:
                distal_reads.add(qname)
                continue

            classified = False
            for pos in positions:
                pos0 = pos - 1  # convert to 0-based
                # Distance check: is read NEAR this junction at all?
                if ref_end < pos0 - junction_window or ref_start > pos0 + junction_window:
                    continue

                # INSERTION diagnostic: soft-clip end touches junction
                left_clip = cigar[0][1] if cigar[0][0] == 4 else 0
                right_clip = cigar[-1][1] if cigar[-1][0] == 4 else 0

                # Right soft-clip → break point is at ref_end → matches 5p junction
                is_insertion = False
                if right_clip >= min_clip and abs(ref_end - pos0) <= junction_window:
                    is_insertion = True
                # Left soft-clip → break point is at ref_start → matches 3p junction
                if left_clip >= min_clip and abs(ref_start - pos0) <= junction_window:
                    is_insertion = True

                if is_insertion:
                    insertion_reads.add(qname)
                    classified = True
                    break

                # WT diagnostic: spans junction (start before & end after) WITHOUT clip
                # at that boundary
                if ref_start <= pos0 - junction_window and ref_end >= pos0 + junction_window:
                    # No clip at the junction boundary
                    if right_clip < min_clip and left_clip < min_clip:
                        wt_reads.add(qname)
                        classified = True
                        break

            if not classified:
                distal_reads.add(qname)

    # Mate propagation: insertion wins over WT.
    # samtools view -N fetches BOTH mates by name, so we only need a name set.
    # A pair is excluded if ANY mate has WT evidence (spans junction w/o clip)
    # AND no mate is insertion-classified. Distal does NOT rescue WT: the
    # 3' mate of a pair is ~300bp away from the junction and gets classified
    # as distal, but that doesn't contradict the 5' mate's WT observation.
    all_names = insertion_reads | wt_reads | distal_reads
    exclude_names = wt_reads - insertion_reads
    keep_names = all_names - exclude_names

    log(f"    Phasing: {len(insertion_reads):,} insertion-allele, "
        f"{len(wt_reads):,} WT-allele, {len(distal_reads):,} distal/ambiguous")
    log(f"    Excluded {len(exclude_names):,} pure-WT pairs")

    namelist = tmp_dir / f"_{site.site_id}_names.txt"
    with open(namelist, "w") as fh:
        for n in sorted(keep_names):
            fh.write(n + "\n")
    log(f"    Junction region reads (phased): {len(keep_names):,} read names")

    # Extract both mates from full BAM
    both_mates_bam = tmp_dir / f"_{site.site_id}_mates.bam"
    with open(both_mates_bam, "wb") as bam_fh:
        subprocess.run(
            ["samtools", "view", "-b", "-N", str(namelist),
             "-@", str(threads), str(host_bam)],
            stdout=bam_fh,
            stderr=subprocess.DEVNULL, check=True,
        )

    # Also get all-unmapped pairs (both mates unmapped, flag 12)
    unmapped_bam = tmp_dir / f"_{site.site_id}_unmapped.bam"
    with open(unmapped_bam, "wb") as bam_fh:
        subprocess.run(
            ["samtools", "view", "-b", "-f", "12", "-@", str(threads), str(host_bam)],
            stdout=bam_fh,
            stderr=subprocess.DEVNULL, check=True,
        )

    # Merge
    merged_bam = tmp_dir / f"_{site.site_id}_merged.bam"
    subprocess.run(
        ["samtools", "merge", "-f", "-@", str(threads), str(merged_bam),
         str(both_mates_bam), str(unmapped_bam)],
        stderr=subprocess.DEVNULL, check=True,
    )

    # Name-sort and convert to paired FASTQ
    nsort_bam = tmp_dir / f"_{site.site_id}_nsort.bam"
    subprocess.run(
        ["samtools", "sort", "-n", "-@", str(threads),
         str(merged_bam), "-o", str(nsort_bam)],
        stderr=subprocess.DEVNULL, check=True,
    )
    subprocess.run(
        ["samtools", "fastq", "-@", str(threads),
         "-1", str(out_r1), "-2", str(out_r2),
         "-s", "/dev/null", "-0", "/dev/null", str(nsort_bam)],
        stderr=subprocess.DEVNULL, check=True,
    )

    # Count reads
    count = 0
    opener = gzip.open if str(out_r1).endswith(".gz") else open
    with opener(str(out_r1), "rt") as fh:
        for _ in fh:
            count += 1
    count = count // 4  # FASTQ: 4 lines per read

    # Cleanup temp files
    for p in [region_bam, namelist, both_mates_bam, unmapped_bam,
              merged_bam, nsort_bam]:
        p.unlink(missing_ok=True)

    return count


def extract_unmapped_paired(
    host_bam: Path,
    workdir: Path,
    threads: int = 4,
) -> tuple[Path, Path]:
    """Extract unmapped reads from host BAM as paired FASTQ (cached)."""
    r1 = workdir / "_unmapped_R1.fq.gz"
    r2 = workdir / "_unmapped_R2.fq.gz"
    if r1.exists() and r2.exists():
        n = 0
        with gzip.open(r1, "rt") as fh:
            for _ in fh:
                n += 1
        log(f"  Unmapped reads (cached): {n // 4:,} pairs")
        return r1, r2

    log(f"  Extracting unmapped reads from host BAM...")
    subprocess.run(
        f"samtools view -f 4 -@ {threads} -b {host_bam} "
        f"| samtools sort -n -@ {threads} - "
        f"| samtools fastq -1 {r1} -2 {r2} -s /dev/null -0 /dev/null - "
        f"2>/dev/null",
        shell=True, check=True,
    )
    n = 0
    with gzip.open(r1, "rt") as fh:
        for _ in fh:
            n += 1
    log(f"  Unmapped reads: {n // 4:,} pairs")
    return r1, r2


# ---------------------------------------------------------------------------
# Phase 3: K-mer Extension + Pilon Gap Fill
# ---------------------------------------------------------------------------

class StrandAwareSeedExtender:
    """K-mer based extension engine with PE strand awareness.

    INDEX RULES:
    - R1 sequences: indexed as-is (forward strand)
    - R2 sequences: reverse-complemented then indexed (fragment strand)
    - Raw R2 (as sequenced): NEVER indexed, NEVER used for extension

    EXTENSION RULES:
    - Track used reads to prevent infinite loops
    - Stop at branch points (base ratio < min_ratio)
    - Majority vote consensus for each new base
    """

    def __init__(
        self,
        k: int = 15,
        min_overlap: int = 20,
        min_depth: int = 2,
        min_ratio: float = 0.7,
    ):
        self.k = k
        self.min_overlap = min_overlap
        self.min_depth = min_depth
        self.min_ratio = min_ratio
        self.seqs: list[str] = []
        self.kmer_index: dict[str, list[tuple[int, int]]] = defaultdict(list)

    def load_paired_reads(self, r1_path: Path, r2_path: Path) -> int:
        """Load R1 as forward, R2 as reverse-complement. Return pair count."""
        r1_seqs = _read_fq_seqs(r1_path)
        r2_seqs = _read_fq_seqs(r2_path)
        k = self.k
        n = 0
        for r1, r2 in zip(r1_seqs, r2_seqs):
            if len(r1) < self.min_overlap or len(r2) < self.min_overlap:
                continue
            if "N" in r1 or "N" in r2:
                continue
            r2_rc = revcomp(r2)

            for seq in (r1, r2_rc):
                idx = len(self.seqs)
                self.seqs.append(seq)
                for p in range(len(seq) - k + 1):
                    self.kmer_index[seq[p:p + k]].append((idx, p))
            n += 1
        return n

    def add_seqs(self, seqs: list[str]) -> int:
        """Add pre-processed sequences to pool and index."""
        existing = set(self.seqs)
        k = self.k
        n_new = 0
        for seq in seqs:
            seq = seq.upper()
            if seq in existing or len(seq) < self.min_overlap or "N" in seq:
                continue
            existing.add(seq)
            idx = len(self.seqs)
            self.seqs.append(seq)
            for p in range(len(seq) - k + 1):
                self.kmer_index[seq[p:p + k]].append((idx, p))
            n_new += 1
        return n_new

    def _extend_right(self, seed: str, used: set[int]) -> str:
        """Extend seed to the right with used-read tracking."""
        k = self.k
        min_ovl = self.min_overlap
        current = seed

        while True:
            if len(current) < k:
                break
            tail_kmer = current[-k:]
            hits = self.kmer_index.get(tail_kmer)
            if not hits:
                break

            votes: dict[int, Counter] = defaultdict(Counter)
            for seq_idx, kpos in hits:
                if seq_idx in used:
                    continue
                overlap = kpos + k
                if overlap < min_ovl:
                    continue
                read_seq = self.seqs[seq_idx]
                if current[-overlap:] == read_seq[:overlap]:
                    for j in range(overlap, len(read_seq)):
                        votes[j - overlap][read_seq[j]] += 1
                    used.add(seq_idx)

            if not votes:
                break

            ext = []
            for p in range(max(votes.keys()) + 1):
                if p not in votes:
                    break
                v = votes[p]
                total = sum(v.values())
                if total < self.min_depth:
                    break
                best_base, best_count = v.most_common(1)[0]
                if best_count / total < self.min_ratio:
                    break
                ext.append(best_base)

            if not ext:
                break
            current = current + "".join(ext)

        return current

    def extend(self, seed: str, max_iterations: int = 100) -> str:
        """Extend seed bidirectionally. Return extended sequence."""
        current = seed
        used: set[int] = set()
        for i in range(max_iterations):
            # Right extension
            extended = self._extend_right(current, used)
            # Left extension (revcomp -> extend right -> revcomp)
            rc_extended = self._extend_right(revcomp(extended), used)
            extended = revcomp(rc_extended)

            growth = len(extended) - len(current)
            if growth == 0:
                log(f"      Iter {i + 1}: {len(current):,}bp → converged")
                break
            log(f"      Iter {i + 1}: {len(current):,} → {len(extended):,}bp "
                f"(+{growth})")
            current = extended
        return current


# ---------------------------------------------------------------------------
# K-mer recruitment
# ---------------------------------------------------------------------------

def recruit_by_kmer(
    contig_seq: str,
    unmapped_r1: Path,
    unmapped_r2: Path,
    k: int = 25,
) -> tuple[list[str], list[str]]:
    """K-mer based paired recruitment from unmapped reads.

    Builds k-mer set from contig (both strands), scans R1 forward + R2 RC.
    Raw R2 is NOT used (head-to-head safety).
    Returns (r1_seqs, r2_rc_seqs).
    """
    contig_kmers: set[str] = set()
    for seq in (contig_seq, revcomp(contig_seq)):
        for i in range(len(seq) - k + 1):
            contig_kmers.add(seq[i:i + k])

    r1_seqs = _read_fq_seqs(unmapped_r1)
    r2_seqs = _read_fq_seqs(unmapped_r2)

    recruited_r1: list[str] = []
    recruited_r2rc: list[str] = []
    stride = max(1, k // 2)

    for r1, r2 in zip(r1_seqs, r2_seqs):
        if "N" in r1 or "N" in r2:
            continue
        r2_rc = revcomp(r2)

        hit = False
        for seq in (r1, r2_rc):
            for i in range(0, len(seq) - k + 1, stride):
                if seq[i:i + k] in contig_kmers:
                    hit = True
                    break
            if hit:
                break

        if hit:
            recruited_r1.append(r1)
            recruited_r2rc.append(r2_rc)

    return recruited_r1, recruited_r2rc


# ---------------------------------------------------------------------------
# Pilon gap fill (CRITICAL: [contig_5p] + NNN + [contig_3p] as ONE scaffold)
# ---------------------------------------------------------------------------

def pilon_fill(
    contig_5p: str,
    contig_3p: str,
    r1_fq: Path,
    r2_fq: Path,
    workdir: Path,
    gap_size: int = 1000,
    threads: int = 4,
) -> tuple[str, str, bool]:
    """Build scaffold [contig_5p]+NNN+[contig_3p] and run Pilon to fill gap.

    KEY INSIGHT: Pilon only fills INTERNAL gaps between two known sequences.
    Both contigs must be in ONE scaffold for Pilon to work.

    Returns (updated_5p, updated_3p, gap_filled).
    gap_filled=True means no N remains (complete assembly).
    """
    workdir.mkdir(parents=True, exist_ok=True)

    if not contig_5p or not contig_3p:
        return contig_5p, contig_3p, False

    # Build scaffold: [contig_5p] + NNN + [contig_3p]
    scaffold = contig_5p + ("N" * gap_size) + contig_3p
    ref_fa = workdir / "pilon_ref.fasta"
    write_fasta(ref_fa, "scaffold", scaffold)

    # Map reads with minimap2
    bam = workdir / "pilon_mapped.bam"
    subprocess.run(
        f"minimap2 -ax sr -t {threads} --secondary=no "
        f"{ref_fa} {r1_fq} {r2_fq} 2>/dev/null "
        f"| samtools sort -@ {threads} -o {bam}",
        shell=True, check=True,
    )
    subprocess.run(["samtools", "index", str(bam)],
                   stderr=subprocess.DEVNULL, check=True)

    # Count mapped reads
    r = subprocess.run(
        ["samtools", "view", "-c", "-F", "4", str(bam)],
        stdout=subprocess.PIPE, text=True, stderr=subprocess.DEVNULL,
    )
    n_mapped = int(r.stdout.strip()) if r.stdout.strip() else 0

    if n_mapped < 5:
        log(f"    Pilon: only {n_mapped} mapped reads, skipping")
        return contig_5p, contig_3p, False

    log(f"    Pilon: {n_mapped:,} reads mapped to scaffold "
        f"({len(contig_5p)}+{gap_size}N+{len(contig_3p)}={len(scaffold)}bp)")

    # Run Pilon
    pilon_prefix = workdir / "pilon_out"
    pilon_log_file = workdir / "pilon.log"
    with open(pilon_log_file, "w") as logfh:
        subprocess.run(
            ["pilon", "--genome", str(ref_fa), "--frags", str(bam),
             "--output", str(pilon_prefix), "--fix", "all",
             "--mindepth", "2", "--gapmargin", "100000"],
            stdout=logfh, stderr=subprocess.STDOUT,
        )

    if pilon_log_file.exists():
        with open(pilon_log_file) as fh:
            for line in fh:
                line = line.strip()
                if any(k in line for k in ["Gap", "fix", "Total", "Confirmed"]):
                    log(f"    [Pilon] {line}")

    pilon_fa = Path(f"{pilon_prefix}.fasta")
    if not pilon_fa.exists():
        return contig_5p, contig_3p, False

    seqs = read_fasta(pilon_fa)
    if not seqs:
        return contig_5p, contig_3p, False

    filled = list(seqs.values())[0]

    # Check if gap is fully filled (no N remaining)
    n_match = re.search(r"N{10,}", filled)
    if n_match is None:
        # Gap fully filled!
        log(f"    Pilon: gap completely filled! ({len(filled):,}bp)")
        return filled, "", True

    # Split at longest remaining N stretch
    n_runs = [(m.start(), m.end()) for m in re.finditer(r"N+", filled)]
    if not n_runs:
        return filled, "", True

    # Find longest N run
    longest = max(n_runs, key=lambda x: x[1] - x[0])
    new_5p = filled[:longest[0]]
    new_3p = filled[longest[1]:]

    growth_5p = len(new_5p) - len(contig_5p)
    growth_3p = len(new_3p) - len(contig_3p)
    remaining_gap = longest[1] - longest[0]
    log(f"    Pilon: 5p {len(contig_5p):,}→{len(new_5p):,}bp (+{growth_5p}), "
        f"3p {len(contig_3p):,}→{len(new_3p):,}bp (+{growth_3p}), "
        f"gap {gap_size}→{remaining_gap}N")

    return new_5p, new_3p, False


# ---------------------------------------------------------------------------
# Host genome termination check
# ---------------------------------------------------------------------------

def check_host_termination(
    contig_5p: str,
    contig_3p: str,
    host_ref: Path,
    workdir: Path,
    check_len: int = 100,
    min_identity: float = 0.95,
    min_match: int = 50,
) -> tuple[bool, bool]:
    """Check if contig ends have reached back to host genome.

    5' contig growing end = right end (extending into insert then through to host)
    3' contig growing end = left end (extending into insert then through to host)

    Returns (reached_5p, reached_3p).
    """
    if len(contig_5p) < check_len * 2 and len(contig_3p) < check_len * 2:
        return False, False

    query_fa = workdir / "_host_term.fa"
    with open(query_fa, "w") as fh:
        if len(contig_5p) >= check_len * 2:
            fh.write(f">growing_5p\n{contig_5p[-check_len:]}\n")
        if len(contig_3p) >= check_len * 2:
            fh.write(f">growing_3p\n{contig_3p[:check_len]}\n")

    paf = workdir / "_host_term.paf"
    subprocess.run(
        ["minimap2", "-c", "--secondary=no", "-t", "1",
         str(host_ref), str(query_fa)],
        stdout=open(paf, "w"), stderr=subprocess.DEVNULL,
    )

    reached_5p = False
    reached_3p = False
    if paf.exists():
        with open(paf) as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 12:
                    continue
                qname = cols[0]
                match_bp = int(cols[9])
                block_len = int(cols[10])
                if match_bp < min_match or block_len < min_match:
                    continue
                identity = match_bp / block_len
                if identity < min_identity:
                    continue
                if qname == "growing_5p":
                    reached_5p = True
                elif qname == "growing_3p":
                    reached_3p = True

    for f in (query_fa, paf):
        f.unlink(missing_ok=True)

    return reached_5p, reached_3p


# ---------------------------------------------------------------------------
# Foreign read refinement (minimap2 chaining + Pilon)
# ---------------------------------------------------------------------------

def extract_foreign_reads(
    s03_r1: Path,
    s03_r2: Path,
    host_ref: Path,
    workdir: Path,
    threads: int = 4,
) -> tuple[Path, Path]:
    """Extract s03 reads that DON'T map to host genome (construct-only reads).

    These 'foreign' reads contain transgene sequences not present in the host.
    By filtering out host-mapping reads, we enrich for construct-internal
    elements (hLF1, G6 EPSPS, Pepc, etc.) that enable Pilon to extend
    assemblies through regions the k-mer assembler couldn't resolve.
    """
    foreign_r1 = workdir / "_foreign_R1.fq.gz"
    foreign_r2 = workdir / "_foreign_R2.fq.gz"
    if foreign_r1.exists() and foreign_r2.exists():
        n = 0
        with gzip.open(foreign_r1, "rt") as fh:
            for _ in fh:
                n += 1
        log(f"  Foreign reads (cached): {n // 4:,} pairs")
        return foreign_r1, foreign_r2

    log(f"  Extracting foreign reads (s03 reads unmapped to host)...")

    # Map s03 reads to host
    host_bam = workdir / "_s03_to_host.bam"
    subprocess.run(
        f"bwa mem -t {threads} {host_ref} {s03_r1} {s03_r2} 2>/dev/null "
        f"| samtools sort -@ 4 -o {host_bam}",
        shell=True, check=True,
    )
    subprocess.run(["samtools", "index", str(host_bam)],
                   stderr=subprocess.DEVNULL, check=True)

    # Get properly mapped read names (MAPQ >= 20, both reads mapped)
    host_mapped = set()
    all_reads = set()
    with pysam.AlignmentFile(str(host_bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            all_reads.add(read.query_name)
            if (not read.is_unmapped and not read.mate_is_unmapped
                    and read.mapping_quality >= 20):
                host_mapped.add(read.query_name)

    foreign_names = all_reads - host_mapped
    log(f"  Total: {len(all_reads):,} pairs, "
        f"host: {len(host_mapped):,}, foreign: {len(foreign_names):,}")

    # Extract foreign reads as FASTQ
    readname_file = workdir / "_foreign_names.txt"
    with open(readname_file, "w") as fh:
        for name in foreign_names:
            fh.write(f"{name}\n")

    subprocess.run(
        f"samtools view -h -N {readname_file} {host_bam} "
        f"| samtools sort -n -@ 4 - "
        f"| samtools fastq -1 {foreign_r1} -2 {foreign_r2} "
        f"-s /dev/null -0 /dev/null - 2>/dev/null",
        shell=True, check=True,
    )

    # Cleanup
    host_bam.unlink(missing_ok=True)
    Path(f"{host_bam}.bai").unlink(missing_ok=True)
    readname_file.unlink(missing_ok=True)

    n = 0
    with gzip.open(foreign_r1, "rt") as fh:
        for _ in fh:
            n += 1
    log(f"  Foreign reads extracted: {n // 4:,} pairs")
    return foreign_r1, foreign_r2


def refine_with_foreign_reads(
    insert_fasta: Path,
    s03_r1: Path,
    s03_r2: Path,
    host_ref: Path,
    workdir: Path,
    threads: int = 4,
    max_rounds: int = 10,
) -> Path:
    """Refine assembled insert using minimap2 chaining + Pilon with foreign reads.

    Approach:
      1. Extract foreign reads (s03 reads unmapped to host)
      2. Iterative: minimap2 map → Pilon --fix all → check convergence
      3. Pilon's local reassembly can extend construct regions using
         foreign reads that the k-mer assembler couldn't resolve
         (especially palindromic/inverted repeat regions in head-to-head T-DNA)

    Returns path to refined FASTA.
    """
    refine_dir = workdir / "_foreign_refine"
    refine_dir.mkdir(parents=True, exist_ok=True)

    # Extract foreign reads
    foreign_r1, foreign_r2 = extract_foreign_reads(
        s03_r1, s03_r2, host_ref, refine_dir, threads=threads,
    )

    # Check if we have enough foreign reads
    n_foreign = 0
    with gzip.open(foreign_r1, "rt") as fh:
        for _ in fh:
            n_foreign += 1
    n_foreign //= 4
    if n_foreign < 10:
        log(f"  Too few foreign reads ({n_foreign}), skipping refinement")
        return insert_fasta

    # Start with current assembly
    seqs = read_fasta(insert_fasta)
    if not seqs:
        return insert_fasta
    seq_name = list(seqs.keys())[0]
    current_seq = list(seqs.values())[0]
    prev_len = len(current_seq)

    log(f"  Foreign read refinement: {prev_len:,}bp scaffold, "
        f"{n_foreign:,} foreign read pairs")

    for rnd in range(1, max_rounds + 1):
        scaffold_fa = refine_dir / f"scaffold_r{rnd}.fa"
        write_fasta(scaffold_fa, "insert_scaffold", current_seq)

        # Map foreign reads with minimap2 (k-mer chaining handles repeats)
        bam = refine_dir / f"r{rnd}.bam"
        subprocess.run(
            f"minimap2 -ax sr -t {threads} {scaffold_fa} "
            f"{foreign_r1} {foreign_r2} 2>/dev/null "
            f"| samtools sort -@ 4 -o {bam}",
            shell=True, check=True,
        )
        subprocess.run(["samtools", "index", str(bam)],
                       stderr=subprocess.DEVNULL, check=True)

        # Count mapped
        r = subprocess.run(
            ["samtools", "view", "-c", "-F", "4", str(bam)],
            stdout=subprocess.PIPE, text=True, stderr=subprocess.DEVNULL,
        )
        n_mapped = int(r.stdout.strip()) if r.stdout.strip() else 0

        if n_mapped < 5:
            log(f"    Round {rnd}: only {n_mapped} mapped, stopping")
            break

        # Run Pilon with --fix all
        pilon_prefix = refine_dir / f"pilon_r{rnd}"
        pilon_log = refine_dir / f"pilon_r{rnd}.log"
        with open(pilon_log, "w") as logfh:
            subprocess.run(
                ["pilon", "--genome", str(scaffold_fa),
                 "--frags", str(bam),
                 "--output", str(pilon_prefix),
                 "--fix", "all", "--mindepth", "1"],
                stdout=logfh, stderr=subprocess.STDOUT,
            )

        pilon_fa = Path(f"{pilon_prefix}.fasta")
        if not pilon_fa.exists():
            log(f"    Round {rnd}: Pilon failed")
            break

        new_seqs = read_fasta(pilon_fa)
        if not new_seqs:
            break

        new_seq = list(new_seqs.values())[0]
        new_len = len(new_seq)
        n_ns = new_seq.upper().count("N")

        log(f"    Round {rnd}: {prev_len:,}→{new_len:,}bp "
            f"(Ns={n_ns}, mapped={n_mapped})")

        # Check convergence
        if new_len == prev_len and new_seq == current_seq:
            log(f"    Converged at round {rnd}")
            break

        current_seq = new_seq
        prev_len = new_len

        # Cleanup intermediate files
        scaffold_fa.unlink(missing_ok=True)
        bam.unlink(missing_ok=True)
        Path(f"{bam}.bai").unlink(missing_ok=True)
        pilon_fa.unlink(missing_ok=True)
        pilon_log.unlink(missing_ok=True)

    # Write refined result
    refined_fa = workdir / insert_fasta.name
    write_fasta(refined_fa, seq_name, current_seq)

    orig_seqs = read_fasta(insert_fasta)
    orig_len = len(list(orig_seqs.values())[0])
    final_len = len(current_seq)
    final_ns = current_seq.upper().count("N")
    log(f"  Refinement: {orig_len:,}→{final_len:,}bp, "
        f"Ns: {list(orig_seqs.values())[0].upper().count('N')}→{final_ns}")

    # Cleanup refine directory
    shutil.rmtree(refine_dir, ignore_errors=True)

    # Update the original file
    write_fasta(insert_fasta, seq_name, current_seq)
    return insert_fasta


# ---------------------------------------------------------------------------
# Overlap / merge detection
# ---------------------------------------------------------------------------

def _minimap2_extend(
    contig: str,
    r1_fq: Path,
    r2_fq: Path,
    workdir: Path,
    threads: int = 4,
    min_clip: int = 20,
    min_depth: int = 2,
    min_ratio: float = 0.6,
) -> str:
    """Extend contig using minimap2 soft-clip consensus.

    Maps reads to contig, collects soft-clipped tails beyond contig ends,
    and extends by majority vote. Handles repeats that k-mer extension
    cannot resolve because minimap2 uses minimizer chaining with gap
    penalties for placement, not exact k-mer lookup.

    Returns extended contig (may be unchanged if no extension possible).
    """
    if len(contig) < 100:
        return contig

    ref_fa = workdir / "_mm2ext_ref.fa"
    write_fasta(ref_fa, "contig", contig)

    bam_path = workdir / "_mm2ext.bam"
    subprocess.run(
        f"minimap2 -ax sr -t {threads} --secondary=no "
        f"{ref_fa} {r1_fq} {r2_fq} 2>/dev/null "
        f"| samtools sort -@ 2 -o {bam_path}",
        shell=True, check=True,
    )
    subprocess.run(["samtools", "index", str(bam_path)],
                   stderr=subprocess.DEVNULL, check=True)

    contig_len = len(contig)

    # Collect right-side soft-clip extensions (reads extending past contig end)
    right_tails: list[str] = []
    # Collect left-side soft-clip extensions (reads extending before contig start)
    left_tails: list[str] = []

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            cigar = read.cigartuples
            if not cigar:
                continue
            seq = read.query_sequence
            if not seq:
                continue

            # Right extension: read maps near contig end, has right soft-clip
            # cigar[-1] == (4, clip_len) means soft-clip at right end
            if (cigar[-1][0] == 4 and cigar[-1][1] >= min_clip
                    and read.reference_end is not None
                    and read.reference_end >= contig_len - 10):
                clip_len = cigar[-1][1]
                tail = seq[-clip_len:]
                right_tails.append(tail)

            # Left extension: read maps near contig start, has left soft-clip
            # cigar[0] == (4, clip_len) means soft-clip at left end
            if (cigar[0][0] == 4 and cigar[0][1] >= min_clip
                    and read.reference_start is not None
                    and read.reference_start <= 10):
                clip_len = cigar[0][1]
                tail = seq[:clip_len]
                left_tails.append(tail)

    # Cleanup
    for f in [ref_fa, bam_path, Path(f"{bam_path}.bai")]:
        f.unlink(missing_ok=True)

    extended = contig

    # Right extension via majority vote on soft-clip tails
    if len(right_tails) >= min_depth:
        ext = _vote_extension(right_tails, min_depth, min_ratio)
        if ext:
            extended = extended + ext

    # Left extension via majority vote on soft-clip tails (reverse direction)
    if len(left_tails) >= min_depth:
        # Reverse tails so we vote from the clip boundary outward
        rev_tails = [t[::-1] for t in left_tails]
        ext = _vote_extension(rev_tails, min_depth, min_ratio)
        if ext:
            extended = ext[::-1] + extended

    return extended


def _vote_extension(
    tails: list[str], min_depth: int = 2, min_ratio: float = 0.6,
) -> str:
    """Majority-vote consensus from a list of tail sequences.

    Each tail starts at the extension boundary and grows outward.
    Returns the consensus extension string (may be empty).
    """
    ext_bases: list[str] = []
    pos = 0
    while True:
        votes: Counter = Counter()
        for t in tails:
            if pos < len(t):
                votes[t[pos]] += 1
        total = sum(votes.values())
        if total < min_depth:
            break
        best_base, best_count = votes.most_common(1)[0]
        if best_count / total < min_ratio:
            break
        ext_bases.append(best_base)
        pos += 1
    return "".join(ext_bases)


def _ssake_extend(
    contig: str,
    r1_fq: Path,
    r2_fq: Path,
    workdir: Path,
    min_overlap: int = 20,
) -> str:
    """Extend contig using SSAKE's overlap-layout-consensus.

    SSAKE uses a prefix tree for rapid read overlap detection,
    complementing k-mer extension (greedy, exact-overlap) and
    minimap2 extension (chain-based, soft-clip). SSAKE can resolve
    moderately repetitive regions where k-mer extension stalls due
    to its all-at-once prefix search over the full read set.

    Runs via 'micromamba run -n ssake SSAKE' (separate Python 2.7 env).
    Returns extended contig (unchanged if SSAKE fails or no growth).
    """
    if len(contig) < 50:
        return contig

    ssake_dir = workdir / "_ssake"
    ssake_dir.mkdir(parents=True, exist_ok=True)

    # SSAKE input: seed FASTA (-s) + reads FASTA (-f)
    seed_fa = ssake_dir / "seed.fa"
    with open(seed_fa, "w") as fh:
        fh.write(f">seed\n{contig}\n")

    # Convert paired FASTQ to interleaved FASTA for SSAKE
    reads_fa = ssake_dir / "reads.fa"
    n_reads = 0
    with open(reads_fa, "w") as out:
        for fq_path in (r1_fq, r2_fq):
            if not fq_path.exists():
                continue
            with open(fq_path) as fq:
                while True:
                    header = fq.readline()
                    if not header:
                        break
                    seq = fq.readline().strip()
                    fq.readline()  # +
                    fq.readline()  # qual
                    if seq and len(seq) >= min_overlap:
                        out.write(f">r{n_reads}\n{seq}\n")
                        n_reads += 1

    if n_reads < 10:
        shutil.rmtree(ssake_dir, ignore_errors=True)
        return contig

    # Run SSAKE: -s seed, -f reads, -m min_overlap, -o 1 (1 pass), -w 1
    ssake_prefix = ssake_dir / "ssake_out"
    try:
        subprocess.run(
            ["micromamba", "run", "-n", "ssake", "SSAKE",
             "-f", str(reads_fa),
             "-s", str(seed_fa),
             "-b", str(ssake_prefix),
             "-m", str(min_overlap),
             "-o", "1",        # single pass
             "-w", "1",        # min reads for base call
             ],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            timeout=120,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError):
        shutil.rmtree(ssake_dir, ignore_errors=True)
        return contig

    # Read SSAKE contigs — pick longest that contains original seed
    ssake_contigs = Path(f"{ssake_prefix}.contigs")
    best = contig
    if ssake_contigs.exists():
        seqs = read_fasta(ssake_contigs)
        for name, seq in seqs.items():
            if len(seq) > len(best):
                # Check the SSAKE contig extends the seed (contains it or overlaps)
                if contig[:50] in seq or contig[-50:] in seq:
                    best = seq

    shutil.rmtree(ssake_dir, ignore_errors=True)
    return best


def _check_merge(contig_5p: str, contig_3p: str, min_overlap: int = 30) -> str | None:
    """Check if contig_5p and contig_3p overlap. Return merged seq or None."""
    if not contig_5p or not contig_3p:
        return None

    # Check if 3' end of 5p overlaps with 5' end of 3p
    max_check = min(len(contig_5p), len(contig_3p), 500)
    for ovl in range(max_check, min_overlap - 1, -1):
        if contig_5p[-ovl:] == contig_3p[:ovl]:
            merged = contig_5p + contig_3p[ovl:]
            log(f"    MERGE: {ovl}bp overlap detected! "
                f"→ {len(merged):,}bp merged contig")
            return merged
    return None


# ---------------------------------------------------------------------------
# Write pool FASTQ for Pilon
# ---------------------------------------------------------------------------

def _write_pool_fastq(extender: StrandAwareSeedExtender,
                      r1_out: Path, r2_out: Path) -> None:
    """Write extender's read pool as pseudo-paired FASTQ for Pilon mapping."""
    with open(r1_out, "w") as f1, open(r2_out, "w") as f2:
        for i in range(0, len(extender.seqs) - 1, 2):
            r1 = extender.seqs[i]
            r2_rc = extender.seqs[i + 1]
            r2_raw = revcomp(r2_rc)
            name = f"read_{i // 2}"
            f1.write(f"@{name}/1\n{r1}\n+\n{'I' * len(r1)}\n")
            f2.write(f"@{name}/2\n{r2_raw}\n+\n{'I' * len(r2_raw)}\n")


# ---------------------------------------------------------------------------
# Main assembly loop
# ---------------------------------------------------------------------------

def assemble_insert(
    site: InsertionSite,
    candidate_r1: Path,
    candidate_r2: Path,
    host_bam: Path,
    host_ref: Path,
    element_db: Path,
    workdir: Path,
    threads: int = 4,
    max_rounds: int = 8,
    ext_k: int = 15,
    recruit_k: int = 25,
    gap_size: int = 1000,
    s03_r1: Path | None = None,
    s03_r2: Path | None = None,
) -> tuple[Path | None, int, str]:
    """Main loop alternating k-mer extension and Pilon gap fill.

    Returns (fasta_path, n_rounds, status).
    Status: 'complete', 'partial', 'no_seeds', 'no_assembly'.
    """
    workdir.mkdir(parents=True, exist_ok=True)

    seed_5p = site.seed_5p
    seed_3p = site.seed_3p

    if not seed_5p and not seed_3p:
        log(f"  No seeds for {site.site_id}")
        return None, 0, "no_seeds"

    log(f"  Seeds: 5p={len(seed_5p)}bp, 3p={len(seed_3p)}bp")

    # Initialize extender
    extender = StrandAwareSeedExtender(
        k=ext_k, min_overlap=20, min_depth=2, min_ratio=0.7,
    )

    # Load construct-extracted reads (s03) as PRIMARY pool — these are the
    # reads that actually hit the construct, much more targeted than host BAM
    if s03_r1 and s03_r2 and s03_r1.exists() and s03_r2.exists():
        log(f"  Loading construct-extracted reads (s03)...")
        n_s03 = extender.load_paired_reads(s03_r1, s03_r2)
        log(f"  S03 reads: {n_s03:,} pairs")
    else:
        n_s03 = 0

    # Load junction-region candidate reads (supplementary)
    log(f"  Loading junction-region reads...")
    n_pairs = extender.load_paired_reads(candidate_r1, candidate_r2)
    log(f"  Total pool: {len(extender.seqs):,} seqs, "
        f"{len(extender.kmer_index):,} unique {ext_k}-mers")

    # Extract unmapped reads (cached)
    unmapped_r1, unmapped_r2 = extract_unmapped_paired(
        host_bam, workdir, threads=threads,
    )

    # Initial k-mer extension of both seeds
    contig_5p = seed_5p
    contig_3p = seed_3p

    if contig_5p:
        log(f"  Initial extension of 5p seed ({len(contig_5p)}bp)...")
        contig_5p = extender.extend(contig_5p, max_iterations=100)
        log(f"  5p: {len(seed_5p)} → {len(contig_5p):,}bp")

    if contig_3p:
        log(f"  Initial extension of 3p seed ({len(contig_3p)}bp)...")
        contig_3p = extender.extend(contig_3p, max_iterations=100)
        log(f"  3p: {len(seed_3p)} → {len(contig_3p):,}bp")

    # Check immediate merge
    merged = _check_merge(contig_5p, contig_3p)
    if merged:
        final_fa = workdir / f"{site.site_id}_insert.fasta"
        write_fasta(final_fa, f"{site.site_id}_assembled_insert", merged)
        log(f"  COMPLETE after initial extension: {len(merged):,}bp")
        return final_fa, 0, "complete"

    # Main iterative loop
    total_rounds = 0
    growth_history: list[int] = []
    status = "partial"

    for rnd in range(1, max_rounds + 1):
        prev_5p_len = len(contig_5p)
        prev_3p_len = len(contig_3p)
        log(f"\n  --- Round {rnd} ---")

        # Step A: K-mer recruit from unmapped pool
        recruit_contig = contig_5p + contig_3p
        r1_new, r2rc_new = recruit_by_kmer(
            recruit_contig, unmapped_r1, unmapped_r2, k=recruit_k,
        )
        if r1_new:
            n_added = extender.add_seqs(r1_new + r2rc_new)
            log(f"    Recruit: {len(r1_new)} pairs → {n_added} new seqs")
        else:
            log(f"    Recruit: 0 pairs")

        # Step B: K-mer extension (strand-aware)
        if contig_5p:
            contig_5p = extender.extend(contig_5p, max_iterations=100)
        if contig_3p:
            contig_3p = extender.extend(contig_3p, max_iterations=100)

        kmer_growth = (len(contig_5p) - prev_5p_len) + (len(contig_3p) - prev_3p_len)
        log(f"    K-mer ext: 5p={prev_5p_len:,}→{len(contig_5p):,}, "
            f"3p={prev_3p_len:,}→{len(contig_3p):,} (+{kmer_growth})")

        # Step B2: minimap2 soft-clip extension (handles repeats k-mer can't)
        mm2_dir = workdir / f"_mm2_r{rnd}"
        mm2_dir.mkdir(parents=True, exist_ok=True)
        pool_r1_mm2 = mm2_dir / "pool_R1.fq"
        pool_r2_mm2 = mm2_dir / "pool_R2.fq"
        _write_pool_fastq(extender, pool_r1_mm2, pool_r2_mm2)

        pre_mm2_5p = len(contig_5p)
        pre_mm2_3p = len(contig_3p)
        if contig_5p:
            contig_5p = _minimap2_extend(
                contig_5p, pool_r1_mm2, pool_r2_mm2, mm2_dir,
                threads=threads,
            )
        if contig_3p:
            contig_3p = _minimap2_extend(
                contig_3p, pool_r1_mm2, pool_r2_mm2, mm2_dir,
                threads=threads,
            )
        mm2_growth = (len(contig_5p) - pre_mm2_5p) + (len(contig_3p) - pre_mm2_3p)
        if mm2_growth > 0:
            log(f"    mm2 ext:  5p={pre_mm2_5p:,}→{len(contig_5p):,}, "
                f"3p={pre_mm2_3p:,}→{len(contig_3p):,} (+{mm2_growth})")
            # Recruit new reads matching the mm2-extended region
            r1_mm2, r2rc_mm2 = recruit_by_kmer(
                contig_5p + contig_3p, unmapped_r1, unmapped_r2, k=recruit_k,
            )
            if r1_mm2:
                n_mm2 = extender.add_seqs(r1_mm2 + r2rc_mm2)
                if n_mm2:
                    log(f"    mm2 recruit: {n_mm2} new seqs")
        shutil.rmtree(mm2_dir, ignore_errors=True)

        # Step C: Check merge
        merged = _check_merge(contig_5p, contig_3p)
        if merged:
            final_fa = workdir / f"{site.site_id}_insert.fasta"
            write_fasta(final_fa, f"{site.site_id}_assembled_insert", merged)
            log(f"  COMPLETE: contigs merged → {len(merged):,}bp (round {rnd})")
            return final_fa, rnd, "complete"

        # Step D: Pilon gap fill
        pre_pilon_5p = len(contig_5p)
        pre_pilon_3p = len(contig_3p)
        pilon_growth = 0
        if contig_5p and contig_3p:
            pilon_dir = workdir / f"_pilon_r{rnd}"
            pilon_dir.mkdir(parents=True, exist_ok=True)

            pool_r1 = pilon_dir / "pool_R1.fq"
            pool_r2 = pilon_dir / "pool_R2.fq"
            _write_pool_fastq(extender, pool_r1, pool_r2)

            contig_5p, contig_3p, gap_filled = pilon_fill(
                contig_5p, contig_3p,
                pool_r1, pool_r2,
                pilon_dir,
                gap_size=gap_size,
                threads=threads,
            )

            shutil.rmtree(pilon_dir, ignore_errors=True)
            pilon_growth = (len(contig_5p) - pre_pilon_5p) + (len(contig_3p) - pre_pilon_3p)

            if gap_filled:
                final_seq = contig_5p
                final_fa = workdir / f"{site.site_id}_insert.fasta"
                write_fasta(final_fa, f"{site.site_id}_assembled_insert", final_seq)
                log(f"  COMPLETE: Pilon filled gap → {len(final_seq):,}bp (round {rnd})")
                return final_fa, rnd, "complete"

        # Step D2: SSAKE overlap-layout-consensus extension
        ssake_dir = workdir / f"_ssake_r{rnd}"
        ssake_dir.mkdir(parents=True, exist_ok=True)
        ssake_r1 = ssake_dir / "pool_R1.fq"
        ssake_r2 = ssake_dir / "pool_R2.fq"
        _write_pool_fastq(extender, ssake_r1, ssake_r2)

        pre_ssake_5p = len(contig_5p)
        pre_ssake_3p = len(contig_3p)
        if contig_5p:
            contig_5p = _ssake_extend(contig_5p, ssake_r1, ssake_r2, ssake_dir)
        if contig_3p:
            contig_3p = _ssake_extend(contig_3p, ssake_r1, ssake_r2, ssake_dir)
        ssake_growth = (len(contig_5p) - pre_ssake_5p) + (len(contig_3p) - pre_ssake_3p)
        if ssake_growth > 0:
            log(f"    SSAKE:   5p={pre_ssake_5p:,}→{len(contig_5p):,}, "
                f"3p={pre_ssake_3p:,}→{len(contig_3p):,} (+{ssake_growth})")
            # Check merge after SSAKE
            merged = _check_merge(contig_5p, contig_3p)
            if merged:
                final_fa = workdir / f"{site.site_id}_insert.fasta"
                write_fasta(final_fa, f"{site.site_id}_assembled_insert", merged)
                log(f"  COMPLETE: SSAKE+merge → {len(merged):,}bp (round {rnd})")
                shutil.rmtree(ssake_dir, ignore_errors=True)
                return final_fa, rnd, "complete"
        shutil.rmtree(ssake_dir, ignore_errors=True)

        # Step E: Check termination — only converge when ALL 4 steps show zero growth
        total_rounds = rnd
        log(f"    Round {rnd} growth: kmer={kmer_growth}, mm2={mm2_growth}, "
            f"pilon={pilon_growth}, ssake={ssake_growth}")

        if kmer_growth == 0 and mm2_growth == 0 and pilon_growth == 0 and ssake_growth == 0:
            log(f"    All 4 assemblers show zero growth → converged")
            status = "converged"
            break

        # Check if contig ends reached host genome
        reached_5p, reached_3p = check_host_termination(
            contig_5p, contig_3p, host_ref, workdir,
        )
        if reached_5p:
            log(f"    5' contig reached host genome!")
        if reached_3p:
            log(f"    3' contig reached host genome!")
        if reached_5p and reached_3p:
            log(f"    >>> Both ends at host genome — insert fully traversed!")
            status = "complete"
            break

        # Cycle detection — per-step growth pattern over 3 rounds
        total_growth = kmer_growth + mm2_growth + pilon_growth + ssake_growth
        growth_history.append(total_growth)
        if len(growth_history) >= 3:
            last3 = growth_history[-3:]
            if last3[0] == last3[1] == last3[2] and last3[0] > 0:
                log(f"    Same total growth ({last3[0]}) for 3 rounds → cycling")
                status = "converged"
                break
        if len(growth_history) >= 6:
            recent = growth_history[-3:]
            earlier = growth_history[-6:-3]
            if recent == earlier:
                log(f"    Cyclic growth pattern detected → stopping")
                status = "converged"
                break

    # Write final result
    if contig_5p and contig_3p:
        # Output both contigs separately if not merged
        final_seq = contig_5p + ("N" * 100) + contig_3p
        final_fa = workdir / f"{site.site_id}_insert.fasta"
        write_fasta(final_fa, f"{site.site_id}_assembled_insert", final_seq)
        log(f"  Final: 5p={len(contig_5p):,}bp + 3p={len(contig_3p):,}bp "
            f"= {len(final_seq):,}bp after {total_rounds} rounds → {status}")
    elif contig_5p:
        final_fa = workdir / f"{site.site_id}_insert.fasta"
        write_fasta(final_fa, f"{site.site_id}_assembled_insert", contig_5p)
        log(f"  Final: {len(contig_5p):,}bp (5p only) after {total_rounds} rounds")
    elif contig_3p:
        final_fa = workdir / f"{site.site_id}_insert.fasta"
        write_fasta(final_fa, f"{site.site_id}_assembled_insert", contig_3p)
        log(f"  Final: {len(contig_3p):,}bp (3p only) after {total_rounds} rounds")
    else:
        return None, total_rounds, "no_assembly"

    # Phase 3b: Foreign read refinement (minimap2 chaining + Pilon)
    # After k-mer assembly converges, use foreign reads (s03 unmapped to host)
    # with iterative minimap2+Pilon to extend construct regions.
    # This resolves palindromic/inverted repeat regions that k-mer extension
    # cannot handle (e.g., head-to-head T-DNA structures).
    if s03_r1 and s03_r2 and s03_r1.exists() and s03_r2.exists():
        log(f"\n  === Foreign read refinement (minimap2 + Pilon) ===")
        final_fa = refine_with_foreign_reads(
            insert_fasta=final_fa,
            s03_r1=s03_r1,
            s03_r2=s03_r2,
            host_ref=host_ref,
            workdir=workdir,
            threads=threads,
            max_rounds=10,
        )

    return final_fa, total_rounds, status


# ---------------------------------------------------------------------------
# Phase 4: Annotation
# ---------------------------------------------------------------------------

def _parse_blast6(path: Path, min_len: int = 30) -> list[dict]:
    """Parse BLAST outfmt-6 with 10+ columns into list of dicts."""
    hits: list[dict] = []
    if not path.exists():
        return hits
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip().split("\t")
            if len(cols) < 10:
                continue
            aln_len = int(cols[3])
            if aln_len < min_len:
                continue
            hits.append({
                "query": cols[0], "subject": cols[1],
                "identity": float(cols[2]), "length": aln_len,
                "q_start": int(cols[4]), "q_end": int(cols[5]),
                "s_start": int(cols[6]), "s_end": int(cols[7]),
                "evalue": float(cols[8]), "bitscore": float(cols[9]),
            })
    return hits


def _run_local_blast(
    insert_fasta: Path, element_db: Path, output_dir: Path,
    tag: str | None = None,
) -> list[dict]:
    """Local BLAST vs element_db. Fast, covers known GMO elements.

    ``tag`` disambiguates intermediate filenames when the function is
    called multiple times in the same ``output_dir`` (e.g. once for the
    primary element_db and again for each extra DB). Defaults to the
    stem of ``element_db`` so callers can omit it.
    """
    suffix = tag if tag is not None else element_db.stem
    db_prefix = output_dir / f"_element_blastdb_{suffix}"
    subprocess.run(
        ["makeblastdb", "-in", str(element_db), "-dbtype", "nucl",
         "-out", str(db_prefix)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
    )
    blast_out = output_dir / f"_local_blast_{suffix}.tsv"
    subprocess.run(
        ["blastn", "-query", str(insert_fasta), "-db", str(db_prefix),
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
         "-evalue", "1e-5", "-max_target_seqs", "50",
         "-out", str(blast_out)],
        stderr=subprocess.DEVNULL, check=True,
    )
    hits = _parse_blast6(blast_out)
    # Cleanup
    blast_out.unlink(missing_ok=True)
    for ext in [".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto", ".njs"]:
        Path(f"{db_prefix}{ext}").unlink(missing_ok=True)
    log(f"  Local BLAST ({suffix}): {len(hits)} hits")
    return hits


def _run_remote_blast(
    insert_fasta: Path, output_dir: Path,
    timeout: int = 600, max_retries: int = 2,
) -> list[dict]:
    """Remote BLAST vs NCBI nt. Slower, but annotates unknown regions."""
    blast_out = output_dir / "_remote_blast.tsv"
    for attempt in range(1, max_retries + 1):
        log(f"  Remote BLAST vs NCBI nt (attempt {attempt}/{max_retries})...")
        try:
            proc = subprocess.run(
                ["blastn", "-query", str(insert_fasta), "-db", "nt",
                 "-remote",
                 "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle",
                 "-evalue", "1e-10", "-max_target_seqs", "10",
                 "-out", str(blast_out)],
                stderr=subprocess.PIPE, timeout=timeout,
            )
            if proc.returncode == 0 and blast_out.exists():
                break
            log(f"    Remote BLAST returned code {proc.returncode}")
        except subprocess.TimeoutExpired:
            log(f"    Remote BLAST timed out ({timeout}s)")
        except FileNotFoundError:
            log("    blastn not found — skipping remote BLAST")
            return []
    else:
        log("  Remote BLAST failed after all retries — skipping")
        return []

    # Parse — outfmt has 11 columns (extra stitle), but _parse_blast6 needs 10+
    hits: list[dict] = []
    if blast_out.exists():
        with open(blast_out) as fh:
            for line in fh:
                cols = line.rstrip().split("\t")
                if len(cols) < 10:
                    continue
                aln_len = int(cols[3])
                if aln_len < 30:
                    continue
                # Use stitle (col 10) as a readable subject description
                stitle = cols[10] if len(cols) > 10 else cols[1]
                hits.append({
                    "query": cols[0],
                    "subject": f"{cols[1]}|{stitle}",
                    "identity": float(cols[2]), "length": aln_len,
                    "q_start": int(cols[4]), "q_end": int(cols[5]),
                    "s_start": int(cols[6]), "s_end": int(cols[7]),
                    "evalue": float(cols[8]), "bitscore": float(cols[9]),
                })
        blast_out.unlink(missing_ok=True)

    log(f"  Remote BLAST: {len(hits)} hits from NCBI nt")
    return hits


def _merge_annotations(
    local_hits: list[dict], remote_hits: list[dict],
) -> list[dict]:
    """Merge local + remote hits, keeping best bitscore per query region.

    For each position along the insert, prefer the hit with highest bitscore.
    Local hits from element_db are preferred at equal score (more specific names).
    """
    # Tag source
    for h in local_hits:
        h["source"] = "element_db"
    for h in remote_hits:
        h["source"] = "ncbi_nt"

    all_hits = local_hits + remote_hits
    if not all_hits:
        return []

    # Group by query sequence
    from collections import defaultdict as _dd
    by_query: dict[str, list[dict]] = _dd(list)
    for h in all_hits:
        by_query[h["query"]].append(h)

    merged: list[dict] = []
    for qname, hits in by_query.items():
        # Sort by bitscore descending, prefer element_db on tie
        hits.sort(key=lambda h: (-h["bitscore"],
                                  0 if h["source"] == "element_db" else 1))

        # Greedy interval selection: pick best non-overlapping hits
        # (allow 80% reciprocal overlap to keep significant alternatives)
        selected: list[dict] = []
        covered: list[tuple[int, int]] = []

        for h in hits:
            qs, qe = min(h["q_start"], h["q_end"]), max(h["q_start"], h["q_end"])
            h_len = qe - qs + 1

            # Check overlap with already-selected hits
            dominated = False
            for cs, ce in covered:
                overlap = max(0, min(qe, ce) - max(qs, cs) + 1)
                if overlap > 0.80 * h_len:
                    dominated = True
                    break
            if not dominated:
                selected.append(h)
                covered.append((qs, qe))

        merged.extend(selected)

    merged.sort(key=lambda h: (h["query"], h["q_start"]))
    return merged


def annotate_insert(
    insert_fasta: Path,
    element_db: Path,
    output_dir: Path,
    sample_name: str,
    no_remote_blast: bool = False,
    extra_dbs: list[Path] | None = None,
) -> tuple[Path, Path]:
    """Annotate insert with local element_db BLAST + remote NCBI nt BLAST.

    Runs both in sequence (local is fast, remote may take 1-5 min),
    then merges by best bitscore per region. Output format is unchanged
    for downstream report generation.

    ``extra_dbs`` mirrors the Phase 1.5 plumbing (common_payload.fa and/or
    per-sample s04b SPAdes contigs): each extra DB is BLASTed separately
    and its hits are concatenated into the local-hit stream before the
    best-bitscore merge, so sample-specific payloads (e.g. bar, AtYUCCA6,
    full T-DNA backbone contigs) surface in ``element_annotation.tsv``
    alongside the shared EUginius catalogue.
    """
    annotation_tsv = output_dir / "element_annotation.tsv"
    border_tsv = output_dir / "border_hits.tsv"

    # ---- Local BLAST (element_db + any extra DBs) ----
    local_hits = _run_local_blast(
        insert_fasta, element_db, output_dir, tag="primary",
    )
    for i, edb in enumerate(extra_dbs or []):
        if edb is None or not edb.exists() or edb.stat().st_size == 0:
            continue
        extra_hits = _run_local_blast(
            insert_fasta, edb, output_dir, tag=f"extra{i}_{edb.stem}",
        )
        local_hits.extend(extra_hits)

    # ---- Remote BLAST (NCBI nt) ----
    if no_remote_blast:
        log("  Remote BLAST skipped (--no-remote-blast)")
        remote_hits = []
    else:
        remote_hits = _run_remote_blast(insert_fasta, output_dir)

    # ---- Merge ----
    merged = _merge_annotations(local_hits, remote_hits)
    log(f"  Merged annotation: {len(merged)} regions "
        f"({sum(1 for h in merged if h['source'] == 'element_db')} local, "
        f"{sum(1 for h in merged if h['source'] == 'ncbi_nt')} remote)")

    # ---- Write annotation TSV ----
    with open(annotation_tsv, "w") as fout:
        fout.write("query\telement\tidentity\tlength\tq_start\tq_end\t"
                   "s_start\ts_end\tevalue\tsource\n")
        for h in merged:
            fout.write(f"{h['query']}\t{h['subject']}\t{h['identity']}\t"
                       f"{h['length']}\t{h['q_start']}\t{h['q_end']}\t"
                       f"{h['s_start']}\t{h['s_end']}\t{h['evalue']}\t"
                       f"{h['source']}\n")

    if merged:
        for h in merged[:15]:
            src = "L" if h["source"] == "element_db" else "R"
            log(f"    {h['q_start']:>6}-{h['q_end']:<6} [{src}] "
                f"{h['subject'][:60]} ({h['identity']:.1f}%, {h['length']}bp)")
    else:
        log("  No BLAST hits found")

    # ---- T-DNA border motif search ----
    border_fa = output_dir / "_borders.fa"
    with open(border_fa, "w") as fh:
        fh.write(">RB_consensus\nTGGCAGGATATATTGTGGTGTAAAC\n")
        fh.write(">LB_consensus\nTGGCAGGATATATTGTGGTGTAAAC\n")

    subprocess.run(
        ["blastn", "-query", str(border_fa), "-subject", str(insert_fasta),
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send",
         "-evalue", "1", "-word_size", "7",
         "-out", str(border_tsv)],
        stderr=subprocess.DEVNULL,
    )

    if border_tsv.exists() and border_tsv.stat().st_size > 0:
        log("  T-DNA border motifs found")
    else:
        log("  No border motifs found (may need manual inspection)")

    border_fa.unlink(missing_ok=True)
    return annotation_tsv, border_tsv


# ---------------------------------------------------------------------------
# Post-assembly host-fraction filter
# ---------------------------------------------------------------------------

def _find_construct_flanking_regions(
    construct_ref: Path,
    host_ref: Path,
    workdir: Path,
    threads: int = 4,
) -> list[tuple[str, int, int]]:
    """BLAST construct reference vs host to find host-flanking regions.

    Some construct references include host genomic flanking DNA at their ends
    (e.g., rice_G281 has 201bp of Chr11 at positions 8758-8958).  Soft-clip
    sites at these host coordinates are false detections.

    Returns list of (host_chr, start, end) tuples for flanking regions.
    """
    if not construct_ref.exists():
        return []

    blast_out = workdir / "_construct_vs_host_flanking.tsv"
    result = subprocess.run(
        ["blastn", "-task", "megablast",
         "-query", str(construct_ref), "-db", str(host_ref),
         "-outfmt", "6 qseqid qlen qstart qend sseqid sstart send pident length",
         "-evalue", "1e-20", "-max_target_seqs", "5",
         "-num_threads", str(threads),
         "-out", str(blast_out)],
        stderr=subprocess.DEVNULL,
    )
    if result.returncode != 0 or not blast_out.exists():
        return []

    flanking: list[tuple[str, int, int]] = []
    with open(blast_out) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            qlen = int(cols[1])
            q_start = int(cols[2])
            q_end = int(cols[3])
            s_chr = cols[4]
            s_start = int(cols[5])
            s_end = int(cols[6])
            pident = float(cols[7])
            aln_len = int(cols[8])

            if pident < CONSTRUCT_FLANK_PIDENT or aln_len < CONSTRUCT_FLANK_MIN_LEN:
                continue

            # Only count hits near the ends of a construct entry (flanking),
            # not internal matches (e.g., host-derived promoters handled
            # separately by _filter_host_endogenous).
            near_start = q_start <= 100
            near_end = q_end >= qlen - 100
            if not (near_start or near_end):
                continue

            lo = min(s_start, s_end)
            hi = max(s_start, s_end)
            flanking.append((s_chr, lo, hi))

    if flanking:
        log(f"  Construct-flanking regions found: {len(flanking)}")
        for chrom, lo, hi in flanking:
            log(f"    {chrom}:{lo:,}-{hi:,} ({hi - lo + 1}bp)")

    return flanking


def _site_overlaps_flanking(
    site_chr: str,
    site_pos: int,
    flanking_regions: list[tuple[str, int, int]],
    slop: int = CONSTRUCT_FLANK_SLOP,
) -> tuple[bool, str]:
    """Check if a detection site overlaps a construct-flanking host region."""
    for chrom, lo, hi in flanking_regions:
        if chrom == site_chr and lo - slop <= site_pos <= hi + slop:
            return True, f"{chrom}:{lo:,}-{hi:,}"
    return False, ""


def _check_chimeric_assembly(
    insert_fasta: Path,
    host_ref: Path,
    site_chr: str,
    workdir: Path,
    threads: int = 4,
) -> tuple[bool, list[tuple[str, int]]]:
    """Check if assembled insert contains DNA from multiple host chromosomes.

    Returns (is_chimeric, off_target_hits) where off_target_hits is a list of
    (chromosome, aligned_bp) for chromosomes other than site_chr.
    """
    blast_out = workdir / f"_{insert_fasta.stem}_vs_host_chrom.tsv"
    # Reuse existing BLAST output if available (from _blast_insert_vs_host)
    if not blast_out.exists():
        result = subprocess.run(
            ["blastn", "-task", "megablast",
             "-query", str(insert_fasta), "-db", str(host_ref),
             "-outfmt", "6 qseqid qstart qend sseqid pident length",
             "-evalue", "1e-10", "-max_target_seqs", "10",
             "-num_threads", str(threads),
             "-out", str(blast_out)],
            stderr=subprocess.DEVNULL,
        )
        if result.returncode != 0 or not blast_out.exists():
            return False, []

    # Accumulate aligned bp per chromosome (strict identity to avoid
    # counting element-level homologies as chimeric evidence)
    chr_bp: dict[str, int] = defaultdict(int)
    with open(blast_out) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue
            s_chr = cols[3]
            pident = float(cols[4])
            aln_len = int(cols[5])
            if pident >= CHIMERIC_MIN_PIDENT:
                chr_bp[s_chr] += aln_len

    # Find off-target chromosomes with significant coverage
    off_target: list[tuple[str, int]] = []
    for chrom, bp in sorted(chr_bp.items(), key=lambda x: -x[1]):
        if chrom != site_chr and bp >= CHIMERIC_MIN_OFFTARGET_BP:
            off_target.append((chrom, bp))

    is_chimeric = len(off_target) >= 2
    return is_chimeric, off_target


def _check_construct_host_coverage(
    insert_fasta: Path,
    element_db: Path,
    host_fraction: float,
    host_bp: int,
    insert_len: int,
    n_count: int,
    workdir: Path,
    threads: int = 4,
) -> tuple[bool, float, float, float]:
    """Check if insert is fully explained by construct + host coverage.

    Real T-DNA inserts have low construct coverage (~10%, only border regions).
    False positives from construct-element homology have high construct
    coverage (~50%, the insert IS the construct fragment).

    Returns (is_fp, construct_frac, host_frac, combined_frac).
    """
    eff_len = insert_len - n_count
    if eff_len <= 0:
        return False, 0.0, 0.0, 0.0

    # BLAST insert vs construct/element_db
    blast_out = workdir / f"_{insert_fasta.stem}_vs_construct.tsv"
    result = subprocess.run(
        ["blastn", "-task", "megablast",
         "-query", str(insert_fasta), "-subject", str(element_db),
         "-outfmt", "6 qstart qend pident",
         "-evalue", "1e-5"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    if result.returncode != 0:
        return False, 0.0, 0.0, 0.0

    # Merge overlapping intervals at >= threshold identity
    intervals: list[tuple[int, int]] = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        cols = line.split("\t")
        qs, qe, pid = int(cols[0]), int(cols[1]), float(cols[2])
        if pid >= CONSTRUCT_HOST_MIN_PIDENT:
            intervals.append((min(qs, qe), max(qs, qe)))

    if not intervals:
        return False, 0.0, host_fraction, host_fraction

    intervals.sort()
    merged: list[tuple[int, int]] = []
    for s, e in intervals:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))

    construct_bp = sum(e - s + 1 for s, e in merged)
    construct_frac = construct_bp / eff_len
    combined_frac = min(1.0, (construct_bp + host_bp) / eff_len)

    is_fp = (construct_frac >= CONSTRUCT_MIN_FRACTION
             and combined_frac >= CONSTRUCT_HOST_MIN_COMBINED)

    return is_fp, construct_frac, host_fraction, combined_frac


def _blast_insert_vs_host(
    insert_fasta: Path,
    host_ref: Path,
    workdir: Path,
    threads: int = 4,
) -> tuple[float, int, int, int]:
    """BLAST assembled insert against host genome to measure host-fraction.

    Returns (host_fraction, host_covered_bp, insert_length, largest_foreign_gap).
    host_fraction = fraction of insert positions covered by host alignments.
    largest_foreign_gap = longest contiguous stretch NOT covered by host.
    """
    seqs = read_fasta(insert_fasta)
    if not seqs:
        return 0.0, 0, 0, 0
    insert_seq = list(seqs.values())[0]
    insert_len = len(insert_seq)
    if insert_len == 0:
        return 0.0, 0, 0, 0

    # Use _chrom variant filename so _check_chimeric_assembly can reuse it
    blast_out = workdir / f"_{insert_fasta.stem}_vs_host_chrom.tsv"
    if not blast_out.exists():
        result = subprocess.run(
            ["blastn", "-task", "megablast",
             "-query", str(insert_fasta), "-db", str(host_ref),
             "-outfmt", "6 qseqid qstart qend sseqid pident length",
             "-evalue", "1e-10", "-max_target_seqs", "10",
             "-num_threads", str(threads),
             "-out", str(blast_out)],
            stderr=subprocess.DEVNULL,
        )
    if not blast_out.exists():
        return 0.0, 0, insert_len, insert_len

    # Merge overlapping host-aligned intervals to compute non-redundant coverage
    intervals: list[tuple[int, int]] = []
    with open(blast_out) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue
            pident = float(cols[4])
            if pident < INSERT_HOST_MIN_PIDENT:
                continue
            q_start = int(cols[1])
            q_end = int(cols[2])
            lo, hi = min(q_start, q_end), max(q_start, q_end)
            intervals.append((lo, hi))

    if not intervals:
        return 0.0, 0, insert_len, insert_len

    # Merge overlapping intervals
    intervals.sort()
    merged: list[tuple[int, int]] = [intervals[0]]
    for lo, hi in intervals[1:]:
        if lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))

    host_bp = sum(hi - lo + 1 for lo, hi in merged)

    # Compute largest gap between host-aligned regions (= foreign/T-DNA region)
    # Include gaps at start and end of insert
    gaps: list[int] = []
    gaps.append(merged[0][0] - 1)                    # gap before first host hit
    for i in range(1, len(merged)):
        gaps.append(merged[i][0] - merged[i - 1][1] - 1)
    gaps.append(insert_len - merged[-1][1])           # gap after last host hit
    largest_gap = max(gaps) if gaps else 0

    # Exclude N-runs from denominator (Ns are gap-fill placeholders, not real sequence)
    n_count = insert_seq.upper().count("N")
    effective_len = insert_len - n_count
    if effective_len <= 0:
        return 0.0, host_bp, insert_len, largest_gap

    return host_bp / effective_len, host_bp, insert_len, largest_gap


def generate_report(
    insert_fasta: Path,
    annotation_tsv: Path,
    border_tsv: Path,
    site: InsertionSite,
    # NOTE: verdict logic (CANDIDATE/FALSE_POSITIVE/UNKNOWN) is embedded here.
    # If a verdict-only mode is needed later, extract into compute_verdict().
    n_rounds: int,
    status: str,
    output_dir: Path,
    host_endogenous_ids: set[str] | None = None,
    host_ref: Path | None = None,
    element_db: Path | None = None,
    threads: int = 4,
    construct_flanking: list[tuple[str, int, int]] | None = None,
) -> tuple[Path, str]:
    """Generate human-readable linear map report."""
    report_path = output_dir / f"{site.site_id}_report.txt"

    # Read insert sequence
    seqs = read_fasta(insert_fasta)
    if not seqs:
        with open(report_path, "w") as fh:
            fh.write("No insert assembled.\n")
        return report_path, "NO_ASSEMBLY"

    insert_name = list(seqs.keys())[0]
    insert_seq = seqs[insert_name]
    insert_len = len(insert_seq)
    n_count = insert_seq.upper().count("N")

    # Read annotation hits (with orientation from s_start/s_end)
    # elements: (q_start, q_end, element_name, s_strand, source)
    elements: list[tuple[int, int, str, str, str]] = []
    orientation_by_elem: dict[str, set[str]] = defaultdict(set)
    if annotation_tsv.exists():
        with open(annotation_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    q_start = int(row["q_start"])
                    q_end = int(row["q_end"])
                    s_start = int(row["s_start"])
                    s_end = int(row["s_end"])
                    elem = row["element"]
                    source = row.get("source", "element_db")
                    # Filter to current insert only
                    if row["query"] != insert_name:
                        continue
                    strand = "+" if s_start < s_end else "-"
                    elements.append((q_start, q_end, elem, strand, source))
                    orientation_by_elem[elem].add(strand)
                except (ValueError, KeyError):
                    continue
    elements.sort(key=lambda x: x[0])

    # Deduplicate overlapping elements (keep longest per region)
    deduped: list[tuple[int, int, str, str, str]] = []
    for start, end, elem, strand, source in elements:
        overlap = False
        for i, (ds, de, dn, _ds, _src) in enumerate(deduped):
            if start < de and end > ds:
                if (end - start) > (de - ds):
                    deduped[i] = (start, end, elem, strand, source)
                overlap = True
                break
        if not overlap:
            deduped.append((start, end, elem, strand, source))
    elements = sorted(deduped, key=lambda x: x[0])

    # Count element occurrences for multi-construct detection
    elem_counts: Counter = Counter()
    for _, _, elem, _, _ in elements:
        elem_counts[elem] += 1

    # ---- Phase 4 post-filter: host-endogenous verdict ----
    # If every annotated element was excluded at the DB-level host BLAST
    # (Tier 1 or Tier 2), the whole assembly is host-endogenous noise.
    # A single genuinely foreign element is enough to keep the site.
    host_endo = host_endogenous_ids or set()
    unique_elems = set(elem_counts)
    foreign_elems: set[str] = set()
    endo_elems: set[str] = set()
    for elem in unique_elems:
        if elem in host_endo:
            endo_elems.add(elem)
        else:
            foreign_elems.add(elem)
    # Need at least one annotation AND at least one foreign element to keep it.
    if unique_elems and not foreign_elems:
        verdict = "FALSE_POSITIVE"
        verdict_reason = (
            "all annotated elements are host-endogenous "
            "(matched host genome at Tier 1/2 BLAST threshold)"
        )
    elif not unique_elems:
        verdict = "UNKNOWN"
        verdict_reason = "no element annotations"
    else:
        verdict = "CANDIDATE"
        verdict_reason = ""

    # ---- Post-assembly false-positive filters (applied to CANDIDATE only) ----
    # Filter A: Host-fraction + small gap → chimeric assembly artifact
    host_fraction = 0.0
    host_bp = 0
    largest_gap = 0
    # Filter B: Construct-flanking overlap
    flanking_hit = ""
    # Filter C: Multi-locus chimeric assembly
    is_chimeric = False
    off_target_chrs: list[tuple[str, int]] = []
    # Filter D: Construct-host coverage
    construct_frac = 0.0
    combined_frac = 0.0

    if verdict == "CANDIDATE" and host_ref is not None:
        # Filter A: host-fraction with non-host gap check
        host_fraction, host_bp, _, largest_gap = _blast_insert_vs_host(
            insert_fasta, host_ref, output_dir, threads=threads,
        )
        if (host_fraction >= INSERT_HOST_FRACTION
                and largest_gap < INSERT_MIN_FOREIGN_GAP):
            verdict = "FALSE_POSITIVE"
            verdict_reason = (
                f"assembled insert is {host_fraction:.0%} host genome "
                f"({host_bp:,}/{insert_len - n_count:,}bp) with only "
                f"{largest_gap}bp non-host gap "
                f"(need ≥{INSERT_MIN_FOREIGN_GAP}bp for real T-DNA)"
            )

    if verdict == "CANDIDATE" and construct_flanking:
        # Filter B: site coordinates overlap with construct→host flanking
        site_pos = site.pos_5p
        overlaps, flanking_hit = _site_overlaps_flanking(
            site.host_chr, site_pos, construct_flanking,
        )
        if overlaps:
            verdict = "FALSE_POSITIVE"
            verdict_reason = (
                f"site {site.host_chr}:{site_pos:,} overlaps construct-flanking "
                f"region {flanking_hit} — host DNA in construct reference"
            )

    if verdict == "CANDIDATE" and host_ref is not None:
        # Filter C: assembled insert spans ≥2 off-target chromosomes
        is_chimeric, off_target_chrs = _check_chimeric_assembly(
            insert_fasta, host_ref, site.host_chr, output_dir, threads=threads,
        )
        if is_chimeric:
            off_str = ", ".join(f"{c}:{bp}bp" for c, bp in off_target_chrs)
            verdict = "FALSE_POSITIVE"
            verdict_reason = (
                f"chimeric assembly — host-aligned portions span "
                f"{len(off_target_chrs)} off-target chromosomes ({off_str})"
            )

    if verdict == "CANDIDATE" and host_ref is not None and element_db is not None:
        # Filter D: construct + host coverage fully explains insert
        is_construct_fp, construct_frac, _, combined_frac = \
            _check_construct_host_coverage(
                insert_fasta, element_db,
                host_fraction, host_bp, insert_len, n_count,
                output_dir, threads=threads,
            )
        if is_construct_fp:
            verdict = "FALSE_POSITIVE"
            verdict_reason = (
                f"insert fully explained by construct ({construct_frac:.0%}) "
                f"+ host ({host_fraction:.0%}) = {combined_frac:.0%} combined "
                f"coverage — host genomic DNA with construct-element homology"
            )

    # ---- UNKNOWN reclassification: host-only insert ----
    # If no element annotations but insert is mostly host DNA with
    # negligible construct match, reclassify as FALSE_POSITIVE.
    if verdict == "UNKNOWN" and host_ref is not None:
        log(f"  UNKNOWN reclassification: BLASTing {insert_fasta.name} vs host")
        host_fraction, host_bp, _, largest_gap = _blast_insert_vs_host(
            insert_fasta, host_ref, output_dir, threads=threads,
        )
        construct_frac_unk = 0.0
        if element_db is not None:
            _, construct_frac_unk, _, _ = _check_construct_host_coverage(
                insert_fasta, element_db,
                host_fraction, host_bp, insert_len, n_count,
                output_dir, threads=threads,
            )
        if (host_fraction >= UNKNOWN_HOST_MIN_FRACTION
                and construct_frac_unk <= UNKNOWN_MAX_CONSTRUCT_FRAC):
            verdict = "FALSE_POSITIVE"
            construct_frac = construct_frac_unk
            verdict_reason = (
                f"no element annotations; insert is {host_fraction:.0%} host "
                f"genome ({host_bp:,}bp) with {construct_frac_unk:.0%} construct "
                f"— host genomic DNA"
            )

    # Build linear map with orientation arrows
    linear_parts = []
    for _, _, elem, strand, _ in elements:
        short_name = elem.split("|")[0] if "|" in elem else (
            elem.split("_")[0] if "_" in elem else elem)
        arrow = "→" if strand == "+" else "←"
        linear_parts.append(f"[{short_name}{arrow}]")

    # Detect head-to-head / tandem arrangement
    # Head-to-head: same element found in BOTH orientations
    structure = "single-copy"
    bidirectional_elems = {e for e, strands in orientation_by_elem.items()
                          if "+" in strands and "-" in strands
                          and len(e) > 30}  # skip short amplicons
    if bidirectional_elems:
        structure = "head-to-head 2-copy T-DNA"
    elif len(elements) > 3:
        if elem_counts.most_common(1)[0][1] >= 2:
            structure = f"multi-copy (≥{elem_counts.most_common(1)[0][1]} copies)"

    # Read border hits
    n_borders = 0
    if border_tsv.exists():
        with open(border_tsv) as fh:
            n_borders = sum(1 for line in fh if line.strip())

    # Determine deletion size
    deletion_size = abs(site.pos_3p - site.pos_5p) if site.pos_3p > 0 else 0

    with open(report_path, "w") as fh:
        fh.write("=" * 70 + "\n")
        fh.write("RedGene Insert Assembly & Annotation Report\n")
        fh.write("=" * 70 + "\n")
        pos_str = f"{site.host_chr}:{site.pos_5p:,}"
        if site.pos_3p > 0 and site.pos_3p != site.pos_5p:
            pos_str += f"-{site.pos_3p:,}"
        fh.write(f"Insertion site: {pos_str}")
        if deletion_size > 0:
            fh.write(f" ({deletion_size}bp deletion)")
        fh.write("\n")
        fh.write(f"Insert length: {insert_len:,} bp")
        if n_count > 0:
            fh.write(f" ({n_count} unresolved N's)")
        fh.write("\n")
        fh.write(f"Assembly status: {status.upper()} (round {n_rounds})\n")
        fh.write(f"Structure: {structure}\n")
        fh.write(f"Verdict: {verdict}")
        if verdict_reason:
            fh.write(f" — {verdict_reason}")
        fh.write("\n")
        if n_borders > 0:
            fh.write(f"T-DNA borders found: {n_borders}\n")
        fh.write("\n")

        if linear_parts:
            fh.write("--- Linear Map ---\n")
            fh.write("5' host --")
            line = ""
            for part in linear_parts:
                if len(line) + len(part) > 60:
                    fh.write(line + "\n         ")
                    line = ""
                line += part + "--"
            fh.write(line + " 3' host\n")
        else:
            fh.write("--- No element annotations found ---\n")

        fh.write("=" * 70 + "\n")

        # Detailed element list
        if elements:
            fh.write("\nDetailed element positions:\n")
            fh.write(f"{'Start':>8}  {'End':>8}  {'Dir':>3}  {'Src':>5}  {'Element'}\n")
            fh.write("-" * 70 + "\n")
            for start, end, elem, strand, source in elements:
                src_tag = "local" if source == "element_db" else "NCBI"
                fh.write(f"{start:>8}  {end:>8}  {strand:>3}  {src_tag:>5}  {elem}\n")

        # Multi-copy detection
        multi = {e: c for e, c in elem_counts.items() if c >= 2}
        if multi:
            fh.write(f"\nMulti-copy elements detected:\n")
            for elem, count in multi.items():
                fh.write(f"  {elem}: {count} copies\n")

        # Head-to-head evidence
        if bidirectional_elems:
            fh.write(f"\nHead-to-head evidence (same element in both orientations):\n")
            for elem in sorted(bidirectional_elems):
                fh.write(f"  {elem}\n")

        # Foreign vs host-endogenous element breakdown
        if host_endo and unique_elems:
            fh.write(f"\nForeign elements (not in host genome): {len(foreign_elems)}\n")
            for elem in sorted(foreign_elems):
                fh.write(f"  + {elem}\n")
            if endo_elems:
                fh.write(f"Host-endogenous elements "
                         f"(excluded at DB-level BLAST): {len(endo_elems)}\n")
                for elem in sorted(endo_elems):
                    fh.write(f"  - {elem}\n")

        # Post-assembly filter diagnostics
        if host_fraction > 0:
            fh.write(f"\nHost-fraction analysis: {host_fraction:.1%} host "
                     f"({host_bp:,}/{insert_len - n_count:,}bp at "
                     f"≥{INSERT_HOST_MIN_PIDENT}% identity), "
                     f"largest non-host gap: {largest_gap:,}bp\n")
        if flanking_hit:
            fh.write(f"Construct-flanking overlap: site overlaps {flanking_hit}\n")
        if off_target_chrs:
            off_str = ", ".join(f"{c}:{bp}bp" for c, bp in off_target_chrs)
            tag = " → chimeric" if is_chimeric else ""
            fh.write(f"Off-target chromosomes in assembly: {off_str}{tag}\n")
        if construct_frac > 0:
            fh.write(f"Construct coverage: {construct_frac:.1%}, "
                     f"combined (construct+host): {combined_frac:.1%}\n")

    log(f"  Report written: {report_path} [{verdict}]")
    return report_path, verdict


# ---------------------------------------------------------------------------
# Stats output
# ---------------------------------------------------------------------------

def write_stats(
    stats_path: Path,
    sample_name: str,
    num_sites: int,
    candidate_reads: int,
    results: list[dict],
) -> None:
    with open(stats_path, "w") as fh:
        fh.write(f"sample\t{sample_name}\n")
        fh.write(f"insertion_sites\t{num_sites}\n")
        fh.write(f"candidate_read_pairs\t{candidate_reads}\n")
        for r in results:
            prefix = r["site_id"]
            fh.write(f"{prefix}_insert_length\t{r['insert_length']}\n")
            fh.write(f"{prefix}_remaining_ns\t{r['remaining_ns']}\n")
            fh.write(f"{prefix}_assembly_rounds\t{r['rounds']}\n")
            fh.write(f"{prefix}_status\t{r['status']}\n")
            fh.write(f"{prefix}_verdict\t{r.get('verdict', 'NO_ASSEMBLY')}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 5: Targeted insert assembly "
                    "(soft-clip detection + k-mer extension + Pilon gap fill)")
    parser.add_argument("--junctions", default=None,
                        help="junctions.tsv from step 6 (fallback if soft-clip "
                             "detection finds nothing)")
    parser.add_argument("--host-bam", required=True,
                        help="Host-mapped BAM from step 7")
    parser.add_argument("--host-ref", required=True,
                        help="Host reference FASTA")
    parser.add_argument("--element-db", required=True,
                        help="Element database FASTA for annotation")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--max-rounds", type=int, default=8,
                        help="Max alternating extension+Pilon rounds")
    parser.add_argument("--seed-k", type=int, default=15,
                        help="k-mer size for seed extension")
    parser.add_argument("--recruit-k", type=int, default=25,
                        help="k-mer size for unmapped read recruitment")
    parser.add_argument("--gap-size", type=int, default=1000,
                        help="Initial N gap size for Pilon scaffold")
    parser.add_argument("--flank", type=int, default=5000,
                        help="Flank size for candidate read extraction (bp)")
    parser.add_argument("--min-clip", type=int, default=20,
                        help="Minimum soft-clip length for junction detection")
    parser.add_argument("--extra-element-db", type=Path, default=None,
                        help="Optional second FASTA for element annotation "
                             "(e.g., per-sample SPAdes construct contigs).")
    parser.add_argument("--common-payload-db", type=Path, default=None,
                        help="Always-on supplementary FASTA of common transgene "
                             "payload genes (bar, nptII, hpt, gusA, gfp, ...).")
    parser.add_argument("--junction-window", type=int, default=20,
                        help="Window (bp) around junction for allele phasing "
                             "(reads within this distance classified as junction-proximal)")
    parser.add_argument("--construct-ref", default=None,
                        help="Construct reference FASTA (for flanking detection)")
    parser.add_argument("--no-remote-blast", action="store_true",
                        help="Skip remote NCBI nt BLAST (use local element_db only)")
    parser.add_argument("--s03-r1", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--s03-r2", default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()

    step_dir = Path(args.outdir) / args.sample_name / STEP
    step_dir.mkdir(parents=True, exist_ok=True)

    host_bam = Path(args.host_bam)
    host_ref = Path(args.host_ref)
    element_db = Path(args.element_db)
    extra_db = args.extra_element_db
    common_payload_db = args.common_payload_db

    # Order: common_payload (shared) first, then per-sample extra.
    # A None entry is dropped naturally.
    extra_dbs = [p for p in [common_payload_db, extra_db] if p is not None]

    log(f"=== Step 5: Targeted Insert Assembly for {args.sample_name} ===")
    log(f"  Host BAM: {host_bam}")
    log(f"  Host ref: {host_ref}")
    log(f"  Element DB: {element_db}")
    if common_payload_db is not None:
        log(f"  Common payload DB: {common_payload_db}")
    if extra_db is not None:
        log(f"  Extra element DB: {extra_db}")

    # ---- Phase 1: Soft-clip junction detection ----
    sites = find_softclip_junctions(
        host_bam, host_ref, element_db, step_dir,
        min_clip=args.min_clip,
        extra_dbs=extra_dbs,
    )

    # Fallback to step 6 junctions if no sites found
    if not sites and args.junctions:
        junctions_path = Path(args.junctions)
        if junctions_path.exists() and junctions_path.stat().st_size > 0:
            log("No soft-clip sites found, falling back to legacy junctions file...")
            legacy_juncs = parse_legacy_junctions(junctions_path)
            if legacy_juncs:
                sites = legacy_junctions_to_sites(legacy_juncs, host_bam)
                log(f"  Converted {len(legacy_juncs)} legacy junctions → "
                    f"{len(sites)} sites")

    if not sites:
        log("No insertion sites found. Nothing to assemble.")
        write_stats(step_dir / "s05_stats.txt", args.sample_name, 0, 0, [])
        return

    log(f"\nPhase 1 found {len(sites)} insertion site(s):")
    for site in sites:
        log(f"  {site.site_id}: {site.host_chr}:{site.pos_5p}"
            f"{f'-{site.pos_3p}' if site.pos_3p else ''} "
            f"({site.confidence}, 5p={len(site.seed_5p)}bp, "
            f"3p={len(site.seed_3p)}bp)")

    # ---- Phase 1.5: Transgene-positive identification ----
    log("\n--- Phase 1.5: Transgene-positive identification ---")
    assembly_sites, skip_sites, tier_results, host_endo_ids = classify_site_tiers(
        sites, element_db, host_ref, step_dir,
        threads=args.threads,
        extra_transgene_dbs=extra_dbs,
    )
    write_tier_classification(
        tier_results, step_dir / "site_tier_classification.tsv",
    )

    if not assembly_sites:
        log("No transgene-positive sites found. Nothing to assemble.")
        write_stats(step_dir / "s05_stats.txt", args.sample_name, len(sites), 0, [])
        return

    sites = assembly_sites  # Only assemble transgene-positive sites

    # ---- Process each insertion site ----
    all_results = []
    for site in sites:
        log(f"\n{'=' * 60}")
        log(f"=== Processing {site.site_id}: "
            f"{site.host_chr}:{site.pos_5p} ===")

        # Phase 2: Extract candidate reads
        cand_r1 = step_dir / f"{site.site_id}_candidate_R1.fastq.gz"
        cand_r2 = step_dir / f"{site.site_id}_candidate_R2.fastq.gz"
        if not cand_r1.exists():
            n_cand = extract_candidate_reads(
                host_bam, site, cand_r1, cand_r2,
                flank=args.flank, threads=args.threads,
                min_clip=args.min_clip,
                junction_window=args.junction_window,
            )
            log(f"  Candidate reads: {n_cand:,} pairs")
        else:
            n_cand = 0
            opener = gzip.open if str(cand_r1).endswith(".gz") else open
            with opener(str(cand_r1), "rt") as fh:
                for _ in fh:
                    n_cand += 1
            n_cand = n_cand // 4
            log(f"  Candidate reads (cached): {n_cand:,} pairs")

        # Phase 3: Iterative assembly
        s03_r1 = Path(args.s03_r1) if args.s03_r1 else None
        s03_r2 = Path(args.s03_r2) if args.s03_r2 else None
        assembly_fa, rounds, status = assemble_insert(
            site=site,
            candidate_r1=cand_r1,
            candidate_r2=cand_r2,
            host_bam=host_bam,
            host_ref=host_ref,
            element_db=element_db,
            workdir=step_dir,
            threads=args.threads,
            max_rounds=args.max_rounds,
            ext_k=args.seed_k,
            recruit_k=args.recruit_k,
            gap_size=args.gap_size,
            s03_r1=s03_r1,
            s03_r2=s03_r2,
        )

        # Record result
        if assembly_fa and assembly_fa.exists():
            contigs = read_fasta(assembly_fa)
            if contigs:
                longest_name = max(contigs, key=lambda k: len(contigs[k]))
                longest_seq = contigs[longest_name]
                insert_len = len(longest_seq)
                insert_ns = longest_seq.upper().count("N")
            else:
                insert_len, insert_ns = 0, 0
                status = "no_assembly"
        else:
            insert_len, insert_ns = 0, 0
            status = "no_assembly"

        all_results.append({
            "site_id": site.site_id,
            "insert_length": insert_len,
            "remaining_ns": insert_ns,
            "rounds": rounds,
            "status": status,
        })

    # ---- Combine inserts ----
    combined_insert = step_dir / "insert_only.fasta"
    with open(combined_insert, "w") as fout:
        for site in sites:
            site_fa = step_dir / f"{site.site_id}_insert.fasta"
            if site_fa.exists():
                fout.write(open(site_fa).read())

    # ---- Phase 4: Annotation & Report ----
    if combined_insert.exists() and combined_insert.stat().st_size > 0:
        log(f"\n{'=' * 60}")
        log("=== Phase 4: Annotation ===")
        ann_tsv, border_tsv = annotate_insert(
            combined_insert, element_db, step_dir, args.sample_name,
            no_remote_blast=args.no_remote_blast,
            extra_dbs=extra_dbs,
        )

        # host_endo_ids returned directly from classify_site_tiers (no file I/O)
        log(f"  Using {len(host_endo_ids)} host-endogenous element IDs "
            f"for post-filter verdict")

        # Detect construct-flanking regions (host DNA in construct reference)
        construct_flanking: list[tuple[str, int, int]] = []
        if args.construct_ref:
            construct_path = Path(args.construct_ref)
            if construct_path.exists():
                log("  Detecting construct-flanking host regions...")
                construct_flanking = _find_construct_flanking_regions(
                    construct_path, host_ref, step_dir, threads=args.threads,
                )

        # Generate report for each site
        verdict_counts: Counter = Counter()
        verdict_by_site: dict[str, str] = {}
        for site in sites:
            site_fa = step_dir / f"{site.site_id}_insert.fasta"
            if site_fa.exists():
                result_info = next(
                    (r for r in all_results if r["site_id"] == site.site_id), None)
                if result_info:
                    _, verdict = generate_report(
                        site_fa, ann_tsv, border_tsv,
                        site, result_info["rounds"], result_info["status"],
                        step_dir,
                        host_endogenous_ids=host_endo_ids,
                        host_ref=host_ref,
                        element_db=element_db,
                        threads=args.threads,
                        construct_flanking=construct_flanking,
                    )
                    verdict_counts[verdict] += 1
                    verdict_by_site[site.site_id] = verdict

        # Attach verdict onto the per-site results so the final summary reflects it.
        for r in all_results:
            r["verdict"] = verdict_by_site.get(r["site_id"], "NO_ASSEMBLY")
        log(f"  Phase 4 verdicts: "
            f"{dict(verdict_counts) if verdict_counts else '(none)'}")

    # ---- Write stats ----
    write_stats(
        step_dir / "s05_stats.txt", args.sample_name,
        len(sites), 0, all_results,
    )

    log(f"\n{'=' * 60}")
    log(f"=== Step 5 complete ===")
    log(f"Output: {step_dir}")
    for r in all_results:
        v = r.get("verdict", "NO_ASSEMBLY")
        log(f"  {r['site_id']}: {r['insert_length']:,}bp, "
            f"{r['remaining_ns']} Ns, {r['rounds']} rounds → {r['status']} [{v}]")


if __name__ == "__main__":
    main()
