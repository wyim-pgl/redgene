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
import json
import pickle
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

# Local sub-modules (Issue #3 full wire-in: compute_verdict + FilterEvidence).
# Session 1 of Issue #4 also re-imports primitives (log, revcomp, dataclasses,
# FASTA/FASTQ helpers) from scripts.s05.primitives so downstream submodules
# share a single definition without cycling back through this entrypoint.
try:
    from s05.verdict import VerdictRules, FilterEvidence, compute_verdict, _apply_canonical_override
    from s05.config_loader import load_verdict_rules
    from s05.primitives import (
        STEP,
        log,
        revcomp,
        read_fasta,
        write_fasta,
        _parse_src_tag,
        _read_fq_seqs,
        JunctionCluster,
        InsertionSite,
        LegacyJunction,
        TierResult,
    )
    from s05.filter_b_flank import (
        CONSTRUCT_FLANK_PIDENT,
        CONSTRUCT_FLANK_MIN_LEN,
        CONSTRUCT_FLANK_SLOP,
        _find_construct_flanking_regions,
        _site_overlaps_flanking,
    )
    from s05.filter_c_chimeric import (
        CHIMERIC_MIN_PIDENT,
        CHIMERIC_MIN_OFFTARGET_BP,
        _check_chimeric_assembly,
    )
    from s05.filter_a_host import (
        INSERT_HOST_MIN_PIDENT,
        _blast_insert_vs_host,
    )
    from s05.filter_d_altlocus import (
        CONSTRUCT_HOST_MIN_COMBINED,
        CONSTRUCT_MIN_FRACTION,
        CONSTRUCT_HOST_MIN_PIDENT,
        _check_construct_host_coverage,
    )
    from s05.annotation import (
        _parse_blast6,
        _run_local_blast,
        _run_remote_blast,
        _merge_annotations,
        annotate_insert,
    )
    from s05.assembly import (
        StrandAwareSeedExtender,
        _check_merge,
        _minimap2_extend,
        _ssake_extend,
        _vote_extension,
        _write_pool_fastq,
        assemble_insert,
        check_host_termination,
        extract_foreign_reads,
        pilon_fill,
        recruit_by_kmer,
        refine_with_foreign_reads,
    )
except ImportError:  # pragma: no cover - allow standalone `python scripts/s05_insert_assembly.py`
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from s05.verdict import VerdictRules, FilterEvidence, compute_verdict, _apply_canonical_override
    from s05.config_loader import load_verdict_rules
    from s05.primitives import (
        STEP,
        log,
        revcomp,
        read_fasta,
        write_fasta,
        _parse_src_tag,
        _read_fq_seqs,
        JunctionCluster,
        InsertionSite,
        LegacyJunction,
        TierResult,
    )
    from s05.filter_b_flank import (
        CONSTRUCT_FLANK_PIDENT,
        CONSTRUCT_FLANK_MIN_LEN,
        CONSTRUCT_FLANK_SLOP,
        _find_construct_flanking_regions,
        _site_overlaps_flanking,
    )
    from s05.filter_c_chimeric import (
        CHIMERIC_MIN_PIDENT,
        CHIMERIC_MIN_OFFTARGET_BP,
        _check_chimeric_assembly,
    )
    from s05.filter_a_host import (
        INSERT_HOST_MIN_PIDENT,
        _blast_insert_vs_host,
    )
    from s05.filter_d_altlocus import (
        CONSTRUCT_HOST_MIN_COMBINED,
        CONSTRUCT_MIN_FRACTION,
        CONSTRUCT_HOST_MIN_PIDENT,
        _check_construct_host_coverage,
    )
    from s05.annotation import (
        _parse_blast6,
        _run_local_blast,
        _run_remote_blast,
        _merge_annotations,
        annotate_insert,
    )
    from s05.assembly import (
        StrandAwareSeedExtender,
        _check_merge,
        _minimap2_extend,
        _ssake_extend,
        _vote_extension,
        _write_pool_fastq,
        assemble_insert,
        check_host_termination,
        extract_foreign_reads,
        pilon_fill,
        recruit_by_kmer,
        refine_with_foreign_reads,
    )

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
# INSERT_HOST_MIN_PIDENT moved to scripts/s05/filter_a_host.py (Issue #4
# Session 2) and re-imported at the top of this module for backward compat.

# Construct-flanking filter constants (CONSTRUCT_FLANK_PIDENT, _MIN_LEN, _SLOP)
# moved to scripts/s05/filter_b_flank.py (Issue #4 Session 1) and re-imported
# at the top of this module for backward compatibility.

# Multi-locus chimeric filter constants (CHIMERIC_MIN_PIDENT,
# CHIMERIC_MIN_OFFTARGET_BP) moved to scripts/s05/filter_c_chimeric.py
# (Issue #4 Session 1) and re-imported at the top of this module for
# backward compatibility.

# (Filter D constants CONSTRUCT_HOST_MIN_COMBINED, CONSTRUCT_MIN_FRACTION,
# CONSTRUCT_HOST_MIN_PIDENT moved to scripts/s05/filter_d_altlocus.py,
# Issue #4 Session 2, and re-imported at the top of this module.)

# UNKNOWN → FALSE_POSITIVE auto-reclassification: if insert has no element
# annotations but is mostly host DNA with negligible construct match,
# it is host genomic DNA, not a real transgene insertion.
UNKNOWN_HOST_MIN_FRACTION = 0.85   # min host fraction to classify as host-only
UNKNOWN_MAX_CONSTRUCT_FRAC = 0.05  # max construct fraction for host-only


# _apply_canonical_override moved to scripts/s05/verdict.py and re-imported
# here for backward compatibility with tests/test_canonical_override.py.


# ---------------------------------------------------------------------------
# Primitives moved to scripts/s05/primitives.py (Issue #4 Session 1).
# ``revcomp``, ``read_fasta``, ``write_fasta``, ``_read_fq_seqs``, and the four
# dataclasses (JunctionCluster, InsertionSite, LegacyJunction, TierResult) are
# imported at the top of this module and re-exported below for test back-compat.
# ---------------------------------------------------------------------------


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
# T10: Pre-mask BED intersect (host-endogenous ortholog regions)
# ---------------------------------------------------------------------------

#: Tier value written onto sites that intersect a T9 host-endogenous mask
#: region.  Centralised here so downstream consumers (generate_report audit
#: appender, tests, future modules) can import the tag instead of duplicating
#: the string literal.  Changing this constant is a wire-format break — update
#: the audit TSVs and regression fixtures in lock-step.
MASKED_SOURCE_TAG = "FALSE_NEGATIVE_MASKED"


def _apply_mask_bed(sites: list, bed_path: Path) -> list:
    """Tag sites whose ``(chr, pos)`` falls inside a BED mask region.

    Sites are NOT dropped — their ``tier`` is changed to
    :data:`MASKED_SOURCE_TAG` (``FALSE_NEGATIVE_MASKED``) and ``mask_element``
    holds the element_name of the overlapping BED region.  This preserves the
    audit trail for regulatory review while preventing CANDIDATE promotion on
    host-endogenous ortholog regions (T9 E-01/E-02).

    Interval semantics
    ------------------
    Regions are treated as **half-open** ``[start, end)`` — i.e. ``start`` is
    inclusive and ``end`` is exclusive.  This matches BED convention and the
    ``start <= pos < end`` check below; a site at ``pos == end`` is NOT
    considered masked.

    BED format (from ``scripts/build_element_mask_bed.sh``):
    ``chr<TAB>start<TAB>end<TAB>name[<TAB>score]``.  An empty or missing BED
    file (including a BED with zero parseable regions) is a no-op.

    Parameters
    ----------
    sites
        A list of site records.  The list itself is returned (not copied); the
        records inside are mutated in place — see *Mutation contract* below.
    bed_path
        Path to the pre-computed mask BED.

    Mutation contract
    -----------------
    Accepts heterogeneous site records and mutates the *original* objects:

    * ``dict`` with ``chr``/``pos`` keys (test + pipeline JSON form) — sets
      ``s["tier"]`` and ``s["mask_element"]``.
    * Any object with ``host_chr`` + ``pos_5p`` attributes (e.g. the
      ``InsertionSite`` dataclass) — sets ``tier`` and ``mask_element`` via
      ``setattr``.  When ``pos_5p`` is falsy the helper falls back to
      ``pos_3p``.

    Callers that need the prior value should snapshot it before invoking this
    function; the helper does not save or restore previous ``tier`` values.

    Complexity
    ----------
    Regions are grouped by chromosome into a list sorted by ``start`` plus a
    parallel ``start`` array; per-site lookup uses :func:`bisect.bisect_right`
    and scans only intervals whose ``start <= pos``.  Runtime is therefore
    ``O(R log R + S log R)`` for ``R`` regions and ``S`` sites, instead of the
    ``O(S * R)`` worst case of a linear scan.
    """
    import bisect

    if not bed_path.exists() or bed_path.stat().st_size == 0:
        return sites
    raw_regions: dict[str, list[tuple[int, int, str]]] = {}
    for raw in bed_path.read_text().splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        cols = raw.split("\t")
        if len(cols) < 3:
            continue
        try:
            start = int(cols[1])
            end = int(cols[2])
        except ValueError:
            continue
        chrom = cols[0]
        name = cols[3] if len(cols) > 3 else "masked"
        raw_regions.setdefault(chrom, []).append((start, end, name))
    if not raw_regions:
        return sites
    # Sort regions by start and stash a parallel starts array for bisect.
    # NOTE: the intervals may still overlap each other; bisect narrows the
    # scan to candidates with start <= pos, then we check end > pos.
    index: dict[str, tuple[list[int], list[tuple[int, int, str]]]] = {}
    for chrom, ivs in raw_regions.items():
        ivs.sort(key=lambda t: t[0])
        index[chrom] = ([iv[0] for iv in ivs], ivs)
    for s in sites:
        if isinstance(s, dict):
            chrom = s.get("chr")
            pos = s.get("pos")
        else:
            chrom = getattr(s, "host_chr", None)
            pos = getattr(s, "pos_5p", None)
            if not pos:
                pos = getattr(s, "pos_3p", None)
        if chrom is None or pos is None:
            continue
        entry = index.get(chrom)
        if entry is None:
            continue
        starts, ivs = entry
        hi = bisect.bisect_right(starts, pos)
        # Walk backwards through candidates whose start <= pos; because the
        # regions in the mask BED are short (<= a few kb) and only a handful
        # overlap any given pos, this inner loop is effectively O(1).
        for idx in range(hi - 1, -1, -1):
            start, end, name = ivs[idx]
            if end <= pos:
                # Half-open: pos >= end means outside; keep scanning earlier
                # intervals because an earlier start may still extend past pos.
                continue
            # start <= pos < end -> match.
            if isinstance(s, dict):
                s["tier"] = MASKED_SOURCE_TAG
                s["mask_element"] = name
            else:
                setattr(s, "tier", MASKED_SOURCE_TAG)
                setattr(s, "mask_element", name)
            break
    return sites


# ---------------------------------------------------------------------------
# Transgene-positive identification
# ---------------------------------------------------------------------------
# Strategy: Instead of asking "is this clip in the host?" (fails with wrong
# cultivar reference), ask "is this clip from a transgene?"
# Uses blastn-short against transgene_db (element DB + UniVec vectors).
# Any transgene hit = assembly target. No hit = skip.
#
# NOTE: TierResult dataclass moved to scripts/s05/primitives.py
# (Issue #4 Session 1) and re-imported at the top of this module.


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


_SRC_TIER = {
    "element_db":    2,
    "payload":       2,
    "sample_contig": 2,
    "univec":        1,
}

# Tier-2 sources collectively count as "element hits" for transgene_positive
# classification. Derived from _SRC_TIER so adding a new tier-2 source (e.g.,
# future 'crispr_ref') stays in sync automatically. Required because cd-hit
# -c 0.95 clustering of gmo_combined_db_v2.fa absorbed multiple element_db
# amplicons into payload-tagged representatives (cluster 0 P-CaMV35S, 35
# nptII, 54 hpt, 71 T-ocs, 77 P-nos in element_db/gmo_combined_db_v2.fa.clstr);
# a literal 'element_db' equality check at classify_site_tiers would drop
# rice_G281 Chr3:16,439,674 (CaMV35S-containing) -> AC-1 regression.
_TIER2_SRCS = frozenset(k for k, v in _SRC_TIER.items() if v >= 2)

# Module-level cache so we only warn once per unknown tag per run.
_UNKNOWN_SRC_WARNED: set[str] = set()


# _parse_src_tag moved to scripts/s05/primitives.py
# (Issue #4 Session 3) and re-imported at the top of this module for
# backward compatibility with existing callers.


def _should_replace(existing: dict | None, new_src: str, new_bit: float) -> bool:
    """Tag-tiered merge priority for classify_site_tiers `hits` map (Task T5).

    Tier-2 sources (element_db / payload / sample_contig) collectively beat
    tier-1 (univec). Within a tier, a strictly greater bitscore wins;
    ties go to the incumbent (matches the original `>` semantics and the
    BUG-3 regression guard in test_extra_element_db).

    The `existing` dict may use either the new `src`/`bit` schema or the
    legacy `source`/`bitscore` schema - both are supported so the historic
    callers in classify_site_tiers keep working while migration proceeds.
    """
    # I-1: warn once per run when an unknown tag shows up. Treated as tier 0
    # so it will never beat a known tag, which is the safe default.
    if new_src and new_src not in _SRC_TIER and new_src not in _UNKNOWN_SRC_WARNED:
        _UNKNOWN_SRC_WARNED.add(new_src)
        print(f"[s05] warn: unknown src tag {new_src!r}, treating as tier 0",
              file=sys.stderr)
    if existing is None:
        return True
    old_src = existing.get("src", existing.get("source", ""))
    old_bit = existing.get("bit", existing.get("bitscore", 0.0))
    new_tier = _SRC_TIER.get(new_src, 0)
    old_tier = _SRC_TIER.get(old_src, 0)
    if new_tier > old_tier:
        return True
    if new_tier < old_tier:
        return False
    # same tier -> strict '>' bitscore
    return new_bit > old_bit


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

                sseqid = cols[1]
                # T5: 4-way source tag. v2 DB headers carry a |src=<tag>
                # suffix (payload / element_db / sample_contig); legacy
                # univec|... headers get mapped to 'univec', everything
                # else to 'element_db'.
                default_main = (
                    "univec" if sseqid.startswith("univec|") else "element_db"
                )
                element_id, src_tag = _parse_src_tag(sseqid, default_main)
                if _should_replace(hits.get(qname), src_tag, bitscore):
                    hits[qname] = {
                        "element": element_id,
                        "identity": pident,
                        "aln_length": aln_len,
                        "bitscore": bitscore,
                        "source": src_tag,
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
        # T5 fallback tag when a header lacks |src=<tag>: per-sample SPAdes
        # contig FASTAs from s04b are labelled 'sample_contig'; everything
        # else (legacy common_payload, ad-hoc FASTAs) defaults to
        # 'element_db' so tier-2 priority still holds.
        default_tag = "sample_contig" if "contigs" in edb.name else "element_db"
        log(f"  BLASTing clips against extra transgene DB "
            f"({edb.name}, default_src={default_tag}, blastn-short)...")
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
                    sseqid = cols[1]
                    element_id, src_tag = _parse_src_tag(sseqid, default_tag)
                    if _should_replace(hits.get(qname), src_tag, bitscore):
                        hits[qname] = {
                            "element": element_id,
                            "identity": pident,
                            "aln_length": aln_len,
                            "bitscore": bitscore,
                            "source": src_tag,
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

        # Require at least one tier-2 element hit (not univec-only).
        # UniVec-only matches at this step are typically short (20-31bp)
        # alignments that match native plant DNA by chance. Real T-DNA
        # insertions always have at least one characteristic element
        # (promoter, selection marker, terminator) matching a tier-2 source.
        #
        # C-1 fix: tier-2 sources (element_db / payload / sample_contig) all
        # count. After cd-hit -c 0.95 clustering, former element_db amplicons
        # may be absorbed into payload-tagged representatives (e.g., CaMV35S
        # cluster 0 in gmo_combined_db_v2.fa.clstr). A literal "element_db"
        # check would miss rice_G281 Chr3:16,439,674 CaMV35S hits -> AC-1
        # regression. Use the _TIER2_SRCS module constant derived from
        # _SRC_TIER for a single source of truth.
        has_element_hit = (hit_5p.get("source") in _TIER2_SRCS) or \
                           (hit_3p.get("source") in _TIER2_SRCS)
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
# Phase 3 (Assembly) functions moved to scripts/s05/assembly.py
# (Issue #4 Session 4) and re-imported at the top of this module for
# backward compatibility with existing callers (Phase 2 → Phase 3 →
# Phase 4 → classify_site_tiers chain remains bit-identical).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Phase 4 (Annotation) functions moved to scripts/s05/annotation.py
# (Issue #4 Session 3) and re-imported at the top of this module for
# backward compatibility with existing callers. The monolith's
# ``classify_site_tiers`` continues to consume the resulting hits.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Post-assembly host-fraction filter
# ---------------------------------------------------------------------------
# Filter B helpers (_find_construct_flanking_regions, _site_overlaps_flanking)
# moved to scripts/s05/filter_b_flank.py (Issue #4 Session 1) and re-imported
# at the top of this module for backward compatibility with existing callers.


# _check_chimeric_assembly moved to scripts/s05/filter_c_chimeric.py
# (Issue #4 Session 1) and re-imported at the top of this module for
# backward compatibility with existing callers.


# _check_construct_host_coverage moved to scripts/s05/filter_d_altlocus.py
# (Issue #4 Session 2) and re-imported at the top of this module for
# backward compatibility with existing callers.

# _blast_insert_vs_host moved to scripts/s05/filter_a_host.py
# (Issue #4 Session 2) and re-imported at the top of this module for
# backward compatibility with existing callers.


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
    rules: VerdictRules | None = None,
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

    # ---- Build element-level evidence (no I/O) ----
    host_endo = host_endogenous_ids or set()
    unique_elems = set(elem_counts)
    foreign_elems: set[str] = set()
    endo_elems: set[str] = set()
    for elem in unique_elems:
        if elem in host_endo:
            endo_elems.add(elem)
        else:
            foreign_elems.add(elem)
    # sources_by_element: last-seen source tag per element name.
    sources_by_element: dict[str, str] = {
        elem: source for _, _, elem, _, source in elements
    }

    # ---- Gather BLAST filter evidence (lazy: only run when useful) ----
    # Evidence is collected for CANDIDATE sites and UNKNOWN sites (host-only
    # reclassification). All four filter values default to 0/None so that
    # compute_verdict receives a complete FilterEvidence for every path.
    host_fraction = 0.0
    host_bp = 0
    largest_gap = 0
    flanking_hit_str = ""           # human-readable for report diagnostics
    flanking_hit_tup: tuple[str, int, int] | None = None  # structured for verdict
    is_chimeric = False
    off_target_chrs: list[tuple[str, int]] = []
    construct_frac = 0.0
    combined_frac = 0.0

    # Determine whether this site is worth running BLASTs on.
    # Sites that are all-host-endogenous skip BLAST, matching old behaviour;
    # canonical_triplet override (highest priority in compute_verdict) still
    # works because ev.host_fraction defaults to 0.0 < cand_host_fraction_max.
    _needs_blast = bool(unique_elems and foreign_elems) or not unique_elems
    # "not unique_elems" → UNKNOWN path, needs BLAST for host-only reclassification.
    # "foreign_elems" → CANDIDATE path, needs BLAST for all four FP filters.

    if _needs_blast and host_ref is not None:
        # Filter A evidence: host-fraction + largest non-host gap
        host_fraction, host_bp, _, largest_gap = _blast_insert_vs_host(
            insert_fasta, host_ref, output_dir, threads=threads,
        )
        if not unique_elems:
            log(f"  UNKNOWN reclassification: BLASTed {insert_fasta.name} vs host")

    if _needs_blast and construct_flanking:
        # Filter B evidence: find first overlapping flanking region
        site_pos = site.pos_5p
        _, flanking_hit_str = _site_overlaps_flanking(
            site.host_chr, site_pos, construct_flanking,
        )
        for chrom, lo, hi in construct_flanking:
            slop = CONSTRUCT_FLANK_SLOP
            if chrom == site.host_chr and lo - slop <= site_pos <= hi + slop:
                # Store slop-expanded coords so compute_verdict's simple
                # lo <= site_pos <= hi boundary check is equivalent.
                flanking_hit_tup = (chrom, lo - slop, hi + slop)
                break

    if _needs_blast and host_ref is not None:
        # Filter C evidence: off-target chromosome spans
        is_chimeric, off_target_chrs = _check_chimeric_assembly(
            insert_fasta, host_ref, site.host_chr, output_dir, threads=threads,
        )

    if _needs_blast and host_ref is not None and element_db is not None:
        # Filter D evidence: construct + host combined coverage
        _, construct_frac, _, combined_frac = _check_construct_host_coverage(
            insert_fasta, element_db,
            host_fraction, host_bp, insert_len, n_count,
            output_dir, threads=threads,
        )

    # ---- Delegate verdict to compute_verdict (Issue #3 full wire-in) ----
    # Build a FilterEvidence bundle from all gathered evidence and call the
    # pure function.  This replaces the old inline if/else verdict chain plus
    # the _apply_canonical_override post-hoc call.
    _ev = FilterEvidence(
        elements=list(unique_elems),
        host_bp=host_bp,
        host_fraction=host_fraction,
        largest_gap=largest_gap,
        flanking_hit=flanking_hit_tup,
        off_target_chrs=off_target_chrs,
        construct_frac=construct_frac,
        combined_frac=combined_frac,
        is_chimeric=is_chimeric,
        site_chr=site.host_chr,
        site_pos=site.pos_5p,
        matched_canonical=unique_elems,   # canonical triplet check uses all elements
        sources_by_element=sources_by_element,
        host_endogenous_elements=endo_elems,
    )
    _verdict_rules = rules if rules is not None else VerdictRules()
    verdict, verdict_reason = compute_verdict(_ev, _verdict_rules)

    # For report diagnostics: restore the human-readable flanking string if
    # compute_verdict returned a flanking FP (flanking_hit_str may differ from
    # the structured tuple representation used internally).
    # The reporting block below uses flanking_hit_str for display; verdict_reason
    # already contains the full human-readable text from compute_verdict.

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
        if flanking_hit_str:
            fh.write(f"Construct-flanking overlap: site overlaps {flanking_hit_str}\n")
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
    parser.add_argument("--max-rounds", type=int, default=5,
                        help="Max alternating extension+Pilon rounds (T3: reduced from 8 per max_rounds_study.md)")
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
    parser.add_argument(
        "--mask-bed",
        default=None,
        help="BED file of host-endogenous ortholog regions (T9/T10). "
             "Sites overlapping a region are tagged FALSE_NEGATIVE_MASKED "
             "(not dropped) after Phase 1 discovery, preventing CANDIDATE "
             "promotion while preserving the audit trail.",
    )
    parser.add_argument("--s03-r1", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--s03-r2", default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        "--config",
        default=None,
        help="Path to config.yaml for verdict rules / canonical triplets "
             "(Issue #3 wire-in). Missing file falls back to DEFAULT_TRIPLETS.",
    )
    # --- T8: per-site SLURM array fan-out (v1.0 임시안) -----------------------
    # v1.1 에서 full module split 완료 후 이 flag 제거 예정.
    parser.add_argument(
        "--phase",
        choices=["all", "1_1.5", "2_3", "4"],
        default="all",
        help="Run specific phase(s). all=default sequential (backwards compat). "
             "1_1.5=site discovery + transgene-positive classify, dumps "
             "positive_sites.json and exits. 2_3=per-site candidate read "
             "extraction + iterative assembly (requires --site-id). "
             "4=annotate + per-site report (consumes per-site "
             "<site_id>_insert.fasta dropped by 2_3 runs).",
    )
    parser.add_argument(
        "--site-id",
        default=None,
        help="For --phase=2_3: run assembly for only the given site_id "
             "(as stored in positive_sites.json).",
    )
    args = parser.parse_args()

    if args.phase == "2_3" and not args.site_id:
        parser.error("--site-id is required when --phase=2_3")

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
    log(f"  Phase: {args.phase}"
        + (f" (site_id={args.site_id})" if args.site_id else ""))
    log(f"  Host BAM: {host_bam}")
    log(f"  Host ref: {host_ref}")
    log(f"  Element DB: {element_db}")
    if common_payload_db is not None:
        log(f"  Common payload DB: {common_payload_db}")
    if extra_db is not None:
        log(f"  Extra element DB: {extra_db}")

    # Paths used by T8 fan-out contract between phases.
    positive_sites_json = step_dir / "positive_sites.json"
    positive_sites_pkl = step_dir / "positive_sites.pkl"

    # ---- Phase 1 + 1.5 (site discovery + transgene-positive classify) ----
    # Runs when --phase is 'all' or '1_1.5'. For 2_3 / 4, we load the
    # pickled state that a previous 1_1.5 (or all) invocation wrote.
    sites: list[InsertionSite] = []
    host_endo_ids: set[str] = set()

    if args.phase in ("all", "1_1.5"):
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
            # Still emit an empty positive_sites.json so the array wrapper
            # can detect 'nothing to fan out' cleanly.
            positive_sites_json.write_text("[]\n")
            positive_sites_pkl.write_bytes(pickle.dumps(
                {"sites": [], "host_endo_ids": set()}
            ))
            return

        log(f"\nPhase 1 found {len(sites)} insertion site(s):")
        for site in sites:
            log(f"  {site.site_id}: {site.host_chr}:{site.pos_5p}"
                f"{f'-{site.pos_3p}' if site.pos_3p else ''} "
                f"({site.confidence}, 5p={len(site.seed_5p)}bp, "
                f"3p={len(site.seed_3p)}bp)")

        # ---- T10: pre-mask BED intersect (host-endogenous ortholog regions) ----
        # Sites inside a T9 mask region are tagged FALSE_NEGATIVE_MASKED and
        # excluded from Phase 1.5 CANDIDATE classification. They still appear
        # in site_tier_classification.tsv for regulatory audit trail.
        masked_sites: list[InsertionSite] = []
        if args.mask_bed:
            bed_path = Path(args.mask_bed)
            _apply_mask_bed(sites, bed_path)
            masked_sites = [
                s for s in sites
                if getattr(s, "tier", None) == MASKED_SOURCE_TAG
            ]
            sites = [
                s for s in sites
                if getattr(s, "tier", None) != MASKED_SOURCE_TAG
            ]
            print(
                f"[Phase 1] mask-bed {bed_path.name}: "
                f"{len(masked_sites)}/{len(masked_sites) + len(sites)} "
                f"sites tagged {MASKED_SOURCE_TAG}",
                file=sys.stderr,
            )
            for s in masked_sites:
                log(f"  MASKED: {s.site_id} {s.host_chr}:{s.pos_5p} "
                    f"→ {getattr(s, 'mask_element', 'masked')}")

        # ---- Phase 1.5: Transgene-positive identification ----
        log("\n--- Phase 1.5: Transgene-positive identification ---")
        assembly_sites, skip_sites, tier_results, host_endo_ids = classify_site_tiers(
            sites, element_db, host_ref, step_dir,
            threads=args.threads,
            extra_transgene_dbs=extra_dbs,
        )
        # T10: append FALSE_NEGATIVE_MASKED entries to tier_results for audit
        for s in masked_sites:
            tier_results.append(TierResult(
                site_id=s.site_id,
                chrom=s.host_chr,
                pos=s.pos_5p or s.pos_3p or 0,
                transgene_positive=False,
                clip_5p_len=len(s.seed_5p) if s.seed_5p else 0,
                clip_3p_len=len(s.seed_3p) if s.seed_3p else 0,
                hit_5p=getattr(s, "mask_element", ""),
                hit_5p_source=MASKED_SOURCE_TAG,
                hit_3p=getattr(s, "mask_element", ""),
                hit_3p_source=MASKED_SOURCE_TAG,
            ))
        write_tier_classification(
            tier_results, step_dir / "site_tier_classification.tsv",
        )

        if not assembly_sites:
            log("No transgene-positive sites found. Nothing to assemble.")
            write_stats(step_dir / "s05_stats.txt", args.sample_name, len(sites), 0, [])
            positive_sites_json.write_text("[]\n")
            positive_sites_pkl.write_bytes(pickle.dumps(
                {"sites": [], "host_endo_ids": host_endo_ids}
            ))
            return

        sites = assembly_sites  # Only assemble transgene-positive sites

        # ---- T8: persist Phase 1/1.5 state for array fan-out -----------------
        # positive_sites.json is the *public* contract for submit_s05_array.sh
        # to enumerate site_ids for --array indexing.  positive_sites.pkl is
        # the *internal* rehydration file for Phase 2_3 / 4 workers so they
        # can run against the exact same InsertionSite objects (including
        # seed_5p/3p consensus clips and JunctionCluster detail).
        positive_sites_json.write_text(json.dumps(
            [
                {
                    "site_id": s.site_id,
                    "chr": s.host_chr,
                    "pos": s.pos_5p,
                    "pos_3p": s.pos_3p,
                    "confidence": s.confidence,
                }
                for s in sites
            ],
            indent=2,
        ) + "\n")
        # SECURITY: pickle is consumed only from the same sample's step_dir
        # (same user writes and reads). Never rehydrate from untrusted sources.
        # v1.1 module split will replace with typed JSON + dataclasses loader.
        positive_sites_pkl.write_bytes(pickle.dumps(
            {"sites": sites, "host_endo_ids": host_endo_ids}
        ))
        log(f"  [T8] wrote {len(sites)} positive sites → "
            f"{positive_sites_json.name} (+ .pkl)")

        if args.phase == "1_1.5":
            log(f"\n[Phase 1.5] early return — "
                f"{len(sites)} positive sites saved for per-site fan-out.")
            return

    else:
        # phase in {"2_3", "4"} — rehydrate from pickle written by a prior
        # 1_1.5 / all run.  Fail loud if the pickle is missing so the user
        # notices they skipped Phase 1.5.
        if not positive_sites_pkl.exists():
            sys.exit(
                f"ERROR: --phase {args.phase} requires a prior "
                f"--phase 1_1.5 (or --phase all) run to have produced "
                f"{positive_sites_pkl} — not found."
            )
        loaded = pickle.loads(positive_sites_pkl.read_bytes())
        sites = loaded["sites"]
        host_endo_ids = loaded["host_endo_ids"]
        log(f"  [T8] loaded {len(sites)} positive sites + "
            f"{len(host_endo_ids)} host-endogenous IDs from {positive_sites_pkl.name}")

    # ---- Phase 2 + 3: Per-site candidate extraction + iterative assembly ----
    all_results: list[dict] = []
    if args.phase in ("all", "2_3"):
        if args.phase == "2_3":
            # Filter to single site for array task.
            matching = [s for s in sites if s.site_id == args.site_id]
            if not matching:
                sys.exit(
                    f"ERROR: --site-id {args.site_id!r} not in "
                    f"positive_sites.json (have: "
                    f"{[s.site_id for s in sites]})"
                )
            sites_to_run = matching
        else:
            sites_to_run = sites

        for site in sites_to_run:
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

            result = {
                "site_id": site.site_id,
                "insert_length": insert_len,
                "remaining_ns": insert_ns,
                "rounds": rounds,
                "status": status,
            }
            all_results.append(result)

            # T8: per-site runs drop a small JSON sidecar so Phase 4
            # can rebuild all_results without re-running assembly.
            if args.phase == "2_3":
                (step_dir / f"{site.site_id}_result.json").write_text(
                    json.dumps(result, indent=2) + "\n"
                )

        if args.phase == "2_3":
            log(f"\n[Phase 2_3] done for site {args.site_id}. "
                f"Phase 4 collects results after the full array completes.")
            return

    # ---- Phase 4: Annotation & Report ----
    # For --phase=4 we rehydrate all_results from per-site sidecars dropped
    # by each --phase=2_3 worker.  For --phase=all we already built
    # all_results above.
    if args.phase == "4":
        for site in sites:
            sidecar = step_dir / f"{site.site_id}_result.json"
            if sidecar.exists():
                all_results.append(json.loads(sidecar.read_text()))
            else:
                log(f"  [T8 Phase 4] missing sidecar for {site.site_id}, "
                    f"marking as no_assembly")
                all_results.append({
                    "site_id": site.site_id,
                    "insert_length": 0,
                    "remaining_ns": 0,
                    "rounds": 0,
                    "status": "no_assembly",
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
        # Issue #3 wire-in: load canonical-triplet + threshold rules from config
        # (falls back to DEFAULT_TRIPLETS when --config is omitted / missing).
        rules = load_verdict_rules(
            Path(args.config) if args.config else Path("config.yaml"),
            args.sample_name,
        )
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
                        rules=rules,
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
