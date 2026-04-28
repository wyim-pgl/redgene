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
    from s05.classify import (
        _batch_check_element_hits,
        _filter_host_endogenous,
        _should_replace,
        _SRC_TIER,
        _TIER2_SRCS,
        _UNKNOWN_SRC_WARNED,
        classify_site_tiers,
        write_tier_classification,
    )
    from s05.site_discovery import (
        _build_consensus,
        _batch_check_maps_to_host,
        find_softclip_junctions,
        _apply_mask_bed,
        MASKED_SOURCE_TAG,
        parse_legacy_junctions,
        _extract_seeds_at_positions,
        legacy_junctions_to_sites,
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
    from s05.classify import (
        _batch_check_element_hits,
        _filter_host_endogenous,
        _should_replace,
        _SRC_TIER,
        _TIER2_SRCS,
        _UNKNOWN_SRC_WARNED,
        classify_site_tiers,
        write_tier_classification,
    )
    from s05.site_discovery import (
        _build_consensus,
        _batch_check_maps_to_host,
        find_softclip_junctions,
        _apply_mask_bed,
        MASKED_SOURCE_TAG,
        parse_legacy_junctions,
        _extract_seeds_at_positions,
        legacy_junctions_to_sites,
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
# Phase 1 (Site Discovery) functions moved to scripts/s05/site_discovery.py
# (Issue #4 Session 5) and re-imported at the top of this module for
# backward compatibility with existing callers.
# Note: _extract_seeds_at_positions is included in site_discovery.py
# (not read_extraction.py) because it is solely called by
# legacy_junctions_to_sites (Phase 1). Moving it to read_extraction
# would require a stage-1→stage-2 import — a DAG violation.
# ---------------------------------------------------------------------------
# Phase 1.5 (Classify) functions moved to scripts/s05/classify.py
# (Issue #4 Session 4) and re-imported at the top of this module for
# backward compatibility. _TIER2_SRCS is derived from _SRC_TIER to stay
# in sync with gmo_combined_db_v2.fa cd-hit clustering; see
# scripts/s05/classify.py for the full audit trail comment.
# ---------------------------------------------------------------------------

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
