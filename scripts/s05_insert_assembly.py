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
    from s05.read_extraction import (
        extract_candidate_reads,
        extract_unmapped_paired,
    )
    from s05.report import (
        generate_report,
        write_stats,
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
    from s05.read_extraction import (
        extract_candidate_reads,
        extract_unmapped_paired,
    )
    from s05.report import (
        generate_report,
        write_stats,
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
# Phase 2 (Read Extraction) functions moved to scripts/s05/read_extraction.py
# (Issue #4 Session 5) and re-imported at the top of this module.
# ---------------------------------------------------------------------------

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



# ---------------------------------------------------------------------------
# Phase 5 (Report + Stats) functions moved to scripts/s05/report.py
# (Issue #4 Session 6) and re-imported at the top of this module.
# ---------------------------------------------------------------------------

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
