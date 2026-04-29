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
    from s05.fanout_orchestrator import main as _fanout_main
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
    from s05.fanout_orchestrator import main as _fanout_main

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
# Main entrypoint (orchestration moved to scripts/s05/fanout_orchestrator.py)
# ---------------------------------------------------------------------------
# T8 fan-out state files (positive_sites.json + positive_sites.pkl) are
# written/read by fanout_orchestrator._run_phase_1_1_5.  The pickle is
# consumed only from the same sample's step_dir (trusted path, same user).
# v1.1 module split will replace pickle with typed JSON + dataclasses loader.
# ---------------------------------------------------------------------------

def main() -> None:
    """Thin entrypoint — full orchestration in scripts/s05/fanout_orchestrator."""
    from s05.fanout_orchestrator import main as _orchestrator_main
    _orchestrator_main()

if __name__ == "__main__":
    main()
