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

After Issue #4 Session 6 this file is a thin entrypoint shim. All
orchestration lives in ``scripts/s05/fanout_orchestrator.py``. Phase-specific
logic lives in the ``scripts/s05/`` package. The re-export block below
preserves backward compatibility for callers that do:
  ``from scripts.s05_insert_assembly import <symbol>``
"""

from __future__ import annotations

import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Sub-module re-exports (backward compatibility shim for external callers).
# All symbols are sourced from the scripts/s05/ package. The try/except
# pattern allows ``python scripts/s05_insert_assembly.py --help`` to work
# even when the package is not on sys.path (standalone invocation).
# ---------------------------------------------------------------------------
try:
    from s05.verdict import VerdictRules, FilterEvidence, compute_verdict, _apply_canonical_override  # noqa: F401
    from s05.config_loader import load_verdict_rules  # noqa: F401
    from s05.primitives import (  # noqa: F401
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
    from s05.filter_b_flank import (  # noqa: F401
        CONSTRUCT_FLANK_PIDENT,
        CONSTRUCT_FLANK_MIN_LEN,
        CONSTRUCT_FLANK_SLOP,
        _find_construct_flanking_regions,
        _site_overlaps_flanking,
    )
    from s05.filter_c_chimeric import (  # noqa: F401
        CHIMERIC_MIN_PIDENT,
        CHIMERIC_MIN_OFFTARGET_BP,
        _check_chimeric_assembly,
    )
    from s05.filter_a_host import (  # noqa: F401
        INSERT_HOST_MIN_PIDENT,
        _blast_insert_vs_host,
    )
    from s05.filter_d_altlocus import (  # noqa: F401
        CONSTRUCT_HOST_MIN_COMBINED,
        CONSTRUCT_MIN_FRACTION,
        CONSTRUCT_HOST_MIN_PIDENT,
        _check_construct_host_coverage,
    )
    from s05.annotation import (  # noqa: F401
        _parse_blast6,
        _run_local_blast,
        _run_remote_blast,
        _merge_annotations,
        annotate_insert,
    )
    from s05.assembly import (  # noqa: F401
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
    from s05.classify import (  # noqa: F401
        _batch_check_element_hits,
        _filter_host_endogenous,
        _should_replace,
        _SRC_TIER,
        _TIER2_SRCS,
        _UNKNOWN_SRC_WARNED,
        classify_site_tiers,
        write_tier_classification,
    )
    from s05.site_discovery import (  # noqa: F401
        _build_consensus,
        _batch_check_maps_to_host,
        find_softclip_junctions,
        _apply_mask_bed,
        MASKED_SOURCE_TAG,
        parse_legacy_junctions,
        _extract_seeds_at_positions,
        legacy_junctions_to_sites,
    )
    from s05.read_extraction import (  # noqa: F401
        extract_candidate_reads,
        extract_unmapped_paired,
    )
    from s05.report import (  # noqa: F401
        generate_report,
        write_stats,
    )
    from s05.fanout_orchestrator import main as _fanout_main  # noqa: F401
except ImportError:  # pragma: no cover - allow standalone `python scripts/s05_insert_assembly.py`
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from s05.verdict import VerdictRules, FilterEvidence, compute_verdict, _apply_canonical_override  # noqa: F401
    from s05.config_loader import load_verdict_rules  # noqa: F401
    from s05.primitives import (  # noqa: F401
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
    from s05.filter_b_flank import (  # noqa: F401
        CONSTRUCT_FLANK_PIDENT,
        CONSTRUCT_FLANK_MIN_LEN,
        CONSTRUCT_FLANK_SLOP,
        _find_construct_flanking_regions,
        _site_overlaps_flanking,
    )
    from s05.filter_c_chimeric import (  # noqa: F401
        CHIMERIC_MIN_PIDENT,
        CHIMERIC_MIN_OFFTARGET_BP,
        _check_chimeric_assembly,
    )
    from s05.filter_a_host import (  # noqa: F401
        INSERT_HOST_MIN_PIDENT,
        _blast_insert_vs_host,
    )
    from s05.filter_d_altlocus import (  # noqa: F401
        CONSTRUCT_HOST_MIN_COMBINED,
        CONSTRUCT_MIN_FRACTION,
        CONSTRUCT_HOST_MIN_PIDENT,
        _check_construct_host_coverage,
    )
    from s05.annotation import (  # noqa: F401
        _parse_blast6,
        _run_local_blast,
        _run_remote_blast,
        _merge_annotations,
        annotate_insert,
    )
    from s05.assembly import (  # noqa: F401
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
    from s05.classify import (  # noqa: F401
        _batch_check_element_hits,
        _filter_host_endogenous,
        _should_replace,
        _SRC_TIER,
        _TIER2_SRCS,
        _UNKNOWN_SRC_WARNED,
        classify_site_tiers,
        write_tier_classification,
    )
    from s05.site_discovery import (  # noqa: F401
        _build_consensus,
        _batch_check_maps_to_host,
        find_softclip_junctions,
        _apply_mask_bed,
        MASKED_SOURCE_TAG,
        parse_legacy_junctions,
        _extract_seeds_at_positions,
        legacy_junctions_to_sites,
    )
    from s05.read_extraction import (  # noqa: F401
        extract_candidate_reads,
        extract_unmapped_paired,
    )
    from s05.report import (  # noqa: F401
        generate_report,
        write_stats,
    )
    from s05.fanout_orchestrator import main as _fanout_main  # noqa: F401

# ---------------------------------------------------------------------------
# Module-level constants (kept here for backward compat; originals live in
# scripts/s05/classify.py and scripts/s05/filter_*.py).
# ---------------------------------------------------------------------------
# Host-endogenous exclusion thresholds (used in classify_site_tiers).
# Two tiers handle cultivar drift — P-Act1 hits Nipponbare at ~78%/39%
# due to Xiushui vs Nipponbare divergence.
HOST_ENDO_T1_PIDENT = 90.0   # Tier 1 (clean host): min % identity
HOST_ENDO_T1_QCOVS = 50.0    # Tier 1: min query coverage %
HOST_ENDO_T2_PIDENT = 75.0   # Tier 2 (divergent host): min % identity
HOST_ENDO_T2_QCOVS = 30.0    # Tier 2: min query coverage %
CLIP_HOST_PIDENT = 95.0       # Per-clip host check: min % identity
CLIP_HOST_MIN_LEN = 30        # Per-clip host check: min alignment length
# Post-assembly host-fraction filter thresholds.
INSERT_HOST_FRACTION = 0.80     # host coverage threshold
INSERT_MIN_FOREIGN_GAP = 500    # non-host gap must be ≥ this to be real T-DNA
# UNKNOWN → FALSE_POSITIVE auto-reclassification thresholds.
UNKNOWN_HOST_MIN_FRACTION = 0.85   # min host fraction to classify as host-only
UNKNOWN_MAX_CONSTRUCT_FRAC = 0.05  # max construct fraction for host-only

# Phase 1.5 (_TIER2_SRCS) is derived from _SRC_TIER to stay in sync with
# gmo_combined_db_v2.fa cd-hit clustering; see scripts/s05/classify.py for
# the full audit trail comment.

# T8 fan-out state files (positive_sites.json + positive_sites.pkl) are
# written/read by fanout_orchestrator._run_phase_1_1_5.  The pickle is
# consumed only from the same sample's step_dir (trusted path, same user).
# v1.1 module split will replace pickle with typed JSON + dataclasses loader.

# ---------------------------------------------------------------------------
# Main entrypoint (orchestration moved to scripts/s05/fanout_orchestrator.py)
# ---------------------------------------------------------------------------

def main() -> None:
    """Thin entrypoint — full orchestration in scripts/s05/fanout_orchestrator."""
    from s05.fanout_orchestrator import main as _orchestrator_main
    _orchestrator_main()


if __name__ == "__main__":
    main()
