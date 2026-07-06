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
    from s05.annotation import annotate_insert  # noqa: F401
    from s05.classify import (  # noqa: F401
        _batch_check_element_hits,
        _should_replace,
        _SRC_TIER,
        _TIER2_SRCS,
    )
except ImportError:  # pragma: no cover - allow standalone `python scripts/s05_insert_assembly.py`
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from s05.annotation import annotate_insert  # noqa: F401
    from s05.classify import (  # noqa: F401
        _batch_check_element_hits,
        _should_replace,
        _SRC_TIER,
        _TIER2_SRCS,
    )

# NOTE: the post-Issue-#4 authoritative thresholds live in scripts/s05/
# (classify.py HOST_ENDO_*/CLIP_HOST_*, filter_a_host.py INSERT_HOST_*,
# verdict.py / config_loader.py for the UNKNOWN-reclass values). A dead
# duplicate block that used to live here — with zero readers — was removed
# to keep a single source of truth.
#
# Audit trail: Phase 1.5 tier classification (_TIER2_SRCS derived from
# _SRC_TIER in scripts/s05/classify.py) is keyed to the cd-hit clustering of
# element_db/gmo_combined_db_v2.fa, whose |src=<tag> headers drive the 4-way
# source tiering. This reference links Phase 1.5 results back to that build
# artefact (guarded by tests/test_element_db_build.py::test_m1_*).

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
