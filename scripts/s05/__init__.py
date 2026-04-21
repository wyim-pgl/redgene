"""RedGene s05 package (Issue #4 multi-session split).

Pure-function modules extracted from monolithic ``s05_insert_assembly.py``.
Session 1 (2026-04-19) added ``primitives``, ``filter_b_flank``, and
``filter_c_chimeric`` alongside the already-extracted ``verdict`` and
``config_loader``. Remaining decomposition (site_discovery, classify,
read_extraction, assembly, annotation, report, fanout_orchestrator) is
tracked in ``docs/superpowers/specs/2026-04-19-s05-module-split-design.md``.
"""
from .primitives import (  # noqa: F401
    STEP,
    log,
    revcomp,
    read_fasta,
    write_fasta,
    _read_fq_seqs,
    JunctionCluster,
    InsertionSite,
    LegacyJunction,
    TierResult,
)
from .verdict import compute_verdict, FilterEvidence, VerdictRules  # noqa: F401
from .config_loader import load_verdict_rules  # noqa: F401
