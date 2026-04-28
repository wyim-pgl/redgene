"""Phase 2-3 per-site re-exports (shim, v1.0).

Session 4 (2026-04-28): ``assemble_insert`` and ``refine_with_foreign_reads``
now live in ``scripts/s05/assembly.py``; this shim sources them from there.
Session 5 (2026-04-28): all 4 re-exports now come from native modules
(.assembly, .read_extraction). The shim itself can be retired in Session 7.
"""
from .assembly import (  # noqa: F401
    assemble_insert,
    refine_with_foreign_reads,
)
from .read_extraction import (  # noqa: F401
    extract_candidate_reads,
    extract_unmapped_paired,
)
