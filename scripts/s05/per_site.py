"""Phase 2-3 per-site re-exports (shim, v1.0).

Session 4 (2026-04-28): ``assemble_insert`` and ``refine_with_foreign_reads``
now live in ``scripts/s05/assembly.py``; this shim sources them from there
while ``extract_candidate_reads`` / ``extract_unmapped_paired`` continue
to ride on the monolith pending Session 5 (read_extraction.py).
"""
from .assembly import (  # noqa: F401
    assemble_insert,
    refine_with_foreign_reads,
)
from scripts.s05_insert_assembly import (  # noqa: F401
    extract_candidate_reads,
    extract_unmapped_paired,
)
