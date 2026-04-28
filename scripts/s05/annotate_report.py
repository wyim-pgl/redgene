"""Phase 4 annotate + report re-exports (shim, v1.0).

Session 3 (2026-04-28): ``annotate_insert`` now lives in
``scripts/s05/annotation.py``; this shim sources it from there while
``generate_report`` and ``write_stats`` continue to ride on the monolith
until Session 6 carves out ``report.py``.
"""
from .annotation import annotate_insert  # noqa: F401
from scripts.s05_insert_assembly import (  # noqa: F401
    generate_report,
    write_stats,
)
