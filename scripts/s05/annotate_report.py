"""Phase 4-5 annotate + report re-exports (shim, v1.0).

Session 6 (2026-04-29): generate_report and write_stats now live in
``scripts/s05/report.py``; this shim sources them from there. Session
3 already migrated annotate_insert to ``scripts/s05/annotation.py``.
The shim itself can be retired in Session 7.
"""
from .annotation import annotate_insert  # noqa: F401
from .report import generate_report, write_stats  # noqa: F401
