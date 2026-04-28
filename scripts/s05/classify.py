"""Phase 1.5 classify re-exports (shim, v1.0).

Session 3 (2026-04-28): added ``_parse_src_tag`` to the re-export list so
external callers can discover it via the ``classify`` surface. The function
itself was moved to ``scripts/s05/primitives.py`` (stage 0) to avoid a
runtime circular import: annotation.py → classify.py → monolith → annotation
(partially initialised). Full classify extraction is scheduled for Session 4
per ``docs/superpowers/specs/2026-04-19-s05-module-split-design.md`` §4.
"""
from .primitives import _parse_src_tag  # noqa: F401
from scripts.s05_insert_assembly import (  # noqa: F401
    classify_site_tiers,
    _batch_check_element_hits,
    _filter_host_endogenous,
    _should_replace,
)
