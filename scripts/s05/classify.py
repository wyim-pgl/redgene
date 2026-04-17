"""Phase 1.5 classify re-exports (shim, v1.0)."""
from scripts.s05_insert_assembly import (  # noqa: F401
    classify_site_tiers,
    _batch_check_element_hits,
    _should_replace,
    _filter_host_endogenous,
)
