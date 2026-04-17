"""Phase 1 site discovery re-exports (shim, v1.0).

Full decomposition to independent module scheduled for v1.1. For now,
this shim exposes the monolithic script's public entry points at the
scripts.s05 import path so T8 per-site fan-out can use stable imports.
"""
from scripts.s05_insert_assembly import (  # noqa: F401
    find_softclip_junctions,
    _build_consensus,
    _batch_check_maps_to_host,
    _apply_mask_bed,    # T10
)
