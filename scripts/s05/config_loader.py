"""Load VerdictRules from config.yaml.

T6 (v1.0 MVP) — parses the ``verdict_rules`` + ``canonical_triplets`` global
sections and returns a VerdictRules dataclass. Defaults mirror the compiled-in
constants in s05_insert_assembly.py so missing config keys do not change
behaviour.
"""
from __future__ import annotations

from dataclasses import fields
from pathlib import Path
from typing import Any

import yaml

from .verdict import VerdictRules


# Canonical triplets shipped as fallback when config.yaml omits the section.
# Aligns with docs/team-review/work_implementation_plan.md Task 6 Step 7.
DEFAULT_TRIPLETS: dict[str, set[str]] = {
    "default": {"bar", "P-CaMV35S", "T-ocs"},
    "rice_G281": {"hLF1", "P-Gt1", "T-nos"},
    "soybean_AtYUCCA6": {"bar", "P-CaMV35S", "T-ocs"},
    "soybean_UGT72E3": {"bar", "P-CaMV35S", "T-nos"},
    "tomato_Cas9_A2_3": {"bar", "SpCas9", "sgRNA_scaffold_generic"},
}


def _coerce_triplets(raw: Any) -> dict[str, set[str]]:
    """Convert YAML-decoded triplet config (dict[str, list[str]]) into sets."""
    if not isinstance(raw, dict):
        return {}
    out: dict[str, set[str]] = {}
    for key, val in raw.items():
        if val is None:
            continue
        if isinstance(val, (list, tuple, set)):
            out[str(key)] = {str(x) for x in val}
    return out


def load_verdict_rules(config_path: Path, sample: str) -> VerdictRules:
    """Parse config.yaml and build a VerdictRules for the given sample.

    Resolution order:
      * threshold fields: ``verdict_rules:`` block at the top level overrides
        compiled-in defaults; missing keys stay at their default.
      * canonical_triplets: merge DEFAULT_TRIPLETS with any top-level
        ``canonical_triplets:`` block. Per-sample keys (e.g. ``rice_G281``)
        remain accessible on the returned mapping so callers can opt in via
        ``rules.canonical_triplets[sample]``.

    The ``sample`` argument is accepted for future per-sample threshold
    overrides (v1.1+); today it is used only to help the caller select the
    matching triplet entry.
    """
    path = Path(config_path)
    if not path.exists():
        return VerdictRules(canonical_triplets=dict(DEFAULT_TRIPLETS))

    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}

    # Threshold overrides.
    threshold_fields = {
        f.name for f in fields(VerdictRules) if f.name != "canonical_triplets"
    }
    thresholds_cfg = raw.get("verdict_rules") or {}
    if not isinstance(thresholds_cfg, dict):
        thresholds_cfg = {}
    overrides: dict[str, Any] = {
        k: v for k, v in thresholds_cfg.items() if k in threshold_fields
    }

    # Canonical-triplet merging: start from DEFAULT_TRIPLETS so missing keys
    # still give callers a sane fallback, then overlay anything explicit.
    triplets: dict[str, set[str]] = {k: set(v) for k, v in DEFAULT_TRIPLETS.items()}
    triplets.update(_coerce_triplets(raw.get("canonical_triplets")))

    # `sample` currently does not select threshold overrides but is reserved
    # for v1.1 per-sample `samples: { foo: { verdict_rules: {...} } }` blocks.
    del sample  # suppress unused-argument lint; kept for forward compatibility

    return VerdictRules(canonical_triplets=triplets, **overrides)
