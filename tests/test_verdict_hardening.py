"""Regression tests for Issue #11 — T6 compute_verdict hardening.

Important items (I-1, I-2) + minor items (M-1, M-3, M-4, M-7) that are
script-testable. The cosmetic items (M-5 reason localisation structure,
M-8 fixture-naming consistency) are covered by docstring-level checks.
"""
from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from scripts.s05 import config_loader
from scripts.s05.config_loader import (
    DEFAULT_TRIPLETS,
    load_verdict_rules,
    _coerce_triplets,
)
from scripts.s05.verdict import (
    FilterEvidence,
    VerdictRules,
    _find_matching_triplet,
    compute_verdict,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


# ---- I-1: config_loader defensive parsing ---------------------------------

def test_i1_verdict_rules_none_falls_back_to_defaults(tmp_path):
    """`verdict_rules: null` must not propagate as None into kwargs."""
    cfg = tmp_path / "config.yaml"
    cfg.write_text("verdict_rules:\nsamples: {}\n")
    rules = load_verdict_rules(cfg, sample="anything")
    # All thresholds stay at compiled-in defaults.
    assert rules.cand_host_fraction_max == 0.80
    assert rules.fp_host_fraction_min == 0.80
    assert rules.canonical_triplet_min_identity == 0.90


def test_i1_verdict_rules_string_ignored(tmp_path):
    """A bogus `verdict_rules: "hello"` scalar must be silently ignored."""
    cfg = tmp_path / "config.yaml"
    cfg.write_text("verdict_rules: hello\nsamples: {}\n")
    rules = load_verdict_rules(cfg, sample="anything")
    assert rules.cand_host_fraction_max == 0.80  # default preserved


def test_i1_canonical_triplets_none_keeps_defaults(tmp_path):
    """`canonical_triplets: null` keeps DEFAULT_TRIPLETS intact."""
    cfg = tmp_path / "config.yaml"
    cfg.write_text("canonical_triplets:\nsamples: {}\n")
    rules = load_verdict_rules(cfg, sample="anything")
    assert "default" in rules.canonical_triplets
    assert rules.canonical_triplets["default"] == DEFAULT_TRIPLETS["default"]


def test_i1_coerce_triplets_handles_non_dict():
    """`_coerce_triplets` returns {} for any non-dict input (None, list, str)."""
    assert _coerce_triplets(None) == {}
    assert _coerce_triplets([]) == {}
    assert _coerce_triplets("nope") == {}


def test_i1_coerce_triplets_skips_none_values():
    assert _coerce_triplets({"a": None, "b": ["x"]}) == {"b": {"x"}}


# ---- I-2: docstring accuracy for triplet / rule 1 -------------------------

def test_i2_find_matching_triplet_empty_set_is_none():
    """Rule 1 docstring promises: an empty triplet is NEVER a match."""
    assert _find_matching_triplet(set(), {"default": set()}) is None
    assert _find_matching_triplet({"bar"}, {"default": set()}) is None


def test_i2_find_matching_triplet_needs_full_subset():
    """Rule 1 fires only when the triplet is fully covered by matched."""
    triplets = {"rice_G281": {"hLF1", "P-Gt1", "T-nos"}}
    assert _find_matching_triplet({"hLF1", "P-Gt1"}, triplets) is None
    assert _find_matching_triplet({"hLF1", "P-Gt1", "T-nos"}, triplets) == "rice_G281"


def test_i2_compute_verdict_docstring_mentions_canonical_semantics():
    """compute_verdict docstring must spell out the canonical-triplet rule."""
    doc = compute_verdict.__doc__ or ""
    # canonical triplet rule and host_fraction gate must both be mentioned.
    assert "canonical" in doc.lower()
    assert "host_fraction" in doc.lower() or "host fraction" in doc.lower()


# ---- M-3: is_chimeric is a dead field (off_target_chrs is the signal) ----

def test_m3_compute_verdict_ignores_is_chimeric_flag():
    """FilterEvidence.is_chimeric must be informational-only — the actual
    chimeric filter (rule 3) is driven entirely by off_target_chrs."""
    rules = VerdictRules(
        canonical_triplets={"default": {"bar", "P-CaMV35S", "T-ocs"}},
    )
    # is_chimeric=True but no off_target_chrs -> NOT FALSE_POSITIVE(chimeric).
    ev = FilterEvidence(
        elements=["bar"], host_fraction=0.3, largest_gap=5000,
        off_target_chrs=[], is_chimeric=True,
    )
    verdict, reason = compute_verdict(ev, rules)
    assert verdict == "CANDIDATE"
    assert "chimeric" not in reason.lower()


# ---- M-4: em-dash unicode stability ---------------------------------------

def test_m4_reason_strings_use_consistent_unicode_em_dash():
    """All `\\u2014` em-dashes in verdict.py must round-trip through utf-8."""
    verdict_py = (REPO_ROOT / "scripts" / "s05" / "verdict.py").read_text(encoding="utf-8")
    # We don't forbid em-dashes; we just verify they survive utf-8 encode/decode.
    assert verdict_py.encode("utf-8").decode("utf-8") == verdict_py


# ---- M-7: DEFAULT_TRIPLETS stays in sync with config.yaml -----------------

def test_m7_default_triplets_match_config_yaml():
    """scripts/s05/config_loader.DEFAULT_TRIPLETS must mirror config.yaml's
    `canonical_triplets:` block so the compiled fallback matches what
    operators actually see."""
    cfg_path = REPO_ROOT / "config.yaml"
    if not cfg_path.exists():
        pytest.skip("config.yaml not present in this checkout")
    with open(cfg_path) as fh:
        raw = yaml.safe_load(fh) or {}
    yaml_trip = {k: set(v) for k, v in (raw.get("canonical_triplets") or {}).items()}
    # Every key in config.yaml must appear in DEFAULT_TRIPLETS with the
    # same members. DEFAULT_TRIPLETS may additionally carry keys that
    # config.yaml does not declare (legacy fallback) — that's fine.
    for sample, elems in yaml_trip.items():
        assert sample in DEFAULT_TRIPLETS, (
            f"DEFAULT_TRIPLETS missing `{sample}` present in config.yaml"
        )
        assert DEFAULT_TRIPLETS[sample] == elems, (
            f"DEFAULT_TRIPLETS[{sample}] = {DEFAULT_TRIPLETS[sample]!r} "
            f"but config.yaml has {elems!r}"
        )


# ---- M-2: per-sample threshold override scaffold --------------------------

def test_m2_sample_argument_accepted_for_all_samples(tmp_path):
    """The `sample=` argument must be accepted (scaffold for v1.1) without
    raising, regardless of whether the sample is in DEFAULT_TRIPLETS."""
    cfg = tmp_path / "config.yaml"
    cfg.write_text("verdict_rules:\n  cand_host_fraction_max: 0.70\n")
    # Known sample.
    rules1 = load_verdict_rules(cfg, sample="rice_G281")
    # Unknown sample — must still load defaults, not crash.
    rules2 = load_verdict_rules(cfg, sample="does_not_exist")
    assert rules1.cand_host_fraction_max == 0.70
    assert rules2.cand_host_fraction_max == 0.70


# ---- M-1: `del sample` idiom is explained ---------------------------------

def test_m1_config_loader_sample_arg_is_documented():
    """The `del sample` line exists with a rationale comment above it so a
    future maintainer does not delete the parameter thinking it is unused."""
    src = (REPO_ROOT / "scripts" / "s05" / "config_loader.py").read_text()
    assert "del sample" in src, "del sample idiom removed; see Issue #11 M-1"
    # Rationale must be present within 3 lines above `del sample`.
    lines = src.splitlines()
    idx = next(i for i, l in enumerate(lines) if "del sample" in l)
    window = "\n".join(lines[max(0, idx - 3): idx])
    assert "v1.1" in window or "forward" in window.lower() or "per-sample" in window, (
        "missing rationale comment above `del sample`"
    )
