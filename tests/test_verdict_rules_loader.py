"""Tests for scripts.s05.config_loader.load_verdict_rules.

T6 pytest scenarios — 3 cases covering default fallback, global override, and
per-sample override.
"""
from pathlib import Path

from scripts.s05.config_loader import load_verdict_rules


def test_loader_default_when_key_missing(tmp_path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text("samples:\n  foo:\n    host_reference: /x\n")
    rules = load_verdict_rules(cfg, sample="foo")
    assert rules.cand_host_fraction_max == 0.80  # default
    assert "bar" in rules.canonical_triplets["default"]


def test_loader_reads_global_verdict_rules(tmp_path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "verdict_rules:\n"
        "  cand_host_fraction_max: 0.70\n"
        "canonical_triplets:\n"
        "  default: [bar, P-CaMV35S, T-ocs]\n"
        "samples:\n  foo: {}\n"
    )
    rules = load_verdict_rules(cfg, sample="foo")
    assert rules.cand_host_fraction_max == 0.70


def test_loader_sample_override(tmp_path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "canonical_triplets:\n"
        "  default: [bar, P-CaMV35S, T-ocs]\n"
        "  rice_G281: [hLF1, P-Gt1, T-nos]\n"
        "samples:\n  rice_G281: {}\n"
    )
    rules = load_verdict_rules(cfg, sample="rice_G281")
    assert rules.canonical_triplets["rice_G281"] == {"hLF1", "P-Gt1", "T-nos"}
