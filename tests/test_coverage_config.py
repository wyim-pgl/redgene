"""Tests for Issue #2 [AC-7] coverage-sensitivity config.yaml entries.

Verifies that the 12 (3 hosts x 4 coverage) sample keys exist in config.yaml
with correct structure so run_pipeline.py can resolve them when invoked by
run_coverage_sensitivity.sh.
"""
from __future__ import annotations

from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")


REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG = REPO_ROOT / "config.yaml"

HOSTS = ("rice_G281", "tomato_Cas9_A2_3", "cucumber_line225")
COVS = ("5x", "10x", "15x", "20x")


def _load_config() -> dict:
    assert CONFIG.exists(), CONFIG
    with CONFIG.open() as fh:
        return yaml.safe_load(fh)


def test_config_yaml_parses():
    cfg = _load_config()
    assert "samples" in cfg, "config.yaml missing top-level 'samples' key"


def test_twelve_coverage_entries_exist():
    cfg = _load_config()
    samples = cfg["samples"]
    missing = []
    for host in HOSTS:
        for cov in COVS:
            key = f"{host}_cov{cov}"
            if key not in samples:
                missing.append(key)
    assert not missing, f"config.yaml missing coverage entries: {missing}"


def test_coverage_entries_point_at_results_subsample_dir():
    """Each entry's FASTQ paths must live under results/<sample>_cov<tag>/.

    That is where run_coverage_sensitivity.sh writes the seqkit output, so
    run_pipeline.py will pick them up automatically once the subsample step
    completes.
    """
    cfg = _load_config()
    samples = cfg["samples"]
    for host in HOSTS:
        for cov in COVS:
            key = f"{host}_cov{cov}"
            entry = samples[key]
            assert "reads" in entry, f"{key} missing 'reads'"
            r1 = entry["reads"]["r1"]
            r2 = entry["reads"]["r2"]
            expected_prefix = f"results/{key}/"
            assert r1.startswith(expected_prefix), (
                f"{key}.reads.r1 = {r1!r} (want prefix {expected_prefix!r})"
            )
            assert r2.startswith(expected_prefix), (
                f"{key}.reads.r2 = {r2!r} (want prefix {expected_prefix!r})"
            )


def test_coverage_entries_inherit_host_reference():
    """Sanity check: each entry inherits its parent's host_reference."""
    cfg = _load_config()
    samples = cfg["samples"]
    for host in HOSTS:
        parent = samples[host]
        for cov in COVS:
            key = f"{host}_cov{cov}"
            child = samples[key]
            assert child["host_reference"] == parent["host_reference"], (
                f"{key} host_reference mismatches parent {host}"
            )
            assert child["construct_reference"] == parent["construct_reference"], (
                f"{key} construct_reference mismatches parent {host}"
            )
