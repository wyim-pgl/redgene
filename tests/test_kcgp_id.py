"""Tests for scripts/util/kcgp_id.py — KCGP tentative_v0 ID format."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "util" / "kcgp_id.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("kcgp_id", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["kcgp_id"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def test_build_basic():
    mod = _load_module()
    out = mod.build_kcgp_id("db/Osativa_323_v7.0.fa", "G281", 2026, 1)
    assert out == "KCGP-OSA-G281-26-0001"


def test_build_all_5_hosts():
    mod = _load_module()
    cases = [
        ("db/Osativa_323_v7.0.fa", "OSA"),
        ("db/SLM_r2.0.pmol.fasta", "SLE"),
        ("db/CucSat_B10v3.fa",     "CSA"),
        ("db/Zm_B73_v5.fa",        "ZMA"),
        ("db/Gmax_v4.0.fa",        "GMA"),
    ]
    for host_path, expected_code in cases:
        out = mod.build_kcgp_id(host_path, "X", 2026, 1)
        assert f"KCGP-{expected_code}-X-26-0001" == out


def test_build_unknown_host_raises():
    mod = _load_module()
    with pytest.raises(KeyError) as exc:
        mod.build_kcgp_id("db/unknown_host.fa", "G281", 2026, 1)
    msg = str(exc.value)
    assert "unknown_host" in msg
    assert "OSA" in msg  # known host list shown for debuggability


def test_round_trip():
    mod = _load_module()
    cases = [
        ("db/Osativa_323_v7.0.fa", "G281", 2026, 1),
        ("db/SLM_r2.0.pmol.fasta", "A2_3", 2026, 42),
        ("db/CucSat_B10v3.fa",     "line225", 2025, 9999),
    ]
    for host, event, year, site_ord in cases:
        kid = mod.build_kcgp_id(host, event, year, site_ord)
        parsed = mod.parse_kcgp_id(kid)
        assert parsed["host_code"] == mod.HOST_CODES[mod.normalize_host_path(host)]
        assert parsed["event"] == event
        assert parsed["year_2digit"] == year % 100
        assert parsed["site_ord"] == site_ord


def test_parse_invalid_format():
    mod = _load_module()
    with pytest.raises(ValueError):
        mod.parse_kcgp_id("KCGP-XX")
    with pytest.raises(ValueError):
        mod.parse_kcgp_id("not-a-kcgp-id")


def test_event_with_underscore():
    mod = _load_module()
    out = mod.build_kcgp_id("db/SLM_r2.0.pmol.fasta", "A2_3", 2026, 1)
    assert out == "KCGP-SLE-A2_3-26-0001"
    parsed = mod.parse_kcgp_id(out)
    assert parsed["event"] == "A2_3"


def test_normalize_host_path():
    mod = _load_module()
    assert mod.normalize_host_path("db/Osativa_323_v7.0.fa") == "Osativa_323_v7.0"
    assert mod.normalize_host_path("db/SLM_r2.0.pmol.fasta") == "SLM_r2.0.pmol"
    assert mod.normalize_host_path("Osativa_323_v7.0") == "Osativa_323_v7.0"


def test_site_overflow():
    mod = _load_module()
    mod.build_kcgp_id("db/Osativa_323_v7.0.fa", "X", 2026, 9999)
    with pytest.raises(ValueError):
        mod.build_kcgp_id("db/Osativa_323_v7.0.fa", "X", 2026, 10000)


def test_build_rejects_event_with_hyphen():
    mod = _load_module()
    with pytest.raises(ValueError) as exc:
        mod.build_kcgp_id("db/Osativa_323_v7.0.fa", "A2-3", 2026, 1)
    assert "event" in str(exc.value).lower()


def test_build_rejects_empty_event():
    mod = _load_module()
    with pytest.raises(ValueError):
        mod.build_kcgp_id("db/Osativa_323_v7.0.fa", "", 2026, 1)
