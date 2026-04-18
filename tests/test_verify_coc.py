"""Tests for Issue #8 — tools/verify_coc.py integrity checker.

Covers:
  (a) Monotone timestamp verification
  (b) Pre/post pairing check
  (c) Hash continuity check (output of step N == input of step N+1)
  (d) Overall exit code behaviour (0 = clean, 1 = error)
  (e) Missing file returns exit code 2
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone, timedelta
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import importlib.util

VERIFIER_PATH = REPO / "tools" / "verify_coc.py"


def _load_verifier():
    assert VERIFIER_PATH.exists(), VERIFIER_PATH
    spec = importlib.util.spec_from_file_location("verify_coc", VERIFIER_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("verify_coc", mod)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ts(offset_sec: float = 0.0) -> str:
    base = datetime(2026, 4, 18, 10, 0, 0, tzinfo=timezone.utc)
    return (base + timedelta(seconds=offset_sec)).isoformat(timespec="microseconds")


def _write_jsonl(path: Path, records: list[dict]) -> None:
    path.write_text("\n".join(json.dumps(r) for r in records) + "\n", encoding="utf-8")


def _make_step_records(
    step: str,
    ts_start: float = 0.0,
    input_sha: dict | None = None,
    output_sha: dict | None = None,
    exit_code: int = 0,
) -> list[dict]:
    """Build [start, pre, post, end] records for a step."""
    return [
        {"ts": _ts(ts_start),     "sample": "S1", "step": step, "event": "start", "user": "u"},
        {"ts": _ts(ts_start + 1), "sample": "S1", "step": step, "event": "pre",
         "cmd": f"python s{step}.py", "input_sha256": input_sha or {}, "user": "u"},
        {"ts": _ts(ts_start + 2), "sample": "S1", "step": step, "event": "post",
         "exit_code": exit_code, "output_sha256": output_sha or {}, "user": "u"},
        {"ts": _ts(ts_start + 3), "sample": "S1", "step": step, "event": "end",
         "exit_code": 0, "wall_time_sec": 3.0, "user": "u"},
    ]


# ---------------------------------------------------------------------------
# Tests for verifier module existence
# ---------------------------------------------------------------------------

def test_verifier_script_exists():
    assert VERIFIER_PATH.exists(), VERIFIER_PATH


# ---------------------------------------------------------------------------
# Tests for monotone timestamps
# ---------------------------------------------------------------------------

def test_verify_monotone_clean(tmp_path):
    """Strictly increasing timestamps pass verification."""
    mod = _load_verifier()
    records = _make_step_records("s03", ts_start=0.0)
    errors = mod.verify_monotone_timestamps(records)
    assert errors == []


def test_verify_monotone_detects_backward_jump(tmp_path):
    """A timestamp that goes backward must produce an error."""
    mod = _load_verifier()
    records = [
        {"ts": _ts(10.0), "step": "s03", "event": "start", "sample": "S1"},
        {"ts": _ts(5.0),  "step": "s03", "event": "end",   "sample": "S1"},
    ]
    errors = mod.verify_monotone_timestamps(records)
    assert len(errors) == 1
    assert "earlier" in errors[0].lower() or "timestamp" in errors[0].lower()


def test_verify_monotone_equal_timestamps_ok(tmp_path):
    """Equal timestamps are acceptable (same-millisecond writes)."""
    mod = _load_verifier()
    records = [
        {"ts": _ts(5.0), "step": "s03", "event": "start", "sample": "S1"},
        {"ts": _ts(5.0), "step": "s03", "event": "end",   "sample": "S1"},
    ]
    errors = mod.verify_monotone_timestamps(records)
    assert errors == []


# ---------------------------------------------------------------------------
# Tests for pre/post pairing
# ---------------------------------------------------------------------------

def test_verify_pre_post_pairs_clean(tmp_path):
    """A complete pre→post sequence for one step is valid."""
    mod = _load_verifier()
    records = _make_step_records("s03", ts_start=0.0)
    errors = mod.verify_pre_post_pairs(records)
    assert errors == []


def test_verify_pre_post_missing_post(tmp_path):
    """A 'pre' without a matching 'post' must be flagged."""
    mod = _load_verifier()
    records = [
        {"ts": _ts(0), "step": "s03", "event": "start", "sample": "S1"},
        {"ts": _ts(1), "step": "s03", "event": "pre",   "sample": "S1", "input_sha256": {}},
        # missing 'post'
        {"ts": _ts(2), "step": "s03", "event": "end",   "sample": "S1"},
    ]
    errors = mod.verify_pre_post_pairs(records)
    assert any("pre" in e.lower() for e in errors)


def test_verify_pre_post_missing_pre(tmp_path):
    """A 'post' without a preceding 'pre' must be flagged."""
    mod = _load_verifier()
    records = [
        {"ts": _ts(0), "step": "s04", "event": "start", "sample": "S1"},
        # missing 'pre'
        {"ts": _ts(1), "step": "s04", "event": "post",   "sample": "S1", "exit_code": 0, "output_sha256": {}},
        {"ts": _ts(2), "step": "s04", "event": "end",    "sample": "S1"},
    ]
    errors = mod.verify_pre_post_pairs(records)
    assert any("post" in e.lower() for e in errors)


def test_verify_pre_post_two_steps_clean(tmp_path):
    """Two consecutive steps each with proper pre/post pairs pass."""
    mod = _load_verifier()
    records = (
        _make_step_records("s03", ts_start=0.0)
        + _make_step_records("s04", ts_start=10.0)
    )
    errors = mod.verify_pre_post_pairs(records)
    assert errors == []


# ---------------------------------------------------------------------------
# Tests for hash continuity
# ---------------------------------------------------------------------------

def test_verify_hash_continuity_clean(tmp_path):
    """When step N's output SHA matches step N+1's input SHA, no error."""
    mod = _load_verifier()
    bam_sha = "a" * 64
    records = (
        _make_step_records("s03", ts_start=0.0, output_sha={"host.bam": bam_sha})
        + _make_step_records("s04", ts_start=10.0, input_sha={"host.bam": bam_sha})
    )
    errors = mod.verify_hash_continuity(records)
    assert errors == []


def test_verify_hash_continuity_mismatch_detected(tmp_path):
    """When step N's output SHA differs from step N+1's input SHA, error reported."""
    mod = _load_verifier()
    records = (
        _make_step_records("s03", ts_start=0.0, output_sha={"host.bam": "a" * 64})
        + _make_step_records("s04", ts_start=10.0, input_sha={"host.bam": "b" * 64})
    )
    errors = mod.verify_hash_continuity(records)
    assert len(errors) == 1
    assert "host.bam" in errors[0]


def test_verify_hash_continuity_different_filenames_ignored(tmp_path):
    """Files present in output but not in next input (or vice versa) are not errors."""
    mod = _load_verifier()
    records = (
        _make_step_records("s03", ts_start=0.0, output_sha={"file_A.bam": "a" * 64})
        + _make_step_records("s04", ts_start=10.0, input_sha={"file_B.bam": "b" * 64})
    )
    errors = mod.verify_hash_continuity(records)
    assert errors == []


# ---------------------------------------------------------------------------
# Integration: verify() end-to-end
# ---------------------------------------------------------------------------

def test_verify_returns_true_for_clean_jsonl(tmp_path):
    """verify() returns True when the JSONL passes all checks."""
    mod = _load_verifier()
    jsonl = tmp_path / "chain_of_custody.jsonl"
    sha = "c" * 64
    records = (
        _make_step_records("s03", ts_start=0.0, output_sha={"x.bam": sha})
        + _make_step_records("s04", ts_start=10.0, input_sha={"x.bam": sha})
    )
    _write_jsonl(jsonl, records)
    assert mod.verify(jsonl, verbose=False) is True


def test_verify_returns_false_for_broken_jsonl(tmp_path):
    """verify() returns False when any check fails."""
    mod = _load_verifier()
    jsonl = tmp_path / "chain_of_custody.jsonl"
    # Backward timestamp
    records = [
        {"ts": _ts(10.0), "step": "s03", "event": "start", "sample": "S1", "user": "u"},
        {"ts": _ts(5.0),  "step": "s03", "event": "end",   "sample": "S1", "user": "u"},
    ]
    _write_jsonl(jsonl, records)
    assert mod.verify(jsonl, verbose=False) is False


# ---------------------------------------------------------------------------
# CLI exit codes
# ---------------------------------------------------------------------------

def test_main_exits_0_for_clean(tmp_path):
    """CLI exits 0 for a clean JSONL."""
    mod = _load_verifier()
    jsonl = tmp_path / "chain_of_custody.jsonl"
    _write_jsonl(jsonl, _make_step_records("s03", ts_start=0.0))
    rc = mod.main([str(jsonl)])
    assert rc == 0


def test_main_exits_1_for_errors(tmp_path):
    """CLI exits 1 when verification fails."""
    mod = _load_verifier()
    jsonl = tmp_path / "chain_of_custody.jsonl"
    # Unpaired pre
    _write_jsonl(jsonl, [
        {"ts": _ts(0), "step": "s03", "event": "pre", "sample": "S1",
         "input_sha256": {}, "user": "u"},
    ])
    rc = mod.main([str(jsonl)])
    assert rc == 1


def test_main_exits_2_for_missing_file(tmp_path):
    """CLI exits 2 when the JSONL file does not exist."""
    mod = _load_verifier()
    missing = tmp_path / "no_such.jsonl"
    rc = mod.main([str(missing)])
    assert rc == 2
