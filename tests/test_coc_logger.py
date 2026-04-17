"""Tests for Issue #8 — scripts/util/coc_logger.py scaffold.

Chain-of-custody logger. Every pipeline step should eventually wrap itself
in a ``CocLogger`` context manager so the resulting
``results/<sample>/chain_of_custody.jsonl`` is an append-only audit trail.

Scope of THIS scaffold: the logger itself (context manager + helpers). The
wire-in into ``run_pipeline.py`` is intentionally deferred to a separate
v1.1 PR (see docs/architecture/chain_of_custody.md).
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "util" / "coc_logger.py"


def _load_module():
    assert SCRIPT.exists(), SCRIPT
    spec = importlib.util.spec_from_file_location("coc_logger", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["coc_logger"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_script_exists():
    assert SCRIPT.exists(), SCRIPT


def test_context_manager_writes_start_and_end(tmp_path):
    """Entering/exiting the logger emits a `start` and `end` event."""
    mod = _load_module()
    jsonl = tmp_path / "chain_of_custody.jsonl"
    with mod.CocLogger(path=jsonl, sample="S1", step="s03") as logger:
        assert logger is not None

    lines = jsonl.read_text().splitlines()
    assert len(lines) == 2
    events = [json.loads(line) for line in lines]
    assert events[0]["event"] == "start"
    assert events[1]["event"] == "end"
    assert all(e["sample"] == "S1" and e["step"] == "s03" for e in events)
    # Wall-time present on `end`
    assert isinstance(events[1].get("wall_time_sec"), (int, float))


def test_append_only_does_not_truncate(tmp_path):
    """Second context manager run must APPEND (not overwrite) prior events."""
    mod = _load_module()
    jsonl = tmp_path / "chain_of_custody.jsonl"
    with mod.CocLogger(path=jsonl, sample="S1", step="s03"):
        pass
    with mod.CocLogger(path=jsonl, sample="S1", step="s04"):
        pass

    lines = jsonl.read_text().splitlines()
    assert len(lines) == 4, lines
    assert [json.loads(l)["step"] for l in lines] == ["s03", "s03", "s04", "s04"]


def test_log_custom_event(tmp_path):
    """`logger.log(event, payload)` writes a custom-event line between start/end."""
    mod = _load_module()
    jsonl = tmp_path / "chain_of_custody.jsonl"
    with mod.CocLogger(path=jsonl, sample="S1", step="s03") as logger:
        logger.log("progress", {"pct": 42, "note": "halfway"})

    events = [json.loads(l) for l in jsonl.read_text().splitlines()]
    assert events[0]["event"] == "start"
    assert events[1]["event"] == "progress"
    assert events[1]["pct"] == 42
    assert events[2]["event"] == "end"


def test_error_event_captured_on_exception(tmp_path):
    """If the context manager body raises, an `error` event must be emitted."""
    mod = _load_module()
    jsonl = tmp_path / "chain_of_custody.jsonl"

    import pytest
    with pytest.raises(RuntimeError):
        with mod.CocLogger(path=jsonl, sample="S1", step="s04"):
            raise RuntimeError("boom")

    events = [json.loads(l) for l in jsonl.read_text().splitlines()]
    assert events[0]["event"] == "start"
    assert events[-1]["event"] == "error"
    # exit_code != 0 on error
    assert events[-1].get("exit_code") != 0
    # Error message preserved
    assert "boom" in events[-1].get("cmd", "") or "boom" in events[-1].get("error", "")


def test_sha256_helper(tmp_path):
    """`sha256_file()` returns the hex digest of a file's contents."""
    mod = _load_module()
    f = tmp_path / "sample.txt"
    f.write_bytes(b"hello world")
    expected = hashlib.sha256(b"hello world").hexdigest()
    assert mod.sha256_file(f) == expected


def test_log_with_input_and_output_sha(tmp_path):
    """Inputs/outputs can be logged with their SHA256 digests."""
    mod = _load_module()
    jsonl = tmp_path / "chain_of_custody.jsonl"
    inp = tmp_path / "in.bin"
    inp.write_bytes(b"x" * 100)
    outp = tmp_path / "out.bin"
    outp.write_bytes(b"y" * 200)

    with mod.CocLogger(path=jsonl, sample="S1", step="s05") as logger:
        logger.log("checksum", {
            "input_sha256": mod.sha256_file(inp),
            "output_sha256": mod.sha256_file(outp),
        })

    events = [json.loads(l) for l in jsonl.read_text().splitlines()]
    middle = events[1]
    assert middle["input_sha256"] == hashlib.sha256(b"x" * 100).hexdigest()
    assert middle["output_sha256"] == hashlib.sha256(b"y" * 200).hexdigest()


def test_iso8601_timestamp_format(tmp_path):
    """`ts` field must be ISO-8601 (`fromisoformat` round-trips)."""
    mod = _load_module()
    from datetime import datetime
    jsonl = tmp_path / "chain_of_custody.jsonl"
    with mod.CocLogger(path=jsonl, sample="S1", step="s03"):
        pass
    events = [json.loads(l) for l in jsonl.read_text().splitlines()]
    for e in events:
        # datetime.fromisoformat parses standard ISO-8601 with timezone
        datetime.fromisoformat(e["ts"])


def test_sidecar_user_field_populated(tmp_path):
    """`user` field should be populated from $USER or getpass (for provenance)."""
    mod = _load_module()
    jsonl = tmp_path / "chain_of_custody.jsonl"
    with mod.CocLogger(path=jsonl, sample="S1", step="s03"):
        pass
    events = [json.loads(l) for l in jsonl.read_text().splitlines()]
    for e in events:
        assert e.get("user"), "user field must be populated for provenance"


def test_docs_chain_of_custody_exists():
    """Schema + usage doc is shipped alongside the scaffold."""
    doc = REPO_ROOT / "docs" / "architecture" / "chain_of_custody.md"
    assert doc.exists(), doc
    text = doc.read_text()
    # Key schema fields must be documented
    for key in ("input_sha256", "output_sha256", "wall_time_sec",
                "exit_code", "coc_chain_id"):
        assert key in text, f"chain_of_custody.md missing `{key}`"
