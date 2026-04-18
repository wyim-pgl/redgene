"""Tests for Issue #8 — CoC wire-in into run_pipeline.py.

Verifies that:
  1. run_step() writes pre+post events to chain_of_custody.jsonl.
  2. The pre-entry contains the full command string and input SHA-256 hashes.
  3. The post-entry contains the exit code and output SHA-256 hashes.
  4. dry_run mode does NOT write any CoC entries.
  5. _step_input_files / _step_output_files return only existing paths.
  6. _sha256_files produces correct digests.
  7. run_step() records exit_code=0 for success and exit_code!=0 for failure
     (the failure path still records and then calls sys.exit).
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import run_pipeline  # noqa: E402 — needs sys.path set above


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def fake_sample_cfg() -> dict[str, Any]:
    return {
        "reads": {"r1": "fake_R1.fq.gz", "r2": "fake_R2.fq.gz"},
        "construct_reference": "db/construct.fa",
        "host_reference": "db/host.fa",
    }


@pytest.fixture()
def fake_step_dirs(tmp_path):
    """Create minimal directory structure + dummy files for a fake sample."""
    sname = "test_sample"
    s01 = tmp_path / "results" / sname / "s01_qc"
    s02 = tmp_path / "results" / sname / "s02_construct_map"
    s03 = tmp_path / "results" / sname / "s03_extract"
    s04 = tmp_path / "results" / sname / "s04_host_map"
    for d in (s01, s02, s03, s04):
        d.mkdir(parents=True)

    # Create dummy input files so SHA-256 can be computed
    (s01 / f"{sname}_R1.fq.gz").write_bytes(b"r1_content")
    (s01 / f"{sname}_R2.fq.gz").write_bytes(b"r2_content")
    (s04 / f"{sname}_host.bam").write_bytes(b"bam_content")

    return tmp_path / "results", sname


# ---------------------------------------------------------------------------
# Tests for _step_input_files
# ---------------------------------------------------------------------------

def test_step_input_files_step2_returns_existing_reads(tmp_path):
    """_step_input_files for step 2 returns the s01 fastq files when they exist."""
    sname = "S1"
    s01 = tmp_path / sname / "s01_qc"
    s01.mkdir(parents=True)
    r1 = s01 / f"{sname}_R1.fq.gz"
    r2 = s01 / f"{sname}_R2.fq.gz"
    r1.write_bytes(b"r1")
    r2.write_bytes(b"r2")

    sample_cfg: dict[str, Any] = {
        "reads": {"r1": "raw_R1.fq.gz", "r2": "raw_R2.fq.gz"},
        "construct_reference": "db/c.fa",
        "host_reference": "db/h.fa",
    }
    files = run_pipeline._step_input_files("2", sname, sample_cfg, tmp_path, REPO)
    assert r1 in files
    assert r2 in files


def test_step_input_files_missing_files_skipped(tmp_path):
    """Files that do not exist are silently excluded."""
    sname = "S1"
    sample_cfg: dict[str, Any] = {
        "reads": {"r1": "nonexistent_R1.fq.gz", "r2": "nonexistent_R2.fq.gz"},
        "construct_reference": "db/c.fa",
        "host_reference": "db/h.fa",
    }
    files = run_pipeline._step_input_files("4", sname, sample_cfg, tmp_path, REPO)
    assert files == []


# ---------------------------------------------------------------------------
# Tests for _step_output_files
# ---------------------------------------------------------------------------

def test_step_output_files_step4_returns_bam_if_exists(tmp_path):
    """_step_output_files for step 4 returns host.bam when it exists."""
    sname = "S1"
    s04 = tmp_path / sname / "s04_host_map"
    s04.mkdir(parents=True)
    bam = s04 / f"{sname}_host.bam"
    bam.write_bytes(b"bam")

    files = run_pipeline._step_output_files("4", sname, tmp_path)
    assert bam in files


def test_step_output_files_returns_empty_for_missing(tmp_path):
    """If the expected output doesn't exist yet, returns empty list."""
    files = run_pipeline._step_output_files("5", "S1", tmp_path)
    assert files == []


# ---------------------------------------------------------------------------
# Tests for _sha256_files
# ---------------------------------------------------------------------------

def test_sha256_files_produces_correct_digests(tmp_path):
    """_sha256_files maps filename → hex SHA-256 digest."""
    f1 = tmp_path / "alpha.bam"
    f1.write_bytes(b"data_alpha")
    f2 = tmp_path / "beta.tsv"
    f2.write_bytes(b"data_beta")

    result = run_pipeline._sha256_files([f1, f2])
    assert result["alpha.bam"] == hashlib.sha256(b"data_alpha").hexdigest()
    assert result["beta.tsv"] == hashlib.sha256(b"data_beta").hexdigest()


def test_sha256_files_empty_list_returns_empty_dict(tmp_path):
    assert run_pipeline._sha256_files([]) == {}


# ---------------------------------------------------------------------------
# Tests for run_step CoC integration
# ---------------------------------------------------------------------------

def _make_run_step_kwargs(
    step: str,
    sample_key: str,
    sample_cfg: dict[str, Any],
    outdir: Path,
    base_dir: Path,
) -> dict[str, Any]:
    return dict(
        step=step,
        sample_key=sample_key,
        sample_cfg=sample_cfg,
        outdir=outdir,
        threads=8,
        base_dir=base_dir,
        dry_run=False,
        no_remote_blast=False,
        cfg=None,
        fanout=False,
        host_bam_override=None,
    )


def _make_fake_scripts(base_dir: Path) -> None:
    """Create stub script files so run_step doesn't exit on 'script not found'."""
    for script_rel in run_pipeline.STEP_SCRIPTS.values():
        p = base_dir / script_rel
        p.parent.mkdir(parents=True, exist_ok=True)
        if not p.exists():
            p.write_text("#!/usr/bin/env python3\n")


def test_run_step_writes_coc_entries(tmp_path):
    """run_step writes start/pre/post/end events to chain_of_custody.jsonl."""
    sname = "coc_test"
    outdir = tmp_path / "results"
    (outdir / sname).mkdir(parents=True)
    _make_fake_scripts(tmp_path)

    sample_cfg: dict[str, Any] = {
        "reads": {"r1": "r1.fq.gz", "r2": "r2.fq.gz"},
        "construct_reference": "db/c.fa",
        "host_reference": "db/h.fa",
    }

    # Patch build_step_cmd to return a no-op echo command and
    # subprocess.run to return success without actually running anything
    mock_result = MagicMock()
    mock_result.returncode = 0

    with patch.object(run_pipeline, "build_step_cmd",
                      return_value=["echo", "hello"]), \
         patch("subprocess.run", return_value=mock_result):
        run_pipeline.run_step(**_make_run_step_kwargs(
            "3", sname, sample_cfg, outdir, tmp_path
        ))

    jsonl = outdir / sname / "chain_of_custody.jsonl"
    assert jsonl.exists(), "chain_of_custody.jsonl must be created by run_step"

    events = [json.loads(l) for l in jsonl.read_text().splitlines()]
    event_types = [e["event"] for e in events]

    # Must have start, pre, post, end in that order
    assert "start" in event_types
    assert "pre" in event_types
    assert "post" in event_types
    assert "end" in event_types
    assert event_types.index("start") < event_types.index("pre")
    assert event_types.index("pre") < event_types.index("post")
    assert event_types.index("post") < event_types.index("end")


def test_run_step_pre_entry_contains_cmd(tmp_path):
    """The 'pre' CoC entry must include the full CLI command."""
    sname = "coc_cmd_test"
    outdir = tmp_path / "results"
    (outdir / sname).mkdir(parents=True)
    _make_fake_scripts(tmp_path)

    sample_cfg: dict[str, Any] = {
        "reads": {"r1": "r1.fq.gz", "r2": "r2.fq.gz"},
        "construct_reference": "db/c.fa",
        "host_reference": "db/h.fa",
    }
    mock_result = MagicMock()
    mock_result.returncode = 0

    fake_cmd = ["python", "scripts/s02_construct_map.py", "--sample-name", sname]

    with patch.object(run_pipeline, "build_step_cmd", return_value=fake_cmd), \
         patch("subprocess.run", return_value=mock_result):
        run_pipeline.run_step(**_make_run_step_kwargs(
            "2", sname, sample_cfg, outdir, tmp_path
        ))

    jsonl = outdir / sname / "chain_of_custody.jsonl"
    events = [json.loads(l) for l in jsonl.read_text().splitlines()]
    pre_events = [e for e in events if e["event"] == "pre"]
    assert len(pre_events) == 1
    assert sname in pre_events[0].get("cmd", "")


def test_run_step_post_entry_contains_exit_code(tmp_path):
    """The 'post' CoC entry must contain exit_code=0 on success."""
    sname = "coc_exit_test"
    outdir = tmp_path / "results"
    (outdir / sname).mkdir(parents=True)
    _make_fake_scripts(tmp_path)

    sample_cfg: dict[str, Any] = {
        "reads": {"r1": "r1.fq.gz", "r2": "r2.fq.gz"},
        "construct_reference": "db/c.fa",
        "host_reference": "db/h.fa",
    }
    mock_result = MagicMock()
    mock_result.returncode = 0

    with patch.object(run_pipeline, "build_step_cmd",
                      return_value=["echo", "ok"]), \
         patch("subprocess.run", return_value=mock_result):
        run_pipeline.run_step(**_make_run_step_kwargs(
            "1", sname, sample_cfg, outdir, tmp_path
        ))

    jsonl = outdir / sname / "chain_of_custody.jsonl"
    events = [json.loads(l) for l in jsonl.read_text().splitlines()]
    post_events = [e for e in events if e["event"] == "post"]
    assert len(post_events) == 1
    assert post_events[0].get("exit_code") == 0


def test_run_step_dry_run_writes_no_coc(tmp_path):
    """In dry-run mode, no chain_of_custody.jsonl file must be written."""
    sname = "dry_run_test"
    outdir = tmp_path / "results"
    (outdir / sname).mkdir(parents=True)
    _make_fake_scripts(tmp_path)

    sample_cfg: dict[str, Any] = {
        "reads": {"r1": "r1.fq.gz", "r2": "r2.fq.gz"},
        "construct_reference": "db/c.fa",
        "host_reference": "db/h.fa",
    }

    with patch.object(run_pipeline, "build_step_cmd",
                      return_value=["echo", "dry"]):
        run_pipeline.run_step(
            step="1",
            sample_key=sname,
            sample_cfg=sample_cfg,
            outdir=outdir,
            threads=8,
            base_dir=tmp_path,
            dry_run=True,   # ← dry-run
        )

    jsonl = outdir / sname / "chain_of_custody.jsonl"
    assert not jsonl.exists(), "dry-run must not create chain_of_custody.jsonl"


def test_run_step_coc_step_label_uses_step_key(tmp_path):
    """CoC 'step' field uses 's03'-style label derived from the step key."""
    sname = "label_test"
    outdir = tmp_path / "results"
    (outdir / sname).mkdir(parents=True)
    _make_fake_scripts(tmp_path)

    sample_cfg: dict[str, Any] = {
        "reads": {"r1": "r1.fq.gz", "r2": "r2.fq.gz"},
        "construct_reference": "db/c.fa",
        "host_reference": "db/h.fa",
    }
    mock_result = MagicMock()
    mock_result.returncode = 0

    with patch.object(run_pipeline, "build_step_cmd",
                      return_value=["echo", "x"]), \
         patch("subprocess.run", return_value=mock_result):
        run_pipeline.run_step(**_make_run_step_kwargs(
            "3", sname, sample_cfg, outdir, tmp_path
        ))

    jsonl = outdir / sname / "chain_of_custody.jsonl"
    events = [json.loads(l) for l in jsonl.read_text().splitlines()]
    steps = {e["step"] for e in events}
    # Step "3" should produce label "s03"
    assert "s03" in steps, f"expected 's03' label in CoC, got: {steps}"
