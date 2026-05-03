"""Tests for scripts/util/assign_kcgp_ids.py — KCGP mapping CLI."""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "util" / "assign_kcgp_ids.py"


def _write_config(path: Path, sample: str, host: str, event: str | None) -> None:
    body = f"""samples:
  {sample}:
    host_reference: {host}
"""
    if event is not None:
        body += f"    event: {event}\n"
    path.write_text(body)


def _write_report(s05_dir: Path, host_chr: str, pos: int, verdict: str) -> None:
    s05_dir.mkdir(parents=True, exist_ok=True)
    path = s05_dir / f"insertion_{host_chr}_{pos}_report.txt"
    path.write_text(
        "Header line\n"
        "More header\n"
        "\n"
        "Soft-clip junction stats\n"
        "  ...\n"
        "  ...\n"
        "  ...\n"
        f"Verdict: {verdict} — explanation text\n"
    )


def _run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True, text=True,
    )


def _read_mapping(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_no_s05_dir(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/Osativa_323_v7.0.fa", "EVT")
    sample_dir = tmp_path / "results" / "sampX"
    sample_dir.mkdir(parents=True)

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode == 0, res.stderr
    rows = _read_mapping(sample_dir / "kcgp_mapping.tsv")
    assert rows == []


def test_no_reports(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/Osativa_323_v7.0.fa", "EVT")
    sample_dir = tmp_path / "results" / "sampX"
    (sample_dir / "s05_insert_assembly").mkdir(parents=True)

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode == 0, res.stderr
    assert _read_mapping(sample_dir / "kcgp_mapping.tsv") == []


def test_assigns_ordinal_to_candidate_only(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/Osativa_323_v7.0.fa", "EVT")
    sample_dir = tmp_path / "results" / "sampX"
    s05 = sample_dir / "s05_insert_assembly"
    _write_report(s05, "Chr1", 100, "CANDIDATE")
    _write_report(s05, "Chr2", 200, "FALSE_POSITIVE")
    _write_report(s05, "Chr3", 300, "UNKNOWN")

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode == 0, res.stderr
    rows = _read_mapping(sample_dir / "kcgp_mapping.tsv")
    by_chr = {r["host_chr"]: r for r in rows}
    assert by_chr["Chr1"]["kcgp_id"] == "KCGP-OSA-EVT-26-0001"
    assert by_chr["Chr1"]["verdict"] == "CANDIDATE"
    assert by_chr["Chr2"]["kcgp_id"] == "-"
    assert by_chr["Chr2"]["verdict"] == "FALSE_POSITIVE"
    assert by_chr["Chr3"]["kcgp_id"] == "-"
    assert by_chr["Chr3"]["verdict"] == "UNKNOWN"
    for r in rows:
        assert r["spec_version"] == "tentative_v0"


def test_sort_order_natural(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/Osativa_323_v7.0.fa", "EVT")
    sample_dir = tmp_path / "results" / "sampX"
    s05 = sample_dir / "s05_insert_assembly"
    _write_report(s05, "Chr2", 100, "CANDIDATE")
    _write_report(s05, "Chr10", 200, "CANDIDATE")
    _write_report(s05, "Chr1", 50, "CANDIDATE")

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode == 0, res.stderr
    rows = _read_mapping(sample_dir / "kcgp_mapping.tsv")
    candidate_rows = [r for r in rows if r["verdict"] == "CANDIDATE"]
    assert [r["host_chr"] for r in candidate_rows] == ["Chr1", "Chr2", "Chr10"]
    assert [r["kcgp_id"] for r in candidate_rows] == [
        "KCGP-OSA-EVT-26-0001",
        "KCGP-OSA-EVT-26-0002",
        "KCGP-OSA-EVT-26-0003",
    ]


def test_missing_event_in_config(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/Osativa_323_v7.0.fa", event=None)
    sample_dir = tmp_path / "results" / "sampX"
    sample_dir.mkdir(parents=True)

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode != 0
    assert "event" in res.stderr.lower()


def test_unknown_host_in_config(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/unknown_genome.fa", "EVT")
    sample_dir = tmp_path / "results" / "sampX"
    sample_dir.mkdir(parents=True)

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode != 0
    err = res.stderr.lower()
    assert "host" in err or "osa" in err
