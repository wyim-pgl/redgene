"""Tests for Issue #6 — scripts/reports/insertion_pdf.py scaffold.

v1.1 PDF report generator. Scaffold only: validates that the PDF generator
produces a non-trivial multi-page PDF with the seven required sections, using
matplotlib's PdfPages backend (no reportlab dependency added).
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "reports" / "insertion_pdf.py"


def _load_module():
    assert SCRIPT.exists(), SCRIPT
    spec = importlib.util.spec_from_file_location("insertion_pdf", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["insertion_pdf"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


_AUDIT_HEADER = {
    "sample": "fixtureX",
    "generated_at": "2026-04-17T00:00:00+00:00",
    "input_sha256": {"r1": "a" * 64, "r2": "b" * 64},
    "pipeline_commit": "deadbeefcafe",
    "pipeline_dirty": False,
    "db_manifest": [
        {"name": "element_db", "md5": "e" * 32, "build_date": "2026-04-16",
         "seq_count": "106"}
    ],
    "software_versions": {"python": "Python 3.11.15"},
}

_REPORT_CAND = """\
======================================================================
RedGene Insert Assembly & Annotation Report
======================================================================
Insertion site: chrX:100-150 (50bp deletion)
Insert length: 200 bp
Verdict: CANDIDATE
T-DNA borders found: 4
======================================================================
"""

_REPORT_TRUE = """\
======================================================================
RedGene Insert Assembly & Annotation Report
======================================================================
Insertion site: chrY:200-250 (50bp deletion)
Verdict: TRUE_INSERTION
======================================================================
"""


def _make_results(tmp_path: Path) -> Path:
    sample_dir = tmp_path / "fixtureX"
    (sample_dir / "s05_insert_assembly").mkdir(parents=True)
    (sample_dir / "audit_header.json").write_text(json.dumps(_AUDIT_HEADER))
    (sample_dir / "s05_insert_assembly" / "insertion_chrX_100_report.txt").write_text(_REPORT_CAND)
    (sample_dir / "s05_insert_assembly" / "insertion_chrY_200_report.txt").write_text(_REPORT_TRUE)
    return sample_dir


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_script_exists():
    assert SCRIPT.exists(), SCRIPT


def test_load_audit_header(tmp_path):
    mod = _load_module()
    sample_dir = _make_results(tmp_path)
    hdr = mod.load_audit_header(sample_dir)
    assert hdr["sample"] == "fixtureX"
    assert hdr["pipeline_commit"] == "deadbeefcafe"


def test_parse_insertion_reports(tmp_path):
    mod = _load_module()
    sample_dir = _make_results(tmp_path)
    rows = mod.parse_insertion_reports(sample_dir / "s05_insert_assembly")
    assert len(rows) == 2
    verdicts = sorted(r["verdict"] for r in rows)
    assert verdicts == ["CANDIDATE", "TRUE_INSERTION"]
    # Each row must carry a site string and a source file
    for r in rows:
        assert "site" in r and r["site"]
        assert "source" in r


def test_summary_counts(tmp_path):
    mod = _load_module()
    sample_dir = _make_results(tmp_path)
    rows = mod.parse_insertion_reports(sample_dir / "s05_insert_assembly")
    summary = mod.summarize(rows)
    assert summary["n_sites"] == 2
    assert summary["by_verdict"]["CANDIDATE"] == 1
    assert summary["by_verdict"]["TRUE_INSERTION"] == 1


def test_generate_pdf_produces_nonzero_output(tmp_path):
    """Smoke: generate_pdf emits a multi-page PDF (non-zero size, starts with %PDF)."""
    mod = _load_module()
    sample_dir = _make_results(tmp_path)
    out = tmp_path / "out.pdf"
    mod.generate_pdf(sample_dir=sample_dir, sample_name="fixtureX", out_pdf=out)
    assert out.exists(), out
    assert out.stat().st_size > 1000, "PDF suspiciously small"
    # Minimal sanity check: PDF magic bytes
    with open(out, "rb") as fh:
        assert fh.read(4) == b"%PDF"


def test_cli_generate(tmp_path):
    sample_dir = _make_results(tmp_path)
    out = tmp_path / "cli.pdf"
    result = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--sample", "fixtureX",
            "--sample-dir", str(sample_dir),
            "--out", str(out),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out.exists() and out.stat().st_size > 1000


def test_generate_pdf_handles_missing_audit_header(tmp_path):
    """If audit_header.json is missing, the PDF still generates (with a placeholder)."""
    mod = _load_module()
    sample_dir = tmp_path / "no_audit"
    (sample_dir / "s05_insert_assembly").mkdir(parents=True)
    (sample_dir / "s05_insert_assembly" / "insertion_chrA_1_report.txt").write_text(_REPORT_CAND)
    out = tmp_path / "no_audit.pdf"
    mod.generate_pdf(sample_dir=sample_dir, sample_name="no_audit", out_pdf=out)
    assert out.exists() and out.stat().st_size > 500


def test_generate_pdf_writes_seven_pages(tmp_path):
    """The scaffold spec requires 7 sections → at least 7 PDF pages."""
    mod = _load_module()
    sample_dir = _make_results(tmp_path)
    out = tmp_path / "seven.pdf"
    n_pages = mod.generate_pdf(
        sample_dir=sample_dir, sample_name="fixtureX", out_pdf=out,
    )
    assert n_pages >= 7, f"expected ≥7 sections/pages, got {n_pages}"
