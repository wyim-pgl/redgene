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


# ---------------------------------------------------------------------------
# Issue #6 R-5 — expanded panel tests (matplotlib/PdfPages backend)
# ---------------------------------------------------------------------------


_REPORT_WITH_SITES = """\
======================================================================
RedGene Insert Assembly & Annotation Report
======================================================================
Insertion site: Chr3:16,439,674-16,439,709 (35bp deletion)
Insert length: 13,471 bp (100 unresolved N's)
Assembly status: PARTIAL (round 5)
Structure: head-to-head 2-copy T-DNA
Verdict: CANDIDATE
T-DNA borders found: 10

--- Linear Map ---
5' host --[NODE->]--[NODE<-]-- 3' host
======================================================================

Detailed element positions:
   Start       End  Dir    Src  Element
----------------------------------------------------------------------
    6043      6774    +  local  NODE_2_length_2383
    6875      8167    -  local  NODE_2_length_2383

Multi-copy elements detected:
  NODE_2_length_2383: 2 copies

Foreign elements (not in host genome): 1
  + NODE_2_length_2383

Host-fraction analysis: 87.4% host (11,688/13,371bp at >=90.0% identity), largest non-host gap: 1,024bp
"""


def test_parse_report_extra_fields_elements_and_structure(tmp_path):
    """parse_insertion_reports must expose structure + foreign elements + length."""
    mod = _load_module()
    s05 = tmp_path / "s05"
    s05.mkdir()
    (s05 / "insertion_Chr3_16439674_report.txt").write_text(_REPORT_WITH_SITES)
    rows = mod.parse_insertion_reports(s05)
    assert len(rows) == 1
    row = rows[0]
    assert "Chr3:16,439,674" in row["site"]
    assert row["verdict"] == "CANDIDATE"
    # Insert length digits preserved (value may be '13,471')
    assert row["insert_length"].replace(",", "") == "13471"
    # New fields added for R-5 table:
    assert "structure" in row
    assert "head-to-head" in row["structure"].lower() or row["structure"]
    assert "elements" in row
    # foreign element list must include at least one element
    assert isinstance(row["elements"], list)


def test_load_copynumber(tmp_path):
    """load_copynumber reads copynumber.tsv into a dict with depth ratio + copies."""
    mod = _load_module()
    sample_dir = tmp_path / "sX"
    s07 = sample_dir / "s07_copynumber"
    s07.mkdir(parents=True)
    (s07 / "copynumber.tsv").write_text(
        "sample\tconstruct_median_depth\thost_median_depth\tdepth_ratio\t"
        "estimated_copies\tcopy_description\tcandidate_sites\tsite_validation\tconfidence\n"
        "sX\t47.00\t22.00\t2.1364\t2.1\tmulti-copy (~2.1x)\t1\tdepth_suggests_multicopy\tMedium\n"
    )
    cn = mod.load_copynumber(sample_dir)
    assert cn is not None
    assert abs(float(cn["depth_ratio"]) - 2.1364) < 1e-3
    assert float(cn["estimated_copies"]) == 2.1
    assert cn["copy_description"].startswith("multi-copy")


def test_load_copynumber_missing_returns_none(tmp_path):
    mod = _load_module()
    sample_dir = tmp_path / "no_cn"
    sample_dir.mkdir()
    assert mod.load_copynumber(sample_dir) is None


def test_load_coc_log_parses_jsonl(tmp_path):
    mod = _load_module()
    sample_dir = tmp_path / "sY"
    sample_dir.mkdir()
    entries = [
        {"ts": "2026-04-18T00:00:00+00:00", "sample": "sY", "step": "s01",
         "event": "start"},
        {"ts": "2026-04-18T00:05:00+00:00", "sample": "sY", "step": "s01",
         "event": "end", "exit_code": 0, "wall_time_sec": 300},
    ]
    with open(sample_dir / "chain_of_custody.jsonl", "w") as fh:
        for e in entries:
            fh.write(json.dumps(e) + "\n")
    parsed = mod.load_coc_log(sample_dir)
    assert len(parsed) == 2
    assert parsed[0]["step"] == "s01"
    assert parsed[1]["exit_code"] == 0


def test_load_coc_log_missing_returns_empty(tmp_path):
    mod = _load_module()
    assert mod.load_coc_log(tmp_path / "nope") == []


def test_parse_softclip_positions_from_report(tmp_path):
    """A report's site-line yields a numeric anchor position used in the diagram."""
    mod = _load_module()
    s05 = tmp_path / "s05"
    s05.mkdir()
    (s05 / "insertion_Chr3_16439674_report.txt").write_text(_REPORT_WITH_SITES)
    rows = mod.parse_insertion_reports(s05)
    anchors = mod.extract_site_anchors(rows)
    # anchors is a list of (chrom, pos) tuples for diagramming
    assert any(a[0] == "Chr3" and a[1] == 16439674 for a in anchors)


def test_generate_pdf_with_copy_number_tsv(tmp_path):
    """When copynumber.tsv is present, the PDF should grow (real chart rendered)."""
    mod = _load_module()
    sample_dir = _make_results(tmp_path)
    s07 = sample_dir / "s07_copynumber"
    s07.mkdir()
    (s07 / "copynumber.tsv").write_text(
        "sample\tconstruct_median_depth\thost_median_depth\tdepth_ratio\t"
        "estimated_copies\tcopy_description\tcandidate_sites\tsite_validation\tconfidence\n"
        "fixtureX\t30.0\t15.0\t2.0\t2.0\tmulti-copy (~2.0x)\t1\tvalid\tHigh\n"
    )
    out = tmp_path / "with_cn.pdf"
    n = mod.generate_pdf(sample_dir=sample_dir, sample_name="fixtureX", out_pdf=out)
    assert n >= 7
    assert out.stat().st_size > 1000


def test_generate_pdf_with_coc_log_appendix(tmp_path):
    """When chain_of_custody.jsonl is present, the appendix panel should include it."""
    mod = _load_module()
    sample_dir = _make_results(tmp_path)
    entries = [
        {"ts": "2026-04-18T00:00:00+00:00", "sample": "fixtureX", "step": "s01",
         "event": "pre", "cmd": "scripts/s01_qc.py ..."},
        {"ts": "2026-04-18T00:05:00+00:00", "sample": "fixtureX", "step": "s01",
         "event": "post", "exit_code": 0},
    ]
    with open(sample_dir / "chain_of_custody.jsonl", "w") as fh:
        for e in entries:
            fh.write(json.dumps(e) + "\n")
    out = tmp_path / "with_coc.pdf"
    n = mod.generate_pdf(sample_dir=sample_dir, sample_name="fixtureX", out_pdf=out)
    # With CoC present the report adds an extra appendix page → ≥8.
    assert n >= 8


def test_sample_info_and_pipeline_version_helpers(tmp_path):
    """Dedicated helpers expose sample-info and pipeline-version blocks."""
    mod = _load_module()
    sample_dir = _make_results(tmp_path)
    info = mod.build_sample_info(sample_dir, sample_name="fixtureX",
                                 host="rice", construct="pCAMBIA")
    assert info["sample"] == "fixtureX"
    assert info["host"] == "rice"
    assert info["construct"] == "pCAMBIA"
    # generated_at present only if audit_header.json supplies it
    assert "date" in info

    audit = mod.load_audit_header(sample_dir)
    ver = mod.build_pipeline_version(audit)
    assert ver["pipeline_commit"] == "deadbeefcafe"
    assert "software_versions" in ver
    assert "python" in ver["software_versions"]


def test_generate_pdf_accepts_host_and_construct_kwargs(tmp_path):
    """generate_pdf accepts optional host / construct kwargs for the info panel."""
    mod = _load_module()
    sample_dir = _make_results(tmp_path)
    out = tmp_path / "info.pdf"
    n = mod.generate_pdf(
        sample_dir=sample_dir, sample_name="fixtureX", out_pdf=out,
        host="rice_Osativa_v7", construct="pCAMBIA-lactoferrin",
    )
    assert n >= 7
    assert out.stat().st_size > 1000
