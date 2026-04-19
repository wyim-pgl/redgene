"""Tests for scripts/util/analyze_coverage_sensitivity.py GT anchor matching.

Covers the Issue #2 v1.1 extension (2026-04-18) that resolves gt_anchor_hit
against ground_truth_baseline.tsv. The analyzer must:

1. Return 'HIT:<verdict>' when insertion_<chrom>_<pos>_report.txt exists for
   a GT coordinate (chrom matching with NCBI version suffix tolerance).
2. Return 'MISS' when the sample has a GT row but no matching report.
3. Return 'no_gt' for samples absent from the baseline TSV.
4. Distinguish 'pipeline complete, 0 insertion reports' from 'no s05 output
   (pipeline incomplete?)' based on whether s05 dir exists.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "analyze_coverage_sensitivity",
    REPO_ROOT / "scripts" / "util" / "analyze_coverage_sensitivity.py",
)
acs = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(acs)


def _write_report(dst: Path, verdict: str) -> Path:
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(
        "======================================================================\n"
        "RedGene Insert Assembly & Annotation Report\n"
        "======================================================================\n"
        f"Verdict: {verdict} — stub\n"
    )
    return dst


def test_chrom_equal_ncbi_version_suffix():
    assert acs._chrom_equal("LKUO03001451", "LKUO03001451.1")
    assert acs._chrom_equal("LKUO03001451.1", "LKUO03001451")
    assert acs._chrom_equal("Chr3", "Chr3")
    assert not acs._chrom_equal("Chr3", "Chr4")
    assert not acs._chrom_equal("LKUO03001451", "LKUO03001452")


def test_load_ground_truth_parses_tsv(tmp_path: Path):
    tsv = tmp_path / "gt.tsv"
    tsv.write_text(
        "sample\tchrom\tposition\tdescription\n"
        "rice_G281\tChr3\t16439674\tT-DNA 2-copy\n"
        "cucumber_line225\tLKUO03001451\t6501\tline 225\n"
    )
    gt = acs._load_ground_truth(tsv)
    assert gt["rice_G281"] == [("Chr3", 16439674)]
    assert gt["cucumber_line225"] == [("LKUO03001451", 6501)]


def test_load_ground_truth_missing_file_returns_empty(tmp_path: Path):
    assert acs._load_ground_truth(tmp_path / "does_not_exist.tsv") == {}
    assert acs._load_ground_truth(None) == {}


def test_match_gt_report_hit_exact(tmp_path: Path):
    s05 = tmp_path / "s05_insert_assembly"
    _write_report(s05 / "insertion_Chr3_16439674_report.txt", "CANDIDATE")
    _write_report(s05 / "insertion_Chr1_1000_report.txt", "UNKNOWN")
    status, path = acs._match_gt_report(s05, [("Chr3", 16439674)])
    assert status == "HIT:CANDIDATE"
    assert path.name == "insertion_Chr3_16439674_report.txt"


def test_match_gt_report_hit_ncbi_suffix(tmp_path: Path):
    """Baseline TSV uses LKUO03001451, filename uses LKUO03001451.1."""
    s05 = tmp_path / "s05_insert_assembly"
    _write_report(s05 / "insertion_LKUO03001451.1_6501_report.txt", "CANDIDATE")
    status, path = acs._match_gt_report(s05, [("LKUO03001451", 6501)])
    assert status == "HIT:CANDIDATE"


def test_match_gt_report_hit_chrom_with_underscore(tmp_path: Path):
    """Tomato chrom SLM_r2.0ch01 contains both underscore and dot."""
    s05 = tmp_path / "s05_insert_assembly"
    _write_report(s05 / "insertion_SLM_r2.0ch01_91002744_report.txt", "CANDIDATE")
    status, path = acs._match_gt_report(s05, [("SLM_r2.0ch01", 91002744)])
    assert status == "HIT:CANDIDATE"


def test_match_gt_report_miss_when_no_matching_report(tmp_path: Path):
    s05 = tmp_path / "s05_insert_assembly"
    _write_report(s05 / "insertion_Chr1_1000_report.txt", "UNKNOWN")
    status, path = acs._match_gt_report(s05, [("Chr3", 16439674)])
    assert status == "MISS"
    assert path is None


def test_match_gt_report_miss_when_dir_missing(tmp_path: Path):
    s05 = tmp_path / "does_not_exist"
    status, path = acs._match_gt_report(s05, [("Chr3", 16439674)])
    assert status == "MISS"


def test_analyze_sample_cov_no_gt_for_sample(tmp_path: Path):
    s05 = tmp_path / "results" / "sample_cov5x" / "s05_insert_assembly"
    _write_report(s05 / "insertion_X_1_report.txt", "CANDIDATE")
    row = acs.analyze_sample_cov(
        tmp_path / "results", "sample", "5x", ground_truth={"other_sample": [("X", 1)]}
    )
    assert row["gt_anchor_hit"] == "no_gt"


def test_analyze_sample_cov_note_distinguishes_empty_vs_missing(tmp_path: Path):
    """s05 dir exists but empty → 'pipeline complete, 0 insertion reports'."""
    s05 = tmp_path / "results" / "sample_cov5x" / "s05_insert_assembly"
    s05.mkdir(parents=True)
    row = acs.analyze_sample_cov(
        tmp_path / "results", "sample", "5x", ground_truth={}
    )
    assert row["n_sites_total"] == 0
    assert row["notes"] == "pipeline complete, 0 insertion reports"

    # No s05 dir → "pipeline incomplete"
    row2 = acs.analyze_sample_cov(
        tmp_path / "results", "sample", "99x", ground_truth={}
    )
    assert row2["notes"] == "no s05 output (pipeline incomplete?)"
