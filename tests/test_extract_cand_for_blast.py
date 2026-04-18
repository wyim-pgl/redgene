"""Tests for Issue #1 — scripts/util/extract_cand_for_blast.py.

Covers extraction of CANDIDATE insertion reports and merging of their
insert FASTA files into a single per-sample query FASTA used by the
remote-BLAST FP verification workflow (AC-2).
"""
from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "util" / "extract_cand_for_blast.py"


def _load_module():
    """Import the helper as a fresh module (script has no package path)."""
    assert SCRIPT.exists(), SCRIPT
    spec = importlib.util.spec_from_file_location(
        "extract_cand_for_blast", SCRIPT
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["extract_cand_for_blast"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


_REPORT_CAND = """\
======================================================================
RedGene Insert Assembly & Annotation Report
======================================================================
Insertion site: chrX:100-150 (50bp deletion)
Insert length: 200 bp
Assembly status: CONVERGED (round 2)
Structure: single-copy
Verdict: CANDIDATE
T-DNA borders found: 4
======================================================================
"""

_REPORT_UNKNOWN = """\
======================================================================
RedGene Insert Assembly & Annotation Report
======================================================================
Insertion site: chrY:200-250 (50bp deletion)
Verdict: UNKNOWN — no element annotations
======================================================================
"""

_REPORT_TRUE = """\
======================================================================
RedGene Insert Assembly & Annotation Report
======================================================================
Insertion site: chrZ:300-350 (50bp deletion)
Verdict: TRUE_INSERTION
======================================================================
"""


def _write_site(dir_: Path, site: str, verdict_text: str, seq: str) -> None:
    (dir_ / f"insertion_{site}_report.txt").write_text(verdict_text)
    (dir_ / f"insertion_{site}_insert.fasta").write_text(
        f">insert_{site}\n{seq}\n"
    )


def _make_sample(tmp_path: Path) -> Path:
    sample = tmp_path / "sample_X"
    s05 = sample / "s05_insert_assembly"
    s05.mkdir(parents=True)
    _write_site(s05, "chrX_100", _REPORT_CAND, "ACGT" * 10)
    _write_site(s05, "chrZ_300", _REPORT_TRUE, "GGGG" * 10)
    _write_site(s05, "chrY_200", _REPORT_UNKNOWN, "TTTT" * 10)
    # Second CAND to exercise multi-record merge
    _write_site(s05, "chrX_900", _REPORT_CAND, "CCCC" * 10)
    return sample


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_script_exists():
    assert SCRIPT.exists(), SCRIPT


def test_find_candidate_reports_returns_only_cand(tmp_path):
    mod = _load_module()
    sample = _make_sample(tmp_path)
    reports = mod.find_candidate_reports(sample / "s05_insert_assembly")
    names = sorted(p.name for p in reports)
    assert names == [
        "insertion_chrX_100_report.txt",
        "insertion_chrX_900_report.txt",
    ], names


def test_report_to_insert_fasta_path(tmp_path):
    mod = _load_module()
    sample = _make_sample(tmp_path)
    s05 = sample / "s05_insert_assembly"
    report = s05 / "insertion_chrX_100_report.txt"
    fasta = mod.report_to_insert_fasta(report)
    assert fasta.name == "insertion_chrX_100_insert.fasta", fasta
    assert fasta.exists()


def test_merge_inserts_produces_concatenated_fasta(tmp_path):
    mod = _load_module()
    sample = _make_sample(tmp_path)
    out = tmp_path / "out_cand.fa"
    count = mod.merge_cand_inserts(
        sample_dir=sample, out_fasta=out, sample_name="sample_X"
    )
    assert count == 2, count
    text = out.read_text()
    # Each record re-headered with sample prefix so BLAST hit TSV is traceable
    assert text.count(">") == 2
    assert "sample_X" in text
    # Body preserved (at minimum one of the inserts)
    assert "ACGT" in text


def test_merge_inserts_handles_no_candidates(tmp_path):
    """If no CAND reports exist we must write a valid empty FASTA and return 0."""
    mod = _load_module()
    sample = tmp_path / "sample_empty"
    (sample / "s05_insert_assembly").mkdir(parents=True)
    out = tmp_path / "empty.fa"
    count = mod.merge_cand_inserts(
        sample_dir=sample, out_fasta=out, sample_name="sample_empty"
    )
    assert count == 0
    assert out.exists() and out.stat().st_size == 0


def test_prepare_remote_blast_has_sbatch_headers():
    """Companion SLURM template must carry explicit headers (no login-node runs)."""
    script = REPO_ROOT / "scripts" / "util" / "prepare_remote_blast.sh"
    assert script.exists(), script
    text = script.read_text()
    for key in ("--partition=", "--account=", "--time=", "--mem="):
        assert key in text, f"prepare_remote_blast.sh missing SBATCH header {key}"
    # Must NOT *call* sbatch itself (submission is operator-driven). Allow
    # comment mentions like `sbatch prepare_remote_blast.sh` in docstrings
    # — only flag uncommented lines that actually invoke sbatch.
    for line in text.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("#"):
            continue
        assert not stripped.startswith("sbatch"), (
            f"prepare_remote_blast.sh must not auto-submit: `{line}`"
        )
    # Must reference blastn -remote (the whole point)
    assert "blastn" in text and "-remote" in text


def test_script_has_cli_main(tmp_path):
    """CLI entry point: --sample-dir / --out / --sample-name plumbing."""
    import subprocess
    sample = _make_sample(tmp_path)
    out = tmp_path / "cli_cand.fa"
    result = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--sample-dir", str(sample),
            "--out", str(out),
            "--sample-name", "sample_X",
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out.exists()
    assert out.read_text().count(">") == 2
