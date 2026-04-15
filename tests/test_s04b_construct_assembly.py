"""Tests for scripts/s04b_construct_assembly.py"""
import importlib.util
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _load_s04b():
    spec = importlib.util.spec_from_file_location(
        "s04b", REPO / "scripts/s04b_construct_assembly.py",
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["s04b"] = mod
    spec.loader.exec_module(mod)
    return mod


# High-complexity 50bp seed (same style as test_extra_element_db); repeated 8x
# gives a 400bp contig that survives BLAST's default dust filter while being
# long enough to clear the 200bp length cutoff.
MARKER_SEED = "ATCGATCGATCGAAGCTTGGATCCAAGCTAGCTAGCTAGAACCGGTTAACC"
MARKER_400 = MARKER_SEED * 8
# Distinct high-complexity seed, 400bp, so contig B won't match MARKER_400.
NONMARKER_SEED = "TGCAAGTTCGATCGTACGTAGCTAGCATCGATCGTTAAGGCCTAGCTTGCA"
NONMARKER_400 = NONMARKER_SEED * 8


def test_spades_failure_on_tiny_input_yields_empty_contigs(tmp_path):
    """3 read pairs is below SPAdes' coverage threshold; wrapper should swallow
    the failure, emit an empty contigs.fasta, and exit 0."""
    r1 = REPO / "tests/fixtures/mini_R1.fq.gz"
    r2 = REPO / "tests/fixtures/mini_R2.fq.gz"
    outdir = tmp_path / "out"
    result = subprocess.run(
        ["python", str(REPO / "scripts/s04b_construct_assembly.py"),
         "--r1", str(r1), "--r2", str(r2),
         "--outdir", str(outdir),
         "--sample-name", "mini",
         "--threads", "2"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    contigs = outdir / "mini" / "s04b_construct_asm" / "contigs.fasta"
    assert contigs.exists()
    assert contigs.stat().st_size == 0, "SPAdes-failure branch should emit empty fasta"


def test_handles_empty_reads_gracefully(tmp_path):
    import gzip
    r1 = tmp_path / "empty_R1.fq.gz"
    r2 = tmp_path / "empty_R2.fq.gz"
    for p in (r1, r2):
        with gzip.open(p, "wt") as fh:
            pass
    outdir = tmp_path / "out"
    result = subprocess.run(
        ["python", str(REPO / "scripts/s04b_construct_assembly.py"),
         "--r1", str(r1), "--r2", str(r2),
         "--outdir", str(outdir),
         "--sample-name", "empty",
         "--threads", "2"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    contigs = outdir / "empty" / "s04b_construct_asm" / "contigs.fasta"
    assert contigs.exists()
    assert contigs.stat().st_size == 0


def test_filter_contigs_keeps_marker_positive_only(tmp_path):
    """_filter_contigs_by_markers must keep a contig that has a >=200bp
    >=90% blastn hit and drop one that has no marker hit at all."""
    s04b = _load_s04b()

    contigs = tmp_path / "all.fa"
    contigs.write_text(
        f">node_A\n{MARKER_400}\n"
        f">node_B\n{NONMARKER_400}\n"
    )
    marker = tmp_path / "marker.fa"
    marker.write_text(f">bar_like\n{MARKER_400}\n")
    out = tmp_path / "filtered.fa"

    kept, total = s04b._filter_contigs_by_markers(
        contigs, [marker], out,
        min_identity=90.0, min_aln_len=200,
    )
    assert total == 2
    assert kept == 1
    out_text = out.read_text()
    assert ">node_A" in out_text, out_text
    assert ">node_B" not in out_text, out_text


def test_filter_contigs_handles_missing_marker_db(tmp_path):
    """A missing marker DB is a warning, not a failure — the other DBs still
    get consulted. With zero valid DBs, nothing passes (expected)."""
    s04b = _load_s04b()

    contigs = tmp_path / "all.fa"
    contigs.write_text(f">node_A\n{MARKER_400}\n")
    out = tmp_path / "filtered.fa"

    # Only DB is a missing path
    kept, total = s04b._filter_contigs_by_markers(
        contigs, [tmp_path / "does_not_exist.fa"], out,
        min_identity=90.0, min_aln_len=200,
    )
    assert total == 1
    assert kept == 0
    assert out.read_text() == ""


def test_no_filter_keeps_raw_contigs(tmp_path):
    """--no-filter must skip the filter step; contigs_all.fasta should NOT
    be produced, and contigs.fasta is whatever SPAdes emitted (here: empty
    because the tiny fixture is below SPAdes' coverage threshold)."""
    r1 = REPO / "tests/fixtures/mini_R1.fq.gz"
    r2 = REPO / "tests/fixtures/mini_R2.fq.gz"
    outdir = tmp_path / "out"
    result = subprocess.run(
        ["python", str(REPO / "scripts/s04b_construct_assembly.py"),
         "--r1", str(r1), "--r2", str(r2),
         "--outdir", str(outdir),
         "--sample-name", "nofilt",
         "--threads", "2",
         "--no-filter"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    step_dir = outdir / "nofilt" / "s04b_construct_asm"
    assert (step_dir / "contigs.fasta").exists()
    assert not (step_dir / "contigs_all.fasta").exists(), \
        "--no-filter must not create contigs_all.fasta"
