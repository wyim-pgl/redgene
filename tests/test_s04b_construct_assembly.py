"""Tests for scripts/s04b_construct_assembly.py"""
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def test_emits_contigs_fasta(tmp_path):
    r1 = REPO / "tests/fixtures/mini_R1.fq.gz"
    r2 = REPO / "tests/fixtures/mini_R2.fq.gz"
    outdir = tmp_path / "out"
    subprocess.run(
        ["python", str(REPO / "scripts/s04b_construct_assembly.py"),
         "--r1", str(r1), "--r2", str(r2),
         "--outdir", str(outdir),
         "--sample-name", "mini",
         "--threads", "2"],
        check=True,
    )
    contigs = outdir / "mini" / "s04b_construct_asm" / "contigs.fasta"
    assert contigs.exists()


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
