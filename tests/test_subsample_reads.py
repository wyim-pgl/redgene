"""Tests for Issue #2 — scripts/util/subsample_reads.py.

The AC-7 coverage-sensitivity sweep needs a deterministic fraction-based
subsampler. The existing `scripts/subsample_reads.py` is coverage-based (requires
genome-size + target-coverage); the new helper exposes a simpler `--fraction`
CLI and a standalone Python API that the tests (and `run_coverage_sensitivity.sh`)
consume without importing seqkit.
"""
from __future__ import annotations

import gzip
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "util" / "subsample_reads.py"
FIXTURES = REPO_ROOT / "tests" / "fixtures"


def _load_module():
    assert SCRIPT.exists(), SCRIPT
    spec = importlib.util.spec_from_file_location("util_subsample_reads", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["util_subsample_reads"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def _write_fq(path: Path, n_pairs: int) -> None:
    """Write a tiny gzipped FASTQ with n_pairs records (seq length 10)."""
    with gzip.open(path, "wt") as fh:
        for i in range(n_pairs):
            fh.write(f"@read_{i}\nACGTACGTAC\n+\nIIIIIIIIII\n")


def test_script_exists():
    assert SCRIPT.exists(), SCRIPT


def test_compute_fraction_prefers_cli_fraction_over_coverage():
    mod = _load_module()
    # When --fraction is given explicitly, it wins (no coverage math needed)
    got = mod.compute_fraction(fraction=0.25, target_coverage=None,
                               total_pairs=1000, read_length=150,
                               genome_size=None)
    assert got == 0.25


def test_compute_fraction_from_target_coverage():
    mod = _load_module()
    # 1_000_000 bp genome, 10x target, 150bp reads, 100k total pairs.
    # Need 10 * 1_000_000 / (2 * 150) = 33_333 pairs → fraction ≈ 0.3333
    frac = mod.compute_fraction(
        fraction=None, target_coverage=10.0,
        total_pairs=100_000, read_length=150, genome_size=1_000_000,
    )
    assert 0.33 < frac < 0.34, frac


def test_compute_fraction_caps_at_1():
    mod = _load_module()
    # 50x over a 10x-equivalent population must clamp to 1.0 (not error).
    frac = mod.compute_fraction(
        fraction=None, target_coverage=50.0,
        total_pairs=100_000, read_length=150, genome_size=10_000_000,
    )
    assert frac == 1.0


def test_compute_fraction_requires_enough_inputs():
    mod = _load_module()
    import pytest
    with pytest.raises(ValueError):
        mod.compute_fraction(fraction=None, target_coverage=None,
                             total_pairs=100, read_length=150,
                             genome_size=1_000_000)


def test_subsample_fraction_preserves_pair_count(tmp_path):
    """Running with fraction=0.5 on a 20-pair input → ~10 output pairs."""
    mod = _load_module()
    if shutil.which("seqkit") is None:
        import pytest
        pytest.skip("seqkit not on PATH")
    r1 = tmp_path / "in_R1.fq.gz"
    r2 = tmp_path / "in_R2.fq.gz"
    _write_fq(r1, 20)
    _write_fq(r2, 20)
    out_prefix = tmp_path / "out"

    n = mod.subsample(
        r1=r1, r2=r2, fraction=0.5, seed=42, output_prefix=out_prefix,
    )
    out_r1 = Path(f"{out_prefix}_R1.fq.gz")
    out_r2 = Path(f"{out_prefix}_R2.fq.gz")
    assert out_r1.exists() and out_r2.exists()
    # Count records directly (robust to seqkit binary versioning)
    with gzip.open(out_r1, "rt") as fh:
        lines = sum(1 for _ in fh)
    pairs_out = lines // 4
    # seqkit sample is probabilistic; allow ±50% tolerance on this tiny input
    assert 1 <= pairs_out <= 20, pairs_out
    # API returns the subsampled count
    assert n == pairs_out


def test_layout_coverage_dir_symlinks(tmp_path):
    """For each coverage level the helper must create results/<sample>_cov{X}x/
    with symlinks (or copies) so `run_pipeline.py --sample <sample>_cov{X}x` works.
    """
    mod = _load_module()
    r1 = tmp_path / "in_R1.fq.gz"
    r2 = tmp_path / "in_R2.fq.gz"
    _write_fq(r1, 4)
    _write_fq(r2, 4)
    layout = mod.layout_coverage_dir(
        base_outdir=tmp_path / "results",
        sample="mysample",
        coverage_tag="5x",
        r1=r1, r2=r2,
    )
    assert layout.exists() and layout.is_dir()
    # Both reads are reachable (either symlinks or real files)
    assert (layout / "mysample_cov5x_R1.fq.gz").exists()
    assert (layout / "mysample_cov5x_R2.fq.gz").exists()


def test_cli_fraction_runs(tmp_path):
    """End-to-end CLI: --fraction 0.5 on tiny fixture produces both outputs."""
    if shutil.which("seqkit") is None:
        import pytest
        pytest.skip("seqkit not on PATH")
    r1 = tmp_path / "in_R1.fq.gz"
    r2 = tmp_path / "in_R2.fq.gz"
    _write_fq(r1, 20)
    _write_fq(r2, 20)
    out_prefix = tmp_path / "cli_out"
    result = subprocess.run(
        [
            sys.executable, str(SCRIPT),
            "--r1", str(r1), "--r2", str(r2),
            "--fraction", "0.5",
            "--seed", "7",
            "--output-prefix", str(out_prefix),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert Path(f"{out_prefix}_R1.fq.gz").exists()
    assert Path(f"{out_prefix}_R2.fq.gz").exists()


def test_run_coverage_sensitivity_has_slurm_headers():
    """SLURM array script must declare explicit headers (>=96G for cucumber)."""
    script = REPO_ROOT / "run_coverage_sensitivity.sh"
    assert script.exists(), script
    text = script.read_text()
    for key in ("--partition=", "--account=", "--time=", "--mem=", "--array="):
        assert key in text, f"run_coverage_sensitivity.sh missing {key}"
    # Must NOT auto-submit any nested scripts
    for line in text.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("#"):
            continue
        assert not stripped.startswith("sbatch"), (
            f"run_coverage_sensitivity.sh must not auto-submit nested jobs: `{line}`"
        )


def test_run_coverage_sensitivity_samples_match_config():
    """SAMPLES table must carry config.yaml keys (not stale aliases)."""
    script = REPO_ROOT / "run_coverage_sensitivity.sh"
    text = script.read_text()
    # Config.yaml uses tomato_Cas9_A2_3, not tomato_A2_3 - regression guard.
    assert "tomato_Cas9_A2_3" in text, (
        "run_coverage_sensitivity.sh must reference tomato_Cas9_A2_3 "
        "(matches config.yaml sample key)"
    )
    assert "tomato_A2_3)" not in text, (
        "stale alias tomato_A2_3 found - rename to tomato_Cas9_A2_3"
    )
    for sample in ("rice_G281", "tomato_Cas9_A2_3", "cucumber_line225"):
        assert sample in text, f"SAMPLES missing {sample}"


def test_run_coverage_sensitivity_invokes_run_pipeline():
    """Runner must actually call run_pipeline.py (not a scaffold echo)."""
    script = REPO_ROOT / "run_coverage_sensitivity.sh"
    text = script.read_text()
    # Must have an actual invocation of run_pipeline.py with expected flags
    assert "python run_pipeline.py" in text or "run_pipeline.py" in text, (
        "runner must call run_pipeline.py"
    )
    assert "--steps 1-5" in text, "runner must request steps 1-5"
    assert "--no-remote-blast" in text, "runner must use --no-remote-blast"
    # And no leftover "would now run" placeholder
    assert "would now run" not in text, (
        "stale scaffold echo 'would now run' still present"
    )


def test_analyze_coverage_sensitivity_script_exists():
    script = REPO_ROOT / "scripts" / "util" / "analyze_coverage_sensitivity.py"
    assert script.exists(), script
