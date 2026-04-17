"""Regression test for Issue #15 — run_rerun_w1_batch.sh SBATCH headers.

BUG-18: cucumber_line212 (MaxRSS 67.1G) and line224 (65.7G) OOM'd under the
cluster-default memory limit because run_rerun_w1_batch.sh lacked explicit
SBATCH headers. Fix: add explicit `--mem=96G` (covers cucumber), plus
partition/account/cpus/time headers for reproducibility.
"""
from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[1]


def _headers(path: Path) -> dict[str, str]:
    """Extract #SBATCH --key=value lines into a dict."""
    out: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"#SBATCH\s+--([a-z-]+)=?(.*)", line.strip())
            if m:
                out[m.group(1)] = m.group(2).strip()
            if line.startswith(("set ", "eval ", "cd ")):
                break  # end of SBATCH preamble
    return out


def test_rerun_w1_batch_has_mem_ge_96g():
    """BUG-18 fix: explicit --mem must be >= 96G to accommodate cucumber."""
    script = REPO_ROOT / "run_rerun_w1_batch.sh"
    assert script.exists(), script
    headers = _headers(script)
    assert "mem" in headers, "missing --mem header (root cause of BUG-18)"
    mem = headers["mem"].rstrip("Gg")
    assert int(mem) >= 96, f"--mem must be >=96G for cucumber, got {headers['mem']}"


def test_rerun_w1_batch_has_cpus_16():
    script = REPO_ROOT / "run_rerun_w1_batch.sh"
    headers = _headers(script)
    assert "cpus-per-task" in headers
    assert int(headers["cpus-per-task"]) == 16


def test_rerun_w1_batch_has_partition():
    script = REPO_ROOT / "run_rerun_w1_batch.sh"
    headers = _headers(script)
    assert "partition" in headers, "missing --partition (reproducibility)"


def test_rerun_w1_batch_has_time():
    script = REPO_ROOT / "run_rerun_w1_batch.sh"
    headers = _headers(script)
    assert "time" in headers, "missing --time (reproducibility)"
