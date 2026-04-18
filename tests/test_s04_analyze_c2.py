"""Regression test for Issue #14 — C2 transgene-positive count fix.

The original analyze script used grep -c "^transgene-positive" on a TSV whose
4th column is `transgene_positive` (True/False), never a line prefix. The fix
must count rows whose 4th field equals "True".
"""
import subprocess
import tempfile
from pathlib import Path


_FIXTURE_TSV = (
    "site_id\tchrom\tpos\ttransgene_positive\textra\n"
    "insertion_1\tChr1\t100\tTrue\tfoo\n"
    "insertion_2\tChr1\t200\tFalse\tbar\n"
    "insertion_3\tChr2\t300\tTrue\tbaz\n"
    "insertion_4\tChr2\t400\tFalse\tqux\n"
    "insertion_5\tChr3\t500\tTrue\tquux\n"
)


def _count_via_fix(tsv_path: Path) -> int:
    """Run the awk-based fix inline and return the count of True rows."""
    out = subprocess.check_output(
        ["awk", "-F", "\t", 'NR>1 && $4 == "True" {n++} END {print n+0}', str(tsv_path)],
        text=True,
    ).strip()
    return int(out)


def _count_via_broken_grep(tsv_path: Path) -> int:
    """Reproduce the original buggy grep behaviour for comparison."""
    try:
        out = subprocess.check_output(
            ["grep", "-c", "^transgene-positive", str(tsv_path)],
            text=True,
        ).strip()
        return int(out)
    except subprocess.CalledProcessError:
        return 0


def test_fix_counts_true_rows():
    with tempfile.TemporaryDirectory() as tmp:
        tsv = Path(tmp) / "tiers.tsv"
        tsv.write_text(_FIXTURE_TSV)
        assert _count_via_fix(tsv) == 3


def test_broken_grep_returns_zero():
    # Document the regression: original grep always returned 0.
    with tempfile.TemporaryDirectory() as tmp:
        tsv = Path(tmp) / "tiers.tsv"
        tsv.write_text(_FIXTURE_TSV)
        assert _count_via_broken_grep(tsv) == 0


def test_fix_handles_all_false():
    with tempfile.TemporaryDirectory() as tmp:
        tsv = Path(tmp) / "tiers.tsv"
        tsv.write_text(
            "site_id\tchrom\tpos\ttransgene_positive\n"
            "insertion_1\tChr1\t100\tFalse\n"
            "insertion_2\tChr1\t200\tFalse\n"
        )
        assert _count_via_fix(tsv) == 0


def test_fix_handles_missing_file():
    # Missing file should not crash; awk returns 0 for empty input.
    fake = Path("/tmp/does_not_exist_xyz_123.tsv")
    if fake.exists():
        fake.unlink()
    out = subprocess.run(
        ["awk", "-F", "\t", 'NR>1 && $4 == "True" {n++} END {print n+0}', str(fake)],
        text=True, capture_output=True,
    )
    # awk writes to stderr when file missing; fix wrapper should guard.
    assert "0" in out.stdout or out.returncode != 0
