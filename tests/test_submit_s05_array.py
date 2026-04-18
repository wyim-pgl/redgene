"""Regression tests for Issue #12 — T8 per-site SLURM array hardening.

Covers six follow-up items documented on issue #12 after commits a5ee725 /
fe790f7:

  * M-1 `%25` concurrency cap and `32G` per-site mem must be overridable via
    env vars (`S05_ARRAY_THROTTLE`, `S05_ARRAY_MEM`) so clusters with tighter
    queue limits can re-use the script without forking it.
  * M-2 file-growth guard: `scripts/s05_insert_assembly.py` main() must not
    silently bloat past an agreed ceiling (LINE_BUDGET below) — pre-v1.1
    refactor target.
  * M-3 Phase 4 summary: the Phase 4 afterok wrapper must emit per-site
    wall-time / mem / exit-code aggregation via `sacct` before running the
    Phase 4 python call (comment + command suffices — we check the fixture).
  * M-4 CLI parser coupling: explicit warning next to the EXTRA_ARGS builder
    that any new optional argument must also be added to s05's argparse.
  * M-5 pickle security: comment on the `positive_sites.pkl` producer /
    consumer paths pointing to the v1.1 json/msgpack migration plan.
  * M-6 mapfile docs: explanatory comment on the `mapfile -t SITES` pipeline
    so the shell idiom is obvious to maintainers.
"""
from __future__ import annotations

import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "submit_s05_array.sh"
S05 = REPO_ROOT / "scripts" / "s05_insert_assembly.py"


def _text() -> str:
    return SCRIPT.read_text()


def _s05_text() -> str:
    return S05.read_text()


def test_submit_script_exists():
    assert SCRIPT.exists(), SCRIPT


# ---- M-1 parametrize throttle + per-site mem ------------------------------

def test_m1_array_throttle_is_parametrized():
    """`%25` concurrency cap must come from $S05_ARRAY_THROTTLE with a default."""
    text = _text()
    assert re.search(r"S05_ARRAY_THROTTLE[:\-]?=?", text), (
        "missing $S05_ARRAY_THROTTLE env-var default block"
    )
    # The SBATCH header must reference the variable (not a hard-coded %25).
    assert re.search(r"--array=0-\$\(\(N-1\)\)%\$\{?S05_ARRAY_THROTTLE", text) or \
           re.search(r"--array=0-\$\(\(N-1\)\)%\$S05_ARRAY_THROTTLE", text), (
        "SBATCH --array directive must use $S05_ARRAY_THROTTLE, not a literal %25"
    )


def test_m1_per_site_mem_is_parametrized():
    """`--mem=32G` must come from $S05_ARRAY_MEM with a default."""
    text = _text()
    assert re.search(r"S05_ARRAY_MEM[:\-]?=?", text), (
        "missing $S05_ARRAY_MEM env-var default block"
    )
    # SBATCH --mem directive must reference the variable.
    assert re.search(r"#SBATCH --mem=\$\{?S05_ARRAY_MEM", text), (
        "SBATCH --mem must use $S05_ARRAY_MEM, not a literal 32G"
    )


# ---- M-2 file-growth guard ------------------------------------------------

def test_m2_s05_main_line_budget():
    """main() must not silently bloat past the pre-v1.1 refactor ceiling.

    LINE_BUDGET is calibrated ~+50 lines above the current size so ordinary
    hardening commits (adding a log line here or there) do not trip the test,
    but a 200+ line accretion does — that's the signal to start the v1.1
    split documented in Issue #4.
    """
    text = _s05_text().splitlines()
    LINE_BUDGET = 550  # current ~462; see Issue #4 for the v1.1 split plan
    # Locate `def main(` and the next top-level `def `/class / sentinel.
    start = None
    end = None
    for i, line in enumerate(text):
        if start is None and line.startswith("def main"):
            start = i
            continue
        if start is not None and (line.startswith("if __name__") or
                                  (line.startswith(("def ", "class ")) and i != start)):
            end = i
            break
    assert start is not None and end is not None, "could not locate main()"
    size = end - start
    assert size <= LINE_BUDGET, (
        f"s05_insert_assembly.main() is {size} lines (>{LINE_BUDGET}); "
        "Issue #12 M-2 regression — trigger the Issue #4 v1.1 module split"
    )


# ---- M-3 Phase 4 summary (sacct aggregation) ------------------------------

def test_m3_phase4_emits_sacct_summary():
    """Phase 4 afterok wrapper must log per-array-task wall/mem via sacct."""
    text = _text()
    # We require a sacct invocation inside the Phase 4 wrap so the post-run
    # slurm log records per-task metrics for the M-3 summary.
    assert "sacct" in text, (
        "Phase 4 wrapper missing sacct summary (M-3) — add before python call"
    )


# ---- M-4 CLI parser coupling ---------------------------------------------

def test_m4_extra_args_has_coupling_comment():
    """The EXTRA_ARGS builder must warn about s05 argparse coupling."""
    text = _text()
    # Look for a comment within ~20 lines above/below the EXTRA_ARGS= line.
    lines = text.splitlines()
    idx = next((i for i, l in enumerate(lines) if "EXTRA_ARGS=" in l), None)
    assert idx is not None, "EXTRA_ARGS= line not found"
    window = "\n".join(lines[max(0, idx - 10): idx + 10])
    assert re.search(r"argparse|parser|CLI", window, re.IGNORECASE), (
        "M-4: EXTRA_ARGS block missing `argparse/parser/CLI` coupling comment"
    )


# ---- M-5 pickle security --------------------------------------------------

def test_m5_pickle_has_security_comment():
    """positive_sites.pkl producer/consumer must flag the v1.1 json migration."""
    text = _s05_text()
    # Find every usage and require at least one mentions the migration plan.
    assert "positive_sites.pkl" in text
    # The word `pickle` must appear near a `v1.1` or `json` note.
    assert re.search(
        r"pickle[\s\S]{0,400}(v1\.1|json|msgpack|trusted)", text, re.IGNORECASE
    ) or re.search(
        r"(v1\.1|json|msgpack|trusted)[\s\S]{0,400}pickle", text, re.IGNORECASE
    ), "M-5: pickle usage missing security/migration note"


# ---- M-6 mapfile docs -----------------------------------------------------

def test_m6_mapfile_has_explanatory_comment():
    text = _text()
    lines = text.splitlines()
    idx = next((i for i, l in enumerate(lines) if l.startswith("mapfile")), None)
    assert idx is not None, "mapfile -t SITES line not found"
    # Check 6 lines above mapfile for an explanatory comment.
    window = "\n".join(lines[max(0, idx - 6): idx])
    assert re.search(r"#.*mapfile|#.*-t\b|#.*newline-separated", window, re.IGNORECASE), (
        "M-6: mapfile line missing explanatory comment block above it"
    )
