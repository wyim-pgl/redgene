"""Tests for Issue #6 Section 8 — CRISPR Editing panel.

Validates load_editing_data() state machine (5 state tests) and the
_page_crispr_summary text-page renderer. Tasks 3-4 will add tests for
the remaining page renderers.
"""
from __future__ import annotations

import importlib.util
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
# Fixture builders
# ---------------------------------------------------------------------------

_GRNA_HEADER = (
    "grna_idx\tgrna_seq\tchrom\tstart\tend\tstrand\tmismatches\thas_pam\t"
    "pam_seq\tsite_type\tcut_pos"
)
_EDIT_HEADER = (
    "chrom\tpos\tref\talt\ttype\tsize\tindel_seq\tqual\tdp\tcount\tfreq\t"
    "gt\tzygosity\tgrna_idx\tsite_type\tmismatches\tdistance_to_cut"
)


def _write_grna(path: Path, rows: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_GRNA_HEADER + "\n" + "\n".join(rows) + ("\n" if rows else ""))


def _write_edits(path: Path, rows: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_EDIT_HEADER + "\n" + "\n".join(rows) + ("\n" if rows else ""))


# Two on-target rows + one off-target row, matching tomato_Cas9_A2_3 shape
_ON_TARGET_ROWS = [
    "1\tCATTAGCCTATGGTGAGCCA\tSLM_r2.0ch04\t2635428\t2635448\t+\t0\tTrue\tTGG\ton-target\t2635445",
    "2\tGGTGTGCTATAAGTACTGAA\tSLM_r2.0ch08\t53314226\t53314246\t-\t0\tTrue\tCCG\ton-target\t53314229",
]
_OFF_TARGET_ROW = (
    "1\tCATTAGCCTATGGTGAGCCA\tSLM_r2.0ch11\t54104526\t54104544\t-\t3\tFalse\tTTA\toff-target\t54104529"
)
_EDIT_ROW_OK = (
    "SLM_r2.0ch04\t2635445\tG\tGT\tinsertion\t1\tT\t60.0\t13\t6\t0.4615\t0/1\theterozygous"
    "\t1\ton-target\t0\t0"
)


# ---------------------------------------------------------------------------
# load_editing_data state machine
# ---------------------------------------------------------------------------

def test_state_missing_no_s06_dir(tmp_path):
    mod = _load_module()
    state, on_targets, sites = mod.load_editing_data(tmp_path)
    assert state == "missing"
    assert on_targets is None
    assert sites is None


def test_state_missing_no_targets(tmp_path):
    mod = _load_module()
    (tmp_path / "s06_indel").mkdir()
    state, on_targets, sites = mod.load_editing_data(tmp_path)
    assert state == "missing"
    assert on_targets is None
    assert sites is None


def test_state_missing_no_on_targets(tmp_path):
    mod = _load_module()
    _write_grna(tmp_path / "s06_indel" / "grna_targets.tsv", [_OFF_TARGET_ROW])
    state, on_targets, sites = mod.load_editing_data(tmp_path)
    assert state == "missing"
    assert on_targets is None
    assert sites is None


def test_state_no_edits(tmp_path):
    mod = _load_module()
    _write_grna(tmp_path / "s06_indel" / "grna_targets.tsv",
                _ON_TARGET_ROWS + [_OFF_TARGET_ROW])
    _write_edits(tmp_path / "s06_indel" / "editing_sites.tsv", [])
    state, on_targets, sites = mod.load_editing_data(tmp_path)
    assert state == "no_edits"
    assert on_targets is not None
    assert len(on_targets) == 2
    assert sites == []


def test_state_ok(tmp_path):
    mod = _load_module()
    _write_grna(tmp_path / "s06_indel" / "grna_targets.tsv",
                _ON_TARGET_ROWS + [_OFF_TARGET_ROW])
    _write_edits(tmp_path / "s06_indel" / "editing_sites.tsv", [_EDIT_ROW_OK])
    state, on_targets, sites = mod.load_editing_data(tmp_path)
    assert state == "ok"
    assert on_targets is not None
    assert len(on_targets) == 2
    assert sites is not None
    assert len(sites) == 1
    assert sites[0]["grna_idx"] == "1"
    assert sites[0]["type"] == "insertion"


# ---------------------------------------------------------------------------
# _page_crispr_summary
# ---------------------------------------------------------------------------

def test_summary_renders_one_page(tmp_path):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    _write_grna(tmp_path / "s06_indel" / "grna_targets.tsv", _ON_TARGET_ROWS)
    _write_edits(tmp_path / "s06_indel" / "editing_sites.tsv", [_EDIT_ROW_OK])
    _state, on_targets, sites = mod.load_editing_data(tmp_path)

    out_pdf = tmp_path / "summary.pdf"
    with PdfPages(out_pdf) as pdf:
        mod._page_crispr_summary(pdf, on_targets, sites)

    assert out_pdf.exists()
    assert out_pdf.stat().st_size > 0
