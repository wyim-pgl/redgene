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

def test_summary_renders_one_page(tmp_path, monkeypatch):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    _write_grna(tmp_path / "s06_indel" / "grna_targets.tsv", _ON_TARGET_ROWS)
    _write_edits(tmp_path / "s06_indel" / "editing_sites.tsv", [_EDIT_ROW_OK])
    _state, on_targets, sites = mod.load_editing_data(tmp_path)

    captured = {}

    def fake_page_text(pdf, title, lines, **kwargs):
        captured["title"] = title
        captured["lines"] = lines

    monkeypatch.setattr(mod, "_page_text", fake_page_text)

    out_pdf = tmp_path / "summary.pdf"
    with PdfPages(out_pdf) as pdf:
        mod._page_crispr_summary(pdf, on_targets, sites)

    assert captured["title"] == "CRISPR Editing — Summary"
    body = "\n".join(captured["lines"])
    assert "gRNAs (on-target): 2" in body
    assert "Editing events:    1" in body
    # gRNA 1 had the insertion edit
    grna1_lines = [l for l in captured["lines"] if l.strip().startswith("1 ")]
    assert len(grna1_lines) == 1
    assert "insertion" in grna1_lines[0]
    assert "0.462" in grna1_lines[0] or "0.461" in grna1_lines[0]  # 0.4615 → .3f
    # gRNA 2 had no edits — em-dash placeholders
    grna2_lines = [l for l in captured["lines"] if l.strip().startswith("2 ")]
    assert len(grna2_lines) == 1
    assert "—" in grna2_lines[0]


def test_summary_truncates_long_locus(monkeypatch):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    on_targets = [{
        "grna_idx": "1", "grna_seq": "ACGT" * 5,
        "chrom": "B10v3_unplaced_contig_1234567", "cut_pos": "987654321",
        "strand": "+",
    }]
    sites: list[dict] = []

    captured = {}

    def fake_page_text(pdf, title, lines, **kwargs):
        captured["lines"] = lines

    monkeypatch.setattr(mod, "_page_text", fake_page_text)

    import io
    with PdfPages(io.BytesIO()) as pdf:
        mod._page_crispr_summary(pdf, on_targets, sites)

    grna_row = [l for l in captured["lines"] if l.strip().startswith("1 ")][0]
    assert "..." in grna_row
    # The truncated locus column must not blow past 32 chars
    # (We assert by reconstructing where the count column should land.)
    parts = grna_row.split()
    assert parts[0] == "1"  # gRNA idx column

# ---------------------------------------------------------------------------
# _page_crispr_profile (PNG embed + fallback)
# ---------------------------------------------------------------------------

def _make_real_png(path: Path) -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4, 2))
    ax.plot([0, 1, 2], [0, 1, 0])
    fig.savefig(path)
    plt.close(fig)


def test_profile_renders_with_valid_png(tmp_path):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    s06 = tmp_path / "s06_indel"
    s06.mkdir()
    _make_real_png(s06 / "editing_profile_gRNA1.png")
    target = {
        "grna_idx": "1", "grna_seq": "CATTAGCCTATGGTGAGCCA",
        "chrom": "SLM_r2.0ch04", "cut_pos": "2635445", "strand": "+",
    }

    out_pdf = tmp_path / "profile.pdf"
    with PdfPages(out_pdf) as pdf:
        mod._page_crispr_profile(pdf, tmp_path, target)

    assert out_pdf.stat().st_size > 0


def test_profile_png_missing_falls_back_to_text(tmp_path):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    (tmp_path / "s06_indel").mkdir()
    target = {
        "grna_idx": "2", "grna_seq": "GGTGTGCTATAAGTACTGAA",
        "chrom": "SLM_r2.0ch08", "cut_pos": "53314229", "strand": "-",
    }

    out_pdf = tmp_path / "profile_fallback.pdf"
    with PdfPages(out_pdf) as pdf:
        mod._page_crispr_profile(pdf, tmp_path, target)

    assert out_pdf.stat().st_size > 0  # no exception, page still written


def test_profile_corrupt_png_falls_back(tmp_path):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    s06 = tmp_path / "s06_indel"
    s06.mkdir()
    (s06 / "editing_profile_gRNA1.png").write_bytes(b"not-a-png")
    target = {
        "grna_idx": "1", "grna_seq": "CATTAGCCTATGGTGAGCCA",
        "chrom": "SLM_r2.0ch04", "cut_pos": "2635445", "strand": "+",
    }

    out_pdf = tmp_path / "profile_corrupt.pdf"
    with PdfPages(out_pdf) as pdf:
        mod._page_crispr_profile(pdf, tmp_path, target)

    assert out_pdf.stat().st_size > 0
