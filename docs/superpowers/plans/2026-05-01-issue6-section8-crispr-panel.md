# Issue #6 Section 8 — CRISPR Editing Panel Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the `_page_placeholder` call at `scripts/reports/insertion_pdf.py:599` with a real CRISPR Editing panel that embeds per-gRNA editing-profile PNGs already produced by `scripts/viz/plot_editing_profile.py`.

**Architecture:** A 3-state dispatcher (`missing` / `no_edits` / `ok`) selects between an informational page and a multi-page block (1 summary + 1 page per on-target gRNA). Each gRNA page embeds the existing `editing_profile_gRNA{idx}.png` via `matplotlib.image.imread` + `ax.imshow`. No new dependencies; PDF reads existing files only (no subprocess to viz scripts).

**Tech Stack:** Python 3.11, matplotlib (PdfPages, image), pytest. No new deps.

**Spec:** [`docs/superpowers/specs/2026-05-01-issue6-section8-crispr-panel-design.md`](../specs/2026-05-01-issue6-section8-crispr-panel-design.md) (commit `c3be1f5`).

**File Structure:**
- Modify: `scripts/reports/insertion_pdf.py` — add 4 functions, swap 1 line in `generate_pdf`, add 1 import
- Create: `tests/test_insertion_pdf_section8.py` — 9 tests for the new state machine + page rendering
- No other files touched. The 242-test baseline (235 pre-existing + 7 from Issue #4) must remain green.

---

## Task 1: Add `load_editing_data` state machine

**Files:**
- Create: `tests/test_insertion_pdf_section8.py`
- Modify: `scripts/reports/insertion_pdf.py` (append new function near other `load_*` helpers, around line 220)

- [ ] **Step 1.1: Create test file with module loader and 5 state-machine tests**

Write `tests/test_insertion_pdf_section8.py`:

```python
"""Tests for Issue #6 Section 8 — CRISPR Editing panel.

Validates load_editing_data() state machine and the four new page renderers
(_page_crispr_editing dispatcher, _page_crispr_summary, _page_crispr_profile).
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


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
```

- [ ] **Step 1.2: Run new tests and confirm they fail with `AttributeError`**

Run:
```bash
cd /data/gpfs/assoc/pgl/develop/redgene
pytest tests/test_insertion_pdf_section8.py -v
```
Expected: 5 FAIL, all with `AttributeError: module 'insertion_pdf' has no attribute 'load_editing_data'`.

- [ ] **Step 1.3: Implement `load_editing_data` in `scripts/reports/insertion_pdf.py`**

Insert immediately after the existing `load_coc_log` function (around line 220, before `build_sample_info`):

```python
def load_editing_data(
    sample_dir: Path,
) -> tuple[str, list[dict[str, Any]] | None, list[dict[str, Any]] | None]:
    """Load step-6 CRISPR outputs and return a 3-state tuple.

    State machine:
      "missing"  → s06_indel/ absent, grna_targets.tsv absent, no on-target rows,
                   or any TSV parse failure (treat as not-applicable)
      "no_edits" → on-target gRNAs present, editing_sites.tsv missing or empty
      "ok"       → on-target gRNAs present and at least one edit row
    """
    s06 = Path(sample_dir) / "s06_indel"
    targets_path = s06 / "grna_targets.tsv"
    sites_path = s06 / "editing_sites.tsv"

    if not targets_path.exists():
        return ("missing", None, None)

    try:
        with open(targets_path, newline="") as fh:
            targets = list(csv.DictReader(fh, delimiter="\t"))
    except (OSError, csv.Error):
        return ("missing", None, None)

    on_targets = [t for t in targets if t.get("site_type") == "on-target"]
    if not on_targets:
        return ("missing", None, None)

    sites: list[dict[str, Any]] = []
    if sites_path.exists():
        try:
            with open(sites_path, newline="") as fh:
                sites = list(csv.DictReader(fh, delimiter="\t"))
        except (OSError, csv.Error):
            sites = []

    if not sites:
        return ("no_edits", on_targets, [])
    return ("ok", on_targets, sites)
```

- [ ] **Step 1.4: Run the 5 state-machine tests and confirm PASS**

Run:
```bash
pytest tests/test_insertion_pdf_section8.py -v -k state_
```
Expected: 5 PASS.

- [ ] **Step 1.5: Commit**

```bash
git add tests/test_insertion_pdf_section8.py scripts/reports/insertion_pdf.py
git commit -m "$(cat <<'EOF'
Issue #6 Section 8 [1/4]: add load_editing_data state machine

Three-state classifier (missing / no_edits / ok) for the upcoming
CRISPR Editing panel. Soft failure: any IO or csv.Error is treated as
"missing" so the PDF still renders.
EOF
)"
```

---

## Task 2: Add `_page_crispr_summary` text page

**Files:**
- Modify: `scripts/reports/insertion_pdf.py` (append new function near `_page_placeholder`, around line 510)
- Modify: `tests/test_insertion_pdf_section8.py` (append summary-render test)

- [ ] **Step 2.1: Add summary rendering test**

Append to `tests/test_insertion_pdf_section8.py`:

```python
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
```

- [ ] **Step 2.2: Run and confirm fail**

Run:
```bash
pytest tests/test_insertion_pdf_section8.py::test_summary_renders_one_page -v
```
Expected: FAIL with `AttributeError: module 'insertion_pdf' has no attribute '_page_crispr_summary'`.

- [ ] **Step 2.3: Implement `_page_crispr_summary`**

Insert in `scripts/reports/insertion_pdf.py` immediately after `_page_placeholder` (around line 511):

```python
def _page_crispr_summary(
    pdf: PdfPages,
    on_targets: list[dict[str, Any]],
    sites: list[dict[str, Any]],
) -> None:
    """Section 8 summary page: gRNA inventory + per-gRNA edit highlights."""
    n_grna = len(on_targets)
    n_edits = len(sites)
    lines = [
        f"gRNAs (on-target): {n_grna}",
        f"Editing events:    {n_edits}",
        "",
        f"  {'gRNA':<5} {'Locus':<32} {'Edits':>6} {'Top freq':>10} {'Top type':<12}",
        f"  {'-'*5} {'-'*32} {'-'*6} {'-'*10} {'-'*12}",
    ]
    for t in on_targets:
        idx = t.get("grna_idx", "?")
        chrom = t.get("chrom", "?")
        cut = t.get("cut_pos", "?")
        strand = t.get("strand", "?")
        locus = f"{chrom}:{cut} ({strand})"
        edits = [s for s in sites if s.get("grna_idx") == idx]
        if edits:
            top = max(edits, key=lambda r: float(r.get("freq", "0") or 0))
            try:
                freq_str = f"{float(top.get('freq', '0') or 0):.3f}"
            except ValueError:
                freq_str = "—"
            top_type = top.get("type", "—")
        else:
            freq_str = "—"
            top_type = "—"
        lines.append(
            f"  {str(idx):<5} {locus:<32} {len(edits):>6d} {freq_str:>10} {top_type:<12}"
        )
    _page_text(pdf, "CRISPR Editing — Summary", lines)
```

- [ ] **Step 2.4: Run and confirm PASS**

Run:
```bash
pytest tests/test_insertion_pdf_section8.py::test_summary_renders_one_page -v
```
Expected: PASS.

- [ ] **Step 2.5: Commit**

```bash
git add tests/test_insertion_pdf_section8.py scripts/reports/insertion_pdf.py
git commit -m "$(cat <<'EOF'
Issue #6 Section 8 [2/4]: add _page_crispr_summary text page

Single-page summary: gRNA on-target count, total edit count, and a
per-gRNA table with top-frequency highlight. Used by the dispatcher
when load_editing_data returns state="ok".
EOF
)"
```

---

## Task 3: Add `_page_crispr_profile` PNG embed + fallback

**Files:**
- Modify: `scripts/reports/insertion_pdf.py` (add `import matplotlib.image as mpimg` near other matplotlib imports, append new function after `_page_crispr_summary`)
- Modify: `tests/test_insertion_pdf_section8.py` (append PNG-missing and corrupt-PNG tests)

- [ ] **Step 3.1: Add 2 fallback tests**

Append to `tests/test_insertion_pdf_section8.py`:

```python
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
```

- [ ] **Step 3.2: Run and confirm fail**

Run:
```bash
pytest tests/test_insertion_pdf_section8.py::test_profile_renders_with_valid_png \
       tests/test_insertion_pdf_section8.py::test_profile_png_missing_falls_back_to_text \
       tests/test_insertion_pdf_section8.py::test_profile_corrupt_png_falls_back -v
```
Expected: 3 FAIL with `AttributeError: ... '_page_crispr_profile'`.

- [ ] **Step 3.3: Add `mpimg` import and implement `_page_crispr_profile`**

In `scripts/reports/insertion_pdf.py`, add the import on the line after `import matplotlib.pyplot as plt` (around line 56):

```python
import matplotlib.image as mpimg
```

Then append immediately after `_page_crispr_summary`:

```python
def _page_crispr_profile(
    pdf: PdfPages,
    sample_dir: Path,
    target: dict[str, Any],
) -> None:
    """Render one PDF page embedding editing_profile_gRNA{idx}.png.

    Falls back to a text page if the PNG is missing or unreadable.
    """
    idx = str(target.get("grna_idx", "?"))
    seq = target.get("grna_seq", "?")
    chrom = target.get("chrom", "?")
    cut = target.get("cut_pos", "?")
    title = f"CRISPR Editing — gRNA {idx}"
    header = f"gRNA {idx} — {chrom}:{cut} ({seq})"

    png_path = Path(sample_dir) / "s06_indel" / f"editing_profile_gRNA{idx}.png"
    if not png_path.exists():
        _page_text(pdf, title, [
            header, "",
            f"(editing_profile_gRNA{idx}.png missing — "
            "scripts/viz/plot_editing_profile.py needs to run first)",
        ])
        return

    try:
        img = mpimg.imread(str(png_path))
    except (OSError, ValueError, SyntaxError):
        _page_text(pdf, title, [
            header, "",
            f"(editing_profile_gRNA{idx}.png unreadable — file may be corrupt)",
        ])
        return

    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis("off")
    ax.text(0.05, 0.96, title, fontsize=18, fontweight="bold",
            transform=ax.transAxes, verticalalignment="top")
    ax.text(0.05, 0.92, header, fontsize=10, family="monospace",
            transform=ax.transAxes, verticalalignment="top")
    img_ax = fig.add_axes((0.05, 0.05, 0.90, 0.83))
    img_ax.imshow(img)
    img_ax.axis("off")
    pdf.savefig(fig)
    plt.close(fig)
```

- [ ] **Step 3.4: Run and confirm PASS**

Run:
```bash
pytest tests/test_insertion_pdf_section8.py -v
```
Expected: 8 PASS (5 state + 1 summary + 3 profile tests).

- [ ] **Step 3.5: Commit**

```bash
git add tests/test_insertion_pdf_section8.py scripts/reports/insertion_pdf.py
git commit -m "$(cat <<'EOF'
Issue #6 Section 8 [3/4]: add _page_crispr_profile PNG embed

Per-gRNA page that embeds editing_profile_gRNA{idx}.png via
matplotlib.image.imread + imshow. Falls back to a text page on
missing or unreadable PNG so PDF generation never crashes.
EOF
)"
```

---

## Task 4: Wire `_page_crispr_editing` dispatcher into `generate_pdf`

**Files:**
- Modify: `scripts/reports/insertion_pdf.py:599` (swap `_page_placeholder` call) and append `_page_crispr_editing` after `_page_crispr_profile`
- Modify: `tests/test_insertion_pdf_section8.py` (append 2 dispatcher tests)

- [ ] **Step 4.1: Add 2 dispatcher tests**

Append to `tests/test_insertion_pdf_section8.py`:

```python
# ---------------------------------------------------------------------------
# _page_crispr_editing dispatcher
# ---------------------------------------------------------------------------

def test_dispatcher_missing_returns_one_page(tmp_path):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    out_pdf = tmp_path / "missing.pdf"
    with PdfPages(out_pdf) as pdf:
        n = mod._page_crispr_editing(pdf, tmp_path)

    assert n == 1
    assert out_pdf.stat().st_size > 0


def test_dispatcher_ok_returns_summary_plus_n_profiles(tmp_path):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    _write_grna(tmp_path / "s06_indel" / "grna_targets.tsv",
                _ON_TARGET_ROWS + [_OFF_TARGET_ROW])
    _write_edits(tmp_path / "s06_indel" / "editing_sites.tsv", [_EDIT_ROW_OK])
    _make_real_png(tmp_path / "s06_indel" / "editing_profile_gRNA1.png")
    _make_real_png(tmp_path / "s06_indel" / "editing_profile_gRNA2.png")

    out_pdf = tmp_path / "ok.pdf"
    with PdfPages(out_pdf) as pdf:
        n = mod._page_crispr_editing(pdf, tmp_path)

    # 1 summary + 2 on-target gRNA profile pages
    assert n == 3
```

- [ ] **Step 4.2: Run and confirm fail**

Run:
```bash
pytest tests/test_insertion_pdf_section8.py::test_dispatcher_missing_returns_one_page \
       tests/test_insertion_pdf_section8.py::test_dispatcher_ok_returns_summary_plus_n_profiles -v
```
Expected: 2 FAIL with `AttributeError: ... '_page_crispr_editing'`.

- [ ] **Step 4.3: Implement `_page_crispr_editing` dispatcher**

Append in `scripts/reports/insertion_pdf.py` immediately after `_page_crispr_profile`:

```python
def _page_crispr_editing(pdf: PdfPages, sample_dir: Path) -> int:
    """Section 8 dispatcher. Returns the number of pages written (>=1)."""
    state, on_targets, sites = load_editing_data(sample_dir)

    if state == "missing":
        _page_text(pdf, "CRISPR Editing", [
            "(s06_indel CRISPR analysis was not run for this sample,",
            " or no on-target gRNA was defined.)",
        ])
        return 1

    if state == "no_edits":
        n = len(on_targets) if on_targets else 0
        _page_text(pdf, "CRISPR Editing", [
            f"On-target gRNAs analyzed: {n}",
            "Editing events detected:  0",
            "",
            "(no editing events passed step-6 thresholds)",
        ])
        return 1

    # state == "ok"
    pages = 0
    _page_crispr_summary(pdf, on_targets or [], sites or [])
    pages += 1
    for target in on_targets or []:
        _page_crispr_profile(pdf, sample_dir, target)
        pages += 1
    return pages
```

- [ ] **Step 4.4: Swap the placeholder call in `generate_pdf`**

In `scripts/reports/insertion_pdf.py`, locate the only `_page_placeholder` call in `generate_pdf()` (line numbers will have shifted from the original 599 after Tasks 1-3 added new code; the call is unique). Replace this exact two-line block:

```python
        _page_placeholder(pdf, "CRISPR Editing",
                          "editing panel pending v1.1"); pages += 1     # (8)
```

with:

```python
        pages += _page_crispr_editing(pdf, sample_dir)                  # (8)
```

Verify nothing else still references `_page_placeholder` (the helper itself stays defined for the future appendix-style placeholders, but no caller should remain in `generate_pdf`):

```bash
grep -n "_page_placeholder" scripts/reports/insertion_pdf.py
```
Expected: only the function definition (one line), no call sites.

- [ ] **Step 4.5: Run new tests + full suite to confirm no regression**

Run:
```bash
pytest tests/test_insertion_pdf_section8.py -v
```
Expected: 10 PASS (5 state + 1 summary + 3 profile + 2 dispatcher).

Run:
```bash
pytest tests/ -q
```
Expected: 252 PASS + 1 skipped (242 baseline + 10 new). Zero failures.

- [ ] **Step 4.6: Commit**

```bash
git add tests/test_insertion_pdf_section8.py scripts/reports/insertion_pdf.py
git commit -m "$(cat <<'EOF'
Issue #6 Section 8 [4/4]: wire dispatcher into generate_pdf

Replaces the _page_placeholder call at line 599 with the new
_page_crispr_editing dispatcher. Section 8 now writes 1 informational
page when CRISPR analysis is N/A, 1 summary page when there are no
edits, or 1 summary + N gRNA-profile pages when state="ok".

Test count grows from 242 to 252 PASS (+10 new in
tests/test_insertion_pdf_section8.py); zero regressions.
EOF
)"
```

---

## Task 5: 3-sample acceptance validation

**Files:**
- No code changes. Run the pipeline against three real samples and inspect the resulting PDFs.

- [ ] **Step 5.1: Pre-generate editing-profile PNGs for tomato_Cas9_A2_3**

Run:
```bash
cd /data/gpfs/assoc/pgl/develop/redgene
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

python scripts/viz/plot_editing_profile.py \
  --treatment-bam results/tomato_Cas9_A2_3/s04_host_map/tomato_Cas9_A2_3_host.bam \
  --wt-bam results/tomato_WT/s04_host_map/tomato_WT_host.bam \
  --host-ref db/SLM_r2.0.pmol.fasta \
  --grna-targets results/tomato_Cas9_A2_3/s06_indel/grna_targets.tsv \
  --editing-sites results/tomato_Cas9_A2_3/s06_indel/editing_sites.tsv \
  --sample-name tomato_Cas9_A2_3 \
  --outdir results
```

Expected: `results/tomato_Cas9_A2_3/s06_indel/editing_profile_gRNA1.png` and `editing_profile_gRNA2.png` exist (≥10 KB each).

Verify:
```bash
ls -la results/tomato_Cas9_A2_3/s06_indel/editing_profile_gRNA*.png
```

- [ ] **Step 5.2: Generate PDF for tomato_Cas9_A2_3 (state="ok")**

Run:
```bash
python run_pipeline.py --sample tomato_Cas9_A2_3 --steps 8
```

Inspect:
```bash
pdfinfo results/tomato_Cas9_A2_3/tomato_Cas9_A2_3_insertion_report.pdf | grep -i pages
```

Expected: page count is 2 higher than the pre-fix baseline (the placeholder was 1 page; new output is 1 summary + 2 gRNA profile pages = 3 pages; net +2).

Open the PDF (`evince` / `okular` / scp to local machine) and confirm:
- Section 8 summary page lists 2 on-target gRNAs (chr04, chr08) and 1 detected edit (chr04:2635445, insertion, freq ≈ 0.46).
- gRNA 1 page embeds the chr04 editing-profile PNG (visible nucleotide quilt + indel highlighting).
- gRNA 2 page embeds the chr08 editing-profile PNG.

- [ ] **Step 5.3: Generate PDF for rice_G281 (state="missing")**

Run:
```bash
python run_pipeline.py --sample rice_G281 --steps 8
```

Inspect:
```bash
pdfinfo results/rice_G281/rice_G281_insertion_report.pdf | grep -i pages
```

Expected: page count unchanged from pre-fix baseline (Section 8 was 1 placeholder page → 1 informational page).

Open and confirm Section 8 reads "(s06_indel CRISPR analysis was not run for this sample, or no on-target gRNA was defined.)".

- [ ] **Step 5.4: Generate PDF for cucumber_line225 (state="missing")**

Run:
```bash
python run_pipeline.py --sample cucumber_line225 --steps 8
```

Inspect:
```bash
pdfinfo results/cucumber_line225/cucumber_line225_insertion_report.pdf | grep -i pages
```

Expected: page count unchanged from pre-fix baseline. Section 8 informational page identical to rice_G281.

- [ ] **Step 5.5: Comment on Issue #6 with PDF artifacts**

Use `gh` (note: `gh` is at `~/micromamba/bin/gh`, NOT in the redgene env per CLAUDE.md):

```bash
~/micromamba/bin/gh issue comment 6 --body "$(cat <<'EOF'
Section 8 (CRISPR Editing) panel implemented in commits [1-4]/[1-4]:
- tomato_Cas9_A2_3: state="ok", +2 pages (1 summary + 2 gRNA profile pages with embedded editing-profile PNGs)
- rice_G281: state="missing", 1 informational page (no s06_indel/grna_targets.tsv)
- cucumber_line225: state="missing", 1 informational page

Test suite: 252 PASS + 1 skipped (was 242 + 1; +10 new in tests/test_insertion_pdf_section8.py).
EOF
)"
```

- [ ] **Step 5.6: Close Issue #6 if all 3 PDFs check out visually**

Only after the user confirms the PDFs look correct:

```bash
~/micromamba/bin/gh issue close 6
```

---

## Spec Coverage Matrix

| Spec section | Implemented in |
|--------------|---------------|
| `load_editing_data` 3-state | Task 1 |
| `_page_crispr_summary` | Task 2 |
| `_page_crispr_profile` (PNG embed + fallback) | Task 3 |
| `_page_crispr_editing` dispatcher | Task 4 |
| `generate_pdf` swap | Task 4.4 |
| 5 state-machine tests | Task 1.1 |
| 1 summary test | Task 2.1 |
| PNG-missing + corrupt-PNG fallback tests | Task 3.1 |
| Dispatcher page-count tests | Task 4.1 |
| 3-sample acceptance | Task 5 |
| Page-count invariant (always ≥1) | Task 4.3 dispatcher; tested in 4.1 |

## Notes for the Implementer

- **Read the spec first.** Path: `docs/superpowers/specs/2026-05-01-issue6-section8-crispr-panel-design.md`. The state machine and page contract are the most important parts.
- **Verbatim style.** This codebase prefers no docstrings/comments unless WHY is non-obvious. Keep new code spare; the spec already documents intent.
- **No subprocess to viz scripts.** PDF generation reads existing files. PNGs must be present already (Task 5.1 generates them for the validation sample).
- **Do not regenerate verdict snapshots.** Section 8 changes do not touch verdict computation; the 26 snapshot tests should stay untouched.
- **Test runtime.** Each new test creates a small tmp PDF (~50 KB). Full suite runtime should stay near the 20s baseline.
- **Commit after every task.** Each task ships independently and leaves the repo in a green state.
