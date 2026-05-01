# Issue #6 Section 8 — CRISPR Editing Panel Design

**Date**: 2026-05-01
**Status**: Spec — pending implementation plan
**GitHub issue**: [wyim-pgl/redgene#6](https://github.com/wyim-pgl/redgene/issues/6)
**Implements**: R-5 (PDF-ready insertion report)

## Goal

Replace the `_page_placeholder` call at `scripts/reports/insertion_pdf.py:599` (Section 8 of `generate_pdf()`) with a real CRISPR Editing panel that embeds the per-gRNA editing-profile PNGs already produced by `scripts/viz/plot_editing_profile.py`.

## Non-Goals

- Generating the editing-profile PNGs themselves. PDF generation reads existing files only — consistent with the step-8 contract documented in `CLAUDE.md` ("opt-in; reads existing s05/s07 output").
- Re-running step 6 indel detection. Inputs are `s06_indel/grna_targets.tsv` and `s06_indel/editing_sites.tsv`.
- BAM-level pileup rendering inside the PDF. The wide-format CRISPResso2-style figure stays in `scripts/viz/plot_editing_profile.py`; we embed its raster output.
- Off-target editing visualization. Only `site_type == "on-target"` rows drive page count.

## Architecture

Three new functions are added to `scripts/reports/insertion_pdf.py`. No other module is touched.

```
generate_pdf()
  └─ _page_crispr_editing(pdf, sample_dir)        ← replaces _page_placeholder
       ├─ load_editing_data(sample_dir)
       │     returns (state, on_target_rows, edit_rows)
       │     state ∈ {"missing", "no_edits", "ok"}
       │
       ├─ state == "missing":  _page_text(...)                        # 1 page
       ├─ state == "no_edits": _page_text(... gRNA N / 0 edits)       # 1 page
       └─ state == "ok":
             _page_crispr_summary(pdf, on_targets, sites)              # 1 page
             for target in on_targets:
                 _page_crispr_profile(pdf, sample_dir, target)         # 1 page each
```

### Invariants

- Section 8 always writes ≥1 page so downstream sections (audit, CoC) keep their page numbers.
- All PDF figures use `figsize=(8.5, 11)` to match existing pages.
- PNG embedding uses `matplotlib.image.imread` + `ax.imshow` (no Pillow / pypdf dependency added).
- `editing_sites.tsv` schema is consumed as-emitted by `s06_indel.py` (chrom / pos / ref / alt / type / size / freq / zygosity / grna_idx / site_type / mismatches / distance_to_cut). No re-parsing into typed records.

## Components

### `load_editing_data(sample_dir: Path) -> tuple[str, list[dict]|None, list[dict]|None]`

Reads `s06_indel/grna_targets.tsv` and `s06_indel/editing_sites.tsv`. Returns one of three states:

| State | Inputs | Meaning |
|-------|--------|---------|
| `"missing"` | `s06_indel/` absent OR `grna_targets.tsv` absent OR no on-target rows | Step 6 was not run, or no gRNAs were defined for this sample |
| `"no_edits"` | on-target rows present, `editing_sites.tsv` absent or empty | gRNAs analyzed, zero edits detected |
| `"ok"` | on-target rows + ≥1 edit row | Both available |

TSV parsing follows the soft-failure pattern of `load_audit_header`, `load_copynumber`, `load_coc_log`: catch `OSError` and `csv.Error`, return `("missing", None, None)`.

### `_page_crispr_editing(pdf, sample_dir) -> int`

Section dispatcher. Returns the number of pages written.

- `state == "missing"` → 1 page via `_page_text("CRISPR Editing", [explanation])`
- `state == "no_edits"` → 1 page via `_page_text` listing `len(on_targets)` and the 0-edit verdict
- `state == "ok"` → calls `_page_crispr_summary` (1 page) then `_page_crispr_profile` once per on-target gRNA

### `_page_crispr_summary(pdf, on_targets, sites)`

Single text page (uses `_page_text`) listing:

- Counts: `gRNAs (on-target): N`, `Editing events: M`
- Table (monospace) per gRNA: index, locus (`chrom:cut_pos (strand)`), edit count, top-frequency value (`max(freq)` across that gRNA's `editing_sites` rows), and the type column of that top-frequency row. Rows with zero edits show `—` placeholders.

### `_page_crispr_profile(pdf, sample_dir, target)`

Single page per on-target gRNA.

- Header at top (0.92 → 0.98 in axes coords): `gRNA {idx} — {chrom}:{cut_pos} ({grna_seq})`
- Body (0.05 → 0.90): `ax.imshow(mpimg.imread(png_path))`, `ax.axis("off")`
- PNG path: `sample_dir / "s06_indel" / f"editing_profile_gRNA{idx}.png"`
- PNG missing or unreadable: fall back to `_page_text(title, ["(editing_profile_gRNA{idx}.png missing — scripts/viz/plot_editing_profile.py needs to run first)"])`

## Data Flow

```
sample_dir/
├── s06_indel/
│   ├── grna_targets.tsv        ← required for state != "missing"
│   ├── editing_sites.tsv       ← required for state == "ok"
│   ├── editing_profile_gRNA1.png   ← embedded by _page_crispr_profile
│   └── editing_profile_gRNA2.png
└── ...

generate_pdf()
  └─ _page_crispr_editing(pdf, sample_dir)
       └─ load_editing_data → reads two TSVs
       └─ _page_crispr_profile → reads one PNG per on-target gRNA
```

No writes to `sample_dir`. No subprocess calls.

## Error Handling

| Case | Behavior |
|------|----------|
| `s06_indel/` directory absent | `state = "missing"`, 1 informational page |
| `grna_targets.tsv` absent | `state = "missing"` |
| `grna_targets.tsv` has only off-target rows | `state = "missing"` (no on-target gRNA = treat as not-applicable) |
| TSV parse failure (`OSError`, `csv.Error`) | `state = "missing"` |
| `editing_sites.tsv` absent or 0 rows | `state = "no_edits"` |
| PNG file absent for a given gRNA | per-page text fallback; other gRNA pages and the summary still render |
| PNG file unreadable (`OSError`, `ValueError` from `imread`) | per-page text fallback |

## Testing

New file: `tests/test_insertion_pdf_section8.py`. The 242-test baseline must remain green.

| Test | Fixture | Assertion |
|------|---------|-----------|
| `test_state_missing_no_s06_dir` | empty `sample_dir` | `state == "missing"` |
| `test_state_missing_no_targets` | empty `s06_indel/` | `state == "missing"` |
| `test_state_missing_no_on_targets` | `grna_targets.tsv` with only off-target rows | `state == "missing"` |
| `test_state_no_edits` | on-target gRNA + empty `editing_sites.tsv` | `state == "no_edits"` |
| `test_state_ok` | tomato-shaped fixture (2 on-target gRNAs, 1 edit) | `state == "ok"`; lengths match |
| `test_pdf_section8_missing_one_page` | empty fixture | `_page_crispr_editing` returns 1 |
| `test_pdf_section8_ok_n_plus_one_pages` | 2-gRNA fixture | returns 3 (1 summary + 2 profiles) |
| `test_pdf_section8_png_missing_fallback` | state=ok but no PNG files | returns 3, no exceptions |
| `test_pdf_section8_corrupt_png` | state=ok with malformed PNG bytes | returns 3, no exceptions |

The existing snapshot suite (`tests/test_verdict_snapshots.py`, 26 PASSED) is unaffected — Section 8 changes do not touch verdict computation.

## Acceptance Criteria

After implementation, three real-sample validation runs are required (per Issue #6):

1. **tomato_Cas9_A2_3** — state must be `"ok"`. Pre-step: run `scripts/viz/plot_editing_profile.py` to produce 2 PNGs. Then `python run_pipeline.py --sample tomato_Cas9_A2_3 --steps 8`. Resulting PDF must show:
   - 1 summary page listing the 2 gRNAs (chr04 / chr08) and the 1 detected edit (chr04:2635445, 1 bp insertion, freq ≈ 0.46)
   - 2 gRNA-profile pages with PNGs embedded
2. **rice_G281** — state must be `"missing"` (non-CRISPR sample, no `s06_indel/grna_targets.tsv`). 1 informational page in Section 8.
3. **cucumber_line225** — state must be `"missing"`. 1 informational page in Section 8.

Verification: `pdfinfo` page counts must increase relative to the pre-fix baseline by exactly the expected delta (tomato: +2, rice/cucumber: +0).

## Out of Scope (Future Issues)

- Element-annotation appendix for samples with GFP/bar/nptII hits — tracked separately under Issue #6 follow-up.
- Off-target editing review pages — would need a different layout and new TSV joins; not part of R-5.
- KCGP nomenclature labels (Issue #7) — orthogonal to Section 8.

## Open Questions

None at spec time. Implementation may surface PNG aspect-ratio details (the editing-profile figure is wide); if `imshow` produces excessive whitespace or distortion at A4-portrait, the implementer may rotate the figure or letterbox it without changing the spec.
