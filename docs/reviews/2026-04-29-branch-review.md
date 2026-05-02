# Code Review — feature branch changes since `6828421`

**Scope:** Commits `46aaf40..39a123f` (7 commits): Issue #4 Session 1 (primitives / filter_b_flank / filter_c_chimeric / DAG test), Issue #6 PDF report, .gitignore + cleanup.
**Reviewer:** Claude (critical pass, no edits).
**Baseline:** `main @ 39a123f`, pytest 235 PASS + 1 skipped.

---

## Summary verdict

Overall the branch is **solid, production-acceptable**. The s05 module split executes the design spec cleanly, the PDF report delivers 3 validated sample outputs, and the DAG guard is a forward-compatible piece of test infrastructure. Issues below are **polish / consistency**, not blockers.

- **Architecture: A-** — DAG design is clean and forward-compatible; one design-spec claim (`_apply_canonical_override` "moved to verdict.py") is **not yet satisfied**; the function still lives in the monolith.
- **Code Quality: B+** — mostly tight; style drift between filter modules (`Optional` vs `| None`); several silent-failure paths.
- **Maintainability: A** — every new module docstring links back to the design spec; commit messages consistent and specific.
- **UX: B** — step 8 wired correctly but **not documented in `--help`**; PDF silently drops overflow content.

---

## 1. Architecture & Design

### 1.1 ✅ DAG enforcement is excellent
`tests/test_s05_import_dag.py` uses an explicit `STAGE` dict pre-populated with all 14 planned modules. When Session 4 adds `assembly.py` (stage 3), the test automatically validates it with **zero test edits required**. This is the kind of forward-compatible infrastructure that pays dividends across a 7-session refactor.

The edge case (`config_loader` at stage 7 > `verdict` at 6) is handled correctly: the test allows later→earlier imports (config_loader depends on VerdictRules, which is the intended direction).

### 1.2 ⚠️ Design-spec claim vs. code reality (medium)
**Finding:** `docs/superpowers/specs/2026-04-19-s05-module-split-design.md` §2.3 lists `_apply_canonical_override` as **"moved from monolith line 127"** — implying it now lives in `scripts/s05/verdict.py`.

**Reality:** `_apply_canonical_override` is still defined in `scripts/s05_insert_assembly.py:161` and **not present** in `scripts/s05/verdict.py` (which only exports `compute_verdict`, `FilterEvidence`, `VerdictRules`). The monolith's `compute_verdict` rule already covers the canonical-triplet path, so the legacy `_apply_canonical_override` is effectively dead / post-hoc duplicate logic.

**Impact:** low today (not a correctness bug; compute_verdict's Rule 1 covers canonical-triplet), but it means Session 1's "priority reconcile" commit (`352a59b`, referenced in v1.0 PR) left an unreferenced legacy helper in the monolith. Future reviewers will be confused about which is authoritative.

**Suggested action (NOT applying):** Either (a) delete `_apply_canonical_override` from the monolith if truly unused (grep shows only one self-reference on line 3391), or (b) update the design spec to state "canonical override logic is fully absorbed into `compute_verdict`; legacy helper retained in monolith pending Session 4 classify extraction."

### 1.3 ✅ Module responsibilities respected
Each extracted module owns exactly what the spec §2.3 table allocates:
- `primitives.py`: log/I-O/dataclasses — no behavior
- `filter_b_flank.py`: BLAST + slop-overlap — one pure sub-helper
- `filter_c_chimeric.py`: BLAST + offtarget aggregate — one pure sub-helper

Naming is consistent (`filter_b_flank`, `filter_c_chimeric` — matches spec §2.1).

### 1.4 ⚠️ `filter_c_chimeric.filter_c_offtarget_chrs` is dead weight (low)
The "pure filter-C check" at `scripts/s05/filter_c_chimeric.py:95-108` accepts `is_chimeric: bool` and immediately `del is_chimeric`s it. The docstring says it's kept "for call-site parity" with a deprecation note (Issue #11 M-3). This is acceptable technical debt but: **no caller is yet wired to use this function** (grep confirms it's only exported, not invoked). Adding API surface before any caller exists is premature abstraction by the "don't design for hypothetical future requirements" principle.

**Suggested action (NOT applying):** Either delete `filter_c_offtarget_chrs` until Session 5+ wires it, or add a TODO comment pointing to the exact Session where the wire-in happens.

---

## 2. Code Quality

### 2.1 Style drift between sibling modules (low)

| File | Optional type hint style |
|---|---|
| `filter_b_flank.py` | `Optional[tuple[str, int, int]]` (imports `from typing import Optional`) |
| `filter_c_chimeric.py` | (no optional types used) |
| `primitives.py` | `JunctionCluster | None` (PEP 604) |

Pick one. Repo is Python 3.11 → prefer PEP 604 (`X | None`) everywhere. This will bite readability as more modules land.

### 2.2 Silent BLAST failures (medium)

**`filter_b_flank.py:70-74`:**
```python
result = subprocess.run(
    ["blastn", ...],
    stderr=subprocess.DEVNULL,
)
if result.returncode != 0 or not blast_out.exists():
    return []
```

Same pattern in `filter_c_chimeric.py:60-65`. Both:
1. Discard stderr completely.
2. On failure (nonzero returncode OR missing output), return empty with no log message.

In the monolith these checks were inline and surrounded by context; extracting them into standalone modules amplifies the invisibility. If BLAST silently fails during production step-5 runs, the Filter B / Filter C evidence will be empty and verdict will be falsely CANDIDATE — **a silent false-positive risk**.

**Suggested action:** replace `stderr=subprocess.DEVNULL` with `stderr=subprocess.PIPE`, and on failure emit `log(f"[filter_b] blast failed rc={result.returncode}: {result.stderr[:200]}")` via the shared `primitives.log`. Zero behavior change on success; debuggability ↑.

### 2.3 `test_s05_import_dag.py` misnamed test function (low)

```python
def test_stage_ranking_is_sorted_and_unique_or_tied() -> None:
    """Sanity check: every ranked module file exists under scripts/s05/."""
```

The function name claims it tests sorting/uniqueness, but the body only checks "required Session-1 modules exist". **Docstring is right, name is wrong.** Rename to `test_required_session1_modules_present` or similar.

### 2.4 `scripts/reports/insertion_pdf.py:254-259` silent content truncation (low)

```python
for line in body_lines:
    ax.text(0.05, y, line, fontsize=line_fontsize, family=family, ...)
    y -= step
    if y < 0.04:
        break  # ← silently drops remaining lines
```

A page with 45+ lines (e.g., a sample with many elements in the foreign-elements block) would silently truncate. The paginator `_paginate` at line 264 handles known-bounded lists but not single-page text pages with variable body length.

**Suggested action:** either paginate text pages that overflow, or log a warning when truncation happens. At minimum, a comment saying "max ~35 lines per page" next to `step = 0.022`.

### 2.5 `extract_site_anchors` edge case (low)

`scripts/reports/insertion_pdf.py:139-157` parses sites like `"Chr3:16,439,674-16,439,709"`. The regex `r"^([^\s:]+):\s*([\d,]+)"` will also match `"Chr3:16,439,674 (35bp deletion)"` but **not** a negative coordinate or a site with an alternate format like `"chr:NODE_5_length_1234:500"`. Today no such sites exist, but as s05 evolves this will silently drop anchors.

**Suggested action:** add one unit test covering a malformed site string and confirm it's skipped (not crash-inducing). Current `extract_site_anchors` already handles this via `try/except ValueError`, so defensive coverage is fine; the test is the gap.

### 2.6 Hardcoded verdict set in junction diagram (low)

`scripts/reports/insertion_pdf.py:391`:
```python
interesting = {"CANDIDATE", "TRUE_INSERTION"}
```

`docs/measurements/verdict_priority.md` (Issue #3) defined the verdict lattice as: `canonical_triplet > host_endogenous > B > C > D > A > elements_present → CAND > UNK_FP > UNK`. Multiple "interesting" verdicts exist (e.g., `UNK_FP` with elements is arguably informative). Hardcoding the filter inside the PDF renderer creates a drift risk with verdict.py.

**Suggested action:** import the "interesting" set from `scripts/s05/verdict.py` or a shared config, or accept as a CLI argument. Mid-term: add to `VerdictRules` as `interesting_for_reports: set[str]`.

---

## 3. Maintainability

### 3.1 ✅ Every new module docstring links back to the spec
`primitives.py`, `filter_b_flank.py`, `filter_c_chimeric.py`, `test_s05_import_dag.py` all reference
`docs/superpowers/specs/2026-04-19-s05-module-split-design.md`. A future agent landing Session 2 will have a direct pointer to the plan. Excellent.

### 3.2 ✅ Commit messages are scannable
`Issue #4 Session 1 [1/4]`, `[2/4]`, etc. — one-glance visibility into which session / slot a commit occupies. Bodies reference "pytest 215→225 PASS" explicitly. Good grep-ability.

### 3.3 ⚠️ `run_pipeline.py:649` — help text omits step 8 (medium, UX)

```python
"--steps", type=str, default="1-5",
help=(
    "Steps to run (1, 2, 3, 4, 4b, 5, 6, 7). "  # ← missing 8!
    ...
),
```

`STEP_ORDER = [..., "7", "8"]` and `STEP_NAMES["8"] = "PDF insertion report ..."` but the `--help` output tells users step 8 doesn't exist. Users running `python run_pipeline.py --help` will not discover the feature they just paid for.

**Also:** example epilog (line 631-635) doesn't mention `--steps 8` or `--steps 1-8`.

**Suggested fix (one-line edit, NOT applying):**
```python
"Steps to run (1, 2, 3, 4, 4b, 5, 6, 7, 8). "
```
Plus one new example: `python run_pipeline.py --steps 8           # PDF report only`.

### 3.4 Test file duplication alert (low)

`tests/test_insertion_pdf.py` has two import-style patterns:
- Lines 9-26: direct `importlib.util.spec_from_file_location` into a module object.
- No `conftest.py` fixture for the same.

Same loader helper is used 14 times inside this file. Mild DRY violation but keeps each test self-contained, which is defensible. Monitor — don't refactor yet.

### 3.5 Inline comments remain synced
The monolith retains *pointer* comments (`filter_b_flank.py` at `scripts/s05_insert_assembly.py:134-136`) where constants used to live:
```python
# Construct-flanking filter constants (CONSTRUCT_FLANK_PIDENT, _MIN_LEN, _SLOP)
# moved to scripts/s05/filter_b_flank.py (Issue #4 Session 1) ...
```
This is the right call for an in-progress refactor: it orients readers during the transition. Remove in Session 6 or 7 when the monolith becomes a thin entrypoint.

---

## 4. User Experience

### 4.1 ⚠️ Step 8 discoverability (see 3.3 above)

### 4.2 PDF copy-number chart axis handling (low)

`scripts/reports/insertion_pdf.py:475-481` plots "Construct median", "Host median", "Depth ratio" as three bars on the **same y-axis**. Depth values are typically 20-50 reads, ratio is 0.5-3.0 — three orders of magnitude apart. The depth-ratio bar will be invisible when plotted against median-depth bars.

**Suggested improvement:** two-axis plot (primary y for depth, secondary for ratio), or two separate subplots. Current output is technically correct but will be misread.

### 4.3 Missing CRISPR panel (deferred)
`_page_placeholder(pdf, "CRISPR Editing", "editing panel pending v1.1")` is explicit. Acceptable placeholder, but note: **Issue #6 text says the PDF should include junction diagrams per site** — Agent #2 interpreted this as a genome-scale scatter (line 382-430). This is a *useful* view but **not the same as a per-site soft-clip pileup**, which regulators typically want for insertion characterization. Worth a conversation with the user about which view is the regulatory requirement.

### 4.4 CoC appendix ordering (low)
`_page_appendix_coc` renders entries in **file order** (as they appear in `chain_of_custody.jsonl`). Entries are typically chronological by `ts`, but nothing enforces it — malformed logs (e.g., out-of-order `pre`/`post`) would render misleadingly. Sort by `ts` before render, or add a CoC integrity check (which exists: `tools/verify_coc.py` from Issue #8 `5ed4f66`).

### 4.5 ✅ Graceful YAML fallback
`_resolve_host_construct_from_config` catches `ImportError` and malformed YAML, returning `(None, None)` without crashing. Good defensive code for optional dependencies.

---

## 5. Minor notes

- **`scripts/s05/__init__.py`** re-exports underscored names (`_read_fq_seqs`, `_find_construct_flanking_regions`, `_site_overlaps_flanking`, `_check_chimeric_assembly`). This is unusual — normally underscored = private. Here it's deliberate for back-compat with monolith callers. Note with a one-line comment: `# Underscored names are re-exported for monolith backward compatibility; Session 7 cleanup will drop leading underscores.`

- **`.gitignore`** correctly anchors new patterns with `/` (root-only). Good.

- **`atyucca6_gap_blast/` and `rice_G281_foreign/`** are now explicitly gitignored. These are audit evidence — if a regulator asks "how did you validate rice_G281 ground truth", these dirs are the answer. Gitignoring them means they won't travel with the repo. **Might be worth an `element_db/`-style `db/`-style allowlist instead**, or committing a summary .md pointing to where they're archived.

- **`tests/test_s05_import_dag.py:58-75` STAGE dict** is declared at module scope (fine for pytest parametrize) — but each entry has a comment about its future home. Consider freezing this after Session 7 into a simple `list[str]` ordering.

---

## Takeaways (priority-ordered, for follow-up)

1. **(medium)** Fix `run_pipeline.py --help` to mention step 8 (§3.3).
2. **(medium)** Un-silence BLAST failures in filter_b / filter_c (§2.2).
3. **(medium)** Decide fate of `_apply_canonical_override` in monolith — delete or document why it stays (§1.2).
4. **(low)** Unify `Optional` vs `| None` style across s05 modules (§2.1).
5. **(low)** Rename `test_stage_ranking_is_sorted_and_unique_or_tied` (§2.3).
6. **(low)** Add two-axis layout for copy-number bars in PDF (§4.2).
7. **(low)** Centralize "interesting verdict set" (§2.6).
8. **(low)** Add test for malformed site strings in `extract_site_anchors` (§2.5).

None of these block merge; all are polish. Branch is ready to ship as-is for an internal milestone; address items 1–3 before an external stakeholder review.
