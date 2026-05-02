# s05 Session 6 — Report + Fanout + Thin Entrypoint (Issue #4)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`).

**Goal:** Extract `generate_report` + `write_stats` into `scripts/s05/report.py`; extract `main()` orchestration into `scripts/s05/fanout_orchestrator.py`; reduce monolith to thin entrypoint shim. Three commits per spec §4 Session 6.

**Architecture:** report = stage 8 (already in DAG STAGE), fanout_orchestrator = stage 9 (top of DAG). Verbatim move + monolith thin re-export. After Session 6: monolith is < 200 lines (target < 120 per spec); fanout_orchestrator owns the CLI.

## Pre-state at HEAD `14889ad`
- monolith: 1127 lines
- pytest: 240 PASS + 1 skipped
- DAG: 15 PASSED
- Snapshots: 26 PASSED, zero drift

## Target post-state (after all 3 commits)
- monolith: ~150 lines (close to but possibly above the 120 spec target — relax if line-budget guard test isn't trivially achievable)
- new `report.py`: ~325 lines (generate_report ~300 + write_stats ~25)
- new `fanout_orchestrator.py`: ~470 lines (main + helpers)
- pytest: 242 PASS + 1 skipped (DAG +2 for report + fanout_orchestrator)
- DAG: 17 PASSED

---

## Commit 1: Extract `scripts/s05/report.py`

**Files:**
- Create: `scripts/s05/report.py`
- Modify: `scripts/s05/annotate_report.py` shim (source `generate_report`/`write_stats` from `.report` instead of monolith)
- Modify: `scripts/s05/__init__.py` (re-export new symbols)
- Modify: `scripts/s05_insert_assembly.py` (BOTH ImportError branches + delete originals + tombstone)

### Items to extract (verbatim, monolith pre-Session-6)

- `def generate_report` — lines 336–635 (~300 lines, the verdict/CANDIDATE/FALSE_POSITIVE classification + report writing)
- `def write_stats` — lines 640–659 (~20 lines, sample-level stats TSV writer)

### Steps

- [ ] **Step 1: Pre-state**
  ```bash
  cd /data/gpfs/assoc/pgl/develop/redgene
  export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"
  git log --oneline -3                       # HEAD = 14889ad
  pytest tests/ -q                            # 240 + 1 skipped
  pytest tests/test_s05_import_dag.py -v      # 15 PASSED
  pytest tests/test_verdict_snapshots.py -v   # 26 PASSED
  git show HEAD:scripts/s05_insert_assembly.py > /tmp/s05_pre_session6.py
  ```

- [ ] **Step 2: Create `scripts/s05/report.py`**

  Module docstring:
  ```python
  """Phase 5 — Report + Stats (Issue #4, Session 6).

  Extracted verbatim from ``scripts/s05_insert_assembly.py``. Owns the
  per-site report generation (insertion_report.tsv with verdict labels)
  and the sample-level stats TSV. Stage 8 in the s05 DAG.

  Functions
  ---------
  - ``generate_report`` — writes per-site insertion_report.tsv with
    embedded CANDIDATE/FALSE_POSITIVE/UNKNOWN verdict logic. NOTE: the
    verdict logic is duplicated with ``scripts/s05/verdict.py::compute_verdict``;
    Session 7 cleanup may consolidate.
  - ``write_stats`` — sample-level summary TSV (insertion_sites count,
    candidate_read_pairs, per-site result rows).

  Replaces the v1.0 ``annotate_report.py`` shim's monolith re-export of
  these two symbols.
  """
  from __future__ import annotations
  ```

  Imports — copy whatever generate_report and write_stats use. From s05:
  ```python
  from .primitives import (log, InsertionSite)
  from .verdict import (FilterEvidence, compute_verdict, ...)  # if used
  ```
  (Grep the bodies for unqualified names; import everything they need from primitives/verdict.)

  Verbatim function bodies. Use:
  ```bash
  awk 'NR>=336 && NR<=635' /tmp/s05_pre_session6.py > /tmp/generate_report.py
  awk 'NR>=640 && NR<=659' /tmp/s05_pre_session6.py > /tmp/write_stats.py
  ```
  Diff against new module bodies for byte-identity.

- [ ] **Step 3: Update `annotate_report.py` shim**

  Replace monolith import with `.report` import:
  ```python
  """Phase 4-5 annotate + report re-exports (shim, v1.0).

  Session 6 (2026-04-29): generate_report and write_stats now live in
  ``scripts/s05/report.py``; this shim sources them from there. Session
  3 already migrated annotate_insert to ``scripts/s05/annotation.py``.
  The shim itself can be retired in Session 7.
  """
  from .annotation import annotate_insert  # noqa: F401
  from .report import generate_report, write_stats  # noqa: F401
  ```

- [ ] **Step 4: Update `__init__.py`**

  Re-export `generate_report`, `write_stats`. Add Session 6 line to docstring.

- [ ] **Step 5: Wire monolith ImportError branches**

  Add `from s05.report import (generate_report, write_stats)` to BOTH branches.

- [ ] **Step 6: Delete originals + tombstone**

  Replace lines 336–659 (generate_report, write_stats, surrounding section dividers) with:
  ```python
  # ---------------------------------------------------------------------------
  # Phase 5 (Report + Stats) functions moved to scripts/s05/report.py
  # (Issue #4 Session 6) and re-imported at the top of this module.
  # ---------------------------------------------------------------------------
  ```

- [ ] **Step 7: Tests**
  ```bash
  pytest tests/test_s05_import_dag.py -v   # 16 PASSED ([report.py] new)
  pytest tests/ -q                          # 241 PASS + 1 skipped (DAG +1)
  pytest tests/test_verdict_snapshots.py -v # 26 PASSED, ZERO DRIFT
  python scripts/s05_insert_assembly.py --help > /dev/null && echo OK
  ```

- [ ] **Step 8: Commit 1**
  ```bash
  git add scripts/s05/report.py scripts/s05/annotate_report.py scripts/s05/__init__.py scripts/s05_insert_assembly.py
  git commit -m "$(cat <<'EOF'
  Issue #4 Session 6 [1/3]: extract report.py from s05 monolith

  Phase 5 report + stats (generate_report, write_stats) move verbatim
  from monolith into scripts/s05/report.py. annotate_report.py shim
  now sources both from .report. Monolith retains thin re-imports.

  DAG: report pre-registered at stage 8. 16 PASSED including [report.py].
  Tests: 240 PASS + 1 skipped → 241 PASS + 1 skipped (DAG +1).
  26 verdict snapshots green; zero drift.

  Per docs/superpowers/specs/2026-04-19-s05-module-split-design.md §4 Session 6.

  Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
  EOF
  )"
  ```

---

## Commit 2: Extract `scripts/s05/fanout_orchestrator.py`

**Files:**
- Create: `scripts/s05/fanout_orchestrator.py`
- Modify: `scripts/s05/__init__.py` (re-export `main`)
- Modify: `scripts/s05_insert_assembly.py` (delete `main` + replace with thin call to fanout_orchestrator.main)

### What to extract

`main()` (lines 664–~1127, currently ~460 lines including the four-phase dispatch). The spec's line budget is `fanout_orchestrator::main < 250 lines` — to meet it, split monolithic `main()` into:

- `main()` — argparse + dispatch (target < 100 lines)
- `_run_phase_1_1_5(args, ...)` — Phase 1 + 1.5 logic
- `_run_phase_2_3(args, ...)` — Phase 2 + 3 logic
- `_run_phase_4(args, ...)` — Phase 4 logic

The split is mechanical (existing `if args.phase in (...)` blocks already define the boundaries). Keep variable names identical to source. The split is the ONE non-verbatim change in this session — document it in the commit message.

**If the split is too risky for one commit**, fall back to verbatim move (main() over 460 lines) and note that line budget is missed. The spec has fallback language: "If violated, split dispatch into separate run_phase_* functions" — do the split if reasonable.

### Steps

- [ ] **Step 1: Pre-state (after Commit 1)**

- [ ] **Step 2: Create `scripts/s05/fanout_orchestrator.py`**

  Module docstring + imports + helper functions + main(). Imports include EVERYTHING `main()` references — pull from primitives, site_discovery, classify, read_extraction, assembly, annotation, verdict, report, config_loader, plus stdlib + filter modules (filter_a/b/c/d for FilterEvidence usage).

  Verbatim function bodies for the four `_run_phase_*` helpers (each carved from the corresponding `if args.phase in (...)` block in original main()). Variable names + control flow IDENTICAL to source.

  ```python
  """Phase 6 — Fanout Orchestrator (Issue #4, Session 6).

  Extracted from ``scripts/s05_insert_assembly.py``. Owns the CLI dispatch
  for s05's four phases (1_1.5 site discovery, 2_3 per-site assembly,
  4 annotation, 'all' end-to-end). Stage 9 in the s05 DAG (top of pipeline).

  main() splits the original monolith ``main()`` into:
  - main: argparse + phase dispatch (~80 lines)
  - _run_phase_1_1_5: Phase 1 + 1.5 (site discovery + classify)
  - _run_phase_2_3: Phase 2 + 3 (read extraction + assembly per site)
  - _run_phase_4: Phase 4 (annotation, FP filter eval, report writing)

  The split is mechanical — variables and control flow preserved verbatim
  from the monolith. Snapshot fixtures + DAG test verify no behavioral drift.
  """
  ```

- [ ] **Step 3: Replace monolith `main()` body with thin call**

  ```python
  # ---------------------------------------------------------------------------
  # Main entrypoint (orchestration moved to scripts/s05/fanout_orchestrator.py)
  # ---------------------------------------------------------------------------

  def main() -> None:
      """Thin entrypoint — full orchestration in scripts/s05/fanout_orchestrator."""
      from s05.fanout_orchestrator import main as _orchestrator_main
      _orchestrator_main()


  if __name__ == "__main__":
      main()
  ```

- [ ] **Step 4: Update `__init__.py`**

  Re-export `main` (or `fanout_main`) from `.fanout_orchestrator`. Update docstring.

- [ ] **Step 5: Tests**
  ```bash
  pytest tests/test_s05_import_dag.py -v   # 17 PASSED ([fanout_orchestrator.py] new)
  pytest tests/ -q                          # 242 PASS + 1 skipped (DAG +1)
  pytest tests/test_verdict_snapshots.py -v # 26 PASSED
  python scripts/s05_insert_assembly.py --help > /dev/null && echo OK
  python -m scripts.s05.fanout_orchestrator --help > /dev/null && echo "fanout OK"
  ```

  **Fanout smoke** (per spec exit gate) — if a 5x rice sample is available, do:
  ```bash
  # See spec §3.4 Step 4. SKIP if no fixture available locally — defer to Session 7.
  ```

- [ ] **Step 6: Commit 2**
  ```bash
  git add scripts/s05/fanout_orchestrator.py scripts/s05/__init__.py scripts/s05_insert_assembly.py
  git commit -m "$(cat <<'EOF'
  Issue #4 Session 6 [2/3]: extract fanout_orchestrator.py from s05 monolith

  main() (~460 lines) extracted into scripts/s05/fanout_orchestrator.py
  and split into 4 functions:
    main: argparse + phase dispatch
    _run_phase_1_1_5: site discovery + classify
    _run_phase_2_3: per-site read extraction + assembly
    _run_phase_4: annotation + FP filter eval + report

  Split is mechanical — variable names and control flow preserved
  verbatim from monolith if-blocks. Monolith's main() reduced to a
  thin call into fanout_orchestrator.main.

  DAG: fanout_orchestrator pre-registered at stage 9. 17 PASSED.
  Tests: 241 → 242 PASS + 1 skipped. 26 verdict snapshots green.

  Per docs/superpowers/specs/2026-04-19-s05-module-split-design.md §4 Session 6.

  Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
  EOF
  )"
  ```

---

## Commit 3: Thin entrypoint + line-budget guard

**Files:**
- Modify: `scripts/s05_insert_assembly.py` (clean up to thin entrypoint < 200 lines)
- Optionally update: `tests/test_submit_s05_array.py` to add line-budget assertion (per spec exit gate)

### Steps

- [ ] **Step 1: Audit monolith**

  After Commits 1+2, monolith has:
  - shebang + docstring
  - `try: ... except ImportError:` block (still needed for backward-compat re-exports)
  - tombstone comments
  - Module-level constants that callers may import
  - thin `main()` (just calls fanout_orchestrator)

  Goal: bring total < 200 lines (relax from spec's 120 if needed). Possible simplifications:
  - Consolidate the try/except block — many of the imports may be redundant (e.g., if no caller imports `_minimap2_extend` from the monolith, drop it from the import group)
  - Drop module-level constants that no caller uses
  - Drop forward `from s05.X import` groups for symbols nothing imports from monolith

  **Audit method**: for each name imported in the try/except block, run:
  ```bash
  grep -RnE "from scripts.s05_insert_assembly import.*\b<NAME>\b" --include="*.py"
  ```
  If nothing matches OUTSIDE the monolith, the import is dead — drop it.

- [ ] **Step 2: Update `tests/test_submit_s05_array.py` (or create `tests/test_line_budget.py`)**

  Add an assertion:
  ```python
  def test_monolith_line_budget():
      monolith = REPO_ROOT / "scripts/s05_insert_assembly.py"
      lines = monolith.read_text().splitlines()
      assert len(lines) < 200, (
          f"monolith should be a thin entrypoint after Session 6; "
          f"got {len(lines)} lines"
      )

  def test_fanout_main_line_budget():
      import ast
      fanout = REPO_ROOT / "scripts/s05/fanout_orchestrator.py"
      tree = ast.parse(fanout.read_text())
      for node in ast.walk(tree):
          if isinstance(node, ast.FunctionDef) and node.name == "main":
              n_lines = node.end_lineno - node.lineno
              assert n_lines < 250, f"main() too large: {n_lines} lines"
              return
      raise AssertionError("main() not found in fanout_orchestrator.py")
  ```

  (Adjust threshold to whatever the actual achievable size is — relax if necessary.)

- [ ] **Step 3: Tests**
  ```bash
  pytest tests/ -q
  pytest tests/test_s05_import_dag.py -v
  python scripts/s05_insert_assembly.py --help > /dev/null && echo OK
  ```

- [ ] **Step 4: Commit 3**
  ```bash
  git add scripts/s05_insert_assembly.py tests/test_submit_s05_array.py
  git commit -m "$(cat <<'EOF'
  Issue #4 Session 6 [3/3]: thin entrypoint + line-budget guard

  Monolith reduced to NN lines (target < 200). Dead re-exports dropped
  after grep audit confirmed no external callers. Line-budget tests
  added to lock in the slim shape:
    - scripts/s05_insert_assembly.py: < 200 lines
    - scripts/s05/fanout_orchestrator.py::main: < 250 lines

  Closes Session 6. Session 7 remaining: docs + retire 4 legacy shims
  (annotate_report, classify, per_site, site_discovery already-renamed).

  Per docs/superpowers/specs/2026-04-19-s05-module-split-design.md §4 Session 6.

  Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
  EOF
  )"
  ```

---

## Risks

- `main()` split: this is the ONE non-verbatim change in the entire refactor. If snapshot tests drift after Commit 2, revert immediately and try a verbatim move (no split).
- The `try: from s05.X import (...)` pattern was specifically designed to make the monolith load standalone. After thinning, verify that dropping unused imports doesn't break standalone invocation: `python scripts/s05_insert_assembly.py --help`.
- Line budget < 120 is aggressive. If after audit the monolith stabilizes at 150–200 lines, document and proceed; the spec's number is a target, not a hard requirement.

Rollback per commit.
