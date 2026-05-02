# s05 Session 5 — Site Discovery + Read Extraction (Issue #4)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`).

**Goal:** Expand `scripts/s05/site_discovery.py` from shim to native module (Phase 1) and create `scripts/s05/read_extraction.py` from Phase 2 (currently in monolith). Two commits per spec §4 Session 5.

**Architecture:** site_discovery = stage 1 (DAG STAGE pre-registered), read_extraction = stage 2. Verbatim move + monolith re-imports + snapshot regression test as primary safety net. After Session 5: import migration is half-complete (per spec exit gate). Lazy import in `assembly.py::assemble_insert` for `extract_unmapped_paired` should be cleaned up to use `from .read_extraction import` once read_extraction lands.

## Pre-state at HEAD `ac4265d`
- monolith: 2026 lines
- pytest: 239 PASS + 1 skipped
- DAG: 14 PASSED
- Snapshots: 26 PASSED, zero drift

## Target post-state
- monolith: ~1100 lines (drop ~900 from Phase 1 + Phase 2)
- new `site_discovery.py`: ~470 lines (was 13-line shim)
- new `read_extraction.py`: ~450 lines
- pytest: 240 PASS + 1 skipped (DAG +1 for read_extraction)
- DAG: 15 PASSED

---

## Commit 1: Expand `scripts/s05/site_discovery.py`

### Phase 1 items to extract (verbatim, in monolith order)

From `scripts/s05_insert_assembly.py` lines ~252–910 (Phase 1 section through end of `legacy_junctions_to_sites`):

1. `# Phase 1: Soft-clip Junction Detection` section header (line 252)
2. `def _build_consensus` — line 255
3. `def _batch_check_maps_to_host` — line 307
4. `def find_softclip_junctions` — line 354 (large, ~320 lines)
5. `def _apply_mask_bed` — line 677
6. `MASKED_SOURCE_TAG` constant (search via `grep -n "MASKED_SOURCE_TAG ="` to find definition line)
7. `def parse_legacy_junctions` — line 815
8. `def legacy_junctions_to_sites` — line 846

`_extract_seeds_at_positions` (line 912) is a Phase 2 helper — DO NOT move it here; it goes with read_extraction in Commit 2.

### Steps

- [ ] **Step 1: Pre-state**
  ```bash
  cd /data/gpfs/assoc/pgl/develop/redgene
  export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"
  git log --oneline -3                      # HEAD = ac4265d
  pytest tests/ -q                           # 239 PASS + 1 skipped
  pytest tests/test_s05_import_dag.py -v     # 14 PASSED
  pytest tests/test_verdict_snapshots.py -v  # 26 PASSED
  git show HEAD:scripts/s05_insert_assembly.py > /tmp/s05_pre_session5.py
  ```

- [ ] **Step 2: Replace `scripts/s05/site_discovery.py` with native module**

  Module structure (mirrors Session 4 classify.py pattern):
  ```python
  """Phase 1 — Site Discovery (Issue #4, Session 5).

  Extracted verbatim from ``scripts/s05_insert_assembly.py``. Owns the
  soft-clip junction detection pipeline that turns a host-mapped BAM
  into a list of ``InsertionSite`` candidates. Stage 1 in the s05 DAG;
  imports only from ``.primitives`` (stage 0).

  Pipeline overview
  -----------------
  ``find_softclip_junctions`` is the public entry point. It scans the
  host BAM for soft-clipped reads, clusters their clips into per-position
  ``JunctionCluster`` objects, and builds ``InsertionSite`` records for
  downstream Phase 1.5 classify and Phase 2/3 read-extraction +
  assembly.

  Helpers
  -------
  - ``_build_consensus`` — majority-vote consensus over a clip stack
  - ``_batch_check_maps_to_host`` — BLAST clip seqs to host, mark host-mapping
  - ``_apply_mask_bed`` — pre-filter sites by user-supplied BED mask (T10)
  - ``parse_legacy_junctions`` / ``legacy_junctions_to_sites`` — fallback
    path that ingests step-6 junctions.tsv when running step 5 alone

  Replaces the v1.0 shim that re-exported these names from the monolith.
  """
  from __future__ import annotations
  ```

  Then imports — copy whatever the monolith body uses (`subprocess`, `pathlib.Path`, `pysam`, `collections.Counter`/`defaultdict`, etc.). From s05:
  ```python
  from .primitives import (
      log, revcomp, read_fasta, write_fasta,
      InsertionSite, JunctionCluster,
  )
  ```
  (Plus any other primitives symbols used inside Phase 1 — grep to confirm.)

  Then verbatim copies of all 8 items in monolith order. Each function body byte-identical to source. Use:
  ```bash
  awk 'NR>=252 && NR<=910' /tmp/s05_pre_session5.py > /tmp/phase1_block.py
  ```
  Diff the new module's bodies against this for verification.

- [ ] **Step 3: Update monolith ImportError branches**

  Currently the monolith imports `from s05.site_discovery import ...` is implicit (only the package's namespace re-exports). After Commit 1, the monolith needs explicit imports. Either:
  - Update existing `try: ... except ImportError: ...` block in the monolith to add `from s05.site_discovery import (...)` with all 8 names, OR
  - Reuse existing pattern: just add to each branch alongside the assembly/classify/annotation/filter import groups.

  The 8 names: `_build_consensus`, `_batch_check_maps_to_host`, `find_softclip_junctions`, `_apply_mask_bed`, `MASKED_SOURCE_TAG`, `parse_legacy_junctions`, `legacy_junctions_to_sites`. The 8th (already covered) is `MASKED_SOURCE_TAG` — count rechecks.

  Add to BOTH branches.

- [ ] **Step 4: Delete original Phase 1 block + tombstone**

  Delete lines ~252–910 from monolith. Replace with:
  ```python
  # ---------------------------------------------------------------------------
  # Phase 1 (Site Discovery) functions moved to scripts/s05/site_discovery.py
  # (Issue #4 Session 5) and re-imported at the top of this module for
  # backward compatibility with existing callers.
  # ---------------------------------------------------------------------------
  ```

- [ ] **Step 5: DAG + tests**
  ```bash
  pytest tests/test_s05_import_dag.py -v   # 14 PASSED (no parametrize change; site_discovery transitions shim → native)
  pytest tests/ -q                          # 239 PASS + 1 skipped (no count change)
  pytest tests/test_verdict_snapshots.py -v # 26 PASSED, ZERO DRIFT
  python scripts/s05_insert_assembly.py --help > /dev/null && echo OK
  ```

- [ ] **Step 6: Commit 1**
  ```bash
  git add scripts/s05/site_discovery.py scripts/s05/__init__.py scripts/s05_insert_assembly.py
  git commit -m "$(cat <<'EOF'
  Issue #4 Session 5 [1/2]: expand site_discovery.py from s05 monolith

  Phase 1 site discovery (8 items: _build_consensus, _batch_check_maps_to_host,
  find_softclip_junctions, _apply_mask_bed, MASKED_SOURCE_TAG,
  parse_legacy_junctions, legacy_junctions_to_sites) moves verbatim from
  monolith into scripts/s05/site_discovery.py, replacing the v1.0 re-export
  shim.

  DAG: site_discovery already at stage 1; transitions from shim to native.
  Tests: 239 PASS + 1 skipped (no count change). 26 verdict snapshots
  green; zero drift. CLI OK.

  Per docs/superpowers/specs/2026-04-19-s05-module-split-design.md §4 Session 5.

  Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
  EOF
  )"
  ```

---

## Commit 2: Create `scripts/s05/read_extraction.py`

### Phase 2 items to extract (verbatim, in monolith order)

From `scripts/s05_insert_assembly.py` (post-Commit-1) Phase 2 region:

1. `def _extract_seeds_at_positions` — was line 912, now ~250-ish after Commit 1's deletion
2. `# Phase 2: Candidate Read Extraction` section header
3. `def extract_candidate_reads` — was line 974
4. `def extract_unmapped_paired` — was line 1167

Plus any other helpers in this region (grep `^def ` between Phase 1.5 tombstone and Phase 3 tombstone after Commit 1 lands).

### Steps

- [ ] **Step 1: Pre-state (after Commit 1)**

- [ ] **Step 2: Create `scripts/s05/read_extraction.py`**

  Module structure:
  ```python
  """Phase 2 — Read Extraction (Issue #4, Session 5).

  Extracted verbatim from ``scripts/s05_insert_assembly.py``. Owns the
  candidate read extraction pipeline that pulls clip-anchored reads + their
  mates from the host BAM for downstream Phase 3 assembly. Stage 2 in the
  s05 DAG; imports from ``.primitives`` (stage 0) and ``.site_discovery``
  (stage 1) — same-stage import (`.classify`) also allowed but not used.

  Functions
  ---------
  - ``_extract_seeds_at_positions`` — Phase 2 helper, fills ``InsertionSite``
    seed_5p / seed_3p attributes from BAM at known junction positions
  - ``extract_candidate_reads`` — public entry; extracts reads + mates near
    insertion sites, writes paired FASTQs for assembly
  - ``extract_unmapped_paired`` — pulls unmapped read pairs (recovery for
    sites where one mate failed to map)

  Replaces the per-site shim re-exports from Sessions 1–4.
  """
  from __future__ import annotations
  ```

  Imports: stdlib + pysam + from primitives + possibly site_discovery (if Phase 2 calls `_build_consensus` directly).

  Verbatim function bodies. Use:
  ```bash
  awk 'NR>=BEFORE_PHASE2_TOMBSTONE && NR<=AFTER_PHASE2_TOMBSTONE' /tmp/s05_pre_session5.py > /tmp/phase2_block.py
  ```
  (The boundaries need recomputing after Commit 1 — use grep for `# Phase 2:` and `# Phase 3` tombstones in the post-Commit-1 monolith.)

- [ ] **Step 3: Update `scripts/s05/per_site.py` shim**

  Source `extract_candidate_reads` / `extract_unmapped_paired` from `.read_extraction`. Keep `assemble_insert` / `refine_with_foreign_reads` from `.assembly` (unchanged from Session 4):
  ```python
  """Phase 2-3 per-site re-exports (shim, v1.0).

  Session 5 (2026-04-28): all 4 re-exports now come from native modules
  (.assembly, .read_extraction). The shim itself can be retired in Session 7.
  """
  from .assembly import (  # noqa: F401
      assemble_insert,
      refine_with_foreign_reads,
  )
  from .read_extraction import (  # noqa: F401
      extract_candidate_reads,
      extract_unmapped_paired,
  )
  ```

- [ ] **Step 4: Update `scripts/s05/assembly.py` lazy import**

  In `assemble_insert`, the lazy import added in Session 4 should now resolve to `.read_extraction` instead of the monolith. Replace:
  ```python
  # OLD (Session 4):
  from scripts.s05_insert_assembly import extract_unmapped_paired
  ```
  with:
  ```python
  # NEW (Session 5):
  from .read_extraction import extract_unmapped_paired
  ```

  This eliminates the cycle entirely (no lazy needed); could even be hoisted to module top. **However**, hoisting requires read_extraction to be at stage ≤ 2; it is. So hoisting IS valid:
  ```python
  # At module top of assembly.py:
  from .read_extraction import extract_unmapped_paired
  ```
  Either approach is acceptable. The implementer's choice. If hoisting, remove the lazy block from `assemble_insert` body.

  **Risk:** the lazy form was introduced specifically to avoid a cycle. Confirm cycle is broken now: assembly (stage 3) imports from read_extraction (stage 2) — DAG strict `>` rule allows this (stage 2 < stage 3). Cycle-free.

- [ ] **Step 5: Update `scripts/s05/__init__.py`**

  Re-export read_extraction symbols. Update docstring with Session 5 line.

- [ ] **Step 6: Update monolith ImportError branches**

  Add `from s05.read_extraction import (...)` with all 3 symbols (`_extract_seeds_at_positions`, `extract_candidate_reads`, `extract_unmapped_paired`) to BOTH branches.

- [ ] **Step 7: Delete original Phase 2 block + tombstone**

  Replace with:
  ```python
  # ---------------------------------------------------------------------------
  # Phase 2 (Read Extraction) functions moved to scripts/s05/read_extraction.py
  # (Issue #4 Session 5) and re-imported at the top of this module.
  # ---------------------------------------------------------------------------
  ```

- [ ] **Step 8: DAG + tests**
  ```bash
  pytest tests/test_s05_import_dag.py -v
  pytest tests/ -q
  pytest tests/test_verdict_snapshots.py -v
  pytest tests/test_extra_element_db.py -v
  python scripts/s05_insert_assembly.py --help > /dev/null && echo OK
  ```
  Expected: DAG 15 PASSED (new `[read_extraction.py]` parametrize), pytest 240 PASS + 1 skipped, snapshots green.

  **DAG STAGE update may be needed:** if `tests/test_s05_import_dag.py::STAGE` doesn't already have `read_extraction: 2`, add it. (The forward-looking entries in the dict from Session 1 should already include this — verify.)

- [ ] **Step 9: Commit 2**
  ```bash
  git add scripts/s05/read_extraction.py scripts/s05/__init__.py \
          scripts/s05/per_site.py scripts/s05/assembly.py \
          scripts/s05_insert_assembly.py
  git commit -m "$(cat <<'EOF'
  Issue #4 Session 5 [2/2]: extract read_extraction.py from s05 monolith

  Phase 2 read extraction (3 functions: _extract_seeds_at_positions,
  extract_candidate_reads, extract_unmapped_paired) moves verbatim from
  monolith into scripts/s05/read_extraction.py.

  per_site.py shim now sources all 4 symbols from native modules
  (.assembly + .read_extraction). assembly.py's Session-4 lazy import for
  extract_unmapped_paired is replaced by a regular `.read_extraction`
  import (cycle gone).

  Closes Session 5: import migration half-complete per spec exit gate.
  Sessions 6 + 7 remaining: report.py + fanout_orchestrator + thin
  entrypoint (Session 6); shim retirement + docs (Session 7).

  Tests: 239 PASS + 1 skipped → 240 PASS + 1 skipped (DAG +1 for
  read_extraction.py). 26 verdict snapshots green; zero drift.

  Per docs/superpowers/specs/2026-04-19-s05-module-split-design.md §4 Session 5.

  Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
  EOF
  )"
  ```

---

## Risks

- Phase 1 / Phase 2 functions reference monolith-level constants — grep audit before deletion.
- `find_softclip_junctions` is large (~320 lines) — extra care on byte-identity diff.
- The Session 4 lazy import in assembly.py needs cleanup; if missed, assembly continues to work but carries dead code.
- `_extract_seeds_at_positions` is between Phase 1.5 tombstone and Phase 2 in current monolith — boundaries shift after Commit 1 deletion; recompute line ranges before Commit 2.

Rollback per commit: `git revert <sha>`.
