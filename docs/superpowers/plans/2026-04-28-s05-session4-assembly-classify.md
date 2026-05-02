# s05 Session 4 — Assembly + Classify Extraction (Issue #4)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Extract Phase 3 (Assembly) into `scripts/s05/assembly.py` and fully extract Phase 1.5 (Classify) into `scripts/s05/classify.py` (replacing the existing shim). Two commits per spec §4 Session 4.

**Architecture:** Per `tests/test_s05_import_dag.py` STAGE: assembly = stage 3, classify = stage 1. Both modules become standalone (verbatim function bodies); the monolith re-imports their public symbols for backward compatibility. `assemble_insert` is the core of the pipeline — verbatim move + snapshot regression tests are the primary safety net. The 3-host 15x e2e gate from the design spec is **deferred** (hours of HPC compute); snapshot fixtures + full pytest are the substitute regression signal.

**Tech Stack:** Python 3.11, pytest 9, pysam. No new deps.

---

## Pre-state (after Session 3 cleanup `bf0c929`)

- monolith: 3618 lines
- pytest: 238 PASS + 1 skipped
- DAG: 13 PASSED
- BLAST smoke: green
- Snapshots: 26 PASSED, zero drift

---

## Commit 1: Extract `scripts/s05/assembly.py`

**Files:**
- Create: `scripts/s05/assembly.py` (~1100 lines)
- Modify: `scripts/s05/per_site.py` (shim — source `assemble_insert` and `refine_with_foreign_reads` from `.assembly`; keep `extract_candidate_reads`/`extract_unmapped_paired` from monolith pending Session 5)
- Modify: `scripts/s05/__init__.py` (re-export new symbols)
- Modify: `scripts/s05_insert_assembly.py` (add imports to BOTH ImportError branches, delete original definitions, add tombstone)

### Phase 3 functions to extract (verbatim)

All functions live in the monolith between `# Phase 3: K-mer Extension + Pilon Gap Fill` (line 1683) and `# Phase 4 (Annotation)` tombstone (~line 2798). Order in the new module mirrors the monolith order:

1. `class StrandAwareSeedExtender` — line 1686, ~136 lines
2. `def recruit_by_kmer` — line 1826
3. `def pilon_fill` — line 1875
4. `def check_host_termination` — line 1986
5. `def extract_foreign_reads` — line 2050
6. `def refine_with_foreign_reads` — line 2127
7. `def _minimap2_extend` — line 2268
8. `def _vote_extension` — line 2362
9. `def _ssake_extend` — line 2388
10. `def _check_merge` — line 2474
11. `def _write_pool_fastq` — line 2494
12. `def assemble_insert` — line 2511 (public entry point, ~280 lines)

Also extract any section header comments between phases (`# ---` lines) verbatim into appropriate spots in the new module.

### Steps

- [ ] **Step 1: Pre-state check**
  ```bash
  cd /data/gpfs/assoc/pgl/develop/redgene
  export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"
  git log --oneline -3                    # HEAD = bf0c929
  pytest tests/ -q                         # 238 PASS + 1 skipped
  pytest tests/test_s05_import_dag.py -v   # 13 PASSED
  pytest tests/test_verdict_snapshots.py -v   # 26 PASSED
  git show HEAD:scripts/s05_insert_assembly.py > /tmp/s05_pre_session4.py
  ```

- [ ] **Step 2: Create `scripts/s05/assembly.py`**

  Module docstring:
  ```python
  """Phase 3 — Assembly (Issue #4, Session 4).

  Extracted verbatim from ``scripts/s05_insert_assembly.py``. Owns the
  k-mer extension + Pilon gap fill + minimap2/SSAKE refinement pipeline
  that turns soft-clip seeds into a contiguous insert assembly.

  Pipeline overview
  -----------------
  ``assemble_insert`` is the public entry point invoked once per
  ``InsertionSite`` after Phase 1.5 classify. It orchestrates:

  1. ``StrandAwareSeedExtender`` — k-mer index over candidate reads
  2. ``recruit_by_kmer`` — pull reads matching seed k-mers
  3. ``_minimap2_extend`` + ``_vote_extension`` — extend contigs by tail voting
  4. ``_ssake_extend`` — fallback assembler for stalled extensions
  5. ``pilon_fill`` — gap-fill via Pilon polishing
  6. ``_check_merge`` — try to merge 5'/3' contigs into a single insert
  7. ``check_host_termination`` — confirm host flanking reached
  8. ``extract_foreign_reads`` + ``refine_with_foreign_reads`` — recruit
     the read pool relevant to the assembled insert and re-polish

  This module is the highest-risk extraction in the s05 split: it owns
  the core assembly logic that produces the FASTA consumed by Phase 4
  annotation and the four post-assembly FP filters. Verbatim move +
  ``tests/fixtures/verdict_snapshots/`` are the primary regression
  signal. The 3-host 15x e2e gate from the design spec is deferred
  (hours of HPC compute); snapshot fixtures cover the verdict-level
  surface that downstream consumers care about.
  """
  from __future__ import annotations
  ```

  Then imports: stdlib (`subprocess`, `tempfile`, `shutil`, `pathlib.Path`, `collections.Counter`, etc. — copy what the monolith needs at top), `pysam`, and from the s05 package:

  ```python
  from .primitives import (
      log, revcomp, read_fasta, write_fasta, _read_fq_seqs,
      InsertionSite, JunctionCluster,
  )
  ```

  (Other dataclasses or constants used inside Phase 3 — add them here. The implementer should grep each function body for unqualified names that resolve to monolith-level state, and import them all from `.primitives`.)

  Then verbatim copies of all 12 items in monolith order. Use
  ```bash
  awk 'NR>=1683 && NR<=2797' /tmp/s05_pre_session4.py > /tmp/phase3_block.py
  ```
  to capture the source range and `diff` against your new module's body for byte-identity verification.

  **Module-level constants** that Phase 3 functions depend on (e.g., `CLIP_HOST_PIDENT`, `CLIP_HOST_MIN_LEN` if used in assembly — verify via `grep`) must move with the functions. Mirror them at the top of `assembly.py`. If a constant is also used by classify or elsewhere, leave a copy in the monolith too (it stays a cheap import surface).

- [ ] **Step 3: Update `scripts/s05/per_site.py` shim**

  Current shim re-exports `extract_candidate_reads`, `extract_unmapped_paired`, `assemble_insert`, `refine_with_foreign_reads`. After Commit 1, the latter two move to `.assembly`. The former two stay in monolith pending Session 5.

  ```python
  """Phase 2-3 per-site re-exports (shim, v1.0).

  Session 4 (2026-04-28): ``assemble_insert`` and ``refine_with_foreign_reads``
  now live in ``scripts/s05/assembly.py``; this shim sources them from there
  while ``extract_candidate_reads`` / ``extract_unmapped_paired`` continue
  to ride on the monolith pending Session 5 (read_extraction.py).
  """
  from .assembly import (  # noqa: F401
      assemble_insert,
      refine_with_foreign_reads,
  )
  from scripts.s05_insert_assembly import (  # noqa: F401
      extract_candidate_reads,
      extract_unmapped_paired,
  )
  ```

- [ ] **Step 4: Update `scripts/s05/__init__.py`**

  After existing re-export blocks (filters, annotation), add:

  ```python
  from .assembly import (  # noqa: F401
      StrandAwareSeedExtender,
      assemble_insert,
      check_host_termination,
      extract_foreign_reads,
      pilon_fill,
      recruit_by_kmer,
      refine_with_foreign_reads,
  )
  ```

  Update the docstring with a Session 4 line.

- [ ] **Step 5: Wire monolith — add imports to BOTH ImportError branches**

  After the `from s05.annotation import (...)` block, add a new group in BOTH branches:

  ```python
      from s05.assembly import (
          StrandAwareSeedExtender,
          _check_merge,
          _minimap2_extend,
          _ssake_extend,
          _vote_extension,
          _write_pool_fastq,
          assemble_insert,
          check_host_termination,
          extract_foreign_reads,
          pilon_fill,
          recruit_by_kmer,
          refine_with_foreign_reads,
      )
  ```

  Include the underscore-private helpers because the monolith may still call them as bare names (from `assemble_insert`'s body? No — those calls move with the function. But `extract_foreign_reads` calls `_minimap2_extend`? Verify by grep.) **Run `grep -n "_minimap2_extend\|_vote_extension\|_ssake_extend\|_check_merge\|_write_pool_fastq" scripts/s05_insert_assembly.py` AFTER deleting the function bodies.** If any monolith call sites remain, the import must include those names.

- [ ] **Step 6: Delete original Phase 3 block + add tombstone**

  Delete monolith lines 1683–2797 (Phase 3 section header through end of `assemble_insert`). Replace with:

  ```python
  # ---------------------------------------------------------------------------
  # Phase 3 (Assembly) functions moved to scripts/s05/assembly.py
  # (Issue #4 Session 4) and re-imported at the top of this module for
  # backward compatibility with existing callers (Phase 2 → Phase 3 →
  # Phase 4 → classify_site_tiers chain remains bit-identical).
  # ---------------------------------------------------------------------------
  ```

- [ ] **Step 7: DAG test**
  ```bash
  pytest tests/test_s05_import_dag.py -v
  ```
  Expected: **14 PASSED** (new `[assembly.py]` parametrize entry).

- [ ] **Step 8: Full pytest + snapshots + CLI smoke**
  ```bash
  pytest tests/ -q                          # 239 passed, 1 skipped (was 238 + 1)
  pytest tests/test_verdict_snapshots.py -v # 26 PASSED, ZERO DRIFT
  pytest tests/test_extra_element_db.py -v  # all green
  python scripts/s05_insert_assembly.py --help > /dev/null && echo OK
  ```
  Snapshot drift = revert immediately. Verbatim move must produce zero drift.

- [ ] **Step 9: Commit 1**
  ```bash
  git add scripts/s05/assembly.py scripts/s05/per_site.py scripts/s05/__init__.py scripts/s05_insert_assembly.py
  git commit -m "$(cat <<'EOF'
  Issue #4 Session 4 [1/2]: extract assembly.py from s05 monolith

  Phase 3 assembly (StrandAwareSeedExtender + 11 helpers + assemble_insert
  public entry point, ~1100 lines) moves verbatim into scripts/s05/assembly.py.
  Highest-risk extraction in the s05 split — assemble_insert is the core
  pipeline function that produces the contig consumed by Phase 4 annotation
  and the four post-assembly FP filters.

  per_site.py shim now sources assemble_insert / refine_with_foreign_reads
  from .assembly; Phase 2 helpers (extract_candidate_reads,
  extract_unmapped_paired) continue to ride on the monolith until
  Session 5.

  DAG: assembly pre-registered at stage 3. 14 PASSED including
  [assembly.py] parametrize entry.

  Tests: 238 PASS + 1 skipped → 239 PASS + 1 skipped (DAG +1). 26
  verdict snapshots green with zero drift; BLAST smoke green; CLI OK.

  3-host 15x e2e gate from spec §4 deferred (hours of HPC compute);
  snapshot fixtures + full pytest are the regression signal.

  Per docs/superpowers/specs/2026-04-19-s05-module-split-design.md §4 Session 4.

  Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
  EOF
  )"
  ```

---

## Commit 2: Expand `scripts/s05/classify.py`

**Files:**
- Modify: `scripts/s05/classify.py` (replace shim with native definitions; keep `_parse_src_tag` re-export from `.primitives`)
- Modify: `scripts/s05/__init__.py` (no symbol changes; comment update)
- Modify: `scripts/s05_insert_assembly.py` (add imports to BOTH ImportError branches, delete original definitions, add tombstone)

### Functions and constants to extract (verbatim)

All from monolith Phase 1.5 region:

1. `_batch_check_element_hits` — line 302 (in Phase 1 area but used by classify; per spec §2.3 this lives in classify alongside the tier logic)
2. `_filter_host_endogenous` — line 801, ~133 lines
3. `_SRC_TIER` constant — line 908
4. `_TIER2_SRCS` constant — line 923
5. `_UNKNOWN_SRC_WARNED` module-level set — line 926
6. `_should_replace` — line 934
7. `classify_site_tiers` — line 966, ~310 lines
8. `write_tier_classification` — line 1276

### Steps

- [ ] **Step 1: Pre-state check (after Commit 1)**
  ```bash
  git log --oneline -1   # HEAD = Commit 1 SHA
  pytest tests/ -q        # 239 PASS + 1 skipped
  ```

- [ ] **Step 2: Replace `scripts/s05/classify.py` with native module**

  Keep the existing `_parse_src_tag` re-export from `.primitives` at the top. Then add module-level constants and all 5 function bodies verbatim. Final structure:

  ```python
  """Phase 1.5 — Classify (Issue #4, Session 4).

  Extracted verbatim from ``scripts/s05_insert_assembly.py``. Owns the
  tier-based classification logic that turns Phase 1 soft-clip junctions
  + Phase 4 annotation hits into TierResult / InsertionSite verdicts.

  Stage 1 in the s05 DAG (depends only on ``primitives``). Imports
  ``_parse_src_tag`` from ``.primitives`` (Session 3 migrated it there
  to break a runtime cycle with ``annotation``).

  Replaces the v1.0 shim that re-exported from the monolith.
  """
  from __future__ import annotations

  # imports per the monolith body — typing, dataclasses, primitives.

  from .primitives import (
      log, _parse_src_tag,
      InsertionSite, TierResult,  # whatever the bodies actually use
  )
  # plus any other s05 imports the bodies need (limited to stage 0)
  ```

  Then 8 verbatim items in this order: `_batch_check_element_hits`, `_filter_host_endogenous`, `_SRC_TIER`, `_TIER2_SRCS`, `_UNKNOWN_SRC_WARNED`, `_should_replace`, `classify_site_tiers`, `write_tier_classification`.

  **Critical:** preserve every comment, including the BUG-3 / cd-hit / cluster references in `_TIER2_SRCS`'s docstring. They are load-bearing audit trail.

- [ ] **Step 3: Update monolith ImportError branches**

  Change the existing `from s05.classify import ...` (which currently imports nothing because the shim only had re-exports) to:

  ```python
      from s05.classify import (
          _batch_check_element_hits,
          _filter_host_endogenous,
          _should_replace,
          _SRC_TIER,
          _TIER2_SRCS,
          _UNKNOWN_SRC_WARNED,
          classify_site_tiers,
          write_tier_classification,
      )
  ```

  Both branches.

- [ ] **Step 4: Delete original definitions from monolith + tombstone**

  Delete:
  - `_batch_check_element_hits` body (line 302)
  - `_filter_host_endogenous` body (line 801)
  - `_SRC_TIER` / `_TIER2_SRCS` / `_UNKNOWN_SRC_WARNED` (lines 908–926)
  - `_should_replace` body (line 934)
  - `classify_site_tiers` body (line 966)
  - `write_tier_classification` body (line 1276)

  Replace each block with a small tombstone (or one umbrella tombstone at the start of the deleted region):

  ```python
  # ---------------------------------------------------------------------------
  # Phase 1.5 (Classify) functions moved to scripts/s05/classify.py
  # (Issue #4 Session 4) and re-imported at the top of this module for
  # backward compatibility. _SRC_TIER / _TIER2_SRCS / _UNKNOWN_SRC_WARNED
  # module state migrated with them.
  # ---------------------------------------------------------------------------
  ```

  Note that `_batch_check_element_hits` lives in the Phase 1 region (line 302), so it gets a separate tombstone there.

- [ ] **Step 5: DAG test**
  ```bash
  pytest tests/test_s05_import_dag.py -v
  ```
  Expected: **14 PASSED** (no parametrize change — `classify.py` was already at stage 1; just transitions from shim to native).

  Critically: `[classify.py]` parametrize must still PASS. classify.py's only allowed s05 sibling import is `from .primitives import ...`. If the implementer accidentally imports from another sibling (e.g., `assembly` or `annotation`), DAG fails.

- [ ] **Step 6: Full pytest + snapshots**
  ```bash
  pytest tests/ -q                          # still 239 passed, 1 skipped
  pytest tests/test_verdict_snapshots.py -v # 26 PASSED, ZERO DRIFT
  pytest tests/test_extra_element_db.py -v  # green
  python scripts/s05_insert_assembly.py --help > /dev/null && echo OK
  ```

- [ ] **Step 7: Commit 2**
  ```bash
  git commit -am "$(cat <<'EOF'
  Issue #4 Session 4 [2/2]: expand classify.py from s05 monolith

  Phase 1.5 classify (8 items: _batch_check_element_hits,
  _filter_host_endogenous, _SRC_TIER + _TIER2_SRCS +
  _UNKNOWN_SRC_WARNED module state, _should_replace, classify_site_tiers,
  write_tier_classification) moves verbatim from monolith into
  scripts/s05/classify.py, replacing the v1.0 re-export shim.

  Closes Session 4. Sessions 5/6/7 remaining: read_extraction +
  site_discovery (Session 5), report + fanout + thin entrypoint
  (Session 6), docs cleanup + shim retirement (Session 7).

  Tests: 239 PASS + 1 skipped (no count change; classify.py transitions
  from shim to native, parametrize entry was already there). 26 verdict
  snapshots green; BLAST smoke green; CLI OK.

  Per docs/superpowers/specs/2026-04-19-s05-module-split-design.md §4 Session 4.

  Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
  EOF
  )"
  ```

---

## Risks and rollback

| Risk | Mitigation |
|------|------------|
| `assemble_insert` byte-drift causing snapshot regression | Verbatim policy; if snapshots drift, revert and diff line-by-line |
| Hidden cross-call inside Phase 3 (e.g., `extract_foreign_reads` calls `_check_merge`) breaking after extraction | Step 5/Step 3 grep audit ensures all referenced names are imported |
| Module-level constants left orphaned in monolith | Step 4 scans the deleted region; constants used by both old + new code get mirrored, not deleted |
| `pysam` import cost | `assembly.py` imports pysam at module top — same as the monolith. No regression. |
| `classify_site_tiers` runtime state (`_UNKNOWN_SRC_WARNED`) — is it shared? | The set is a module-level singleton. After extraction, it lives only in `classify.py`. Monolith's call sites resolve to the imported reference, so the set is shared correctly. Verify by checking that the warning behavior in `classify_site_tiers` still de-duplicates across calls. |

Rollback: each commit is independent. `git revert` or `git reset --hard HEAD~N`.
