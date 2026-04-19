# s05 Module Split — Design Spec (Issue #4)

**Date:** 2026-04-19
**Scope:** v1.1 multi-day refactor of `scripts/s05_insert_assembly.py` (4234 lines) into a 10-module package.
**Issue:** [#4 [v1.1] s05_insert_assembly.py 전체 8-10 모듈 분할](https://github.com/wyim-pgl/redgene/issues/4)
**Related commits:** `7a3224a` (T7 shim package, v1.0), `a5ee725` (T8 fanout, v1.0), `352a59b` (T6 compute_verdict split, v1.0).

## 1. Goals and non-goals

### Goals

- Turn the 4234-line monolith into **10 focused modules** with clear single-responsibility boundaries.
- Keep the existing CLI (`python scripts/s05_insert_assembly.py --host-bam ...`) fully backward-compatible so `run_pipeline.py`, SLURM array scripts, and external operators need no changes.
- Preserve every current behavior: snapshot tests (`test_verdict_snapshots.py` rice/cucumber golden), filter priority order, CoC `run_step()` context, `--phase` fan-out semantics.
- Enable clean per-module unit tests and make future refactors (Issue #4 follow-ups) cheaper.
- Leave a documented architecture (`docs/architecture/s05-modules.md`) for v2 contributors.

### Non-goals

- **No behavior change.** Every snapshot test passes unchanged. No verdict, no filter threshold, no assembly parameter is modified.
- **No `run_pipeline.py` refactor.** The orchestrator keeps invoking `scripts/s05_insert_assembly.py` exactly as before.
- **No new features.** PDF report (#6), KCGP (#7), and other v1.1 items stay out of this refactor's scope.
- **No deletion of the monolith file.** It survives as a thin entrypoint. Full removal (if desired) is a v2 decision.
- **No test-behavior changes.** Tests that currently use `importlib.util.spec_from_file_location` against `scripts/s05_insert_assembly.py` keep working; we do not force migration of that pattern.

## 2. Architecture (final state)

### 2.1 Package layout

```
scripts/
├── s05_insert_assembly.py           ~100 lines — thin entrypoint (argparse + dispatch only)
└── s05/
    ├── __init__.py                  ~30 lines — public API re-exports for b/c
    ├── primitives.py          NEW   ~150 lines — dataclasses + FASTA/FASTQ I/O + log + revcomp
    ├── site_discovery.py      EXPAND ~700 lines — soft-clip junctions, batch host/element checks, mask BED
    ├── classify.py            EXPAND ~500 lines — tier classification, host_endogenous filter
    ├── read_extraction.py     RENAMED ~400 lines — legacy junctions, candidate reads, unmapped mates (was per_site.py)
    ├── assembly.py            NEW   ~1200 lines — SeedExtender, k-mer recruit, Pilon gap fill, assemble_insert
    ├── annotation.py          NEW   ~450 lines — local/remote BLAST, merge_annotations, annotate_insert
    ├── filter_a_host.py       NEW   ~100 lines — Filter A (host_fraction + small gap)
    ├── filter_b_flank.py      NEW   ~150 lines — Filter B (construct-flanking)
    ├── filter_c_chimeric.py   NEW   ~150 lines — Filter C (chimeric multi-chrom)
    ├── filter_d_altlocus.py   NEW   ~200 lines — Filter D (construct+host coverage, alt-locus)
    ├── verdict.py             EXPAND ~360 lines — compute_verdict + _apply_canonical_override (moved from monolith)
    ├── report.py              NEW   ~330 lines — generate_report, write_stats
    ├── fanout_orchestrator.py NEW   ~200 lines — --phase dispatch (1_1.5 / 2_3 / 4 / all)
    └── config_loader.py       UNCHANGED 85 lines — load_verdict_rules, DEFAULT_TRIPLETS
```

Total: **15 files in `scripts/s05/` + 1 thin entrypoint = 16 files**, of which:

- **9 modules match Issue #4 spec verbatim**: `assembly.py`, `annotation.py`, `filter_a_host.py`, `filter_b_flank.py`, `filter_c_chimeric.py`, `filter_d_altlocus.py`, `verdict.py` (already extracted in T6), `report.py`, `fanout_orchestrator.py`.
- **5 supporting modules added for DAG clarity** (not in Issue #4 text but agreed during 2026-04-19 brainstorming): `primitives.py`, `site_discovery.py` (promoted from shim), `classify.py` (promoted from shim), `read_extraction.py` (renamed from `per_site.py` shim), `config_loader.py` (already extracted in T6).
- **1 package marker**: `__init__.py`.

User-facing count = **14 modules**; Issue #4's "8-10 modules" text is interpreted as a coarse target and the richer decomposition is justified by the single-responsibility-per-file principle (see §2.3).

### 2.2 Dependency DAG (no cycles)

```
                      primitives.py
                     /             \
           site_discovery.py   classify.py
                    |               |
                    └───────┬───────┘
                            v
                    read_extraction.py
                            |
                            v
                       assembly.py
                            |
                            v
                      annotation.py
                            |
         ┌──────┬───────────┼──────────┬──────┐
         v      v           v          v      v
     filter_a filter_b  filter_c  filter_d
         └──────┴───────────┴──────────┘
                         |
                         v
                      verdict.py ─── config_loader.py
                         |
                         v
                       report.py
                         |
                         v
                fanout_orchestrator.py
                         |
                         v
              s05_insert_assembly.py (thin entrypoint)
```

Cycle prevention is enforced by a new AST-walking unit test (`tests/test_s05_import_dag.py`) that fails if any module imports from a later-stage module.

### 2.3 Module responsibilities

| Module | Symbols |
|--------|---------|
| `primitives.py` | `log`, `revcomp`, `read_fasta`, `write_fasta`, `_read_fq_seqs`, `JunctionCluster`, `InsertionSite`, `LegacyJunction`, `TierResult` |
| `site_discovery.py` | `find_softclip_junctions`, `_build_consensus`, `_batch_check_maps_to_host`, `_batch_check_element_hits`, `_apply_mask_bed`, `MASKED_SOURCE_TAG` |
| `classify.py` | `_filter_host_endogenous`, `_parse_src_tag`, `_should_replace`, `classify_site_tiers`, `write_tier_classification` |
| `read_extraction.py` | `parse_legacy_junctions`, `legacy_junctions_to_sites`, `_extract_seeds_at_positions`, `extract_candidate_reads`, `extract_unmapped_paired` |
| `assembly.py` | `StrandAwareSeedExtender`, `recruit_by_kmer`, `pilon_fill`, `check_host_termination`, `extract_foreign_reads`, `refine_with_foreign_reads`, `_minimap2_extend`, `_vote_extension`, `_ssake_extend`, `_check_merge`, `_write_pool_fastq`, `assemble_insert` |
| `annotation.py` | `_parse_blast6`, `_run_local_blast`, `_run_remote_blast`, `_merge_annotations`, `annotate_insert` |
| `filter_a_host.py` | Pure filter-A check function; consumed by verdict.py |
| `filter_b_flank.py` | `_find_construct_flanking_regions`, `_site_overlaps_flanking`, pure filter-B check |
| `filter_c_chimeric.py` | `_check_chimeric_assembly`, pure filter-C check |
| `filter_d_altlocus.py` | `_check_construct_host_coverage`, `_blast_insert_vs_host`, pure filter-D check |
| `verdict.py` (expand) | `compute_verdict`, `FilterEvidence`, `VerdictRules`, `_find_matching_triplet`, `_apply_canonical_override` (moved from monolith line 127) |
| `report.py` | `generate_report`, `write_stats` |
| `fanout_orchestrator.py` | `run_phase_1_1_5`, `run_phase_2_3`, `run_phase_4`, `run_all`, top-level dispatch |

## 3. Migration strategy

### 3.1 Import surface (3 tiers)

**Tier 1 — stable, public (no change):**
- `from scripts.s05.verdict import compute_verdict, FilterEvidence, VerdictRules`
- `from scripts.s05.config_loader import load_verdict_rules`

**Tier 2 — test imports that need migration:**
Tests using `from scripts.s05_insert_assembly import X` are rewritten to `from scripts.s05.<new_module> import X`. Affected files (confirmed by grep):
- `tests/test_canonical_override.py`
- `tests/test_source_tag_priority.py`
- `tests/test_extra_element_db.py`
- `tests/test_element_db_build.py` (source-text scan, only path string updates)
- `tests/test_submit_s05_array.py` (M-2 guard: update target to new entry)
- `tests/test_mask_bed_intersect.py` (already uses shim, adjust to new module path)
- Any other file flagged by `grep -l "scripts.s05_insert_assembly" tests/`

Migration is mechanical (grep/sed driven) and covered by Session 7's cleanup commit.

**Tier 3 — importlib.util fixture pattern (keep working, no migration needed):**
Tests that do:
```python
_spec = importlib.util.spec_from_file_location("s05", REPO / "scripts/s05_insert_assembly.py")
s05 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(s05)
```
continue to point at the thin entrypoint. The entrypoint re-exports names via `from scripts.s05 import *` so that `s05.classify_site_tiers` etc. resolve transparently.

### 3.2 Public package API (`scripts/s05/__init__.py`)

```python
# Re-export everything currently callable via `from scripts.s05_insert_assembly import X`
from scripts.s05.primitives import (  # noqa: F401
    log, revcomp, read_fasta, write_fasta,
    JunctionCluster, InsertionSite, LegacyJunction, TierResult,
)
from scripts.s05.site_discovery import (  # noqa: F401
    find_softclip_junctions, _batch_check_maps_to_host, _batch_check_element_hits,
    _apply_mask_bed, MASKED_SOURCE_TAG, _build_consensus,
)
from scripts.s05.classify import (  # noqa: F401
    classify_site_tiers, write_tier_classification,
    _filter_host_endogenous, _parse_src_tag, _should_replace,
)
from scripts.s05.read_extraction import (  # noqa: F401
    parse_legacy_junctions, legacy_junctions_to_sites,
    extract_candidate_reads, extract_unmapped_paired,
)
from scripts.s05.assembly import (  # noqa: F401
    StrandAwareSeedExtender, assemble_insert,
    recruit_by_kmer, pilon_fill, refine_with_foreign_reads,
)
from scripts.s05.annotation import (  # noqa: F401
    annotate_insert, _run_local_blast, _run_remote_blast,
)
from scripts.s05.verdict import (  # noqa: F401
    compute_verdict, FilterEvidence, VerdictRules,
    _apply_canonical_override, _find_matching_triplet,
)
from scripts.s05.report import generate_report, write_stats  # noqa: F401
```

### 3.3 Thin entrypoint (`scripts/s05_insert_assembly.py` after split)

```python
#!/usr/bin/env python3
"""Step 5 thin entrypoint.

Delegates to scripts.s05.fanout_orchestrator. Preserves the CLI surface
run_pipeline.py and SLURM array scripts depend on. All logic lives in
scripts/s05/ submodules; see docs/architecture/s05-modules.md.
"""
from __future__ import annotations

# Public re-exports for backwards-compat (tests using importlib.util
# spec_from_file_location against this file still work).
from scripts.s05 import *  # noqa: F401,F403

from scripts.s05.fanout_orchestrator import main

if __name__ == "__main__":
    main()
```

### 3.4 Regression-guard rules

1. **Snapshot immutability** — `pytest tests/test_verdict_snapshots.py` green after every commit.
2. **Full pytest run** at the end of every session: 215 PASS + 1 skipped (baseline; may grow as new module-unit tests land).
3. **CLI smoke test** after every commit:
   ```bash
   python scripts/s05_insert_assembly.py --help > /dev/null && echo OK
   ```
4. **Fanout smoke test** after Session 6: run a short (5x) rice sample through `--phase 1_1.5` then `--phase 2_3 --site-id X` then `--phase 4`; diff insertion reports against a pre-refactor baseline.
5. **Line budget (M-2 guard)** — update `tests/test_submit_s05_array.py` to assert:
   - `scripts/s05_insert_assembly.py` total < 120 lines
   - `scripts/s05/fanout_orchestrator.py::main` < 250 lines
6. **No-cycle test** — new `tests/test_s05_import_dag.py` walks AST imports across all `scripts/s05/*.py` and asserts the DAG in §2.2 (fails CI if violated).
7. **3-host e2e** (Session 4 only, before merging assembly.py) — rerun rice_G281, tomato_Cas9_A2_3, cucumber_line225 at 15x and confirm verdict/GT anchor unchanged from `docs/measurements/coverage_sensitivity_matrix.tsv`.

## 4. Session-by-session plan

| Session | Focus | Modules | Commits (est.) | Exit gate |
|---------|-------|---------|----------------|-----------|
| **1** (today, 착수) | Design + low-risk first extractions | `primitives.py`, `filter_b_flank.py`, `filter_c_chimeric.py` | 4 (design commit + 3 module commits) | 215 PASS; CLI smoke; no snapshot change |
| 2 | Remaining filters | `filter_a_host.py`, `filter_d_altlocus.py` | 2 | 215+ PASS |
| 3 | Annotation | `annotation.py` | 1 | 215+ PASS; BLAST path smoke |
| 4 | Assembly (highest risk) | `assembly.py` + `classify.py` expand | 2 | 215+ PASS; **3-host 15x e2e** |
| 5 | Site discovery + read extraction | `site_discovery.py` expand, `read_extraction.py` rename | 2 | import migration half-complete |
| 6 | Report + fanout + thin entrypoint | `report.py`, `fanout_orchestrator.py`, thin entrypoint | 3 | LINE_BUDGET guard; fanout smoke |
| 7 | Documentation + cleanup | `docs/architecture/s05-modules.md`, delete 4 old shims, final test migration | 1-2 | pytest clean, docs complete, Issue #4 closable |

**Session 1 deliverables (this session only):**
1. This design spec (committed) → `docs/superpowers/specs/2026-04-19-s05-module-split-design.md`
2. Writing-plans output (to be created after spec approval) → `docs/superpowers/plans/2026-04-19-s05-module-split-plan.md`
3. `scripts/s05/primitives.py` — extracted, `__init__.py` re-exports updated, pytest green
4. `scripts/s05/filter_b_flank.py` — extracted, pytest green
5. `scripts/s05/filter_c_chimeric.py` — extracted, pytest green
6. New test: `tests/test_s05_import_dag.py` — DAG no-cycle guard

**Explicitly out of session 1:** `assembly.py`, `classify.py` expand, `annotation.py`, `report.py`, `fanout_orchestrator.py`, thin entrypoint rewrite, `--phase` dispatch migration, shim deletion.

## 5. Risks and mitigations

| # | Risk | Likelihood | Impact | Mitigation |
|---|------|-----------|--------|-----------|
| R1 | Circular import between extracted modules | Medium | High | DAG test enforced from Session 1; imports go through `__init__.py` re-exports rather than sibling paths |
| R2 | Snapshot test failure after assembly.py extraction | Low | High | Session 4 gate requires 3-host 15x e2e before commit is kept; full revert plan per commit |
| R3 | BUG-7 microhomology MAPQ regression (classify/site_discovery changes) | Low | Medium | Session 5 adds explicit replay of `rice_G281_foreign/` test fixtures (MAPQ 0-19 path) |
| R4 | M-2 line budget violation if fanout_orchestrator grows | Medium | Low | New guard at < 250 lines; if violated, split dispatch into separate `run_phase_*` functions |
| R5 | `run_pipeline.py` subprocess call breaks (CLI drift) | Low | High | Thin entrypoint preserves argparse verbatim; CLI smoke after every commit |
| R6 | Test migration churn breaks unrelated tests | Medium | Medium | Migration confined to Session 7; grep-audited file list; each rewrite is mechanical |
| R7 | Shim deletion breaks an external caller I didn't find | Low | Medium | Shim deletion only in Session 7, after package `__init__.py` re-exports verified; deprecation warning in shim bodies one session earlier |
| R8 | Session-gap drift (graph of changes grows between sessions) | Medium | Low | Each session starts by confirming git state clean and pytest green; design spec linked from every session commit |

## 6. Success criteria (Issue #4 closable)

- [ ] All 10 modules extracted with line counts within ±20% of estimates in §2.1
- [ ] `scripts/s05_insert_assembly.py` ≤ 120 lines
- [ ] `pytest tests/` passes ≥ current baseline (215 + 1 skipped), no snapshots changed
- [ ] DAG no-cycle test passes
- [ ] `docs/architecture/s05-modules.md` written and committed
- [ ] 3-host 15x e2e matches `coverage_sensitivity_matrix.tsv` baseline exactly
- [ ] CLI help output diff-clean vs pre-refactor snapshot
- [ ] Old shim files replaced: `classify.py` and `site_discovery.py` retain the same name but become real modules (no longer re-export shims); `per_site.py` file is removed, its symbols live in new `read_extraction.py`; `annotate_report.py` file is removed, its symbols split across new `annotation.py` and `report.py`.
- [ ] `gh issue close 4` with session summary comments on the issue

## 7. Open questions (for writing-plans)

- Exact ordering of module unit tests added per session (TDD or post-hoc?) — writing-plans to decide per TDD skill.
- Whether to add a deprecation shim for the old 4 re-export files between Sessions 6 and 7 (one-session overlap with `DeprecationWarning`), vs hard-delete in Session 7.
- Whether `primitives.py` should live in `scripts/util/` instead (cross-pipeline reuse). **Decision deferred to Session 1 implementation; default: keep inside `scripts/s05/` to avoid scope creep.**
