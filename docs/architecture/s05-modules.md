# s05 Module Architecture

*This document supersedes `docs/superpowers/specs/2026-04-19-s05-module-split-design.md`
as the authoritative architecture reference for the `scripts/s05/` package.
Issue #4 (7-session refactor) is complete as of 2026-04-29.*

## Overview

`scripts/s05/` is a 11-module pure-Python package that implements the core insert
assembly pipeline (Step 5). It was extracted from a 4003-line monolith
(`scripts/s05_insert_assembly.py`) across Sessions 1-7 of Issue #4.

The monolith remains as a thin entrypoint (102 lines) that:
- Delegates to `fanout_orchestrator.main()` for execution
- Re-exports a small set of symbols needed by importlib-based test fixtures
- Is executable directly (`python scripts/s05_insert_assembly.py --help`)

## Modules

| Module | Stage | Role |
|--------|-------|------|
| `primitives.py` | 0 | `log`, `revcomp`, FASTA/FASTQ I/O, dataclasses (`JunctionCluster`, `InsertionSite`, `LegacyJunction`, `TierResult`) |
| `site_discovery.py` | 1 | Phase 1: soft-clip junction detection, `find_softclip_junctions`, `_apply_mask_bed`, `MASKED_SOURCE_TAG` |
| `classify.py` | 1 | Tier classification (`classify_site_tiers`, `_SRC_TIER`, `_TIER2_SRCS`, `_batch_check_element_hits`) |
| `read_extraction.py` | 2 | Phase 2: `extract_candidate_reads`, `extract_unmapped_paired` |
| `assembly.py` | 3 | Phase 3: k-mer extension, Pilon gap fill, minimap2/SSAKE refinement (`assemble_insert`, `refine_with_foreign_reads`) |
| `annotation.py` | 4 | Phase 4: BLAST-driven element annotation (`annotate_insert`, `_run_local_blast`, `_run_remote_blast`) |
| `filter_a_host.py` | 5 | Post-assembly host-fraction false positive filter |
| `filter_b_flank.py` | 5 | Construct-flanking false positive filter (BLAST construct vs host, slop-overlap) |
| `filter_c_chimeric.py` | 5 | Multi-locus chimeric false positive filter (strict-identity off-target chrom) |
| `filter_d_altlocus.py` | 5 | Alternative-locus false positive filter (minimap2, host-derived elements) |
| `verdict.py` | 6 | Pure `compute_verdict` + `FilterEvidence` + `VerdictRules` + `_apply_canonical_override` + `REPORT_INTERESTING_VERDICTS` |
| `config_loader.py` | 7 | `load_verdict_rules` (YAML -> `VerdictRules`) |
| `report.py` | 8 | Phase 5: per-site report generation (`generate_report`, `write_stats`) |
| `fanout_orchestrator.py` | 9 | Phase 6: main() CLI dispatch + phase helpers; sole entry point |

## DAG Stage Map

Strict layering enforced by `tests/test_s05_import_dag.py`. A module at stage N
may only import from stages <= N; reverse imports are forbidden.

```
[0] primitives
       |
  +----+----+
  v         v
[1] site_discovery   [1] classify
       |
       v
[2] read_extraction
       |
       v
[3] assembly
       |
       v
[4] annotation
       |
  +----+----+----+----+
  v    v    v    v
[5] filter_a  filter_b  filter_c  filter_d
  +----+----+----+----+
       |
       v
[6] verdict  <-- [7] config_loader
       |
       v
[8] report
       |
       v
[9] fanout_orchestrator
```

## Public Surface (`__init__.py` re-exports)

`scripts/s05/__init__.py` re-exports all public symbols from every module so
callers can do `from scripts.s05 import <symbol>` without knowing the source
module. Key groups:

- **Primitives**: `STEP`, `log`, `revcomp`, `read_fasta`, `write_fasta`, `JunctionCluster`, `InsertionSite`, `LegacyJunction`, `TierResult`
- **Filters**: All `CONSTRUCT_*`, `CHIMERIC_*`, `INSERT_HOST_*`, `CONSTRUCT_HOST_*` constants + their check functions
- **Annotation**: `annotate_insert`, `_parse_blast6`, `_run_local_blast`, `_run_remote_blast`, `_merge_annotations`
- **Verdict**: `compute_verdict`, `FilterEvidence`, `VerdictRules`, `_apply_canonical_override`, `REPORT_INTERESTING_VERDICTS`
- **Site discovery**: `find_softclip_junctions`, `_apply_mask_bed`, `MASKED_SOURCE_TAG`, `parse_legacy_junctions`, `legacy_junctions_to_sites`
- **Assembly**: `assemble_insert`, `refine_with_foreign_reads`, `pilon_fill`, `recruit_by_kmer`, `StrandAwareSeedExtender`
- **Orchestration**: `main` (from `fanout_orchestrator`)

The monolith entrypoint (`scripts/s05_insert_assembly.py`) re-exports a smaller
subset: only symbols actively referenced by importlib-based test fixtures
(`annotate_insert`, `_batch_check_element_hits`, `_should_replace`, `_SRC_TIER`,
`_TIER2_SRCS`).

## Verdict Priority Order

```
canonical_triplet > host_endogenous > Filter B > Filter C > Filter D > Filter A
  > elements_present -> CAND > UNK_FP > UNK
```

Documented in `docs/measurements/verdict_priority.md`; implemented in
`scripts/s05/verdict.py:compute_verdict`. Changes require snapshot regeneration.

## Refactor History (Issue #4)

| Session | Date | Commits | What was extracted |
|---------|------|---------|-------------------|
| Session 1 | 2026-04-19 | `39da513`, `a9230d7`, `81f22fe`, `65b1c6e` | `primitives`, `filter_b_flank`, `filter_c_chimeric`; DAG test added |
| Session 2 | 2026-04-28 | `2c59f1e`, `6658fb8` | `filter_a_host`, `filter_d_altlocus` |
| Session 3 | 2026-04-28 | `6225d42`, `bf0c929` | `annotation` (Phase 4 BLAST-driven element annotation) |
| Session 4 | 2026-04-28 | `aeee876`, `ac4265d` | `assembly` (Phase 3 k-mer/Pilon/minimap2), `classify` promoted |
| Session 5 | 2026-04-28 | `a11e6ea`, `14889ad` | `site_discovery` (Phase 1), `read_extraction` (Phase 2) |
| Session 6 | 2026-04-29 | `972e0d7`, `90a4343`, `96d16fd` | `report` (Phase 5), `fanout_orchestrator` (Phase 6); monolith thinned to 262 lines |
| Session 7 | 2026-04-29 | *(this commit)* | Retired `annotate_report.py` + `per_site.py` shims; thinned monolith to 102 lines; migrated test imports |

**Net result**: 4003-line monolith -> 102-line thin entrypoint + 11 native modules.
Pytest: 235 PASS (pre-refactor) -> 242 PASS (post-Session 7; +7 net new tests, -2
DAG entries for retired shims vs Session 6). Zero snapshot drift across all sessions.

## Thin Entrypoint

`scripts/s05_insert_assembly.py` is the subprocess entry point wired in
`run_pipeline.py` (step 5) and SLURM batch arrays. It must remain at this path.
Its only logic is:

1. A try/except re-export block (5 symbols for importlib-based test fixtures)
2. A `main()` shim that delegates to `fanout_orchestrator.main()`
3. `if __name__ == "__main__": main()`

Do not add logic here. All new functionality goes into a `scripts/s05/` module.
