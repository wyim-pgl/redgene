# Code Review: Pipeline Reorganization (Step 9 → Step 5 centric)

**Date**: 2026-04-10
**Scope**: run_pipeline.py rewrite, 4 new scripts (s04-s07), filter_univec.py plant mode, s02 UniVec integration, CLAUDE.md update

---

## Critical Issues

### 1. CRASH: s11_multiqc.py still reads from removed steps (s04_assembly, s06_junction)

**File**: `scripts/s11_multiqc.py:40, 74`
```python
stats_file = outdir / sample / "s04_assembly" / "assembly_stats.txt"   # removed step!
junc_file = outdir / sample / "s06_junction" / "junctions.tsv"         # removed step!
```
`s04_assembly` and `s06_junction` no longer exist in the pipeline. These paths will silently fail (files won't exist), but the multiqc report will be missing sections without warning. Should reference `s05_insert_assembly/s05_stats.txt` instead.

### 2. CRASH: plot_insert_structure.py + plot_sample_summary.py reference `s09_stats.txt` (renamed to `s05_stats.txt`)

**File**: `scripts/viz/plot_insert_structure.py`, `scripts/viz/plot_sample_summary.py`
```
--stats results/{sample}/s05_insert_assembly/s09_stats.txt   # file is now s05_stats.txt!
stats = parse_stats(s09_dir / "s09_stats.txt")               # will FileNotFoundError
```
The stats file was renamed to `s05_stats.txt` inside `s05_insert_assembly.py`, but the viz scripts still look for `s09_stats.txt`. **This will crash at runtime.**

### 3. BUG: s02 `_build_combined_ref()` reads construct_ref 3 times

**File**: `scripts/s02_construct_map.py:64-78`
```python
out.write(construct_ref.read_text())              # read 1
if not construct_ref.read_text().endswith("\n"):   # read 2
    out.write("\n")
n_construct = sum(1 for line in construct_ref.read_text().splitlines()  # read 3
```
`read_text()` is called 3 times on the same file. For a typical construct ref (~10KB) this is harmless, but it's wasteful and breaks if the file changes between reads. Read once into a variable.

---

## Stale References (docstrings/help text)

### 4. WARNING: New scripts have old step numbers in docstrings

| File | Line | Stale Text |
|------|------|------------|
| `scripts/s04_host_map.py:2` | `"""Step 7: Map all trimmed reads...` | Should be Step 4 |
| `scripts/s05_insert_assembly.py:2` | `"""Step 9 — Targeted insert assembly...` | Should be Step 5 |
| `scripts/s05_insert_assembly.py:3269` | `description="Step 9: Targeted insert assembly"` | argparse |
| `scripts/s05_insert_assembly.py:3311` | `log(f"=== Step 9: Targeted Insert Assembly...` | runtime log |
| `scripts/s05_insert_assembly.py:3496` | `log(f"=== Step 9 complete ===")` | runtime log |
| `scripts/s06_indel.py:2` | `"""Step 8: CRISPR editing...` | Should be Step 6 |
| `scripts/s06_indel.py:20-30` | Docstring examples reference `s07_host_map`, `s08_indel_detection.py` | |
| `scripts/s06_indel.py:536` | `description="Step 8:..."` | argparse |
| `scripts/s07_copynumber.py:2` | `"""Step 10: Estimate transgene...` | Should be Step 7 |
| `scripts/s07_copynumber.py:15` | `{outdir}/{sample}/s10_copynumber/copynumber.tsv` | docstring path |

Only the log prefix (`[s04_host_map]`) and directory names were updated via `sed`; docstrings, argparse descriptions, and example commands were not. Users running `--help` will see old step numbers. Log messages at runtime will say "Step 9" when this is Step 5.

### 5. WARNING: CLAUDE.md still references "step 8" for gRNA config

**File**: `CLAUDE.md:96`
```
optional `grna` (path to gRNA file for step 8)
```
Should say "step 6".

### 6. WARNING: `--s09-dir` flag name in plot_sample_summary.py

**File**: `scripts/viz/plot_sample_summary.py`
The CLI flag is still `--s09-dir` while the directory is now `s05_insert_assembly`. The flag name is confusing but functional. Consider renaming to `--insert-dir` or `--s05-dir`.

---

## Architecture Concerns

### 7. Old script files not removed

Both old and new scripts coexist:
```
scripts/s07_host_map.py          # old
scripts/s04_host_map.py          # new (copy)
scripts/s09_targeted_assembly.py # old  
scripts/s05_insert_assembly.py   # new (copy)
scripts/s08_indel_detection.py   # old
scripts/s06_indel.py             # new (copy)
scripts/s10_copynumber.py        # old
scripts/s07_copynumber.py        # new (copy)
```
`run_pipeline.py` only references the new scripts, but having both is confusing. Old scripts are still tracked by git. Either:
- Delete old scripts, or
- Move to `scripts/archive/` with a note

### 8. s05 `--junctions` fallback is now dead code

`run_pipeline.py` step 5 no longer passes `--junctions`. The fallback in `s05_insert_assembly.py:3323` (`if not sites and args.junctions`) can never trigger through the pipeline. It still works if someone calls the script directly with `--junctions`, but this is an undocumented escape hatch. Either document it or remove.

### 9. Step 5 still passes `--construct-ref` from run_pipeline.py but this was a new addition

`run_pipeline.py:196` passes `--construct-ref` to step 5. This was added in the s09 refactor for construct-flanking detection but the old committed s09 does not have this flag. Since s05 is a copy of the uncommitted s09, this works, but the dependency is implicit.

---

## Minor

### 10. s07_copynumber `count_candidate_sites` parses by checking `endswith("\tCANDIDATE")`

This is fragile — any field value containing "CANDIDATE" at end of line would match. The stats format is `key\tvalue`, so `insertion_X_verdict\tCANDIDATE` works, but a more precise check would split on tab and check the value column.

---

## Summary

| # | Severity | Issue |
|---|----------|-------|
| 1 | WARNING | s11_multiqc reads from removed step directories |
| 2 | **CRASH** | viz scripts reference `s09_stats.txt` (renamed to `s05_stats.txt`) |
| 3 | BUG | s02 reads construct_ref file 3 times |
| 4 | WARNING | 10+ stale step numbers in docstrings/argparse/logs |
| 5 | WARNING | CLAUDE.md "step 8" → "step 6" |
| 6 | COSMETIC | `--s09-dir` flag name stale |
| 7 | CLEANUP | Old + new scripts both exist |
| 8 | DEAD CODE | s05 `--junctions` fallback unreachable from pipeline |
| 9 | IMPLICIT | construct-ref dependency not in committed code |
| 10 | MINOR | Fragile verdict parsing in s07_copynumber |

**Recommendation**: Fix #2 (crash) and #3 (bug) before any test run. Fix #4 docstrings in bulk with sed. Remove or archive old scripts (#7).
