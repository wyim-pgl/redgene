# s04 minimap2 PoC - rice_G281 (T12)

**Date submitted:** 2026-04-16
**Operator:** wyim
**Branch:** `feature/v1.0-mvp-2026-04-16`
**HEAD at submit:** <fill after commit>
**SLURM JOBID:** <fill after submit>
**Partition/Account:** cpu-s2-core-0 / cpu-s2-pgl-0

## Goal

Validate whether s04 host mapping can be migrated from BWA to minimap2 `-ax sr`
without regressing the v1.0 acceptance criteria, per won_yim's conditional
approval (team-consensus.md §5 "s04 migration timing" and §8 Unresolved Item #1).

The canonical BWA output tree in `results/rice_G281/` is preserved;
the PoC writes to `results/rice_G281_mm2_poc/` so both can be compared.

## Experiment design

1. Build minimap2 host index (`db/Osativa_323_v7.0.fa.mmi`, ~2-3 min, cached).
2. Run `minimap2 -ax sr -t 16` on the same s01-QC reads that feed the BWA path.
3. Pipe through `samtools sort`, index, record wall-time to
   `results/rice_G281_mm2_poc/rice_G281/s04_host_map/wall_time_seconds.txt`.
4. Symlink s03 outputs from the BWA baseline into the PoC tree so s04b
   (SPAdes) and s05 (insert assembly) run against identical upstream state.
5. Invoke `run_pipeline.py --steps 4b,5 --host-bam-override ... --outdir-override ...`
   to run the downstream steps with the minimap2 BAM swapped in.

Only s04 differs between the two runs; every other input is byte-identical.

## Acceptance criteria

All four must PASS for v1.1 adoption. Any FAIL defers minimap2 to v2.0 grant.

| # | Criterion | Target | Result | Status |
|---|---|---|---|---|
| 1 | Chr3:16,439,674 verdict preserved | CANDIDATE or CONFIRMED (matching BWA) | <fill> | <PASS/FAIL> |
| 2 | Phase 1 transgene-positive site count | ±20% of BWA count | <fill> | <PASS/FAIL> |
| 3 | s04 wall time reduction | ≥30% vs BWA baseline (5h lower bound) | <fill> | <PASS/FAIL> |
| 4 | MAPQ<20 soft-clip read ratio | ≤ BWA ratio + 0.05 (BUG-7 guard) | <fill> | <PASS/FAIL> |

## Decision matrix

- **4/4 PASS:** minimap2 `-ax sr` approved for rice_G281 in v1.1. Remaining
  four hosts (tomato / cucumber / corn / soybean) still require individual
  PoC runs before full migration — run this script with the host reference
  swapped and a different `OUT_ROOT`.
- **Any FAIL:** BWA path retained for v1.0 and v1.1. minimap2 migration
  re-evaluated in the v2.0 grant proposal (haibao Level 2/3 scope).

## Raw artifacts

- `results/rice_G281_mm2_poc/rice_G281/s04_host_map/rice_G281_host.bam` - minimap2 BAM
- `results/rice_G281_mm2_poc/rice_G281/s04_host_map/mm2.log` - minimap2 stderr
- `results/rice_G281_mm2_poc/rice_G281/s04_host_map/wall_time_seconds.txt` - s04 wall time
- `results/rice_G281_mm2_poc/rice_G281/s04b_construct_asm/` - SPAdes per-sample DB
- `results/rice_G281_mm2_poc/rice_G281/s05_insert_assembly/` - s05 verdicts + tier TSV
- `results/rice_G281_mm2_poc/pipeline.log` - full run_pipeline.py stdout
- `results/rice_G281/` - BWA baseline (for comparison)
- `results/rg_s04_mm2_poc_<JOBID>.out` - SLURM log
- `results/rg_s04_mm2_poc_<JOBID>.err` - SLURM stderr

## Running the analysis

After SLURM completion (~4-6h):

```bash
bash docs/measurements/s04_minimap2_poc_analyze.sh <JOBID>
```

This prints all four criterion results side-by-side with BWA.
Transcribe the numbers into the table above.

## Decision

**Decided by:** <operator name>
**Decision date:** <fill>
**Verdict:** <ACCEPT-v1.1 / DEFER-v2.0>
**Rationale:** <1-2 sentences>
