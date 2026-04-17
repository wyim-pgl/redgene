# AC-7: Coverage sensitivity sweep (3 hosts x 4 coverages)

Issue tracker: [#2](https://github.com/wyim-pgl/redgene/issues/2)

## Purpose

Empirically validate the rule-of-thumb coverage thresholds documented in
`CLAUDE.md` (rice ≥15x, tomato ≥10x, cucumber ≥10x) by running the core
pipeline at 5x / 10x / 15x / 20x for each of three flagship samples and
measuring (a) GT anchor recall and (b) false-positive count.

## Operator SOP

```bash
# (1) Ensure master FASTQs are staged at input/<sample>_R[12].fq.gz
ls input/rice_G281_R1.fq.gz input/tomato_A2_3_R1.fq.gz input/cucumber_line225_R1.fq.gz

# (2) Edit run_coverage_sensitivity.sh if you need different coverage fractions
#     (the defaults assume ~20x master reads → 0.25/0.50/0.75/1.00).

# (3) Submit the 12-task SLURM array — operator only, no auto-submit here
sbatch run_coverage_sensitivity.sh

# (4) After all 12 tasks succeed, collect the summary TSV
python scripts/util/analyze_coverage_sensitivity.py \
    --results-dir results \
    --samples rice_G281 tomato_A2_3 cucumber_line225 \
    --coverages 5x 10x 15x 20x \
    --out docs/measurements/coverage_sensitivity_summary.tsv
```

## Blank results matrix (fill in post-sweep)

Leave this table as-is until the 12 tasks complete. The analyzer script
produces a one-row-per-cell TSV that can be pasted into the matrix below.

| sample           |  5x  |  10x |  15x |  20x |
|------------------|------|------|------|------|
| rice_G281        | `—`  | `—`  | `—`  | `—`  |
| tomato_A2_3      | `—`  | `—`  | `—`  | `—`  |
| cucumber_line225 | `—`  | `—`  | `—`  | `—`  |

Cell format: `n_true / n_cand / n_unknown (gt_hit?)`

## Scaffold status (as of this commit)

- [x] `scripts/util/subsample_reads.py` — fraction-first, also supports
  `--target-coverage` + genome-size math.
- [x] `scripts/util/analyze_coverage_sensitivity.py` — tallies verdict mix
  per (sample, coverage) cell.
- [x] `run_coverage_sensitivity.sh` — SLURM array template (12 tasks,
  `--mem=96G --cpus-per-task=16 --time=24:00:00`, explicit headers).
- [ ] **Pending operator work:** add `<sample>_cov<tag>` entries to
  `config.yaml` so `run_pipeline.py` can find the subsampled inputs. This is
  explicitly out of scope for the weekend scaffold (user request: no
  config.yaml mutation).

## Risks / notes

- cucumber at 5x is expected to fail entirely — the literature shows
  assembly-based junction detection breaks below ~8x for complex genomes.
  The sweep will still run and log the zero-result cell to quantify it.
- Seqkit random sampling is seeded (`--seed 42`) for reproducibility but
  two runs with different seeds would need to be averaged before publishing
  the threshold curves (v1.1 enhancement).
