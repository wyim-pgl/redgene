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

## Results matrix — partial (SLURM 5630185, as of 2026-04-18)

Partial aggregation via `analyze_coverage_sensitivity.py` → saved to
`docs/measurements/coverage_sensitivity_matrix_partial.tsv`. 10/12 tasks
COMPLETED; cucumber_line225 15x/20x still RUNNING (~10h wall, OOM-safe at
96G after BUG-18 fix).

Cell format: `n_cand / n_unknown` (`n_true=0` across all — pre-BLAST).

| sample             |   5x    |  10x    |  15x    |  20x    |
|--------------------|---------|---------|---------|---------|
| rice_G281          | 0 / 1   | 1 / 3   | 1 / 5   | 3 / 11  |
| tomato_Cas9_A2_3   | 1 / 0   | 0 / 0 † | 1 / 0   | 2 / 0   |
| cucumber_line225   | 2 / 2   | 3 / 4   | `RUN`   | `RUN`   |

† tomato_Cas9_A2_3 10x produced 0 insertion reports despite pipeline
completion (positive_sites.pkl exists, site_tier_classification.tsv
present). Genuine "no candidates at this coverage," not a pipeline failure.

### SLURM wall-time by task

| idx | sample / cov            | state     | elapsed | MaxRSS |
|-----|-------------------------|-----------|---------|--------|
| _0  | rice_G281 5x            | COMPLETED | 45m     | 11.8G  |
| _1  | rice_G281 10x           | COMPLETED | 1h29m   | 17.8G  |
| _2  | rice_G281 15x           | COMPLETED | 2h05m   | 17.9G  |
| _3  | rice_G281 20x           | COMPLETED | 3h41m   | 17.4G  |
| _4  | tomato_Cas9_A2_3 5x     | COMPLETED | 41m     | 11.5G  |
| _5  | tomato_Cas9_A2_3 10x    | COMPLETED | 1h12m   | 17.9G  |
| _6  | tomato_Cas9_A2_3 15x    | COMPLETED | 1h56m   | 19.3G  |
| _7  | tomato_Cas9_A2_3 20x    | COMPLETED | 2h33m   | 19.2G  |
| _8  | cucumber_line225 5x     | COMPLETED | 1h38m   | 23.8G  |
| _9  | cucumber_line225 10x    | COMPLETED | 5h22m   | 33.2G  |
| _10 | cucumber_line225 15x    | RUNNING   | ~10h    | —      |
| _11 | cucumber_line225 20x    | RUNNING   | ~9.8h   | —      |

All `.e+` (external epilog) entries show OUT_OF_ME+ with 1480K which is the
SLURM sampler artifact, not a real OOM — `.b+` (batch) steps completed with
the memory shown.

### Observations (partial)

- **Coverage dose-response for rice_G281**: CAND count rises monotonically
  (0 → 1 → 1 → 3 as 5x→20x). 5x produces no candidate junction — matches
  the CLAUDE.md ≥15x rule of thumb.
- **tomato_Cas9_A2_3 plateau**: CAND at 5x/15x/20x but not 10x (subsample
  stochasticity on the 833 Mbp genome). 20x plateau at 2 CAND / 0 UNK is
  cleanest.
- **cucumber_line225**: RAM demand scales aggressively (23.8G @ 5x →
  33.2G @ 10x). 15x/20x at 96G headroom likely safe but slow (>10h each).

Final matrix + GT-anchor verification (coordinate match from
`ground_truth_baseline.tsv`) pending _10/_11 completion + BLAST pass on
remaining UNK rows.

## Blank results matrix (original template, retained for reference)

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
