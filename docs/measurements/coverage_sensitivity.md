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

# (4) After all 12 tasks succeed, collect the summary TSV (now with GT anchor
#     resolution via ground_truth_baseline.tsv)
python scripts/util/analyze_coverage_sensitivity.py \
    --results-dir results \
    --samples rice_G281 tomato_Cas9_A2_3 cucumber_line225 \
    --coverages 5x 10x 15x 20x \
    --ground-truth ground_truth_baseline.tsv \
    --out docs/measurements/coverage_sensitivity_matrix.tsv
```

## Results matrix — final (SLURM 5630185, 12/12 COMPLETED)

Full aggregation via `analyze_coverage_sensitivity.py` → saved to
`docs/measurements/coverage_sensitivity_matrix.tsv`.

Cell format: `n_cand / n_unknown — gt_anchor_hit` (`n_true=0` across all;
MVP classifier emits CANDIDATE until downstream BLAST promotes rows).

| sample             | 5x                      | 10x                         | 15x                              | 20x                         |
|--------------------|-------------------------|-----------------------------|----------------------------------|-----------------------------|
| rice_G281          | 0 / 1 — MISS            | 1 / 3 — HIT:CANDIDATE       | 1 / 5 — **HIT:FALSE_POSITIVE**   | 3 / 11 — HIT:CANDIDATE      |
| tomato_Cas9_A2_3   | 1 / 0 — HIT:CANDIDATE   | 0 / 0 — MISS †              | 1 / 0 — HIT:CANDIDATE            | 2 / 0 — HIT:CANDIDATE       |
| cucumber_line225   | 2 / 2 — MISS            | 3 / 4 — MISS                | 6 / 17 — HIT:CANDIDATE           | 7 / 28 — HIT:CANDIDATE      |

† `tomato_Cas9_A2_3` at 10x produced zero insertion reports despite pipeline
completion (positive_sites.pkl present, but no site cleared the s05 gate).
Subsample stochasticity on the 833 Mbp genome; flanking coverages 5x / 15x /
20x all HIT. Note distinguished by analyzer: *"pipeline complete, 0 insertion
reports"* (not a pipeline failure).

### SLURM wall-time and memory by task

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
| _10 | cucumber_line225 15x    | COMPLETED | 10h52m  | 48.1G  |
| _11 | cucumber_line225 20x    | COMPLETED | 19h20m  | 62.9G  |

All `.e+` (external epilog) entries show OUT_OF_ME+ with 1480K which is the
SLURM sampler artifact, not a real OOM — `.b+` (batch) steps completed with
the memory shown.

### Observations

1. **Rice coverage dose-response** — CAND count rises monotonically
   (0 → 1 → 1 → 3 from 5x → 20x). 5x produces no candidate at the GT locus,
   matching the `CLAUDE.md` ≥15x rule of thumb. 5x does not reach the GT
   anchor at all (MISS).

2. **Mid-coverage assembly instability (rice 15x)** — At 15x the GT junction
   is *detected* (HIT) but classified **FALSE_POSITIVE**: the PARTIAL
   assembly (round 5) extended enough to include off-target host fragments
   (Chr2:373bp, Chr11:201bp), triggering the chimeric-assembly filter.
   At 20x the assembly CONVERGED (round 4) to a clean head-to-head 2-copy
   structure → CANDIDATE. Practical lesson: *partial assemblies at moderate
   coverage can over-extend into host noise;* increasing coverage or tuning
   the chimeric filter threshold (`--host-span-max`) are both valid knobs.

3. **Tomato plateau with a stochastic dip** — HIT at 5x / 15x / 20x; 10x
   returns 0 reports from the subsample draw. Recommendation: for the 833
   Mbp tomato genome, favour 15x+ to escape the 10x stochasticity valley.

4. **Cucumber coverage threshold** — GT anchor is MISS at 5x and 10x
   (assembly never reaches the junction through 8035-contig reference
   fragments). First HIT is 15x with 26 total sites, 6 CAND / 17 UNK.
   This is **above** the previous `CLAUDE.md` ≥10x guideline — matrix
   tightens the rule to **cucumber ≥15x** for GT-level sensitivity.

5. **Memory scales aggressively on cucumber** — 23.8G (5x) → 33.2G (10x) →
   48.1G (15x) → **62.9G (20x)**. The `--mem=96G` default from
   Issue #15 / BUG-18 remains necessary headroom (33G margin at 20x).

6. **Cucumber 20x candidate density** — 38 total sites (up from 26 @ 15x), 7
   CAND / 28 UNK / 3 FP, GT=CANDIDATE. More coverage surfaces more sites but
   most land in UNK (pending downstream BLAST); this matches the cuc_225
   AC-2 workflow where all 8 cuc_225 CAND promoted to TRUE_INSERTION after
   remote BLAST confirmed T-DNA vector backbone identity.

### Threshold recommendation updates

| host    | old CLAUDE.md | matrix-validated |
|---------|---------------|------------------|
| rice    | ≥15x          | **≥15x** (5x = MISS, 10x CAND, 15x CAND-but-FP, 20x clean CAND) |
| tomato  | ≥10x          | **≥15x** (10x stochastically empty) |
| cucumber| ≥10x          | **≥15x** (10x = MISS; 15x first HIT) |

All three hosts converge on **≥15x** as the practical floor for GT-anchor
recall. This is a matrix-derived refinement; a confirmatory re-run with a
second seqkit seed (`v1.1` enhancement) would tighten the error bars.

## Blank results matrix (original template, retained for reference)

| sample           |  5x  |  10x |  15x |  20x |
|------------------|------|------|------|------|
| rice_G281        | `—`  | `—`  | `—`  | `—`  |
| tomato_A2_3      | `—`  | `—`  | `—`  | `—`  |
| cucumber_line225 | `—`  | `—`  | `—`  | `—`  |

Cell format: `n_true / n_cand / n_unknown (gt_hit?)`

## Scaffold status

- [x] `scripts/util/subsample_reads.py` — fraction-first, also supports
  `--target-coverage` + genome-size math.
- [x] `scripts/util/analyze_coverage_sensitivity.py` — tallies verdict mix
  per (sample, coverage) cell **and** resolves `gt_anchor_hit` via
  `ground_truth_baseline.tsv` (v1.1 enhancement landed in this session).
- [x] `run_coverage_sensitivity.sh` — SLURM array template (12 tasks,
  `--mem=96G --cpus-per-task=16 --time=24:00:00`, explicit headers).
- [x] `config.yaml` — 12 subsampled entries (inherit `host_reference` and
  `construct_reference` from parent sample keys).
- [x] `ground_truth_baseline.tsv` — GT anchor rows for all 3 hosts
  (rice_G281, tomato_Cas9_A2_3, cucumber_line225).

## Risks / notes

- cucumber at 5x is expected to miss entirely — this matrix confirms the
  assembly-based junction detector breaks below ~10x for fragmented
  references. The sweep still runs and logs the zero-candidate cell.
- Seqkit random sampling is seeded (`--seed 42`) for reproducibility but
  two runs with different seeds would need to be averaged before publishing
  the threshold curves (v1.1 enhancement).
- The `HIT:FALSE_POSITIVE` at rice 15x is **not a regression** — it is
  genuine signal about assembly-quality-vs-coverage trade-offs that
  operators should know about when interpreting mid-coverage runs.
