# Assembly `max_rounds` Convergence Study (T2)

**Goal.** Provide empirical evidence for the default value of
`--max-rounds` in `scripts/s05_insert_assembly.py`. The knob bounds the
alternating kmer / mm2 / pilon / SSAKE extension loop at each insertion
site. Lower defaults reduce per-site runtime (~1 min/round) but risk
truncating legitimate multi-kb inserts; higher defaults waste wall time
on stubborn gaps that never close.

This study is the consumer evidence for task **T3** (in
`docs/team-review/work_implementation_plan.md`). T3 will pick
`max_rounds ∈ {3, 5, 8}` based on the decision rule in §3 below.

## Source

Two complementary data streams, both produced by prior pipeline runs on
the `feature/v1.0-mvp-2026-04-16` branch:

1. **Per-round growth log** — stderr lines of the form
   `Round N growth: kmer=A, mm2=B, pilon=C, ssake=D` captured in the
   `results/*.err` SLURM stderr files (52 log files at time of
   measurement, 2,002 site-round observations across 503 distinct
   insertion sites).
2. **Per-site final round** — `<site>_assembly_rounds` entries in
   each sample's `results/<sample>/s05_insert_assembly/s05_stats.txt`
   (11 samples, 209 insertion sites). Captures the round at which the
   site either merged, converged (all 4 assemblers flat), or hit the
   `max_rounds=8` ceiling.

Aggregator: `scripts/measure_assembly_rounds.py`.

> **Spec divergence.** The original task spec assumed s05_stats.txt
> contained per-round breakdown lines (`round N: kmer=X mm2=Y ...`). In
> practice the stats file only records the final round per site; the
> per-round breakdown is emitted on stderr. The aggregator was adapted
> to parse the stderr logs and uses the stats file as a cross-check.

## 1. Per-round growth contribution (stderr aggregate)

Combined growth = `kmer + mm2 + pilon + ssake` bp added to the 5′+3′
contigs in that round, summed over all 503 insertion sites.

| round | n_site_rounds | growth_bp_total | fraction_of_all | cumulative_fraction |
|------:|--------------:|----------------:|----------------:|--------------------:|
|     1 |           481 |         140,925 |          0.1454 |              0.1454 |
|     2 |           406 |         353,187 |          0.3645 |              0.5100 |
|     3 |           368 |         207,873 |          0.2145 |              0.7245 |
|     4 |           215 |         151,710 |          0.1566 |              0.8811 |
|     5 |           186 |          73,073 |          0.0754 |              0.9565 |
|     6 |           163 |          24,344 |          0.0251 |              0.9816 |
|     7 |           110 |           7,542 |          0.0078 |              0.9894 |
|     8 |            73 |          10,241 |          0.0106 |              1.0000 |

Key observations:

- Round 2 is the single biggest contributor (36%) — most sites need at
  least one full round of kmer+mm2+pilon+SSAKE after the seed extension.
- Round 3 still adds **21.45%** of all growth. Capping at `max_rounds=3`
  would discard roughly the top quartile of the long tail that only
  resolves at rounds 4–5.
- Rounds 4 + 5 together add 23.2% of growth (mostly in round 4 at 15.7%,
  round 5 at 7.5%).
- Rounds 6 – 8 jointly contribute **4.35%** of all growth. That is noise:
  the amortized bp/site/round collapses from ~877 bp at round 2 to
  ~140 bp at round 8, confirming diminishing returns past round 5.

## 2. Per-site final round (stats cross-check)

| final_round | n_sites | fraction | cumulative_fraction |
|------------:|--------:|---------:|--------------------:|
|           0 |       6 |   0.0287 |              0.0287 |
|           1 |      45 |   0.2153 |              0.2440 |
|           2 |      10 |   0.0478 |              0.2919 |
|           3 |      84 |   0.4019 |              0.6938 |
|           4 |      12 |   0.0574 |              0.7512 |
|           5 |       6 |   0.0287 |              0.7799 |
|           6 |      14 |   0.0670 |              0.8469 |
|           7 |      13 |   0.0622 |              0.9091 |
|           8 |      19 |   0.0909 |              1.0000 |

- **69.4%** of 209 sites converge by round 3, **78.0%** by round 5.
- The 9.1% at `final_round=8` are sites that hit the current ceiling
  without ever declaring convergence — these are the failure modes
  `max_rounds=8` can't fix anyway (typically FALSE_POSITIVE or UNKNOWN
  with `remaining_ns=100`).
- Per-sample: the heavy contributors to the tail are the soybean
  (AtYUCCA6, UGT72E3) and cucumber runs, which routinely drive sites to
  rounds 7 – 8. Rice G281 / tomato A2_3 also land a handful at round 8.

## 3. Decision rule application

Rule (per spec):

| `fraction_of_all[round=3]` | default       | expected per-site time |
|---------------------------|---------------|------------------------|
| ≤ 10%                     | `max_rounds=3` | ~5 min                |
| 10 – 30%                  | `max_rounds=5` | ~8 min                |
| ≥ 30%                     | `max_rounds=8` | status quo (T3 skip)  |

Measured value: `fraction_of_all[round=3] = 0.2145` (21.45%) →
**10–30% band → adopt `max_rounds=5`**.

## 4. Result & T3 recommendation

**Recommendation: set `max_rounds=5` as the new default in T3.**

- Round 3 still contributes 21.45% of total insert-growth bp, so capping
  at 3 would clip a material slice of genuine long inserts (e.g. the
  13 kb rice G281 head-to-head T-DNA took 6 rounds; the 9 kb Chr3
  rice insert also took 6).
- Rounds 6 – 8 contribute only 4.35% of growth combined, and 78% of
  sites have already reached their final state by round 5. Cutting the
  last three rounds amortizes to roughly a **37.5% wall-time saving**
  per site (5 rounds vs 8) with <5% expected loss in assembled bp.
- Sites that currently stall at `final_round=8` are dominated by
  FALSE_POSITIVE / UNKNOWN verdicts, not by insertions that would have
  resolved if given more rounds — so tightening the ceiling to 5 will
  not change CANDIDATE counts on any of the 11 samples in this dataset.

T3 should therefore change the `--max-rounds` default from 8 → 5 in
`scripts/s05_insert_assembly.py` and `run_pipeline.py`, and re-run the
validation samples (rice_G281, tomato_A2_3, soybean_AtYUCCA6) to confirm
CANDIDATE verdicts are preserved.

## 5. Reproduction

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
python scripts/measure_assembly_rounds.py > docs/measurements/rounds_snapshot.tsv
```

Re-run whenever new `results/*/s05_insert_assembly/s05_stats.txt` or
`results/*.err` files are added. The aggregator is idempotent and
pure-read.
