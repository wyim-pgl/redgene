# Resume — RedGene: algorithm hardening on `algo-improvements-2026-07`

**Date:** 2026-07-09
**Branch:** `algo-improvements-2026-07` — **not merged, not pushed**. Base `main` @ `f3cb0c0`;
run `git log --oneline main..HEAD` for the current commit list (a hard-coded
count here goes stale the moment this file is amended).
**Working dir:** `/data/gpfs/assoc/pgl/develop/redgene`
**Previous resume:** 2026-04-29 (Issue #4 s05 module split CLOSED)

---

## What this session did

Re-reviewed `work_implementation_plan.md` (2026-04-11) against the codebase, then
fixed the correctness and algorithm defects that a second, adversarial pass over
the Phase 1 / Phase 4 hot path turned up.

**`work_implementation_plan.md` is superseded as a code plan.** Its two code
tasks already shipped:

| Plan task | Status | Where |
|---|---|---|
| Task 1 — Filter D v2 | DONE | `scripts/s05/filter_d_altlocus.py`, `verdict.py` Rule 4 |
| Task 2 — UNKNOWN → FALSE_POSITIVE (host-only) | DONE (renamed constants) | `verdict.py` Rule 7, `unknown_to_fp_*` in `config.yaml` |

What remains in that plan is **SLURM re-runs only** (Phases 1–5, ~3–4 days wall
time). No code is blocking them. The rice reports on disk (dated 2026-04-16,
3 CANDIDATE / 32 FALSE_POSITIVE / 33 UNKNOWN) predate the Rule 7 code and would
collapse most of those 33 UNKNOWN on a re-run.

Open GitHub issue: **#7 only** (KCGP nomenclature, blocked on a spec).

---

## Commits

Code + tests (later commits amend this file and the wiki):

```
57245d4 Tests: end-to-end regression for find_softclip_junctions on a synthetic BAM
67125aa Follow-up from code review: keep the clip when the junction base is ambiguous
84ce853 Site discovery: N-free seeds, deterministic pairing, read-quality filter
44ed5ad Border scan: query both T-DNA repeats, collapse HSPs into distinct loci
2e183b7 Verdict: make canonical_triplet_min_identity actually gate Rule 1
a76f97a Filter D: log BLAST failures instead of failing silently toward CANDIDATE
```

Plan + measurements: `docs/superpowers/plans/2026-07-09-algorithm-hardening.md`.

pytest **322 + 1 skipped → 387 + 1 skipped**. Verdict snapshots: 26, zero drift,
fixtures untouched.

---

## The four defects, with evidence

1. **Border scan queried one 25-mer twice.** `RB_consensus` and `LB_consensus`
   both held `TGGCAGGATATATTGTGGTGTAAAC`, so every border was counted twice and
   the labels meant nothing. Both constructs actually carry **two distinct**
   repeats (`db/G281_construct.fa` `>rice_G281` at 6,556 / 279;
   `db/Cas9_construct.fa` `>tomato_A2_3` at 0 / 10,141). Measured: 4 rows → 2 loci
   per construct.

2. **`canonical_triplet_min_identity` was a dead knob.** Rule 1 outranks Filters
   B/C/D and was fed *every* annotated element, at any identity — while
   annotation runs at a 0.70 floor.

3. **Filter D swallowed BLAST failures.** `stderr=DEVNULL` + a bare
   `return False, …`. Note the verdict was already unaffected (Rule 4 gates on
   `construct_frac`, which is 0.0 either way) — the fix is that the run log now
   records the crash.

4. **Left-clip consensus leaked internal `N` into assembly seeds.** 2,327 `N`
   bases across rice Chr1+Chr2. `StrandAwareSeedExtender` skips `N`-containing
   reads and `_extend_right` needs an exact overlap, so those seeds could not
   extend in the direction they existed to extend.

Plus: order-dependent greedy junction pairing → globally distance-sorted
(bisect-windowed); no MAPQ/duplicate/QC-fail filter on clip anchors → added;
`cluster_window` / cluster depths hard-coded → new `--min-mapq`,
`--cluster-window`, `--min-cluster-depth`, `--min-single-depth`.

---

## RB vs LB — an OPEN question, deliberately not answered

The two repeats are labelled `TDNA_border_A` / `TDNA_border_B` in
`border_hits.tsv`, **not** RB / LB. Nothing in this repo establishes the
assignment, and the `canonical_v1` `RB-TDNA` / `LB-TDNA` entries in `db/*.fa`
(`TGAGCGTCGCAAAGGCGCTCGGTCT` / `GGCCTCGGCCTGAGAGCCAAAACAC`) appear in **no**
construct and contradict the motif the scan has always used — they look wrong.

**Next session:** settle this against an external reference (the Ti plasmid
accession the constructs derive from), then either relabel or fix the DB entries.
Do not infer T-DNA orientation from the current labels.

---

## Verification performed (no HPC runs)

- `_run_border_blast` on both real constructs → 2 loci each, 100% identity.
- Old vs new Phase 1 on `results/rice_G281/s04_host_map/rice_G281_host.bam`:
  - Chr3 ground-truth region: **identical** — 25 sites, identical seeds,
    `Chr3:16,439,674-16,439,709` recovered by both.
  - Chr1+Chr2: sites 5,864 → 5,706; `N` in seeds 2,327 → **0**. Every dropped
    site's old seed only cleared `min_clip=20` because `N`s padded it
    (`Chr1:457,108`: 24 bp seed, 8 `N`, unambiguous run 0 bp).
  - Defaults are byte-identical on real data: 0 clipped reads discarded
    (BAMs are not duplicate-marked).

---

## Next steps

1. **Merge decision** on `algo-improvements-2026-07` (branch, unpushed).
2. **Re-run s05** for rice + all species — the reports on disk predate Rule 7
   *and* every fix above. This is the plan's Phases 1–3.
3. Settle the **RB/LB assignment** (see above).
4. Deferred, highest detection ROI: **PE-discordancy in site discovery**.
   `find_softclip_junctions` uses soft clips only; discordant-pair primitives
   already exist in `s06b_junction_verify.py` but only to *count* support at
   known junctions, never to *discover*. Needs real-BAM validation.
5. Deferred: seed k-mer 15 → 21 for large plant genomes (gated on a sweep);
   a unified `run_blastn` wrapper over the 12 scattered call sites with
   inconsistent stderr handling.
6. When one side of a paired cluster now yields a sub-`min_clip` consensus the
   site is dropped (as before). Demoting it to a single-direction site would
   recover sensitivity — needs a validation run.

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"
git log --oneline main..HEAD
pytest tests/ -q                                  # 387 PASS + 1 skipped
pytest tests/test_verdict_snapshots.py -q         # 26 PASS, 0 drift
pytest tests/test_find_softclip_junctions.py -v   # Phase 1 e2e guard (needs minimap2)
less docs/superpowers/plans/2026-07-09-algorithm-hardening.md
```
