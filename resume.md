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
time). No code is blocking them.

> **Corrected 2026-07-10.** This file previously said the rice reports on disk
> were "3 CANDIDATE / 32 FALSE_POSITIVE / 33 UNKNOWN" and that a re-run would
> collapse most of the 33 UNKNOWN. That count came from globbing the 68
> `insertion_*_report.txt` files in `results/rice_G281/s05_insert_assembly/`,
> which s05 never cleans: they are three runs' worth (13 from 2026-04-11, 34 from
> 2026-04-14, 21 from 2026-04-16). The 2026-04-16 run itself found **21 sites:
> 3 CANDIDATE / 9 FALSE_POSITIVE / 9 UNKNOWN** — read `s05_stats.txt`, which is
> rewritten (`mode="w"`) each run. Rule 7's effect was therefore already present
> in April, not pending.

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

pytest **322 + 1 skipped → 401 + 1 skipped** (verified at HEAD 2026-07-10).
Verdict snapshots: 26, zero drift, fixtures untouched.

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

## RB vs LB — RESOLVED (2026-07-09)

| repeat | assignment | primary source |
|---|---|---|
| `TGACAGGATATATTGGCGGGTAAAC` | **RB** (`TDNA_RB`) | GenBank **J01826**, DEFINITION "T-DNA 3' (right) border"; COMMENT "right border at base 158" |
| `TGGCAGGATATATTGTGGTGTAAAC` | **LB** (`TDNA_LB`) | GenBank **J01825**, DEFINITION "T-DNA 5' (left) border" |

Yadav et al. 1982 PNAS 79:6322 (nopaline strain T37). Reproduced independently by
pBIN19 (U09365) and pCAMBIA-1300 (AF234296). Octopine pTi15955 (X00493) carries
neither. The scan now reports rice_G281 as RB@279 + LB@6557 and tomato_A2_3 as
LB@1 + RB@10,142, matching AF234296's feature table coordinate for coordinate.

Integration is **precise at RB, imprecise at LB** (J01825/J01826 COMMENT;
Shaw et al. 1984) — worth weighting RB junctions higher in confidence scoring.

### Two findings that fell out of this

1. **`db/G281_construct.fa`'s first record is byte-identical to AF234296** (md5
   match) — the **empty** pCAMBIA-1300 vector, not the G281 construct. Its T-DNA
   (lacZα MCS, CaMV35S → hptII → CaMV polyA) spans the numbering origin between
   LB(6,557) and RB(298); the 6,303 bp between RB and LB in increasing
   coordinates is backbone (pVS1 STA/REP, pBR322 bom/ori, aadA).
   `benchmark.md` and `results.md` do say "pCAMBIA-1300 backbone", but
   **`manuscript.md` claims rice G281 used an "exact construct sequence"** — the
   file does not support that. Sample-specific payload reaches step 5 only via
   the step-4b SPAdes contigs (`--extra-element-db`). **Needs a decision.**
2. **The `canonical_v1` RB/LB entries in the tracked element DBs are wrong**
   (`TGAGCGTCGCAAAGGCGCTCGGTCT` / `GGCCTCGGCCTGAGAGCCAAAACAC`): no
   `CAGGATATAT` core, present in no construct. In `element_db/gmo_combined_db_v2.fa`,
   `db/gmo_all_combined_db.fa`, `db/gmo_corn_combined_db.fa`. The border scan
   does not use them. **Curating them out is still OPEN** — deferred because
   job 5799803 reads those DBs at runtime.

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

## 2026-07-10/11 full re-run result (SLURM array 5803578, all 11 samples)

All 11 samples now done: 10 completed in-array; `soybean_UGT72E3` timed out at
16 h sequential, then completed via the fixed fan-out (array 5811471, 106/106).
No ground-truth regression (`compare_s05_runs.py` exit 0). **Two ground truths
recovered that April missed: soybean_UGT72E3 (CANDIDATE, see 4g) — and corn is
found but wrongly rejected (see 4e/#19).**

`soybean_UGT72E3` was resubmitted as per-site fan-out (array 5811054) — which
**exposed a concurrency bug** ([#18](https://github.com/wyim-pgl/redgene/issues/18)):
`assemble_insert`'s scratch dirs (`_mm2_r{rnd}`, `_pilon_r{rnd}`, `_ssake_r{rnd}`,
`check_host_termination`'s `_host_term.*`, `refine_with_foreign_reads`'s
`_foreign_refine`) were scoped to the shared `step_dir`, not per-site, so 25
concurrent tasks collided on the same `_mm2ext.bam` — 26 of 106 crashed (incl.
task 0, the ground-truth site), and the "completed" ones may have mapped reads
against another site's contig. **FIXED** (working tree, uncommitted): a per-site
`_scratch_{site_id}` root, plus atomic writes for the shared `_unmapped_R*.fq.gz`
cache. Guarded by `tests/test_s05_fanout_scratch_isolation.py` (3 tests). pytest
**404 + 1 skipped**. UGT72E3 re-run resubmitted with the fix as array **5811471**
(Phase 4 = **5811472**); array 5811054's output was discarded.

**The hardening did exactly what it was for.** Across the 9 samples with a clean
baseline (excluding **both** `soybean_AtYUCCA6` and `soybean_UGT72E3`, whose
April runs were truncated — 45 of 100 and 53 of 105 positives assembled before
the run stopped — so their site/UNKNOWN counts are not comparable; see 4f/4g):

| | sites | CANDIDATE | FALSE_POSITIVE | UNKNOWN |
|---|---|---|---|---|
| baseline | 113 | 17 | 19 | 77 |
| hardened | 81 | 18 | 19 | 44 |
| **delta** | **−32** | **+1** | **0** | **−33** |

Noise UNKNOWNs fall by ~40% (the N-free-seed fix), CANDIDATE rises by one
(`cucumber_line224`), FALSE_POSITIVE is flat. No CANDIDATE was lost to a code
change: `cucumber_line225` swapped one CANDIDATE for a better assembly of the
same payload (4f), and `rice_G281`'s single verdict move was UNKNOWN→FP. On the
two truncated soybean samples the re-run assembled the full positive set for the
first time, which is how UGT72E3's ground truth surfaced (4g).

**Three findings that outrank the hardening itself, all recorded above:**
- **4e** — Filter B structurally rejects `corn_ND207`'s true insertion
  (discovered 5 bp off the documented site, then killed FALSE_POSITIVE by its
  own event's border-flank). Highest-value fix on the branch.
- **4c** — 3 of 4 bp-resolution ground truths in `benchmark.md` are not recovered
  (tomato_Cas9_A2_3, corn — as CAND, soybean_UGT72E3 pending); this contradicts
  `benchmark.md` line 113.
- **4d** — the `tomato_WT` negative control now yields a transgene-positive site
  (UNKNOWN) off a 26 bp `P-Act1-rice` hit.

**Not committed yet.** Working tree carries the doc corrections plus 4 new
untracked helper scripts (`stage_hardened_rerun.sh`, `run_s05_hardened_array.sh`,
`run_s05_ugt72e3_fanout.sh`, `compare_s05_runs.py`). Decide what to commit after
UGT72E3's fan-out lands.

---

## Next steps

1. **Merge decision** on `algo-improvements-2026-07` → PR
   [#17](https://github.com/wyim-pgl/redgene/pull/17).
2. **Re-run s05 for every production sample** — the reports on disk predate the
   fixes above. This is the plan's Phases 2–3. **IN FLIGHT: SLURM array 5803578**
   (11 tasks, `%4` throttle, 16 CPU / 96G / 16 h each) writing into a fresh dated
   root `results/hardened_20260710/<sample>/`. Submitted from
   `run_s05_hardened_array.sh`; inputs staged by `stage_hardened_rerun.sh`, which
   symlinks each sample's s03 reads, s04b `contigs.fasta` and s04 host BAM into
   the override tree (an `--outdir-override` run otherwise sees an empty tree).
   `tomato_WT`, `tomato_Cas9_A2_1` and `corn_ND207` have no s04b assembly, so they
   run without `--extra-element-db`, exactly as their April runs did.

   **Rice's earlier `results/rice_G281_algo_v2/` (SLURM 5799803) is superseded**:
   it finished 4 minutes before the per-insert border-count fix (`0914323`) and
   12 minutes before the element-DB border correction (`496c328`/`6ba1f8b`). It
   is kept as historical evidence, not as the branch's result.
   What it did establish, against the correct 2026-04-16 baseline (21 sites,
   3 CANDIDATE / 9 FALSE_POSITIVE / 9 UNKNOWN — see the correction above):
   **the same 21 sites with identical site IDs**, the same 3 CANDIDATEs, and
   exactly one verdict change, `Chr11_28363144` UNKNOWN → FALSE_POSITIVE.
   Ground truth `Chr3:16,439,674` is CANDIDATE in both. Phase 1 logged
   `0 clipped reads discarded` — production confirmation that the new read
   filters are a no-op at default settings.
   The re-run exposed and fixed two more defects: `generate_report` counted the
   global `border_hits.tsv` rather than this insert's rows (the ground-truth
   report claimed 10 borders it did not own), and rice has in fact **never** had
   a genuine T-DNA border detected — zero HSPs ≥ 20 bp across all 21 inserts even
   under permissive BLAST.
3. ~~Curate the bogus `canonical_v1` RB/LB entries~~ — **DONE.** They were
   mis-sliced from AF485783.1 (pBI121): the records held bases 1–25 and
   ~14,040–14,065 instead of the annotated `complement(2454..2478)` (RB) and
   `complement(8621..8646)` (LB). Corrected in every DB plus the gitignored
   build sources, all indices rebuilt, manifest md5 regenerated, guarded by
   `tests/test_element_db_borders.py`. A clip spanning a real RB now hits
   `RB-TDNA` 100%/25 bp in classify's tier check; before, no clip could ever
   match. **Results on disk predate this — they need a re-run to benefit.**
4. **Decide what to do about `db/G281_construct.fa`** ([#23](https://github.com/wyim-pgl/redgene/issues/23)) being the empty
   pCAMBIA-1300 vector while `manuscript.md` calls it an exact construct.
   Corroborated independently: the stale `.fai` still named that record
   `vector|pCAMBIA-1300|T-DNA_backbone|AF234296.1|binary_vector|v1`.
4b. **OPEN ([#22](https://github.com/wyim-pgl/redgene/issues/22)): s05 leaves stale `insertion_<site>_report.txt` behind.** Re-running
   into a directory that already holds reports leaves the previous run's files
   in place, name-indistinguishable from the live ones.
   `results/rice_G281/s05_insert_assembly/` carries 47 such orphans across three
   April runs; the other ten samples are clean (report count == `insertion_sites`).
   Only `s05_stats.txt` is authoritative — it is rewritten each run. This already
   produced one wrong claim in this branch's own docs. Same misattribution class
   as the `border_hits.tsv` global-count bug. Not fixed yet because array 5803578
   must observe one fixed working tree; it writes into fresh directories, where
   the bug cannot bite. Fix after the array drains: clear prior
   `insertion_*_report.txt` at the start of Phase 5, or emit a run manifest.

   The per-insert border fix is **confirmed in production**: all 21 of April's
   2026-04-16 rice reports printed *"T-DNA borders found: 10"*, and the 34 from
   2026-04-14 printed *"2"* — the same global count stamped on every site. The
   re-run's `border_hits.tsv` is empty and no report claims a border. While
   there: `report.py:320` guards that line with `if n_borders > 0`, so "no border
   detected" and "border never checked" read identically in the report. A
   regulatory deliverable should print the explicit `0`.
4c. **OPEN ([#20](https://github.com/wyim-pgl/redgene/issues/20)): three of the four
   bp-resolution ground truths in `benchmark.md` are not recovered at all.**
   Measured against the April baselines (`compare_s05_runs.py`, ±1 kb window):

   | sample | documented site | April baseline | 2026-07-10 re-run |
   |---|---|---|---|
   | `rice_G281` | `Chr3:16,439,674` | **CANDIDATE**, 0 bp off | CANDIDATE |
   | `tomato_Cas9_A2_3` | `SLM_r2.0ch08:65,107,378` | **not discovered** — 0 of its 2 sites are on ch08 (they are ch01, ch06) | still not discovered |
   | `corn_ND207` | `NC_050098.1:181,367,276` (chr3) | **not discovered** — its only site is on chr3 but 53.97 Mb away | **discovered at 181,367,271 (5 bp off) — and rejected FALSE_POSITIVE, see 4e** |
   | `soybean_UGT72E3` | `NC_038254.2:51,882,860` (Gm18) | **not discovered** — 0 of its 53 sites are on Gm18 | **RECOVERED as CANDIDATE** at `51,882,867` (7 bp off) — 30,496 bp multi-copy insert, host 73%, see 4g |

   The soybean coordinate is confirmed to be in *our* reference's coordinate
   system, not the source paper's: `db/Gmax_v4.0.gff3.gz` puts `LOC102669442`
   (= Glyma.18g226800) at `NC_038254.2:51,789,501-51,906,722`, and 51,882,860
   falls inside it. So "not discovered" there is solid. For corn and tomato the
   coordinates' assembly provenance is *not* independently verified (no corn GFF
   in `db/`), so treat those two as strongly suggestive rather than settled.

   This **contradicts `benchmark.md` line 113**, which tabulates tomato_Cas9_A2_3
   as RedGene *"Yes | chr08:65,107,378 (T-DNA junction) | TP 1 | FP 0"*. Note
   65,107,378 is not a gRNA on-target either — gRNA 2's on-target is
   `SLM_r2.0ch08:53,314,226`. Same class of unsupported claim as the
   `manuscript.md` "exact construct sequence" problem in item 4.

   Separately, `tomato_Cas9_A2_1` yields **zero** transgene-positive sites in both
   April and the 2026-07-10 re-run: Phase 1 finds 6,293 junctions, 0 with any
   element hit (6,261 no hit, 32 `univec`). A Cas9 line with no detectable T-DNA
   needs explaining — segregant, mislabelled sample, or a real recall failure.

   None of this is caused by the hardening — the baselines already showed it.
   It does mean the branch must not be reported as validated on more than rice
   until these are resolved. Re-check with `python compare_s05_runs.py` once
   array 5803578 drains; it fails (exit 1) only on a ground truth that *was*
   CANDIDATE and no longer is.
4d. **WATCH ([#21](https://github.com/wyim-pgl/redgene/issues/21)): the WT negative control now yields a transgene-positive site.**
   In the 2026-07-10 re-run, `tomato_WT` produces one assembled site,
   `SLM_r2.0ch08:48,294,203`, verdict **UNKNOWN**. Its 5' clip (74 bp) hits
   `P-Act1-rice` (rice actin1 promoter, S44221.1) at **96.2% over 26 bp** in
   classify's `blastn-short` tier check. April's `tomato_WT` run assembled
   nothing. This is the project's top documented pitfall — host-derived promoters
   in the element DB — reappearing in a wild-type control.

   Two mitigations exist and neither applies here:
   - The **host-endogenous exclusion** BLASTs the whole element against the host
     genome; a 26 bp spurious hit is far below the threshold that would retire
     `P-Act1-rice` for tomato.
   - The **mask BED** (`docs/host_masks/tomato_slm_r2.bed`) is a **0-byte file**,
     so `--mask-bed` is a no-op for tomato. That is *consistent* with
     `host_masked_rationale.tsv`, which has **no tomato rows at all**: mask
     regions carry code `X-01` (or `E-01`), while `E-0*-EMPTY` records an element
     analysed with nothing to mask. rice has 7 `X-01` → 7 BED lines; corn 12
     `X-01` + 1 `E-01` → 13; soybean has only `E-03/E-04-EMPTY` → empty BED, by
     design. Whether tomato was ever analysed is simply not recorded.
     Note `P-Act1-rice` has a rationale row for the **rice** host only
     (`E-02-EMPTY`) — it was never evaluated against the tomato genome.

   **Do not yet blame this branch.** April's `tomato_WT` s05 ran on 2026-04-13
   code, well before `--mask-bed`, the tier check and the common-payload DB, so
   it is not a clean control for attribution. The verdict layer did contain the
   site (UNKNOWN, not CANDIDATE). Confirm against the other re-run samples, then
   decide whether tomato needs a mask BED.
4g. **NEW DETECTION WIN: soybean_UGT72E3's ground truth is recovered.** With the
   fan-out fix ([#18](https://github.com/wyim-pgl/redgene/issues/18)) the full
   106-site assembly completed (array 5811471, 106/106, 0 fail; Phase 4 re-run
   separately as job 5811947 after the wrapper's micromamba activation failed —
   fixed in `submit_s05_array.sh`). The site `NC_038254.2:51,882,867` — 7 bp from
   the documented `51,882,860`, inside `LOC102669442` (= Glyma.18g226800) — is a
   **30,496 bp multi-copy insert, host 73%, all FP filters passed → CANDIDATE**.
   April found nothing on Gm18. This is a direct recall gain from the hardening,
   and it means UGT72E3 must NOT be listed under #20's recall gaps.

4f. **The hardening does what it was for: it strips noise UNKNOWNs while keeping
   the CANDIDATE signal.** Per-sample re-run vs the true April baselines
   (`compare_s05_runs.py`, `s05_stats.txt` on both sides):

   | sample | sites | UNKNOWN | CANDIDATE |
   |---|---|---|---|
   | `rice_G281` | 21 → 21 | 9 → 8 | 3 → 3 (same) |
   | `tomato_Cas9_A2_2` | 2 → 2 | 0 → 0 | 2 → 2 (same) |
   | `tomato_Cas9_A2_3` | 2 → 2 | 0 → 0 | 2 → 2 (same) |
   | `cucumber_line212` | 31 → 17 | 26 → 13 | 1 → 1 (same site) |
   | `cucumber_line224` | 18 → 10 | 15 → 6 | 1 → 2 |
   | `cucumber_line225` | 38 → 26 | 26 → 15 | 8 → 8 (7 same, 1 swapped) |

   The cucumber drops are the N-free-seed fix measured in production: the lost
   sites were UNKNOWNs whose seeds only cleared `min_clip` because `N`s padded
   them (the same mechanism attributed on rice Chr1+Chr2). No CANDIDATE was lost
   to the seed change.

   `cucumber_line224` **gains** a CANDIDATE (`LKUO03001518.1_2091613`,
   UNKNOWN → CANDIDATE). `cucumber_line225` keeps 8 CANDIDATEs but swaps one:
   `LKUO03000939.1_800850` (baseline, 601 bp insert, 464 bp payload) is no longer
   discovered, and `LKUO03001728.1_511594` (new, **11,600 bp insert, 11,277 bp
   foreign gap, host 1.5%**) appears. Different contigs — not a relocation — but
   both annotate the same payload `NODE_1_length_5568`, consistent with the
   documented 2-copy insertion; the new one is the far more complete assembly.
   These are net detection *gains*, produced by the deterministic pairing + N-free
   seed changes, not regressions.

   **`soybean_AtYUCCA6` is the one sample that grows (45 → 99 sites, CAND 0 → 3),
   and it is a baseline artefact, not a hardening effect.** Both runs identify
   ~the same transgene-positive set — baseline `positive_sites.json` = 100,
   hardened = 99 — but the April run wrote only **45** reports / stats rows before
   stopping (all dated 2026-04-14; the other 55 positives were never assembled,
   i.e. the April run was truncated). The re-run assembled all 99. So the "growth"
   is the April baseline having been partial, the same accumulation-vs-truncation
   trap as the rice 68-file count. The 3 new CANDIDATEs are large clean inserts
   (host 2% / 30% / 85%; foreign gaps 13 kb / 26 kb / 1 kb). One,
   `NC_038254.2:52,710,589`, is on Gm18 but inside `LOC100813760` (RIC1), **not**
   the documented Site I gene Glyma.18g235800 — so it is not a ground-truth
   confirmation. `soybean_UGT72E3` (still running) must be read the same way:
   check `positive_sites.json`, not just the report count.

4e. **HIGHEST-VALUE BUG ([#19](https://github.com/wyim-pgl/redgene/issues/19)): Filter B is guaranteed to reject
   corn_ND207's true insertion, because the construct DB contains that event's
   own genomic flanks.**

   The re-run newly discovers `NC_050098.1:181,367,271` — **5 bp** from
   `benchmark.md`'s documented `181,367,276`. Its report is a textbook true
   positive: 833 bp insert, `5' host --[QL-CON-00-015 (Bt event)]--[P-CaMV35S]--
   [ND207-LB]-- 3' host`, host fraction **18%**, largest non-host gap **683 bp**.
   Filter A would pass it. It is called **FALSE_POSITIVE** by Filter B:
   *"overlaps construct-flanking region NC_050098.1:181,366,769-181,367,847 —
   host DNA in construct reference"*.

   Mechanism, verified by BLAST:

   | record in `db/gmo_corn_combined_db.fa` | length | aligns to maize | identity |
   |---|---|---|---|
   | `ND207-LB` | 139 bp | `NC_050098.1:181,367,269-181,367,347` (q 1-79) | 98.7% |
   | `ND207-RB` | 183 bp | `NC_050098.1:181,367,141-181,367,242` (q 82-183) | 99.0% |

   Both border records carry ~50-60% **native maize flanking taken from the
   ND207 junction itself** (Sci Rep 2025 Table S1). Filter B derives its
   "construct reference contains host DNA" region by BLASTing the construct DB
   against the host — so that region *is* the insertion locus, and any site
   landing there is rejected by construction. CLAUDE.md already warns "Border
   sequences contain ~100bp native flanking"; what was missed is that this makes
   Filter B **anti-correlated with truth** for such records. Verdict priority
   puts Filter B above Filter A, so the strong host-fraction evidence cannot
   rescue it, and the `[bar, P-CaMV35S, T-ocs]` canonical triplet does not
   complete (only `P-CaMV35S` is present), so Rule 1 cannot either.

   The same flanking also causes `ND207-LB` to be retired as
   *"Host-endogenous (excluded at DB-level BLAST)"*, so the annotation loses the
   border evidence too.

   Fix direction (needs a decision, not yet implemented): trim the native flank
   off the 62 corn LB/RB border records, or exclude records tagged `corn_border`
   from Filter B's construct-vs-host BLAST and from the host-endogenous
   exclusion. Either way, re-run corn and re-check that
   `NC_050098.1:181,367,271` promotes to CANDIDATE. This is the single clearest
   detection win available on this branch.
5. Deferred, highest detection ROI: **PE-discordancy in site discovery**.
   `find_softclip_junctions` uses soft clips only; discordant-pair primitives
   already exist in `s06b_junction_verify.py` but only to *count* support at
   known junctions, never to *discover*. Needs real-BAM validation.
6. Deferred: seed k-mer 15 → 21 for large plant genomes (gated on a sweep);
   a unified `run_blastn` wrapper over the 12 scattered call sites with
   inconsistent stderr handling.
7. When one side of a paired cluster now yields a sub-`min_clip` consensus the
   site is dropped (as before). Demoting it to a single-direction site would
   recover sensitivity — needs a validation run.

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"
git log --oneline main..HEAD
pytest tests/ -q                                  # 401 PASS + 1 skipped
pytest tests/test_verdict_snapshots.py -q         # 26 PASS, 0 drift
pytest tests/test_find_softclip_junctions.py -v   # Phase 1 e2e guard (needs minimap2)
less docs/superpowers/plans/2026-07-09-algorithm-hardening.md
```
