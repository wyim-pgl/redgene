# RedGene — Algorithm Hardening (2026-07-09)

**Context.** `work_implementation_plan.md` (2026-04-11) is largely superseded: its
code tasks (Filter D v2, UNKNOWN→FALSE_POSITIVE reclassification) landed as
`scripts/s05/filter_d_altlocus.py` and `compute_verdict` Rule 4 / Rule 7. What
remains there is HPC re-runs (days of SLURM), not code.

This plan covers the code work that a second, adversarial read of the hot path
turned up: four correctness defects and three algorithmic improvements, all
unit-testable without HPC.

Baseline: `pytest tests/ -q` → **322 passed, 1 skipped**.

---

## A. Correctness defects

### A1 — T-DNA border scan queries the same 25-mer twice
`scripts/s05/annotation.py:323-324` writes `RB_consensus` and `LB_consensus`
with the byte-identical sequence `TGGCAGGATATATTGTGGTGTAAAC`.

Measured evidence (this repo, record-scoped scan of the first FASTA record):

| construct | length | border repeats found |
|---|---|---|
| `db/G281_construct.fa` (`>rice_G281`) | 8,958 bp | pos 278 `TGACAGGATATATTGGCGGGTAAAC`, pos 6,556 `TGGCAGGATATATTGTGGTGTAAAC` |
| `db/Cas9_construct.fa` (`>tomato_A2_3`) | 16,419 bp | pos 0 `TGGCAGGATATATTGTGGTGTAAAC`, pos 10,141 `TGACAGGATATATTGGCGGGTAAAC` |

So each construct carries **two distinct** 25-bp imperfect direct repeats, and
the scan queries only one of them — twice. Consequences:

1. every genuine border yields two identical BLAST rows → `T-DNA borders found:`
   (`scripts/s05/report.py:253`, a raw line count) is inflated ≈2×;
2. the `RB_consensus` / `LB_consensus` labels carry zero information;
3. the second repeat variant is never queried under its own identity.

**Fix.** Query both distinct repeats; collapse overlapping HSPs on the insert
into distinct border loci before writing `border_hits.tsv`. Extract the scan
into `_run_border_blast` (also closes the deferred item in `resume.md`).

**RB/LB assignment — RESOLVED 2026-07-09 against primary deposits.** Labels are
`TDNA_RB` / `TDNA_LB`:

| repeat | assignment | source |
|---|---|---|
| `TGACAGGATATATTGGCGGGTAAAC` | **RB** | GenBank J01826, DEFINITION "T-DNA 3' (right) border"; COMMENT "right border at base 158" (the 25-mer sits at 158–182) |
| `TGGCAGGATATATTGTGGTGTAAAC` | **LB** | GenBank J01825, DEFINITION "T-DNA 5' (left) border" |

Both from Yadav, Vanderleyden, Bennett, Barnes & Chilton (1982) PNAS 79:6322,
nopaline strain T37. Independently reproduced by pBIN19 (U09365) and
pCAMBIA-1300 (AF234296), which annotate `left border repeat from C58 T-DNA` at
6557..6582 and `right border T-DNA repeat` at 298..323. The octopine T-DNA
(X00493) contains neither repeat.

Two consequences fell out of establishing this:

1. **`db/G281_construct.fa`'s first record is byte-identical to AF234296** —
   the *empty* pCAMBIA-1300 vector, not the G281 construct. Its T-DNA (lacZα
   MCS, CaMV35S → hptII → CaMV polyA) spans the numbering origin between
   LB(6,557) and RB(298); the 6,303 bp between RB and LB in increasing
   coordinates is the backbone (pVS1 STA/REP, pBR322 bom/ori, aadA) — which is
   why a BLAST of that interval against the element DB returns zero hits.
2. **The `canonical_v1` RB/LB entries in the element DBs are wrong**
   (`TGAGCGTCGCAAAGGCGCTCGGTCT` / `GGCCTCGGCCTGAGAGCCAAAACAC`): no
   `CAGGATATAT` core, present in no construct. Curating them out of the tracked
   `element_db/gmo_combined_db_v2.fa` etc. is still open.

### A2 — Filter D swallows BLAST failure and fails toward CANDIDATE
`scripts/s05/filter_d_altlocus.py:70-78`: `stderr=subprocess.DEVNULL` and a bare
`return False, 0.0, 0.0, 0.0` on `returncode != 0`. `False` means "not a
construct+host FP", so a crashed blastn silently promotes the site. Filters A and
B already log (`filter_a_host.py:74`, `filter_b_flank.py:73`); D was missed. The
declared `threads` parameter is also unused, and `blast_out` is dead.

### A3 — `canonical_triplet_min_identity` is a dead config knob
Loaded (`config_loader.py`), documented (`config.yaml:37`), asserted in
`tests/test_verdict_hardening.py:41` — and never applied. `report.py:216` passes
`matched_canonical=unique_elems`, i.e. *every* annotated element regardless of
BLAST identity. Rule 1 (canonical triplet) is the **highest-priority** rule and
overrides Filters B/C/D, so a low-identity spurious triplet forces CANDIDATE.

**Fix.** Gate `matched_canonical` on each element's best BLAST identity.

### A4 — left-clip consensus leaks internal `N` into assembly seeds
`_build_consensus(..., "left")` emits `N` at ambiguous positions and then only
`.strip("N")`s the ends, so internal `N`s survive into `seed_3p`. Downstream,
`StrandAwareSeedExtender.add_seqs` skips any sequence containing `N`
(`assembly.py:122`) and `_extend_right` matches overlaps exactly, so an internal
`N` blocks extension in exactly the direction the 3' seed exists to extend.

**Fix.** Return the maximal `N`-free suffix (the junction-adjacent end), which
makes the `left` branch symmetric with `right` (which already stops at the first
ambiguous position).

---

## B. Algorithmic improvements

### B1 — no read-quality filtering in junction discovery
`find_softclip_junctions` accepts soft clips from PCR duplicates, QC-fail reads,
and MAPQ-0 multimappers. Elsewhere the codebase requires MAPQ ≥ 20
(`assembly.py:466`), and MAPQ-driven false positives are the project's
top documented pitfall.

**Fix.** Skip `is_duplicate` / `is_qcfail`; add a `min_mapq` parameter
(default 0 → byte-identical behaviour on today's un-deduplicated BAMs) exposed
as `--min-mapq`.

### B2 — junction pairing is order-dependent greedy
Right clusters are walked in positional order and each grabs the nearest unused
left cluster, so the result depends on iteration order and can mis-pair
(`r=[100,110], l=[108]` pairs 100↔108 and drops 110↔108).

**Fix.** Enumerate all candidate pairs within the window, sort by
`(distance, r_pos, l_pos)`, and accept greedily — deterministic, and it never
picks a worse pair when a closer one exists.

### B3 — clustering constants are unreachable from the CLI
`cluster_window=50`, `MIN_CLUSTER_DEPTH=3` and the single-direction depth `5` are
hard-coded; `--junction-window` only feeds read extraction. This is the knob the
coverage-sensitivity work needs (`docs/measurements/coverage_sensitivity.md`
records detection failing below 15×).

**Fix.** Expose `--cluster-window`, `--min-cluster-depth`, `--min-single-depth`.

---

## Out of scope (recorded, not done)

- **PE-discordancy site discovery** (the single biggest detection gap; primitives
  exist in `s06b_junction_verify.py`). Needs real-BAM validation runs.
- **k=15 → 21 seed k-mer** for large plant genomes. Gated on a measurement sweep.
- **Unified `run_blastn` wrapper** over the 12 scattered `subprocess.run(["blastn"…])`
  call sites with inconsistent stderr handling.
- Every SLURM re-run in `work_implementation_plan.md` Phases 1–5.

---

## Task list

- [x] T1 `select_canonical_elements` pure helper + identity gate wired in `report.py` (A3)
- [x] T2 `_build_consensus` left-branch N-free suffix (A4)
- [x] T3 `_pair_clusters` deterministic global pairing (B2)
- [x] T4 read-quality filter + `min_mapq` (B1)
- [x] T5 configurable cluster depth/window, CLI wiring (B3)
- [x] T6 `_run_border_blast` + both repeats + HSP dedupe (A1)
- [x] T7 Filter D BLAST failure logging + dead-param cleanup (A2)
- [x] T8 full suite green, no verdict-snapshot drift

---

## Measured outcomes

`pytest tests/ -q` → **371 passed, 1 skipped** (from 322 + 1). Verdict snapshots:
26 fixtures, zero drift. DAG no-cycle: 15 modules. Line budgets: unchanged.

### A1 — border scan, on the real constructs

`_run_border_blast` against the first record of each construct FASTA:

| construct | old rows | new rows | new loci |
|---|---|---|---|
| `rice_G281` (8,958 bp) | 4 | 2 | `TDNA_border_B@279` (100%), `TDNA_border_A@6557` (100%) |
| `tomato_A2_3` (16,419 bp) | 4 | 2 | `TDNA_border_A@1` (100%), `TDNA_border_B@10142` (100%) |

The old scan reported each locus twice and matched the second repeat only as an
84%-identity degenerate hit of the first. `BORDER_MIN_ALIGN_LEN = 20` was needed
after the fix: querying two similar 25-mers at `-word_size 7` surfaces ~9 bp
scraps (rice_G281 positions 4,321 and 7,747) that are not borders.

### A4 / B2 — Phase 1 on the real rice host BAM

Old vs new clustering + consensus + pairing over `results/rice_G281/s04_host_map/`
(`Chr1`+`Chr2`, 78,697 right-clips / 77,069 left-clips):

| metric | old | new |
|---|---|---|
| sites | 5,864 | 5,706 |
| **`N` bases embedded in assembly seeds** | **2,327** | **0** |
| total seed bp | 327,325 | 317,777 |
| sites only in old | — | 169 |
| sites only in new | — | 11 |

Every dropped site lost a seed whose old length only cleared `min_clip=20`
because `N` characters padded it — `Chr1:457,108`'s 3' seed was 24 bp of which
**8 were `N`** (longest unambiguous run: 0 bp); `Chr1:1,726,127`'s was 23 bp of
which **14 were `N`**. Those were never usable seeds. A further 5 paired sites
moved position because the new pairing picks the closer partner
(`Chr1:2,014,692↔2,014,723`, 31 bp apart → `2,014,723↔2,014,743`, 20 bp apart).

The consensus rule is "longest unambiguous run, ties toward the junction", not
"longest `N`-free suffix". The suffix rule was the first attempt and code review
caught that it throws the whole clip away whenever the deepest,
junction-adjacent column happens to tie — a sensitivity regression the old
`.strip("N")` did not have. Correcting it recovered 48 sites (5,658 → 5,706)
while still embedding zero `N`.

On the `Chr3:16,380,000-16,500,000` region containing the rice ground truth, old
and new agree exactly: 25 sites, identical seeds, and
`Chr3:16,439,674-16,439,709` recovered by both.

### Full rice s05 re-run on this branch (SLURM 5799803, 1 h 00 m, COMPLETED)

Run into a **fresh** `results/rice_G281_algo_v2/` so the 2026-04-16 reports stay
intact. `--outdir-override` + `--host-bam-override`, `--no-remote-blast`.

**Superseded as the branch's result.** This run finished at 23:17 PDT, four
minutes before the per-insert border-count fix (`0914323`) and twelve before the
element-DB border correction (`496c328` / `6ba1f8b`). It is retained as the
evidence that produced those two fixes. The branch's measured result is the
2026-07-10 array (SLURM 5803578) over all 11 production samples, into
`results/hardened_20260710/`.

> **Corrected 2026-07-10.** An earlier revision of this section reported the old
> side of the table as 68 sites / 3 CANDIDATE / 32 FALSE_POSITIVE / 33 UNKNOWN,
> read off the 68 `insertion_*_report.txt` files in
> `results/rice_G281/s05_insert_assembly/`. **s05 never cleans its output
> directory**, so those 68 files are the union of three separate runs — 13 dated
> 2026-04-11, 34 dated 2026-04-14, 21 dated 2026-04-16 — of which 47 no longer
> correspond to any site the pipeline currently finds. `write_stats` opens
> `s05_stats.txt` with mode `"w"`, so that file alone describes one run, and its
> 21 site IDs are exactly the 21 reports dated 2026-04-16. Comparing report-file
> counts therefore measured accumulated debris, not a change in behaviour.

| | old (2026-04-16 run) | new (this branch) |
|---|---|---|
| sites reported | 21 | 21 |
| CANDIDATE | 3 | 3 |
| FALSE_POSITIVE | 9 | 10 |
| UNKNOWN | 9 | 8 |

Both runs discover the **same 21 sites, with identical site IDs**. Exactly one
verdict moves: `Chr11_28363144` UNKNOWN → FALSE_POSITIVE. Phase 1's site set for
rice is unchanged end to end, which is what the Chr3-region check predicted.

**The CANDIDATE set is identical**: `Chr3_16439674`, `Chr3_29074002`,
`Chr11_2877409`. The ground truth `Chr3:16,439,674` is CANDIDATE in both, and the
new run states the golden-case reason verbatim — *"1 element(s) annotated;
host_fraction=87% (all FP filters passed)"* — i.e. Filter A correctly did not
fire on the 1,024 bp foreign gap.

Phase 1 logged `0 clipped reads discarded`, confirming on a production BAM what
the unit tests assert: the new duplicate / QC-fail / MAPQ filters are a no-op at
default settings (these BAMs are not duplicate-marked).

**Three further defects the re-run exposed:**

1. **FIXED.** `generate_report` counted **every** line of `border_hits.tsv`, which
   is written once for all inserts. Every old report printed the same global
   count. `insertion_Chr3_16439674_report.txt` said *"T-DNA borders found: 10"*
   while all ten hits belonged to three other inserts.
2. RedGene has **never** detected a genuine T-DNA border in rice_G281. Under
   permissive settings (`-word_size 7 -evalue 1`), RB+LB against all 21 newly
   assembled inserts (57,285 bp) return **zero** HSPs ≥ 20 bp — the best is 13 bp
   at 92%. The old counts were sub-motif scraps. Consistent with
   `db/G281_construct.fa` being the empty pCAMBIA-1300 vector.
3. **OPEN — stale per-site reports accumulate.** s05 writes
   `insertion_<site>_report.txt` per site and never removes reports from a prior
   run into the same directory. `results/rice_G281/s05_insert_assembly/` holds 47
   reports for sites the pipeline no longer finds, indistinguishable by name from
   the 21 live ones. Only `s05_stats.txt` (truncated per run) is authoritative.
   This is the same misattribution class as defect 1 — a reader, or a report
   generator globbing `insertion_*_report.txt`, silently mixes runs. Rice is the
   only sample on disk affected: every other sample's report count equals its
   `insertion_sites` stat. Fix deferred so the 2026-07-10 re-run array observes
   one fixed tree; the array writes into fresh directories, where the bug cannot
   manifest.

### The full 11-sample re-run (2026-07-10/11) and what it surfaced

SLURM array `5803578` re-ran s05 for every production sample into
`results/hardened_20260710/`. 10 of 11 finished in-array; `soybean_UGT72E3`
timed out at 16 h assembling 106 sites sequentially. Result across the 9 samples
with a clean baseline (excluding the truncated `soybean_AtYUCCA6` April run):
sites −32, UNKNOWN −33, CANDIDATE +1, FALSE_POSITIVE +0 — the hardening strips
noise UNKNOWNs (N-free seeds) while preserving the CANDIDATE signal, no
ground-truth regression. Comparator: `compare_s05_runs.py` (reads `s05_stats.txt`,
never report globs). Six issues were filed from this run:

- **#18 (FIXED here)** — the per-site fan-out (`submit_s05_array.sh`) corrupted
  assemblies: `assemble_insert`'s scratch dirs (`_mm2_r{rnd}`, `_pilon_r{rnd}`,
  `_ssake_r{rnd}`, `_host_term.*`, `_foreign_refine`) were scoped to the shared
  `step_dir`, so concurrent sites raced on identical paths (`soybean_UGT72E3`
  array 5811054: 26/106 crashed on `_mm2ext.bam`). Fixed with a per-site
  `_scratch_{site_id}` root + atomic `_unmapped_R*.fq.gz` cache writes; guarded
  by `tests/test_s05_fanout_scratch_isolation.py`. UGT72E3 re-run as array 5811471.
- **#19** — Filter B rejects `corn_ND207`'s true insertion (5 bp off the
  documented site) because the corn border records carry the event's own maize
  flank; verdict priority puts Filter B above Filter A so the site is
  FALSE_POSITIVE by construction. Highest-value detection fix outstanding.
- **#20** — 3 of 4 bp-resolution `benchmark.md` ground truths not recovered;
  line 113 (tomato_Cas9_A2_3 "TP 1") contradicts the pipeline.
- **#21** — `tomato_WT` negative control now yields a transgene-positive site off
  a 26 bp `P-Act1-rice` tier hit; tomato/cucumber/soybean mask BEDs are empty.
- **#22** — s05 never cleans its output dir; stale `insertion_*_report.txt`
  accumulate (rice: 47 orphans across 3 runs) and already produced a wrong doc
  claim. Also `report.py:320` hides "T-DNA borders found: 0".
- **#23** — `db/G281_construct.fa` is the empty pCAMBIA-1300 vector while
  `manuscript.md` claims an exact construct sequence.

### Follow-up surfaced by this work

When one side of a paired cluster now yields a sub-`min_clip` consensus the whole
site is dropped, as it was before. Demoting it to a single-direction site using
the surviving seed would recover sensitivity, but needs a validation run.
