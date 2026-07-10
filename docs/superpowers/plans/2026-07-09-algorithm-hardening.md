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

### Follow-up surfaced by this work

When one side of a paired cluster now yields a sub-`min_clip` consensus the whole
site is dropped, as it was before. Demoting it to a single-direction site using
the surviving seed would recover sensitivity, but needs a validation run.
