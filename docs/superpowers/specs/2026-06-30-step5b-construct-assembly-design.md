# Step 5b — Construct Assembly: Design Spec

**Date:** 2026-06-30
**Status:** Approved design → ready for implementation plan
**Goal:** Add a new pipeline sub-step `5b` that reconstructs and characterizes the
*complete inserted transgene construct* from Illumina reads, at two levels:
**global** (sample-level construct/foreign contig inventory) and **per-site**
(which contig corresponds to each s05 CANDIDATE insertion site).

This fulfills the long-standing "step 9 / construct build" goal (orphan
`run_s09_*.sh` / `run_foreign_assembly.sh` exploratory scripts from 2026-04),
now integrated as a first-class, tested pipeline step.

---

## 1. Motivation

s05 finds insertion **sites** (bp resolution) and assembles a **per-site insert**
(`insertion_{site}_insert.fasta`), but those fragments are junction-anchored. They
do not give a consolidated answer to "what construct was inserted into this
sample, and which elements does it contain?" — especially when the payload
(e.g. AtYUCCA6, UGT72E3) is not in the shared `element_db`.

`s04b_construct_assembly.py` already de-novo assembles the construct-hitting reads
(from s03) with SPAdes, but it **filters contigs to known marker DBs** and discards
everything else — so novel/unknown payload contigs are thrown away.

Step 5b reuses s04b's *pre-filter* assembly (`contigs_all.fasta`) and instead of
discarding non-marker contigs, **classifies and annotates all of them**, preserving
unknown payload and junction-spanning (chimeric) contigs, then **links them back to
s05 sites**.

### Relationship to existing steps (no overlap, clean division of labor)

| Step | Produces | Keeps unknown payload? | Per-site? |
|------|----------|------------------------|-----------|
| s04b | marker-filtered contigs → element DB for s05 | No (filtered out) | No |
| s05  | insertion sites + per-site junction-anchored insert | n/a | Yes (sites) |
| **5b** | **classified + annotated construct inventory, linked to sites** | **Yes** | **Yes (links)** |

---

## 2. Pipeline position

`STEP_ORDER = ["1", "2", "3", "4", "4b", "5", "5b", "6", "7", "8"]`

- Runs **after s05** (depends on s05 site reports for linkage).
- **Opt-in.** Because `"5b"` sits *after* the `"5"` boundary in `STEP_ORDER`, the
  default core range `--steps 1-5` does **not** include it (contrast: `4b` sits
  between `4` and `5`, so `1-5` *does* include `4b`). Run via `--steps 5b` or
  `--steps 1-5b`. This is intentional: 5b is a characterization extra, like
  steps 6/7/8.
- Follows the existing `4b` "sub-step" naming convention; avoids the name
  collision with the unrelated dev tool `scripts/s09_ground_truth_baseline.py`
  (which is moved to `scripts/util/` as part of this work — see §9).

---

## 3. Inputs

| Input | Source | Required? |
|-------|--------|-----------|
| `results/{sample}/s04b_construct_asm/contigs_all.fasta` | s04b pre-filter SPAdes contigs (reused) | Preferred |
| s03 construct reads `{sample}_construct_R1/R2.fq.gz` | for SPAdes fallback when `contigs_all.fasta` is absent/empty | Fallback |
| host reference FASTA | config `host_reference` | For host classification |
| `element_db/gmo_combined_db.fa` + `element_db/common_payload.fa` | local element annotation | For element classification |
| `results/{sample}/s05_insert_assembly/insertion_*_report.txt` | CANDIDATE site coords + verdict (parsed) | For per-site linkage |
| `--remote-blast` flag (opt-in) | NCBI nt on UNKNOWN contigs only | Optional |

**s05 site source of truth:** there is no consolidated `insertion_report.tsv`.
Sites are derived by scanning `insertion_*_report.txt` filenames for the
`{chrom}_{pos}` site_id and parsing each file's `Verdict:` line. Step 5b links
to **CANDIDATE** sites by default (configurable to include all verdicts).

---

## 4. Algorithm (assemble-then-classify)

### 4.1 Obtain contigs
1. If `s04b/contigs_all.fasta` exists and is non-empty → use it.
2. Else if `s04b/contigs.fasta` exists (filtered, when `--no-filter` was used the
   raw is here) → use it.
3. Else run SPAdes on s03 construct reads using the same invocation as s04b
   (`spades.py --only-assembler --careful`). On SPAdes failure / too few reads →
   empty output, exit 0 (see §6).
4. Length filter: discard contigs `< MIN_CONTIG_LEN` (200 bp).

### 4.2 Classify each contig
For each surviving contig, compute two coverage fractions:

- **`host_frac`** = fraction of contig length covered by `blastn` hits vs the host
  reference with `pident ≥ HOST_PIDENT` (90.0). Overlapping hit intervals are
  merged before summing covered bp.
- **`element_frac`** = fraction of contig length covered by `blastn` hits vs the
  element DBs with `pident ≥ ELEMENT_PIDENT` (80.0; element DB matches are
  lower-identity by design — see CLAUDE.md). Hits from **both** element DBs
  (`gmo_combined_db.fa` + `common_payload.fa`) are pooled, then intervals merged.

**Classification (first match wins — priority order):**

```
1. host_frac ≥ CHIMERIC_MIN  AND  element_frac ≥ CHIMERIC_MIN   → "chimeric"
2. element_frac ≥ ELEMENT_MIN                                    → "construct"
3. host_frac   ≥ HOST_MIN                                        → "host"        (discarded from FASTA)
4. otherwise                                                     → "foreign_unknown"
```

**Concrete thresholds (hard-coded module constants):**

```python
MIN_CONTIG_LEN  = 200     # bp
HOST_PIDENT     = 90.0    # % identity for a host blastn hit to count
ELEMENT_PIDENT  = 80.0    # % identity for an element blastn hit to count
CHIMERIC_MIN    = 0.20    # both host_frac AND element_frac ≥ this → junction-spanning
ELEMENT_MIN     = 0.50    # element-dominant → construct
HOST_MIN        = 0.70    # host-dominant → host (discarded)
LINK_SLOP       = 100     # bp tolerance for contig↔site coordinate match
```

Rationale: `chimeric` is checked first because a junction-spanning contig has
both host and element signal; classifying it as plain `host` or `construct` would
lose the junction. `construct` before `host` so an element-dominant contig with
minor host noise is not mislabeled host.

### 4.3 Annotate
- **Always (local):** record per-contig top element hit(s) from the element-DB
  blastn already run in 4.2 (sseqid, pident, aligned length, e-value).
- **`--remote-blast` (opt-in):** run `blastn -remote -db nt` **only on
  `foreign_unknown` contigs** (the candidate novel payload) to identify them.
  Bounded `max_target_seqs 5`, `evalue 1e-10`. Failure/timeout → leave
  `remote_hit` empty, warn, continue (never fails the step).

### 4.4 Per-site linkage
1. `minimap2 -x asm10` align **all non-host contigs** to the host reference.
2. For each s05 CANDIDATE site `(chrom, pos)`: a contig links to the site if it has
   a host alignment on `chrom` whose start or end lies within `LINK_SLOP` (100 bp)
   of `pos`.
3. `junction_side` recorded as `upstream` / `downstream` from alignment orientation
   relative to `pos`.
4. `link_confidence`:
   - `high` — linked contig is `chimeric` (physically spans the junction)
   - `medium` — linked contig is `construct`/`foreign_unknown` abutting within slop
   - (a site with no qualifying contig produces no link row)

---

## 5. Outputs (`results/{sample}/s05b_construct_assembly/`)

| File | Contents |
|------|----------|
| `construct_contigs.fasta` | all non-host contigs (`construct`, `chimeric`, `foreign_unknown`). **Chimeric contigs are written whole** (not trimmed) to preserve junction context. FASTA headers annotated: `>{contig_id} class={class} top_element={elem} host_frac={..} element_frac={..}` |
| `contig_classification.tsv` | `contig_id, length, class, host_frac, element_frac, top_element, element_pident, element_alnlen, remote_hit`. Includes **all** classified contigs (host rows too, for transparency); only non-host contigs are written to the FASTA. |
| `site_construct_links.tsv` | `site_id, verdict, contig_id, junction_side, link_confidence` |
| `construct_summary.txt` | sample-level: total construct bp; element inventory (distinct elements detected with best identity); contig counts by class; number of CANDIDATE sites linked; presence of novel/`foreign_unknown` payload |
| `_classification_blast/` | raw blastn TSVs (host, element) + minimap2 PAF, for debugging |

All TSVs have a header row. Empty runs still write headers (see §6).

---

## 6. Error / edge handling (mirrors s04b "always exit 0 + always produce expected files" contract)

| Condition | Behavior |
|-----------|----------|
| `contigs_all.fasta` absent **and** s03 reads empty/absent | write empty `construct_contigs.fasta` + header-only TSVs, **exit 0** |
| SPAdes fallback fails (too few reads) | empty outputs, **exit 0**, log to stderr |
| 0 contigs after length filter | header-only TSVs, empty FASTA, exit 0 |
| host reference missing | skip host classification (host_frac=0 for all), stderr warning, continue |
| element DB missing | skip element classification (element_frac=0), stderr warning, continue |
| s05 output dir missing / no CANDIDATE sites | skip linkage (still produce global outputs), `site_construct_links.tsv` header-only, warn |
| `--remote-blast` failure / timeout | proceed with local-only annotation, warn, **do not fail the step** |

The step **never** raises an unhandled exception for missing/empty inputs;
downstream consumers (step 8 PDF) can always find the expected output paths.

---

## 7. `run_pipeline.py` wiring

- Add `"5b"` to `STEP_ORDER` after `"5"`.
- Add `STEP_MAP["5b"] = "scripts/s05b_construct_assembly.py"`.
- Add `s05b = outdir / sname / "s05b_construct_assembly"` path var.
- Add `elif step == "5b":` branch in `build_step_cmd()` passing:
  `--contigs-all {s04b}/contigs_all.fasta`, `--s03-r1/--s03-r2` (SPAdes fallback),
  `--s05-dir {s05}`, `--host-ref {host_reference}`,
  `--element-db element_db/gmo_combined_db.fa --element-db element_db/common_payload.fa`,
  `--outdir --sample-name --threads`, and `--remote-blast` **only when the global
  `--no-remote-blast` flag is NOT set** (matches s05 remote semantics).
- Update the module docstring step list and the `--steps` help text.

---

## 8. Step 8 PDF integration (included in this plan)

Extend `scripts/reports/insertion_pdf.py` with a **"Construct Inventory"** section
that, **when `s05b_construct_assembly/` output is present**, renders:
- element inventory table (from `construct_summary.txt` / `contig_classification.tsv`)
- per-CANDIDATE-site linked contig + class (from `site_construct_links.tsv`)
- a flag when `foreign_unknown` payload was detected

When 5b output is absent, the section is omitted (the PDF remains valid without
5b having run). A dispatcher/guard test asserts both presence and absence paths.

---

## 9. Name-collision cleanup (in scope)

`scripts/s09_ground_truth_baseline.py` is a **benchmark/dev tool** (extracts
ground-truth soft-clip metrics; NOT a pipeline step, NOT in `STEP_MAP`), paired
with `scripts/util/analyze_coverage_sensitivity.py`. To free the `s0Nb`/step
namespace and reduce confusion:

- Move `scripts/s09_ground_truth_baseline.py` → `scripts/util/ground_truth_baseline.py`.
- Update references in `tests/` and `docs/measurements/coverage_sensitivity.md`
  that invoke it (historical plan docs under `docs/superpowers/plans/` are
  immutable records — left as-is).

---

## 10. Module structure

Single standalone script `scripts/s05b_construct_assembly.py`:
- argparse CLI, `pathlib.Path`, stderr logging, exit-0 contract — matches s04b.
- Reuse the pure helper `read_fasta` (and `revcomp` if needed) from
  `scripts/s05/primitives.py`; keep BLAST/minimap2 runners self-contained
  (small, like s04b's `_filter_contigs_by_markers`).
- Pure, testable helpers: `_merge_intervals`, `_coverage_fraction`,
  `_classify_contig`, `_parse_s05_sites`, `_link_contigs_to_sites`.

---

## 11. Testing

New test file `tests/test_s05b_construct_assembly.py`:

- **Classification unit tests** — synthetic contigs vs a tiny host ref + tiny
  element DB fixture exercising all four classes (`host`, `construct`, `chimeric`,
  `foreign_unknown`) at the exact threshold boundaries.
- **Interval/coverage helpers** — `_merge_intervals`, `_coverage_fraction` with
  overlapping/adjacent/disjoint cases.
- **s05 site parsing** — `_parse_s05_sites` from synthetic `insertion_*_report.txt`
  (CANDIDATE filtering correctness).
- **Linkage** — synthetic contigs + minimap2 PAF (or mocked) → expected links,
  slop boundary, `link_confidence` levels.
- **Edge/exit-0 contract** — empty contigs, missing host ref, missing element DB,
  missing s05 dir, remote-blast failure → all produce expected files + exit 0.
- **`--help` smoke test.**
- **run_pipeline dispatch test** — `build_step_cmd("5b", ...)` emits the expected
  argv (extend the existing dispatch test pattern).
- **PDF section test** — Construct Inventory renders with 5b output present;
  omitted when absent.

No verdict snapshot changes (5b does not touch s05) — the existing 26 snapshot
fixtures must remain byte-identical.

---

## 12. Out of scope (deferred to the separate "code-advancement" track)

Per the agreed plan, broad code-quality work runs as a **separate spec** after
5b lands: deferred s05 cleanups (`HOST_ENDO_*`/`CLIP_HOST_*` constant duplication,
`__all__` definitions, `_run_border_blast` extraction, RB==LB consensus
verification), static-analysis sweep, and coverage gaps. Not part of this spec.
