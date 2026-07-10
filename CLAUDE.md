# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RedGene is an assembly-based Illumina WGS pipeline for characterizing transgenic plant insertions. It finds transgene insertion sites (bp resolution), estimates copy number, determines zygosity, and detects CRISPR/Cas9 editing events. Designed for Korean quarantine GMO/LMO detection assays.

**Repository**: [wyim-pgl/redgene](https://github.com/wyim-pgl/redgene)

## Running the Pipeline

```bash
# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

# Run core pipeline (QC → construct map → extract → host map → insert assembly)
python run_pipeline.py --sample rice_G281 --steps 1-5 --threads 16

# Run full pipeline (+ CRISPR indel + copy number)
python run_pipeline.py --sample rice_G281 --steps 1-7 --threads 16

# Re-run insert assembly only (skip slow NCBI remote BLAST)
python run_pipeline.py --sample rice_G281 --steps 5 --threads 16 --no-remote-blast

# Dry run (preview commands)
python run_pipeline.py --sample rice_G281 --steps 1-5 --dry-run

# PDF insertion report (Issue #6, opt-in; reads existing results/ without re-running)
python run_pipeline.py --sample rice_G281 --steps 8

# Run individual step scripts directly
python scripts/s01_qc.py --r1 reads_R1.fq.gz --r2 reads_R2.fq.gz \
  --outdir results --sample-name my_sample --threads 16
```

### Testing

```bash
# Full pytest suite (baseline: 387 PASS + 1 skipped, ~20s)
pytest tests/ -q

# Single file / single test
pytest tests/test_verdict_snapshots.py -v
pytest tests/test_compute_verdict.py::test_rule_1_canonical_triplet_promotes -v

# DAG no-cycle guard for scripts/s05/ (Issue #4 invariant — 15 modules registered)
pytest tests/test_s05_import_dag.py -v

# Line-budget guards added in Session 6 (monolith < 300, fanout::main < 250)
pytest tests/test_submit_s05_array.py -v

# Phase 1 site discovery: consensus, clustering, cluster pairing, read filters
pytest tests/test_site_discovery.py -v

# Phase 1 end-to-end on a synthetic BAM (needs minimap2) — the ONLY guard on
# which sites find_softclip_junctions emits; verdict snapshots sit downstream
pytest tests/test_find_softclip_junctions.py -v

# T-DNA border scan (includes a real-blastn check against db/G281_construct.fa)
pytest tests/test_border_scan.py -v
```

Snapshot fixtures live in `tests/fixtures/verdict_snapshots/` (rice_G281, cucumber_line225). Never edit these without a deliberate behavior change — snapshot drift is the primary regression signal for `compute_verdict`. Note the snapshots cover `compute_verdict` only — they cannot detect a change in which sites Phase 1 discovers; `tests/test_find_softclip_junctions.py` is that guard.

### SLURM batch run
```bash
sbatch run_clean.sh   # Runs rice + all tomato samples end-to-end including visualization
```

### Visualization scripts (run after pipeline steps complete)
```bash
# CRISPResso2-style editing profile
python scripts/viz/plot_editing_profile.py \
  --treatment-bam results/{sample}/s04_host_map/{sample}_host.bam \
  --wt-bam results/tomato_WT/s04_host_map/tomato_WT_host.bam \
  --host-ref db/SLM_r2.0.pmol.fasta \
  --grna-targets results/{sample}/s06_indel/grna_targets.tsv \
  --editing-sites results/{sample}/s06_indel/editing_sites.tsv \
  --sample-name {sample} --outdir results

# Variant effect annotation (easyGWAS-style)
python scripts/viz/plot_editing_effects.py \
  --editing-sites results/{sample}/s06_indel/editing_sites.tsv \
  --gff db/SLM_r2.0.gff3.gz --host-ref db/SLM_r2.0.pmol.fasta \
  --sample-name {sample} --outdir results
```

## Architecture

`run_pipeline.py` orchestrates 7 core + 2 optional (4b, 8) analysis steps by calling standalone scripts in `scripts/` via subprocess. Config is loaded from `config.yaml` (YAML with per-sample settings). Each step script accepts `--outdir`, `--sample-name`, and step-specific arguments. Inter-step dependencies are wired in `build_step_cmd()` in run_pipeline.py.

`STEP_ORDER` in run_pipeline.py is the single source of truth for step sequencing: `["1", "2", "3", "4", "4b", "5", "5b", "6", "7", "8"]`. Steps 5b (construct assembly) and 8 (PDF report) sit after the `5` boundary and are opt-in — never part of default `--steps 1-5`.

### Pipeline flow and step dependencies
```
[1] fastp QC → [2] bwa → construct+UniVec BAM → [3] extract reads + mates
  → [3b] WT homology filter (optional)
  → [4] bwa → host BAM (bottleneck: ~5-7h)
  → [4b] SPAdes → sample-specific construct contigs
  → [5] Targeted insert assembly + FP filtering (uses 4b contigs + common_payload)
  → [5b] Construct assembly (opt-in: global construct inventory + per-site links)
  → [6] CRISPR indel detection (optional, needs WT)
  → [7] copy number (depth ratio)
  → [8] PDF insertion report (opt-in; reads existing s05/s07 output)
```

Steps 1-3 are fast (<30 min each). Step 4 is the bottleneck (~5-7h per sample). Step 4b is optional but recommended: it de novo assembles the construct-hitting reads from s03 with SPAdes to produce a per-sample FASTA that s05 consumes via `--extra-element-db`, ensuring sample-specific payload sequences (e.g., AtYUCCA6, bar marker) get annotated rather than reported as UNKNOWN. Step 5 is the core: finds insertion sites from host BAM soft-clips, assembles inserts, annotates elements, and applies 4 post-assembly false positive filters (host-fraction, construct-flanking, chimeric, alternative-locus). Steps 6 and 7 are optional downstream analyses.

### Key scripts
| Script | Purpose |
|--------|---------|
| `scripts/s03_extract_reads.py` | Extracts construct-hitting reads + mates |
| `scripts/s03b_homology_filter.py` | WT-based filtering of host-derived false positives |
| `scripts/s04b_construct_assembly.py` | De novo SPAdes assembly of construct-hitting reads for per-sample element DB |
| `scripts/s05_insert_assembly.py` | **Thin entrypoint** (102 lines, post-Issue #4). Re-exports backward-compat symbols from `scripts/s05/` package; orchestration lives in `scripts/s05/fanout_orchestrator.py::main`. Step 5 is the core: targeted assembly + 4 FP filters (host-fraction, construct-flanking, chimeric, alt-locus). |
| `scripts/s05b_construct_assembly.py` | Step 5b (opt-in): reconstruct + characterize the inserted construct — global classified/annotated contig inventory (assemble-then-classify, reuses s04b `contigs_all.fasta`) + per-site links to s05 CANDIDATE sites. Self-contained; local-first BLAST with opt-in remote nt. |
| `scripts/s06_indel.py` | CRISPR editing detection (pileup-based, NOT bcftools) |
| `scripts/s07_copynumber.py` | Depth-ratio copy number estimation |
| `scripts/reports/insertion_pdf.py` | Step 8 PDF report (matplotlib/PdfPages; no reportlab dep) |
| `scripts/viz/plot_editing_profile.py` | CRISPResso2-style nucleotide quilt |
| `scripts/viz/plot_editing_effects.py` | Variant effect annotation (frameshift/synonymous/etc.) |
| `scripts/viz/plot_sample_summary.py` | 6-panel publication summary figure |
| `tools/verify_coc.py` | Chain-of-custody log integrity checker (monotone ts, pre/post pairing, hash chain) |

### scripts/s05/ package (Issue #4 — 7-session refactor, COMPLETED 2026-04-29)

The 4003-line `s05_insert_assembly.py` monolith was split into 11 single-responsibility modules. The monolith remains the CLI entrypoint (run_pipeline.py and SLURM arrays still invoke `python scripts/s05_insert_assembly.py`), but its `main()` is a thin shim that delegates to `scripts/s05/fanout_orchestrator.py::main`. All function bodies were extracted byte-identical (verbatim policy); the ONE non-verbatim change was Session 6's split of `main()` into 4 phase helpers.

**Authoritative reference:** `docs/architecture/s05-modules.md`. Original design spec: `docs/superpowers/specs/2026-04-19-s05-module-split-design.md`.

DAG stage ordering enforced by `tests/test_s05_import_dag.py`; later-stage modules may import from earlier stages, never the reverse.

| Module | Stage | Role |
|--------|-------|------|
| `primitives.py` | 0 | `log`, `revcomp`, FASTA/FASTQ I/O, dataclasses (`JunctionCluster`, `InsertionSite`, `LegacyJunction`, `TierResult`), `_parse_src_tag` |
| `site_discovery.py` | 1 | Phase 1: `find_softclip_junctions` + `_read_donates_clip` + `_cluster_positions` + `_pair_clusters` + `_build_consensus` + `_apply_mask_bed` + `MASKED_SOURCE_TAG` + `parse_legacy_junctions` + `_extract_seeds_at_positions` |
| `classify.py` | 1 | Phase 1.5: `classify_site_tiers` + `_filter_host_endogenous` + `_should_replace` + `_SRC_TIER`/`_TIER2_SRCS`/`_UNKNOWN_SRC_WARNED` state |
| `read_extraction.py` | 2 | Phase 2: `extract_candidate_reads` + `extract_unmapped_paired` |
| `assembly.py` | 3 | Phase 3: `StrandAwareSeedExtender` + `recruit_by_kmer` + `pilon_fill` + `assemble_insert` (1178 lines, the core pipeline function) |
| `annotation.py` | 4 | Phase 4: `_run_local_blast` + `_run_remote_blast` + `_merge_annotations` + `annotate_insert` + `_run_border_blast` / `TDNA_BORDER_REPEATS` |
| `filter_a_host.py` | 5 | `_blast_insert_vs_host` (host-fraction + foreign-gap measurement) |
| `filter_b_flank.py` | 5 | Construct-flanking FP check (BLAST construct vs host, slop-overlap test) |
| `filter_c_chimeric.py` | 5 | Multi-locus chimeric FP check (strict-identity off-target chrom aggregation) |
| `filter_d_altlocus.py` | 5 | Alt-locus FP check (construct + host combined coverage) |
| `verdict.py` | 6 | Pure `compute_verdict` + `select_canonical_elements` + `FilterEvidence` + `VerdictRules` + `_apply_canonical_override` |
| `config_loader.py` | 7 | `load_verdict_rules` (YAML → `VerdictRules`) |
| `report.py` | 8 | Phase 5: `generate_report` + `write_stats` |
| `fanout_orchestrator.py` | 9 | `main()` + `_run_phase_1_1_5` + `_run_phase_2_3` + `_run_phase_4` (CLI dispatch) |

When adding a new `scripts/s05/*.py`, add its stage to the `STAGE` dict in `tests/test_s05_import_dag.py`. The 2 v1.0 shims (`annotate_report.py`, `per_site.py`) were retired in Session 7.

### Verdict priority order

Canonical priority (documented in `docs/measurements/verdict_priority.md`, implemented in `scripts/s05/verdict.py:compute_verdict`):

```
canonical_triplet > host_endogenous > Filter B > Filter C > Filter D > Filter A
  > elements_present → CAND > UNK_FP > UNK
```

Changes to this ordering require a snapshot regeneration pass.

## Critical Design Decisions

1. **Assembly-based junction detection**, not depth-gap (depth-gap requires long reads)
2. **Direct pileup parsing for CRISPR** (`samtools mpileup -Q 0`), not `bcftools call`. bcftools decomposes complex indels at low depth (e.g., splits 9bp deletion into 1bp events)
3. **-Q 0 base quality** in gRNA-guided mode: CRISPR indel anchor bases sometimes have Q18, below the Q20 default threshold
4. **WT-based homology filtering** is essential: plant T-DNA constructs contain host-derived promoters (Ubi1, Act1, TA29) that create MAPQ=60 false positives
5. **Element database** (131 EUginius elements) instead of requiring exact vector sequence — needs lower identity threshold (0.70 vs 0.90)
6. **Alternative-locus filter** (Filter D): minimap2-based check that catches host genomic DNA mis-assembled via construct-element homology (e.g., CaMV 35S promoter copies in host genome)

## Config Format

`config.yaml` defines samples with: `host_reference`, `construct_reference`, `reads.r1/r2`, optional `wt_control` (sample key for WT BAM), optional `grna` (path to gRNA file for step 6). See existing entries for examples.

## Reference Data (in db/, gitignored)

| File | Source |
|------|--------|
| `Osativa_323_v7.0.fa` | Rice genome (Phytozome, 374 Mbp) |
| `Osativa_323_v7.0.gene_exons.gff3` | MSU/RGAP v7.0 annotation (55,986 genes) |
| `SLM_r2.0.pmol.fasta` | Tomato Micro-Tom genome (Kazusa, 833 Mbp) |
| `SLM_r2.0.gff3.gz` | NCBI Gnomon annotation (chr names remapped to SLM_r2.0ch\*) |
| `CucSat_B10v3.fa` | Cucumber B10v3 genome (GCA_001483825.3, 332 Mbp, 8035 contigs) |
| `Zm_B73_v5.fa` | Corn B73 RefGen_v5 (GCF_902167145.1, 2.18 Gbp, 10 chr) |
| `Gmax_v4.0.fa` | Soybean Wm82.a4.v1 (GCF_000004515.6, 1.1 Gbp) |
| `element_db/gmo_combined_db.fa` | 131 GMO elements from EUginius (incl. thaumatin II) |
| `gmo_corn_combined_db.fa` | 192 seqs: 130 EUginius + 62 corn LB/RB borders (Sci Rep 2025) |
| `element_db/common_payload.fa` | 9 canonical bacterial/viral transgene markers (bar, nptII, hpt, gusA, gfp, egfp, P-CaMV35S, P-nos, T-ocs). Built by `element_db/build_common_payload.sh` via NCBI efetch. Wired as always-on `--common-payload-db` for s05. |

## Coding Conventions

- Python 3.11, type hints preferred
- Each step script is standalone with argparse CLI
- All file paths use `pathlib.Path`
- Logging to stderr (`print(..., file=sys.stderr)` or `logging` module)
- Output directories: `results/{sample}/s{NN}_{step_name}/`

## Known Pitfalls

- **MAPQ=60 false positives**: In plant genomes, unique mapping to a host-derived promoter region is still a false positive. Only WT-based filtering (s03b) reliably removes these.
- **Assembly stochasticity**: SPAdes contig extension direction affects junction detection. Contigs extending toward host-derived elements (TA29, pinII) produce overlapping alignments → false positive. Extending toward bacterial-origin elements (nptII, nos) → clean junction.
- **Coverage requirements**: **≥15x for reliable detection across all three hosts** (rice 374 Mbp, cucumber 332 Mbp, tomato 833 Mbp). The `docs/measurements/coverage_sensitivity.md` matrix (SLURM 5630185, 12/12 tasks) tightened the earlier "cucumber/tomato ≥10x" guidance: rice 5x = MISS, tomato 10x stochastically empty, cucumber ≤10x = MISS with first HIT at 15x. At 10-15x detection is partial (often one junction side); at ≤5x it fails. Note the mid-coverage failure mode: at rice 15x the GT junction is detected but the partial assembly over-extends into off-target host and is mis-called FALSE_POSITIVE by the chimeric filter.
- **Identity threshold for element_db**: Default `--min-identity 0.90` silently filters genuine junctions when using element_db (minimap2 alignment identity ~0.84). `run_pipeline.py` auto-detects element_db/combined_db and uses 0.70.
- **Maize-specific false positives**: When host IS maize, endogenous genes (Ubi1, zSSIIb, wx012) match construct elements. Border sequences contain ~100bp native flanking → 2.25M extracted reads vs 6K for rice.
- **BWA threading**: Earlier versions used `-t 2` due to futex deadlock on Pronghorn GPFS. Now uses `-t 16` successfully.
- **T-DNA border repeats — RB/LB assignment is NOT established in this repo.** Both constructs carry two distinct 25-bp repeats: `TGGCAGGATATATTGTGGTGTAAAC` (`TDNA_border_A`) and `TGACAGGATATATTGGCGGGTAAAC` (`TDNA_border_B`) — rice_G281 at 6,557/279, tomato_A2_3 at 1/10,142. The `canonical_v1` `RB-TDNA`/`LB-TDNA` entries in `db/*.fa` appear in *no* construct and contradict both. Until an external reference settles it, `border_hits.tsv` labels are sequence-derived, and code must not infer strand or T-DNA orientation from them.
- **`min_clip` counts consensus bases, not `N`.** `_build_consensus` returns the longest unambiguous run (ties resolve toward the junction). A seed padded with `N` used to clear `min_clip` and then silently fail k-mer extension (`add_seqs` skips `N`-containing reads; `_extend_right` requires an exact overlap).
- **Junction discovery now rejects `is_duplicate` / `is_qcfail` reads even at `--min-mapq 0`.** This is a no-op on today's BAMs — `bwa mem | samtools sort` never marks duplicates and sets no QC-fail flag. If a `samtools markdup` step is ever added upstream, those reads will start being dropped from clip clusters. That is intended (a duplicate inflates cluster depth without adding independent evidence), but it is a silent change in cluster depth, so re-baseline `--min-cluster-depth` if you add markdup.
- **Canonical triplet promotion is identity-gated.** Rule 1 outranks Filters B/C/D, so an element only joins `matched_canonical` when its best BLAST hit clears `canonical_triplet_min_identity` (default 0.90). Annotation itself runs at a 0.70 floor.

## SLURM Settings

- Partition: `cpu-s1-pgl-0`, Account: `cpu-s1-pgl-0`
- Resources: 16 CPUs, 64G RAM, 24h (for full pipeline with step 4 host mapping)
- Allocate 4 GB RAM per thread.
- See `run_clean.sh` for sbatch configuration
- **`unset SBATCH_PARTITION SBATCH_ACCOUNT SBATCH_TIMELIMIT` at the top of every batch script.** The login shell exports them, and SLURM precedence is CLI > env > script — so `SBATCH_TIMELIMIT=13-23:00:00` silently overrides your `#SBATCH --time=16:00:00` (observed on job 5799803). Verify with `sacct -j <id> --format=Timelimit,ReqCPUS,ReqMem`. Full write-up: lab wiki `guide/pronghorn.md`.

## Environment

```bash
micromamba create -n redgene -c conda-forge -c bioconda \
  python=3.11 bwa minimap2 samtools fastp spades bcftools \
  sra-tools pysam matplotlib biopython seqkit pigz
```

### Env path gotcha

Despite the create command above, `micromamba activate redgene` on this machine resolves to `/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/`, not `~/micromamba/envs/redgene/` (the user's `.bashrc` chains miniconda3 + micromamba such that shared pgl envs take precedence). For non-interactive shells / subprocesses where `micromamba` is not on PATH, prepend the env bin directly:

```bash
export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"
```

`gh` CLI is NOT in the redgene env; it lives at `~/micromamba/bin/gh`.
