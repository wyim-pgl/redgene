# RedGene — GMO/LMO Transgene Characterization Pipeline

Assembly-based Illumina WGS pipeline for characterizing transgenic plant
insertions. Designed for Korean quarantine GMO/LMO detection assays and
CRISPR/Cas9 editing verification.

**Repository**: [wyim-pgl/redgene](https://github.com/wyim-pgl/redgene)

## What It Does

Given paired-end Illumina reads from a transgenic plant:

1. **Finds transgene insertion sites** at base-pair resolution (assembly-based junctions)
2. **Estimates copy number** from read depth ratios
3. **Determines zygosity** (homozygous vs heterozygous insertion)
4. **Detects CRISPR editing** (treatment-specific indels vs WT)
5. **Filters false positives** with four post-assembly filters and a rule-based verdict
6. **Reports** each site with an interpretable verdict (CANDIDATE / FALSE_POSITIVE / UNKNOWN)

## Background

Plant T-DNA constructs often contain host-derived promoters (e.g., rice Ubi1,
maize Act1, tobacco TA29) that create false positive chimeric reads when mapped
against construct databases. This pipeline distinguishes true transgene
junctions from homologous-sequence artifacts using WT-based homology filtering
and four assembly-level false-positive filters (see below).

The approach is based on the modified TranSeq method (Bae et al. 2022,
Communications Biology) adapted for generic GMO element databases rather than
specific vector sequences.

## Pipeline Architecture

`run_pipeline.py` orchestrates the steps below by invoking standalone scripts in
`scripts/` via subprocess. `STEP_ORDER` is the single source of truth for
sequencing: `["1", "2", "3", "4", "4b", "5", "6", "7", "8"]`.

```
                    Illumina PE150 Reads
                            |
                     [1] fastp QC + trim
                            |
             [2] bwa mem → construct + UniVec BAM
                            |
          [3] Extract construct-hitting reads + mates
                            |
       [3b] WT-based homology filter (optional, recommended)
                            |
        [4] bwa mem → host BAM (full genome)   ← bottleneck (~5-7h)
                            |
   [4b] SPAdes → per-sample construct contigs (element DB for s05)
                            |
        [5] Targeted insert assembly + 4 FP filters + verdict   ← core
              (uses 4b contigs + common_payload DB)
             |                    |                     |
   [6] CRISPR indel      [7] Copy number       [8] PDF insertion report
   detection (vs WT)      (depth ratio)         (opt-in; reads s05/s07)
```

Steps 1-3 are fast (<30 min each). Step 4 (host mapping) is the bottleneck.
Steps 4b, 6, 7, 8 are optional. The core run is `--steps 1-5`; a full run is
`--steps 1-7`; the PDF report (step 8) is opt-in and never part of `1-5`.

### Step Details

| Step | Script | Description |
|------|--------|-------------|
| 1 | `s01_qc.py` | QC + adapter trimming (fastp) |
| 2 | `s02_construct_map.py` | Map reads to construct + UniVec (bwa mem) |
| 3 | `s03_extract_reads.py` | Extract construct-hitting reads + mates |
| 3b | `s03b_homology_filter.py` | WT-based filter of host-derived false positives (optional) |
| 4 | `s04_host_map.py` | Map all reads to host genome (bwa mem) — **bottleneck** |
| 4b | `s04b_construct_assembly.py` | De novo SPAdes assembly of construct-hitting reads → per-sample element DB |
| 5 | `s05_insert_assembly.py` | **Core**: targeted insert assembly + 4 FP filters + verdict |
| 6 | `s06_indel.py` | CRISPR indel detection (pileup-based, treatment vs WT) |
| 7 | `s07_copynumber.py` | Copy number estimation (depth ratio) |
| 8 | `reports/insertion_pdf.py` | PDF insertion report (opt-in; reads existing s05/s07 output) |

**Step 4b — sample-specific construct assembly**: Many transgenic samples carry
payload sequences (e.g., AtYUCCA6, bar marker, custom RNAi hairpins) not present
in the shared `element_db/gmo_combined_db.fa`. Step 4b runs SPAdes on the
construct-hitting reads from s03 to build a per-sample FASTA, which s05 loads via
`--extra-element-db` so the payload sequences are annotated instead of reported
as UNKNOWN. `--steps 1-5` includes 4b by default (it slots between 4 and 5);
bypass with `--steps 1-4,5`.

### Step 5 — the core: four false-positive filters + verdict

Step 5 finds insertion sites from host-BAM soft-clips, assembles the insert, and
applies four post-assembly false-positive filters, then computes a rule-based
verdict:

| Filter | Module | Catches |
|--------|--------|---------|
| A — host-fraction | `s05/filter_a_host.py` | inserts that are mostly host genomic DNA |
| B — construct-flanking | `s05/filter_b_flank.py` | junctions inside construct-flanking host homology |
| C — chimeric | `s05/filter_c_chimeric.py` | multi-locus chimeric mis-assemblies |
| D — alternative-locus | `s05/filter_d_altlocus.py` | host DNA mis-assembled via construct-element homology (e.g. CaMV 35S copies) |

**Verdict priority** (see `docs/measurements/verdict_priority.md`, implemented in
`scripts/s05/verdict.py`):

```
canonical_triplet > host_endogenous > Filter B > Filter C > Filter D > Filter A
  > elements_present → CANDIDATE > UNKNOWN_FP > UNKNOWN
```

Step 5 was refactored (Issue #4) from a 4003-line monolith into 11
single-responsibility modules under `scripts/s05/`; the monolith
`s05_insert_assembly.py` is now a thin CLI entrypoint. See
`docs/architecture/s05-modules.md`.

### Planned — Step 5b: construct assembly (spec written, not yet implemented)

An opt-in sub-step that reconstructs and characterizes the *complete inserted
construct*: a sample-level classified/annotated contig inventory plus per-site
links to s05 CANDIDATE sites. It reuses s04b's pre-filter SPAdes contigs
(assemble-then-classify), preserving unknown payload and junction-spanning
(chimeric) contigs that s04b's marker filter discards. Design:
`docs/superpowers/specs/2026-06-30-step5b-construct-assembly-design.md`;
implementation plan: `docs/superpowers/plans/2026-06-30-step5b-construct-assembly.md`.

## Quick Start

```bash
# 1. Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

# 2. Configure the sample in config.yaml (host_reference, construct_reference,
#    reads.r1/r2, optional wt_control, optional grna). See existing entries.

# 3. Core pipeline (QC → construct map → extract → host map → insert assembly)
python run_pipeline.py --sample rice_G281 --steps 1-5 --threads 16

# Full pipeline (+ CRISPR indel + copy number)
python run_pipeline.py --sample rice_G281 --steps 1-7 --threads 16

# Re-run insert assembly only, skipping slow remote NCBI BLAST
python run_pipeline.py --sample rice_G281 --steps 5 --threads 16 --no-remote-blast

# PDF insertion report (opt-in; reads existing results/ without re-running)
python run_pipeline.py --sample rice_G281 --steps 8

# Preview commands without running
python run_pipeline.py --sample rice_G281 --steps 1-5 --dry-run
```

Individual step scripts are standalone with their own argparse CLI, e.g.:

```bash
python scripts/s01_qc.py --r1 reads_R1.fq.gz --r2 reads_R2.fq.gz \
  --outdir results --sample-name my_sample --threads 16
```

### SLURM

```bash
sbatch run_clean.sh   # rice + all tomato samples end-to-end incl. visualization
```

Partition/account `cpu-s1-pgl-0`; 16 CPUs, 64G RAM, 24h for a full run with step 4.

## Element Database

`element_db/gmo_combined_db.fa` contains 131 GMO-related sequences scraped from
the EUginius database (https://euginius.eu/), including:

- **Promoters**: CaMV 35S, nos, Ubi1 (maize), Act1 (rice), TA29 (tobacco), etc.
- **Terminators**: nos, OCS, CaMV 35S, pinII, etc.
- **Selection markers**: nptII, bar, hpt, pat, etc.
- **Construct-specific**: pCAMBIA vectors, event-specific amplicons
- **Full-length**: reference constructs from NCBI

When s05 uses an element DB, the alignment-identity floor is auto-relaxed from
0.90 to 0.70 (element_db/combined_db matches are ~0.84 identity by minimap2);
`run_pipeline.py` auto-detects this. See "Known Pitfalls" in `CLAUDE.md`.

### Common payload DB

`element_db/common_payload.fa` provides 9 canonical bacterial/viral transgene
markers (bar, nptII, hpt, gusA, gfp, egfp, P-CaMV35S, P-nos, T-ocs), wired as an
always-on `--common-payload-db` argument in s05 so the most frequently used
selection and reporter genes are annotated even when the per-sample or
event-specific construct DB does not list them. Built by
`element_db/build_common_payload.sh` via NCBI efetch.

## False Positive Filtering

Plant genomes contain sequences homologous to common construct elements, creating
false positive junction calls that require filtering:

| Source | Example | MAPQ | Solution |
|--------|---------|------|----------|
| Host-derived promoters | Rice Ubi1 in pCAMBIA | Low (0-10) | element/WT filter |
| Unique homologous regions | Chr2:8.4M in rice G281 | High (60) | WT filter (s03b) |
| Repetitive elements | Multiple chromosomes | Low (0-5) | MAPQ / Filter C |
| Element-homology mis-assembly | CaMV 35S copies in host | High | Filter D (alt-locus) |

**Key insight — MAPQ alone is insufficient.** In rice G281 testing,
Chr2:8,432,860 had MAPQ=60 (unique mapping) but was a false positive from pCAMBIA
vector homology. WT-based filtering (s03b) is the only reliable method for these
cases; the four assembly-level filters (A-D) then catch the rest at the insert
level.

## CRISPR Editing Detection (Step 6)

Two modes for detecting CRISPR editing-induced indels:

### Mode 1: gRNA-guided (recommended)
1. **BLAST gRNA** to host genome to find on-target and off-target sites
2. **Direct pileup parsing** at predicted cut sites (±50bp window)
3. **Treatment vs WT subtraction**: keep only treatment-specific indels
4. **Allele frequency** from read counts

### Mode 2: de novo (no gRNA)
1. **Genome-wide variant calling** in treatment and WT
2. **Subtraction**: keep only treatment-specific indels (≥ 2bp)
3. **PAM motif check**: NGG within 3-8bp of the indel
4. **NHEJ signature**: deletions 1-20bp, insertions 1-5bp

### Key design choice
gRNA-guided mode parses `samtools mpileup -Q 0` directly (no base-quality
filter) to avoid losing CRISPR indels whose anchor base has borderline quality.
`bcftools call` decomposes complex indels at low depth (e.g. splitting a 9bp
deletion into 1bp events) — direct pileup parsing preserves the original indel
representation from reads.

## Environment

```bash
micromamba create -n redgene -c conda-forge -c bioconda \
  python=3.11 bwa minimap2 samtools fastp spades bcftools \
  sra-tools pysam matplotlib biopython seqkit pigz
```

> Env-path gotcha: on Pronghorn, `micromamba activate redgene` resolves to
> `/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/`. For non-interactive
> subprocesses, prepend the env bin directly:
> `export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"`.

### Key tools
| Tool | Version | Purpose |
|------|---------|---------|
| bwa | 0.7.18 | Read mapping to construct/host |
| minimap2 | 2.28 | Contig mapping + homology detection |
| samtools | 1.21 | BAM processing + pileup |
| bcftools | 1.23 | Variant calling (de novo indel mode) |
| fastp | 0.23.4 | QC + adapter trimming |
| SPAdes | 4.0.0 | Assembly of extracted reads |
| BLAST+ | 2.x | Element annotation (local + remote nt) |
| sra-tools | 3.1.1 | SRA download (prefetch + fasterq-dump) |
| pysam | 0.22.1 | BAM parsing (Python) |
| matplotlib | 3.9 | Visualization + PDF report |
| biopython | 1.84 | Sequence handling |

## Testing

```bash
pytest tests/ -q                              # full suite
pytest tests/test_verdict_snapshots.py -v     # verdict regression (snapshot fixtures)
pytest tests/test_s05_import_dag.py -v         # scripts/s05 DAG no-cycle guard
```

Snapshot fixtures in `tests/fixtures/verdict_snapshots/` (rice_G281,
cucumber_line225) are the primary regression signal for `compute_verdict` — never
edit them without a deliberate behavior change.

## Test Results

### Rice G281 (2-copy head-to-head T-DNA, ~29x)
- **True insertion**: Chr3:16,439,674 (lactoferrin RNAi, 36bp deletion) → CANDIDATE
- **False positives**: pCAMBIA vector-host homology → FALSE_POSITIVE (Filters A-D)
- **Coverage minimum**: ≥15x for reliable detection

### Tomato Micro-Tom Cas9 (CRISPR editing, PRJNA692070)
- **T-DNA insertion**: SLM_r2.0ch08:65,107,378 detected (A2_3, ~10x)
- **Cas9 construct**: CaMV 35S + nos + nptII + OCS identified in assembled contigs
- **CRISPR editing**: 9bp deletion at SlPHD_MS1 (chr04); 6bp deletion at SlAMS (chr08)
- **Coverage minimum**: ≥10x for junction detection with element_db

### Coverage sensitivity
| Coverage | Junction detection | Notes |
|----------|-------------------|-------|
| ≥ 15x | Reliable, high confidence | Recommended (rice) |
| 10-15x | Usually works | ≥10x for cucumber/tomato |
| 5-10x | Construct presence confirmed | One junction side only |
| < 5x | Unreliable | Insufficient for assembly |

## Config Format

`config.yaml` defines samples with `host_reference`, `construct_reference`,
`reads.r1/r2`, optional `wt_control` (sample key for the WT BAM used by s03b/s06),
and optional `grna` (path to a gRNA file for step 6). See existing entries.

## Output Structure

```
results/<sample>/
  s01_qc/                # Trimmed reads + fastp HTML/JSON
  s02_construct_map/     # Construct + UniVec BAM + stats
  s03_extract/           # Extracted read pairs (+ filtered if s03b)
  s04_host_map/          # Host BAM + flagstat + depth
  s04b_construct_asm/    # Per-sample SPAdes construct contigs
  s05_insert_assembly/   # Per-site insert FASTA + report.txt (verdicts) + annotation TSVs
  s06_indel/             # editing_sites.tsv (CRISPR samples)
  s07_copynumber/        # Copy number estimates
```

## References

- Bae S, Park YC, et al. (2022). Characterization of transgene insertions by
  resequencing. *Communications Biology* 5:671.
- Bae S, et al. (2022). CRISPR/Cas9-mediated editing of SlFAD2.
  *Horticultural Science and Technology* 40:81-92.
- EUginius GMO detection methods: https://euginius.eu/

## License

MIT
