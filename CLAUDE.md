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

# Run individual step scripts directly
python scripts/s01_qc.py --r1 reads_R1.fq.gz --r2 reads_R2.fq.gz \
  --outdir results --sample-name my_sample --threads 16
```

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

`run_pipeline.py` orchestrates 7 analysis steps by calling standalone scripts in `scripts/` via subprocess. Config is loaded from `config.yaml` (YAML with per-sample settings). Each step script accepts `--outdir`, `--sample-name`, and step-specific arguments. Inter-step dependencies are wired in `build_step_cmd()` in run_pipeline.py.

### Pipeline flow and step dependencies
```
[1] fastp QC → [2] bwa → construct+UniVec BAM → [3] extract reads + mates
  → [3b] WT homology filter (optional)
  → [4] bwa → host BAM (bottleneck: ~5-7h)
  → [5] Targeted insert assembly + FP filtering  ★ CORE STEP
  → [6] CRISPR indel detection (optional, needs WT)
  → [7] copy number (depth ratio)
```

Steps 1-3 are fast (<30 min each). Step 4 is the bottleneck (~5-7h per sample). Step 5 is the core: finds insertion sites from host BAM soft-clips, assembles inserts, annotates elements, and applies 4 post-assembly false positive filters (host-fraction, construct-flanking, chimeric, alternative-locus). Steps 6 and 7 are optional downstream analyses.

### Key scripts
| Script | Purpose |
|--------|---------|
| `scripts/s03_extract_reads.py` | Extracts construct-hitting reads + mates |
| `scripts/s03b_homology_filter.py` | WT-based filtering of host-derived false positives |
| `scripts/s05_insert_assembly.py` | **CORE** — targeted assembly + 4 FP filters (host-fraction, construct-flanking, chimeric, alt-locus) |
| `scripts/s06_indel.py` | CRISPR editing detection (pileup-based, NOT bcftools) |
| `scripts/s07_copynumber.py` | Depth-ratio copy number estimation |
| `scripts/viz/plot_editing_profile.py` | CRISPResso2-style nucleotide quilt |
| `scripts/viz/plot_editing_effects.py` | Variant effect annotation (frameshift/synonymous/etc.) |
| `scripts/viz/plot_sample_summary.py` | 6-panel publication summary figure |

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

## Coding Conventions

- Python 3.11, type hints preferred
- Each step script is standalone with argparse CLI
- All file paths use `pathlib.Path`
- Logging to stderr (`print(..., file=sys.stderr)` or `logging` module)
- Output directories: `results/{sample}/s{NN}_{step_name}/`

## Known Pitfalls

- **MAPQ=60 false positives**: In plant genomes, unique mapping to a host-derived promoter region is still a false positive. Only WT-based filtering (s03b) reliably removes these.
- **Assembly stochasticity**: SPAdes contig extension direction affects junction detection. Contigs extending toward host-derived elements (TA29, pinII) produce overlapping alignments → false positive. Extending toward bacterial-origin elements (nptII, nos) → clean junction.
- **Coverage requirements**: ≥10x for cucumber (332 Mbp), ≥15x for rice (374 Mbp), ≥10x for tomato (833 Mbp). At 5x only partial detection (one junction side). At 3x complete failure.
- **Identity threshold for element_db**: Default `--min-identity 0.90` silently filters genuine junctions when using element_db (minimap2 alignment identity ~0.84). `run_pipeline.py` auto-detects element_db/combined_db and uses 0.70.
- **Maize-specific false positives**: When host IS maize, endogenous genes (Ubi1, zSSIIb, wx012) match construct elements. Border sequences contain ~100bp native flanking → 2.25M extracted reads vs 6K for rice.
- **BWA threading**: Earlier versions used `-t 2` due to futex deadlock on Pronghorn GPFS. Now uses `-t 16` successfully.

## SLURM Settings

- Partition: `cpu-s1-pgl-0`, Account: `cpu-s1-pgl-0`
- Resources: 16 CPUs, 64G RAM, 24h (for full pipeline with step 4 host mapping)
- See `run_clean.sh` for sbatch configuration

## Environment

```bash
micromamba create -n redgene -c conda-forge -c bioconda \
  python=3.11 bwa minimap2 samtools fastp spades bcftools \
  sra-tools pysam matplotlib biopython seqkit pigz
```
