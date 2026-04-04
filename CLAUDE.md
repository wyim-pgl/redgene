# CLAUDE.md

## Project: GMO Positive Control Characterization Pipeline

Read DESIGN.md first. It contains the full architecture, all design decisions,
and rationale from extensive discussion.

## What this project does

Illumina PE150 sequencing-based pipeline to characterize transgenic plant
positive control materials for Korean quarantine GMO detection assays.

Input: Illumina PE150 reads (10-15x WGS) from transgenic plants
Output: Insertion site coordinates, copy number, event-specific primers,
        validation reports, and publication-ready plots

## Tech stack

- Python 3.11 wrapper with SLURM job submission
- bwa mem (read mapping), minimap2 (contig mapping), samtools (BAM processing)
- fastp (QC), SPAdes (assembly), Kraken2 (annotation), Primer3 (primers)
- matplotlib (visualization), pysam (BAM parsing), biopython (sequence handling)
- Conda/Mamba for environment management
- Runs on Pronghorn HPC (SLURM scheduler)

## Architecture

13-step pipeline. See DESIGN.md Section 2 for the full step list and
dependency graph. Key design principles:

1. Assembly-based junction detection (NOT depth-gap; depth-gap is for long reads)
2. Extract construct-hitting reads + mates -> local assembly -> BLAST contigs to host
3. Kraken2 is used as ANNOTATION (label reads), never as FILTER (would remove bacterial-origin transgene CDS)
4. Copy number from read depth ratio only (no Southern blot, no ddPCR)
5. Three matplotlib visualizations: construct coverage, junction structure, copy number

## File organization

```
gmo-positive-control/
├── CLAUDE.md              <- you are here
├── DESIGN.md              <- full design document (READ THIS)
├── environment.yml        <- conda environment
├── config.yaml            <- sample sheet
├── run_pipeline.py        <- main entry point
├── scripts/               <- per-step scripts (s01 through s11)
│   ├── viz/               <- matplotlib visualization scripts
│   └── utils/             <- SLURM helper, BAM/BLAST parsers
├── db/                    <- reference databases
└── results/               <- per-sample output
```

## Coding conventions

- Python 3.11, type hints preferred
- Each step is a standalone script callable with --sample and --config args
- SLURM jobs submitted via scripts/utils/slurm.py helper
- All file paths use pathlib.Path
- Logging to stderr, results to stdout or files
- Config in YAML format
- Output directories created automatically per sample

## SLURM settings

- Partition: cpu-s2-core-0
- Account: cpu-s5-bch709-0
- Default: 8 CPUs, 32G RAM, 4 hours
- Phase 3 (unmapped analysis): 8 CPUs, 64G RAM, 8 hours

## Key implementation details

### Step 3 (read extraction) -- most critical step

```bash
samtools view -F 4 construct.bam | cut -f1 > names.txt
samtools view -f 2048 construct.bam | cut -f1 >> names.txt
sort -u names.txt > hit_names.txt
samtools sort -n construct.bam -o nsort.bam
samtools view -N hit_names.txt nsort.bam \
  | samtools fastq -1 R1.fq -2 R2.fq -s singles.fq -n
```

This handles all read name formats and guarantees mate extraction.

### Step 6 (junction detection)

Parse chimeric contigs: contigs that have partial alignment to host AND
partial alignment to construct. The boundary between the two alignments
is the junction coordinate.

### Step 10 (copy number)

```python
copy_number = construct_median_depth / genome_median_depth
# ~0.5 = hemizygous single-copy
# ~1.0 = homozygous single-copy
```

Cross-validate with junction contig count (2 = single-copy, >2 = multi-copy).

### Confidence scoring

- High: both LB and RB junction contigs found, bp-level resolution
- Medium: one junction contig only
- Low: discordant pair evidence but no junction contig

## Visualization (matplotlib)

Three plots per sample:
1. **Construct coverage profile**: depth along construct, color-coded by element
2. **Junction structure diagram**: linear diagram of chimeric contig with host/transgene annotation
3. **Copy number depth plot**: normalized depth comparison, reference line at 1.0x

## Things to watch out for

- manual_sequences.fa contains PLACEHOLDER sequences for some Bt genes
  (cry2Ab2, cry34Ab1, etc.) -- these must be replaced with actual construct sequences
- element_master_list.tsv coordinates for CaMV V00141.1 need verification
- SPAdes on ~3000 read pairs finishes in <1 min, do not over-allocate resources
- Kraken2 MUST NOT be used as a filter (see DESIGN.md Section 3.4)
- Host-derived elements (Ubi1, Act1) cause multi-mapping; use unique bacterial
  markers (nptII, bar) for reliable copy number estimation

## Development priority

Round 1 (core): Steps 1-6 + Plot 1 + Plot 2
  - Get one Arabidopsis event through the pipeline end-to-end first

Round 2 (backbone): Steps 7-9
  - Run on Round 1 data, assess necessity

Round 3 (validation): Steps 10-13 + Plot 3
  - Wet lab coordination needed for Steps 12-13

## Related resources

- GMO element DB builder: see element_db/ directory
  (gmo_db.py + element.list + manual_sequences.fa)
- TgIF reference: https://github.com/jhuapl-bio/TgIF
  (ONT-based tool; we borrowed visualization and primer3 integration concepts)
