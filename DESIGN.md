# GMO Positive Control Characterization Pipeline
# Project Design Document

## 1. Project Overview

**Goal**: Develop a cost-effective, sequencing-only pipeline to characterize
transgenic positive control materials for Korean quarantine agency GMO
detection assays.

**Context**:
- Korean quarantine agency currently operates multiplex amplicon PCR-based
  GMO detection (Option B)
- Won Yim lab (UNR) is creating transgenic plants as positive control
  reference materials for that amplicon assay
- Need to characterize each event: confirm construct integrity, find
  insertion site, estimate copy number, design event-specific primers
- Must be low-cost: Illumina PE150, 10-15x WGS, no Southern blot or ddPCR

**Host species**: Arabidopsis thaliana, Oryza sativa (rice), Glycine max (soybean)
**Transformation methods**: Agrobacterium-mediated + biolistic (both)
**Constructs**: Multiple promoter/terminator combinations to cover the
quarantine amplicon panel (CaMV 35S, FMV, Ubi1, Act1, NOS, OCS, nptII,
bar, hpt, cry1Ab, cp4-epsps, pat, etc.)

## 2. Pipeline Architecture

### Tool Stack (final, confirmed)
- **bwa mem**: all read-to-reference mapping (construct, host)
- **minimap2**: contig-to-reference mapping only (assembly contigs -> host)
- **samtools**: BAM processing, read extraction
- **fastp**: QC and adapter trimming
- **SPAdes**: local assembly (extracted reads only, very fast)
- **Kraken2**: taxonomic annotation (annotation, NOT filtering)
- **Primer3**: junction-spanning primer design (TgIF-style automation)
- **Python**: wrapper, junction parsing, copy number, visualization (matplotlib)

### Pipeline Steps

```
Step 1:   QC + trim (fastp)
Step 2:   Map reads to construct reference (bwa mem)
Step 2v:  Visualization: Construct coverage profile plot
Step 3:   Extract construct-hitting reads + all mate reads (samtools)
Step 4:   Local assembly of extracted reads (SPAdes --careful, multi-k)
Step 5:   Map assembled contigs to host reference (minimap2 -x asm5)
Step 6:   Chimeric contig detection + junction coordinate extraction
Step 6v:  Visualization: Junction structure diagram
Step 7:   Map all reads to host reference (bwa mem) [PARALLEL with Steps 2-6]
Step 8:   Collect unmapped reads, Kraken2 annotation [PARALLEL]
Step 9:   De novo assembly of unmapped reads, BLAST vs backbone [PARALLEL]
Step 10:  Copy number estimation (read depth ratio + junction count)
Step 10v: Visualization: Copy number depth comparison plot
Step 11:  Event-specific primer design (Primer3, TgIF-style)
Step 12:  Amplicon panel validation (expected vs observed concordance)
Step 13:  Sensitivity test (0.1%-5% GMO DNA mixing series)
```

### Dependency Graph

```
Step 1 (QC)
  |
  +---> Step 2 (map to construct) ---> Step 3 (extract) ---> Step 4 (assembly)
  |         |                              ---> Step 5 (contigs to host)
  |         |                              ---> Step 6 (junction extraction)
  |         |                              ---> Step 11 (primer design)
  |         |
  |         +--- (construct depth) -------> Step 10 (copy number, needs Step 7 too)
  |
  +---> Step 7 (map to host) --+---> Step 8 (Kraken2) ---> Step 9 (de novo, backbone)
                                |
                                +---> Step 10 (copy number, needs Step 2 too)
```

### Execution Environment
- Pronghorn HPC (SLURM)
- Conda/Mamba for software management
- Python wrapper submitting SLURM jobs with dependency chaining

## 3. Critical Design Decisions (from discussion)

### 3.1 Why assembly, not split-read or depth-gap
- TgIF (JHU-APL) uses depth-gap method because it targets ONT long reads
  where single reads span entire transgene + flanking host
- Illumina 150bp reads cannot span junctions reliably enough for depth-gap
- Assembly combines multiple overlapping reads to build junction-spanning
  contigs, which is the correct approach for short reads
- Assembly input is only ~2,000-6,000 read pairs (construct hits + mates),
  so SPAdes finishes in under 1 minute

### 3.2 Read extraction strategy
- Map reads to construct with bwa mem
- Extract mapped read NAMES (samtools view -F 4 | cut -f1)
- Include supplementary alignments (samtools view -f 2048) for soft-clipped reads
- Name-sort the BAM, then use samtools view -N + samtools fastq to extract
  both reads of each pair
- This avoids FASTQ name format issues (/1 /2 vs space-delimited)
- BWA mem outputs unmapped mates in the BAM, so mate reads are automatically included

Key commands for Step 3:
```bash
# Get names of reads hitting construct
samtools view -F 4 construct_mapped.bam | cut -f1 > names.txt
samtools view -f 2048 construct_mapped.bam | cut -f1 >> names.txt
sort -u names.txt > hit_names.txt

# Name-sort BAM and extract pairs
samtools sort -n construct_mapped.bam -o nsort.bam
samtools view -N hit_names.txt nsort.bam \
  | samtools fastq -1 R1.fq -2 R2.fq -s singles.fq -n
```

### 3.3 Why bwa (not minimap2) for read mapping
- bwa mem is the de facto standard for Illumina short read mapping
- No reviewer justification needed
- minimap2 -x sr works but edge sensitivity on short references (~10kb
  construct) may be slightly lower
- minimap2 is used for contig-to-host mapping (Step 5) where it excels
  (-x asm5 for assembly-to-reference)

### 3.4 Kraken2 as annotation, not filter
- Unmapped reads from host mapping contain: microbial contamination,
  organellar divergence, host PAV regions, AND transgene sequences
- Transgene CDS of bacterial origin (nptII, cry1Ab, cp4-epsps) would be
  classified as bacterial by Kraken2 and incorrectly filtered out
- Solution: use Kraken2 to LABEL reads taxonomically, but do NOT remove
  them before de novo assembly
- Post-assembly, use multi-layer annotation (Kraken + BLAST + ORF) to
  distinguish genuine transgene contigs from contamination

### 3.5 Copy number from sequencing only
- Construct median depth / genome-wide median depth = copy number estimate
- ~0.5 = hemizygous single-copy, ~1.0 = homozygous or 2-copy hemizygous
- Cross-validate with junction contig count:
  - 2 junction contigs (LB + RB) = single-copy
  - Additional internal junctions = tandem multi-copy
  - Different host chromosomes = multi-locus
- Limitation: 10-15x coverage has high variance; single vs 2-copy may overlap
- Defense: position this as "sequencing-only pipeline" methodological novelty
- Validate accuracy using known single-copy Arabidopsis T-DNA lines (ABRC)

### 3.6 Host-only read pairs are not a problem
- If both reads of a pair map entirely to host (no transgene sequence),
  they are invisible to Step 2 construct mapping
- This is fine: these reads contain no transgene information anyway
- Junction detection relies on reads where at least ONE mate contains
  transgene sequence, and those reads ARE captured
- Assembly then extends into host using mate information
- Host extension is limited to ~150bp (one mate length), which is usually
  sufficient for unique locus identification
- If not sufficient: iterative assembly (use round-1 contigs as bait for
  round-2 read extraction) or targeted PCR + Sanger

### 3.7 Confidence scoring (inspired by TgIF)
- **High**: chimeric contig with bp-level junction + both LB and RB junctions found
- **Medium**: one junction contig found (one border only)
- **Low**: discordant pairs suggest insertion region but no precise junction sequence

## 4. Visualization Requirements (inspired by TgIF)

All plots generated with Python matplotlib.

### Plot 1: Construct coverage profile (Step 2v)
- X-axis: position along construct (bp)
- Y-axis: read depth
- Color-coded regions for each element (promoter, CDS, terminator)
- Annotated element boundaries
- Purpose: verify construct integrity, detect truncation/rearrangement

### Plot 2: Junction structure diagram (Step 6v)
- Linear diagram showing chimeric contig structure
- Left: host sequence with chr:position
- Center: junction breakpoint
- Right: transgene element with name
- Below: number of supporting reads, confidence level
- Both LB and RB junctions shown if found

### Plot 3: Copy number depth comparison (Step 10v)
- X-axis: genomic position or element name
- Y-axis: normalized read depth
- Horizontal reference line at genome-wide median (1.0x)
- Construct region depth shown relative to reference
- Visual distinction of 0.5x (hemizygous) vs 1.0x (homozygous)

## 5. GMO Element Database

A separately maintained curated FASTA database of 43 transgenic elements:
- 12 promoters, 7 terminators, 19 CDS, 5 regulatory elements
- 29 elements downloadable from NCBI GenBank (accession + coordinates defined)
- 14 elements manually curated (synthetic/construct-specific variants)
- Standardized header format: class|name|organism|accession|boundary|version
- Build script: build_gmo_element_db.py (Biopython + NCBI Entrez)

This database is used for:
- Step 2: construct reference (primary use case, construct FASTA from own lab)
- Future extension: multi-element screening for unknown samples

The database itself is a novel contribution (no existing open-access GMO
element FASTA database exists) and should be deposited on GitHub.

## 6. File Structure

```
gmo-positive-control/
├── CLAUDE.md                    # Claude Code instructions
├── DESIGN.md                    # This document
├── environment.yml              # Conda environment
├── config.yaml                  # Sample sheet + paths
├── run_pipeline.py              # Main Python wrapper (SLURM job submission)
├── scripts/
│   ├── s01_qc.py               # fastp wrapper
│   ├── s02_map_construct.py    # bwa mem to construct
│   ├── s03_extract_reads.py    # samtools read extraction
│   ├── s04_local_assembly.py   # SPAdes
│   ├── s05_map_contigs.py      # minimap2 contigs to host
│   ├── s06_junction.py         # chimeric contig parser
│   ├── s07_map_host.py         # bwa mem to host
│   ├── s08_kraken.py           # Kraken2 annotation
│   ├── s09_denovo.py           # de novo assembly + BLAST backbone
│   ├── s10_copynumber.py       # read depth copy number
│   ├── s11_primer.py           # Primer3 (TgIF-style)
│   ├── viz/
│   │   ├── plot_construct_coverage.py   # Plot 1
│   │   ├── plot_junction_structure.py   # Plot 2
│   │   └── plot_copynumber_depth.py     # Plot 3
│   └── utils/
│       ├── slurm.py            # SLURM job submission helper
│       └── parse.py            # BAM/BLAST parsers
├── db/
│   ├── constructs/             # User construct FASTA refs
│   ├── host/                   # Host genome refs (TAIR10, IRGSP, Wm82)
│   ├── kraken2/                # Kraken2 custom DB
│   └── gmo_elements/           # GMO element DB (build_gmo_element_db.py output)
└── results/                    # Per-sample output
    └── {sample_name}/
        ├── qc/
        ├── construct_mapping/
        ├── extracted_reads/
        ├── assembly/
        ├── junction/
        ├── host_mapping/
        ├── unmapped/
        ├── copynumber/
        ├── primers/
        └── plots/
```

## 7. Config File Format

```yaml
samples:
  - name: At_35S_nptII_event1
    r1: /path/to/reads_R1.fastq.gz
    r2: /path/to/reads_R2.fastq.gz
    host: arabidopsis
    construct: db/constructs/35S_nptII.fa
    method: agrobacterium

hosts:
  arabidopsis:
    fasta: db/host/TAIR10.fa
    size_mb: 135
  rice:
    fasta: db/host/IRGSP-1.0.fa
    size_mb: 400
  soybean:
    fasta: db/host/Wm82.a4.v1.fa
    size_mb: 1100

kraken2_db: db/kraken2/custom_db

slurm:
  partition: cpu-s2-core-0
  account: cpu-s5-bch709-0
  mail: wyim@unr.edu
```

## 8. Conda Environment

```yaml
name: gmo-ctrl
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - fastp=0.23.4
  - bwa=0.7.18
  - minimap2=2.28
  - samtools=1.19
  - spades=3.15.5
  - blast=2.15.0
  - kraken2=2.1.3
  - primer3-py=2.0.3
  - bedtools=2.31.1
  - pandas
  - numpy
  - matplotlib
  - pyyaml
  - biopython
  - pysam
```

## 9. Key References

- Debode et al. 2019, Sci Rep 9:15595 (NGS + enrichment for GMO detection,
  39-element database design)
- TgIF: https://github.com/jhuapl-bio/TgIF (transgene insertion finder,
  ONT-based, inspiration for visualization and primer3 integration)
- Genomics & Informatics 2024 (Korean context NGS-based GMO characterization review)
- EUginius database (element nomenclature)
- JRC GMO-Amplicons (EU Open Data Portal, amplicon FASTA download)

## 10. Known Limitations

1. Illumina 150bp cannot resolve complex insertion structures (tandem,
   inverted repeat, multi-locus) as well as long-read sequencing
2. Copy number estimation at 10-15x has high variance; single vs 2-copy
   distinction may be unreliable
3. Host-derived elements in constructs (Ubi1, Act1) cause multi-mapping
   ambiguity that can affect coverage estimation
4. Phase 3 (unknown element discovery via unmapped reads) is exploratory,
   not guaranteed -- sensitivity for truly novel elements is limited
5. Manual sequence curation needed for some elements in GMO element DB
   (synthetic/codon-optimized variants differ from wild-type GenBank entries)
