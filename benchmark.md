# RedGene Benchmark — T-DNA Insertion Site Detection Tools

## Overview

Comparative benchmark of four T-DNA insertion site detection tools using validated transgenic plant WGS datasets across five crop species.

### Tools Compared

| Tool | Year | Paper | Method | Input | Language |
|------|------|-------|--------|-------|----------|
| **RedGene** | 2025 | — | Assembly-based (SPAdes → minimap2 chimeric contig) | WGS PE reads | Python 3.11 |
| **TDNAscan** | 2019 | Front Genet 10:685 | Read mapping + split-read (BWA) | WGS PE reads | Python 3 |
| **T-DNAreader** | 2025 | Genome Biol 26:36 | Read mapping (Bowtie2) + junction extraction | WGS PE reads | Python 3 |
| **TC-hunter** | 2022 | BMC Genomics 23:717 | Read mapping + discordant pairs (BWA) + Nextflow | WGS PE reads | Python 2.7/Nextflow |

### Construct Availability

A critical differentiator is construct/vector sequence availability. TDNAscan, T-DNAreader, and TC-hunter all require the **exact construct FASTA**. RedGene can use an **element database** (131 EUginius elements) when the full vector sequence is unavailable.

| Sample | Construct | Public FASTA Available | Source |
|--------|-----------|----------------------|--------|
| rice_G281 | pCAMBIA-1300 derivative (cry1Ac + bar) | Yes (pCAMBIA-1300 backbone, 8,958 bp) | GenBank |
| tomato_Cas9_A2_3 | pBAtC-SlAMS/SlPHD_MS1 (Cas9 + gRNAs) | Yes (Cas9 construct, lab-provided) | PRJNA692070 |
| cucumber (224/212/225) | pRUR528s (thaumatin II) | **No** — not deposited in NCBI/GenBank | Park et al. 2020 Sci Rep |
| corn_ND207 | ND207 (Bt cry1Ab) | **No** — proprietary (Syngenta) | Ji et al. 2024 |
| soybean_AtYUCCA6 | pB7WG2-AtYUCCA6 (bar + AtYUCCA6) | **No** — not deposited | Kim et al. 2021 Plant Biotechnol Rep |
| soybean_UGT72E3 | pB7WG2-UGT72E3 (bar + UGT72E3) | **No** — not deposited | Kim et al. 2021 Plant Biotechnol Rep |

**Implication**: For 4 of 6 test datasets, only RedGene (via element_db) can attempt detection. Other tools are **not applicable** when construct FASTA is unavailable.

### Test Datasets

| Sample | Species | Genome Size | Coverage | Expected Insertion | Source |
|--------|---------|-------------|----------|-------------------|--------|
| rice_G281 | *Oryza sativa* | 374 Mbp | ~30x | Chr3:16,439,674 (2 copies, head-to-head) | Phytozome |
| tomato_Cas9_A2_3 | *S. lycopersicum* Micro-Tom | 833 Mbp | ~25x | chr08:65,107,378 (T-DNA + CRISPR edits) | PRJNA692070 |
| cucumber_line224 | *Cucumis sativus* | 332 Mbp | ~37x | Chr2, promoter of G6936, 1 T-DNA copy | PRJNA638559 |
| cucumber_line212 | *Cucumis sativus* | 332 Mbp | ~37x | Chr6, intergenic, 1 T-DNA copy | PRJNA638559 |
| cucumber_line225 | *Cucumis sativus* | 332 Mbp | ~37x | Chr2, intergenic, 2 T-DNA copies + backbone | PRJNA638559 |
| corn_ND207 | *Zea mays* | 2.18 Gbp | ~5x | Chr3:181,367,276 (Bt construct) | CRA026358 |
| soybean_AtYUCCA6 | *Glycine max* | 1.1 Gbp | ~28x | Site I: Chr18 Glyma.18g235800 (5+ tandem copies), Site II: Chr19 Glyma.19g245800 | PRJNA627303 |
| soybean_UGT72E3 | *Glycine max* | 1.1 Gbp | ~29x | Chr18:51,882,860 (Glyma.18g226800 CDS) | PRJNA627303 |

---

## Environment Setup

### Env 1: TDNAscan
```bash
micromamba create -n env_tdnascan -c conda-forge -c bioconda python=3.10 bwa samtools -y
micromamba activate env_tdnascan
cd ~/tools/TDNAscan
```
- **Repo**: https://github.com/BCH-RC/TDNAscan
- **Clone**: `~/tools/TDNAscan/`
- **Status**: Ready. Example test: 20 insertions detected in mt4_chr1 (19.5s)
- **Output**: `{prefix}/5.{prefix}_insertion.bed`
- **Note**: Input reads must be uncompressed FASTQ (not .gz)

### Env 2: T-DNAreader
```bash
micromamba create -n env_tdnareader -c conda-forge -c bioconda \
  python=3.10 bowtie2 star samtools pysam pandas numpy -y
micromamba activate env_tdnareader
cd ~/tools/TDNAreader
```
- **Repo**: https://github.com/CDL-HongwooLee/TDNAreader
- **Clone**: `~/tools/TDNAreader/`
- **Status**: Ready. Example test: 20 insertion sites (40 TIS entries with L/R junctions)
- **Output**: `{prefix}.TDNA.txt`
- **Note**: bowtie2 index prefix must match FASTA path; `--gzip` for .gz input; `-g` is bowtie2 index prefix. Requires pre-built bowtie2 indices for **both** construct AND host genome.

### Env 3: TC-hunter
```bash
micromamba create -n env_tchunter -c conda-forge -c bioconda -c r \
  python=2.7 bwa samtools nextflow pandas \
  r-base r-circlize r-dplyr r-data.table -y
micromamba activate env_tchunter
cd ~/tools/TC_hunter
```
- **Repo**: https://github.com/vborjesson/TC_hunter
- **Clone**: `~/tools/TC_hunter/`
- **Status**: Ready (Python 2.7, Nextflow). Requires `NXF_VER=22.10.8` (DSL1 incompatible with nextflow ≥24.x).
- **Workflows**: `TC_hunter.nf` (BAM input), `TC_hunter_BWA.nf` (FASTQ input)
- **Note**: Python 2.7 required; `createOutput.py` has python3 shebang (may need fix); config uses `bwa_threads`. Requires exact construct FASTA.

### Env 4: RedGene (reference)
```bash
micromamba activate redgene
cd /data/gpfs/assoc/pgl/develop/redgene
```
- **Status**: Validated (rice, tomato, cucumber, corn, soybean)

---

## Benchmark Results

### Rice G281 (Chr3:16,439,674, 2 copies head-to-head)

| Tool | Detected | Position(s) | TP | FP | Notes | Runtime |
|------|----------|-------------|----|----|-------|---------|
| **RedGene** | Yes | Chr3:16,439,719, Chr3:31,443,557, Chr2:8,432,860 | 2 | 2 | Head-to-head confirmed; Chr2/Chr3:31M are host-element FPs (TA29, Act1) | ~15 min |
| **TDNAscan** | Yes | Chr3:16,439,710, Chr3:31,443,564, Chr10:13,500,286, Chr2:~8,433,216 | 2 | 2 | Chr10/Chr2 FPs; Chr3 positions match ground truth ±36 bp | ~40 min |
| **T-DNAreader** | Yes | Chr3:16,439,674 (R), Chr3:16,439,710 (L), Chr10:13,500,285, Chr11:12,123,278, Chr3:29,074,003 | 2 | 3 | Chr10/Chr11/Chr3:29M are FPs; best positional accuracy (0 bp offset) | ~30 min |
| **TC-hunter** | Yes | Chr3:16,439,675 (3 links), Chr3:16,439,710 (3 links) | 2 | ~4 | Supplementary reads also hit Chr10, Chr2:8,432,860, Chr3:31,443,452, Chr11 | ~35 min |

**Summary**: All 4 tools detect the true insertion site. All produce false positives from host-derived construct elements (Act1, TA29 promoters). RedGene and TDNAscan produce 2 FPs each; T-DNAreader 3 FPs; TC-hunter ~4 FP clusters from supplementary alignments.

### Tomato Cas9 A2_3 (chr08:65,107,378)

| Tool | Detected | Position(s) | TP | FP | Notes | Runtime |
|------|----------|-------------|----|----|-------|---------|
| **RedGene** | Yes | chr08:65,107,378 (T-DNA junction) + CRISPR edits (SlAMS, SlPHD_MS1) | 1 | 0 | Also detects 6bp del (SlAMS) and 9bp del (SlPHD_MS1) via s08 | ~25 min |
| **TDNAscan** | No | 0 hits (header-only output) | 0 | 0 | Failed to detect any insertion | ~45 min |
| **T-DNAreader** | Pending | — | — | — | Job 5623365 running | — |
| **TC-hunter** | No | 0 links, 0 supplementary links | 0 | 0 | Failed — no discordant pairs found | ~40 min |

### Cucumber Line 224 (Chr2, G6936 promoter, 1 copy)

| Tool | Detected | Position(s) | Notes | Runtime |
|------|----------|-------------|-------|---------|
| **RedGene** | Yes | LKUO03001512.1:580,628–581,332 (LB+RB) | 9 junctions, 28 bp resolution | ~20 min |
| **TDNAscan** | N/A | — | Construct unavailable (pRUR528s not deposited) | — |
| **T-DNAreader** | N/A | — | Construct unavailable | — |
| **TC-hunter** | N/A | — | Construct unavailable | — |

### Cucumber Line 212 (Chr6, intergenic, 1 copy)

| Tool | Detected | Position(s) | Notes | Runtime |
|------|----------|-------------|-------|---------|
| **RedGene** | Yes | 3 junctions detected | element_db mode | ~18 min |
| **TDNAscan** | N/A | — | Construct unavailable | — |
| **T-DNAreader** | N/A | — | Construct unavailable | — |
| **TC-hunter** | N/A | — | Construct unavailable | — |

### Cucumber Line 225 (Chr2, intergenic, 2 copies + backbone)

| Tool | Detected | Position(s) | Notes | Runtime |
|------|----------|-------------|-------|---------|
| **RedGene** | Yes | 2 junctions detected | element_db mode | ~18 min |
| **TDNAscan** | N/A | — | Construct unavailable | — |
| **T-DNAreader** | N/A | — | Construct unavailable | — |
| **TC-hunter** | N/A | — | Construct unavailable | — |

### Corn ND207 (~5x, Chr3:181,367,276)

| Tool | Detected | Position(s) | TP | FP | Notes | Runtime |
|------|----------|-------------|----|----|-------|---------|
| **RedGene** | Yes | NC_050098.1:181,367,276 (LB) | 1 | ~6 | Event-specific border match (0 bp offset). FPs from endogenous Ubi1, zSSIIb | ~2.5 h |
| **TDNAscan** | No | 0 hits (header-only output) | 0 | 0 | Uses element_db → construct not recognized | — |
| **T-DNAreader** | N/A | — | — | — | Construct unavailable (ND207 proprietary) | — |
| **TC-hunter** | N/A | — | — | — | Construct unavailable | — |

### Soybean UGT72E3 (Chr18:51,882,860, Glyma.18g226800 CDS)

Ground truth: Insertion in first exon of Glyma.18g226800 (= LOC102669442), Chr18.

| Tool | Detected | Position(s) | TP | FP | Notes | Runtime |
|------|----------|-------------|----|----|-------|---------|
| **RedGene** | Yes | NC_038254.2:51,882,860 + 51,882,903 | 1 | ~7 | Confirmed in Glyma.18g226800 CDS (51,880,693–51,883,761). Required `--min-coverage-frac 0.30` for element_db | ~45 min |
| **TDNAscan** | No | 0 hits (header-only output) | 0 | 0 | element_db → construct not recognized | — |
| **T-DNAreader** | N/A | — | — | — | Construct unavailable (pB7WG2-UGT72E3 not deposited) | — |
| **TC-hunter** | N/A | — | — | — | Construct unavailable | — |

### Soybean AtYUCCA6 (2 sites: Chr18 + Chr19)

Ground truth (Kim et al. 2021):
- **Site I**: Chr18, Glyma.18g235800 — 5+ tandem copies (multi-copy insertion)
- **Site II**: Chr19, Glyma.19g245800 (NC_038255.2:49,789,752)

| Tool | Detected | Position(s) | TP | FP | Notes | Runtime |
|------|----------|-------------|----|----|-------|---------|
| **RedGene** | Partial | NC_038255.2:49,789,752 (Site II) | 1/2 | ~3 | Site I (5+ tandem copies) defeats short-read assembly — SPAdes collapses repeats. Site II confirmed | ~45 min |
| **TDNAscan** | N/A | — | — | — | Construct unavailable | — |
| **T-DNAreader** | N/A | — | — | — | Construct unavailable | — |
| **TC-hunter** | N/A | — | — | — | Construct unavailable | — |

**Note on Site I**: Multi-copy tandem insertions (5+ copies in the same locus) are an inherent limitation of short-read assembly. SPAdes collapses identical tandem repeats into a single contig, preventing junction detection. Indirect detection via depth anomaly at Glyma.18g235800 using s07 BAM is pending.

---

## Coverage Sensitivity Comparison

### Corn ND207 (key junction: NC_050098.1:181,367,276)

| Coverage | RedGene | TDNAscan | T-DNAreader | TC-hunter |
|----------|---------|----------|-------------|-----------|
| ~5x (full) | Yes (7 junctions, High confidence, MAPQ 60) | N/A | N/A | N/A |
| 5x subsampled | Yes (5 junctions, High confidence) | N/A | N/A | N/A |
| 3x subsampled | Yes (11 junctions*, High confidence) | N/A | N/A | N/A |
| 1x subsampled | Yes (5 junctions, Medium confidence, MAPQ 36) | N/A | N/A | N/A |

*3x produces more junctions due to assembly fragmentation (more short chimeric contigs).

All coverages detect the key ND207 junction at NC_050098.1:181,367,276. At 1x, confidence drops to Medium with MAPQ 36.

### Cucumber Line 224

| Coverage | RedGene | TDNAscan | T-DNAreader | TC-hunter |
|----------|---------|----------|-------------|-----------|
| ~37x (full) | Yes (High, 9 junctions) | N/A | N/A | N/A |
| 15x | Yes (High, 14 junctions) | N/A | N/A | N/A |
| 10x | Yes (High, 14 junctions) | N/A | N/A | N/A |
| 5x | Partial (Medium, 1 junction) | N/A | N/A | N/A |
| 3x | Failed (0 junctions) | N/A | N/A | N/A |

### Rice G281

| Coverage | RedGene | TDNAscan | T-DNAreader | TC-hunter |
|----------|---------|----------|-------------|-----------|
| ~30x (full) | Yes | — | — | — |
| 15x | Yes | — | — | — |
| 10x | Partial | — | — | — |
| 5x | Partial | — | — | — |
| 3x | Failed | — | — | — |

---

## Key Findings

### 1. Construct availability is the primary barrier
For 4 of 6 test datasets (cucumber, corn, soybean×2), no public construct FASTA exists. TDNAscan, T-DNAreader, and TC-hunter are **not applicable** for these samples. RedGene's element database approach (131 EUginius elements + 62 corn LB/RB border sequences) is the only viable method.

### 2. All tools produce false positives from host-derived elements
For rice G281 (the one sample where all tools can run), every tool produces 2–4 false positives caused by construct elements derived from host genes (Act1, TA29 promoters). RedGene's WT-based homology filter (s03b) can mitigate this.

### 3. Read-mapping tools fail on tomato
Both TDNAscan and TC-hunter failed to detect the tomato Cas9 T-DNA insertion (0 hits), despite having the exact construct FASTA. RedGene's assembly-based approach detected it correctly.

### 4. Element_db requires lower thresholds
When using element_db instead of exact construct:
- `--min-identity 0.70` (vs default 0.90): minimap2 alignment identity to fragmented elements is ~0.84
- `--min-coverage-frac 0.30` (vs default 0.50): element_db doesn't cover full vector, so chimeric contig combined coverage is 45–49%

These are auto-applied by `run_pipeline.py` when element_db or combined_db is detected.

### 5. Multi-copy tandem insertions defeat short-read assembly
Soybean AtYUCCA6 Site I (5+ tandem copies at Glyma.18g235800) was not detected by any method. SPAdes collapses identical repeats → no chimeric contig formed. This is a fundamental limitation of short-read WGS, not tool-specific.

### 6. Coverage requirements
Assembly-based junction detection (RedGene) requires ≥10x for reliable detection, ≥5x for partial detection. At 1x (corn), the key junction is still detectable but at reduced confidence (MAPQ 36). At 3x, cucumber assembly fails completely.

---

## Notes

- All tools run on the same hardware: Pronghorn HPC, 16 threads, 64GB RAM
- Runtime measured wall-clock for junction detection steps only (excluding host full-genome mapping)
- "Distance" = |detected position − published/expected position|
- RedGene results from steps 1–6 only (assembly-based junction detection)
- TDNAscan, T-DNAreader, TC-hunter require exact construct/vector FASTA (not element DB)
- N/A = construct FASTA unavailable, tool cannot run
- T-DNAreader results for tomato/cucumber/soybean/corn are pending (job 5623365 still running — building host bowtie2 indices). Note: for samples without exact construct (cucumber, soybean), T-DNAreader was run with element_db as a best-effort test, but results may be unreliable.
