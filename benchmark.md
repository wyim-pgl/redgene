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

### Test Datasets

| Sample | Species | Genome Size | Coverage | Expected Insertion | Source |
|--------|---------|-------------|----------|-------------------|--------|
| rice_G281 | *Oryza sativa* | 374 Mbp | ~30x | Chr3:16,439,674 (2 copies, head-to-head) | Phytozome |
| tomato_Cas9_A2_3 | *S. lycopersicum* Micro-Tom | 833 Mbp | ~25x | chr08:65,107,378 (T-DNA + CRISPR edits) | PRJNA692070 |
| cucumber_line224 | *Cucumis sativus* | 332 Mbp | ~37x | Chr2, promoter of G6936, 1 T-DNA copy | PRJNA638559 |
| cucumber_line212 | *Cucumis sativus* | 332 Mbp | ~37x | Chr6, intergenic, 1 T-DNA copy | PRJNA638559 |
| cucumber_line225 | *Cucumis sativus* | 332 Mbp | ~37x | Chr2, intergenic, 2 T-DNA copies + backbone | PRJNA638559 |
| corn_ND207 | *Zea mays* | 2.18 Gbp | ~5x | Chr3:181,367,276 (Bt construct) | CRA026358 |
| soybean_AtYUCCA6 | *Glycine max* | 1.1 Gbp | ~28x | Unknown (bar marker) | PRJNA627303 |
| soybean_UGT72E3 | *Glycine max* | 1.1 Gbp | ~29x | Unknown (bar marker) | PRJNA627303 |

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
- **Note**: bowtie2 index prefix must match FASTA path; `--gzip` for .gz input; `-g` is bowtie2 index prefix

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
- **Status**: Ready (Python 2.7, Nextflow). Test BAM not included (Google Drive download).
- **Workflows**: `TC_hunter.nf` (BAM input), `TC_hunter_BWA.nf` (FASTQ input)
- **Note**: Python 2.7 required; `createOutput.py` has python3 shebang (may need fix); config uses `bwa_threads`

### Env 4: RedGene (reference)
```bash
micromamba activate redgene
cd /data/gpfs/assoc/pgl/develop/redgene
```
- **Status**: Validated (rice, tomato, cucumber, corn)

---

## Benchmark Results

### Example Data Test (TDNAscan mt4_chr1 simulated)

| Tool | Detected | Position | Distance from Truth | Runtime |
|------|----------|----------|---------------------|---------|
| TDNAscan | — | — | — | — |
| T-DNAreader | — | — | — | — |
| TC-hunter | — | — | — | — |
| RedGene | — | — | — | — |

### Rice G281 (Chr3:16,439,674)

| Tool | Detected | Position(s) | Distance (bp) | Copy # | Confidence | Runtime |
|------|----------|-------------|---------------|--------|------------|---------|
| RedGene | Yes | Chr3:16,439,719 | 45 bp | 2 (head-to-head) | Medium | ~15 min |
| TDNAscan | — | — | — | — | — | — |
| T-DNAreader | — | — | — | — | — | — |
| TC-hunter | Yes | Chr3:16,439,675 + Chr3:31,443,564 | 1 bp | 2 (head-to-head) | 3+22 reads | ~35 min (BWA) |

### Cucumber Line 224 (Chr2, G6936 promoter)

| Tool | Detected | Position(s) | Distance (bp) | Confidence | Runtime |
|------|----------|-------------|---------------|------------|---------|
| RedGene | Yes | LKUO03001512.1:580,628–581,332 (LB+RB) | 28 bp | High | ~20 min |
| TDNAscan | — | — | — | — | — |
| T-DNAreader | — | — | — | — | — |
| TC-hunter | — | — | — | — | — |

### Corn ND207 (Chr3:181,367,276)

| Tool | Detected | Position(s) | Distance (bp) | Confidence | Runtime |
|------|----------|-------------|---------------|------------|---------|
| RedGene | Yes | NC_050098.1:181,367,276 (LB) | 0 bp (event-specific border match) | High | ~2.5 h |
| TDNAscan | — | — | — | — | — |
| T-DNAreader | — | — | — | — | — |
| TC-hunter | — | — | — | — | — |

### Tomato Cas9 A2_3 (chr08:65,107,378)

| Tool | Detected | Position(s) | Notes | Runtime |
|------|----------|-------------|-------|---------|
| RedGene | Yes (via s08) | T-DNA on chr08, CRISPR edits at SlAMS/SlPHD_MS1 | Assembly-based + pileup | ~25 min |
| TDNAscan | — | — | — | — |
| T-DNAreader | — | — | — | — |
| TC-hunter | — | — | — | — |

### Soybean AtYUCCA6 / UGT72E3

| Tool | Sample | Detected | Position(s) | Notes | Runtime |
|------|--------|----------|-------------|-------|---------|
| RedGene | AtYUCCA6 | — | — | — | — |
| RedGene | UGT72E3 | — | — | — | — |
| TDNAscan | — | — | — | — | — |
| T-DNAreader | — | — | — | — | — |
| TC-hunter | — | — | — | — | — |

---

## Coverage Sensitivity Comparison

### Corn ND207

| Coverage | RedGene | TDNAscan | T-DNAreader | TC-hunter |
|----------|---------|----------|-------------|-----------|
| ~5x (full) | Yes (High) | — | — | — |
| 5x subsampled | (pending) | — | — | — |
| 3x subsampled | (pending) | — | — | — |
| 1x subsampled | Yes (Medium, 1 true + 1 FP, MAPQ 36) | — | — | — |

### Cucumber Line 224

| Coverage | RedGene | TDNAscan | T-DNAreader | TC-hunter |
|----------|---------|----------|-------------|-----------|
| ~37x (full) | Yes (High, 9 junctions) | — | — | — |
| 15x | Yes (High, 14 junctions) | — | — | — |
| 10x | Yes (High, 14 junctions) | — | — | — |
| 5x | Partial (Medium, 1 junction) | — | — | — |
| 3x | Failed (0 junctions) | — | — | — |

### Rice G281

| Coverage | RedGene | TDNAscan | T-DNAreader | TC-hunter |
|----------|---------|----------|-------------|-----------|
| ~30x (full) | Yes | — | — | — |
| 15x | Yes | — | — | — |
| 10x | Partial | — | — | — |
| 5x | Partial | — | — | — |
| 3x | Failed | — | — | — |

---

## Notes

- All tools run on the same hardware: Pronghorn HPC, 16 threads, 64GB RAM
- Runtime measured wall-clock for junction detection steps only (excluding host full-genome mapping)
- "Distance" = |detected position - published/expected position|
- RedGene results from steps 1-6 only (assembly-based junction detection)
- TDNAscan, T-DNAreader, TC-hunter require exact construct/vector FASTA (not element DB)
