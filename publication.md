# Test Dataset Publications

## 1. GM Rice G281 (RNAi Event)

- **Paper**: Xu et al. 2022, "A paired-end whole-genome sequencing approach enables comprehensive characterization of transgene integration in rice"
- **Journal**: Communications Biology, 5(1):667
- **DOI**: 10.1038/s42003-022-03608-1
- **PMID**: 35790849 / **PMCID**: PMC9256713
- **Institution**: Shanghai Jiao Tong University

### SRA Data

| Sample | Accession | BioProject | Platform | Layout | Depth |
|--------|-----------|------------|----------|--------|-------|
| GM (G281) | SRR18236702 | PRJNA812718 | HiSeq 2000 | PE100 | 28.9x |
| WT (Xiushui 110) | SRR6983858 (SRX3923908) | PRJNA448500 | HiSeq 2000 | PE90 | 29.3x |

### Known Results (Ground Truth)

- **Host**: Oryza sativa cv. Xiushui 110
- **Insertion site**: Chr 3, position 16,439,674
- **Genomic alteration**: 36 bp deletion at insertion site
- **Copy structure**: Head-to-head 2-copy T-DNA
- **Backbone integration**: None detected
- **Construct**: p1300-S450RNAi-hLF-G6
  - hLF (human lactoferrin) with rice Gt1 promoter + PEPC terminator
  - G6 EPSPS (glyphosate tolerance) with maize ubiquitin promoter + PEPC terminator
  - CYP81A6 RNAi cassette with CaMV35S promoter

### Notes

- Construct sequence NOT deposited in GenBank; must reconstruct from paper supplementary
- WT control from different BioProject (sequenced earlier, 2018)
- PE100 shorter than our target PE150, but 29x coverage compensates
- Excellent ground truth for pipeline validation (known site, known copy number, known backbone status)

## 2. GM Tomato Gene-Editing Line (CRISPR/Cas9)

- **Paper**: Bae et al. 2022, "Establishment of a rapid assay for sequencing of carried DNA and edited sites in gene-editing tomato plants"
- **Journal**: Hortic Environ Biotechnol, 63:515-521
- **DOI**: 10.1007/s13580-022-00427-5
- **BioProject**: PRJNA692070
- **Institution**: National Institute of Agricultural Sciences, RDA, Jeonju, Korea

### SRA Data

| Sample | Accession | Platform | Layout | Depth | Size |
|--------|-----------|----------|--------|-------|------|
| Micro-Tom WT | SRR13450615 | HiSeq 3000 | PE151 | ~38x | 36.0 Gb |
| 20_A2_1 (transgenic) | SRR13450618 | HiSeq 3000 | PE151 | ~48x | 45.3 Gb |
| 20_A2_2 (transgenic) | SRR13450617 | HiSeq 3000 | PE151 | ~10x | 9.4 Gb |
| 20_A2_3 (transgenic) | SRR13450616 | HiSeq 3000 | PE151 | ~10x | 9.6 Gb |

### Key Info

- **Host**: Solanum lycopersicum (Micro-Tom, M82)
- **Event**: CRISPR/Cas9 gene-editing tomato (target: MS1-like, MS1035 male sterility genes)
- **Construct**: pAGM4723 (Cas9 + sgRNA), sequence available via Addgene #48015
- **Editing results**: bi-allelic indels at target sites (e.g., 5bp del + 1bp ins in MS1-like)

### Ground Truth Status

- **Insertion site**: NOT reported in paper
- **Copy number**: NOT reported in paper
- Paper focus is rapid screening assay (BioCube) validation, not insertion characterization
- No construct sequence deposited to NCBI; pAGM4723 available from Addgene

### Notes

- Korean research team — directly relevant to Korean quarantine context
- PE151 read length matches our target PE150
- Useful as de novo discovery test (no ground truth to compare)
- Priority: 2nd (after G281 rice which has full ground truth)
