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

### Primary Reference (CRISPR targets)
- **Paper**: Seol et al. 2025, "CRISPR/Cas9-mediated mutagenesis of SlAMS and SlPHD_MS1 for male sterility in tomato"
- **Journal**: Plant Biotechnology Reports
- **DOI**: 10.1007/s11816-025-00952-0
- **PMCID**: PMC11897730
- **CRISPR targets**: SlAMS (Solyc08g062780), SlPHD_MS1 (Solyc04g008420)

### SRA Data Reference
- **Paper**: Bae et al. 2022, "Establishment of a rapid assay for sequencing of carried DNA and edited sites in gene-editing tomato plants"
- **Journal**: Hortic Environ Biotechnol, 63:515-521
- **DOI**: 10.1007/s13580-022-00427-5
- **BioProject**: PRJNA692070
- **Institution**: National Institute of Agricultural Sciences, RDA, Jeonju, Korea

### SRA Data

| Sample | Accession | Library | Platform | Layout | Cas9 | Edit |
|--------|-----------|---------|----------|--------|------|------|
| Micro-Tom WT | SRR13450615 | MT_WT | HiSeq 3000 | PE151 | N/A | None |
| 20_A2_1 | SRR13450618 | 20_A2_1 | HiSeq 3000 | PE151 | Removed | Hetero |
| 20_A2_2 | SRR13450617 | 20_A2_2 | HiSeq 3000 | PE151 | Present | Homo |
| 20_A2_3 | SRR13450616 | 20_A2_3 | HiSeq 3000 | PE151 | Present | Hetero |

### Key Info

- **Host**: Solanum lycopersicum cv. Micro-Tom
- **Reference genome**: SLM_r2.0.pmol.fasta (plantgarden.jp, 832.8 Mbp, 12 chr)
- **Event**: CRISPR/Cas9 gene-editing for male sterility
- **Editing targets**:
  - **SlAMS** (Solyc08g062780) - Aborted Microspores homolog
    - gRNA+PAM: `GGTGTGCTATAAGTACTGAACGG`
    - Known edit: 6bp homozygous deletion "GTACTT"
  - **SlPHD_MS1** (Solyc04g008420) - PHD-finger MS1 homolog
    - gRNA+PAM: `CATTAGCCTATGGTGAGCCATGG`
    - Known edit: 9bp homozygous deletion "GTGAGCCAT"
- **Construct**: pAGM4723 binary vector with Cas9, nptII, CaMV 35S promoter, AtU6 promoter
- **Our finding**: T-DNA insertion at chr08:65,107,378 (A2_3 sample)

### Ground Truth Status

- **Insertion site**: NOT reported in papers; our pipeline detected chr08:65,107,378
- **Copy number**: NOT reported
- **Editing sites**: SlAMS on chr08, SlPHD_MS1 on chr04 — to be verified with step 8

### Notes

- Korean research team — directly relevant to Korean quarantine context
- PE151 read length matches our target PE150
- Useful as de novo discovery test (no ground truth for T-DNA position)
- CRISPR editing sites can be validated against known target genes
