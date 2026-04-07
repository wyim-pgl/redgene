# RedGene GMO Element Database

## Database Files

| File | Sequences | Description |
|------|-----------|-------------|
| `element_db/gmo_combined_db.fa` | 131 | 131 EUginius elements + custom additions |
| `db/corn_border_db.fa` | 62 | 62 corn LB/RB event-specific border sequences |
| `element_db/gm_new_elements.fa` | 3 | 3 newly added GM regulatory elements |
| `db/gmo_all_combined_db.fa` | 196 | 196 total = all above combined (USE THIS) |
| `db/gmo_corn_combined_db.fa` | 192 | 192 = EUginius + corn borders (legacy, corn-specific) |

## Summary (196 total sequences, 131 kb)

### By Source

| Source | Sequences | Total bp | File |
|--------|-----------|----------|------|
| EUginius DB | 131 | 121,190 bp | `element_db/gmo_combined_db.fa` |
| Corn border DB (Sci Rep 2025) | 62 | 7,912 bp | `db/corn_border_db.fa` |
| Custom GenBank additions | 3 | 2,037 bp | `element_db/gm_new_elements.fa` |
| **Total (unified)** | **196** | **131,139 bp** | **`db/gmo_all_combined_db.fa`** |

### By Category (EUginius 131 seqs)

| Category | Count | Examples |
|----------|-------|---------|
| Promoters | 10 | P-35S-CaMV, P-e35S, P-FMV, P-nos, P-ocs, P-Ubi1, P-Act1, P-TA29, P-SSuAra, P-CsVMV |
| Terminators | 4 | T-nos, T-35S-CaMV, T-E9, T-35S-pCAMBIA |
| Coding sequences | 14 | nptII, bar, pat, hpt, aadA, cp4-epsps, cry1Ab/Ac/F, cry3Bb1, gus, barnase, barstar, thaumatin |
| Regulatory | 2 | LB-TDNA, RB-TDNA (border repeats) |
| Construct-specific | 9 | Junction/amplicon targets spanning multiple GM events |
| Element-specific | 52 | qPCR amplicon targets for individual elements |
| Event-specific | 23 | Event-junction qPCR targets (KMD1, DAS1507, MON87419, etc.) |
| Taxon-specific | 17 | Species ID targets (rice, corn, wheat, canola, salmon, etc.) |

### Corn Border Sequences (62 seqs)

31 approved corn events, each with LB + RB border sequences (~128 bp each):
3272, 5307, 59122, BFL4-2, BT506, Bt11, CBH351, CM8101, DAS-40278-9, DBN9501, DBN9858, DBN9936, DP4114, GA21, GAB-3, IE09S034, MIR098, MIR162, MIR604, MON810, MON87411, MON87419, MON87427, MON87460, MON88017, MON89034, ND207, RF125, T25, TC1507, ZZM030

Source: Table S1, Sci Rep 2025 (DOI: 10.1038/s41598-025-18593-8)

### Custom Additions (3 seqs)

| Name | GenBank | Size | Description |
|------|---------|------|-------------|
| P-FMV34S-full | X06166:6300-7743 | 1,444 bp | Full-length FMV 34S promoter |
| CTP4-PhEPSPS | M21084:28-243 | 216 bp | Petunia EPSPS chloroplast transit peptide |
| I-hsp70-maize | X03714:482-858 | 377 bp | Maize hsp70 intron 1 |

## Newly Added Elements (2025-04)

### 1. P-FMV34S (Figwort Mosaic Virus 34S Promoter)
- **GenBank**: X06166 (FMV strain DxS complete genome, 7743 bp)
- **Region**: 6300-7743 (1444 bp, full-length including enhancer)
- **Reference**: Richins et al. 1987 NAR; Sanger et al. 1990 Plant Mol Biol
- **Usage**: MON87705, MON88302, GT200, H7-1 etc. (often as P-FMV/Tsf1 chimeric)
- **Note**: Short amplicon (79bp) already in EUginius; this is the full-length version

### 2. P-Ubi1 (Maize Ubiquitin-1 Promoter)
- **GenBank**: S94464 (1994 bp, promoter + exon1 + intron1)
- **Reference**: Christensen & Quail 1996 Transgenic Res
- **Usage**: Widely used in monocot transformation (Bt corn events)
- **Status**: Already in EUginius DB (canonical_v1)

### 3. T-E9 (Pea rbcS E9 3' Terminator)
- **GenBank**: X00806 (region 1793-2351, 559 bp)
- **Reference**: Coruzzi et al. 1984 EMBO J
- **Usage**: MON89788 soybean (P-FMV/Tsf1 + CP4-EPSPS + T-E9)
- **Status**: Already in EUginius DB (canonical_v1)

### 4. P-Act1 (Rice Actin-1 Promoter)
- **GenBank**: S44221 (1442 bp, promoter + 5'UTR + intron1)
- **Reference**: McElroy et al. 1990 Plant Cell
- **Usage**: TA29 rice (host-derived, causes false positives in GMO detection)
- **Status**: Already in EUginius DB (canonical_v1)

### 5. CTP4 (Petunia EPSPS Chloroplast Transit Peptide)
- **GenBank**: M21084 (sig_peptide 28..243, 216 bp = 72 aa)
- **Reference**: Shah et al. 1986 Science; della-Cioppa et al. 1986 PNAS
- **Usage**: Fused N-terminally to CP4-EPSPS in Roundup Ready events
- **Status**: NEW — added to db

### 6. I-hsp70 (Maize Heat Shock Protein 70 Intron 1)
- **GenBank**: X03714 (intron1 at 482..858, 377 bp)
- **Reference**: Rochester et al. 1986 EMBO J
- **Usage**: MON810 (P-e35S + I-hsp70 + cry1Ab), enhances expression in monocots
- **Status**: NEW — added to db

### CTP2 (Arabidopsis EPSPS Chloroplast Transit Peptide)
- **GenBank**: Correct accession unclear (M14358 is E. coli flagellin, NOT AtEPSPS)
- **Known sequence**: 72 amino acids (216 bp), from Roundup Ready soybean GTS 40-3-2 petition
- **Reference**: Klee et al. 1987 Mol Gen Genet; Padgette et al. 1995
- **Status**: NOT ADDED — accession verification needed. CP4-EPSPS CDS is already in DB.

---

## EUginius Element List (131 sequences)

Source: [EUginius](https://euginius.eu) — EU Database of Reference Sequences for GMO Detection

### Promoters (10)

| # | Name | Organism | Accession | Description |
|---|------|----------|-----------|-------------|
| 1 | P-35S-CaMV | Cauliflower mosaic virus | V00141.1 | Functional 35S promoter region, -421 to +109 relative to TSS |
| 2 | P-e35S-CaMV | Cauliflower mosaic virus | V00141.1 | Enhanced/doubled 35S including upstream enhancer duplication |
| 3 | P-FMV | Figwort mosaic virus | NC_003554.1 | Functional FMV promoter region |
| 4 | P-nos | Agrobacterium tumefaciens | V00087.1 | Nopaline synthase promoter region |
| 5 | P-ocs | Agrobacterium tumefaciens | AJ311874.1 | Octopine synthase promoter region |
| 6 | P-SSuAra | Arabidopsis thaliana | X13611.1 | Small subunit of RuBisCO Ara promoter |
| 7 | P-TA29 | Nicotiana tabacum | X52283.1 | Tapetum-specific TA29 promoter |
| 8 | P-Ubi1-maize | Zea mays | S94464.1 | Promoter plus 5'UTR exon1 and intron1 |
| 9 | P-Act1-rice | Oryza sativa | S44221.1 | Promoter plus 5'UTR leader and intron |
| 10 | P-CsVMV | Cassava vein mosaic virus | U59751.1 | Cassava vein mosaic virus promoter |

### Terminators (4)

| # | Name | Organism | Accession | Description |
|---|------|----------|-----------|-------------|
| 1 | T-nos | Agrobacterium tumefaciens | V00087.1 | nos 3' UTR and transcription termination region |
| 2 | T-35S-CaMV | Cauliflower mosaic virus | V00141.1 | CaMV 35S polyadenylation and termination region |
| 3 | T-E9 | Pisum sativum | X00806.1 | Pea rbcS E9 3' terminator |
| 4 | T-35S-pCAMBIA | Cauliflower mosaic virus | AF234296.1 | pCAMBIA vector-derived 35S terminator variant |

### Coding Sequences (14)

| # | Name | Organism | Accession | Description |
|---|------|----------|-----------|-------------|
| 1 | nptII | Escherichia coli Tn5 | V00618.1 | Neomycin phosphotransferase II (kanamycin resistance) |
| 2 | bar | Streptomyces hygroscopicus | X17220.1 | Phosphinothricin acetyltransferase (glufosinate resistance) |
| 3 | pat | Streptomyces viridochromogenes | M22827.1 | Phosphinothricin acetyltransferase (glufosinate resistance) |
| 4 | hpt | Streptomyces hygroscopicus | U40398.1 | Hygromycin phosphotransferase |
| 5 | aadA | Escherichia coli | X02340.1 | Aminoglycoside adenyltransferase (spectinomycin resistance) |
| 6 | cp4-epsps | Agrobacterium sp. CP4 | L29358.1 | CP4 EPSPS synthase, excluding transit peptide |
| 7 | cry1Ab | Bacillus thuringiensis | M13898.1 | Crystal protein Cry1Ab (Bt insecticidal) |
| 8 | cry1Ac | Bacillus thuringiensis | M11068.1 | Crystal protein Cry1Ac (Bt insecticidal) |
| 9 | cry1F | Bacillus thuringiensis | M73254.1 | Crystal protein Cry1F (Bt insecticidal) |
| 10 | cry3Bb1 | Bacillus thuringiensis | M89794.1 | Crystal protein Cry3Bb1 (rootworm) |
| 11 | gus | Escherichia coli | S69414.1 | Beta-glucuronidase (reporter gene) |
| 12 | barnase | Bacillus amyloliquefaciens | M14442.1 | Ribonuclease barnase (male sterility) |
| 13 | barstar | Bacillus amyloliquefaciens | X15545.1 | Barnase inhibitor barstar (fertility restorer) |
| 14 | thaumatin II | Thaumatococcus daniellii | J01209.1 | Preprothaumatin-2 (sweetness protein, added for cucumber) |

### Regulatory Elements (2)

| # | Name | Organism | Accession | Description |
|---|------|----------|-----------|-------------|
| 1 | LB-TDNA | Agrobacterium tumefaciens | AF485783.1 | Canonical T-DNA left border repeat |
| 2 | RB-TDNA | Agrobacterium tumefaciens | AF485783.1 | Canonical T-DNA right border repeat |

### Construct-Specific (9)

| # | ID | GM Events | Accession | Size |
|---|------|----------|-----------|------|
| 1 | QL-CON-00-008 | GT200, GT73, GTSB77, H7-1, J101, J163, etc. | FN550387.1 | 88 bp |
| 2 | QL-CON-00-012 | 55-1, 63-1 | KR076687.1 | 896 bp |
| 3 | QL-CON-00-013 | DAS1507, DAS59122, KMD1, Kefeng6 | EUginius | 90 bp |
| 4 | QL-CON-00-014 | 1345-4, 16-0-1, 18-2-4, 55-1, 63-1, etc. | MF521566.1 | 3304 bp |
| 5 | QL-CON-00-015 | CBH351, CM8101, DBT418, DLL25, etc. | JX139719.1 | 985 bp |
| 6 | QL-CON-307-F_652-R | KMD1, Kefeng6 | EUginius | 175 bp |
| 7 | QL-CON-35S-F_nptII-R | 23-18-17, 5345, 8338, ATBT04 series | KT184682.1 | 888 bp |
| 8 | QL-CON-CP1_pFMV2 | GT200, GT73, MON87705, MON88302, etc. | BT022026.1 | 1804 bp |
| 9 | QL-CON-gat_T-pinII | 61061, 73496, 98140, DP356043 | KP784699.1 | 917 bp |

### Element-Specific (52)

| # | ID | Target | Accession | Size |
|---|------|--------|-----------|------|
| 1 | QL-ELE-00-001 | P-35S-CaMV | P-35S-CaMV | 195 bp |
| 2 | QL-ELE-00-002 | nptII | nptII | 215 bp |
| 3 | QL-ELE-00-003 | nptII | nptII | 173 bp |
| 4 | QL-ELE-00-004 | P-35S-CaMV | P-35S-CaMV | 123 bp |
| 5 | QL-ELE-00-006 | T-nos | P-nos | 180 bp |
| 6 | QL-ELE-00-008 | P-nos/P-ocs | P-ocs | 94 bp |
| 7 | QL-ELE-00-009 | T-nos | P-nos | 118 bp |
| 8 | QL-ELE-00-010 | P-FMV | P-FMV | 196 bp |
| 9 | QL-ELE-00-011 | T-nos | EUginius | 84 bp |
| 10 | QL-ELE-00-012 | P-35S/E-35S | EUginius | 82 bp |
| 11 | QL-ELE-00-014 | bar | EUginius | 60 bp |
| 12 | QL-ELE-00-015 | P-FMV | P-FMV | 78 bp |
| 13 | QL-ELE-00-016 | cry1Ab/1Ac | Y09787.1 | 1857 bp |
| 14 | QL-ELE-00-017 | P-35S-CaMV | P-35S-CaMV | 75 bp |
| 15 | QL-ELE-00-018 | T-nos | P-nos | 69 bp |
| 16 | QL-ELE-00-022 | bar | X17220.1 | 69 bp |
| 17 | QL-ELE-00-023 | T-35S-CaMV | EUginius | 137 bp |
| 18 | QL-ELE-00-024 | T-E9 | EUginius | 87 bp |
| 19 | QL-ELE-00-025 | pat | JX139722.1 | 493 bp |
| 20 | QL-ELE-AgroBS1 | Ti-nos flanking | EUginius | 76 bp |
| 21 | QL-ELE-AgroBS2 | Ti-ocs flanking | HF950132.1 | 1130 bp |
| 22 | QL-ELE-AINT | Actin intron | X63830.1 | 623 bp |
| 23 | QL-ELE-bar_2-5_2-3 | bar | EUginius | 186 bp |
| 24 | QL-ELE-Bnase | barnase | EUginius | 102 bp |
| 25 | QL-ELE-Bstar | barstar | EUginius | 145 bp |
| 26 | QL-ELE-cat-F1_R1 | cat (chloramphenicol) | NG_047580.1 | 851 bp |
| 27 | QL-ELE-cat-F2_R2 | cat | EUginius | 529 bp |
| 28 | QL-ELE-cat-F_R | cat | EUginius | 96 bp |
| 29 | QL-ELE-cry1A.105 | cry1A.105 synthetic | EUginius | 133 bp |
| 30 | QL-ELE-cry1A_4-5_4-3 | cry1Ab | AF465640.1 | 204 bp |
| 31 | QL-ELE-Cry1Ac | cry1Ac | HQ161057.1 | 1867 bp |
| 32 | QL-ELE-Cry1F | cry1F truncated | EUginius | 90 bp |
| 33 | QL-ELE-cry2Ab2 | cry2Ab2 | EUginius | 121 bp |
| 34 | QL-ELE-cry3A | mcry3A/cry3A | EUginius | 117 bp |
| 35 | QL-ELE-Cry3Bb | cry3Bb1 | EUginius | 232 bp |
| 36 | QL-ELE-epsps_1-5_3-3 | CP4-EPSPS | AF464188.1 | 1529 bp |
| 37 | QL-ELE-epsps2 | CP4-EPSPS | EUginius | 118 bp |
| 38 | QL-ELE-npt | nptII | EUginius | 155 bp |
| 39 | QL-ELE-Pat | pat | EUginius | 144 bp |
| 40 | QL-ELE-PFMV_1-5_1-3 | P-FMV34S | EUginius | 110 bp |
| 41 | QL-ELE-P-FMV | P-FMV34S | EUginius | 79 bp |
| 42 | QL-ELE-P-Rice_actin | P-act1 rice | EUginius | 95 bp |
| 43 | QL-ELE-P-SSuAra | P-rbcS Arabidopsis | XM_002888530.2 | 1170 bp |
| 44 | QL-ELE-P-TA29 | P-TA29 tobacco | X52283.1 | 6254 bp |
| 45 | QL-ELE-P-Ubi | P-Ubi1 maize | EUginius | 76 bp |
| 46 | QL-ELE-RERIO | DsRed2 fluorescent | EUginius | 506 bp |
| 47 | QL-ELE-T-35S | T-35S CaMV | FJ154952.1 | 1632 bp |
| 48 | QL-ELE-T-OCS | T-ocs | EUginius | 85 bp |
| 49 | QL-TAX-AT-001 | nos (Agro) | OP712404.1 | 1133 bp |
| 50 | QT-ELE-00-001 | P-35S-CaMV | P-35S-CaMV | 79 bp |
| 51 | QT-ELE-NOS_ter | T-nos | EUginius | 151 bp |
| 52 | QT-ELE-P35S | P-35S-CaMV | EUginius | 100 bp |

### Event-Specific (23)

| # | ID | GM Event | Accession | Size |
|---|------|----------|-----------|------|
| 1 | QL-EVE-16-0-1RB | 16-0-1 | KU376441.1 | 234 bp |
| 2 | QL-EVE-18-2-4RB | 18-2-4 | KU376439.1 | 309 bp |
| 3 | QL-EVE-319FalF | Falcon GS 40/90 | AY096774.1 | 216 bp |
| 4 | QL-EVE-375LibF | Liberator | EUginius | 96 bp |
| 5 | QL-EVE-558F | e871 B. subtilis | LT622644.1 | 2668 bp |
| 6 | QL-EVE-816AvaF | Falcon GS 40/90 | XM_009128664.3 | 1834 bp |
| 7 | QL-EVE-CP-001 | Huanong No. 1 | EUginius | 285 bp |
| 8 | QL-EVE-Kef6 | Kefeng6 | HM124448.1 | 506 bp |
| 9 | QL-EVE-MDB110 | OXY-235 | EU099579.1 | 1109 bp |
| 10 | QL-EVE-MS-J101 | J101 | EUginius | 102 bp |
| 11 | QL-EVE-MS-J163 | J163 | EUginius | 118 bp |
| 12 | QL-EVE-MS-KK179 | KK179 | HF949564.1 | 1158 bp |
| 13 | QL-EVE-OS-KM2 | KMD1 | EU980363.1 | 556 bp |
| 14 | QL-EVE-Oxy | OXY-235 | EUginius | 124 bp |
| 15 | QL-EVE-SS-AquAd | AquAdvantage Salmon | S65567.1 | 3350 bp |
| 16 | QT-EVE-CP-001 | Huanong No. 1 | EUginius | 175 bp |
| 17 | QT-EVE-DBN9501 | DBN9501 | EUginius | 96 bp |
| 18 | QT-EVE-GH-012 | COT102 | EUginius | 101 bp |
| 19 | QT-EVE-HF-1F | Huafan No. 1 | KC820135.1 | 753 bp |
| 20 | QT-EVE-OS-002 | LLRICE62 | JQ406881.1 | 628 bp |
| 21 | QT-EVE-OS-G6H1 | G6H1 | EUginius | 90 bp |
| 22 | QT-EVE-ZM-010 | DAS1507 | EUginius | 58 bp |
| 23 | QT-EVE-ZM-029 | MON87419 | HF949548.1 | 1133 bp |

### Taxon-Specific (17)

| # | ID | Species | Accession | Size |
|---|------|---------|-----------|------|
| 1 | QL-TAX-Acc1 | Medicago sativa (alfalfa) | L25042.1 | 7175 bp |
| 2 | QL-TAX-CaMVF | Cauliflower mosaic virus | HE978785.1 | 1391 bp |
| 3 | QL-TAX-CMV-001 | Cauliflower mosaic virus | V00141.1 | 68 bp |
| 4 | QL-TAX-CP-001 | Carica papaya | AY803756.1 | 719 bp |
| 5 | QL-TAX-DC-001 | Dianthus caryophyllus (carnation) | AB727362.2 | 4604 bp |
| 6 | QL-TAX-FMVorf7 | Figwort mosaic virus | NC_003554.1 | 7743 bp |
| 7 | QL-TAX-PS-001 | Pisum sativum (pea) | HM346476.1 | 284 bp |
| 8 | QL-TAX-SS-AM5F | Salmo salar (Atlantic salmon) | AY614009.1 | 4673 bp |
| 9 | QT-TAX-BV-001 | Beta vulgaris (sugar beet) | AY026353.1 | 1296 bp |
| 10 | QT-TAX-CP-001 | Carica papaya | EUginius | 74 bp |
| 11 | QT-TAX-FatA | Brassica napus (canola) | XM_010437010.2 | 1476 bp |
| 12 | QT-TAX-OS-002 | Oryza sativa (rice) | XM_052303253.1 | 2740 bp |
| 13 | QT-TAX-SPS | Oryza sativa (rice) | KT225496.1 | 721 bp |
| 14 | QT-TAX-wx012 | Triticum aestivum (wheat) | LC373576.1 | 4431 bp |
| 15 | QT-TAX-ZM-003 | Zea mays (corn) | K03285.1 | 360 bp |
| 16 | QT-TAX-ZM-004 | Zea mays (corn) | EU963452.1 | 928 bp |
| 17 | QT-TAX-zSSIIb | Zea mays (corn) | AF019297.1 | 2480 bp |

---

## Source Databases

- **EUginius** (https://euginius.eu): 131 GMO detection elements (promoters, terminators, CDS, construct-specific, event-specific, taxon-specific)
- **Sci Rep 2025 corn borders**: 62 LB/RB event-specific border sequences from 31 approved corn events (DOI: 10.1038/s41598-025-18593-8)
- **Custom additions**: Full-length P-FMV34S, CTP4, hsp70 intron from GenBank

## Usage

```bash
# Use the unified database for all species
python run_pipeline.py --sample my_sample --steps 1-6 --threads 16
# In config.yaml, set construct_reference to:
#   db/gmo_all_combined_db.fa
```
