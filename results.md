# GMO Positive Control Characterization Pipeline - Results

## Rice G281 (Oryza sativa cv. Xiushui-110)

### Dataset
- **SRA**: SRR18236702 (BioProject PRJNA812718)
- **Paper**: Communications Biology, 2022 - paired-end WGS transgene characterization
- **Construct**: p1300-S450RNAi-hFL-G6 (pCAMBIA-1300 backbone)
- **Platform**: Illumina HiSeq2000, PE100
- **Total reads**: 108.3M reads (54.1M pairs, ~29x coverage)
- **Host reference**: Osativa_323_v7.0.fa (374 Mbp, Phytozome)
- **Known insertion**: Chr3:16,439,674 (two T-DNA copies, head-to-head)

### Pipeline Results (Full Coverage ~29x)

| Step | Status | Key Results |
|------|--------|-------------|
| Step 1 (QC) | DONE | 108.3M -> 105.7M reads (97.6% pass), Q30: 92.8% |
| Step 2 (Construct map) | DONE | 7,128 reads mapped (0.01%), 40.1% construct covered |
| Step 3 (Read extract) | DONE | 6,029 read pairs extracted |
| Step 4 (Assembly) | DONE | 569 contigs, N50=504bp, longest=2,703bp |
| Step 5 (Contig map) | DONE | 553 host contigs, 8 construct contigs |
| Step 6 (Junction) | DONE | 4 junctions detected (1 true, 3 false positive) |
| Step 7 (Host map) | RUNNING | SLURM job 5618730 (bwa mem, 105.7M reads) |
| Step 10 (Copy number) | PENDING | Awaiting step 7 |

### Junction Detection Results

| Contig | Host Position | MAPQ | Element | Type | Verdict |
|--------|--------------|------|---------|------|---------|
| NODE_14 | **Chr3:16,439,719** | 10 | pCAMBIA-1300 backbone | RB | **TRUE** (45bp from known) |
| NODE_14 | Chr3:31,443,557 | 35 | pCAMBIA-1300 backbone | RB | FALSE POSITIVE |
| NODE_6 | Chr3:31,443,557 | 60 | pCAMBIA-1300 backbone | RB | FALSE POSITIVE |
| NODE_6 | Chr2:8,432,860 | 60 | pCAMBIA-1300 backbone | RB | FALSE POSITIVE |

**Critical finding**: MAPQ alone is insufficient for false positive removal.
Chr2:8,432,860 has MAPQ=60 (unique mapping) but is a false positive from
pCAMBIA vector homology with the rice genome. WT-based filtering is essential.

### Homologous Regions Detected (Construct vs Host)
Minimap2-based detection between pCAMBIA-1300/element_db and rice genome:
- Chr01: 17 homologous region pairs
- Chr02: 90 pairs (includes 8,432,860 false positive region)
- Chr03: 9 pairs (includes 31,443,557 false positive region)
- Chr10: 141 pairs (includes 13,500,269 false positive region)

### Coverage Depth Sensitivity

| Coverage | True Junction (Chr3:16.4M) | False Positives | Junction Contigs |
|----------|---------------------------|----------------|------------------|
| ~29x (full) | **DETECTED** (45bp accuracy) | 3 | NODE_6 (1369bp), NODE_14 (349bp) |
| ~15x | **DETECTED** (45bp accuracy) | 6+ | NODE_14 (336bp) |
| ~10x | **NOT DETECTED** | 8 | Assembly fragmentation |
| ~5x | **NOT DETECTED** | 2 | Too few reads |
| ~3x | **NOT DETECTED** | 3 | Unreliable |

**Conclusion**: Minimum ~15x coverage required for reliable assembly-based junction detection.

---

## Tomato Micro-Tom (PRJNA692070)

### Dataset (Bae et al. 2022, Horticultural Science and Technology)
CRISPR/Cas9-edited Micro-Tom tomato targeting SlFAD2 gene:

| Sample | SRA | Library | Cas9 Status | Edit Type | Raw Reads |
|--------|-----|---------|-------------|-----------|-----------|
| WT | SRR13450615 | MT_WT | N/A | None | 238.7M |
| Cas9 A2_1 | SRR13450618 | 20_A2_1 | Removed | Heterozygous | 299.8M |
| Cas9 A2_2 | SRR13450617 | 20_A2_2 | Present | Homozygous | 62.3M |
| Cas9 A2_3 | SRR13450616 | 20_A2_3 | Present | Heterozygous | 120.5M |

- **Host reference**: SLM_r2.0.pmol.fasta (Micro-Tom, plantgarden.jp, 832.8 Mbp, 12 chr)
- **Construct reference**: element_db/gmo_combined_db.fa (NOT vector-specific)
- **Note**: Identity threshold lowered to 0.70 (element_db has only partial homology to actual Cas9 vector)

### Cas9 A2_3 (SRR13450616) - Cas9 present, heterozygous edit

| Step | Key Results |
|------|-------------|
| QC | 120.5M -> 114.4M reads (95.0%), ~10.3x coverage |
| Construct map | 3,658 reads mapped (0.01%) |
| Read extract | 3,436 read pairs |
| Assembly | 501 contigs, N50=347, longest=4,160bp |
| Contig map | 5,774 host + 15 construct alignments (4 contigs) |
| Junction | **1 junction: SLM_r2.0ch08:65,107,378** (RB, MAPQ=60) |

**Junction details:**
- NODE_10 (1031bp, cov=4.1): Chimeric contig spanning host chr08:65,107,378 and construct element QL-CON-00-014 (MF521566.1)
- NODE_1 (4160bp, cov=8.6): Pure construct (CaMV 35S, nos promoter, nptII, OCS terminator) - T-DNA body

### Cas9 A2_2 (SRR13450617) - Cas9 present, homozygous edit

| Step | Key Results |
|------|-------------|
| QC | 62.3M -> 59.4M reads (95.4%), ~5.4x coverage |
| Construct map | 3,119 reads mapped (0.01%) |
| Read extract | 2,951 read pairs |
| Assembly | 355 contigs, N50=320, longest=3,350bp |
| Contig map | 1,281 host + 17 construct alignments (5 contigs) |
| Junction | **No junctions detected** |

**Analysis:**
- 3 chimeric contigs found but all had overlapping host/construct alignments
  (homologous regions: TA29 tobacco promoter naturally present in tomato genome)
- NODE_1 (3350bp) and NODE_6 (1190bp): Pure construct contigs confirming
  Cas9 T-DNA presence but no junction-spanning assembly
- **Root cause**: Low coverage (~5.4x) insufficient for junction assembly
  (consistent with rice finding: minimum 15x needed)

### Cas9 A2_1 (SRR13450618) - Cas9 removed, heterozygous edit
- **Status**: Download pending

### WT (SRR13450615) - Wild type control
- **Status**: Download completing (compression in progress)
- Will be used for WT-based filtering of all transgenic samples

---

## Key Findings

### 1. Assembly-Based Detection Works for T-DNA Insertions
- Successfully detects junction coordinates within 45bp accuracy (rice G281)
- Detects Cas9 T-DNA insertion site in tomato (chr08:65,107,378)
- Requires sufficient coverage (>=15x) for reliable results

### 2. False Positive Mitigation is Critical
- **WT-based filtering** (s03b_homology_filter.py): Most effective method
- **MAPQ filtering**: Necessary but insufficient (MAPQ=60 false positives exist)
- **Structural warnings**: Inter-chromosomal and distant intra-chromosomal junctions flagged
- Plant genomes have extensive construct-host homology (especially host-derived promoters)

### 3. Element_db vs Specific Vector
- Using element_db requires lowered identity threshold (0.70 vs 0.90)
- Advantage: No prior knowledge of vector sequence needed
- Disadvantage: More false positive risk, reduced detection sensitivity

### 4. Coverage Requirements
| Application | Minimum Coverage | Notes |
|------------|-----------------|-------|
| Reliable junction detection | 15x | Assembly-based approach |
| Junction detection possible | 10x | May miss some junctions |
| Construct presence only | 5x | Can confirm T-DNA but not locate insertion |
| Unreliable | 3x | Insufficient for any confident result |

---

## File Locations

### Pipeline Scripts
- Main entry: `run_pipeline.py`
- Step scripts: `scripts/s01_qc.py` through `scripts/s07_host_map.py`
- Homology filter: `scripts/s03b_homology_filter.py` (WT-based filtering)
- Junction verify: `scripts/s06b_junction_verify.py` (multi-evidence scoring)
- Zygosity: `scripts/s06c_zygosity.py` (homo/hetero estimation)
- Copy number: `scripts/s10_copynumber.py`
- Config: `config.yaml`
- Environment: micromamba env `redgene` (`/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene`)

### Reference Data
- Rice genome: `db/Osativa_323_v7.0.fa`
- Tomato genome: `db/SLM_r2.0.pmol.fasta`
- Element DB: `element_db/gmo_combined_db.fa` (131 sequences)

### Test Data
- Rice: `test_data/rice_G281_R{1,2}.fastq.gz`
- Rice subsampled: `test_data/subsampled/rice_G281_{15x,10x,5x,3x}_R{1,2}.fq.gz`
- Tomato: `test_data/tomato/SRR134506{15,16,17,18}_{1,2}.fastq.gz`

### Results
- Rice (full): `results/rice_G281/`
- Rice (coverage tests): `results/rice_G281_{15x,10x,5x,3x}/`
- Tomato: `results/tomato_Cas9_A2_{1,2,3}/`, `results/tomato_WT/`
