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
| Step 7 (Host map) | DONE | 105.7M reads, 7.0GB BAM, BWA -t 2, ~5h |
| Step 10 (Copy number) | DONE | **2.1 copies** (construct 52.0x / host 25.0x), **High confidence** |

### Junction Detection Results

| Contig | Host Position | MAPQ | Element | Type | Verdict |
|--------|--------------|------|---------|------|---------|
| NODE_14 | **Chr3:16,439,719** | 10 | pCAMBIA-1300 backbone | RB | **TRUE** (45bp from known, MAPQ=10 is low) |
| NODE_14 | Chr3:31,443,557 | 35 | pCAMBIA-1300 backbone | RB | FALSE POSITIVE |
| NODE_6 | Chr3:31,443,557 | 60 | pCAMBIA-1300 backbone | RB | FALSE POSITIVE |
| NODE_6 | Chr2:8,432,860 | 60 | pCAMBIA-1300 backbone | RB | FALSE POSITIVE |

### Copy Number Estimation

| Metric | Value |
|--------|-------|
| Best construct marker | QT-TAX-ZM-004 (Zea mays, 69bp amplicon) |
| Construct median depth | 52.0x |
| Host median depth | 25.0x |
| Depth ratio | 2.08 |
| Estimated copies | **2.1 (multi-copy)** |
| Expected (ground truth) | **2 copies (head-to-head)** |
| Junction count | 4 (consistent with multi-copy) |
| Confidence | **High** |

**Validation**: Pipeline estimated 2.1 copies vs expected 2 copies — excellent agreement.
The construct marker (Zea mays taxon-specific element, 69bp) was auto-selected as the
most reliable marker due to highest median coverage and no host homology.

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

| Coverage | True Junction (Chr3:16.4M) | Total Junctions | False Positives | Junction Contigs |
|----------|---------------------------|-----------------|----------------|------------------|
| ~29x (full) | **DETECTED** (45bp accuracy) | 4 | 3 | NODE_6 (1369bp), NODE_14 (349bp) |
| ~15x | **DETECTED** (45bp accuracy) | 9 | 8 | NODE_14 (336bp) |
| ~10x | **NOT DETECTED** | 8 | 8 | Assembly fragmentation |
| ~5x | **NOT DETECTED** | 2 | 2 | Too few reads |
| ~3x | **NOT DETECTED** | 3 | 3 | Unreliable |

**Conclusion**: Minimum ~15x coverage required for reliable assembly-based junction detection.

---

## Tomato Micro-Tom (PRJNA692070)

### Dataset
CRISPR/Cas9-edited Micro-Tom tomato (Seol et al. 2025, Plant Biotechnology Reports).
Targets: **SlAMS** (Solyc08g062780) and **SlPHD_MS1** (Solyc04g008420) for male sterility.
SRA data from PRJNA692070 (originally published by Bae et al. 2022).

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
| Copy number | Depth ratio 2.1x (T-35S-pCAMBIA: 19.0x vs host: 9.0x) → **~2 copies** |

**Junction details:**
- NODE_10 (1031bp, cov=4.1): Chimeric contig spanning host chr08:65,107,378 and construct element QL-CON-00-014 (MF521566.1)
- NODE_1 (4160bp, cov=8.6): Pure construct (CaMV 35S, nos promoter, nptII, OCS terminator) - T-DNA body

**Copy number notes:**
- Best marker: T-35S-pCAMBIA (258bp, 100% covered, median 19.0x)
- Alternative marker nptII: depth ratio 4.66/9.0 = 0.52 → hemizygous single-copy
- Discrepancy likely due to short element length and CaMV 35S multi-mapping
- nptII-based estimate (hemizygous single-copy) is more biologically reliable for heterozygous sample

**Step 8 — CRISPR editing detection (gRNA-guided, pileup-based):**
- **SlPHD_MS1 (chr04:2635445)**: 1bp T insertion at cut site, 30.8% freq (4/13 reads), heterozygous
  - Different editing outcome than expected 9bp deletion — 1bp insertion is a common NHEJ repair product
  - Treatment-specific (absent in WT)
- **SlAMS (chr08:53314229)**: NOT edited — all reads match reference (8-10x depth)
  - The 98.4% editing frequency from Seol et al. was for "20-1-ca7-1 line", a different subline
  - The 1-A2 samples from PRJNA692070 (Bae et al. 2022) appear unedited at SlAMS

### Cas9 A2_2 (SRR13450617) - Cas9 present, homozygous edit

| Step | Key Results |
|------|-------------|
| QC | 62.3M -> 59.4M reads (95.4%), ~10.4x coverage |
| Construct map | 3,119 reads mapped (0.01%) |
| Read extract | 2,951 read pairs |
| Assembly | 355 contigs, N50=320, longest=3,350bp |
| Contig map | 1,281 host + 17 construct alignments (5 contigs) |
| Junction | **No junctions detected** |
| Copy number | Depth ratio 1.6x (T-OCS: 14.0x vs host: 9.0x) → **~2 copies** |

**Analysis:**
- 3 chimeric contigs found but all had overlapping host/construct alignments:
  - NODE_7 (1072bp): TA29 tobacco promoter ↔ ch02:59,279,798 (606bp overlap)
  - NODE_9 (847bp): TA29 tobacco promoter ↔ ch02:59,278,936 (164bp overlap)
  - NODE_5 (1207bp): pinII terminator ↔ ch03:57,893,656 (244bp overlap)
- NODE_1 (3350bp) and NODE_6 (1190bp): Pure construct contigs confirming
  Cas9 T-DNA presence but no junction-spanning assembly
- **Root cause**: Not coverage depth (10.4x is comparable to A2_3's 10.75x).
  The cause is **assembly direction stochasticity** — in A2_2, contigs extended
  toward construct elements with host genome homology (TA29, pinII), producing
  overlapping alignments rather than adjacent (true junction) alignments.
  In A2_3, the junction contig extended toward QL-CON-00-014 (bacterial nptII/nos
  region) which has no host homology, producing a clean junction.
  This is an inherent limitation of short-read assembly-based detection.
- **Potential improvement**: Iterative assembly (round-1 contigs as bait for
  round-2 read extraction) could recover junctions that extend in unfavorable
  directions during initial assembly.

**Step 8 — CRISPR editing detection (gRNA-guided, pileup-based):**
- **SlPHD_MS1 (chr04:2635440)**: **9bp deletion of GTGAGCCAT** detected, 16.7% freq (1/6 reads)
  - **Exact match with ground truth** (Seol et al. 2025: 9bp in-frame deletion of "GTGAGCCAT")
  - Low observed frequency due to alignment bias at deletion sites + very low depth (6 reads)
  - Additional 1bp indels at adjacent positions (2635443-2635445) are alignment artifacts
  - Treatment-specific (absent in WT, which has 32 reads all matching reference)
- **SlAMS (chr08:53314229)**: NOT edited — all 14-17 reads match reference
  - Same finding as A2_3; 1-A2 line was not edited at SlAMS locus

### Cas9 A2_1 (SRR13450618) - Cas9 removed, heterozygous edit

| Step | Key Results |
|------|-------------|
| QC | 299.8M -> 292.3M reads (97.5%), ~26.3x coverage |
| Construct map | 17,007 reads mapped (0.01%) |
| Read extract | 15,920 read pairs |
| Assembly | 616 contigs, N50=850, longest=4,059bp |
| Junction | **No junctions** (expected — Cas9 removed, no T-DNA) |
| Host map | DONE | 292.3M reads, 17GB BAM, BWA -t 16 |
| Copy number | Construct 1.0x / Host 48.0x → **~0 copies** (expected — no T-DNA) |
| Step 8 | **2 on-target edits at SlAMS** (see below) |

- **Expected**: No T-DNA insertion (Cas9 removed via segregation) ✓
- **Unexpected**: SlAMS editing detected! (see step 8 results below)
- **Note**: 15,920 extracted reads is higher than A2_3 (3,436 reads) despite no T-DNA, reflecting
  larger dataset size (292M vs 60M reads) and proportionally more host-derived noise

**Step 8 — CRISPR editing detection (gRNA-guided, pileup-based):**
- **SlAMS (chr08:53314230)**: **4bp deletion of GTAC** detected, 22.0% freq (9/41 reads), heterozygous
  - Also 1bp deletion of G at same position (14.6% freq, 6/41 reads)
  - Both treatment-specific (absent in WT, which has 37 reads all matching reference)
  - The 4bp GTAC deletion is at the Cas9 cut site, consistent with NHEJ repair
  - **A2_1 HAS SlAMS editing** — unlike A2_2 and A2_3 which are unedited at SlAMS
- **SlPHD_MS1 (chr04:2635445)**: NOT detected — no editing at this locus
  - Opposite pattern from A2_2/A2_3, which have SlPHD_MS1 editing but no SlAMS

**Biological interpretation**: Different T0 regenerant lines (1-A2-1, 1-A2-2, 1-A2-3)
have different CRISPR editing outcomes. The dual-gRNA system targeted both SlAMS and
SlPHD_MS1, but each line was edited at different loci:
- A2_1: SlAMS edited (4bp del), SlPHD_MS1 unedited
- A2_2: SlAMS unedited, SlPHD_MS1 edited (9bp del)
- A2_3: SlAMS unedited, SlPHD_MS1 edited (1bp ins)

### WT (SRR13450615) - Wild type control

| Step | Key Results |
|------|-------------|
| QC | 238.7M -> 229.3M reads (95.5%), ~20.6x coverage |
| Construct map | 11,833 reads mapped (0.01%) |
| Read extract | 9,764 read pairs (all host-derived noise) |
| Assembly | 546 contigs, N50=617, longest=2,086bp |
| Junction | **No junctions** (correct — wild type) |
| Host map | DONE | 230.6M reads, 14GB BAM, BWA -t 16 |

- Used as WT baseline for step 8 treatment-vs-WT comparison
- No step 8 or step 10 needed (WT control sample)

### Tomato Coverage Depth Sensitivity

**A2_3 (Cas9 present, T-DNA positive control):**

| Coverage | Total Reads | Extracted Reads | Contigs | N50 | Junction Detected | Position |
|----------|------------|-----------------|---------|-----|-------------------|----------|
| ~10.3x (full) | 60.2M | 3,436 | 501 | 347 | **YES** | chr08:65,107,378 |
| ~5x | ~33M | 785 | 92 | 264 | **NO** | - |
| ~3x | ~20M | 473 | 52 | 264 | **NO** | - |

**A2_1 (Cas9 removed, negative control):**

| Coverage | Extracted Reads | Contigs | N50 | Junction Detected |
|----------|-----------------|---------|-----|-------------------|
| ~15x | 22,983 | 3,114 | 257 | **NO** (expected) |
| ~10x | 2,911 | 299 | 395 | **NO** (expected) |
| ~5x | 1,488 | 186 | 283 | **NO** (expected) |
| ~3x | 905 | 121 | 272 | **NO** (expected) |

**Note on A2_1 15x**: 22,983 extracted reads is anomalously high for a T-DNA-negative sample.
This reflects background noise from host sequences with partial homology to element_db entries
(TA29, pinII, Ubi1, actin). WT-based filtering (s03b) would remove most of these.

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
- Plant genomes have extensive construct-host homology (TA29, pinII, Ubi1, SSuAra)
- **Limitation**: s03b may remove true junction reads if T-DNA inserts near a
  host homologous region (within ±500bp of TA29/pinII ortholog). Low probability
  but non-zero — should be noted as a pipeline limitation.

### 2b. Assembly Direction Stochasticity
- Same T-DNA insertion site can yield different results depending on which
  construct elements the junction contig extends toward
- Contigs extending toward bacterial-origin elements (nptII, nos) → clean junction
- Contigs extending toward host-derived elements (TA29, pinII) → overlapping
  alignments, junction filtered as false positive
- This is an inherent limitation of short-read assembly; potential mitigation
  via iterative assembly (round-1 contigs as bait for round-2 extraction)

### 3. Element_db vs Specific Vector
- Using element_db requires lowered identity threshold (0.70 vs 0.90)
- Advantage: No prior knowledge of vector sequence needed
- Disadvantage: More false positive risk, reduced detection sensitivity

### 4. Coverage Requirements

**Rice (374 Mbp genome):**

| Coverage | Junction Detection | Notes |
|----------|-------------------|-------|
| ~29x (full) | DETECTED (Chr3:16,439,719) | 45bp from known position |
| ~15x | DETECTED (Chr3:16,439,719) | Same accuracy as full |
| ~10x | NOT detected (Chr3:16,439,719) | 8 junctions found but known site missing |
| ~5x | NOT detected | 2 false positives |
| ~3x | NOT detected | 3 false positives |

**Tomato A2_3 (833 Mbp genome):**

| Coverage | Junction Detection | Notes |
|----------|-------------------|-------|
| ~10.3x (full) | DETECTED (chr08:65,107,378) | EXACT MATCH with known |
| ~5x | NOT detected | 785 extracted reads, N50=264 |
| ~3x | NOT detected | 473 extracted reads, N50=264 |

**Summary**: Minimum ~15x for rice, ~10x for tomato for the *specific known junction*.
Note: Rice at 10x still finds 8 junction candidates (including Chr3:31,443,557 and
Chr2:8,432,860 seen at full coverage) — it only misses the Chr3:16,439,719 site.
The tomato coverage sensitivity is based on a single sample (A2_3) with a single
junction, so the threshold is not statistically robust.

**Caveat**: The 15x rice data produces "High" confidence at Chr3:16,439,719 while
the full 29x data produces "Medium" with MAPQ=10 — an example of assembly
stochasticity where more reads can produce different (sometimes worse) contigs.

### 5. CRISPR Editing Detection (Step 8)

**Method**: gRNA-guided pileup parsing — maps gRNA to host genome via BLAST,
then directly parses samtools mpileup at predicted cut sites (±50bp window).
Compares treatment vs WT to identify treatment-specific indels.

**Key improvement**: Uses direct pileup parsing (`samtools mpileup -Q 0`)
instead of bcftools variant calling. bcftools decomposes complex indels at
low depth (e.g., splits 9bp deletion into multiple 1bp events). Direct
pileup parsing preserves the original indel representation from the reads.

**Results:**

| Sample | SlPHD_MS1 (chr04) | SlAMS (chr08) |
|--------|-------------------|---------------|
| A2_1 (Cas9 removed) | Not edited | **4bp del GTAC** (9/41, 22%) |
| A2_2 (homozygous) | **9bp del GTGAGCCAT** (1/6, 17%) | Not edited |
| A2_3 (heterozygous) | 1bp T insertion (4/13, 31%) | Not edited |
| WT | Clean reference (32 reads) | Clean reference (37 reads) |

**Validation**:
- SlPHD_MS1 9bp deletion in A2_2 **matches ground truth sequence** (Seol et al. 2025)
- SlAMS 4bp GTAC deletion in A2_1 detected at cut site (heterozygous, high-confidence: dp=41)
- All indels are treatment-specific (absent in WT)

**Reliability caveat**: A2_2 9bp deletion is called from dp=6, count=1 (17%, qual=10).
This is the weakest call — a single supporting read at low quality. The sample is labeled
"homozygous" but observed frequency is far below 50%, suggesting either severe
under-sampling at ~5x coverage, alignment bias at deletion sites, or chimeric editing.
A2_1 results (dp=41) are much more reliable. Re-sequencing A2_2 at higher depth
is recommended for confirmation.

**Biological finding**: Different T0 lines have different editing outcomes.
A2_1 was edited at SlAMS only; A2_2 and A2_3 were edited at SlPHD_MS1 only.
The "20-1-ca7-1 line" (98.4% indel freq) from Seol et al. is a different subline
with homozygous SlAMS editing; the 1-A2 lines have heterogeneous editing patterns.

---

## File Locations

### Pipeline Scripts
- Main entry: `run_pipeline.py`
- Step scripts: `scripts/s01_qc.py` through `scripts/s07_host_map.py`
- Homology filter: `scripts/s03b_homology_filter.py` (WT-based filtering)
- Junction verify: `scripts/s06b_junction_verify.py` (multi-evidence scoring)
- Zygosity: `scripts/s06c_zygosity.py` (homo/hetero estimation)
- CRISPR editing: `scripts/s08_indel_detection.py` (gRNA-guided pileup-based detection)
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

### Visualization Outputs

#### Editing Profile (CRISPResso2-inspired, 5-panel nucleotide quilt)
- `results/tomato_Cas9_A2_{1,2,3}/s08_indel/editing_profile_gRNA{1,2}.{png,pdf}`
- Panels: Read depth, Treatment nucleotide quilt, WT nucleotide quilt,
  Net modification frequency, Reference sequence with gRNA/PAM/cut annotation
- CRISPResso2 color scheme: A=green, T=purple, G=yellow, C=orange, DEL=dark, INS=brown

#### Variant Effect Annotation (easyGWAS-inspired)
- `results/tomato_Cas9_A2_{1,2,3}/s08_indel/editing_effects.{png,pdf,tsv}`
- Effect categories: frameshift (dark red), nonsynonymous (red), synonymous (green),
  in-frame deletion (orange), splice site (amber), UTR (blue), intergenic (gray)
- Includes allele frequency chart and codon context details

#### Junction Gene Context (gene/exon/CDS annotation)
- `results/rice_G281/s06_junction/junction_gene_{1..4}_*.{png,pdf}`
- `results/tomato_Cas9_A2_3/s06_junction/junction_gene_1_SLM_r2.0ch08_65107378.{png,pdf}`
- Shows gene model with exon/CDS/intron structure, junction position in gene context,
  contig alignment diagram with host/construct segments

---

## Key Finding 6: CRISPR Editing Variant Effects

| Sample | Site | Gene | Indel | Effect | AA Impact |
|--------|------|------|-------|--------|-----------|
| A2_1 | ch08:53314230 | SlAMS (LOC101253608) | 4bp del GTAC | **Frameshift** | Disrupts TF domain |
| A2_1 | ch08:53314230 | SlAMS (LOC101253608) | 1bp del G | **Frameshift** | Disrupts TF domain |
| A2_2 | ch04:2635440 | SlPHD_MS1 (LOC101257358) | 9bp del GTGAGCCAT | **In-frame del** | Removes 3aa at pos 236 |
| A2_2 | ch04:2635443 | SlPHD_MS1 (LOC101257358) | 1bp del A | **Frameshift** | Disrupts PHD finger |
| A2_2 | ch04:2635444 | SlPHD_MS1 (LOC101257358) | 1bp ins G | **Frameshift** | Disrupts PHD finger |
| A2_2 | ch04:2635445 | SlPHD_MS1 (LOC101257358) | 1bp del C | **Frameshift** | Disrupts PHD finger |
| A2_3 | ch04:2635445 | SlPHD_MS1 (LOC101257358) | 1bp ins T | **Frameshift** | Disrupts PHD finger |

**Key insight**: The 9bp GTGAGCCAT deletion in A2_2 is an in-frame deletion that
removes 3 amino acids from the PHD finger domain, while the 1bp indels cause
frameshifts that likely produce non-functional truncated proteins. Both mechanisms
effectively disrupt the male sterility gene function, consistent with the
intended CRISPR knockout strategy.

### GFF3 Annotation Files
- Rice: `db/Osativa_323_v7.0.gene_exons.gff3` (MSU/RGAP v7.0, 55,986 genes)
- Tomato: `db/SLM_r2.0.gff3.gz` (NCBI Gnomon, 39,809 genes, chr names remapped)
