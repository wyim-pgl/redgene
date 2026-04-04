# RedGene - GMO/LMO Transgene Characterization Pipeline

Assembly-based Illumina WGS pipeline for characterizing transgenic plant
insertions. Designed for Korean quarantine GMO/LMO detection assays and
CRISPR/Cas9 editing verification.

## What It Does

Given paired-end Illumina reads from a transgenic plant:

1. **Finds transgene insertion sites** at base-pair resolution
2. **Estimates copy number** from read depth ratios
3. **Determines zygosity** (homozygous vs heterozygous insertion)
4. **Detects CRISPR editing** (treatment-specific indels vs WT)
5. **Filters false positives** from construct-host sequence homology
6. **Warns about structural anomalies** (translocations, large deletions)

## Background

Plant T-DNA constructs often contain host-derived promoters (e.g., rice Ubi1,
maize Act1, tobacco TA29) that create false positive chimeric reads when mapped
against construct databases. This pipeline implements multiple filtering
strategies including WT-based filtering and MAPQ analysis to distinguish true
transgene junctions from homologous sequence artifacts.

The approach is based on the modified TranSeq method (Bae et al. 2022,
Communications Biology) adapted for generic GMO element databases rather than
specific vector sequences.

## Pipeline Architecture

```
                    Illumina PE150 Reads
                           |
                     [Step 1] fastp QC
                           |
                  [Step 2] bwa mem → construct BAM
                           |
              [Step 3] Extract construct-hitting reads + mates
                           |
          [Step 3b] WT-based homology filter (optional, recommended)
                           |
              [Step 4] SPAdes local assembly → contigs
                           |
      [Step 5] minimap2 contigs → host & construct PAF
                           |
      [Step 6] Chimeric contig detection → junction coordinates
           |                              |
  [Step 6b] Junction verify     [Step 6c] Zygosity estimation
  (multi-evidence scoring)      (allele fraction analysis)
                           |
             [Step 7] bwa mem → host BAM (full genome)
                    |                    |
        [Step 8] CRISPR indel     [Step 10] Copy number
        detection (vs WT)          (depth ratio)
```

### Step Details

| Step | Script | Input | Output | Description |
|------|--------|-------|--------|-------------|
| 1 | `s01_qc.py` | Raw FASTQ | Trimmed FASTQ | QC + adapter trimming (fastp) |
| 2 | `s02_construct_map.py` | Trimmed FASTQ | BAM | Map reads to construct/element DB |
| 3 | `s03_extract_reads.py` | Construct BAM | FASTQ pairs | Extract construct-hitting reads + mates |
| 3b | `s03b_homology_filter.py` | Extracted FASTQ | Filtered FASTQ | Remove reads from homologous regions |
| 4 | `s04_assembly.py` | Extracted FASTQ | Contigs FASTA | Local assembly with SPAdes |
| 5 | `s05_contig_map.py` | Contigs | PAF files | Map contigs to host + construct (minimap2) |
| 6 | `s06_junction.py` | PAF files | junctions.tsv | Detect chimeric contigs, extract junction coords |
| 6b | `s06b_junction_verify.py` | junctions.tsv + host BAM | Verified junctions | Multi-evidence scoring |
| 6c | `s06c_zygosity.py` | junctions.tsv + host BAM | Zygosity report | Allele fraction analysis |
| 7 | `s07_host_map.py` | Trimmed FASTQ | Host BAM | Full genome mapping for depth analysis |
| 8 | `s08_indel_detection.py` | Host BAMs (treatment + WT) | editing_sites.tsv | CRISPR editing detection |
| 10 | `s10_copynumber.py` | Host BAM + construct BAM | Copy number | Depth ratio estimation |

## Quick Start

```bash
# 1. Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

# 2. Prepare references
#    - Host genome: db/host.fa (+ BWA index: bwa index db/host.fa)
#    - Element DB: element_db/gmo_combined_db.fa (included, 131 GMO elements)

# 3. Run pipeline (example: rice G281)
SAMPLE=rice_G281
HOST=db/Osativa_323_v7.0.fa
CONSTRUCT=element_db/gmo_combined_db.fa
R1=test_data/rice_G281_R1.fastq.gz
R2=test_data/rice_G281_R2.fastq.gz

# Core pipeline: steps 1-6
python scripts/s01_qc.py --r1 $R1 --r2 $R2 --outdir results --sample-name $SAMPLE
python scripts/s02_construct_map.py \
  --r1 results/$SAMPLE/s01_qc/${SAMPLE}_R1.fq.gz \
  --r2 results/$SAMPLE/s01_qc/${SAMPLE}_R2.fq.gz \
  --construct-ref $CONSTRUCT --outdir results --sample-name $SAMPLE
python scripts/s03_extract_reads.py \
  --bam results/$SAMPLE/s02_construct_map/${SAMPLE}_construct.bam \
  --outdir results --sample-name $SAMPLE
python scripts/s04_assembly.py \
  --r1 results/$SAMPLE/s03_extract/${SAMPLE}_construct_R1.fq.gz \
  --r2 results/$SAMPLE/s03_extract/${SAMPLE}_construct_R2.fq.gz \
  --outdir results --sample-name $SAMPLE
python scripts/s05_contig_map.py \
  --contigs results/$SAMPLE/s04_assembly/contigs.fasta \
  --host-ref $HOST --construct-ref $CONSTRUCT \
  --outdir results --sample-name $SAMPLE
python scripts/s06_junction.py \
  --host-paf results/$SAMPLE/s05_contig_map/${SAMPLE}_contigs_to_host.paf \
  --construct-paf results/$SAMPLE/s05_contig_map/${SAMPLE}_contigs_to_construct.paf \
  --contigs results/$SAMPLE/s04_assembly/contigs.fasta \
  --outdir results --sample-name $SAMPLE

# For CRISPR samples: add step 7 + step 8
python scripts/s07_host_map.py \
  --r1 results/$SAMPLE/s01_qc/${SAMPLE}_R1.fq.gz \
  --r2 results/$SAMPLE/s01_qc/${SAMPLE}_R2.fq.gz \
  --host-ref $HOST --outdir results --sample-name $SAMPLE
python scripts/s08_indel_detection.py \
  --treatment-bam results/$SAMPLE/s07_host_map/${SAMPLE}_host.bam \
  --wt-bam results/WT/s07_host_map/WT_host.bam \
  --host-ref $HOST --outdir results --sample-name $SAMPLE
```

## Element Database

The `element_db/gmo_combined_db.fa` contains 131 GMO-related sequences scraped
from the EUginius database (https://euginius.eu/), including:

- **Promoters**: CaMV 35S, nos, Ubi1 (maize), Act1 (rice), TA29 (tobacco), etc.
- **Terminators**: nos, OCS, CaMV 35S, pinII, etc.
- **Selection markers**: nptII, bar, hpt, pat, etc.
- **Construct-specific**: pCAMBIA vectors, event-specific amplicons
- **Full-length sequences**: Reference constructs from NCBI

To rebuild the database: `python element_db/gmo_db.py`

## False Positive Filtering

Plant genomes contain sequences homologous to common construct elements.
This creates false positive junction calls that require filtering:

### Problem
| Source | Example | MAPQ | Solution |
|--------|---------|------|----------|
| Host-derived promoters | Rice Ubi1 in pCAMBIA | Low (0-10) | MAPQ filter |
| Unique homologous regions | Chr2:8.4M in rice G281 | High (60) | WT filter |
| Repetitive elements | Multiple chromosomes | Low (0-5) | MAPQ filter |

### Solution: Multi-layer filtering
1. **MAPQ filter** (step 6, `--min-host-mapq 10`): Removes multi-mapping artifacts
2. **WT-based filter** (step 3b): Removes reads mapping to construct-host homologous regions
3. **Structural warnings** (step 6): Flags inter-chromosomal and distant junctions
4. **Verification** (step 6b): Multi-evidence scoring with composite verdicts

### Key insight
**MAPQ alone is insufficient.** In rice G281 testing, Chr2:8,432,860 had MAPQ=60
(unique mapping) but was a false positive from pCAMBIA vector homology. WT-based
filtering is the only reliable method for these cases.

## CRISPR Editing Detection (Step 8)

For Cas9/CRISPR-edited samples, step 8 detects editing-induced indels:

1. **Genome-wide variant calling** in both treatment and WT
2. **Subtraction**: Keep only treatment-specific indels (>= 2bp)
3. **PAM motif check**: Look for NGG within 3-8bp of indel
4. **NHEJ signature**: Filter for typical CRISPR editing patterns
   - Deletions 1-20bp (most common)
   - Insertions 1-5bp
5. **Zygosity**: Homozygous, heterozygous, or biallelic

Note: CRISPR editing sites are typically at the guide RNA target gene,
NOT at the T-DNA insertion site. The search is genome-wide by default.

## Environment

```bash
# Create environment
micromamba create -n redgene -c conda-forge -c bioconda \
  python=3.11 bwa minimap2 samtools fastp spades bcftools \
  sra-tools pysam matplotlib biopython seqkit pigz
```

### Key tools
| Tool | Version | Purpose |
|------|---------|---------|
| bwa | 0.7.18 | Read mapping to construct/host |
| minimap2 | 2.28 | Contig mapping + homology detection |
| samtools | 1.21 | BAM processing |
| bcftools | 1.23 | Variant calling (indel detection) |
| fastp | 0.23.4 | QC + adapter trimming |
| SPAdes | 4.0.0 | Local assembly of extracted reads |
| sra-tools | 3.1.1 | SRA download (prefetch + fasterq-dump) |
| pysam | 0.22.1 | BAM parsing (Python) |
| matplotlib | 3.9 | Visualization |
| biopython | 1.84 | Sequence handling |

## Test Results

### Rice G281 (T-DNA insertion, ~29x coverage)
- **True insertion**: Chr3:16,439,719 detected (45bp from known Chr3:16,439,674)
- **False positives**: 3 detected, all from pCAMBIA vector-host homology
- **Coverage minimum**: 15x for reliable detection

### Tomato Micro-Tom Cas9 (CRISPR editing, PRJNA692070)
- **T-DNA insertion**: SLM_r2.0ch08:65,107,378 detected (A2_3 sample, ~10x)
- **Cas9 construct**: CaMV 35S + nos + nptII + OCS identified in assembled contigs
- **Coverage minimum**: 10x for junction detection with element_db

### Coverage Sensitivity

| Coverage | Junction Detection | Notes |
|----------|-------------------|-------|
| >= 15x | Reliable, high confidence | Recommended |
| 10-15x | Usually works | May need relaxed thresholds |
| 5-10x | Construct presence confirmed | Junction assembly may fail |
| < 5x | Unreliable | Insufficient for assembly approach |

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-identity` | 0.90 | Min alignment identity. Lower to 0.70 for element_db |
| `--min-coverage-frac` | 0.50 | Min contig coverage. Lower to 0.40 for element_db |
| `--min-host-mapq` | 10 | Min host MAPQ (uniqueness filter) |
| `--min-host-aln-len` | 50 | Min host alignment length in bp |
| `--min-construct-aln-len` | 50 | Min construct alignment length in bp |
| `--min-indel-size` | 2 | Min indel size for CRISPR detection |

## Output Structure

```
results/<sample>/
  s01_qc/              # Trimmed reads + fastp HTML/JSON reports
  s02_construct_map/   # Construct BAM + flagstat + depth stats
  s03_extract/         # Extracted read pairs (+ filtered if s03b)
  s04_assembly/        # SPAdes contigs + assembly stats
  s05_contig_map/      # PAF alignments (host + construct)
  s06_junction/        # junctions.tsv + junction_contigs.fa
  s07_host_map/        # Host BAM + flagstat + depth stats
  s08_indel/           # editing_sites.tsv (CRISPR samples only)
  s10_copynumber/      # Copy number estimates
```

## References

- Bae S, Park YC, et al. (2022). Characterization of transgene insertions
  by resequencing. *Communications Biology* 5:671.
- Bae S, et al. (2022). CRISPR/Cas9-mediated editing of SlFAD2.
  *Horticultural Science and Technology* 40:81-92.
- EUginius GMO detection methods: https://euginius.eu/

## License

MIT
