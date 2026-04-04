# RedGene - GMO Positive Control Characterization Pipeline

Illumina PE150 sequencing-based pipeline to characterize transgenic plant
positive control materials for Korean quarantine GMO/LMO detection assays.

## Quick Start

```bash
# 1. Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

# 2. Prepare references in db/
#    - Host genome FASTA (+ BWA index)
#    - Element DB: element_db/gmo_combined_db.fa

# 3. Configure samples in config.yaml

# 4. Run pipeline step by step
python scripts/s01_qc.py --r1 data/R1.fq.gz --r2 data/R2.fq.gz --outdir results --sample-name my_sample
python scripts/s02_construct_map.py --r1 results/my_sample/s01_qc/my_sample_R1.fq.gz --r2 results/my_sample/s01_qc/my_sample_R2.fq.gz --construct-ref element_db/gmo_combined_db.fa --outdir results --sample-name my_sample
python scripts/s03_extract_reads.py --bam results/my_sample/s02_construct_map/my_sample_construct.bam --outdir results --sample-name my_sample
# Optional: WT-based homology filter
python scripts/s03b_homology_filter.py --construct-ref element_db/gmo_combined_db.fa --host-ref db/host.fa --r1 results/my_sample/s03_extract/my_sample_construct_R1.fq.gz --r2 results/my_sample/s03_extract/my_sample_construct_R2.fq.gz --outdir results --sample-name my_sample
python scripts/s04_assembly.py --r1 results/my_sample/s03_extract/my_sample_construct_R1.fq.gz --r2 results/my_sample/s03_extract/my_sample_construct_R2.fq.gz --outdir results --sample-name my_sample
python scripts/s05_contig_map.py --contigs results/my_sample/s04_assembly/contigs.fasta --host-ref db/host.fa --construct-ref element_db/gmo_combined_db.fa --outdir results --sample-name my_sample
python scripts/s06_junction.py --host-paf results/my_sample/s05_contig_map/my_sample_contigs_to_host.paf --construct-paf results/my_sample/s05_contig_map/my_sample_contigs_to_construct.paf --contigs results/my_sample/s04_assembly/contigs.fasta --outdir results --sample-name my_sample
```

## Environment

- **Manager**: micromamba 2.5.0
- **Env name**: `redgene`
- **Env path**: `/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene`
- **Python**: 3.11.15

### Key tools
| Tool | Version | Purpose |
|------|---------|---------|
| bwa | 0.7.18 | Read mapping |
| minimap2 | 2.28 | Contig mapping + homology detection |
| samtools | 1.21 | BAM processing |
| fastp | 0.23.4 | QC + trimming |
| SPAdes | 4.0.0 | Local assembly |
| sra-tools | 3.1.1 | SRA download |
| pysam | 0.22.1 | BAM parsing (Python) |
| matplotlib | 3.9 | Visualization |
| biopython | 1.84 | Sequence handling |

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | s01_qc.py | QC + adapter trimming (fastp) |
| 2 | s02_construct_map.py | Map reads to construct/element DB (bwa mem) |
| 3 | s03_extract_reads.py | Extract construct-hitting reads + mates |
| 3b | s03b_homology_filter.py | WT-based homologous sequence filter (minimap2) |
| 4 | s04_assembly.py | Local assembly of extracted reads (SPAdes) |
| 5 | s05_contig_map.py | Map contigs to host + construct (minimap2) |
| 6 | s06_junction.py | Chimeric contig / junction detection |
| 6b | s06b_junction_verify.py | Junction verification (multi-evidence scoring) |
| 6c | s06c_zygosity.py | Zygosity estimation (homo/hetero) |
| 7 | s07_host_map.py | Map all reads to host genome (bwa mem) |
| 10 | s10_copynumber.py | Copy number estimation (depth ratio) |

## Test Data

| Sample | Species | SRA | Coverage | Expected |
|--------|---------|-----|----------|----------|
| rice_G281 | *O. sativa* | SRR18236702 | ~29x | Chr3:16,439,674 (T-DNA) |
| tomato_Cas9_A2_3 | *S. lycopersicum* (Micro-Tom) | SRR13450616 | ~10x | Cas9 T-DNA + CRISPR edit |
| tomato_Cas9_A2_2 | *S. lycopersicum* (Micro-Tom) | SRR13450617 | ~5x | Cas9 T-DNA + CRISPR edit |
| tomato_Cas9_A2_1 | *S. lycopersicum* (Micro-Tom) | SRR13450618 | ~27x | CRISPR edit (Cas9 removed) |
| tomato_WT | *S. lycopersicum* (Micro-Tom) | SRR13450615 | ~22x | WT control |

## Output Structure

```
results/<sample>/
  s01_qc/            # Trimmed reads + QC reports
  s02_construct_map/  # Construct-mapped BAM + flagstat
  s03_extract/       # Extracted read pairs (+ filtered if s03b run)
  s04_assembly/      # SPAdes contigs + stats
  s05_contig_map/    # PAF alignments (host + construct)
  s06_junction/      # junctions.tsv + junction_contigs.fa
  s07_host_map/      # Host-mapped BAM (for copy number, zygosity)
```

## Key Parameters

| Parameter | Default | When to adjust |
|-----------|---------|---------------|
| --min-identity | 0.90 | Lower to 0.70 when using element_db (not specific vector) |
| --min-coverage-frac | 0.50 | Lower to 0.40 for element_db |
| --min-host-mapq | 10 | Set to 0 only for debugging |
| --min-host-aln-len | 50 | Reduce for very short contigs |

## Coverage Requirements

| Coverage | Junction Detection | Confidence |
|----------|-------------------|------------|
| >= 15x | Reliable | High |
| 10x | Possible | Medium |
| 5x | Construct presence only | Low |
| 3x | Unreliable | Very low |
