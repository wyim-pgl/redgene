#!/bin/bash
#SBATCH --job-name=redgene_clean
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=results/clean_rerun_%j.log

# RedGene Pipeline - Clean Re-run Script
# Runs the full pipeline on rice G281 + tomato samples from scratch.
# Steps 1-6 are fast (<30 min each). Step 7 is the bottleneck (~5h per sample).
#
# Usage:
#   sbatch run_clean.sh              # submit to SLURM
#   bash run_clean.sh                # run locally (not recommended - very long)

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene
cd /data/gpfs/assoc/pgl/develop/redgene

THREADS=16
OUTDIR=results

echo "=== Clean Re-run started at $(date) ==="
echo "Threads: $THREADS"
echo ""

# ---------------------------------------------------------------------------
# Phase 1: Rice G281 (steps 1-6, 7, 10)
# ---------------------------------------------------------------------------
echo "=== Phase 1: Rice G281 ==="
SAMPLE=rice_G281

python run_pipeline.py --sample $SAMPLE --steps 1-6 --threads $THREADS 2>&1
echo "Rice steps 1-6 complete at $(date)"

python run_pipeline.py --sample $SAMPLE --steps 7 --threads $THREADS 2>&1
echo "Rice step 7 complete at $(date)"

python run_pipeline.py --sample $SAMPLE --steps 10 --threads $THREADS 2>&1
echo "Rice step 10 complete at $(date)"

# ---------------------------------------------------------------------------
# Phase 2: Tomato WT (steps 1-7)
# ---------------------------------------------------------------------------
echo ""
echo "=== Phase 2: Tomato WT ==="
SAMPLE=tomato_WT

python run_pipeline.py --sample $SAMPLE --steps 1-6 --threads $THREADS 2>&1
echo "WT steps 1-6 complete at $(date)"

python run_pipeline.py --sample $SAMPLE --steps 7 --threads $THREADS 2>&1
echo "WT step 7 complete at $(date)"

# ---------------------------------------------------------------------------
# Phase 3: Tomato Cas9 samples (steps 1-6, 7, 8, 10)
# ---------------------------------------------------------------------------
for SAMPLE in tomato_Cas9_A2_1 tomato_Cas9_A2_2 tomato_Cas9_A2_3; do
    echo ""
    echo "=== Phase 3: $SAMPLE ==="

    python run_pipeline.py --sample $SAMPLE --steps 1-6 --threads $THREADS 2>&1
    echo "$SAMPLE steps 1-6 complete at $(date)"

    python run_pipeline.py --sample $SAMPLE --steps 7 --threads $THREADS 2>&1
    echo "$SAMPLE step 7 complete at $(date)"

    python run_pipeline.py --sample $SAMPLE --steps 8,10 --threads $THREADS 2>&1
    echo "$SAMPLE steps 8+10 complete at $(date)"
done

# ---------------------------------------------------------------------------
# Phase 4: Visualizations
# ---------------------------------------------------------------------------
echo ""
echo "=== Phase 4: Visualizations ==="

# Junction gene context plots
python scripts/viz/plot_junction_gene.py \
    --junctions results/rice_G281/s06_junction/junctions.tsv \
    --gff db/Osativa_323_v7.0.gene_exons.gff3 \
    --contigs results/rice_G281/s04_assembly/contigs.fasta \
    --sample-name rice_G281 --outdir results 2>&1

for SAMPLE in tomato_Cas9_A2_1 tomato_Cas9_A2_2 tomato_Cas9_A2_3; do
    # Editing profile (CRISPResso2-style)
    python scripts/viz/plot_editing_profile.py \
        --treatment-bam results/$SAMPLE/s07_host_map/${SAMPLE}_host.bam \
        --wt-bam results/tomato_WT/s07_host_map/tomato_WT_host.bam \
        --host-ref db/SLM_r2.0.pmol.fasta \
        --grna-targets results/$SAMPLE/s08_indel/grna_targets.tsv \
        --editing-sites results/$SAMPLE/s08_indel/editing_sites.tsv \
        --sample-name $SAMPLE --outdir results 2>&1

    # Variant effect annotation (easyGWAS-style)
    python scripts/viz/plot_editing_effects.py \
        --editing-sites results/$SAMPLE/s08_indel/editing_sites.tsv \
        --gff db/SLM_r2.0.gff3.gz \
        --host-ref db/SLM_r2.0.pmol.fasta \
        --sample-name $SAMPLE --outdir results 2>&1
done

# Junction gene context for tomato (only A2_3 has junctions)
python scripts/viz/plot_junction_gene.py \
    --junctions results/tomato_Cas9_A2_3/s06_junction/junctions.tsv \
    --gff db/SLM_r2.0.gff3.gz \
    --contigs results/tomato_Cas9_A2_3/s04_assembly/contigs.fasta \
    --sample-name tomato_Cas9_A2_3 --outdir results 2>&1

echo "Visualizations complete at $(date)"

# ---------------------------------------------------------------------------
# Phase 5: MultiQC Report
# ---------------------------------------------------------------------------
echo ""
echo "=== Phase 5: MultiQC Report ==="
python scripts/s11_multiqc.py --outdir $OUTDIR --title "RedGene Pipeline Report" 2>&1
echo "MultiQC report complete at $(date)"

echo ""
echo "=== Clean Re-run COMPLETE at $(date) ==="
echo "Check results in: $OUTDIR/"
