#!/bin/bash
# Issue #2 [AC-7] — coverage sensitivity sweep (3 hosts x 4 coverage = 12 tasks).
#
# TEMPLATE ONLY — do NOT auto-submit. Operator must invoke `sbatch` manually:
#
#     sbatch run_coverage_sensitivity.sh
#
# The array index selects one (sample, coverage) pair from the SAMPLES/COVS
# tables below. Each task:
#   1. subsamples the master FASTQ pair with scripts/util/subsample_reads.py
#   2. lays out results/<sample>_cov<tag>/ with symlinks
#   3. runs the core pipeline (steps 1-5) via run_pipeline.py
#
#SBATCH --job-name=rg_cov_sensitivity
#SBATCH --output=results/rg_cov_%A_%a.out
#SBATCH --error=results/rg_cov_%A_%a.err
#SBATCH --array=0-11
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
# --mem=96G is the cucumber-safe ceiling (BUG-18 lesson). rice/tomato peak
# 15-40G so the extra reservation is a no-op penalty for those tasks.

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

# 3 hosts x 4 coverage = 12 tasks. Edit these tables to change the grid.
SAMPLES=(rice_G281 tomato_A2_3 cucumber_line225)
COVS=(5x 10x 15x 20x)
FRACTIONS=(0.25 0.50 0.75 1.00)   # approx 5/10/15/20x for a ~20x master set

N_SAMPLES=${#SAMPLES[@]}
N_COVS=${#COVS[@]}
IDX=${SLURM_ARRAY_TASK_ID:?must be run as SLURM array task}

SAMPLE_IDX=$(( IDX / N_COVS ))
COV_IDX=$(( IDX % N_COVS ))
SAMPLE=${SAMPLES[$SAMPLE_IDX]}
COV=${COVS[$COV_IDX]}
FRAC=${FRACTIONS[$COV_IDX]}

echo "[cov-sensitivity] task=$IDX sample=$SAMPLE cov=$COV fraction=$FRAC"

# --- Step 1: subsample master FASTQ ---------------------------------------
MASTER_R1="input/${SAMPLE}_R1.fq.gz"
MASTER_R2="input/${SAMPLE}_R2.fq.gz"
OUT_PREFIX="results/${SAMPLE}_cov${COV}/${SAMPLE}_cov${COV}"

mkdir -p "results/${SAMPLE}_cov${COV}"

python scripts/util/subsample_reads.py \
    --r1 "$MASTER_R1" --r2 "$MASTER_R2" \
    --fraction "$FRAC" --seed 42 \
    --output-prefix "$OUT_PREFIX"

# --- Step 2: run the pipeline (steps 1-5) ---------------------------------
# Note: this requires config.yaml to contain a <SAMPLE>_cov<COV> entry that
# points at the subsampled paths above. Scaffold only — config wiring is
# an explicit follow-up task (see docs/measurements/coverage_sensitivity.md).
echo "[cov-sensitivity] (scaffold) would now run: python run_pipeline.py --sample ${SAMPLE}_cov${COV} --steps 1-5 --threads 16"

echo "[cov-sensitivity] task=$IDX done"
