#!/bin/bash
# Issue #2 [AC-7] - coverage sensitivity sweep (3 hosts x 4 coverage = 12 tasks).
#
# Operator submits manually (no auto-submit):
#
#     sbatch run_coverage_sensitivity.sh
#
# The array index selects one (sample, coverage) pair from the SAMPLES/COVS
# tables below. Each task:
#   1. subsamples the master FASTQ pair with scripts/util/subsample_reads.py
#   2. lays out results/<sample>_cov<tag>/ with symlinks
#   3. runs the core pipeline (steps 1-5) via run_pipeline.py
#      (needs config.yaml entry <sample>_cov<tag> - added in sibling commit)
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

# 3 hosts x 4 coverage = 12 tasks. Keys MUST match config.yaml samples.
SAMPLES=(rice_G281 tomato_Cas9_A2_3 cucumber_line225)
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

# --- Step 1: map master FASTQ pair per sample (from config.yaml) ----------
# Case-based mapping avoids a yq dependency. Keep in sync with config.yaml.
case "$SAMPLE" in
    rice_G281)
        MASTER_R1="test_data/rice_G281_R1.fastq.gz"
        MASTER_R2="test_data/rice_G281_R2.fastq.gz"
        ;;
    tomato_Cas9_A2_3)
        MASTER_R1="test_data/tomato/SRR13450616_1.fastq.gz"
        MASTER_R2="test_data/tomato/SRR13450616_2.fastq.gz"
        ;;
    cucumber_line225)
        MASTER_R1="data/cucumber/SRR12082195_1.fastq.gz"
        MASTER_R2="data/cucumber/SRR12082195_2.fastq.gz"
        ;;
    *)
        echo "[cov-sensitivity] ERROR: unknown sample $SAMPLE" >&2
        exit 2
        ;;
esac

if [[ ! -f "$MASTER_R1" || ! -f "$MASTER_R2" ]]; then
    echo "[cov-sensitivity] ERROR: master FASTQ not found for $SAMPLE:" >&2
    echo "  R1=$MASTER_R1" >&2
    echo "  R2=$MASTER_R2" >&2
    exit 3
fi

COV_SAMPLE="${SAMPLE}_cov${COV}"
OUT_DIR="results/${COV_SAMPLE}"
OUT_PREFIX="${OUT_DIR}/${COV_SAMPLE}"

mkdir -p "$OUT_DIR"

# --- Step 2: subsample master FASTQ ---------------------------------------
python scripts/util/subsample_reads.py \
    --r1 "$MASTER_R1" --r2 "$MASTER_R2" \
    --fraction "$FRAC" --seed 42 \
    --output-prefix "$OUT_PREFIX"

echo "[cov-sensitivity] subsample done -> ${OUT_PREFIX}_R[12].fq.gz"

# --- Step 3: run the core pipeline (steps 1-5) ----------------------------
# run_pipeline.py uses the default config.yaml and resolves the sample key
# <SAMPLE>_cov<COV> which must already exist in config.yaml (see commit
# "Issue #2 [AC-7] add 12 coverage-sensitivity sample entries to config.yaml").
python run_pipeline.py \
    --sample "$COV_SAMPLE" \
    --steps 1-5 \
    --threads 16 \
    --no-remote-blast

echo "[cov-sensitivity] task=$IDX done (sample=$COV_SAMPLE)"
