#!/bin/bash
#SBATCH --job-name=rg_s05_hard
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=16:00:00
#SBATCH --array=0-10%4
#SBATCH --output=logs/rg_s05_hard_%A_%a.out
#SBATCH --error=logs/rg_s05_hard_%A_%a.err

# Re-run step 5 for every production sample at the current HEAD of
# algo-improvements-2026-07. Every s05 report on disk predates this branch:
# the April 2026 reports predate all seven algorithm fixes, and even
# results/rice_G281_algo_v2 (SLURM 5799803) finished 4 minutes before the
# per-insert border-count fix and 12 minutes before the element-DB border
# correction landed. Output goes to a fresh, dated root so nothing is
# overwritten.
#
# 96G, not the usual 64G: cucumber s05 OOMs at 64G (batch 5628138 needed 96G).
# One allocation for all samples keeps the array homogeneous.
#
# Prereq: bash stage_hardened_rerun.sh

set -euo pipefail

# The login shell exports SBATCH_PARTITION / SBATCH_ACCOUNT / SBATCH_TIMELIMIT,
# and SLURM precedence is CLI > env > script, so those env vars silently
# override the #SBATCH directives above (observed on job 5799803, which got a
# 13-day time limit instead of 16 h). See the lab wiki, guide/pronghorn.md.
unset SBATCH_PARTITION SBATCH_ACCOUNT SBATCH_TIMELIMIT SLURM_ACCOUNT || true

cd /data/gpfs/assoc/pgl/develop/redgene
export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"

SAMPLES=(
    rice_G281
    tomato_WT
    tomato_Cas9_A2_1
    tomato_Cas9_A2_2
    tomato_Cas9_A2_3
    corn_ND207
    cucumber_line212
    cucumber_line224
    cucumber_line225
    soybean_AtYUCCA6
    soybean_UGT72E3
)
S="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
OUTROOT=results/hardened_20260710

echo "host:    $(hostname)"
echo "sample:  $S  (task $SLURM_ARRAY_TASK_ID)"
echo "branch:  $(git rev-parse --abbrev-ref HEAD) @ $(git rev-parse --short HEAD)"
echo "start:   $(date -Is)"

python run_pipeline.py --sample "$S" --steps 5 --threads 16 --no-remote-blast \
    --outdir-override "$OUTROOT" \
    --host-bam-override "results/$S/s04_host_map/${S}_host.bam"

echo "end:     $(date -Is)"
