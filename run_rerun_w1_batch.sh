#!/bin/bash
#SBATCH --job-name=rg_w1_revalidate
#SBATCH --output=results/rg_w1_%A_%a.out
#SBATCH --error=results/rg_w1_%A_%a.err
#SBATCH --array=0-7
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=16:00:00
# Issue #15 / BUG-18: --mem=96G covers cucumber MaxRSS ~67G (64G default was
# the root cause of T11 _3/_4 OOM failures). rice/tomato peak ~15-40G so the
# extra reservation is safe headroom. Per-host memory requirements:
#   rice     ~10G  | tomato ~38G  | cucumber ~67G  | soybean ~4G
# If cluster scheduling pressure grows, split this into per-host arrays.

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene
SAMPLES=(
    rice_G281
    tomato_Cas9_A2_3
    tomato_Cas9_A2_2
    cucumber_line212
    cucumber_line224
    cucumber_line225
    soybean_AtYUCCA6
    soybean_UGT72E3
)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "=== $SAMPLE (task $SLURM_ARRAY_TASK_ID) ==="

# Skip if config entry missing (defensive; all 8 present at submit time)
if ! grep -q "^  ${SAMPLE}:" config.yaml; then
    echo "WARN: $SAMPLE not in config.yaml — skipping"
    exit 0
fi

# soybean 2 samples use --fanout (T8 output; UGT72E3 AC-4 path)
case "$SAMPLE" in
    soybean_*)
        python run_pipeline.py --sample "$SAMPLE" --steps 4b,5 \
            --threads 16 --no-remote-blast --fanout
        ;;
    *)
        python run_pipeline.py --sample "$SAMPLE" --steps 4b,5 \
            --threads 16 --no-remote-blast
        ;;
esac
