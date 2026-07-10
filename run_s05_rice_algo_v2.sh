#!/bin/bash
#SBATCH --job-name=rg_rice_s05_v2
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=logs/rg_rice_s05_v2_%j.out
#SBATCH --error=logs/rg_rice_s05_v2_%j.err

set -euo pipefail

# The login shell exports SBATCH_PARTITION / SBATCH_ACCOUNT / SBATCH_TIMELIMIT,
# and SLURM precedence is CLI > env > script, so those env vars silently
# override the #SBATCH directives above. Unset them before any child sbatch.
# See the lab wiki, guide/pronghorn.md ("SBATCH_PARTITION environment variable
# propagation").
unset SBATCH_PARTITION SBATCH_ACCOUNT SBATCH_TIMELIMIT SLURM_ACCOUNT || true

cd /data/gpfs/assoc/pgl/develop/redgene
export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"

echo "host: $(hostname)"
echo "branch: $(git rev-parse --abbrev-ref HEAD) @ $(git rev-parse --short HEAD)"
echo "start: $(date -Is)"

python run_pipeline.py --sample rice_G281 --steps 5 --threads 16 --no-remote-blast \
  --outdir-override results/rice_G281_algo_v2 \
  --host-bam-override results/rice_G281/s04_host_map/rice_G281_host.bam

echo "end: $(date -Is)"
