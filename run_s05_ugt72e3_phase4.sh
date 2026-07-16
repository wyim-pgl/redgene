#!/bin/bash
#SBATCH --job-name=ugt_phase4
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/ugt_phase4_%j.out
#SBATCH --error=logs/ugt_phase4_%j.err

# Phase 4 (annotate + FP filters + report + s05_stats.txt) for soybean_UGT72E3.
# The fan-out array 5811471 completed all 106 sites, but its afterok Phase 4 job
# (5811472) died at exit 127: submit_s05_array.sh's --wrap uses `micromamba
# activate`, and micromamba is not resolvable in that job's shell on this
# cluster. This re-runs Phase 4 with the env bin prepended directly, the pattern
# CLAUDE.md documents for non-interactive shells.

set -euo pipefail
unset SBATCH_PARTITION SBATCH_ACCOUNT SBATCH_TIMELIMIT SLURM_ACCOUNT || true

cd /data/gpfs/assoc/pgl/develop/redgene
export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"

REPO=/data/gpfs/assoc/pgl/develop/redgene
OUTROOT="$REPO/results/hardened_20260710"

echo "host:   $(hostname)"
echo "branch: $(git rev-parse --abbrev-ref HEAD) @ $(git rev-parse --short HEAD)"
echo "start:  $(date -Is)"

python scripts/s05_insert_assembly.py --phase 4 --threads 8 \
    --host-bam "$REPO/results/soybean_UGT72E3/s04_host_map/soybean_UGT72E3_host.bam" \
    --host-ref "$REPO/db/Gmax_v4.0.fa" \
    --element-db "$REPO/element_db/gmo_combined_db.fa" \
    --construct-ref "$REPO/element_db/gmo_combined_db.fa" \
    --extra-element-db "$OUTROOT/soybean_UGT72E3/s04b_construct_asm/contigs.fasta" \
    --common-payload-db "$REPO/element_db/common_payload.fa" \
    --s03-r1 "$OUTROOT/soybean_UGT72E3/s03_extract/soybean_UGT72E3_construct_R1.fq.gz" \
    --s03-r2 "$OUTROOT/soybean_UGT72E3/s03_extract/soybean_UGT72E3_construct_R2.fq.gz" \
    --outdir "$OUTROOT" \
    --sample-name soybean_UGT72E3 \
    --config "$REPO/config.yaml" \
    --no-remote-blast

echo "end:    $(date -Is)"
