#!/bin/bash
# soybean_UGT72E3 timed out at 16 h in the array 5803578 run: 106 transgene-
# positive sites assembled sequentially, each recruiting ~1.27M candidate read
# pairs (Gmax is 1.1 Gbp with heavy repeat content). Phase 1+1.5 DID finish and
# wrote positive_sites.{json,pkl}, so this re-enters at Phase 2_3 via the per-
# site array fan-out — one SLURM task per site, up to 25 concurrent — and a
# Phase 4 annotate+report as an afterok dependency. Wall time drops from
# 106 x sequential to ~ceil(106/25) x per-site.
#
# All paths are lifted verbatim from the timed-out run's CMD line
# (logs/rg_s05_hard_5803578_10.err) so the fan-out sees byte-identical inputs.

set -euo pipefail

unset SBATCH_PARTITION SBATCH_ACCOUNT SBATCH_TIMELIMIT SLURM_ACCOUNT || true

cd /data/gpfs/assoc/pgl/develop/redgene
export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"

REPO=/data/gpfs/assoc/pgl/develop/redgene
STEPDIR="$REPO/results/hardened_20260710/soybean_UGT72E3/s05_insert_assembly"

# Paths the s05 CLI needs (submit_s05_array.sh reads these from the environment).
export S05_HOST_BAM="$REPO/results/soybean_UGT72E3/s04_host_map/soybean_UGT72E3_host.bam"
export S05_HOST_REF="$REPO/db/Gmax_v4.0.fa"
export S05_ELEMENT_DB="$REPO/element_db/gmo_combined_db.fa"
export S05_CONSTRUCT_REF="$REPO/element_db/gmo_combined_db.fa"
export S05_EXTRA_ELEMENT_DB="$REPO/results/hardened_20260710/soybean_UGT72E3/s04b_construct_asm/contigs.fasta"
export S05_COMMON_PAYLOAD_DB="$REPO/element_db/common_payload.fa"
export S05_S03_R1="$REPO/results/hardened_20260710/soybean_UGT72E3/s03_extract/soybean_UGT72E3_construct_R1.fq.gz"
export S05_S03_R2="$REPO/results/hardened_20260710/soybean_UGT72E3/s03_extract/soybean_UGT72E3_construct_R2.fq.gz"
export S05_NO_REMOTE_BLAST=1

# NB: submit_s05_array.sh does NOT forward --mask-bed. soybean_gmax_v4.bed is a
# 0-byte file (its mask rationale rows are all *-EMPTY), so --mask-bed was a
# no-op in the sequential run too — dropping it changes nothing.

# Run on the same partition/account as the rest of this work, not the script's
# cpu-s2 defaults.
export SLURM_PARTITION=cpu-s1-pgl-0
export SLURM_ACCOUNT=cpu-s1-pgl-0

# Per-task resources: 8 CPU, 64G (candidate-read load is ~1.27M pairs/site),
# 4 h ceiling per site (the sequential run averaged ~9 min/site but a few are
# far heavier). 25 concurrent by default.
export S05_ARRAY_MEM=64G
export S05_ARRAY_TIME=4:00:00
export S05_ARRAY_THROTTLE=25

bash scripts/submit_s05_array.sh soybean_UGT72E3 "$STEPDIR" 8 \
    "$REPO/results/hardened_20260710"
