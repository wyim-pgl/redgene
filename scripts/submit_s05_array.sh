#!/bin/bash
# scripts/submit_s05_array.sh
#
# T8 임시안 (v1.0 MVP): per-site SLURM array fan-out for s05 Phase 2+3,
# then Phase 4 annotate+report as an afterok dependency.
#
# Prereq: run_pipeline.py --fanout (or a manual invocation) must have
# produced <step_dir>/positive_sites.json and positive_sites.pkl via
# `python scripts/s05_insert_assembly.py ... --phase 1_1.5`.
#
# Usage:
#   scripts/submit_s05_array.sh <sample> <step_dir> [threads] [outdir]
#
# Arguments:
#   sample    sample name (e.g. soybean_UGT72E3)
#   step_dir  results/<sample>/s05_insert_assembly path
#   threads   cpus-per-task for each Phase 2+3 task (default 8)
#   outdir    --outdir to pass to s05 (default: dirname(dirname(step_dir)))
#
# v1.1 replaces this wrapper with the full module split in scripts/s05/.

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "usage: $0 <sample> <step_dir> [threads] [outdir]" >&2
    exit 2
fi

SAMPLE="$1"
STEPDIR="$2"
THREADS="${3:-8}"
# Derive outdir from step_dir: <outdir>/<sample>/s05_insert_assembly
if [ "$#" -ge 4 ]; then
    OUTDIR="$4"
else
    OUTDIR="$(dirname "$(dirname "$STEPDIR")")"
fi

SITES_JSON="$STEPDIR/positive_sites.json"
if [ ! -f "$SITES_JSON" ]; then
    echo "ERROR: $SITES_JSON not found — run Phase 1+1.5 first" >&2
    echo "       (e.g. python scripts/s05_insert_assembly.py ... --phase 1_1.5)" >&2
    exit 1
fi

# Load config paths the s05 CLI needs (host_bam, host_ref, element_db, etc.)
# The caller (run_pipeline.py) knows these; submit_s05_array.sh can cheat and
# ask run_pipeline.py to dry-print the s05 command.  But for simplicity we
# require the invoker to set a few env vars.  Path resolution happens inside
# run_pipeline.py when --fanout is used, so these env vars come from there.
: "${S05_HOST_BAM:?S05_HOST_BAM must be set (exported by run_pipeline.py --fanout)}"
: "${S05_HOST_REF:?S05_HOST_REF must be set}"
: "${S05_ELEMENT_DB:?S05_ELEMENT_DB must be set}"
# Optional args — collect into a string passed through to the array tasks.
# Values are single-quoted so paths with whitespace/special chars survive
# the `sbatch --wrap="..."` heredoc intact (M-7, same fix-class as I-2).
EXTRA_ARGS=""
[ -n "${S05_CONSTRUCT_REF:-}" ] && EXTRA_ARGS+=" --construct-ref '${S05_CONSTRUCT_REF}'"
[ -n "${S05_EXTRA_ELEMENT_DB:-}" ] && EXTRA_ARGS+=" --extra-element-db '${S05_EXTRA_ELEMENT_DB}'"
[ -n "${S05_COMMON_PAYLOAD_DB:-}" ] && EXTRA_ARGS+=" --common-payload-db '${S05_COMMON_PAYLOAD_DB}'"
[ -n "${S05_S03_R1:-}" ] && EXTRA_ARGS+=" --s03-r1 '${S05_S03_R1}'"
[ -n "${S05_S03_R2:-}" ] && EXTRA_ARGS+=" --s03-r2 '${S05_S03_R2}'"
if [ "${S05_NO_REMOTE_BLAST:-0}" = "1" ]; then
    EXTRA_ARGS+=" --no-remote-blast"
fi

SLURM_PARTITION="${SLURM_PARTITION:-cpu-s2-core-0}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-cpu-s2-pgl-0}"

# Count positive sites.
N=$(python -c "import json, sys; print(len(json.load(open(sys.argv[1]))))" "$SITES_JSON")
echo "[T8] $SAMPLE: $N positive sites in $SITES_JSON"

mkdir -p "${OUTDIR}/slurm_logs"

# Build the shared CLI args once — they're identical across array tasks
# modulo the --site-id that each task picks.  All paths are single-quoted
# so the `sbatch --wrap="..."` heredoc below passes them through a single
# extra shell layer without word-splitting on whitespace/glob chars (I-2).
COMMON_ARGS="--host-bam '$S05_HOST_BAM' --host-ref '$S05_HOST_REF'"
COMMON_ARGS+=" --element-db '$S05_ELEMENT_DB'"
COMMON_ARGS+=" --outdir '$OUTDIR' --sample-name '$SAMPLE'"
COMMON_ARGS+="${EXTRA_ARGS}"

if [ "$N" -eq 0 ]; then
    echo "[T8] WARNING: 0 positive sites for $SAMPLE; skipping array, running Phase 4 only" >&2
    PHASE4_JOBID=$(sbatch --parsable \
        --partition="$SLURM_PARTITION" --account="$SLURM_ACCOUNT" \
        --time=1:00:00 --mem=16G --cpus-per-task=4 \
        --job-name="s05_${SAMPLE}_phase4" \
        --output="${OUTDIR}/slurm_logs/s05_phase4_${SAMPLE}_%j.out" \
        --wrap="set -euo pipefail; \
eval \"\$(micromamba shell hook --shell bash)\"; \
micromamba activate redgene; \
python scripts/s05_insert_assembly.py --phase 4 --threads 4 $COMMON_ARGS")
    echo "[T8] Phase 4 job (no array): $PHASE4_JOBID"
    exit 0
fi

# Build the SITES array — one site_id per line via mapfile so whitespace
# or special chars in a site_id can't break bash word-splitting (I-3).
mapfile -t SITES < <(python -c "
import json, sys
with open(sys.argv[1]) as f:
    for s in json.load(f):
        print(s['site_id'])
" "$SITES_JSON")
# Shell-safe re-materialization for the sbatch heredoc below: each element
# is individually `printf %q`-escaped so the array survives a second parse
# without word-splitting (I-3).
SITES_QUOTED=$(printf '%q ' "${SITES[@]}")
# Space-separated form for the human log line only; never re-parse.
SITES_STR="${SITES[*]}"

# Array job — one task per positive site, max 25 concurrent.
ARRAY_SCRIPT="${STEPDIR}/s05_array_${SAMPLE}.sbatch"
cat > "$ARRAY_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --job-name=s05_${SAMPLE}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --array=0-$((N-1))%25
#SBATCH --output=${OUTDIR}/slurm_logs/s05_array_${SAMPLE}_%A_%a.out
#SBATCH --error=${OUTDIR}/slurm_logs/s05_array_${SAMPLE}_%A_%a.err
set -euo pipefail
eval "\$(micromamba shell hook --shell bash)"
micromamba activate redgene

SITES=(${SITES_QUOTED})
SITE="\${SITES[\$SLURM_ARRAY_TASK_ID]}"
echo "=== [T8] site \$SITE (task \$SLURM_ARRAY_TASK_ID / $((N-1))) ==="

python scripts/s05_insert_assembly.py \\
    --phase 2_3 --site-id "\$SITE" \\
    --threads ${THREADS} \\
    ${COMMON_ARGS}
EOF
chmod +x "$ARRAY_SCRIPT"

ARRAY_JOBID=$(sbatch --parsable "$ARRAY_SCRIPT")
echo "[T8] Array job: $ARRAY_JOBID  (0-$((N-1)), ${SITES_STR})"

# Phase 4 as afterok dependency on the whole array.
PHASE4_JOBID=$(sbatch --parsable \
    --partition="$SLURM_PARTITION" --account="$SLURM_ACCOUNT" \
    --time=1:00:00 --mem=16G --cpus-per-task=4 \
    --dependency=afterok:${ARRAY_JOBID} \
    --job-name="s05_${SAMPLE}_phase4" \
    --output="${OUTDIR}/slurm_logs/s05_phase4_${SAMPLE}_%j.out" \
    --wrap="set -euo pipefail; \
eval \"\$(micromamba shell hook --shell bash)\"; \
micromamba activate redgene; \
python scripts/s05_insert_assembly.py --phase 4 --threads 4 $COMMON_ARGS")
echo "[T8] Phase 4 job (afterok:${ARRAY_JOBID}): $PHASE4_JOBID"

# Keep the sbatch script on disk for post-hoc audit; no cleanup on success.
echo "[T8] Array script retained at: $ARRAY_SCRIPT"
