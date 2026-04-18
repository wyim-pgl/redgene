#!/bin/bash
#SBATCH --job-name=rg_s04_mm2_poc
#SBATCH --output=results/rg_s04_mm2_poc_%j.out
#SBATCH --error=results/rg_s04_mm2_poc_%j.err
#
# T12: rice_G281 s04 minimap2 -ax sr PoC (conditional approval experiment).
# Per docs/team-review/work_implementation_plan.md Task 12 and team-consensus.md
# §5 + §8 Unresolved Item #1. Target: 4/4 PASS -> v1.1 migration approved.
# Isolates output in results/rice_G281_mm2_poc/ to avoid touching the canonical
# BWA baseline in results/rice_G281/ (which T11 batch also consumes).

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene

SAMPLE=rice_G281
HOST=db/Osativa_323_v7.0.fa

# Extract R1/R2 paths from config.yaml. yq not installed on Pronghorn, so use
# an awk fallback scoped to the rice_G281 block.
R1=$(awk '
    /^  rice_G281:[[:space:]]*$/ { in_block = 1; next }
    in_block && /^  [a-z]/       { in_block = 0 }
    in_block && /r1:/            { print $2; exit }
' config.yaml)
R2=$(awk '
    /^  rice_G281:[[:space:]]*$/ { in_block = 1; next }
    in_block && /^  [a-z]/       { in_block = 0 }
    in_block && /r2:/            { print $2; exit }
' config.yaml)

if [ -z "${R1:-}" ] || [ -z "${R2:-}" ]; then
    echo "ERROR: could not read R1/R2 from config.yaml rice_G281 block" >&2
    exit 1
fi
echo "[T12] R1=$R1"
echo "[T12] R2=$R2"
[ -f "$R1" ] || { echo "ERROR: R1 missing: $R1" >&2; exit 1; }
[ -f "$R2" ] || { echo "ERROR: R2 missing: $R2" >&2; exit 1; }

OUT_ROOT=results/rice_G281_mm2_poc
OUT_SAMPLE=${OUT_ROOT}/${SAMPLE}
OUT_S04=${OUT_SAMPLE}/s04_host_map
mkdir -p "$OUT_S04"
mkdir -p "$OUT_SAMPLE/s03_extract"

# Symlink s03 outputs from the BWA baseline so s04b/s05 have their inputs
# without re-running 1-3. If baseline is absent, bail with a clear message.
BWA_S03=results/${SAMPLE}/s03_extract
if [ ! -f "${BWA_S03}/${SAMPLE}_construct_R1.fq.gz" ]; then
    echo "ERROR: BWA baseline s03 outputs missing at ${BWA_S03}." >&2
    echo "       Run the standard pipeline for rice_G281 first." >&2
    exit 1
fi
for f in "${SAMPLE}_construct_R1.fq.gz" "${SAMPLE}_construct_R2.fq.gz" \
         "${SAMPLE}_singles.fq.gz" hit_names.txt; do
    src="${BWA_S03}/${f}"
    dst="${OUT_SAMPLE}/s03_extract/${f}"
    if [ -f "$src" ] && [ ! -e "$dst" ]; then
        ln -sf "$(readlink -f "$src")" "$dst"
    fi
done

# minimap2 index cached next to the FASTA; built once per host ref (~2-3 min
# for rice 374 Mbp). Reuse across future PoC runs.
IDX="${HOST}.mmi"
if [ ! -f "$IDX" ]; then
    echo "[T12] minimap2 -d $IDX $HOST (one-time, ~2-3 min)"
    minimap2 -d "$IDX" "$HOST"
fi

# Map with timing. -ax sr = short-read preset (paired-end).
echo "[T12] minimap2 -ax sr mapping..."
START=$(date +%s)
/usr/bin/time -f "%E elapsed, %M maxRSS KB" \
    minimap2 -ax sr -t 16 "$IDX" "$R1" "$R2" 2> "$OUT_S04/mm2.log" \
    | samtools sort -@ 8 -o "$OUT_S04/${SAMPLE}_host.bam" -
samtools index "$OUT_S04/${SAMPLE}_host.bam"
END=$(date +%s)
ELAPSED_S04=$((END - START))
echo "[T12] s04 minimap2 wall time: ${ELAPSED_S04}s"
echo "$ELAPSED_S04" > "$OUT_S04/wall_time_seconds.txt"

# Run 4b + 5 via run_pipeline.py with the minimap2 BAM swapped in and an
# isolated outdir. --no-remote-blast keeps the PoC local-only (matches how
# won_yim wants the comparison: identical s05 config, only s04 engine swapped).
echo "[T12] running 4b + 5 with minimap2 BAM..."
python run_pipeline.py --sample ${SAMPLE} --steps 4b,5 \
    --threads 8 --no-remote-blast \
    --host-bam-override "$(readlink -f "$OUT_S04/${SAMPLE}_host.bam")" \
    --outdir-override "$OUT_ROOT" 2>&1 | tee "$OUT_ROOT/pipeline.log"

echo "[T12] PoC complete. Post-analysis:"
echo "      bash docs/measurements/s04_minimap2_poc_analyze.sh"
