#!/bin/bash
# Post-W1-completion analysis script.
# Run after SLURM job <JOBID> completes.
# Usage: bash docs/superpowers/runs/2026-04-24-w1-post-analysis.sh <JOBID>

set -euo pipefail
JOBID="${1:-$(cat docs/superpowers/runs/2026-04-24-w1-jobid.txt)}"
OUT_DIR="docs/superpowers/runs"
POSTFIX="${OUT_DIR}/2026-04-24-post-w1.txt"

echo "=== SLURM state ==="
sacct -j "$JOBID" --format=JobID,JobName,State,Elapsed,MaxRSS -n | head

echo ""
echo "=== Verdict summary per sample ==="
: > "$POSTFIX"
for s in rice_G281 tomato_Cas9_A2_3 tomato_Cas9_A2_2 cucumber_line212 \
         cucumber_line224 cucumber_line225 soybean_AtYUCCA6 soybean_UGT72E3; do
    echo "=== $s ===" | tee -a "$POSTFIX"
    if ls results/${s}/s05_insert_assembly/insertion_*_report.txt > /dev/null 2>&1; then
        grep -h "^Verdict:" results/${s}/s05_insert_assembly/insertion_*_report.txt 2>/dev/null \
            | awk -F' —' '{print $1}' | sort | uniq -c | tee -a "$POSTFIX"
    else
        echo "(no s05 results)" | tee -a "$POSTFIX"
    fi
done

echo ""
echo "=== AC-1 GT anchor check ==="
declare -A GT=(
    [rice_G281]="Chr3_16439674"
    [tomato_Cas9_A2_3]="SLM_r2.0ch01_9100"
    [tomato_Cas9_A2_2]="SLM_r2.0ch08_6510"
    [cucumber_line212]="LKUO03001392_2751687"
    [cucumber_line224]="LKUO03001512_581328"
    [cucumber_line225]="LKUO03001451_6501"
)
ANCHOR_PASS=0
for s in "${!GT[@]}"; do
    prefix=${GT[$s]}
    f=$(ls results/${s}/s05_insert_assembly/insertion_${prefix}*_report.txt 2>/dev/null | head -1)
    if [ -n "$f" ]; then
        v=$(grep -h "^Verdict:" "$f" | head -1)
        echo "$s $prefix -> $v"
        if [[ "$v" == *"CANDIDATE"* ]]; then ANCHOR_PASS=$((ANCHOR_PASS+1)); fi
    else
        echo "$s $prefix -> (no report)"
    fi
done
echo "AC-1 anchor recall: $ANCHOR_PASS / ${#GT[@]}"

echo ""
echo "=== AC-2 FP budget (target <=5 per sample) ==="
for s in rice_G281 tomato_Cas9_A2_3 cucumber_line225 soybean_AtYUCCA6 soybean_UGT72E3; do
    n=$(grep -l "^Verdict: CANDIDATE" results/${s}/s05_insert_assembly/insertion_*_report.txt 2>/dev/null | wc -l)
    echo "$s CAND total: $n"
done

echo ""
echo "=== AC-4 wall time ==="
sacct -j "$JOBID" --format=JobName,Elapsed,MaxRSS -n

echo ""
echo "=== AC-6 audit header check ==="
for s in rice_G281 soybean_AtYUCCA6; do
    if [ -f "results/$s/audit_header.json" ]; then
        n_fields=$(jq -r '[.input_sha256, .pipeline_commit, .db_manifest, .software_versions] | length' "results/$s/audit_header.json" 2>/dev/null)
        echo "$s audit_header.json: $n_fields/4 fields"
    fi
done

echo ""
echo "=== diff vs pre-W1 baseline ==="
diff docs/superpowers/runs/2026-04-24-pre-w1-baseline.txt "$POSTFIX" || true
