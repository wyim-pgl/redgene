#!/bin/bash
# T12 s04 minimap2 PoC post-completion analysis.
# Verifies 4 acceptance criteria for won_yim conditional-approval experiment.
# Usage: bash docs/measurements/s04_minimap2_poc_analyze.sh [JOBID]

set -euo pipefail
JOBID="${1:-$(cat docs/measurements/s04_minimap2_poc_jobid.txt 2>/dev/null || echo unknown)}"

echo "=== SLURM state (JOBID=$JOBID) ==="
if [ "$JOBID" != "unknown" ]; then
    sacct -j "$JOBID" --format=JobID,State,Elapsed,MaxRSS -n 2>/dev/null \
        || echo "(sacct not available or job not found)"
fi

echo ""
echo "=== Criterion 1: Chr3:16,439,674 verdict preserved ==="
MM2_REPORT=$(ls results/rice_G281_mm2_poc/rice_G281/s05_insert_assembly/insertion_Chr3_16439674*_report.txt 2>/dev/null | head -1 || true)
BWA_REPORT=$(ls results/rice_G281/s05_insert_assembly/insertion_Chr3_16439674*_report.txt 2>/dev/null | head -1 || true)
if [ -n "$MM2_REPORT" ]; then
    MM2_V=$(grep "^Verdict:" "$MM2_REPORT" | head -1 || echo "(no verdict line)")
    echo "minimap2 ($MM2_REPORT): $MM2_V"
else
    echo "minimap2: (no report; check results/rice_G281_mm2_poc/rice_G281/s05_insert_assembly/)"
fi
if [ -n "$BWA_REPORT" ]; then
    BWA_V=$(grep "^Verdict:" "$BWA_REPORT" | head -1 || echo "(no verdict line)")
    echo "BWA      ($BWA_REPORT): $BWA_V"
else
    echo "BWA:      (no report; run baseline first)"
fi
echo "(PASS if both contain CANDIDATE/CONFIRMED; FAIL if mm2 < BWA)"

echo ""
echo "=== Criterion 2: Phase 1 transgene-positive site count (+/-20%) ==="
MM2_TSV=results/rice_G281_mm2_poc/rice_G281/s05_insert_assembly/site_tier_classification.tsv
BWA_TSV=results/rice_G281/s05_insert_assembly/site_tier_classification.tsv
# Issue #14 fix: 4th column `transgene_positive` is True/False, not a line prefix.
# Original grep -c "^transgene-positive" always returned 0.
_count_positive() {
    local tsv="$1"
    [ -f "$tsv" ] || { echo 0; return; }
    awk -F'\t' 'NR>1 && $4 == "True" {n++} END {print n+0}' "$tsv"
}
MM2_N=$(_count_positive "$MM2_TSV")
BWA_N=$(_count_positive "$BWA_TSV")
echo "minimap2 transgene-positive: $MM2_N"
echo "BWA      transgene-positive: $BWA_N"
python3 - <<PY
mm2, bwa = $MM2_N, $BWA_N
ratio = mm2 / bwa if bwa > 0 else 0
print(f'ratio = {ratio:.2f}  (+/-20% band: 0.80..1.20)')
print('PASS' if bwa > 0 and 0.80 <= ratio <= 1.20 else 'FAIL or insufficient data')
PY

echo ""
echo "=== Criterion 3: s04 wall time reduction (>=30% vs BWA) ==="
MM2_WALL=$(cat results/rice_G281_mm2_poc/rice_G281/s04_host_map/wall_time_seconds.txt 2>/dev/null || echo "")
# BWA baseline: from audit_header or BWA log timestamps. Approximate with
# CLAUDE.md "5-7h bottleneck" -> 5h = 18000s lower bound.
BWA_WALL=$(python3 -c "
import json, pathlib, sys
p = pathlib.Path('results/rice_G281/audit_header.json')
if p.exists():
    d = json.loads(p.read_text())
    # no wall time recorded in audit_header; fall through
print('')
")
echo "minimap2 s04 wall: ${MM2_WALL:-unknown}s"
echo "BWA      s04 wall: (not logged; CLAUDE.md reports 5-7h bottleneck)"
python3 - <<PY
mm2 = int("$MM2_WALL") if "$MM2_WALL".isdigit() else 0
bwa_est = 5 * 3600  # 5h lower bound of stated 5-7h window
save = (bwa_est - mm2) / bwa_est * 100 if mm2 > 0 else 0
print(f'estimated savings vs 5h BWA baseline: {save:.0f}%  (need >=30%)')
print('PASS (savings >= 30%)' if mm2 > 0 and save >= 30 else 'FAIL or insufficient data')
PY

echo ""
echo "=== Criterion 4: MAPQ<20 soft-clip ratio (BUG-7 guard) ==="
for label_bam in \
    "BWA:results/rice_G281/s04_host_map/rice_G281_host.bam" \
    "minimap2:results/rice_G281_mm2_poc/rice_G281/s04_host_map/rice_G281_host.bam"; do
    label="${label_bam%%:*}"
    bam="${label_bam##*:}"
    if [ -f "$bam" ]; then
        stats=$(samtools view -F 256 -f 2 "$bam" 2>/dev/null \
            | awk '{
                if ($6 ~ /S/) {
                    if ($5 < 20) low++; else high++
                }
            } END {
                low = low + 0; high = high + 0; tot = low + high;
                if (tot > 0) printf "low<20=%d high>=20=%d ratio=%.3f (n=%d)", low, high, low/tot, tot;
                else print "(no soft-clip reads)"
            }')
        echo "$label: $stats"
    else
        echo "$label: (BAM missing: $bam)"
    fi
done
echo "(PASS if minimap2 ratio <= BWA ratio + 0.05)"

echo ""
echo "=== OVERALL ==="
echo "4/4 PASS -> v1.1 minimap2 -ax sr adoption approved (rice only)"
echo "Any FAIL  -> BWA retained for v1.0+v1.1; minimap2 deferred to v2.0 grant"
echo ""
echo "Record the decision in docs/measurements/s04_minimap2_poc.md."
