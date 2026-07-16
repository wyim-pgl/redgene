#!/bin/bash
# Stage the override trees for the 2026-07-10 s05 re-run.
#
# s05 re-runs read s03 reads and the s04b per-sample construct assembly from
# {outdir}/{sample}/, and --outdir-override repoints outdir. So an override run
# sees an empty tree unless those inputs are staged first. This is what the rice
# algo_v2 run (SLURM 5799803) did by hand; this script does it for every sample.
#
# Inputs are symlinked, never copied: the re-run must consume byte-identical
# reads and contigs so the only variable is the s05 code at HEAD.
#
# Usage: bash stage_hardened_rerun.sh [DEST_ROOT]

set -euo pipefail

REPO=/data/gpfs/assoc/pgl/develop/redgene
DEST_ROOT="${1:-$REPO/results/hardened_20260710}"

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

cd "$REPO"
printf '%-20s %-10s %-10s %-8s\n' SAMPLE s03 s04b host_bam

for S in "${SAMPLES[@]}"; do
    SRC="$REPO/results/$S"
    DST="$DEST_ROOT/$S"
    mkdir -p "$DST/s03_extract" "$DST/s04_host_map"

    # s03 reads. run_pipeline prefers {S}_filtered_R1.fq.gz (s03b output) over
    # {S}_construct_R1.fq.gz, so link the whole directory and let it choose the
    # same file it chose in April.
    s03_n=0
    for f in "$SRC"/s03_extract/*; do
        [ -e "$f" ] || continue
        ln -sfn "$f" "$DST/s03_extract/$(basename "$f")"
        s03_n=$((s03_n + 1))
    done

    # s04b construct assembly. Absent for tomato_WT, tomato_Cas9_A2_1 and
    # corn_ND207 — those never ran step 4b, and their April s05 ran without
    # --extra-element-db. Leaving the directory absent reproduces that exactly.
    s04b=missing
    if [ -s "$SRC/s04b_construct_asm/contigs.fasta" ]; then
        mkdir -p "$DST/s04b_construct_asm"
        for f in "$SRC"/s04b_construct_asm/contigs*.fasta; do
            ln -sfn "$f" "$DST/s04b_construct_asm/$(basename "$f")"
        done
        s04b=linked
    fi

    # Host BAM. --host-bam-override already points s05 at the original, but the
    # rice run staged these too; keep the trees identical.
    bam=missing
    if [ -f "$SRC/s04_host_map/${S}_host.bam" ]; then
        ln -sfn "$SRC/s04_host_map/${S}_host.bam" "$DST/s04_host_map/${S}_host.bam"
        [ -f "$SRC/s04_host_map/${S}_host.bam.bai" ] &&
            ln -sfn "$SRC/s04_host_map/${S}_host.bam.bai" "$DST/s04_host_map/${S}_host.bam.bai"
        bam=linked
    fi

    printf '%-20s %-10s %-10s %-8s\n' "$S" "${s03_n} files" "$s04b" "$bam"
done

echo
echo "Staged under: $DEST_ROOT"
