#!/bin/bash
# Issue #1 [AC-2] — Remote BLAST SLURM template for CANDIDATE FP verification.
#
# This script is a TEMPLATE. It is NOT to be executed or submitted automatically
# — see docs/measurements/ac2_cuc225_fp_workflow.md for the operator workflow.
# The SBATCH headers are intentionally declared so `sbatch prepare_remote_blast.sh`
# works once the operator fills in SAMPLE / QUERY / OUT_TSV variables below.
#
# Why SLURM? `blastn -remote` can take minutes-to-hours per MB of query and must
# run from a node allowed to reach NCBI. Login-node submission is forbidden.
#
# Why max_target_seqs=5? We only need "closest host ortholog vs. vector sequence"
# evidence; deeper hit tables inflate runtime without improving FP/TP calls.
#
#SBATCH --job-name=rg_remote_blast
#SBATCH --output=results/rg_remote_blast_%j.out
#SBATCH --error=results/rg_remote_blast_%j.err
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

# ---------------------------------------------------------------------------
# Operator: set these three variables, then `sbatch scripts/util/prepare_remote_blast.sh`
# ---------------------------------------------------------------------------
: "${SAMPLE:?must set SAMPLE (e.g. cucumber_line225)}"
: "${QUERY:?must set QUERY (e.g. cucumber_line225_cand_inserts.fa from extract_cand_for_blast.py)}"
OUT_TSV="${OUT_TSV:-${SAMPLE}_cand_vs_nt.tsv}"

if [[ ! -s "$QUERY" ]]; then
    echo "ERROR: query FASTA $QUERY is empty or missing" >&2
    exit 2
fi

# 6-column outfmt: qseqid, sseqid, pident, length, evalue, stitle
# (stitle is essential for naming the FP taxon in the verdict TSV)
OUTFMT="6 qseqid sseqid pident length evalue stitle"

echo "[remote-blast] sample=$SAMPLE query=$QUERY out=$OUT_TSV"
echo "[remote-blast] NCBI policy: keep max_target_seqs small; avoid parallel submits"

blastn \
    -remote \
    -db nt \
    -query "$QUERY" \
    -out "$OUT_TSV" \
    -outfmt "$OUTFMT" \
    -max_target_seqs 5 \
    -evalue 1e-5

echo "[remote-blast] done → $OUT_TSV"
