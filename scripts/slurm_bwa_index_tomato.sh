#!/bin/bash
#SBATCH -J bwa_idx_mt
#SBATCH -A cpu-s1-pgl-0
#SBATCH -p cpu-s1-pgl-0
#SBATCH -c 4
#SBATCH --mem=16g
#SBATCH -t 2:00:00
#SBATCH -o /data/gpfs/assoc/pgl/develop/redgene/db/bwa_index_tomato_%j.out
#SBATCH -e /data/gpfs/assoc/pgl/develop/redgene/db/bwa_index_tomato_%j.err

eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene

echo "=== Building BWA index for Micro-Tom reference ==="
bwa index db/SLM_r2.0.pmol.fasta
echo "=== Done ==="
