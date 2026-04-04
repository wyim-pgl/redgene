#!/bin/bash
#SBATCH -J redgene_s07
#SBATCH -A cpu-s1-pgl-0
#SBATCH -p cpu-s1-pgl-0
#SBATCH -c 16
#SBATCH --mem=64g
#SBATCH -t 4:00:00
#SBATCH -o results/rice_G281/s07_host_map/slurm_%j.out
#SBATCH -e results/rice_G281/s07_host_map/slurm_%j.err

eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene

python scripts/s07_host_map.py \
  --r1 results/rice_G281/s01_qc/rice_G281_R1.fq.gz \
  --r2 results/rice_G281/s01_qc/rice_G281_R2.fq.gz \
  --host-ref db/Osativa_323_v7.0.fa \
  --outdir results \
  --threads 16 \
  --sample-name rice_G281
