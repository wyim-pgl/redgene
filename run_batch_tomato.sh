#!/bin/bash
#SBATCH --job-name=rg_tomato
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=results/rg_tomato_%j.out
#SBATCH --error=results/rg_tomato_%j.err

eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

# WT needs steps 3-4 only (s01/s02 done, no insert assembly for WT)
echo "=== Tomato WT (steps 3-4) ==="
python run_pipeline.py --sample tomato_WT --steps 3-4 --threads 16

echo "=== Tomato Cas9 A2_1 (steps 3-5) ==="
python run_pipeline.py --sample tomato_Cas9_A2_1 --steps 3-5 --threads 16 --no-remote-blast

echo "=== Tomato Cas9 A2_2 (steps 3-5) ==="
python run_pipeline.py --sample tomato_Cas9_A2_2 --steps 3-5 --threads 16 --no-remote-blast

echo "=== Tomato Cas9 A2_3 (steps 3-5) ==="
python run_pipeline.py --sample tomato_Cas9_A2_3 --steps 3-5 --threads 16 --no-remote-blast
