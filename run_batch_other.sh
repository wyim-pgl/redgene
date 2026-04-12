#!/bin/bash
#SBATCH --job-name=rg_other
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=results/rg_other_%j.out
#SBATCH --error=results/rg_other_%j.err

eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

echo "=== Corn ND207 (steps 3-5) ==="
python run_pipeline.py --sample corn_ND207 --steps 3-5 --threads 16 --no-remote-blast

echo "=== Soybean AtYUCCA6 (steps 3-5) ==="
python run_pipeline.py --sample soybean_AtYUCCA6 --steps 3-5 --threads 16 --no-remote-blast

echo "=== Soybean UGT72E3 (steps 3-5) ==="
python run_pipeline.py --sample soybean_UGT72E3 --steps 3-5 --threads 16 --no-remote-blast

echo "=== Cucumber line212 (steps 3-5) ==="
python run_pipeline.py --sample cucumber_line212 --steps 3-5 --threads 16 --no-remote-blast

echo "=== Cucumber line224 (steps 3-5) ==="
python run_pipeline.py --sample cucumber_line224 --steps 3-5 --threads 16 --no-remote-blast

echo "=== Cucumber line225 (steps 3-5) ==="
python run_pipeline.py --sample cucumber_line225 --steps 3-5 --threads 16 --no-remote-blast
