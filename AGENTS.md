# AGENTS.md

## SLURM Policy

- **NEVER run heavy computational jobs on the login node.** This is a policy violation on Pronghorn HPC. Always use `sbatch` to submit jobs via SLURM.
- Heavy jobs include: BWA, bowtie2, SPAdes, minimap2, samtools sort/merge, fastp, pigz decompression, and any bioinformatics pipeline step.
- Partition: `cpu-s1-pgl-0`, Account: `cpu-s1-pgl-0`
- Light commands (ls, cat, grep, squeue, git, python scripts that just parse files) are OK on the login node.

## Memory Allocation

- Always allocate **4GB RAM per thread** when running bioinformatics tools.
- Examples: 16 threads = 64GB (`--mem=64G`), 8 threads = 32GB (`--mem=32G`), 4 threads = 16GB (`--mem=16G`).
- For SPAdes, also set `--memory` flag to match (e.g., `-m 60` for 16 threads with some headroom).
