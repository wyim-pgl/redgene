#!/bin/bash
#SBATCH -J tomato_dl
#SBATCH -A cpu-s1-pgl-0
#SBATCH -p cpu-s1-pgl-0
#SBATCH -c 8
#SBATCH --mem=32g
#SBATCH -t 8:00:00
#SBATCH -o /data/gpfs/assoc/pgl/develop/redgene/test_data/tomato/download_%j.out
#SBATCH -e /data/gpfs/assoc/pgl/develop/redgene/test_data/tomato/download_%j.err

eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene/test_data/tomato

# Download and convert all 4 samples
for SRR in SRR13450616 SRR13450617 SRR13450615 SRR13450618; do
    echo "=== Downloading $SRR ==="
    
    # Check if already downloaded
    if [ -f "${SRR}_1.fastq.gz" ] && [ -f "${SRR}_2.fastq.gz" ]; then
        echo "$SRR already downloaded, skipping"
        continue
    fi
    
    # Prefetch first (faster than direct fasterq-dump)
    prefetch $SRR --max-size 50G -O . 2>&1
    
    # Convert to FASTQ
    fasterq-dump $SRR/${SRR}.sra --split-3 --threads 8 -O . 2>&1 || \
    fasterq-dump $SRR --split-3 --threads 8 -O . 2>&1
    
    # Compress
    if [ -f "${SRR}_1.fastq" ]; then
        pigz -p 8 ${SRR}_1.fastq ${SRR}_2.fastq
        echo "$SRR: FASTQ compressed"
    fi
    
    # Clean up SRA cache
    rm -rf $SRR
    
    echo "$SRR done"
done

echo "=== All downloads complete ==="
ls -lh *.fastq.gz 2>/dev/null
