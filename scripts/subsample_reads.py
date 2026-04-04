#!/usr/bin/env python3
"""Subsample paired-end FASTQ files to a target coverage depth.

Calculates the number of reads needed based on genome size and read length,
then uses seqkit sample to randomly subsample.

Usage:
  python subsample_reads.py --r1 in_R1.fq.gz --r2 in_R2.fq.gz \
      --genome-size 374000000 --target-coverage 10 \
      --read-length 100 --output-prefix out_10x
"""

import argparse
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Subsample reads to target coverage")
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    parser.add_argument("--genome-size", type=int, required=True,
                        help="Genome size in bp")
    parser.add_argument("--target-coverage", type=float, required=True,
                        help="Target coverage depth (e.g. 10)")
    parser.add_argument("--read-length", type=int, default=100,
                        help="Average read length (default: 100)")
    parser.add_argument("--total-pairs", type=int, default=0,
                        help="Total pairs in input (0=auto-count)")
    parser.add_argument("--output-prefix", type=str, required=True,
                        help="Output prefix (creates PREFIX_R1.fq.gz, PREFIX_R2.fq.gz)")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    # Calculate required pairs
    bases_needed = args.genome_size * args.target_coverage
    pairs_needed = int(bases_needed / (2 * args.read_length))

    # Get total pairs if not provided
    total_pairs = args.total_pairs
    if total_pairs == 0:
        print(f"[subsample] Counting reads in {args.r1}...", file=sys.stderr)
        result = subprocess.run(
            ["seqkit", "stats", "-T", str(args.r1)],
            capture_output=True, text=True,
        )
        # Parse seqkit stats output
        lines = result.stdout.strip().split("\n")
        if len(lines) >= 2:
            total_pairs = int(lines[1].split("\t")[3].replace(",", ""))
        else:
            print("ERROR: Could not count reads", file=sys.stderr)
            sys.exit(1)

    fraction = min(1.0, pairs_needed / total_pairs)

    print(f"[subsample] Genome size: {args.genome_size:,} bp", file=sys.stderr)
    print(f"[subsample] Target coverage: {args.target_coverage}x", file=sys.stderr)
    print(f"[subsample] Pairs needed: {pairs_needed:,}", file=sys.stderr)
    print(f"[subsample] Total pairs: {total_pairs:,}", file=sys.stderr)
    print(f"[subsample] Fraction: {fraction:.4f}", file=sys.stderr)

    if fraction >= 1.0:
        print(f"[subsample] WARNING: Not enough reads for {args.target_coverage}x. "
              f"Using all reads ({total_pairs * 2 * args.read_length / args.genome_size:.1f}x).",
              file=sys.stderr)
        # Just symlink
        out_r1 = Path(f"{args.output_prefix}_R1.fq.gz")
        out_r2 = Path(f"{args.output_prefix}_R2.fq.gz")
        out_r1.symlink_to(args.r1.resolve())
        out_r2.symlink_to(args.r2.resolve())
        return

    out_r1 = f"{args.output_prefix}_R1.fq.gz"
    out_r2 = f"{args.output_prefix}_R2.fq.gz"

    # Subsample R1
    print(f"[subsample] Subsampling R1...", file=sys.stderr)
    subprocess.run([
        "seqkit", "sample", "-p", str(fraction), "-s", str(args.seed),
        "-o", out_r1, str(args.r1),
    ], check=True)

    # Subsample R2 with same seed
    print(f"[subsample] Subsampling R2...", file=sys.stderr)
    subprocess.run([
        "seqkit", "sample", "-p", str(fraction), "-s", str(args.seed),
        "-o", out_r2, str(args.r2),
    ], check=True)

    # Verify
    result = subprocess.run(
        ["seqkit", "stats", "-T", out_r1],
        capture_output=True, text=True,
    )
    lines = result.stdout.strip().split("\n")
    if len(lines) >= 2:
        actual_reads = int(lines[1].split("\t")[3].replace(",", ""))
        actual_cov = actual_reads * 2 * args.read_length / args.genome_size
        print(f"[subsample] Actual reads: {actual_reads:,} pairs", file=sys.stderr)
        print(f"[subsample] Actual coverage: {actual_cov:.1f}x", file=sys.stderr)

    print(f"[subsample] Output: {out_r1}, {out_r2}", file=sys.stderr)


if __name__ == "__main__":
    main()
