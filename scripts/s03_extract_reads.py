#!/usr/bin/env python3
"""Step 3: Extract construct-hitting reads and their mates.

CRITICAL STEP -- extracts all reads that hit the construct reference
(mapped + supplementary/chimeric) plus their mate pairs, then outputs
them as paired FASTQ for downstream local assembly.

Strategy (from DESIGN.md):
  1. Get names of all mapped reads (samtools view -F 4 | cut -f1)
  2. Get names of supplementary alignments (samtools view -f 2048 | cut -f1)
  3. Merge and deduplicate names
  4. Name-sort BAM
  5. Extract read pairs by name list -> paired FASTQ

This guarantees mate extraction regardless of read name format
(/1 /2 vs space-delimited).
"""

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract construct-hitting reads + mates from BAM"
    )
    parser.add_argument("--bam", type=Path, required=True,
                        help="Construct-mapped BAM from s02_construct_map")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="Base output directory")
    parser.add_argument("--sample-name", type=str, required=True,
                        help="Sample name for output file naming")
    parser.add_argument("--threads", type=int, default=8,
                        help="Threads for samtools sort (default: 8)")
    return parser.parse_args()


def extract_reads(bam: Path, outdir: Path, sample_name: str, threads: int = 8) -> None:
    """Extract construct-hitting reads and their mates."""
    # Create output directory
    step_dir = outdir / sample_name / "s03_extract"
    step_dir.mkdir(parents=True, exist_ok=True)

    names_file = step_dir / "names.txt"
    hit_names_file = step_dir / "hit_names.txt"
    nsort_bam = step_dir / f"{sample_name}_nsort.bam"
    out_r1 = step_dir / f"{sample_name}_construct_R1.fq.gz"
    out_r2 = step_dir / f"{sample_name}_construct_R2.fq.gz"
    singles_fq = step_dir / f"{sample_name}_singles.fq.gz"

    # Uncompressed intermediates for samtools fastq piping
    out_r1_raw = step_dir / f"{sample_name}_construct_R1.fq"
    out_r2_raw = step_dir / f"{sample_name}_construct_R2.fq"
    singles_raw = step_dir / f"{sample_name}_singles.fq"

    # Validate input
    if not bam.exists():
        print(f"ERROR: BAM file not found: {bam}", file=sys.stderr)
        sys.exit(1)

    print(f"[s03_extract] Extracting construct-hitting reads: {sample_name}",
          file=sys.stderr)
    print(f"[s03_extract] Input BAM: {bam}", file=sys.stderr)

    # Step 1: Get names of mapped reads (-F 4 = NOT unmapped)
    print(f"[s03_extract] Collecting mapped read names (-F 4)...",
          file=sys.stderr)
    with open(names_file, "w") as fh:
        # Get mapped reads
        mapped_proc = subprocess.Popen(
            ["samtools", "view", "-F", "4", str(bam)],
            stdout=subprocess.PIPE,
        )
        cut_proc = subprocess.Popen(
            ["cut", "-f1"],
            stdin=mapped_proc.stdout,
            stdout=fh,
        )
        mapped_proc.stdout.close()
        cut_proc.communicate()
        mapped_proc.wait()

        if mapped_proc.returncode != 0:
            print("ERROR: samtools view -F 4 failed", file=sys.stderr)
            sys.exit(1)

    # Step 2: Append supplementary alignment names (-f 2048)
    print(f"[s03_extract] Collecting supplementary read names (-f 2048)...",
          file=sys.stderr)
    with open(names_file, "a") as fh:
        supp_proc = subprocess.Popen(
            ["samtools", "view", "-f", "2048", str(bam)],
            stdout=subprocess.PIPE,
        )
        cut_proc = subprocess.Popen(
            ["cut", "-f1"],
            stdin=supp_proc.stdout,
            stdout=fh,
        )
        supp_proc.stdout.close()
        cut_proc.communicate()
        supp_proc.wait()

        if supp_proc.returncode != 0:
            print("ERROR: samtools view -f 2048 failed", file=sys.stderr)
            sys.exit(1)

    # Step 3: Sort and deduplicate names
    print(f"[s03_extract] Deduplicating read names...", file=sys.stderr)
    with open(hit_names_file, "w") as fh:
        sort_proc = subprocess.Popen(
            ["sort", "-u", str(names_file)],
            stdout=fh,
        )
        sort_proc.communicate()
        if sort_proc.returncode != 0:
            print("ERROR: sort -u failed", file=sys.stderr)
            sys.exit(1)

    # Count unique hit names
    wc_result = subprocess.run(
        ["wc", "-l", str(hit_names_file)],
        capture_output=True, text=True, check=True,
    )
    hit_count = int(wc_result.stdout.strip().split()[0])
    print(f"[s03_extract] Unique read names hitting construct: {hit_count:,}",
          file=sys.stderr)

    if hit_count == 0:
        print("ERROR: No reads mapped to construct. Check your construct "
              "reference and input BAM.", file=sys.stderr)
        sys.exit(1)

    # Step 4+5: Filter by name list, then name-sort the SUBSET, then extract FASTQ
    # This is much faster than name-sorting the entire BAM
    print(f"[s03_extract] Filtering + name-sorting + extracting reads...",
          file=sys.stderr)

    # First, extract only the matching reads into a small BAM, then name-sort that
    filtered_bam = step_dir / f"{sample_name}_filtered.bam"

    # samtools view -N names -b bam | samtools sort -n -@ threads > nsort.bam
    view_proc = subprocess.Popen(
        ["samtools", "view", "-N", str(hit_names_file), "-b", str(bam)],
        stdout=subprocess.PIPE,
    )
    sort_proc = subprocess.Popen(
        ["samtools", "sort", "-n", "-@", str(threads), "-o", str(nsort_bam)],
        stdin=view_proc.stdout,
    )
    view_proc.stdout.close()
    sort_proc.communicate()
    view_proc.wait()

    if view_proc.returncode != 0 or sort_proc.returncode != 0:
        print(f"ERROR: Filter+sort failed", file=sys.stderr)
        sys.exit(1)

    # Extract FASTQ from the small name-sorted BAM
    print(f"[s03_extract] Converting to FASTQ...", file=sys.stderr)
    fastq_result = subprocess.run(
        ["samtools", "fastq",
         "-1", str(out_r1_raw),
         "-2", str(out_r2_raw),
         "-s", str(singles_raw),
         "-n",
         str(nsort_bam)],
        capture_output=True,
    )

    if fastq_result.returncode != 0:
        print(f"ERROR: samtools fastq failed (exit {fastq_result.returncode})",
              file=sys.stderr)
        print(fastq_result.stderr.decode(), file=sys.stderr)
        sys.exit(1)

    # Gzip the FASTQ outputs
    print(f"[s03_extract] Compressing FASTQ outputs...", file=sys.stderr)
    for raw_fq, gz_fq in [(out_r1_raw, out_r1), (out_r2_raw, out_r2),
                           (singles_raw, singles_fq)]:
        if raw_fq.exists() and raw_fq.stat().st_size > 0:
            subprocess.run(["gzip", "-f", str(raw_fq)], check=True)
        elif raw_fq.exists():
            # Empty file -- gzip it anyway for consistency
            subprocess.run(["gzip", "-f", str(raw_fq)], check=True)

    # Count extracted read pairs
    r1_count = 0
    if out_r1.exists():
        count_result = subprocess.run(
            ["zcat", str(out_r1)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        # Each read is 4 lines in FASTQ
        line_count = count_result.stdout.count(b"\n")
        r1_count = line_count // 4

    r2_count = 0
    if out_r2.exists():
        count_result = subprocess.run(
            ["zcat", str(out_r2)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        line_count = count_result.stdout.count(b"\n")
        r2_count = line_count // 4

    singles_count = 0
    if singles_fq.exists():
        count_result = subprocess.run(
            ["zcat", str(singles_fq)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        line_count = count_result.stdout.count(b"\n")
        singles_count = line_count // 4

    print(f"[s03_extract] === Extraction Summary for {sample_name} ===",
          file=sys.stderr)
    print(f"[s03_extract]   Unique read names: {hit_count:,}", file=sys.stderr)
    print(f"[s03_extract]   R1 reads extracted: {r1_count:,}", file=sys.stderr)
    print(f"[s03_extract]   R2 reads extracted: {r2_count:,}", file=sys.stderr)
    print(f"[s03_extract]   Read pairs: {min(r1_count, r2_count):,}",
          file=sys.stderr)
    print(f"[s03_extract]   Singleton reads: {singles_count:,}",
          file=sys.stderr)

    # Clean up intermediate files
    for tmp in [names_file, nsort_bam]:
        if tmp.exists():
            tmp.unlink()
    print(f"[s03_extract] Cleaned up intermediate files", file=sys.stderr)

    # Verify critical outputs
    if not out_r1.exists():
        print(f"ERROR: R1 output not created: {out_r1}", file=sys.stderr)
        sys.exit(1)
    if not out_r2.exists():
        print(f"ERROR: R2 output not created: {out_r2}", file=sys.stderr)
        sys.exit(1)

    print(f"[s03_extract] Done. Outputs in {step_dir}", file=sys.stderr)
    print(f"[s03_extract]   R1: {out_r1}", file=sys.stderr)
    print(f"[s03_extract]   R2: {out_r2}", file=sys.stderr)


def main() -> None:
    args = parse_args()
    extract_reads(args.bam, args.outdir, args.sample_name, args.threads)


if __name__ == "__main__":
    main()
