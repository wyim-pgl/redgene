#!/usr/bin/env python3
"""Step 7: Map all trimmed reads to the host reference genome with bwa mem.

Maps paired-end reads to the host genome, sorts and indexes the BAM,
generates flagstat, and calculates genome-wide depth statistics.
These results are used downstream for copy number estimation (Step 10).
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path


def log(msg: str) -> None:
    """Print message to stderr."""
    print(f"[s07_host_map] {msg}", file=sys.stderr, flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map trimmed reads to host reference genome with bwa mem"
    )
    parser.add_argument("--r1", type=Path, required=True,
                        help="Trimmed forward reads (from s01_qc)")
    parser.add_argument("--r2", type=Path, required=True,
                        help="Trimmed reverse reads (from s01_qc)")
    parser.add_argument("--host-ref", type=Path, required=True,
                        help="Host reference genome FASTA")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="Base output directory")
    parser.add_argument("--threads", type=int, default=8,
                        help="Number of threads (default: 8)")
    parser.add_argument("--sample-name", type=str, required=True,
                        help="Sample name for output file naming")
    return parser.parse_args()


def index_reference(ref: Path) -> None:
    """Index host reference with bwa index if not already done."""
    bwt_file = ref.parent / f"{ref.name}.bwt"
    if bwt_file.exists():
        log(f"BWA index exists: {bwt_file}")
        return

    log(f"Indexing reference: {ref}")
    subprocess.run(["bwa", "index", str(ref)], check=True)

    # Also create samtools faidx if missing
    fai_file = ref.parent / f"{ref.name}.fai"
    if not fai_file.exists():
        log(f"Creating samtools faidx: {ref}")
        subprocess.run(["samtools", "faidx", str(ref)], check=True)


def run_mapping(r1: Path, r2: Path, host_ref: Path, outdir: Path,
                threads: int, sample_name: str) -> None:
    """Map reads to host genome, sort, index, and report stats."""
    # Create output directory
    step_dir = outdir / sample_name / "s07_host_map"
    step_dir.mkdir(parents=True, exist_ok=True)

    bam_file = step_dir / f"{sample_name}_host.bam"
    bai_file = step_dir / f"{sample_name}_host.bam.bai"
    flagstat_file = step_dir / f"{sample_name}_host.flagstat"
    depth_stats_file = step_dir / f"{sample_name}_host_depth_stats.txt"

    # Validate inputs
    for label, path in [("R1", r1), ("R2", r2), ("Host ref", host_ref)]:
        if not path.exists():
            log(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # Index reference if needed
    index_reference(host_ref)

    # Map with bwa mem, pipe to samtools sort
    log(f"Mapping reads to host genome: {sample_name}")
    log(f"R1: {r1}")
    log(f"R2: {r2}")
    log(f"Reference: {host_ref}")

    rg_tag = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"

    bwa_cmd = [
        "bwa", "mem",
        "-t", str(threads),
        "-R", rg_tag,
        str(host_ref),
        str(r1),
        str(r2),
    ]
    sort_cmd = [
        "samtools", "sort",
        "-@", str(threads),
        "-o", str(bam_file),
        "-",
    ]

    # Redirect BWA stderr to a file to avoid pipe buffer deadlock
    # on large datasets (>64KB stderr output from progress messages)
    bwa_stderr_file = step_dir / f"{sample_name}_bwa.log"
    log("Running bwa mem | samtools sort")
    with open(bwa_stderr_file, "w") as bwa_err_fh:
        bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE,
                                    stderr=bwa_err_fh)
        sort_proc = subprocess.Popen(sort_cmd, stdin=bwa_proc.stdout,
                                     stderr=subprocess.PIPE)
        # Allow bwa_proc to receive SIGPIPE if sort_proc exits
        bwa_proc.stdout.close()
        sort_stdout, sort_stderr = sort_proc.communicate()
        bwa_proc.wait()

    if bwa_proc.returncode != 0:
        bwa_stderr = bwa_stderr_file.read_text()
        log(f"ERROR: bwa mem failed (exit {bwa_proc.returncode})")
        log(bwa_stderr[-2000:])  # Last 2000 chars of error log
        sys.exit(1)
    if sort_proc.returncode != 0:
        log(f"ERROR: samtools sort failed (exit {sort_proc.returncode})")
        log(sort_stderr.decode())
        sys.exit(1)

    # Index BAM
    log("Indexing BAM")
    subprocess.run(["samtools", "index", str(bam_file)], check=True)

    # Generate flagstat
    log("Generating flagstat")
    with open(flagstat_file, "w") as fh:
        subprocess.run(["samtools", "flagstat", str(bam_file)],
                       stdout=fh, check=True)

    # Print flagstat to stderr
    flagstat_text = flagstat_file.read_text()
    log(f"=== Flagstat for {sample_name} ===")
    for line in flagstat_text.strip().splitlines():
        log(f"  {line}")

    # Extract mapping rate from flagstat
    match = re.search(r"(\d+) \+ \d+ mapped \((\S+)", flagstat_text)
    if match:
        mapped_reads = int(match.group(1))
        mapping_rate = match.group(2)
        log(f"Mapped reads: {mapped_reads:,}")
        log(f"Mapping rate: {mapping_rate}")

    # Calculate genome-wide depth statistics
    # Use samtools depth on the full BAM; for large genomes, pipe and compute
    # running statistics to avoid loading everything into memory
    log("Calculating genome-wide depth statistics...")
    depth_proc = subprocess.Popen(
        ["samtools", "depth", "-a", str(bam_file)],
        stdout=subprocess.PIPE, text=True,
    )

    total_positions = 0
    covered_positions = 0
    total_depth = 0
    depths_sample: list[int] = []
    sample_every = 1  # Will be adjusted after first pass count

    # Stream depth output to compute stats without loading all into memory
    # Collect every Nth position for median estimation
    for line in depth_proc.stdout:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        d = int(parts[2])
        total_positions += 1
        total_depth += d
        if d > 0:
            covered_positions += 1
        # Sample positions for median calculation (keep up to ~1M positions)
        if total_positions <= 1_000_000 or total_positions % sample_every == 0:
            depths_sample.append(d)
            # Adjust sampling rate once we know the scale
            if len(depths_sample) == 1_000_000 and total_positions == 1_000_000:
                # Genome is large; start sub-sampling
                sample_every = 10

    depth_proc.wait()

    if total_positions > 0:
        mean_depth = total_depth / total_positions
        sorted_sample = sorted(depths_sample)
        median_depth = sorted_sample[len(sorted_sample) // 2]
        coverage_pct = 100 * covered_positions / total_positions

        log(f"=== Host Genome Depth Stats ===")
        log(f"  Total positions: {total_positions:,}")
        log(f"  Covered positions (>0x): {covered_positions:,} ({coverage_pct:.1f}%)")
        log(f"  Mean depth: {mean_depth:.1f}x")
        log(f"  Median depth: {median_depth}x")

        # Write depth stats to file
        with open(depth_stats_file, "w") as fh:
            fh.write(f"total_positions\t{total_positions}\n")
            fh.write(f"covered_positions\t{covered_positions}\n")
            fh.write(f"coverage_pct\t{coverage_pct:.2f}\n")
            fh.write(f"mean_depth\t{mean_depth:.2f}\n")
            fh.write(f"median_depth\t{median_depth}\n")
    else:
        log("WARNING: No depth data (no reads mapped?)")
        with open(depth_stats_file, "w") as fh:
            fh.write(f"total_positions\t0\n")
            fh.write(f"covered_positions\t0\n")
            fh.write(f"coverage_pct\t0.00\n")
            fh.write(f"mean_depth\t0.00\n")
            fh.write(f"median_depth\t0\n")

    # Verify outputs
    for f in [bam_file, bai_file, flagstat_file, depth_stats_file]:
        if not f.exists():
            log(f"ERROR: Expected output not created: {f}")
            sys.exit(1)

    log(f"Done. BAM: {bam_file}")


def main() -> None:
    args = parse_args()
    run_mapping(args.r1, args.r2, args.host_ref, args.outdir,
                args.threads, args.sample_name)


if __name__ == "__main__":
    main()
