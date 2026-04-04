#!/usr/bin/env python3
"""Step 4: Local assembly of construct-hitting reads using SPAdes.

Assembles the extracted read pairs (from Step 3) into contigs using SPAdes
in --careful mode with multiple k-mer sizes. These are small assemblies
(~2000-6000 read pairs) that finish in under 1 minute.

Output: {outdir}/{sample}/s04_assembly/contigs.fasta
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def log(msg: str) -> None:
    """Print message to stderr."""
    print(f"[s04_assembly] {msg}", file=sys.stderr, flush=True)


def compute_n50(lengths: list[int]) -> int:
    """Compute N50 from a list of contig lengths."""
    if not lengths:
        return 0
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    cumulative = 0
    for length in sorted_lengths:
        cumulative += length
        if cumulative >= total / 2:
            return length
    return 0


def count_contigs(fasta_path: Path) -> tuple[int, int, list[int]]:
    """Count contigs, total bases, and collect lengths from a FASTA file.

    Returns:
        (num_contigs, total_bases, list_of_lengths)
    """
    num_contigs = 0
    total_bases = 0
    lengths: list[int] = []
    current_len = 0

    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current_len > 0:
                    lengths.append(current_len)
                    total_bases += current_len
                num_contigs += 1
                current_len = 0
            else:
                current_len += len(line)
        # Last contig
        if current_len > 0:
            lengths.append(current_len)
            total_bases += current_len

    return num_contigs, total_bases, lengths


def run_spades(
    r1: Path,
    r2: Path,
    outdir: Path,
    threads: int,
    sample_name: str,
) -> Path:
    """Run SPAdes assembly and return path to output contigs.fasta."""
    step_dir = outdir / sample_name / "s04_assembly"
    step_dir.mkdir(parents=True, exist_ok=True)

    spades_workdir = step_dir / "spades_work"
    output_contigs = step_dir / "contigs.fasta"

    # Validate inputs
    if not r1.exists():
        log(f"ERROR: R1 file not found: {r1}")
        sys.exit(1)
    if not r2.exists():
        log(f"ERROR: R2 file not found: {r2}")
        sys.exit(1)

    log(f"Sample: {sample_name}")
    log(f"R1: {r1}")
    log(f"R2: {r2}")
    log(f"Threads: {threads}")
    log(f"Output directory: {step_dir}")

    # Count input reads to adapt k-mer sizes
    import subprocess as _sp
    count_result = _sp.run(
        ["seqkit", "stats", "-T", str(r1)],
        capture_output=True, text=True,
    )
    read_count = 0
    try:
        lines = count_result.stdout.strip().split("\n")
        if len(lines) >= 2:
            read_count = int(lines[1].split("\t")[3].replace(",", ""))
    except (IndexError, ValueError):
        pass
    log(f"Input reads: {read_count:,} pairs")

    # Adapt k-mer sizes based on read count
    if read_count < 500:
        kmers = "21,33"
        log("Very low read count - using k-mers: 21,33")
    elif read_count < 2000:
        kmers = "21,33,55"
        log("Low read count - using k-mers: 21,33,55")
    else:
        kmers = "21,33,55,77"

    # Build SPAdes command
    cmd = [
        "spades.py",
        "--careful",
        "-1", str(r1),
        "-2", str(r2),
        "-o", str(spades_workdir),
        "-t", str(threads),
        "-k", kmers,
    ]

    log(f"Running: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        log(f"ERROR: SPAdes failed with return code {e.returncode}")
        sys.exit(1)
    except FileNotFoundError:
        log("ERROR: spades.py not found in PATH. Is SPAdes installed?")
        sys.exit(1)

    # Check for SPAdes output
    spades_contigs = spades_workdir / "contigs.fasta"
    if not spades_contigs.exists():
        # SPAdes may produce scaffolds.fasta but no contigs.fasta if
        # there were too few reads. Check for scaffolds as fallback.
        spades_scaffolds = spades_workdir / "scaffolds.fasta"
        if spades_scaffolds.exists():
            log("WARNING: contigs.fasta not found, using scaffolds.fasta")
            spades_contigs = spades_scaffolds
        else:
            log("ERROR: SPAdes produced no contigs.fasta or scaffolds.fasta")
            sys.exit(1)

    # Copy contigs to step output directory
    shutil.copy2(spades_contigs, output_contigs)
    log(f"Contigs copied to: {output_contigs}")

    # Report assembly statistics
    num_contigs, total_bases, lengths = count_contigs(output_contigs)
    n50 = compute_n50(lengths)
    max_len = max(lengths) if lengths else 0
    min_len = min(lengths) if lengths else 0

    log(f"Assembly statistics:")
    log(f"  Number of contigs: {num_contigs}")
    log(f"  Total bases:       {total_bases}")
    log(f"  N50:               {n50}")
    log(f"  Longest contig:    {max_len}")
    log(f"  Shortest contig:   {min_len}")

    # Write stats to file
    stats_file = step_dir / "assembly_stats.txt"
    with open(stats_file, "w") as fh:
        fh.write(f"sample\t{sample_name}\n")
        fh.write(f"num_contigs\t{num_contigs}\n")
        fh.write(f"total_bases\t{total_bases}\n")
        fh.write(f"n50\t{n50}\n")
        fh.write(f"max_length\t{max_len}\n")
        fh.write(f"min_length\t{min_len}\n")

    log(f"Stats written to: {stats_file}")
    return output_contigs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 4: Local assembly of extracted reads with SPAdes",
    )
    parser.add_argument(
        "--r1", required=True, type=Path,
        help="Path to extracted R1 FASTQ (from Step 3)",
    )
    parser.add_argument(
        "--r2", required=True, type=Path,
        help="Path to extracted R2 FASTQ (from Step 3)",
    )
    parser.add_argument(
        "--outdir", required=True, type=Path,
        help="Base output directory (results/)",
    )
    parser.add_argument(
        "--threads", type=int, default=8,
        help="Number of threads for SPAdes (default: 8)",
    )
    parser.add_argument(
        "--sample-name", required=True,
        help="Sample name for output organization",
    )
    args = parser.parse_args()

    run_spades(
        r1=args.r1,
        r2=args.r2,
        outdir=args.outdir,
        threads=args.threads,
        sample_name=args.sample_name,
    )

    log("Done.")


if __name__ == "__main__":
    main()
