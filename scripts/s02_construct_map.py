#!/usr/bin/env python3
"""Step 2: Map trimmed reads to construct reference with bwa mem.

Indexes the construct reference (if needed), maps paired-end reads,
sorts and indexes the BAM, and reports mapping statistics.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map reads to construct reference with bwa mem"
    )
    parser.add_argument("--r1", type=Path, required=True,
                        help="Trimmed forward reads (from s01_qc)")
    parser.add_argument("--r2", type=Path, required=True,
                        help="Trimmed reverse reads (from s01_qc)")
    parser.add_argument("--construct-ref", type=Path, required=True,
                        help="Construct reference FASTA")
    default_univec = Path(__file__).resolve().parent.parent / "db" / "univec_plant_vectors.fa"
    parser.add_argument("--univec", type=Path, default=default_univec,
                        help="UniVec plant vectors FASTA (default: db/univec_plant_vectors.fa, "
                             "use --no-univec to disable)"
    )
    parser.add_argument("--no-univec", action="store_true",
                        help="Disable UniVec plant vector inclusion")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="Base output directory")
    parser.add_argument("--threads", type=int, default=8,
                        help="Number of threads (default: 8)")
    parser.add_argument("--sample-name", type=str, required=True,
                        help="Sample name for output file naming")
    return parser.parse_args()


def index_reference(ref: Path) -> None:
    """Index construct reference with bwa index if not already done."""
    bwt_file = ref.parent / f"{ref.name}.bwt"
    if bwt_file.exists():
        print(f"[s02_construct_map] BWA index exists: {bwt_file}",
              file=sys.stderr)
        return

    print(f"[s02_construct_map] Indexing reference: {ref}", file=sys.stderr)
    subprocess.run(["bwa", "index", str(ref)], check=True)

    # Also create samtools faidx if missing
    fai_file = ref.parent / f"{ref.name}.fai"
    if not fai_file.exists():
        print(f"[s02_construct_map] Creating samtools faidx: {ref}",
              file=sys.stderr)
        subprocess.run(["samtools", "faidx", str(ref)], check=True)


def _build_combined_ref(construct_ref: Path, univec: Path | None,
                        step_dir: Path) -> Path:
    """Combine construct ref + UniVec plant vectors into a single FASTA.

    Returns the path to use for BWA mapping (original if no UniVec,
    combined file if UniVec provided).
    """
    if univec is None or not univec.exists():
        return construct_ref

    combined = step_dir / "construct_plus_univec.fa"
    print(f"[s02_construct_map] Combining construct ref + UniVec plant vectors",
          file=sys.stderr)
    construct_text = construct_ref.read_text()
    univec_text = univec.read_text()
    with open(combined, "w") as out:
        out.write(construct_text)
        if not construct_text.endswith("\n"):
            out.write("\n")
        out.write(univec_text)

    n_construct = sum(1 for line in construct_text.splitlines()
                      if line.startswith(">"))
    n_univec = sum(1 for line in univec_text.splitlines()
                   if line.startswith(">"))
    print(f"[s02_construct_map]   Construct: {n_construct} seqs, "
          f"UniVec plant: {n_univec} seqs, "
          f"Combined: {n_construct + n_univec} seqs", file=sys.stderr)
    return combined


def run_mapping(r1: Path, r2: Path, construct_ref: Path, outdir: Path,
                threads: int, sample_name: str,
                univec: Path | None = None) -> None:
    """Map reads to construct, sort, index, and report stats."""
    # Create output directory
    step_dir = outdir / sample_name / "s02_construct_map"
    step_dir.mkdir(parents=True, exist_ok=True)

    bam_file = step_dir / f"{sample_name}_construct.bam"
    bai_file = step_dir / f"{sample_name}_construct.bam.bai"
    flagstat_file = step_dir / f"{sample_name}_construct.flagstat"
    depth_file = step_dir / f"{sample_name}_construct_depth.txt"

    # Validate inputs
    for label, path in [("R1", r1), ("R2", r2), ("Construct ref", construct_ref)]:
        if not path.exists():
            print(f"ERROR: {label} not found: {path}", file=sys.stderr)
            sys.exit(1)

    # Build combined reference (construct + UniVec plant vectors)
    mapping_ref = _build_combined_ref(construct_ref, univec, step_dir)

    # Index reference if needed
    index_reference(mapping_ref)

    # Map with bwa mem, pipe to samtools sort
    print(f"[s02_construct_map] Mapping reads to construct: {sample_name}",
          file=sys.stderr)
    print(f"[s02_construct_map] R1: {r1}", file=sys.stderr)
    print(f"[s02_construct_map] R2: {r2}", file=sys.stderr)
    print(f"[s02_construct_map] Reference: {mapping_ref}", file=sys.stderr)

    rg_tag = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"

    # bwa mem | samtools sort
    bwa_cmd = [
        "bwa", "mem",
        "-t", str(threads),
        "-R", rg_tag,
        str(mapping_ref),
        str(r1),
        str(r2),
    ]
    sort_cmd = [
        "samtools", "sort",
        "-@", str(threads),
        "-o", str(bam_file),
        "-",
    ]

    # Redirect BWA stderr to file to avoid pipe buffer deadlock on large datasets
    bwa_stderr_file = step_dir / f"{sample_name}_bwa.log"
    print(f"[s02_construct_map] Running bwa mem | samtools sort",
          file=sys.stderr)
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
        print(f"ERROR: bwa mem failed (exit {bwa_proc.returncode})",
              file=sys.stderr)
        print(bwa_stderr[-2000:], file=sys.stderr)
        sys.exit(1)
    if sort_proc.returncode != 0:
        print(f"ERROR: samtools sort failed (exit {sort_proc.returncode})",
              file=sys.stderr)
        print(sort_stderr.decode(), file=sys.stderr)
        sys.exit(1)

    # Index BAM
    print(f"[s02_construct_map] Indexing BAM", file=sys.stderr)
    subprocess.run(["samtools", "index", str(bam_file)], check=True)

    # Generate flagstat
    print(f"[s02_construct_map] Generating flagstat", file=sys.stderr)
    with open(flagstat_file, "w") as fh:
        subprocess.run(["samtools", "flagstat", str(bam_file)],
                       stdout=fh, check=True)

    # Print flagstat to stderr
    flagstat_text = flagstat_file.read_text()
    print(f"[s02_construct_map] === Flagstat for {sample_name} ===",
          file=sys.stderr)
    for line in flagstat_text.strip().splitlines():
        print(f"[s02_construct_map]   {line}", file=sys.stderr)

    # Calculate construct coverage depth
    print(f"[s02_construct_map] Calculating depth", file=sys.stderr)
    depth_result = subprocess.run(
        ["samtools", "depth", "-a", str(bam_file)],
        capture_output=True, text=True, check=True,
    )

    # Write depth file
    with open(depth_file, "w") as fh:
        fh.write(depth_result.stdout)

    # Parse depth to compute summary stats
    depths: list[int] = []
    for line in depth_result.stdout.strip().splitlines():
        if line:
            parts = line.split("\t")
            if len(parts) >= 3:
                depths.append(int(parts[2]))

    if depths:
        total_positions = len(depths)
        covered_positions = sum(1 for d in depths if d > 0)
        mean_depth = sum(depths) / total_positions
        sorted_depths = sorted(depths)
        median_depth = sorted_depths[total_positions // 2]

        print(f"[s02_construct_map] === Construct Coverage ===",
              file=sys.stderr)
        print(f"[s02_construct_map]   Total positions: {total_positions:,}",
              file=sys.stderr)
        print(f"[s02_construct_map]   Covered positions (>0x): "
              f"{covered_positions:,} "
              f"({100 * covered_positions / total_positions:.1f}%)",
              file=sys.stderr)
        print(f"[s02_construct_map]   Mean depth: {mean_depth:.1f}x",
              file=sys.stderr)
        print(f"[s02_construct_map]   Median depth: {median_depth}x",
              file=sys.stderr)
        print(f"[s02_construct_map]   Min depth: {sorted_depths[0]}x",
              file=sys.stderr)
        print(f"[s02_construct_map]   Max depth: {sorted_depths[-1]}x",
              file=sys.stderr)
    else:
        print(f"[s02_construct_map] WARNING: No depth data (no reads mapped?)",
              file=sys.stderr)

    # Extract mapped read count from flagstat
    match = re.search(r"(\d+) \+ \d+ mapped", flagstat_text)
    if match:
        mapped_reads = int(match.group(1))
        print(f"[s02_construct_map]   Mapped reads: {mapped_reads:,}",
              file=sys.stderr)

    # Verify outputs
    for f in [bam_file, bai_file, flagstat_file]:
        if not f.exists():
            print(f"ERROR: Expected output not created: {f}", file=sys.stderr)
            sys.exit(1)

    print(f"[s02_construct_map] Done. BAM: {bam_file}", file=sys.stderr)


def main() -> None:
    args = parse_args()
    univec = None if args.no_univec else args.univec
    run_mapping(args.r1, args.r2, args.construct_ref, args.outdir,
                args.threads, args.sample_name, univec=univec)


if __name__ == "__main__":
    main()
