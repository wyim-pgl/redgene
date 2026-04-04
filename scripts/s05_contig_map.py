#!/usr/bin/env python3
"""Step 5: Map assembled contigs to host and construct references.

Uses minimap2 -x asm5 to align assembled contigs (from Step 4) against
both the host genome reference and the construct reference. The PAF output
files are used by Step 6 to identify chimeric contigs spanning
host-transgene junctions.

Output:
    {outdir}/{sample}/s05_contig_map/{sample}_contigs_to_host.paf
    {outdir}/{sample}/s05_contig_map/{sample}_contigs_to_construct.paf
"""

import argparse
import subprocess
import sys
from pathlib import Path


def log(msg: str) -> None:
    """Print message to stderr."""
    print(f"[s05_contig_map] {msg}", file=sys.stderr, flush=True)


def parse_paf_stats(paf_path: Path) -> dict[str, int]:
    """Parse a PAF file and return basic mapping statistics."""
    total_alignments = 0
    mapped_contigs: set[str] = set()
    total_aligned_bases = 0
    targets: set[str] = set()

    if not paf_path.exists() or paf_path.stat().st_size == 0:
        return {
            "total_alignments": 0,
            "mapped_contigs": 0,
            "total_aligned_bases": 0,
            "target_sequences": 0,
        }

    with open(paf_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 12:
                continue
            total_alignments += 1
            qname = fields[0]
            tname = fields[5]
            aligned_bases = int(fields[9])  # Number of matching bases
            mapped_contigs.add(qname)
            targets.add(tname)
            total_aligned_bases += aligned_bases

    return {
        "total_alignments": total_alignments,
        "mapped_contigs": len(mapped_contigs),
        "total_aligned_bases": total_aligned_bases,
        "target_sequences": len(targets),
    }


def run_minimap2(
    contigs: Path,
    reference: Path,
    output_paf: Path,
    threads: int,
    label: str,
) -> None:
    """Run minimap2 -x asm5 and write PAF output."""
    cmd = [
        "minimap2",
        "-x", "asm5",
        "-t", str(threads),
        "--secondary=yes",
        "-o", str(output_paf),
        str(reference),
        str(contigs),
    ]

    log(f"Mapping contigs to {label}: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        log(f"ERROR: minimap2 ({label}) failed with return code {e.returncode}")
        sys.exit(1)
    except FileNotFoundError:
        log("ERROR: minimap2 not found in PATH. Is minimap2 installed?")
        sys.exit(1)

    if not output_paf.exists():
        log(f"ERROR: minimap2 did not produce output file: {output_paf}")
        sys.exit(1)

    # Report stats
    stats = parse_paf_stats(output_paf)
    log(f"  {label} mapping statistics:")
    log(f"    Total alignments:    {stats['total_alignments']}")
    log(f"    Mapped contigs:      {stats['mapped_contigs']}")
    log(f"    Total aligned bases: {stats['total_aligned_bases']}")
    log(f"    Target sequences:    {stats['target_sequences']}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 5: Map assembled contigs to host and construct references",
    )
    parser.add_argument(
        "--contigs", required=True, type=Path,
        help="Path to assembled contigs FASTA (from Step 4)",
    )
    parser.add_argument(
        "--host-ref", required=True, type=Path,
        help="Path to host genome reference FASTA",
    )
    parser.add_argument(
        "--construct-ref", required=True, type=Path,
        help="Path to construct reference FASTA",
    )
    parser.add_argument(
        "--outdir", required=True, type=Path,
        help="Base output directory (results/)",
    )
    parser.add_argument(
        "--threads", type=int, default=8,
        help="Number of threads for minimap2 (default: 8)",
    )
    parser.add_argument(
        "--sample-name", required=True,
        help="Sample name for output organization",
    )
    args = parser.parse_args()

    # Validate inputs
    for label, path in [
        ("contigs", args.contigs),
        ("host reference", args.host_ref),
        ("construct reference", args.construct_ref),
    ]:
        if not path.exists():
            log(f"ERROR: {label} file not found: {path}")
            sys.exit(1)

    step_dir = args.outdir / args.sample_name / "s05_contig_map"
    step_dir.mkdir(parents=True, exist_ok=True)

    host_paf = step_dir / f"{args.sample_name}_contigs_to_host.paf"
    construct_paf = step_dir / f"{args.sample_name}_contigs_to_construct.paf"

    log(f"Sample: {args.sample_name}")
    log(f"Contigs: {args.contigs}")
    log(f"Host ref: {args.host_ref}")
    log(f"Construct ref: {args.construct_ref}")
    log(f"Output directory: {step_dir}")

    # Map contigs to host
    run_minimap2(args.contigs, args.host_ref, host_paf, args.threads, "host")

    # Map contigs to construct
    run_minimap2(
        args.contigs, args.construct_ref, construct_paf, args.threads, "construct",
    )

    log(f"Host PAF:      {host_paf}")
    log(f"Construct PAF: {construct_paf}")
    log("Done.")


if __name__ == "__main__":
    main()
