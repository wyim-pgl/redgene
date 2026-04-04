#!/usr/bin/env python3
"""Step 6c: Zygosity estimation (homozygous vs heterozygous insertion).

Determines whether a transgene insertion is homozygous or heterozygous
by analyzing read evidence at junction sites.

Methods:
1. Allele fraction: At the junction breakpoint, count reads supporting:
   - Normal allele: reads that span continuously across the junction
   - Insertion allele: reads that are split/clipped at the junction
   Allele fraction = insertion_reads / total_reads
   ~0.5 = heterozygous, ~1.0 = homozygous

2. Depth ratio: Compare depth at the junction vs flanking regions:
   - Heterozygous: junction depth ≈ 50% of flanking (some reads span intact allele)
   - Homozygous: junction depth similar to flanking (all reads come from inserted allele)

3. Mate-pair evidence: Read pairs where one mate maps to host and
   the other maps to construct. In heterozygous, there should be
   similar number of read pairs mapping entirely to host across the
   junction site.

Input: junctions.tsv, host BAM from Step 7
Output: zygosity_report.tsv
"""

import argparse
import subprocess
import sys
from pathlib import Path


def log(msg: str) -> None:
    print(f"[s06c_zygosity] {msg}", file=sys.stderr, flush=True)


def estimate_zygosity(
    host_bam: Path,
    chrom: str,
    junction_pos: int,
    window: int = 200,
    flank: int = 2000,
) -> dict:
    """Estimate zygosity at a junction position.

    Returns dict with:
    - spanning_reads: reads that span continuously across junction (normal allele)
    - clipped_reads: reads that are soft-clipped at junction (insertion allele)
    - allele_fraction: insertion / total
    - depth_at_junction: depth at exact junction position
    - depth_flanking: average depth in flanking regions
    - depth_ratio: junction / flanking
    - zygosity: 'homozygous', 'heterozygous', or 'uncertain'
    """

    # Count spanning vs clipped reads at junction
    spanning = 0
    clipped = 0
    discordant = 0

    result = subprocess.run(
        ["samtools", "view", str(host_bam),
         f"{chrom}:{junction_pos-window}-{junction_pos+window}"],
        capture_output=True, text=True,
    )

    import re

    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue

        flag = int(fields[1])
        pos = int(fields[3])
        cigar = fields[5]
        mapq = int(fields[4])

        if flag & 4 or mapq < 10:  # unmapped or low quality
            continue

        # Calculate alignment end position from CIGAR
        ref_consumed = 0
        for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar):
            length = int(match.group(1))
            op = match.group(2)
            if op in "MDN=X":
                ref_consumed += length

        read_end = pos + ref_consumed

        # Check if read spans the junction
        if pos < junction_pos - 10 and read_end > junction_pos + 10:
            # Read spans across junction - check if it's clipped near junction
            has_clip_at_junction = False
            current_pos = pos
            for match in re.finditer(r"(\d+)([MIDNSHP=X])", cigar):
                length = int(match.group(1))
                op = match.group(2)
                if op == "S" and length >= 15:
                    if abs(current_pos - junction_pos) < 30:
                        has_clip_at_junction = True
                        break
                if op in "MDN=X":
                    current_pos += length

            if has_clip_at_junction:
                clipped += 1
            else:
                spanning += 1

        # Check for soft clips near junction (reads starting/ending at junction)
        elif abs(pos - junction_pos) < 30 or abs(read_end - junction_pos) < 30:
            clips = re.findall(r"(\d+)S", cigar)
            if any(int(c) >= 15 for c in clips):
                clipped += 1

        # Discordant pairs
        if flag & 1 and not (flag & 2) and not (flag & 8):
            mate_chr = fields[6]
            if mate_chr == "=" or mate_chr == chrom:
                mate_pos = int(fields[7])
                if abs(mate_pos - junction_pos) < window:
                    discordant += 1

    # Depth analysis
    depth_at_junction = _get_depth(host_bam, chrom, junction_pos)
    depth_left = _get_avg_depth(host_bam, chrom,
                                max(0, junction_pos - flank),
                                junction_pos - 200)
    depth_right = _get_avg_depth(host_bam, chrom,
                                 junction_pos + 200,
                                 junction_pos + flank)
    depth_flanking = (depth_left + depth_right) / 2 if (depth_left + depth_right) > 0 else 1

    # Calculate allele fraction
    total = spanning + clipped
    allele_fraction = clipped / total if total > 0 else 0

    # Depth ratio
    depth_ratio = depth_at_junction / depth_flanking if depth_flanking > 0 else 0

    # Determine zygosity
    if total < 5:
        zygosity = "uncertain"
        confidence = "Low (< 5 reads at junction)"
    elif allele_fraction >= 0.8:
        zygosity = "homozygous"
        confidence = f"High (allele fraction = {allele_fraction:.2f})"
    elif allele_fraction <= 0.65 and allele_fraction >= 0.3:
        zygosity = "heterozygous"
        confidence = f"High (allele fraction = {allele_fraction:.2f})"
    elif allele_fraction > 0.65:
        zygosity = "likely_homozygous"
        confidence = f"Medium (allele fraction = {allele_fraction:.2f})"
    elif allele_fraction < 0.3 and allele_fraction > 0.1:
        zygosity = "likely_heterozygous"
        confidence = f"Medium (allele fraction = {allele_fraction:.2f})"
    else:
        zygosity = "uncertain"
        confidence = f"Low (allele fraction = {allele_fraction:.2f})"

    return {
        "spanning_reads": spanning,
        "clipped_reads": clipped,
        "discordant_pairs": discordant,
        "allele_fraction": allele_fraction,
        "depth_at_junction": depth_at_junction,
        "depth_flanking_left": depth_left,
        "depth_flanking_right": depth_right,
        "depth_flanking_avg": depth_flanking,
        "depth_ratio": depth_ratio,
        "zygosity": zygosity,
        "confidence": confidence,
    }


def _get_depth(bam: Path, chrom: str, pos: int) -> int:
    """Get depth at a single position."""
    result = subprocess.run(
        ["samtools", "depth", "-r", f"{chrom}:{pos}-{pos+1}", str(bam)],
        capture_output=True, text=True,
    )
    for line in result.stdout.strip().split("\n"):
        if line:
            parts = line.split("\t")
            if len(parts) >= 3:
                return int(parts[2])
    return 0


def _get_avg_depth(bam: Path, chrom: str, start: int, end: int) -> float:
    """Get average depth in a region."""
    if start >= end:
        return 0.0
    result = subprocess.run(
        ["samtools", "depth", "-r", f"{chrom}:{start}-{end}", str(bam)],
        capture_output=True, text=True,
    )
    depths = []
    for line in result.stdout.strip().split("\n"):
        if line:
            parts = line.split("\t")
            if len(parts) >= 3:
                depths.append(int(parts[2]))
    return sum(depths) / len(depths) if depths else 0.0


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 6c: Zygosity estimation at junction sites",
    )
    parser.add_argument(
        "--junctions", required=True, type=Path,
        help="junctions.tsv from Step 6",
    )
    parser.add_argument(
        "--host-bam", required=True, type=Path,
        help="Host BAM from Step 7 (must be indexed)",
    )
    parser.add_argument(
        "--outdir", required=True, type=Path,
        help="Base output directory",
    )
    parser.add_argument(
        "--sample-name", required=True,
        help="Sample name",
    )
    args = parser.parse_args()

    if not args.host_bam.exists():
        log(f"ERROR: Host BAM not found: {args.host_bam}")
        sys.exit(1)

    # Check index
    bai = Path(str(args.host_bam) + ".bai")
    if not bai.exists():
        log("Indexing host BAM...")
        subprocess.run(["samtools", "index", str(args.host_bam)], check=True)

    step_dir = args.outdir / args.sample_name / "s06_junction"
    step_dir.mkdir(parents=True, exist_ok=True)

    # Parse junctions
    junctions = []
    with open(args.junctions) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            junctions.append(dict(zip(header, fields)))

    log(f"Analyzing zygosity for {len(junctions)} junction(s)...")

    # Estimate zygosity for each junction
    results = []
    for j in junctions:
        chrom = j["host_chr"]
        pos = int(j["junction_pos_host"])
        log(f"\n  {chrom}:{pos} ({j['junction_type']}):")

        zyg = estimate_zygosity(args.host_bam, chrom, pos)

        log(f"    Spanning reads (normal allele): {zyg['spanning_reads']}")
        log(f"    Clipped reads (insertion allele): {zyg['clipped_reads']}")
        log(f"    Allele fraction: {zyg['allele_fraction']:.2f}")
        log(f"    Depth (junction/flanking): {zyg['depth_at_junction']}/{zyg['depth_flanking_avg']:.1f}")
        log(f"    Zygosity: {zyg['zygosity']} ({zyg['confidence']})")

        results.append({**j, **zyg})

    # Write output
    output_tsv = step_dir / "zygosity_report.tsv"
    with open(output_tsv, "w") as f:
        out_fields = [
            "contig_name", "host_chr", "junction_pos_host", "junction_type",
            "spanning_reads", "clipped_reads", "allele_fraction",
            "depth_at_junction", "depth_flanking_avg", "depth_ratio",
            "zygosity", "confidence",
        ]
        f.write("\t".join(out_fields) + "\n")
        for r in results:
            f.write("\t".join(str(r.get(k, "")) for k in out_fields) + "\n")

    log(f"\nZygosity report written to: {output_tsv}")
    log("Done.")


if __name__ == "__main__":
    main()
