#!/usr/bin/env python3
"""Step 3b: Homologous sequence filter (modified TranSeq approach).

Identifies and filters reads originating from homologous sequences between
the construct/vector and host genome. This is critical for plant genomes
which contain many duplicated and repetitive sequences.

Approach (from Bae et al. 2022, Communications Biology):
1. BLAST construct reference against host genome to find homologous regions
2. For each extracted read pair (from Step 3), check if the host-mapping mate
   falls within a known homologous region (±insert_size)
3. Discard homologous read pairs, retain only true transgenic candidates

This step is essential because:
- Plant T-DNA constructs often contain host-derived promoters (e.g., rice Ubi1,
  Gt1, Actin1) that create false positive chimeric reads
- Backbone sequences may have partial homology to host genome
- Without filtering, >75% of extracted reads can be false positives

Input: Extracted reads from Step 3, construct BAM from Step 2
Output: Filtered reads (true candidates only), homology report

Usage:
  python s03b_homology_filter.py \
    --construct-ref db/construct.fa \
    --host-ref db/host.fa \
    --construct-bam results/{sample}/s02_construct_map/{sample}_construct.bam \
    --r1 results/{sample}/s03_extract/{sample}_construct_R1.fq.gz \
    --r2 results/{sample}/s03_extract/{sample}_construct_R2.fq.gz \
    --outdir results --sample-name {sample}
"""

import argparse
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def log(msg: str) -> None:
    print(f"[s03b_homfilter] {msg}", file=sys.stderr, flush=True)


def find_homologous_regions(
    construct_ref: Path,
    host_ref: Path,
    min_identity: float = 0.80,
    min_length: int = 50,
) -> dict[str, list[tuple[int, int]]]:
    """Use minimap2 to find homologous regions between construct and host genome.

    Much faster than BLAST and already a pipeline dependency.
    minimap2 -x asm5 maps construct sequences to host genome.

    Returns dict of {host_chr: [(start, end), ...]} for homologous regions.
    """
    log("Finding construct-host homologous regions with minimap2...")

    result = subprocess.run(
        [
            "minimap2",
            "-x", "asm5",      # assembly-to-assembly alignment
            "-N", "100",       # report up to 100 secondary alignments
            "--secondary=yes", # include secondary alignments
            str(host_ref),
            str(construct_ref),
        ],
        capture_output=True, text=True,
    )

    homologous: dict[str, list[tuple[int, int]]] = defaultdict(list)
    total_hits = 0

    for line in result.stdout.strip().split("\n"):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue

        # PAF format: query_name, query_len, query_start, query_end, strand,
        #             target_name, target_len, target_start, target_end,
        #             matching_bases, alignment_block_len, mapping_quality
        aln_len = int(fields[3]) - int(fields[2])  # query alignment length
        matching = int(fields[9])
        block_len = int(fields[10])
        identity = matching / block_len if block_len > 0 else 0

        host_chr = fields[5]
        host_start = int(fields[7])
        host_end = int(fields[8])

        if identity >= min_identity and aln_len >= min_length:
            homologous[host_chr].append((host_start, host_end))
            total_hits += 1

    log(f"Found {total_hits} homologous regions across {len(homologous)} chromosomes")

    # Report per-chromosome
    for chrom in sorted(homologous.keys()):
        regions = homologous[chrom]
        total_bp = sum(e - s for s, e in regions)
        log(f"  {chrom}: {len(regions)} regions, {total_bp:,} bp total")

    return dict(homologous)


def is_in_homologous_region(
    chrom: str, pos: int,
    homologous: dict[str, list[tuple[int, int]]],
    insert_size: int = 500,
) -> bool:
    """Check if a position falls within or near a known homologous region."""
    if chrom not in homologous:
        return False

    for start, end in homologous[chrom]:
        if start - insert_size <= pos <= end + insert_size:
            return True
    return False


def filter_reads_by_homology(
    host_ref: Path,
    r1: Path,
    r2: Path,
    homologous: dict[str, list[tuple[int, int]]],
    outdir: Path,
    sample_name: str,
    insert_size: int = 500,
) -> tuple[int, int, dict[str, int]]:
    """Filter extracted read pairs by checking if they originate from
    homologous regions.

    For each read pair extracted in Step 3:
    - One mate maps to construct (from Step 2)
    - The other mate maps to host genome
    - If the host-mapping mate falls in a homologous region: DISCARD
    - Otherwise: RETAIN as true transgenic candidate

    Returns: (total_pairs, retained_pairs, discarded_by_chr)
    """
    step_dir = outdir / sample_name / "s03_extract"
    step_dir.mkdir(parents=True, exist_ok=True)

    filtered_r1 = step_dir / f"{sample_name}_filtered_R1.fq.gz"
    filtered_r2 = step_dir / f"{sample_name}_filtered_R2.fq.gz"
    discarded_r1 = step_dir / f"{sample_name}_homologous_R1.fq.gz"
    discarded_r2 = step_dir / f"{sample_name}_homologous_R2.fq.gz"

    # Get the host mapping positions for construct-hitting reads
    # We need to check the host BAM, but at this stage we may not have it
    # Instead, we can map the extracted reads to host quickly with bwa mem
    log("Mapping extracted reads to host genome for position check...")

    # Map R1+R2 to host genome
    host_bam = step_dir / f"{sample_name}_extract_to_host.bam"
    cmd = (
        f"bwa mem -t 8 {host_ref} {r1} {r2} 2>/dev/null | "
        f"samtools sort -@ 4 -o {host_bam} -"
    )
    subprocess.run(cmd, shell=True, check=True)
    subprocess.run(["samtools", "index", str(host_bam)], check=True)

    # Read all positions from the host BAM
    log("Checking read positions against homologous regions...")
    result = subprocess.run(
        ["samtools", "view", str(host_bam)],
        capture_output=True, text=True,
    )

    # Collect read names that fall in homologous regions
    discard_names: set[str] = set()
    discard_by_chr: dict[str, int] = defaultdict(int)
    total_reads = 0

    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue

        read_name = fields[0]
        flag = int(fields[1])
        chrom = fields[2]
        pos = int(fields[3])

        if flag & 4:  # unmapped
            continue

        total_reads += 1

        if is_in_homologous_region(chrom, pos, homologous, insert_size):
            discard_names.add(read_name)
            discard_by_chr[chrom] += 1

    total_pairs = total_reads // 2
    discarded_pairs = len(discard_names)
    retained_pairs = total_pairs - discarded_pairs

    log(f"Total read pairs: {total_pairs}")
    log(f"Discarded (homologous): {discarded_pairs} ({discarded_pairs/max(total_pairs,1)*100:.1f}%)")
    log(f"Retained (true candidates): {retained_pairs} ({retained_pairs/max(total_pairs,1)*100:.1f}%)")

    for chrom in sorted(discard_by_chr):
        log(f"  Discarded from {chrom}: {discard_by_chr[chrom]}")

    # Filter FASTQ files
    log("Writing filtered FASTQ files...")
    _filter_fastq(r1, filtered_r1, discard_names)
    _filter_fastq(r2, filtered_r2, discard_names)

    # Write discarded reads too (for transparency)
    _filter_fastq(r1, discarded_r1, discard_names, keep_matching=True)
    _filter_fastq(r2, discarded_r2, discard_names, keep_matching=True)

    log(f"Filtered reads: {filtered_r1}")
    log(f"Discarded reads: {discarded_r1}")

    # Write report
    report_file = step_dir / "homology_filter_report.txt"
    with open(report_file, "w") as f:
        f.write(f"Homology Filter Report - {sample_name}\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total extracted read pairs: {total_pairs}\n")
        f.write(f"Discarded (homologous): {discarded_pairs} ({discarded_pairs/max(total_pairs,1)*100:.1f}%)\n")
        f.write(f"Retained (true candidates): {retained_pairs} ({retained_pairs/max(total_pairs,1)*100:.1f}%)\n\n")
        f.write("Discarded by chromosome:\n")
        for chrom in sorted(discard_by_chr):
            f.write(f"  {chrom}: {discard_by_chr[chrom]}\n")
        f.write("\nHomologous regions detected:\n")
        for chrom in sorted(homologous):
            for start, end in homologous[chrom]:
                f.write(f"  {chrom}:{start}-{end} ({end-start} bp)\n")

    log(f"Report: {report_file}")

    # Clean up
    host_bam.unlink(missing_ok=True)
    Path(str(host_bam) + ".bai").unlink(missing_ok=True)

    return total_pairs, retained_pairs, dict(discard_by_chr)


def _filter_fastq(
    input_fq: Path, output_fq: Path,
    names: set[str], keep_matching: bool = False,
) -> None:
    """Filter a FASTQ file, keeping or discarding reads by name."""
    import gzip

    with gzip.open(input_fq, "rt") as fin, gzip.open(output_fq, "wt") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            name = header.split()[0].lstrip("@").split("/")[0]
            in_set = name in names

            if (keep_matching and in_set) or (not keep_matching and not in_set):
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 3b: Homologous sequence filter",
    )
    parser.add_argument(
        "--construct-ref", required=True, type=Path,
        help="Construct/vector reference FASTA",
    )
    parser.add_argument(
        "--host-ref", required=True, type=Path,
        help="Host reference genome FASTA",
    )
    parser.add_argument(
        "--r1", required=True, type=Path,
        help="Extracted R1 FASTQ from Step 3",
    )
    parser.add_argument(
        "--r2", required=True, type=Path,
        help="Extracted R2 FASTQ from Step 3",
    )
    parser.add_argument(
        "--outdir", required=True, type=Path,
        help="Base output directory",
    )
    parser.add_argument(
        "--sample-name", required=True,
        help="Sample name",
    )
    parser.add_argument(
        "--insert-size", type=int, default=500,
        help="Insert size for proximity check (default: 500)",
    )
    parser.add_argument(
        "--min-identity", type=float, default=0.80,
        help="Minimum identity for homologous regions (default: 0.80 = 80%%)",
    )
    parser.add_argument(
        "--min-aln-length", type=int, default=50,
        help="Minimum alignment length for homologous regions (default: 50)",
    )
    args = parser.parse_args()

    # Step 1: Find homologous regions
    homologous = find_homologous_regions(
        args.construct_ref, args.host_ref,
        min_identity=args.min_identity,
        min_length=args.min_aln_length,
    )

    if not homologous:
        log("No homologous regions found. All reads retained.")
        return

    # Step 2: Filter reads
    total, retained, discarded = filter_reads_by_homology(
        host_ref=args.host_ref,
        r1=args.r1,
        r2=args.r2,
        homologous=homologous,
        outdir=args.outdir,
        sample_name=args.sample_name,
        insert_size=args.insert_size,
    )

    log("Done.")


if __name__ == "__main__":
    main()
