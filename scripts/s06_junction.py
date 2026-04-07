#!/usr/bin/env python3
"""Step 6: Chimeric contig detection and junction coordinate extraction.

Parses minimap2 PAF files from Step 5 to find contigs that align partially
to both host and construct references. These chimeric contigs span the
transgene insertion junction, and the boundary between host and construct
alignments defines the junction coordinate.

Algorithm:
1. For each contig, collect all host alignments and construct alignments.
2. If a contig has alignments to BOTH host and construct, it is chimeric.
3. The junction position is where the host alignment ends and construct
   alignment begins (or vice versa).
4. Determine LB/RB based on which end of the construct is at the junction.
5. Host junction coordinate = the host reference position at the boundary.

Output:
    {outdir}/{sample}/s06_junction/junctions.tsv
    {outdir}/{sample}/s06_junction/junction_contigs.fa
"""

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path


def log(msg: str) -> None:
    """Print message to stderr."""
    print(f"[s06_junction] {msg}", file=sys.stderr, flush=True)


@dataclass
class PafAlignment:
    """A single alignment record from a PAF file."""

    query_name: str
    query_len: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_len: int
    target_start: int
    target_end: int
    matching_bases: int
    alignment_block_len: int
    mapping_quality: int

    @classmethod
    def from_line(cls, line: str) -> "PafAlignment | None":
        """Parse a PAF line into a PafAlignment object."""
        fields = line.strip().split("\t")
        if len(fields) < 12:
            return None
        return cls(
            query_name=fields[0],
            query_len=int(fields[1]),
            query_start=int(fields[2]),
            query_end=int(fields[3]),
            strand=fields[4],
            target_name=fields[5],
            target_len=int(fields[6]),
            target_start=int(fields[7]),
            target_end=int(fields[8]),
            matching_bases=int(fields[9]),
            alignment_block_len=int(fields[10]),
            mapping_quality=int(fields[11]),
        )

    @property
    def query_aligned_frac(self) -> float:
        """Fraction of query sequence covered by this alignment."""
        if self.query_len == 0:
            return 0.0
        return (self.query_end - self.query_start) / self.query_len

    @property
    def identity(self) -> float:
        """Alignment identity (matching bases / block length)."""
        if self.alignment_block_len == 0:
            return 0.0
        return self.matching_bases / self.alignment_block_len


@dataclass
class Junction:
    """A detected host-transgene junction."""

    contig_name: str
    contig_len: int
    host_chr: str
    host_start: int
    host_end: int
    host_strand: str
    construct_element: str
    construct_start: int
    construct_end: int
    junction_pos_host: int
    junction_type: str  # LB or RB
    confidence: str  # High, Medium, Low
    host_mapq: int = 0  # MAPQ of host alignment

    def to_tsv_row(self) -> str:
        """Format as a TSV row."""
        return "\t".join([
            self.contig_name,
            str(self.contig_len),
            self.host_chr,
            str(self.host_start),
            str(self.host_end),
            self.host_strand,
            self.construct_element,
            str(self.construct_start),
            str(self.construct_end),
            str(self.junction_pos_host),
            self.junction_type,
            self.confidence,
            str(self.host_mapq),
        ])


TSV_HEADER = "\t".join([
    "contig_name",
    "contig_len",
    "host_chr",
    "host_start",
    "host_end",
    "host_strand",
    "construct_element",
    "construct_start",
    "construct_end",
    "junction_pos_host",
    "junction_type",
    "confidence",
    "host_mapq",
])


def parse_paf(paf_path: Path) -> dict[str, list[PafAlignment]]:
    """Parse a PAF file and group alignments by query name."""
    alignments: dict[str, list[PafAlignment]] = {}
    with open(paf_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            aln = PafAlignment.from_line(line)
            if aln is None:
                continue
            alignments.setdefault(aln.query_name, []).append(aln)
    return alignments


def find_chimeric_contigs(
    host_alns: dict[str, list[PafAlignment]],
    construct_alns: dict[str, list[PafAlignment]],
    min_aligned_frac: float = 0.05,
    min_identity: float = 0.90,
    min_host_mapq: int = 10,
    min_host_aln_len: int = 50,
    min_construct_aln_len: int = 50,
    min_coverage_frac: float = 0.5,
) -> list[Junction]:
    """Find contigs with alignments to both host and construct.

    For each chimeric contig, determine the junction coordinates and type.
    Applies multiple filters to reduce false positives:
    - MAPQ: host alignment mapping quality (uniqueness)
    - Minimum alignment length for both host and construct
    - Minimum combined coverage of the contig

    Args:
        host_alns: Host alignments grouped by contig name.
        construct_alns: Construct alignments grouped by contig name.
        min_aligned_frac: Minimum fraction of contig aligned to count.
        min_identity: Minimum alignment identity to consider.
        min_host_mapq: Minimum MAPQ for host alignments (default: 10).
        min_host_aln_len: Minimum host alignment length in bp (default: 50).
        min_construct_aln_len: Minimum construct alignment length in bp (default: 50).
        min_coverage_frac: Minimum fraction of contig covered by host+construct (default: 0.5).

    Returns:
        List of Junction objects.
    """
    # Find contig names present in both sets
    chimeric_names = set(host_alns.keys()) & set(construct_alns.keys())
    if not chimeric_names:
        log("WARNING: No chimeric contigs found (no contig maps to both host and construct)")
        return []

    log(f"Found {len(chimeric_names)} contig(s) with alignments to both host and construct")

    junctions: list[Junction] = []

    for contig_name in sorted(chimeric_names):
        h_alns = host_alns[contig_name]
        c_alns = construct_alns[contig_name]

        # Filter by identity
        h_alns = [a for a in h_alns if a.identity >= min_identity]
        c_alns = [a for a in c_alns if a.identity >= min_identity]

        if not h_alns or not c_alns:
            continue

        # Filter host alignments by MAPQ (quality of uniqueness)
        h_alns_filtered = [a for a in h_alns if a.mapping_quality >= min_host_mapq]
        if h_alns_filtered:
            h_alns = h_alns_filtered
        else:
            log(f"  FILTERED: {contig_name} - no host alignments with MAPQ>={min_host_mapq} (max MAPQ={max(a.mapping_quality for a in h_alns)})")
            continue  # Skip contigs with no reliable host alignments

        # Filter by minimum alignment length
        h_alns = [a for a in h_alns if (a.query_end - a.query_start) >= min_host_aln_len]
        c_alns = [a for a in c_alns if (a.query_end - a.query_start) >= min_construct_aln_len]

        if not h_alns or not c_alns:
            continue

        # Filter by minimum aligned fraction
        h_alns = [a for a in h_alns if a.query_aligned_frac >= min_aligned_frac]
        c_alns = [a for a in c_alns if a.query_aligned_frac >= min_aligned_frac]

        if not h_alns or not c_alns:
            continue

        # Compute union of all construct alignment intervals on this contig
        # (element_db has many small reference sequences, so we merge them)
        c_intervals = sorted([(a.query_start, a.query_end) for a in c_alns])
        c_merged: list[tuple[int, int]] = []
        for start, end in c_intervals:
            if c_merged and start <= c_merged[-1][1]:
                c_merged[-1] = (c_merged[-1][0], max(c_merged[-1][1], end))
            else:
                c_merged.append((start, end))
        total_construct_cov = sum(e - s for s, e in c_merged)

        # For each pair of host/construct alignments, check if they cover
        # complementary portions of the contig (i.e., non-overlapping regions)
        for h_aln in h_alns:
            for c_aln in c_alns:
                # Check that the two alignments cover different parts of the contig
                # Allow some overlap (up to 50 bp) at the junction
                overlap = min(h_aln.query_end, c_aln.query_end) - max(
                    h_aln.query_start, c_aln.query_start,
                )
                if overlap > 50:
                    # Too much overlap -- these alignments cover the same region
                    continue

                # Check combined coverage using union of all construct alignments
                # + this specific host alignment (not just the single c_aln)
                h_len = h_aln.query_end - h_aln.query_start
                # Subtract overlap between host and merged construct intervals
                h_construct_overlap = 0
                for cs, ce in c_merged:
                    ov = min(h_aln.query_end, ce) - max(h_aln.query_start, cs)
                    if ov > 0:
                        h_construct_overlap += ov
                combined_coverage = (h_len + total_construct_cov - h_construct_overlap) / h_aln.query_len
                if combined_coverage < min_coverage_frac:
                    log(f"  FILTERED: {contig_name} - low combined coverage "
                        f"({combined_coverage:.1%} < {min_coverage_frac:.0%})")
                    continue

                # Determine which comes first on the contig
                if h_aln.query_start < c_aln.query_start:
                    # Host is on the left, construct on the right
                    # Junction is at the boundary
                    junction_query_pos = (h_aln.query_end + c_aln.query_start) // 2

                    # Host junction coordinate: the host position at the boundary
                    if h_aln.strand == "+":
                        junction_pos_host = h_aln.target_end
                    else:
                        junction_pos_host = h_aln.target_start

                    # Determine LB vs RB: if construct alignment starts near
                    # position 0 of the construct, this is a Left Border (LB)
                    # junction; if it starts near the end, it is RB.
                    construct_len = c_aln.target_len
                    if c_aln.strand == "+":
                        construct_junction_pos = c_aln.target_start
                    else:
                        construct_junction_pos = c_aln.target_end

                    if construct_junction_pos < construct_len * 0.3:
                        junction_type = "LB"
                    elif construct_junction_pos > construct_len * 0.7:
                        junction_type = "RB"
                    else:
                        junction_type = "LB"  # Default; internal junction
                else:
                    # Construct is on the left, host on the right
                    junction_query_pos = (c_aln.query_end + h_aln.query_start) // 2

                    if h_aln.strand == "+":
                        junction_pos_host = h_aln.target_start
                    else:
                        junction_pos_host = h_aln.target_end

                    construct_len = c_aln.target_len
                    if c_aln.strand == "+":
                        construct_junction_pos = c_aln.target_end
                    else:
                        construct_junction_pos = c_aln.target_start

                    if construct_junction_pos > construct_len * 0.7:
                        junction_type = "RB"
                    elif construct_junction_pos < construct_len * 0.3:
                        junction_type = "LB"
                    else:
                        junction_type = "RB"  # Default; internal junction

                junction = Junction(
                    contig_name=contig_name,
                    contig_len=h_aln.query_len,
                    host_chr=h_aln.target_name,
                    host_start=h_aln.target_start,
                    host_end=h_aln.target_end,
                    host_strand=h_aln.strand,
                    construct_element=c_aln.target_name,
                    construct_start=c_aln.target_start,
                    construct_end=c_aln.target_end,
                    junction_pos_host=junction_pos_host,
                    junction_type=junction_type,
                    confidence="Medium",  # Updated later
                    host_mapq=h_aln.mapping_quality,
                )
                junctions.append(junction)

    return junctions


def _check_rearrangement_warnings(junctions: list[Junction]) -> None:
    """Check for patterns that suggest chromosomal rearrangements or
    potential false positives, and emit warnings.

    Checks:
    1. Same contig mapping to different chromosomes → translocation
    2. Same contig mapping to same chromosome but far apart → deletion/inversion
    3. Multiple junctions with low MAPQ → homologous sequence artifacts
    """
    from collections import defaultdict

    # Group junctions by contig
    by_contig: dict[str, list[Junction]] = defaultdict(list)
    for j in junctions:
        by_contig[j.contig_name].append(j)

    warnings: list[str] = []

    for contig, jcts in by_contig.items():
        if len(jcts) < 2:
            continue

        chromosomes = {j.host_chr for j in jcts}
        positions = [(j.host_chr, j.junction_pos_host) for j in jcts]

        # Check 1: Inter-chromosomal (translocation)
        if len(chromosomes) > 1:
            chrom_list = ", ".join(sorted(chromosomes))
            warnings.append(
                f"WARNING: Contig {contig} maps to multiple chromosomes ({chrom_list}).\n"
                f"  Possible explanations:\n"
                f"  1. Inter-chromosomal translocation mediated by T-DNA insertion\n"
                f"  2. Assembly artifact: reads from different loci merged into one contig\n"
                f"  3. Homologous sequence between construct and multiple host locations\n"
                f"  → Verify with WT data filtering or PCR validation"
            )

        # Check 2: Intra-chromosomal distant junctions
        for i in range(len(jcts)):
            for k in range(i + 1, len(jcts)):
                if jcts[i].host_chr == jcts[k].host_chr:
                    distance = abs(jcts[i].junction_pos_host - jcts[k].junction_pos_host)
                    if distance > 100_000:  # > 100kb apart
                        dist_str = f"{distance/1_000_000:.1f} Mb" if distance > 1_000_000 else f"{distance/1_000:.0f} kb"
                        warnings.append(
                            f"WARNING: Contig {contig} has junctions on {jcts[i].host_chr} "
                            f"separated by {dist_str} ({jcts[i].junction_pos_host:,} and {jcts[k].junction_pos_host:,}).\n"
                            f"  Possible explanations:\n"
                            f"  1. Large chromosomal deletion ({dist_str}) at insertion site\n"
                            f"  2. Chromosomal inversion mediated by T-DNA\n"
                            f"  3. Complex insertion with filler DNA from distant genomic region\n"
                            f"  4. Assembly artifact or homologous sequence false positive\n"
                            f"  → Requires WT filtering, depth analysis, or long-read validation"
                        )

    # Check 3: Multiple low-MAPQ junctions (systematic false positive pattern)
    low_mapq = [j for j in junctions if j.host_mapq < 10]
    if len(low_mapq) >= 2:
        warnings.append(
            f"WARNING: {len(low_mapq)} junction(s) have MAPQ < 10.\n"
            f"  This may indicate homologous sequences between construct and host genome.\n"
            f"  Plant genomes often contain duplicated/repetitive regions that share\n"
            f"  partial homology with transformation vector sequences.\n"
            f"  → Essential to run WT-based filtering (Step 3b) to distinguish true insertions"
        )

    # Emit all warnings
    if warnings:
        log("")
        log("=" * 60)
        log("STRUCTURAL WARNINGS")
        log("=" * 60)
        for w in warnings:
            for line in w.split("\n"):
                log(line)
            log("")
        log("NOTE: These warnings do not necessarily mean the junctions are false.")
        log("They flag patterns that require additional verification.")
        log("Use Step 3b (WT filtering) and Step 6b (verification) to validate.")
        log("=" * 60)
        log("")


def assign_confidence(junctions: list[Junction]) -> None:
    """Assign confidence levels based on junction types found.

    - High: both LB and RB junction contigs found (bp-level resolution)
    - Medium: only one junction type found
    - Low: ambiguous or overlapping alignments
    """
    junction_types_found = {j.junction_type for j in junctions}

    if "LB" in junction_types_found and "RB" in junction_types_found:
        confidence = "High"
    elif len(junction_types_found) == 1:
        confidence = "Medium"
    else:
        confidence = "Low"

    for j in junctions:
        j.confidence = confidence


def read_fasta(fasta_path: Path) -> dict[str, str]:
    """Read a FASTA file and return a dict of name -> sequence."""
    sequences: dict[str, str] = {}
    current_name = ""
    current_seq: list[str] = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            sequences[current_name] = "".join(current_seq)

    return sequences


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 6: Chimeric contig detection and junction extraction",
    )
    parser.add_argument(
        "--host-paf", required=True, type=Path,
        help="PAF file from contigs-to-host mapping (Step 5)",
    )
    parser.add_argument(
        "--construct-paf", required=True, type=Path,
        help="PAF file from contigs-to-construct mapping (Step 5)",
    )
    parser.add_argument(
        "--contigs", required=True, type=Path,
        help="Assembled contigs FASTA (from Step 4)",
    )
    parser.add_argument(
        "--outdir", required=True, type=Path,
        help="Base output directory (results/)",
    )
    parser.add_argument(
        "--sample-name", required=True,
        help="Sample name for output organization",
    )
    parser.add_argument(
        "--min-identity", type=float, default=0.90,
        help="Minimum alignment identity to consider (default: 0.90)",
    )
    parser.add_argument(
        "--min-aligned-frac", type=float, default=0.05,
        help="Minimum fraction of contig aligned (default: 0.05)",
    )
    parser.add_argument(
        "--min-host-mapq", type=int, default=10,
        help="Minimum MAPQ for host alignments (default: 10, filters non-unique mappings)",
    )
    parser.add_argument(
        "--min-host-aln-len", type=int, default=50,
        help="Minimum host alignment length in bp (default: 50)",
    )
    parser.add_argument(
        "--min-construct-aln-len", type=int, default=50,
        help="Minimum construct alignment length in bp (default: 50)",
    )
    parser.add_argument(
        "--min-coverage-frac", type=float, default=0.5,
        help="Minimum fraction of contig covered by host+construct (default: 0.5)",
    )
    args = parser.parse_args()

    # Validate inputs
    for label, path in [
        ("host PAF", args.host_paf),
        ("construct PAF", args.construct_paf),
        ("contigs", args.contigs),
    ]:
        if not path.exists():
            log(f"ERROR: {label} file not found: {path}")
            sys.exit(1)

    step_dir = args.outdir / args.sample_name / "s06_junction"
    step_dir.mkdir(parents=True, exist_ok=True)

    junctions_tsv = step_dir / "junctions.tsv"
    junction_contigs_fa = step_dir / "junction_contigs.fa"

    log(f"Sample: {args.sample_name}")
    log(f"Host PAF: {args.host_paf}")
    log(f"Construct PAF: {args.construct_paf}")
    log(f"Contigs: {args.contigs}")

    # Parse PAF files
    log("Parsing host PAF alignments...")
    host_alns = parse_paf(args.host_paf)
    log(f"  {sum(len(v) for v in host_alns.values())} alignments for {len(host_alns)} contigs")

    log("Parsing construct PAF alignments...")
    construct_alns = parse_paf(args.construct_paf)
    log(f"  {sum(len(v) for v in construct_alns.values())} alignments for {len(construct_alns)} contigs")

    # Find chimeric contigs and extract junctions
    log("Detecting chimeric contigs...")
    junctions = find_chimeric_contigs(
        host_alns,
        construct_alns,
        min_aligned_frac=args.min_aligned_frac,
        min_identity=args.min_identity,
        min_host_mapq=args.min_host_mapq,
        min_host_aln_len=args.min_host_aln_len,
        min_construct_aln_len=args.min_construct_aln_len,
        min_coverage_frac=args.min_coverage_frac,
    )

    if not junctions:
        log("WARNING: No junctions detected.")
        log("This could mean:")
        log("  - No contigs span the host-transgene boundary")
        log("  - Assembly did not produce junction-spanning contigs")
        log("  - The construct may not be inserted in the host genome")
        # Write empty output files
        with open(junctions_tsv, "w") as fh:
            fh.write(TSV_HEADER + "\n")
        with open(junction_contigs_fa, "w") as fh:
            pass
        log(f"Empty junctions file written to: {junctions_tsv}")
        log("Done.")
        return

    # Assign confidence
    assign_confidence(junctions)

    # Check for potential chromosomal rearrangements and warn
    _check_rearrangement_warnings(junctions)

    # Write junctions TSV
    with open(junctions_tsv, "w") as fh:
        fh.write(TSV_HEADER + "\n")
        for j in junctions:
            fh.write(j.to_tsv_row() + "\n")
    log(f"Wrote {len(junctions)} junction(s) to: {junctions_tsv}")

    # Extract chimeric contig sequences
    log("Extracting chimeric contig sequences...")
    all_seqs = read_fasta(args.contigs)
    chimeric_names = {j.contig_name for j in junctions}

    with open(junction_contigs_fa, "w") as fh:
        for name in sorted(chimeric_names):
            if name in all_seqs:
                fh.write(f">{name}\n")
                seq = all_seqs[name]
                # Write sequence in 80-character lines
                for i in range(0, len(seq), 80):
                    fh.write(seq[i : i + 80] + "\n")
    log(f"Wrote {len(chimeric_names)} chimeric contig(s) to: {junction_contigs_fa}")

    # Summary report
    junction_types = {}
    for j in junctions:
        junction_types.setdefault(j.junction_type, []).append(j)

    log("Junction summary:")
    for jtype in sorted(junction_types):
        count = len(junction_types[jtype])
        log(f"  {jtype}: {count} junction(s)")
    log(f"  Confidence: {junctions[0].confidence}")

    for j in junctions:
        log(
            f"  {j.contig_name} ({j.contig_len}bp): "
            f"{j.host_chr}:{j.junction_pos_host} <-> "
            f"{j.construct_element}:{j.construct_start}-{j.construct_end} "
            f"[{j.junction_type}]"
        )

    log("Done.")


if __name__ == "__main__":
    main()
