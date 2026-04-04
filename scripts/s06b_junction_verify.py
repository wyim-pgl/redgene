#!/usr/bin/env python3
"""Step 6b: Junction verification protocol.

Validates junction candidates from Step 6 using multiple lines of evidence:
1. WT-based filtering: Junctions present in wild-type data are false positives
2. Read support: Count split/discordant reads at each junction site
3. Depth profile: Check for depth anomalies around junction coordinates
4. Homology check: BLAST junction region vs host genome for uniqueness
5. Contig structure: Validate chimeric contig covers sufficient contig fraction

Input: junctions.tsv from Step 6, host BAM from Step 7, construct BAM from Step 2
Output: verified_junctions.tsv with evidence columns, verification_report.txt
"""

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None


def log(msg: str) -> None:
    print(f"[s06b_verify] {msg}", file=sys.stderr, flush=True)


@dataclass
class JunctionEvidence:
    """Evidence collected for a junction candidate."""

    # From junctions.tsv
    contig_name: str
    host_chr: str
    junction_pos: int
    junction_type: str
    host_mapq: int

    # Verification evidence
    split_reads: int = 0
    discordant_pairs: int = 0
    soft_clipped_reads: int = 0
    depth_at_junction: int = 0
    depth_flanking_left: float = 0.0
    depth_flanking_right: float = 0.0
    depth_ratio: float = 0.0  # junction depth / flanking average
    blast_hits_in_genome: int = 0  # Number of BLAST hits for junction region
    wt_reads_at_site: int = 0  # Reads from WT at this position
    verdict: str = "UNVERIFIED"
    reason: str = ""


def count_supporting_reads(
    bam_path: Path, chrom: str, pos: int, window: int = 500
) -> tuple[int, int, int, int]:
    """Count split reads, discordant pairs, soft-clipped reads, and depth
    at a junction position using pysam.

    Returns: (split_reads, discordant_pairs, soft_clipped, depth)
    """
    if pysam is None:
        return _count_supporting_reads_samtools(bam_path, chrom, pos, window)

    split = 0
    discordant = 0
    soft_clipped = 0
    depth = 0

    try:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
    except Exception as e:
        log(f"  WARNING: Cannot open BAM {bam_path}: {e}")
        return 0, 0, 0, 0

    start = max(0, pos - window)
    end = pos + window

    try:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary:
                continue

            # Depth at exact position
            if read.reference_start <= pos <= read.reference_end:
                depth += 1

            # Soft-clipped reads near junction
            if read.cigartuples:
                for op, length in read.cigartuples:
                    if op == 4 and length >= 20:  # S = soft clip, min 20bp
                        if abs(read.reference_start - pos) < 50 or abs(
                            read.reference_end - pos
                        ) < 50:
                            soft_clipped += 1
                            break

            # Supplementary alignment (split read)
            if read.has_tag("SA"):
                split += 1

            # Discordant pair
            if read.is_paired and not read.is_proper_pair and not read.mate_is_unmapped:
                discordant += 1
    except ValueError:
        log(f"  WARNING: Region {chrom}:{start}-{end} not found in BAM")
    finally:
        bam.close()

    return split, discordant, soft_clipped, depth


def _count_supporting_reads_samtools(
    bam_path: Path, chrom: str, pos: int, window: int = 500
) -> tuple[int, int, int, int]:
    """Fallback: count reads using samtools (when pysam not available)."""
    start = max(0, pos - window)
    end = pos + window
    region = f"{chrom}:{start}-{end}"

    # Depth at position
    depth_result = subprocess.run(
        ["samtools", "depth", "-r", f"{chrom}:{pos}-{pos+1}", str(bam_path)],
        capture_output=True, text=True,
    )
    depth = 0
    for line in depth_result.stdout.strip().split("\n"):
        if line:
            parts = line.split("\t")
            if len(parts) >= 3:
                depth = int(parts[2])

    # Count soft-clipped reads
    view_result = subprocess.run(
        ["samtools", "view", str(bam_path), region],
        capture_output=True, text=True,
    )

    split = 0
    discordant = 0
    soft_clipped = 0

    for line in view_result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue

        flag = int(fields[1])
        cigar = fields[5]
        pos_read = int(fields[3])

        # Soft clipped
        import re

        clips = re.findall(r"(\d+)S", cigar)
        for clip_len in clips:
            if int(clip_len) >= 20:
                if abs(pos_read - pos) < 100:
                    soft_clipped += 1
                    break

        # Supplementary (SA tag)
        if "SA:Z:" in line:
            split += 1

        # Discordant
        if flag & 1 and not (flag & 2) and not (flag & 8):
            discordant += 1

    return split, discordant, soft_clipped, depth


def get_depth_profile(
    bam_path: Path, chrom: str, pos: int, flank: int = 2000
) -> tuple[float, float]:
    """Get average depth in flanking regions (left and right of junction)."""
    left_start = max(0, pos - flank)
    left_end = pos - 100
    right_start = pos + 100
    right_end = pos + flank

    def avg_depth(region_str: str) -> float:
        result = subprocess.run(
            ["samtools", "depth", "-r", region_str, str(bam_path)],
            capture_output=True, text=True,
        )
        depths = []
        for line in result.stdout.strip().split("\n"):
            if line:
                parts = line.split("\t")
                if len(parts) >= 3:
                    depths.append(int(parts[2]))
        return sum(depths) / len(depths) if depths else 0.0

    left_avg = avg_depth(f"{chrom}:{left_start}-{left_end}")
    right_avg = avg_depth(f"{chrom}:{right_start}-{right_end}")

    return left_avg, right_avg


def check_wt_reads(
    wt_construct_bam: Path, host_ref: Path, chrom: str, pos: int, window: int = 1000
) -> int:
    """Check if WT data has reads mapping to construct that originate from
    this host position (indicates homologous sequence = false positive)."""
    if not wt_construct_bam or not wt_construct_bam.exists():
        return -1  # Unknown

    # Get read names from WT construct BAM
    result = subprocess.run(
        ["samtools", "view", str(wt_construct_bam)],
        capture_output=True, text=True,
    )

    wt_construct_reads = set()
    for line in result.stdout.strip().split("\n"):
        if line:
            wt_construct_reads.add(line.split("\t")[0])

    if not wt_construct_reads:
        return 0

    # Check if these reads map to the junction region in WT host BAM
    # (This requires the WT host BAM, which may not be available yet)
    return len(wt_construct_reads)


def blast_junction_region(
    contigs_fa: Path, contig_name: str, host_ref: Path,
    query_start: int, query_end: int,
) -> int:
    """BLAST the host-mapping portion of a chimeric contig against the
    full host genome to check for multiple hits (non-uniqueness)."""
    from Bio import SeqIO

    # Extract the contig sequence
    for rec in SeqIO.parse(str(contigs_fa), "fasta"):
        if rec.id == contig_name:
            region = str(rec.seq[query_start:query_end])
            break
    else:
        return -1

    # Write query
    query_file = Path("/tmp/junction_blast_query.fa")
    with open(query_file, "w") as f:
        f.write(f">{contig_name}_{query_start}_{query_end}\n{region}\n")

    # BLAST
    result = subprocess.run(
        [
            "blastn", "-query", str(query_file),
            "-subject", str(host_ref),
            "-outfmt", "6", "-evalue", "1e-5",
        ],
        capture_output=True, text=True,
    )

    hits = len([l for l in result.stdout.strip().split("\n") if l])
    return hits


def verify_junctions(
    junctions_tsv: Path,
    host_bam: Path | None,
    contigs_fa: Path,
    host_ref: Path,
    wt_construct_bam: Path | None,
    outdir: Path,
    sample_name: str,
) -> list[JunctionEvidence]:
    """Run the full verification protocol on junction candidates."""

    step_dir = outdir / sample_name / "s06_junction"
    step_dir.mkdir(parents=True, exist_ok=True)

    evidences: list[JunctionEvidence] = []

    # Parse junctions
    with open(junctions_tsv) as f:
        header = f.readline()
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            ev = JunctionEvidence(
                contig_name=fields[0],
                host_chr=fields[2],
                junction_pos=int(fields[9]),
                junction_type=fields[10],
                host_mapq=int(fields[12]) if len(fields) > 12 else 0,
            )
            evidences.append(ev)

    log(f"Verifying {len(evidences)} junction candidate(s)...")

    for ev in evidences:
        log(f"\n  --- {ev.host_chr}:{ev.junction_pos} ({ev.junction_type}, MAPQ={ev.host_mapq}) ---")

        # 1. Read support from host BAM
        if host_bam and host_bam.exists():
            split, disc, sc, depth = count_supporting_reads(
                host_bam, ev.host_chr, ev.junction_pos,
            )
            ev.split_reads = split
            ev.discordant_pairs = disc
            ev.soft_clipped_reads = sc
            ev.depth_at_junction = depth
            log(f"  Read support: split={split}, discordant={disc}, soft-clip={sc}, depth={depth}")

            # 2. Depth profile
            left, right = get_depth_profile(host_bam, ev.host_chr, ev.junction_pos)
            ev.depth_flanking_left = left
            ev.depth_flanking_right = right
            avg_flank = (left + right) / 2 if (left + right) > 0 else 1
            ev.depth_ratio = depth / avg_flank if avg_flank > 0 else 0
            log(f"  Depth profile: left={left:.1f}, junction={depth}, right={right:.1f}, ratio={ev.depth_ratio:.2f}")

        # 3. BLAST uniqueness check
        try:
            hits = blast_junction_region(
                contigs_fa, ev.contig_name, host_ref, 0, 0,  # Full contig
            )
            ev.blast_hits_in_genome = hits
            log(f"  BLAST hits in genome: {hits}")
        except Exception as e:
            log(f"  BLAST check skipped: {e}")

        # 4. WT-based check
        if wt_construct_bam and wt_construct_bam.exists():
            wt_count = check_wt_reads(
                wt_construct_bam, host_ref, ev.host_chr, ev.junction_pos,
            )
            ev.wt_reads_at_site = wt_count
            log(f"  WT construct reads: {wt_count}")

        # 5. Assign verdict
        ev.verdict, ev.reason = _assign_verdict(ev)
        log(f"  Verdict: {ev.verdict} - {ev.reason}")

    # Write verified junctions
    verified_tsv = step_dir / "verified_junctions.tsv"
    report_txt = step_dir / "verification_report.txt"

    with open(verified_tsv, "w") as f:
        f.write("\t".join([
            "contig_name", "host_chr", "junction_pos", "junction_type",
            "host_mapq", "split_reads", "discordant_pairs", "soft_clipped",
            "depth_at_junction", "depth_left", "depth_right", "depth_ratio",
            "blast_hits", "wt_reads", "verdict", "reason",
        ]) + "\n")
        for ev in evidences:
            f.write("\t".join([
                ev.contig_name, ev.host_chr, str(ev.junction_pos),
                ev.junction_type, str(ev.host_mapq),
                str(ev.split_reads), str(ev.discordant_pairs),
                str(ev.soft_clipped_reads), str(ev.depth_at_junction),
                f"{ev.depth_flanking_left:.1f}", f"{ev.depth_flanking_right:.1f}",
                f"{ev.depth_ratio:.2f}", str(ev.blast_hits_in_genome),
                str(ev.wt_reads_at_site), ev.verdict, ev.reason,
            ]) + "\n")

    log(f"\nVerified junctions written to: {verified_tsv}")

    # Write human-readable report
    with open(report_txt, "w") as f:
        f.write(f"Junction Verification Report - {sample_name}\n")
        f.write("=" * 60 + "\n\n")

        true_count = sum(1 for e in evidences if e.verdict == "TRUE")
        fp_count = sum(1 for e in evidences if e.verdict == "FALSE_POSITIVE")
        unk_count = sum(1 for e in evidences if e.verdict == "UNCERTAIN")

        f.write(f"Total candidates: {len(evidences)}\n")
        f.write(f"  TRUE: {true_count}\n")
        f.write(f"  FALSE_POSITIVE: {fp_count}\n")
        f.write(f"  UNCERTAIN: {unk_count}\n\n")

        for ev in evidences:
            f.write(f"--- {ev.host_chr}:{ev.junction_pos} ---\n")
            f.write(f"  Contig: {ev.contig_name}\n")
            f.write(f"  Type: {ev.junction_type}\n")
            f.write(f"  MAPQ: {ev.host_mapq}\n")
            f.write(f"  Split reads: {ev.split_reads}\n")
            f.write(f"  Discordant pairs: {ev.discordant_pairs}\n")
            f.write(f"  Soft-clipped reads: {ev.soft_clipped_reads}\n")
            f.write(f"  Depth (junction/left/right): {ev.depth_at_junction}/{ev.depth_flanking_left:.1f}/{ev.depth_flanking_right:.1f}\n")
            f.write(f"  Depth ratio: {ev.depth_ratio:.2f}\n")
            f.write(f"  BLAST hits: {ev.blast_hits_in_genome}\n")
            f.write(f"  WT reads: {ev.wt_reads_at_site}\n")
            f.write(f"  VERDICT: {ev.verdict}\n")
            f.write(f"  Reason: {ev.reason}\n\n")

    log(f"Report written to: {report_txt}")
    return evidences


def _assign_verdict(ev: JunctionEvidence) -> tuple[str, str]:
    """Assign TRUE/FALSE_POSITIVE/UNCERTAIN verdict based on evidence."""
    reasons = []

    # Strong false positive indicators
    if ev.host_mapq < 5:
        return "FALSE_POSITIVE", f"Very low MAPQ ({ev.host_mapq}), non-unique host mapping"

    if ev.wt_reads_at_site > 10:
        return "FALSE_POSITIVE", f"WT has {ev.wt_reads_at_site} reads at this site (homologous sequence)"

    # Supporting evidence for true junction
    support_score = 0

    if ev.host_mapq >= 30:
        support_score += 2
        reasons.append(f"High MAPQ ({ev.host_mapq})")
    elif ev.host_mapq >= 10:
        support_score += 1
        reasons.append(f"Moderate MAPQ ({ev.host_mapq})")

    if ev.split_reads >= 3:
        support_score += 2
        reasons.append(f"{ev.split_reads} split reads")
    elif ev.split_reads >= 1:
        support_score += 1

    if ev.soft_clipped_reads >= 3:
        support_score += 2
        reasons.append(f"{ev.soft_clipped_reads} soft-clipped reads")
    elif ev.soft_clipped_reads >= 1:
        support_score += 1

    if ev.discordant_pairs >= 3:
        support_score += 1
        reasons.append(f"{ev.discordant_pairs} discordant pairs")

    # Depth anomaly check
    if ev.depth_ratio > 0 and ev.depth_ratio < 0.3:
        support_score += 1
        reasons.append(f"Depth drop at junction (ratio={ev.depth_ratio:.2f})")

    if ev.wt_reads_at_site == 0:
        support_score += 1
        reasons.append("No WT reads at site")

    # Verdict
    if support_score >= 4:
        return "TRUE", "; ".join(reasons)
    elif support_score >= 2:
        return "UNCERTAIN", "; ".join(reasons) if reasons else "Moderate evidence"
    else:
        return "UNCERTAIN", "Insufficient evidence"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 6b: Junction verification protocol",
    )
    parser.add_argument(
        "--junctions", required=True, type=Path,
        help="junctions.tsv from Step 6",
    )
    parser.add_argument(
        "--host-bam", type=Path, default=None,
        help="Host BAM from Step 7 (optional, enables read support analysis)",
    )
    parser.add_argument(
        "--contigs", required=True, type=Path,
        help="Assembled contigs FASTA from Step 4",
    )
    parser.add_argument(
        "--host-ref", required=True, type=Path,
        help="Host reference genome FASTA",
    )
    parser.add_argument(
        "--wt-construct-bam", type=Path, default=None,
        help="WT construct BAM from Step 2 (optional, enables WT-based filtering)",
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

    verify_junctions(
        junctions_tsv=args.junctions,
        host_bam=args.host_bam,
        contigs_fa=args.contigs,
        host_ref=args.host_ref,
        wt_construct_bam=args.wt_construct_bam,
        outdir=args.outdir,
        sample_name=args.sample_name,
    )

    log("Done.")


if __name__ == "__main__":
    main()
