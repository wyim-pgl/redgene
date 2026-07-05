#!/usr/bin/env python3
"""Extract soft-clip metrics for known ground truth junction sites.

Scans the host BAM at each ground truth position for soft-clipped reads,
measures clip lengths, read counts, element DB identity, and host mapping.
Outputs ground_truth_baseline.tsv.
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import pysam

# Known ground truth sites (chr, position, sample, description)
GROUND_TRUTH = [
    # Rice G281
    ("Chr3", 16439674, "rice_G281", "T-DNA insertion, 2-copy head-to-head"),
    ("Chr3", 16439710, "rice_G281", "T-DNA insertion, second junction"),
    # Tomato Cas9 A2_3
    ("SLM_r2.0ch08", 65107378, "tomato_Cas9_A2_3", "T-DNA RB junction"),
    # Cucumber line 224
    ("LKUO03001512.1", 580628, "cucumber_line224", "T-DNA LB junction"),
    ("LKUO03001512.1", 581332, "cucumber_line224", "T-DNA RB junction"),
    # Cucumber line 212
    ("LKUO03001392.1", 2751693, "cucumber_line212", "T-DNA junction"),
    # Cucumber line 225
    ("LKUO03002166.1", 547987, "cucumber_line225", "T-DNA junction"),
    # Corn ND207
    ("NC_050098.1", 181367276, "corn_ND207", "T-DNA LB junction"),
    # Soybean UGT72E3
    ("NC_038254.2", 51882860, "soybean_UGT72E3", "T-DNA junction in Glyma.18g226800"),
    # Soybean AtYUCCA6 Site II
    ("NC_038255.2", 49789752, "soybean_AtYUCCA6", "T-DNA Site II junction"),
]


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def extract_clips_at_position(
    bam_path: Path, chrom: str, pos: int, window: int = 100, min_clip: int = 10,
) -> dict:
    """Extract soft-clip stats at a ground truth position."""
    right_clips: list[str] = []
    left_clips: list[str] = []
    total_reads = 0

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(chrom, max(0, pos - window), pos + window):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if not read.cigartuples or not read.query_sequence:
                    continue
                total_reads += 1
                cigar = read.cigartuples
                seq = read.query_sequence

                # Right soft-clip (5' junction seed)
                if cigar[-1][0] == 4 and cigar[-1][1] >= min_clip:
                    if read.reference_end and abs(read.reference_end - pos) <= window:
                        clip_len = cigar[-1][1]
                        right_clips.append(seq[-clip_len:])

                # Left soft-clip (3' junction seed)
                if cigar[0][0] == 4 and cigar[0][1] >= min_clip:
                    if read.reference_start and abs(read.reference_start - pos) <= window:
                        clip_len = cigar[0][1]
                        left_clips.append(seq[:clip_len])
    except ValueError:
        return {"found": False}

    if not right_clips and not left_clips:
        return {"found": False, "total_reads": total_reads}

    return {
        "found": True,
        "total_reads": total_reads,
        "right_clip_count": len(right_clips),
        "left_clip_count": len(left_clips),
        "right_clip_max_len": max((len(c) for c in right_clips), default=0),
        "left_clip_max_len": max((len(c) for c in left_clips), default=0),
        "right_clip_median_len": sorted(len(c) for c in right_clips)[len(right_clips)//2] if right_clips else 0,
        "left_clip_median_len": sorted(len(c) for c in left_clips)[len(left_clips)//2] if left_clips else 0,
        "right_consensus": max(right_clips, key=len) if right_clips else "",
        "left_consensus": max(left_clips, key=len) if left_clips else "",
    }


def check_element_hit(seq: str, element_db: Path, evalue: str = "1e-5") -> dict:
    """BLAST a clip sequence against element_db, return best hit."""
    if not seq or len(seq) < 15:
        return {"hit": False}
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
        fh.write(f">query\n{seq}\n")
        query_fa = Path(fh.name)
    out_path = query_fa.with_suffix(".blast")
    subprocess.run(
        ["blastn", "-query", str(query_fa), "-subject", str(element_db),
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
         "-evalue", evalue, "-max_target_seqs", "1",
         "-out", str(out_path)],
        stderr=subprocess.DEVNULL,
    )
    result = {"hit": False}
    if out_path.exists() and out_path.stat().st_size > 0:
        cols = open(out_path).readline().strip().split("\t")
        if len(cols) >= 6:
            result = {
                "hit": True,
                "element": cols[1],
                "identity": float(cols[2]),
                "aln_length": int(cols[3]),
                "evalue": cols[4],
                "bitscore": float(cols[5]),
            }
    query_fa.unlink(missing_ok=True)
    out_path.unlink(missing_ok=True)
    return result


def check_host_mapping(seq: str, host_ref: Path) -> bool:
    """Check if clip sequence maps to host genome."""
    if not seq or len(seq) < 15:
        return False
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
        fh.write(f">query\n{seq}\n")
        query_fa = Path(fh.name)
    paf_path = query_fa.with_suffix(".paf")
    subprocess.run(
        ["minimap2", "-c", "--secondary=no", str(host_ref), str(query_fa)],
        stdout=open(paf_path, "w"), stderr=subprocess.DEVNULL,
    )
    maps = False
    if paf_path.exists():
        for line in open(paf_path):
            cols = line.strip().split("\t")
            if len(cols) >= 12:
                match_bp = int(cols[9])
                if match_bp >= 30:
                    maps = True
                    break
    query_fa.unlink(missing_ok=True)
    paf_path.unlink(missing_ok=True)
    return maps


def main():
    parser = argparse.ArgumentParser(description="Extract ground truth clip baselines")
    parser.add_argument("--host-bam", required=True)
    parser.add_argument("--host-ref", required=True)
    parser.add_argument("--element-db", required=True)
    parser.add_argument("--sample", required=True,
                        help="Sample name (must match GROUND_TRUTH entries)")
    parser.add_argument("--output", default="ground_truth_baseline.tsv")
    parser.add_argument("--min-clip", type=int, default=10,
                        help="Min clip for scanning (use low value to catch everything)")
    args = parser.parse_args()

    host_bam = Path(args.host_bam)
    host_ref = Path(args.host_ref)
    element_db = Path(args.element_db)

    sample_gt = [gt for gt in GROUND_TRUTH if gt[2] == args.sample]
    if not sample_gt:
        print(f"No ground truth entries for sample '{args.sample}'", file=sys.stderr)
        print(f"Available: {sorted(set(g[2] for g in GROUND_TRUTH))}", file=sys.stderr)
        sys.exit(1)

    write_header = not Path(args.output).exists()
    with open(args.output, "a") as fout:
        if write_header:
            fout.write("sample\tchrom\tposition\tdescription\t"
                       "right_clip_count\tleft_clip_count\t"
                       "right_clip_max_len\tleft_clip_max_len\t"
                       "right_clip_median_len\tleft_clip_median_len\t"
                       "element_hit_5p\telement_id_5p\telement_identity_5p\telement_aln_len_5p\t"
                       "element_hit_3p\telement_id_3p\telement_identity_3p\telement_aln_len_3p\t"
                       "maps_to_host_5p\tmaps_to_host_3p\n")

        for chrom, pos, sample, desc in sample_gt:
            print(f"Checking {sample} {chrom}:{pos} ({desc})...", file=sys.stderr)
            clips = extract_clips_at_position(
                host_bam, chrom, pos, min_clip=args.min_clip)

            if not clips.get("found"):
                print(f"  WARNING: No clips found at {chrom}:{pos}", file=sys.stderr)
                fout.write(f"{sample}\t{chrom}\t{pos}\t{desc}\t"
                           f"0\t0\t0\t0\t0\t0\t"
                           f"False\t-\t0\t0\tFalse\t-\t0\t0\t"
                           f"False\tFalse\n")
                continue

            elem_5p = check_element_hit(clips.get("right_consensus", ""), element_db)
            elem_3p = check_element_hit(clips.get("left_consensus", ""), element_db)
            host_5p = check_host_mapping(clips.get("right_consensus", ""), host_ref)
            host_3p = check_host_mapping(clips.get("left_consensus", ""), host_ref)

            print(f"  5p: {clips['right_clip_count']} clips, "
                  f"max {clips['right_clip_max_len']}bp, "
                  f"element={'YES' if elem_5p['hit'] else 'no'}, "
                  f"host={'YES' if host_5p else 'no'}", file=sys.stderr)
            print(f"  3p: {clips['left_clip_count']} clips, "
                  f"max {clips['left_clip_max_len']}bp, "
                  f"element={'YES' if elem_3p['hit'] else 'no'}, "
                  f"host={'YES' if host_3p else 'no'}", file=sys.stderr)

            fout.write(
                f"{sample}\t{chrom}\t{pos}\t{desc}\t"
                f"{clips['right_clip_count']}\t{clips['left_clip_count']}\t"
                f"{clips['right_clip_max_len']}\t{clips['left_clip_max_len']}\t"
                f"{clips['right_clip_median_len']}\t{clips['left_clip_median_len']}\t"
                f"{elem_5p.get('hit', False)}\t{elem_5p.get('element', '-')}\t"
                f"{elem_5p.get('identity', 0)}\t{elem_5p.get('aln_length', 0)}\t"
                f"{elem_3p.get('hit', False)}\t{elem_3p.get('element', '-')}\t"
                f"{elem_3p.get('identity', 0)}\t{elem_3p.get('aln_length', 0)}\t"
                f"{host_5p}\t{host_3p}\n"
            )

    print(f"Baseline written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
