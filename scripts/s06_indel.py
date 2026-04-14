#!/usr/bin/env python3
"""Step 6: CRISPR editing site detection via indel calling.

Two modes of operation:

1. **With gRNA sequences** (recommended):
   - Map gRNA to host genome to find on-target and off-target sites
   - Call variants (indels + SNVs) at predicted sites in treatment vs WT
   - High specificity, fast execution

2. **Without gRNA** (de novo screening):
   - Genome-wide or region-based indel calling
   - Filter for treatment-specific indels with PAM + NHEJ signature
   - Lower specificity, slower, more false positives

For non-CRISPR samples (standard T-DNA), this step is skipped.

Usage:
  # With gRNA (recommended)
  python s06_indel.py \
    --treatment-bam results/{sample}/s04_host_map/{sample}_host.bam \
    --wt-bam results/WT/s04_host_map/WT_host.bam \
    --host-ref db/host.fa \
    --grna ATCGATCGATCGATCGATCG,GCTAGCTAGCTAGCTAGCTA \
    --outdir results --sample-name {sample}

  # Without gRNA (de novo)
  python s06_indel.py \
    --treatment-bam results/{sample}/s04_host_map/{sample}_host.bam \
    --wt-bam results/WT/s04_host_map/WT_host.bam \
    --host-ref db/host.fa --genome-wide \
    --outdir results --sample-name {sample}
"""

import argparse
from collections import Counter
import re
import subprocess
import sys
from pathlib import Path


def log(msg: str) -> None:
    print(f"[s06_indel] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# gRNA target site prediction
# ---------------------------------------------------------------------------

def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


def find_grna_targets(
    grna_seqs: list[str],
    host_ref: Path,
    max_mismatches: int = 4,
) -> list[dict]:
    """Find on-target and off-target sites for gRNA sequences using BLAST.

    Uses blastn-short to find sites with up to max_mismatches mismatches.
    BLAST is used instead of minimap2 because 20bp gRNA sequences are too
    short for minimap2's default seed parameters.
    """
    import tempfile

    # Ensure BLAST db exists
    nhr = host_ref.parent / f"{host_ref.name}.nhr"
    if not nhr.exists():
        log("Creating BLAST database for host genome...")
        subprocess.run(
            ["makeblastdb", "-in", str(host_ref), "-dbtype", "nucl"],
            check=True, capture_output=True,
        )

    targets = []

    for i, grna in enumerate(grna_seqs):
        grna = grna.upper().strip()
        if len(grna) < 17 or len(grna) > 25:
            log(f"WARNING: gRNA {i+1} length {len(grna)} unusual (expected 17-25bp)")

        # Create temp FASTA with gRNA sequence
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(f">gRNA_{i+1}\n{grna}\n")
            query_fa = f.name

        # Use blastn-short for sensitive mapping of short sequences
        result = subprocess.run(
            ["blastn", "-query", query_fa,
             "-db", str(host_ref),
             "-task", "blastn-short",
             "-evalue", "100",
             "-word_size", "7",
             "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"],
            capture_output=True, text=True,
        )

        Path(query_fa).unlink(missing_ok=True)

        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue

            pident = float(fields[2])
            aln_len = int(fields[3])
            n_mismatch = int(fields[4])
            q_start = int(fields[6])
            q_end = int(fields[7])
            s_start = int(fields[8])
            s_end = int(fields[9])
            t_name = fields[1]

            # Only consider near-full-length alignments
            if aln_len < len(grna) - max_mismatches:
                continue
            if n_mismatch > max_mismatches:
                continue

            # Determine strand
            if s_start < s_end:
                strand = "+"
                t_start = s_start - 1  # Convert to 0-based
                t_end = s_end
            else:
                strand = "-"
                t_start = s_end - 1
                t_end = s_start

            # Check for PAM (NGG) at the expected position
            # For + strand: PAM is immediately after the gRNA (3' end)
            # For - strand: PAM (CCN on reference) is immediately before
            pam_start = t_end if strand == "+" else t_start - 3
            pam_end = pam_start + 3

            pam_check = subprocess.run(
                ["samtools", "faidx", str(host_ref),
                 f"{t_name}:{max(1, pam_start + 1)}-{pam_end}"],
                capture_output=True, text=True,
            )
            pam_seq = "".join(pam_check.stdout.strip().split("\n")[1:]).upper()

            has_pam = False
            if strand == "+" and len(pam_seq) >= 3 and pam_seq[1:3] == "GG":
                has_pam = True
            elif strand == "-" and len(pam_seq) >= 3 and pam_seq[0:2] == "CC":
                has_pam = True

            mismatches = n_mismatch + (len(grna) - aln_len)  # count unaligned bases
            site_type = "on-target" if mismatches == 0 and has_pam else "off-target"
            if mismatches == 0 and not has_pam:
                site_type = "on-target-no-PAM"

            # Cas9 cut site: 3bp upstream of PAM
            if strand == "+":
                cut_pos = t_end - 3
            else:
                cut_pos = t_start + 3

            targets.append({
                "grna_idx": i + 1,
                "grna_seq": grna,
                "chrom": t_name,
                "start": t_start,
                "end": t_end,
                "strand": strand,
                "mismatches": mismatches,
                "mapq": 60,  # BLAST doesn't report MAPQ; use 60 for exact matches
                "has_pam": has_pam,
                "pam_seq": pam_seq,
                "site_type": site_type,
                "cut_pos": cut_pos,
            })

    # Sort: on-target first, then by mismatches
    targets.sort(key=lambda x: (x["mismatches"], x["site_type"] != "on-target"))

    log(f"Found {len(targets)} gRNA target sites:")
    on_target = [t for t in targets if t["site_type"] == "on-target"]
    off_target = [t for t in targets if t["site_type"] == "off-target"]
    log(f"  On-target (0 mismatch + PAM): {len(on_target)}")
    log(f"  Off-target (1-{max_mismatches} mismatches or no PAM): {len(off_target)}")

    return targets


def call_variants_at_sites(
    treatment_bam: Path,
    wt_bam: Path,
    host_ref: Path,
    targets: list[dict],
    window: int = 50,
    min_indel_size: int = 1,
) -> list[dict]:
    """Call variants at predicted gRNA target sites using direct pileup parsing.

    Uses samtools mpileup directly instead of bcftools call to avoid
    variant decomposition issues at low depth. Parses indel strings
    (+N/-N) from pileup to accurately capture CRISPR editing events.
    """
    results = []

    for target in targets:
        chrom = target["chrom"]
        cut = target["cut_pos"]
        region = f"{chrom}:{max(1,cut-window)}-{cut+window}"

        # Parse pileup for treatment and WT
        t_indels = _parse_pileup_indels(treatment_bam, host_ref, region)
        wt_indels = _parse_pileup_indels(wt_bam, host_ref, region)

        # Build WT set for filtering (key = chrom:pos:type:seq)
        wt_set = set()
        for indel in wt_indels:
            key = f"{indel['chrom']}:{indel['pos']}:{indel['type']}:{indel['indel_seq']}"
            wt_set.add(key)

        for indel in t_indels:
            key = f"{indel['chrom']}:{indel['pos']}:{indel['type']}:{indel['indel_seq']}"
            if key in wt_set:
                continue  # Present in WT, skip

            if indel["size"] < min_indel_size:
                continue

            indel["grna_idx"] = target["grna_idx"]
            indel["grna_seq"] = target["grna_seq"]
            indel["site_type"] = target["site_type"]
            indel["mismatches"] = target["mismatches"]
            indel["distance_to_cut"] = abs(indel["pos"] - cut)

            results.append(indel)

    return results


def _parse_pileup_indels(
    bam: Path, host_ref: Path, region: str, min_freq: float = 0.1,
) -> list[dict]:
    """Parse indels directly from samtools mpileup output.

    This avoids bcftools variant decomposition that can split complex
    indels into multiple smaller variants at low depth.
    """
    # Use low base quality threshold for gRNA-guided mode to avoid
    # filtering out indel-supporting reads with borderline quality
    cmd = [
        "samtools", "mpileup", "-f", str(host_ref),
        "-r", region, "-d", "500", "-q", "10", "-Q", "0",
        str(bam),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    indels = []  # list of {chrom, pos, type, size, indel_seq, dp, count, freq, ref}
    seen = {}  # deduplicate: (chrom, pos, type, seq) -> index

    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue

        chrom = fields[0]
        pos = int(fields[1])
        ref_base = fields[2].upper()
        dp = int(fields[3])
        pileup = fields[4]

        if dp == 0:
            continue

        # Parse indels from pileup string
        pos_indels = _extract_indels_from_pileup(pileup)

        # Aggregate counts case-insensitively
        indel_counts = Counter(
            (t, s.upper()) for t, s in pos_indels
        )

        for (indel_type, indel_seq_upper), count in indel_counts.items():
            freq = count / dp

            if freq < min_freq:
                continue

            key = (chrom, pos, indel_type, indel_seq_upper)
            if key in seen:
                continue
            seen[key] = True

            size = len(indel_seq)

            # Determine zygosity from allele frequency
            if freq >= 0.85:
                gt = "1/1"
                zygosity = "homozygous"
            elif freq >= 0.2:
                gt = "0/1"
                zygosity = "heterozygous"
            else:
                gt = "0/0"
                zygosity = "low-frequency"

            # Build ref/alt in VCF-like format
            if indel_type == "deletion":
                ref_allele = ref_base + indel_seq_upper
                alt_allele = ref_base
            else:  # insertion
                ref_allele = ref_base
                alt_allele = ref_base + indel_seq_upper

            indels.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref_allele,
                "alt": alt_allele,
                "type": indel_type,
                "size": size,
                "indel_seq": indel_seq_upper,
                "qual": min(freq * dp * 10, 999),  # Approximate quality
                "dp": dp,
                "count": count,
                "freq": freq,
                "gt": gt,
                "zygosity": zygosity,
            })

    return indels


def _extract_indels_from_pileup(pileup: str) -> list[tuple[str, str]]:
    """Extract all indels from a pileup string.

    Returns list of (type, sequence) tuples where type is 'insertion'
    or 'deletion' and sequence is the indel bases.
    """
    indels = []
    i = 0
    while i < len(pileup):
        c = pileup[i]
        if c == "^":
            i += 2  # skip ^ and mapping quality char
            continue
        elif c == "$":
            i += 1
            continue
        elif c in "+-":
            indel_type = "insertion" if c == "+" else "deletion"
            i += 1
            # Read the length digits
            num_str = ""
            while i < len(pileup) and pileup[i].isdigit():
                num_str += pileup[i]
                i += 1
            if num_str:
                n = int(num_str)
                seq = pileup[i:i+n]
                i += n
                indels.append((indel_type, seq))
            continue
        elif c == "*":
            # Deleted base placeholder
            i += 1
            continue
        else:
            i += 1
            continue
    return indels


# ---------------------------------------------------------------------------
# De novo mode (no gRNA)
# ---------------------------------------------------------------------------

def call_indels(bam: Path, host_ref: Path, region: str | None,
                min_indel_size: int = 2) -> list[dict]:
    """Call indels from a BAM file using bcftools."""
    cmd = ["bcftools", "mpileup", "-f", str(host_ref)]
    if region:
        cmd.extend(["-r", region])
    cmd.extend(["-d", "500", "-q", "10", "-Q", "20", str(bam)])

    mpileup = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    call_cmd = ["bcftools", "call", "-mv", "--ploidy", "2", "-Ov", "--threads", "4"]
    call_proc = subprocess.Popen(call_cmd, stdin=mpileup.stdout,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  text=True)
    mpileup.stdout.close()
    stdout, _ = call_proc.communicate()
    mpileup.wait()

    indels = []
    for line in stdout.strip().split("\n"):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue

        ref_allele = fields[3]
        alt_allele = fields[4]
        if len(ref_allele) == len(alt_allele):
            continue  # SNP

        indel_size = abs(len(alt_allele) - len(ref_allele))
        if indel_size < min_indel_size:
            continue

        info = fields[7]
        dp = 0
        if "DP=" in info:
            dp = int(info.split("DP=")[1].split(";")[0])

        gt_data = fields[9] if len(fields) > 9 else ""

        indels.append({
            "chrom": fields[0],
            "pos": int(fields[1]),
            "ref": ref_allele,
            "alt": alt_allele,
            "qual": float(fields[5]) if fields[5] != "." else 0,
            "type": "deletion" if len(ref_allele) > len(alt_allele) else "insertion",
            "size": indel_size,
            "dp": dp,
            "gt": gt_data.split(":")[0] if gt_data else ".",
        })

    return indels


def find_treatment_specific_indels(
    treatment_bam: Path, wt_bam: Path, host_ref: Path,
    regions: list[str | None], min_indel_size: int = 2,
) -> list[dict]:
    """Find indels present in treatment but absent in WT."""
    treatment_indels = []
    wt_set = set()

    for region in regions:
        log(f"Scanning region: {region or 'genome-wide'}...")
        t_indels = call_indels(treatment_bam, host_ref, region, min_indel_size)
        wt_indels = call_indels(wt_bam, host_ref, region, min_indel_size)
        log(f"  Treatment: {len(t_indels)}, WT: {len(wt_indels)}")

        for w in wt_indels:
            wt_set.add(f"{w['chrom']}:{w['pos']}:{w['ref']}:{w['alt']}")
        treatment_indels.extend(t_indels)

    specific = [i for i in treatment_indels
                if f"{i['chrom']}:{i['pos']}:{i['ref']}:{i['alt']}" not in wt_set]
    log(f"Treatment-specific indels: {len(specific)}")
    return specific


def annotate_pam_and_nhej(indels: list[dict], host_ref: Path) -> None:
    """Annotate indels with PAM proximity and NHEJ signature."""
    for indel in indels:
        indel["has_pam"] = False
        indel["pam_motif"] = ""
        indel["nhej_signature"] = False
        indel["crispr_candidate"] = False

        # PAM check at specific positions (Cas9 cuts 3bp upstream of NGG)
        try:
            region = f"{indel['chrom']}:{max(1,indel['pos']-8)}-{indel['pos']+8}"
            result = subprocess.run(
                ["samtools", "faidx", str(host_ref), region],
                capture_output=True, text=True,
            )
            seq = "".join(result.stdout.strip().split("\n")[1:]).upper()
            mid = 8
            if mid + 6 <= len(seq) and seq[mid+4:mid+6] == "GG":
                indel["has_pam"] = True
                indel["pam_motif"] = seq[mid+3:mid+6] + "(+)"
            elif mid >= 6 and seq[mid-6:mid-4] == "CC":
                indel["has_pam"] = True
                indel["pam_motif"] = seq[mid-6:mid-3] + "(-)"
        except Exception:
            pass

        # NHEJ signature
        if indel["type"] == "deletion" and 1 <= indel["size"] <= 20:
            indel["nhej_signature"] = True
        elif indel["type"] == "insertion" and 1 <= indel["size"] <= 5:
            indel["nhej_signature"] = True

        if indel["has_pam"] and indel["nhej_signature"] and indel["qual"] >= 20:
            indel["crispr_candidate"] = True

        # Zygosity
        gt = indel.get("gt", ".")
        if gt in ("1/1", "1|1"):
            indel["zygosity"] = "homozygous"
        elif gt in ("0/1", "0|1", "1/0", "1|0"):
            indel["zygosity"] = "heterozygous"
        elif "/" in gt or "|" in gt:
            alleles = gt.replace("|", "/").split("/")
            if len(set(alleles)) == 1 and alleles[0] != "0":
                indel["zygosity"] = "homozygous"
            elif "0" not in alleles:
                indel["zygosity"] = "biallelic"
            else:
                indel["zygosity"] = "heterozygous"
        else:
            indel["zygosity"] = "unknown"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_grna_input(grna_arg: str) -> list[str]:
    """Parse gRNA input: comma-separated sequences or file path."""
    p = Path(grna_arg)
    if p.exists():
        seqs = []
        with open(p) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith(">"):
                    continue
                seqs.append(line.upper())
        return seqs
    return [s.strip().upper() for s in grna_arg.split(",") if s.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 6: CRISPR editing site detection",
    )
    parser.add_argument("--treatment-bam", required=True, type=Path)
    parser.add_argument("--wt-bam", required=True, type=Path)
    parser.add_argument("--host-ref", required=True, type=Path)
    parser.add_argument(
        "--grna", type=str, default=None,
        help="gRNA sequences (comma-separated or file path). "
             "If provided, uses gRNA-guided mode for higher specificity.",
    )
    parser.add_argument("--junctions", type=Path, default=None)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--min-indel-size", type=int, default=2)
    parser.add_argument("--search-window", type=int, default=50000)
    parser.add_argument("--genome-wide", action="store_true")
    parser.add_argument(
        "--max-off-target-mismatches", type=int, default=4,
        help="Max mismatches for off-target prediction (default: 4)",
    )
    args = parser.parse_args()

    for label, path in [("Treatment BAM", args.treatment_bam),
                        ("WT BAM", args.wt_bam),
                        ("Host ref", args.host_ref)]:
        if not path.exists():
            log(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    step_dir = args.outdir / args.sample_name / "s06_indel"
    step_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # MODE 1: gRNA-guided analysis
    # -----------------------------------------------------------------------
    if args.grna:
        grna_seqs = parse_grna_input(args.grna)
        log(f"gRNA-guided mode: {len(grna_seqs)} gRNA sequence(s)")
        for i, seq in enumerate(grna_seqs, 1):
            log(f"  gRNA {i}: {seq} ({len(seq)}bp)")

        # Find target sites
        targets = find_grna_targets(
            grna_seqs, args.host_ref,
            max_mismatches=args.max_off_target_mismatches,
        )

        # Write target sites
        targets_tsv = step_dir / "grna_targets.tsv"
        with open(targets_tsv, "w") as f:
            header = ["grna_idx", "grna_seq", "chrom", "start", "end", "strand",
                      "mismatches", "has_pam", "pam_seq", "site_type", "cut_pos"]
            f.write("\t".join(header) + "\n")
            for t in targets:
                f.write("\t".join(str(t.get(k, "")) for k in header) + "\n")
        log(f"Target sites written to: {targets_tsv}")

        # Call variants at target sites
        log("Calling variants at predicted target sites...")
        variants = call_variants_at_sites(
            args.treatment_bam, args.wt_bam, args.host_ref,
            targets, window=50, min_indel_size=1,  # Include 1bp indels for gRNA mode
        )

        # Write results
        output_tsv = step_dir / "editing_sites.tsv"
        with open(output_tsv, "w") as f:
            header = ["chrom", "pos", "ref", "alt", "type", "size",
                      "indel_seq", "qual", "dp", "count", "freq",
                      "gt", "zygosity", "grna_idx", "site_type",
                      "mismatches", "distance_to_cut"]
            f.write("\t".join(header) + "\n")
            for v in variants:
                f.write("\t".join(str(v.get(k, "")) for k in header) + "\n")

        # Summary
        on_target_edits = [v for v in variants if v["site_type"] == "on-target"]
        off_target_edits = [v for v in variants if v["site_type"] == "off-target"]

        log(f"\n=== gRNA-Guided Editing Summary for {args.sample_name} ===")
        log(f"Total treatment-specific variants at target sites: {len(variants)}")
        log(f"  On-target edits: {len(on_target_edits)}")
        log(f"  Off-target edits: {len(off_target_edits)}")

        if on_target_edits:
            log("\nOn-target editing events:")
            for v in on_target_edits:
                log(f"  gRNA {v['grna_idx']}: {v['chrom']}:{v['pos']} "
                    f"{v['type']} {v['size']}bp ({v.get('indel_seq', '')}) "
                    f"freq={v.get('freq', 0):.1%} ({v.get('count', 0)}/{v['dp']}), "
                    f"{v['zygosity']}, dist={v['distance_to_cut']}bp from cut")

        if off_target_edits:
            log(f"\nOff-target edits ({len(off_target_edits)} total):")
            for v in off_target_edits[:10]:
                log(f"  gRNA {v['grna_idx']}: {v['chrom']}:{v['pos']} "
                    f"{v['type']} {v['size']}bp "
                    f"({v['mismatches']} mismatches, "
                    f"dist={v['distance_to_cut']}bp)")

        log(f"\nOutput: {output_tsv}")
        log("Done.")
        return

    # -----------------------------------------------------------------------
    # MODE 2: De novo screening (no gRNA)
    # -----------------------------------------------------------------------
    log("De novo mode (no gRNA provided)")
    log("WARNING: Without gRNA, specificity is limited. Consider providing "
        "--grna sequences for more accurate results.")

    # Parse junction regions
    regions: list[str | None] = []
    junction_positions = {}
    if args.junctions and args.junctions.exists():
        with open(args.junctions) as f:
            header = f.readline().strip().split("\t")
            for line in f:
                if not line.strip():
                    continue
                fields = line.strip().split("\t")
                row = dict(zip(header, fields))
                chrom = row["host_chr"]
                pos = int(row["junction_pos_host"])
                start = max(0, pos - args.search_window)
                end = pos + args.search_window
                region = f"{chrom}:{start}-{end}"
                if region not in regions:
                    regions.append(region)
                junction_positions[f"{chrom}:{pos}"] = row

    if not regions or args.genome_wide:
        log("Scanning genome-wide...")
        regions = [None]

    specific_indels = find_treatment_specific_indels(
        args.treatment_bam, args.wt_bam, args.host_ref,
        regions=regions, min_indel_size=args.min_indel_size,
    )

    specific_indels.sort(key=lambda x: (x["chrom"], x["pos"]))

    # Annotate
    annotate_pam_and_nhej(specific_indels, args.host_ref)

    # Junction distance
    for indel in specific_indels:
        min_dist = float("inf")
        nearest = ""
        for jkey in junction_positions:
            jchrom, jpos = jkey.rsplit(":", 1)
            if indel["chrom"] == jchrom:
                d = abs(indel["pos"] - int(jpos))
                if d < min_dist:
                    min_dist = d
                    nearest = jkey
        indel["nearest_junction"] = nearest
        indel["distance_to_junction"] = min_dist if min_dist < float("inf") else -1

    # Write output
    output_tsv = step_dir / "editing_sites.tsv"
    with open(output_tsv, "w") as f:
        header = ["chrom", "pos", "ref", "alt", "type", "size", "qual",
                  "dp", "gt", "zygosity", "has_pam", "pam_motif",
                  "nhej_signature", "crispr_candidate",
                  "nearest_junction", "distance_to_junction"]
        f.write("\t".join(header) + "\n")
        for indel in specific_indels:
            f.write("\t".join(str(indel.get(k, "")) for k in header) + "\n")

    candidates = [i for i in specific_indels if i.get("crispr_candidate")]

    log(f"\n=== De Novo Editing Summary for {args.sample_name} ===")
    log(f"Treatment-specific indels (>= {args.min_indel_size}bp): {len(specific_indels)}")
    log(f"CRISPR candidates (PAM + NHEJ + QUAL>=20): {len(candidates)}")

    if candidates:
        log("\nTop candidates:")
        candidates.sort(key=lambda x: -x["qual"])
        for indel in candidates[:20]:
            log(f"  {indel['chrom']}:{indel['pos']} "
                f"{indel['type']} {indel['size']}bp "
                f"(QUAL={indel['qual']:.0f}, {indel['zygosity']}, "
                f"PAM={indel['pam_motif']})")

    log(f"\nOutput: {output_tsv}")
    log("Done.")


if __name__ == "__main__":
    main()
