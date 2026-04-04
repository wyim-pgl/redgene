#!/usr/bin/env python3
"""Step 8: CRISPR editing site detection via indel calling.

For Cas9/CRISPR-edited samples, detects editing-induced indels by comparing
treatment vs wild-type BAMs. Finds indels (>=2bp) present in treatment but
absent in WT, focusing on regions near detected T-DNA insertion sites.

For non-CRISPR samples (standard T-DNA), this step is skipped.

Approach:
1. Run bcftools mpileup + call on treatment and WT BAMs
2. Filter for indels >= min_indel_size (default: 2bp)
3. Subtract WT variants to find treatment-specific indels
4. Report indels near T-DNA insertion sites as potential editing events

Input: Host BAMs (step 7), junctions.tsv (step 6)
Output: editing_sites.tsv with indel calls

Usage:
  python s08_indel_detection.py \
    --treatment-bam results/{sample}/s07_host_map/{sample}_host.bam \
    --wt-bam results/tomato_WT/s07_host_map/tomato_WT_host.bam \
    --host-ref db/host.fa \
    --junctions results/{sample}/s06_junction/junctions.tsv \
    --outdir results --sample-name {sample}
"""

import argparse
import subprocess
import sys
from pathlib import Path


def log(msg: str) -> None:
    print(f"[s08_indel] {msg}", file=sys.stderr, flush=True)


def call_indels(bam: Path, host_ref: Path, region: str | None,
                min_indel_size: int = 2) -> list[dict]:
    """Call indels from a BAM file using bcftools."""
    cmd = ["bcftools", "mpileup", "-f", str(host_ref)]
    if region:
        cmd.extend(["-r", region])
    cmd.extend(["-d", "500", "-q", "10", "-Q", "20", str(bam)])

    mpileup = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    call_cmd = ["bcftools", "call", "-mv", "--ploidy", "2",
                "-Ov", "--threads", "4"]
    call_proc = subprocess.Popen(call_cmd, stdin=mpileup.stdout,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  text=True)
    mpileup.stdout.close()
    stdout, stderr = call_proc.communicate()
    mpileup.wait()

    indels = []
    for line in stdout.strip().split("\n"):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue

        chrom = fields[0]
        pos = int(fields[1])
        ref_allele = fields[3]
        alt_allele = fields[4]
        qual = float(fields[5]) if fields[5] != "." else 0
        info = fields[7]

        # Filter: indels only (ref and alt differ in length)
        if len(ref_allele) == len(alt_allele):
            continue  # SNP, not indel

        indel_size = abs(len(alt_allele) - len(ref_allele))
        if indel_size < min_indel_size:
            continue

        indel_type = "deletion" if len(ref_allele) > len(alt_allele) else "insertion"

        # Parse genotype info
        gt_format = fields[8] if len(fields) > 8 else ""
        gt_data = fields[9] if len(fields) > 9 else ""

        # Get depth from INFO or FORMAT
        dp = 0
        if "DP=" in info:
            dp = int(info.split("DP=")[1].split(";")[0])

        indels.append({
            "chrom": chrom,
            "pos": pos,
            "ref": ref_allele,
            "alt": alt_allele,
            "qual": qual,
            "type": indel_type,
            "size": indel_size,
            "dp": dp,
            "gt": gt_data.split(":")[0] if gt_data else ".",
        })

    return indels


def find_treatment_specific_indels(
    treatment_bam: Path,
    wt_bam: Path,
    host_ref: Path,
    regions: list[str] | None = None,
    min_indel_size: int = 2,
    search_window: int = 50000,
) -> list[dict]:
    """Find indels present in treatment but absent in WT.

    Args:
        treatment_bam: Treatment sample BAM
        wt_bam: Wild-type BAM
        host_ref: Host reference FASTA
        regions: List of "chr:start-end" regions to focus on (from junctions)
        min_indel_size: Minimum indel size in bp (default: 2)
        search_window: Window around junction to search (default: 50kb)
    """
    treatment_indels = []
    wt_indels_set = set()

    # If we have junction regions, focus on those + genome-wide scan
    scan_regions = regions if regions else [None]

    for region in scan_regions:
        log(f"Scanning region: {region or 'genome-wide'}...")

        # Call indels in treatment
        t_indels = call_indels(treatment_bam, host_ref, region, min_indel_size)
        log(f"  Treatment: {len(t_indels)} indels >= {min_indel_size}bp")

        # Call indels in WT
        wt_indels = call_indels(wt_bam, host_ref, region, min_indel_size)
        log(f"  WT: {len(wt_indels)} indels >= {min_indel_size}bp")

        # Build WT set for subtraction (chr:pos:ref:alt)
        for w in wt_indels:
            wt_indels_set.add(f"{w['chrom']}:{w['pos']}:{w['ref']}:{w['alt']}")

        treatment_indels.extend(t_indels)

    # Filter: keep only treatment-specific
    specific = []
    for indel in treatment_indels:
        key = f"{indel['chrom']}:{indel['pos']}:{indel['ref']}:{indel['alt']}"
        if key not in wt_indels_set:
            specific.append(indel)

    log(f"Treatment-specific indels: {len(specific)}")
    return specific


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 8: CRISPR editing site detection via indel calling",
    )
    parser.add_argument(
        "--treatment-bam", required=True, type=Path,
        help="Treatment sample host BAM (from step 7)",
    )
    parser.add_argument(
        "--wt-bam", required=True, type=Path,
        help="Wild-type host BAM (from step 7)",
    )
    parser.add_argument(
        "--host-ref", required=True, type=Path,
        help="Host reference genome FASTA",
    )
    parser.add_argument(
        "--junctions", type=Path, default=None,
        help="junctions.tsv from step 6 (optional, to focus search)",
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
        "--min-indel-size", type=int, default=2,
        help="Minimum indel size in bp (default: 2)",
    )
    parser.add_argument(
        "--search-window", type=int, default=50000,
        help="Window around junction to search for indels (default: 50000)",
    )
    parser.add_argument(
        "--genome-wide", action="store_true",
        help="Scan entire genome (slow for large genomes)",
    )
    args = parser.parse_args()

    # Validate inputs
    for label, path in [("Treatment BAM", args.treatment_bam),
                        ("WT BAM", args.wt_bam),
                        ("Host ref", args.host_ref)]:
        if not path.exists():
            log(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    step_dir = args.outdir / args.sample_name / "s08_indel"
    step_dir.mkdir(parents=True, exist_ok=True)

    # Parse junction regions if provided
    regions = []
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
        log(f"Loaded {len(regions)} junction region(s)")

    # Default: genome-wide scan for CRISPR editing sites
    # CRISPR edits can be anywhere (target gene, NOT near T-DNA insertion)
    if not regions or args.genome_wide:
        log("Scanning genome-wide for treatment-specific indels...")
        regions = [None]

    # Find treatment-specific indels
    specific_indels = find_treatment_specific_indels(
        args.treatment_bam,
        args.wt_bam,
        args.host_ref,
        regions=regions,
        min_indel_size=args.min_indel_size,
        search_window=args.search_window,
    )

    # Sort by chromosome and position
    specific_indels.sort(key=lambda x: (x["chrom"], x["pos"]))

    # Annotate distance to nearest junction
    for indel in specific_indels:
        min_dist = float("inf")
        nearest_junction = ""
        for jkey, jrow in junction_positions.items():
            jchrom, jpos = jkey.rsplit(":", 1)
            if indel["chrom"] == jchrom:
                dist = abs(indel["pos"] - int(jpos))
                if dist < min_dist:
                    min_dist = dist
                    nearest_junction = jkey
        indel["nearest_junction"] = nearest_junction
        indel["distance_to_junction"] = min_dist if min_dist < float("inf") else -1

    # PAM motif and NHEJ signature analysis for CRISPR candidate filtering
    log("Checking PAM motifs and NHEJ signatures...")
    for indel in specific_indels:
        indel["has_pam"] = False
        indel["pam_motif"] = ""
        indel["nhej_signature"] = False
        indel["crispr_candidate"] = False

        # Check PAM (NGG) within 3-8bp downstream of indel
        # Extract flanking sequence from reference
        try:
            region = f"{indel['chrom']}:{max(1,indel['pos']-10)}-{indel['pos']+30}"
            result = subprocess.run(
                ["samtools", "faidx", str(args.host_ref), region],
                capture_output=True, text=True,
            )
            seq = "".join(result.stdout.strip().split("\n")[1:]).upper()
            # Look for NGG in the flanking region (both strands)
            # Forward strand PAM: NGG at position +3 to +8 from cut
            # Reverse strand PAM: CCN at position -8 to -3 from cut
            for i in range(len(seq) - 2):
                if seq[i+1:i+3] == "GG":
                    indel["has_pam"] = True
                    indel["pam_motif"] = seq[i:i+3]
                    break
                if seq[i:i+2] == "CC":
                    indel["has_pam"] = True
                    indel["pam_motif"] = seq[i:i+3]
                    break
        except Exception:
            pass

        # NHEJ signature check:
        # - Deletions 1-20bp are most common
        # - Insertions 1-2bp are frequent
        # - Larger deletions possible but rarer
        if indel["type"] == "deletion" and 1 <= indel["size"] <= 20:
            indel["nhej_signature"] = True
        elif indel["type"] == "insertion" and 1 <= indel["size"] <= 5:
            indel["nhej_signature"] = True

        # CRISPR candidate: has PAM + NHEJ signature + good quality
        if indel["has_pam"] and indel["nhej_signature"] and indel["qual"] >= 20:
            indel["crispr_candidate"] = True

        # Check zygosity from genotype
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

    # Summary
    crispr_candidates = [i for i in specific_indels if i.get("crispr_candidate")]
    nhej_indels = [i for i in specific_indels if i.get("nhej_signature")]

    log(f"\n=== Editing Site Summary for {args.sample_name} ===")
    log(f"Treatment-specific indels (>= {args.min_indel_size}bp): {len(specific_indels)}")
    log(f"With NHEJ signature: {len(nhej_indels)}")
    log(f"CRISPR candidates (PAM + NHEJ + QUAL>=20): {len(crispr_candidates)}")

    if crispr_candidates:
        log("\nTop CRISPR editing candidates:")
        # Sort by quality descending
        crispr_candidates.sort(key=lambda x: -x["qual"])
        for indel in crispr_candidates[:20]:
            log(f"  {indel['chrom']}:{indel['pos']} "
                f"{indel['type']} {indel['size']}bp "
                f"(QUAL={indel['qual']:.0f}, DP={indel['dp']}, "
                f"GT={indel['gt']}, {indel['zygosity']}, "
                f"PAM={indel['pam_motif']})")

    log(f"\nOutput: {output_tsv}")
    log("Done.")


if __name__ == "__main__":
    main()
