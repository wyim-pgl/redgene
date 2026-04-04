#!/usr/bin/env python3
"""Step 8: CRISPR editing site detection via indel calling.

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
  python s08_indel_detection.py \
    --treatment-bam results/{sample}/s07_host_map/{sample}_host.bam \
    --wt-bam results/WT/s07_host_map/WT_host.bam \
    --host-ref db/host.fa \
    --grna ATCGATCGATCGATCGATCG,GCTAGCTAGCTAGCTAGCTA \
    --outdir results --sample-name {sample}

  # Without gRNA (de novo)
  python s08_indel_detection.py \
    --treatment-bam results/{sample}/s07_host_map/{sample}_host.bam \
    --wt-bam results/WT/s07_host_map/WT_host.bam \
    --host-ref db/host.fa --genome-wide \
    --outdir results --sample-name {sample}
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path


def log(msg: str) -> None:
    print(f"[s08_indel] {msg}", file=sys.stderr, flush=True)


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
    """Find on-target and off-target sites for gRNA sequences using minimap2.

    Maps gRNA+PAM (23bp) to host genome with relaxed parameters to find
    sites with up to max_mismatches mismatches.
    """
    import tempfile

    targets = []

    for i, grna in enumerate(grna_seqs):
        grna = grna.upper().strip()
        if len(grna) < 17 or len(grna) > 25:
            log(f"WARNING: gRNA {i+1} length {len(grna)} unusual (expected 17-25bp)")

        # Create temp FASTA with gRNA + PAM (NGG)
        # We search for the gRNA itself; PAM will be checked at hit positions
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(f">gRNA_{i+1}\n{grna}\n")
            query_fa = f.name

        # Use minimap2 with short-read preset for sensitive gRNA mapping
        result = subprocess.run(
            ["minimap2", "-x", "sr", "-N", "1000", "--secondary=yes",
             "-c", str(host_ref), query_fa],
            capture_output=True, text=True,
        )

        Path(query_fa).unlink(missing_ok=True)

        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue

            q_len = int(fields[1])
            strand = fields[4]
            t_name = fields[5]
            t_start = int(fields[7])
            t_end = int(fields[8])
            matching = int(fields[9])
            block_len = int(fields[10])
            mapq = int(fields[11])

            mismatches = block_len - matching
            if mismatches > max_mismatches:
                continue

            # Check for PAM (NGG) at the expected position
            # For + strand: PAM is immediately after the gRNA (3' end)
            # For - strand: PAM (CCN on reference) is immediately before
            pam_start = t_end if strand == "+" else t_start - 3
            pam_end = pam_start + 3

            pam_check = subprocess.run(
                ["samtools", "faidx", str(host_ref),
                 f"{t_name}:{max(1,pam_start)}-{pam_end}"],
                capture_output=True, text=True,
            )
            pam_seq = "".join(pam_check.stdout.strip().split("\n")[1:]).upper()

            has_pam = False
            if strand == "+" and len(pam_seq) >= 3 and pam_seq[1:3] == "GG":
                has_pam = True
            elif strand == "-" and len(pam_seq) >= 3 and pam_seq[0:2] == "CC":
                has_pam = True

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
                "mapq": mapq,
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
    """Call variants at predicted gRNA target sites.

    Includes both indels and SNVs for comprehensive editing detection.
    """
    results = []

    for target in targets:
        chrom = target["chrom"]
        cut = target["cut_pos"]
        region = f"{chrom}:{max(1,cut-window)}-{cut+window}"

        # Call variants in treatment
        t_vars = _call_all_variants(treatment_bam, host_ref, region)
        # Call variants in WT
        wt_vars = _call_all_variants(wt_bam, host_ref, region)

        wt_set = {f"{v['chrom']}:{v['pos']}:{v['ref']}:{v['alt']}" for v in wt_vars}

        for var in t_vars:
            key = f"{var['chrom']}:{var['pos']}:{var['ref']}:{var['alt']}"
            if key in wt_set:
                continue  # Present in WT, skip

            # Determine variant type
            ref_len = len(var["ref"])
            alt_len = len(var["alt"])
            if ref_len == alt_len == 1:
                var["type"] = "SNV"
                var["size"] = 1
            elif ref_len > alt_len:
                var["type"] = "deletion"
                var["size"] = ref_len - alt_len
            else:
                var["type"] = "insertion"
                var["size"] = alt_len - ref_len

            if var["type"] != "SNV" and var["size"] < min_indel_size:
                continue

            var["grna_idx"] = target["grna_idx"]
            var["grna_seq"] = target["grna_seq"]
            var["site_type"] = target["site_type"]
            var["mismatches"] = target["mismatches"]
            var["distance_to_cut"] = abs(var["pos"] - cut)

            # Zygosity from genotype
            gt = var.get("gt", ".")
            if gt in ("1/1", "1|1"):
                var["zygosity"] = "homozygous"
            elif gt in ("0/1", "0|1", "1/0", "1|0"):
                var["zygosity"] = "heterozygous"
            else:
                var["zygosity"] = "unknown"

            results.append(var)

    return results


def _call_all_variants(bam: Path, host_ref: Path, region: str) -> list[dict]:
    """Call all variants (indels + SNVs) in a region."""
    cmd = ["bcftools", "mpileup", "-f", str(host_ref),
           "-r", region, "-d", "500", "-q", "10", "-Q", "20", str(bam)]
    mpileup = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    call_cmd = ["bcftools", "call", "-mv", "--ploidy", "2", "-Ov"]
    call_proc = subprocess.Popen(call_cmd, stdin=mpileup.stdout,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  text=True)
    mpileup.stdout.close()
    stdout, _ = call_proc.communicate()
    mpileup.wait()

    variants = []
    for line in stdout.strip().split("\n"):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 8:
            continue

        info = fields[7]
        dp = 0
        if "DP=" in info:
            dp = int(info.split("DP=")[1].split(";")[0])

        gt_data = fields[9] if len(fields) > 9 else ""

        variants.append({
            "chrom": fields[0],
            "pos": int(fields[1]),
            "ref": fields[3],
            "alt": fields[4],
            "qual": float(fields[5]) if fields[5] != "." else 0,
            "dp": dp,
            "gt": gt_data.split(":")[0] if gt_data else ".",
        })

    return variants


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
        description="Step 8: CRISPR editing site detection",
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

    step_dir = args.outdir / args.sample_name / "s08_indel"
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
            header = ["chrom", "pos", "ref", "alt", "type", "size", "qual",
                      "dp", "gt", "zygosity", "grna_idx", "site_type",
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
                    f"{v['type']} {v['size']}bp "
                    f"(QUAL={v['qual']:.0f}, GT={v['gt']}, {v['zygosity']}, "
                    f"dist={v['distance_to_cut']}bp from cut)")

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
