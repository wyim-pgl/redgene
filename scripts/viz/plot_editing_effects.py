#!/usr/bin/env python3
"""CRISPR editing variant effect visualization.

Annotates CRISPR editing sites with functional impact:
  - Synonymous / Nonsynonymous / Frameshift / Splice site / UTR / Intergenic
  - Amino acid changes for in-frame indels and SNPs
  - Codon context around editing site

Generates a publication-ready figure combining:
  1. Gene model with exon/CDS/intron structure
  2. Nucleotide-level view at editing site with codon grid
  3. Variant effect summary (easyGWAS-inspired color coding)
  4. Allele frequency bar chart per editing site

Color scheme (easyGWAS-inspired):
  - Nonsynonymous/Missense: #E53935 (red)
  - Frameshift: #B71C1C (dark red)
  - Synonymous: #43A047 (green)
  - Splice site: #FF6F00 (amber)
  - UTR: #1E88E5 (blue)
  - Intron: #78909C (gray)
  - Intergenic: #BDBDBD (light gray)

Usage:
  python plot_editing_effects.py \
    --editing-sites results/{sample}/s08_indel/editing_sites.tsv \
    --gff db/genome.gff3 \
    --host-ref db/host.fa \
    --sample-name {sample} \
    --outdir results
"""

import argparse
import csv
import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np


def log(msg: str) -> None:
    print(f"[plot_effects] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Color scheme (easyGWAS-inspired)
# ---------------------------------------------------------------------------

EFFECT_COLORS = {
    "frameshift": "#B71C1C",
    "nonsynonymous": "#E53935",
    "missense": "#E53935",
    "synonymous": "#43A047",
    "splice_site": "#FF6F00",
    "UTR": "#1E88E5",
    "intron": "#78909C",
    "intergenic": "#BDBDBD",
    "in-frame_deletion": "#FF7043",
    "in-frame_insertion": "#FF7043",
    "stop_gained": "#880E4F",
    "stop_lost": "#AD1457",
}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA_NAMES = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
    "*": "Stop",
}


# ---------------------------------------------------------------------------
# GFF3 parsing for CDS context
# ---------------------------------------------------------------------------

def find_overlapping_cds(
    gff_path: Path, chrom: str, pos: int,
) -> dict | None:
    """Find CDS region containing a position. Returns gene info + CDS coords."""
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    mode = "rt" if str(gff_path).endswith(".gz") else "r"

    genes = {}
    mrna_to_gene = {}
    gene_cds = defaultdict(list)
    gene_exons = defaultdict(list)

    with opener(gff_path, mode) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[0] != chrom:
                continue

            f_type = fields[2]
            f_start = int(fields[3])
            f_end = int(fields[4])
            f_strand = fields[6]
            attrs = fields[8]

            # Only look within 50kb
            if abs(f_start - pos) > 50000 and abs(f_end - pos) > 50000:
                continue

            attr_dict = {}
            for attr in attrs.split(";"):
                if "=" in attr:
                    k, v = attr.split("=", 1)
                    attr_dict[k] = v

            feat_id = attr_dict.get("ID", "")
            parent = attr_dict.get("Parent", "")

            if f_type == "gene" and f_start <= pos <= f_end:
                name = attr_dict.get("Name", feat_id)
                desc = attr_dict.get("description", "").replace("%2C", ",").replace("%20", " ")
                genes[feat_id] = {
                    "id": feat_id, "name": name, "description": desc,
                    "start": f_start, "end": f_end, "strand": f_strand,
                }
            elif f_type == "mRNA":
                if parent in genes:
                    mrna_to_gene[feat_id] = parent
            elif f_type == "CDS":
                gid = mrna_to_gene.get(parent, parent)
                if gid in genes:
                    frame = int(fields[7]) if fields[7] != "." else 0
                    gene_cds[gid].append((f_start, f_end, frame))
            elif f_type == "exon":
                gid = mrna_to_gene.get(parent, parent)
                if gid in genes:
                    gene_exons[gid].append((f_start, f_end))

    # Find gene where position falls in CDS
    for gid, cds_list in gene_cds.items():
        for cs, ce, frame in sorted(cds_list):
            if cs <= pos <= ce:
                return {
                    **genes[gid],
                    "cds": sorted(set(cds_list)),
                    "exons": sorted(set(gene_exons.get(gid, []))),
                    "hit_cds": (cs, ce, frame),
                }

    # Check if in exon (UTR) or intron
    for gid in genes:
        gene = genes[gid]
        for es, ee in gene_exons.get(gid, []):
            if es <= pos <= ee:
                return {
                    **gene,
                    "cds": sorted(set(gene_cds.get(gid, []))),
                    "exons": sorted(set(gene_exons.get(gid, []))),
                    "hit_cds": None,
                    "location": "UTR",
                }
        if gene["start"] <= pos <= gene["end"]:
            return {
                **gene,
                "cds": sorted(set(gene_cds.get(gid, []))),
                "exons": sorted(set(gene_exons.get(gid, []))),
                "hit_cds": None,
                "location": "intron",
            }

    return None


def get_ref_seq(ref: Path, chrom: str, start: int, end: int) -> str:
    """Get reference sequence for a region."""
    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        ["samtools", "faidx", str(ref), region],
        capture_output=True, text=True,
    )
    return "".join(result.stdout.strip().split("\n")[1:]).upper()


def rev_comp(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


# ---------------------------------------------------------------------------
# Variant effect prediction
# ---------------------------------------------------------------------------

def predict_effect(
    chrom: str, pos: int, ref: str, alt: str,
    indel_type: str, indel_size: int, indel_seq: str,
    gene_info: dict | None, host_ref: Path,
) -> dict:
    """Predict functional effect of a variant.

    Returns dict with: effect, aa_change, codon_ref, codon_alt, aa_pos, details
    """
    result = {
        "effect": "intergenic",
        "aa_change": "",
        "codon_ref": "",
        "codon_alt": "",
        "aa_pos": 0,
        "gene_name": "",
        "gene_desc": "",
        "details": "",
    }

    if gene_info is None:
        return result

    result["gene_name"] = gene_info["name"]
    result["gene_desc"] = gene_info.get("description", "")
    strand = gene_info["strand"]

    # Check location
    if gene_info.get("hit_cds") is None:
        loc = gene_info.get("location", "intergenic")
        result["effect"] = loc
        result["details"] = f"In {loc} of {gene_info['name']}"
        return result

    # Position is in CDS
    cds_list = gene_info["cds"]  # sorted list of (start, end, frame)

    # Build CDS coordinate map
    cds_positions = []  # list of genomic positions in CDS order
    for cs, ce, _ in cds_list:
        cds_positions.extend(range(cs, ce + 1))

    if strand == "-":
        cds_positions = cds_positions[::-1]

    # Adjust for reading frame phase (GFF3 phase field on first CDS exon)
    first_phase = cds_list[0][2] if strand == "+" else cds_list[-1][2]
    if first_phase > 0:
        cds_positions = cds_positions[first_phase:]

    # Find position in CDS
    try:
        cds_idx = cds_positions.index(pos)
    except ValueError:
        result["effect"] = "CDS_boundary"
        result["details"] = f"At CDS boundary of {gene_info['name']}"
        return result

    # Check splice site proximity
    for cs, ce, _ in cds_list:
        if abs(pos - cs) <= 2 or abs(pos - ce) <= 2:
            result["effect"] = "splice_site"
            result["details"] = f"Near splice site of {gene_info['name']}"
            # Don't return yet - also check frameshift

    # For indels
    if indel_type in ("deletion", "insertion"):
        if indel_size % 3 != 0:
            result["effect"] = "frameshift"
            result["details"] = (
                f"Frameshift {indel_type} of {indel_size}bp "
                f"({indel_seq}) in {gene_info['name']}"
            )

            # Calculate amino acid position
            codon_pos = cds_idx // 3 + 1
            result["aa_pos"] = codon_pos
        else:
            # In-frame indel
            n_aa = indel_size // 3
            codon_pos = cds_idx // 3 + 1
            result["aa_pos"] = codon_pos

            if indel_type == "deletion":
                result["effect"] = "in-frame_deletion"
                # Get deleted amino acids
                if len(cds_positions) > cds_idx + indel_size:
                    del_start = min(cds_positions[cds_idx:cds_idx + indel_size + 3])
                    del_end = max(cds_positions[cds_idx:cds_idx + indel_size + 3])
                    del_seq = get_ref_seq(host_ref, chrom, del_start, del_end)
                    if strand == "-":
                        del_seq = rev_comp(del_seq)

                    # Translate deleted codons
                    del_aas = ""
                    for j in range(0, len(del_seq) - 2, 3):
                        codon = del_seq[j:j+3]
                        del_aas += CODON_TABLE.get(codon, "?")

                    result["aa_change"] = f"del {n_aa}aa ({del_aas})"
                else:
                    result["aa_change"] = f"del {n_aa}aa"

                result["details"] = (
                    f"In-frame {indel_size}bp deletion removes {n_aa} amino acid(s) "
                    f"at position {codon_pos} in {gene_info['name']}"
                )
            else:
                result["effect"] = "in-frame_insertion"
                result["aa_change"] = f"ins {n_aa}aa"
                result["details"] = (
                    f"In-frame {indel_size}bp insertion adds {n_aa} amino acid(s) "
                    f"at position {codon_pos} in {gene_info['name']}"
                )

        return result

    # For SNPs (point mutations)
    codon_offset = cds_idx % 3
    codon_start_idx = cds_idx - codon_offset

    if codon_start_idx + 3 > len(cds_positions):
        result["effect"] = "CDS_boundary"
        return result

    # Get codon genomic positions
    codon_genomic = cds_positions[codon_start_idx:codon_start_idx + 3]

    # Get reference codon
    ref_codon = ""
    for gp in codon_genomic:
        base = get_ref_seq(host_ref, chrom, gp, gp)
        ref_codon += base

    if strand == "-":
        ref_codon = rev_comp(ref_codon)

    # Make alt codon
    alt_base = alt[-1] if len(alt) > 0 else "N"
    if strand == "-":
        alt_base = rev_comp(alt_base)
    alt_codon = list(ref_codon)
    alt_codon[codon_offset] = alt_base
    alt_codon = "".join(alt_codon)

    ref_aa = CODON_TABLE.get(ref_codon, "?")
    alt_aa = CODON_TABLE.get(alt_codon, "?")
    aa_pos = codon_start_idx // 3 + 1

    result["codon_ref"] = ref_codon
    result["codon_alt"] = alt_codon
    result["aa_pos"] = aa_pos

    if ref_aa == alt_aa:
        result["effect"] = "synonymous"
        result["aa_change"] = f"p.{AA_NAMES.get(ref_aa, ref_aa)}{aa_pos}="
        result["details"] = (
            f"Synonymous: {ref_codon}>{alt_codon} "
            f"({ref_aa}{aa_pos}{alt_aa}) in {gene_info['name']}"
        )
    elif alt_aa == "*":
        result["effect"] = "stop_gained"
        result["aa_change"] = f"p.{AA_NAMES.get(ref_aa, ref_aa)}{aa_pos}*"
        result["details"] = (
            f"Stop gained: {ref_codon}>{alt_codon} "
            f"({ref_aa}{aa_pos}*) in {gene_info['name']}"
        )
    elif ref_aa == "*":
        result["effect"] = "stop_lost"
        result["aa_change"] = f"p.*{aa_pos}{AA_NAMES.get(alt_aa, alt_aa)}"
        result["details"] = (
            f"Stop lost: {ref_codon}>{alt_codon} "
            f"(*{aa_pos}{alt_aa}) in {gene_info['name']}"
        )
    else:
        result["effect"] = "nonsynonymous"
        result["aa_change"] = (
            f"p.{AA_NAMES.get(ref_aa, ref_aa)}{aa_pos}"
            f"{AA_NAMES.get(alt_aa, alt_aa)}"
        )
        result["details"] = (
            f"Missense: {ref_codon}>{alt_codon} "
            f"({ref_aa}{aa_pos}{alt_aa}) in {gene_info['name']}"
        )

    return result


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_editing_effects(
    sites_with_effects: list[dict],
    sample_name: str,
    output_path: Path,
) -> None:
    """Generate variant effect summary figure."""

    n_sites = len(sites_with_effects)
    if n_sites == 0:
        log("No editing sites to plot")
        return

    fig = plt.figure(figsize=(14, max(4, 2 + n_sites * 1.8)))
    gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1], hspace=0.35, wspace=0.3)

    # ---- Panel 1: Effect summary table (top-left) ----
    ax_table = fig.add_subplot(gs[0, 0])
    ax_table.axis("off")

    table_data = []
    row_colors = []
    for s in sites_with_effects:
        eff = s["effect_info"]
        site = s["site"]
        effect = eff["effect"]
        color = EFFECT_COLORS.get(effect, "#BDBDBD")
        row_colors.append((*matplotlib.colors.to_rgb(color), 0.2))

        gene = eff.get("gene_name", "")
        aa = eff.get("aa_change", "")
        freq = float(site.get("freq", 0))

        table_data.append([
            f"{site.get('chrom', '')}:{site.get('pos', '')}",
            f"{site.get('type', '')} {site.get('size', '')}bp",
            site.get("indel_seq", ""),
            effect.replace("_", " ").title(),
            gene,
            aa if aa else "-",
            f"{freq:.1%}",
            site.get("zygosity", ""),
        ])

    col_labels = ["Position", "Type", "Sequence", "Effect", "Gene",
                   "AA Change", "Freq", "Zygosity"]

    if table_data:
        table = ax_table.table(
            cellText=table_data,
            colLabels=col_labels,
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(7)
        table.scale(1, 1.4)

        # Color rows
        for i, color in enumerate(row_colors):
            for j in range(len(col_labels)):
                table[i + 1, j].set_facecolor(color)

        # Bold header
        for j in range(len(col_labels)):
            table[0, j].set_text_props(fontweight="bold", fontsize=8)
            table[0, j].set_facecolor("#E0E0E0")

    ax_table.set_title(
        f"CRISPR Editing Effects — {sample_name}",
        fontsize=12, fontweight="bold", pad=15,
    )

    # ---- Panel 2: Effect type pie chart (top-right) ----
    ax_pie = fig.add_subplot(gs[0, 1])
    effect_counts = defaultdict(int)
    for s in sites_with_effects:
        effect_counts[s["effect_info"]["effect"]] += 1

    labels = []
    sizes = []
    colors = []
    for eff, count in sorted(effect_counts.items()):
        labels.append(eff.replace("_", " ").title())
        sizes.append(count)
        colors.append(EFFECT_COLORS.get(eff, "#BDBDBD"))

    if sizes:
        wedges, texts, autotexts = ax_pie.pie(
            sizes, labels=labels, colors=colors, autopct="%1.0f%%",
            startangle=90, pctdistance=0.75,
            textprops={"fontsize": 8},
        )
        for t in autotexts:
            t.set_fontsize(7)
            t.set_fontweight("bold")
    ax_pie.set_title("Effect Distribution", fontsize=10, fontweight="bold")

    # ---- Panel 3: Allele frequency + depth bar chart (bottom-left) ----
    ax_freq = fig.add_subplot(gs[1, 0])

    x_pos = np.arange(n_sites)
    freqs = [float(s["site"].get("freq", 0)) for s in sites_with_effects]
    bar_colors = [EFFECT_COLORS.get(s["effect_info"]["effect"], "#BDBDBD")
                  for s in sites_with_effects]

    bars = ax_freq.bar(x_pos, [f * 100 for f in freqs], color=bar_colors,
                       edgecolor="black", linewidth=0.5, width=0.6)

    # Depth as text on bars
    for i, s in enumerate(sites_with_effects):
        dp = s["site"].get("dp", "?")
        count = s["site"].get("count", "?")
        ax_freq.text(i, freqs[i] * 100 + 1, f"{count}/{dp}",
                     ha="center", va="bottom", fontsize=7, fontweight="bold")

    ax_freq.set_xticks(x_pos)
    x_labels = [f"{s['site'].get('chrom', '')[-4:]}:{s['site'].get('pos', '')}"
                for s in sites_with_effects]
    ax_freq.set_xticklabels(x_labels, fontsize=7, rotation=30, ha="right")
    ax_freq.set_ylabel("Allele Frequency (%)", fontsize=9)
    ax_freq.set_title("Allele Frequency (reads/depth)", fontsize=10,
                       fontweight="bold")
    ax_freq.set_ylim(0, max(f * 100 for f in freqs) * 1.3 + 5 if freqs else 50)

    # ---- Panel 4: Codon context (bottom-right) ----
    ax_codon = fig.add_subplot(gs[1, 1])
    ax_codon.axis("off")

    codon_text = "Codon Context:\n\n"
    for s in sites_with_effects:
        eff = s["effect_info"]
        if eff.get("details"):
            codon_text += f"  {eff['details']}\n\n"

    ax_codon.text(0.05, 0.95, codon_text, transform=ax_codon.transAxes,
                  fontsize=8, va="top", fontfamily="monospace",
                  bbox=dict(boxstyle="round,pad=0.5", facecolor="#FAFAFA",
                            edgecolor="#ccc"))

    # Effect color legend
    legend_patches = [
        mpatches.Patch(color=c, label=e.replace("_", " ").title())
        for e, c in EFFECT_COLORS.items()
        if e in effect_counts
    ]
    if legend_patches:
        ax_codon.legend(handles=legend_patches, loc="lower right",
                        fontsize=7, framealpha=0.9)

    plt.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    log(f"Saved: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CRISPR editing variant effect annotation and visualization",
    )
    parser.add_argument("--editing-sites", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True,
                        help="GFF3 annotation file (may be .gz)")
    parser.add_argument("--host-ref", type=Path, required=True)
    parser.add_argument("--sample-name", type=str, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    out_dir = args.outdir / args.sample_name / "s08_indel"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load editing sites
    sites = []
    with open(args.editing_sites) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sites.append(row)

    if not sites:
        log("No editing sites found")
        return

    log(f"Annotating {len(sites)} editing sites")

    # Annotate each site
    annotated = []
    for site in sites:
        chrom = site["chrom"]
        pos = int(site["pos"])
        ref = site.get("ref", "")
        alt = site.get("alt", "")
        indel_type = site.get("type", "")
        indel_size = int(site.get("size", 0))
        indel_seq = site.get("indel_seq", "")

        log(f"  {chrom}:{pos} {indel_type} {indel_size}bp {indel_seq}")

        # Find gene/CDS context
        gene_info = find_overlapping_cds(args.gff, chrom, pos)

        if gene_info:
            gene_name = gene_info["name"]
            gene_desc = gene_info.get("description", "")
            log(f"    Gene: {gene_name} ({gene_desc[:50]})")
            if gene_info.get("hit_cds"):
                cs, ce, frame = gene_info["hit_cds"]
                log(f"    CDS: {cs}-{ce} (frame {frame})")
        else:
            log(f"    Intergenic")

        # Predict effect
        effect = predict_effect(
            chrom, pos, ref, alt,
            indel_type, indel_size, indel_seq,
            gene_info, args.host_ref,
        )

        log(f"    Effect: {effect['effect']} {effect.get('aa_change', '')}")
        log(f"    {effect.get('details', '')}")

        annotated.append({
            "site": site,
            "effect_info": effect,
        })

    # Write annotated TSV
    ann_tsv = out_dir / "editing_effects.tsv"
    with open(ann_tsv, "w") as f:
        headers = [
            "chrom", "pos", "type", "size", "indel_seq", "freq", "dp",
            "count", "zygosity", "effect", "gene", "aa_change", "details",
            "grna_idx", "site_type",
        ]
        f.write("\t".join(headers) + "\n")
        for a in annotated:
            s = a["site"]
            e = a["effect_info"]
            row = [
                s.get("chrom", ""), s.get("pos", ""),
                s.get("type", ""), s.get("size", ""),
                s.get("indel_seq", ""), s.get("freq", ""),
                s.get("dp", ""), s.get("count", ""),
                s.get("zygosity", ""),
                e["effect"], e.get("gene_name", ""),
                e.get("aa_change", ""), e.get("details", ""),
                s.get("grna_idx", ""), s.get("site_type", ""),
            ]
            f.write("\t".join(str(x) for x in row) + "\n")
    log(f"Annotated TSV: {ann_tsv}")

    # Generate plot
    out_png = out_dir / "editing_effects.png"
    out_pdf = out_dir / "editing_effects.pdf"

    plot_editing_effects(annotated, args.sample_name, out_png)
    plot_editing_effects(annotated, args.sample_name, out_pdf)

    log("Done.")


if __name__ == "__main__":
    main()
