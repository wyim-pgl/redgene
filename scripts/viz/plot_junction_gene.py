#!/usr/bin/env python3
"""Enhanced junction visualization with gene/exon/CDS annotation.

Shows each junction site in the context of host genome gene structure,
with the inserted construct sequence annotated alongside.

Panels per junction:
  - Gene model track (exons/CDS/UTR/introns) with junction position
  - Contig alignment diagram showing host and construct segments
  - Inserted sequence annotation

Usage:
  python plot_junction_gene.py \
    --junctions results/{sample}/s06_junction/junctions.tsv \
    --gff db/genome.gff3 \
    --contigs results/{sample}/s04_assembly/contigs.fasta \
    --sample-name {sample} \
    --output results/{sample}/s06_junction/junction_gene_context.png
"""

import argparse
import csv
import gzip
import re
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
    print(f"[plot_junction_gene] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# GFF3 parsing
# ---------------------------------------------------------------------------

def parse_gff3_region(
    gff_path: Path, chrom: str, region_start: int, region_end: int,
    flank: int = 10000,
) -> dict:
    """Parse GFF3 for genes/mRNA/exon/CDS features overlapping a region.

    Returns dict with:
      genes: [{id, name, description, start, end, strand,
               exons: [(s,e),...], cds: [(s,e,frame),...], utrs: [(s,e),...]}]
    """
    query_start = region_start - flank
    query_end = region_end + flank

    opener = gzip.open if str(gff_path).endswith(".gz") else open
    mode = "rt" if str(gff_path).endswith(".gz") else "r"

    # First pass: find genes and their hierarchy
    genes = {}  # gene_id -> gene info
    mrna_to_gene = {}  # mRNA_id -> gene_id
    gene_mrnas = defaultdict(list)  # gene_id -> [mRNA_ids]

    features_by_parent = defaultdict(list)

    with opener(gff_path, mode) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            f_chrom = fields[0]
            if f_chrom != chrom:
                continue

            f_type = fields[2]
            f_start = int(fields[3])
            f_end = int(fields[4])
            f_strand = fields[6]
            attrs = fields[8]

            # Check overlap
            if f_end < query_start or f_start > query_end:
                continue

            # Parse attributes
            attr_dict = {}
            for attr in attrs.split(";"):
                if "=" in attr:
                    k, v = attr.split("=", 1)
                    attr_dict[k] = v

            feat_id = attr_dict.get("ID", "")
            parent = attr_dict.get("Parent", "")

            if f_type == "gene":
                name = attr_dict.get("Name", feat_id)
                desc = attr_dict.get("description", "")
                # URL decode
                desc = desc.replace("%2C", ",").replace("%20", " ")
                genes[feat_id] = {
                    "id": feat_id,
                    "name": name,
                    "description": desc,
                    "start": f_start,
                    "end": f_end,
                    "strand": f_strand,
                    "exons": [],
                    "cds": [],
                    "utrs": [],
                }
            elif f_type == "mRNA":
                if parent in genes:
                    mrna_to_gene[feat_id] = parent
                    gene_mrnas[parent].append(feat_id)
            elif f_type == "exon":
                # Find which gene this belongs to
                if parent in mrna_to_gene:
                    gid = mrna_to_gene[parent]
                    genes[gid]["exons"].append((f_start, f_end))
                elif parent in genes:
                    genes[parent]["exons"].append((f_start, f_end))
            elif f_type == "CDS":
                frame = fields[7] if fields[7] != "." else "0"
                if parent in mrna_to_gene:
                    gid = mrna_to_gene[parent]
                    genes[gid]["cds"].append((f_start, f_end, int(frame)))
                elif parent in genes:
                    genes[parent]["cds"].append((f_start, f_end, int(frame)))
            elif f_type in ("five_prime_UTR", "three_prime_UTR"):
                if parent in mrna_to_gene:
                    gid = mrna_to_gene[parent]
                    genes[gid]["utrs"].append((f_start, f_end))

    # Sort exons and CDS
    for g in genes.values():
        g["exons"] = sorted(set(g["exons"]))
        g["cds"] = sorted(set(g["cds"]))
        g["utrs"] = sorted(set(g["utrs"]))

    return {"genes": list(genes.values())}


# ---------------------------------------------------------------------------
# Contig sequence extraction
# ---------------------------------------------------------------------------

def read_contigs(fasta_path: Path) -> dict[str, str]:
    """Read contigs from FASTA file."""
    contigs = {}
    current = None
    seq_parts = []

    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    mode = "rt" if str(fasta_path).endswith(".gz") else "r"

    with opener(fasta_path, mode) as f:
        for line in f:
            if line.startswith(">"):
                if current:
                    contigs[current] = "".join(seq_parts)
                current = line.strip()[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if current:
            contigs[current] = "".join(seq_parts)

    return contigs


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

# Colors
HOST_COLOR = "#4CAF50"       # green
CONSTRUCT_COLOR = "#F44336"  # red
EXON_COLOR = "#1565C0"       # blue
CDS_COLOR = "#0D47A1"        # dark blue
UTR_COLOR = "#90CAF9"        # light blue
INTRON_COLOR = "#78909C"     # gray
JUNCTION_COLOR = "#FF6F00"   # amber


def plot_junction_with_gene(
    junction: dict,
    gene_info: dict,
    contig_seq: str | None,
    sample_name: str,
    output_path: Path,
) -> None:
    """Create a multi-panel junction diagram with gene context."""

    host_chr = junction["host_chr"]
    host_start = int(junction["host_start"])
    host_end = int(junction["host_end"])
    junction_pos = int(junction["junction_pos_host"])
    junction_type = junction["junction_type"]
    contig_name = junction["contig_name"]
    contig_len = int(junction["contig_len"])
    construct_element = junction["construct_element"]
    confidence = junction["confidence"]
    host_mapq = junction.get("host_mapq", "?")

    # Parse construct element name
    element_parts = construct_element.split("|")
    element_type = element_parts[0] if len(element_parts) > 0 else "unknown"
    element_name = element_parts[1] if len(element_parts) > 1 else construct_element[:30]

    genes = gene_info.get("genes", [])

    # Determine plot window around junction
    flank = 5000
    plot_start = junction_pos - flank
    plot_end = junction_pos + flank

    # Adjust to show full gene if it extends beyond
    for g in genes:
        if g["start"] < plot_start and g["end"] > plot_start - flank:
            plot_start = g["start"] - 500
        if g["end"] > plot_end and g["start"] < plot_end + flank:
            plot_end = g["end"] + 500

    n_genes = len(genes)
    fig_height = max(6, 3 + n_genes * 1.5 + 2)

    fig = plt.figure(figsize=(16, fig_height))
    gs = gridspec.GridSpec(
        3 + n_genes, 1,
        height_ratios=[1.5] + [1.2] * n_genes + [1.5, 0.8],
        hspace=0.15,
    )

    # ---- Panel 1: Genomic context overview ----
    ax_ctx = fig.add_subplot(gs[0])
    ax_ctx.set_xlim(plot_start, plot_end)
    ax_ctx.set_ylim(-0.5, 0.5)

    # Chromosome bar
    ax_ctx.barh(0, plot_end - plot_start, left=plot_start, height=0.15,
                color="#E0E0E0", edgecolor="#999", linewidth=0.5, zorder=1)

    # Junction marker
    ax_ctx.axvline(x=junction_pos, color=JUNCTION_COLOR, linewidth=2.5,
                   linestyle="-", zorder=10)
    ax_ctx.annotate(
        f"Junction: {host_chr}:{junction_pos:,}\n[{junction_type}] MAPQ={host_mapq}",
        xy=(junction_pos, 0.15),
        xytext=(junction_pos, 0.4),
        fontsize=8, fontweight="bold", color=JUNCTION_COLOR,
        ha="center", va="bottom",
        arrowprops=dict(arrowstyle="->", color=JUNCTION_COLOR, lw=1.5),
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF3E0",
                  edgecolor=JUNCTION_COLOR, alpha=0.9),
    )

    # Host alignment region
    ax_ctx.barh(0, host_end - host_start, left=host_start, height=0.3,
                color=HOST_COLOR, alpha=0.3, edgecolor=HOST_COLOR,
                linewidth=1, zorder=2)

    ax_ctx.set_title(
        f"Junction Gene Context — {sample_name}\n"
        f"{host_chr}:{plot_start:,}-{plot_end:,} | "
        f"Construct: {element_name} | Confidence: {confidence}",
        fontsize=11, fontweight="bold", pad=10,
    )
    ax_ctx.set_ylabel("Genome", fontsize=9)
    ax_ctx.set_yticks([])
    ax_ctx.tick_params(labelbottom=False)
    ax_ctx.spines["left"].set_visible(False)

    # ---- Panels 2..N: Gene models ----
    for gi, gene in enumerate(genes):
        ax_gene = fig.add_subplot(gs[1 + gi])
        ax_gene.set_xlim(plot_start, plot_end)
        ax_gene.set_ylim(-0.8, 0.8)

        gene_name = gene["name"]
        gene_desc = gene.get("description", "")
        gene_strand = gene["strand"]
        gene_start = gene["start"]
        gene_end = gene["end"]

        # Intron line (gene body)
        ax_gene.plot([gene_start, gene_end], [0, 0],
                     color=INTRON_COLOR, linewidth=1.5, zorder=1)

        # Strand arrows along intron
        arrow_step = max(200, (gene_end - gene_start) // 20)
        for apos in range(gene_start, gene_end, arrow_step):
            if plot_start <= apos <= plot_end:
                dx = 80 if gene_strand == "+" else -80
                ax_gene.annotate("", xy=(apos + dx, 0), xytext=(apos, 0),
                                 arrowprops=dict(arrowstyle="->",
                                                 color=INTRON_COLOR, lw=0.8))

        # Draw exons (light boxes)
        for es, ee in gene["exons"]:
            if ee < plot_start or es > plot_end:
                continue
            ax_gene.barh(0, ee - es, left=es, height=0.5,
                         color=UTR_COLOR, edgecolor=EXON_COLOR,
                         linewidth=0.8, zorder=3)

        # Draw CDS (filled boxes, darker)
        for cs, ce, frame in gene["cds"]:
            if ce < plot_start or cs > plot_end:
                continue
            ax_gene.barh(0, ce - cs, left=cs, height=0.5,
                         color=CDS_COLOR, edgecolor="black",
                         linewidth=0.8, alpha=0.85, zorder=4)

            # Frame annotation (small text)
            mid = (cs + ce) / 2
            if plot_start <= mid <= plot_end:
                ax_gene.text(mid, -0.35, f"f{frame}",
                             fontsize=5, ha="center", va="top",
                             color="#666")

        # Draw UTRs
        for us, ue in gene["utrs"]:
            if ue < plot_start or us > plot_end:
                continue
            ax_gene.barh(0, ue - us, left=us, height=0.35,
                         color=UTR_COLOR, edgecolor=EXON_COLOR,
                         linewidth=0.5, zorder=3)

        # Junction line through gene
        ax_gene.axvline(x=junction_pos, color=JUNCTION_COLOR, linewidth=2,
                        linestyle="--", alpha=0.7, zorder=10)

        # Check if junction falls in CDS/exon/intron/UTR
        location = "intergenic"
        for cs, ce, _ in gene["cds"]:
            if cs <= junction_pos <= ce:
                location = "CDS"
                break
        if location == "intergenic":
            for es, ee in gene["exons"]:
                if es <= junction_pos <= ee:
                    location = "exon (UTR)"
                    break
        if location == "intergenic":
            if gene_start <= junction_pos <= gene_end:
                location = "intron"

        # Gene label
        label = f"{gene_name} ({gene_strand})"
        if gene_desc:
            label += f"\n{gene_desc[:60]}"
        label += f"\nJunction in: {location}"

        ax_gene.text(plot_start + (plot_end - plot_start) * 0.02, 0.65,
                     label, fontsize=8, va="top",
                     bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                               edgecolor="#ccc", alpha=0.9))

        ax_gene.set_ylabel(f"Gene {gi+1}", fontsize=8)
        ax_gene.set_yticks([])
        ax_gene.tick_params(labelbottom=False)
        ax_gene.spines["left"].set_visible(False)

    # ---- Panel: Contig alignment diagram ----
    ax_contig = fig.add_subplot(gs[1 + n_genes])
    ax_contig.set_xlim(0, contig_len)
    ax_contig.set_ylim(-0.8, 1.2)

    # Contig backbone
    ax_contig.barh(0, contig_len, height=0.4, color="#E0E0E0",
                   edgecolor="black", linewidth=0.8, zorder=1)

    # Host alignment portion (approximate)
    host_aln_len = abs(host_end - host_start)
    construct_start = int(junction.get("construct_start", 0))
    construct_end = int(junction.get("construct_end", 0))
    construct_aln_len = abs(construct_end - construct_start)

    # Determine which side is host vs construct
    if junction_type in ("LB", "5prime"):
        # Host on right, construct on left
        host_contig_start = contig_len - host_aln_len
        host_contig_end = contig_len
        const_contig_start = 0
        const_contig_end = construct_aln_len
    else:
        # Host on left, construct on right
        host_contig_start = 0
        host_contig_end = host_aln_len
        const_contig_start = contig_len - construct_aln_len
        const_contig_end = contig_len

    # Clamp
    host_contig_start = max(0, host_contig_start)
    const_contig_end = min(contig_len, const_contig_end)

    # Draw host region
    ax_contig.barh(0, host_contig_end - host_contig_start,
                   left=host_contig_start, height=0.4,
                   color=HOST_COLOR, alpha=0.6, edgecolor=HOST_COLOR,
                   linewidth=1, zorder=3)
    ax_contig.text((host_contig_start + host_contig_end) / 2, 0,
                   f"Host\n{host_chr}:{host_start:,}-{host_end:,}",
                   ha="center", va="center", fontsize=7,
                   fontweight="bold", color="white")

    # Draw construct region
    ax_contig.barh(0, const_contig_end - const_contig_start,
                   left=const_contig_start, height=0.4,
                   color=CONSTRUCT_COLOR, alpha=0.6, edgecolor=CONSTRUCT_COLOR,
                   linewidth=1, zorder=3)
    ax_contig.text((const_contig_start + const_contig_end) / 2, 0,
                   f"Construct\n{element_name}",
                   ha="center", va="center", fontsize=7,
                   fontweight="bold", color="white")

    # Junction point on contig
    junction_on_contig = host_contig_end if junction_type in ("LB", "5prime") else host_contig_end
    ax_contig.axvline(x=junction_on_contig, color=JUNCTION_COLOR,
                      linewidth=2.5, zorder=10)
    ax_contig.annotate("Junction", xy=(junction_on_contig, 0.25),
                       xytext=(junction_on_contig, 0.9),
                       fontsize=8, fontweight="bold", color=JUNCTION_COLOR,
                       ha="center",
                       arrowprops=dict(arrowstyle="->", color=JUNCTION_COLOR))

    # Show contig sequence around junction (if available)
    if contig_seq:
        # Show ~20bp around junction
        jx = min(junction_on_contig, len(contig_seq) - 1)
        jx = max(0, jx)
        seq_start = max(0, jx - 15)
        seq_end = min(len(contig_seq), jx + 15)
        seq_around = contig_seq[seq_start:seq_end]

        # Color-code: host side green, construct side red
        host_side_len = jx - seq_start
        host_part = seq_around[:host_side_len]
        const_part = seq_around[host_side_len:]

        ax_contig.text(junction_on_contig, -0.5,
                       f"...{host_part}|{const_part}...",
                       ha="center", va="top", fontsize=7,
                       fontfamily="monospace",
                       bbox=dict(facecolor="white", edgecolor="#ccc",
                                 pad=2, boxstyle="round"))

    ax_contig.set_ylabel("Contig", fontsize=8)
    ax_contig.set_xlabel(f"Contig position (bp) — {contig_name}", fontsize=8)
    ax_contig.set_yticks([])
    ax_contig.spines["left"].set_visible(False)

    # ---- Panel: Legend ----
    ax_leg = fig.add_subplot(gs[2 + n_genes])
    ax_leg.axis("off")

    legend_patches = [
        mpatches.Patch(color=HOST_COLOR, alpha=0.6, label="Host genome"),
        mpatches.Patch(color=CONSTRUCT_COLOR, alpha=0.6, label="Construct/T-DNA"),
        mpatches.Patch(color=CDS_COLOR, alpha=0.85, label="CDS"),
        mpatches.Patch(color=UTR_COLOR, label="Exon/UTR"),
        plt.Line2D([0], [0], color=INTRON_COLOR, linewidth=1.5, label="Intron"),
        plt.Line2D([0], [0], color=JUNCTION_COLOR, linewidth=2.5, label="Junction site"),
    ]
    ax_leg.legend(handles=legend_patches, loc="center", ncol=6,
                  fontsize=8, framealpha=0.9)

    plt.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    log(f"Saved: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Junction visualization with gene/exon/CDS context",
    )
    parser.add_argument("--junctions", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True,
                        help="GFF3 annotation file (may be .gz)")
    parser.add_argument("--contigs", type=Path, default=None,
                        help="Contigs FASTA from SPAdes assembly")
    parser.add_argument("--sample-name", type=str, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=10000,
                        help="Flanking region for gene search (default: 10kb)")
    args = parser.parse_args()

    out_dir = args.outdir / args.sample_name / "s06_junction"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Read junctions
    junctions = []
    with open(args.junctions) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            junctions.append(row)

    if not junctions:
        log("No junctions found")
        return

    # Read contigs if available
    contig_seqs = {}
    if args.contigs and args.contigs.exists():
        contig_seqs = read_contigs(args.contigs)
        log(f"Loaded {len(contig_seqs)} contigs")

    log(f"Processing {len(junctions)} junctions")

    for i, junc in enumerate(junctions):
        host_chr = junc["host_chr"]
        junction_pos = int(junc["junction_pos_host"])
        contig_name = junc["contig_name"]

        log(f"  Junction {i+1}: {host_chr}:{junction_pos}")

        # Get gene info from GFF3
        gene_info = parse_gff3_region(
            args.gff, host_chr, junction_pos, junction_pos,
            flank=args.flank,
        )

        n_genes = len(gene_info.get("genes", []))
        log(f"    Found {n_genes} gene(s) in ±{args.flank}bp")

        # Get contig sequence
        contig_seq = contig_seqs.get(contig_name)

        out_png = out_dir / f"junction_gene_{i+1}_{host_chr}_{junction_pos}.png"
        out_pdf = out_dir / f"junction_gene_{i+1}_{host_chr}_{junction_pos}.pdf"

        plot_junction_with_gene(
            junction=junc,
            gene_info=gene_info,
            contig_seq=contig_seq,
            sample_name=args.sample_name,
            output_path=out_png,
        )

        # Also save PDF
        plot_junction_with_gene(
            junction=junc,
            gene_info=gene_info,
            contig_seq=contig_seq,
            sample_name=args.sample_name,
            output_path=out_pdf,
        )

    log("All junction gene context plots generated.")


if __name__ == "__main__":
    main()
