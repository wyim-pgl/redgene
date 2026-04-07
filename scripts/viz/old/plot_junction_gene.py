#!/usr/bin/env python3
"""Junction visualization with gene/exon/CDS annotation.

Redesigned layout:
  Row 1: Gene track — all genes in a single wide row, junction centered
  Row 2: Funnel/wedge contig — host narrows to junction, construct shown
         at true scale with kb annotation
  Row 3: Legend

Usage:
  python plot_junction_gene.py \
    --junctions results/{sample}/s06_junction/junctions.tsv \
    --gff db/genome.gff3 \
    --contigs results/{sample}/s04_assembly/contigs.fasta \
    --sample-name {sample} --outdir results
"""

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Polygon
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
    """Parse GFF3 for genes overlapping a region. Returns dict with genes list."""
    query_start = region_start - flank
    query_end = region_end + flank

    opener = gzip.open if str(gff_path).endswith(".gz") else open
    mode = "rt" if str(gff_path).endswith(".gz") else "r"

    genes = {}
    mrna_to_gene = {}

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

            if f_end < query_start or f_start > query_end:
                continue

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
                desc = desc.replace("%2C", ",").replace("%20", " ")
                genes[feat_id] = {
                    "id": feat_id, "name": name, "description": desc,
                    "start": f_start, "end": f_end, "strand": f_strand,
                    "exons": [], "cds": [], "utrs": [],
                }
            elif f_type == "mRNA":
                if parent in genes:
                    mrna_to_gene[feat_id] = parent
            elif f_type == "exon":
                gid = mrna_to_gene.get(parent, parent)
                if gid in genes:
                    genes[gid]["exons"].append((f_start, f_end))
            elif f_type == "CDS":
                frame = int(fields[7]) if fields[7] != "." else 0
                gid = mrna_to_gene.get(parent, parent)
                if gid in genes:
                    genes[gid]["cds"].append((f_start, f_end, frame))
            elif f_type in ("five_prime_UTR", "three_prime_UTR"):
                gid = mrna_to_gene.get(parent, parent)
                if gid in genes:
                    genes[gid]["utrs"].append((f_start, f_end))

    for g in genes.values():
        g["exons"] = sorted(set(g["exons"]))
        g["cds"] = sorted(set(g["cds"]))
        g["utrs"] = sorted(set(g["utrs"]))

    return {"genes": list(genes.values())}


# ---------------------------------------------------------------------------
# Contig FASTA
# ---------------------------------------------------------------------------

def read_contigs(fasta_path: Path) -> dict[str, str]:
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
# Colors
# ---------------------------------------------------------------------------

HOST_COLOR = "#4CAF50"
CONSTRUCT_COLOR = "#E53935"
EXON_COLOR = "#1565C0"
CDS_COLOR = "#0D47A1"
UTR_COLOR = "#90CAF9"
INTRON_COLOR = "#78909C"
JUNCTION_COLOR = "#FF6F00"
GENE_COLORS = ["#1565C0", "#6A1B9A", "#00838F", "#2E7D32", "#E65100"]


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------

def _classify_junction_location(junction_pos: int, gene: dict) -> str:
    """Classify where junction falls relative to gene features."""
    for cs, ce, _ in gene["cds"]:
        if cs <= junction_pos <= ce:
            return "CDS"
    for es, ee in gene["exons"]:
        if es <= junction_pos <= ee:
            return "UTR"
    if gene["start"] <= junction_pos <= gene["end"]:
        return "intron"
    return "intergenic"


def _fmt_kb(bp: int) -> str:
    """Format bp as human-readable."""
    if abs(bp) >= 1000:
        return f"{bp/1000:.1f} kb"
    return f"{bp} bp"


def plot_junction_with_gene(
    junction: dict,
    gene_info: dict,
    contig_seq: str | None,
    sample_name: str,
    output_path: Path,
) -> None:
    """Create redesigned junction diagram: gene track + funnel contig."""

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

    element_parts = construct_element.split("|")
    element_name = element_parts[1] if len(element_parts) > 1 else construct_element[:30]

    genes = gene_info.get("genes", [])
    # Sort genes by start position
    genes = sorted(genes, key=lambda g: g["start"])
    # Limit to closest 5 genes to junction
    if len(genes) > 5:
        genes = sorted(genes, key=lambda g: abs((g["start"]+g["end"])/2 - junction_pos))[:5]
        genes = sorted(genes, key=lambda g: g["start"])

    # --- Window: center on junction ---
    flank = 10000
    # Expand to include all gene boundaries
    all_positions = [junction_pos - flank, junction_pos + flank]
    for g in genes:
        all_positions.extend([g["start"], g["end"]])
    win_start = min(all_positions) - 500
    win_end = max(all_positions) + 500
    # Ensure junction is roughly centered (at least 40-60% from left)
    win_width = win_end - win_start
    junc_frac = (junction_pos - win_start) / win_width if win_width > 0 else 0.5
    if junc_frac < 0.35:
        win_end = junction_pos + int((junction_pos - win_start) / 0.4 * 0.6)
    elif junc_frac > 0.65:
        win_start = junction_pos - int((win_end - junction_pos) / 0.4 * 0.6)
    win_width = win_end - win_start

    # --- Contig geometry ---
    host_aln_len = abs(host_end - host_start)
    construct_start = int(junction.get("construct_start", 0))
    construct_end = int(junction.get("construct_end", 0))
    construct_aln_len = abs(construct_end - construct_start)

    if junction_type in ("LB", "5prime"):
        host_contig_start = contig_len - host_aln_len
        host_contig_end = contig_len
        const_contig_start = 0
        const_contig_end = construct_aln_len
    else:
        host_contig_start = 0
        host_contig_end = host_aln_len
        const_contig_start = contig_len - construct_aln_len
        const_contig_end = contig_len

    host_contig_start = max(0, host_contig_start)
    const_contig_end = min(contig_len, const_contig_end)
    junction_on_contig = host_contig_start if junction_type in ("LB", "5prime") else host_contig_end

    # --- Figure ---
    fig = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(
        3, 1,
        height_ratios=[2.5, 2.0, 0.3],
        hspace=0.05,
    )

    # ==== ROW 1: Gene track (single wide row, all genes) ====
    ax_gene = fig.add_subplot(gs[0])
    ax_gene.set_xlim(win_start, win_end)

    n_genes = len(genes)
    # If genes overlap vertically, stack them
    gene_y_positions = []
    GENE_HEIGHT = 0.35
    GENE_SPACING = 0.85
    for gi, gene in enumerate(genes):
        # Simple stacking: check overlap with previously placed genes
        y = 0
        for prev_gi, prev_gene in enumerate(genes[:gi]):
            prev_y = gene_y_positions[prev_gi]
            if gene["start"] <= prev_gene["end"] + 200 and gene["end"] >= prev_gene["start"] - 200:
                y = max(y, prev_y + GENE_SPACING)
        gene_y_positions.append(y)

    max_y = max(gene_y_positions) if gene_y_positions else 0
    ax_gene.set_ylim(-0.8, max_y + GENE_SPACING + 0.3)

    # Chromosome bar (background)
    chr_y = -0.5
    ax_gene.barh(chr_y, win_width, left=win_start, height=0.12,
                 color="#E0E0E0", edgecolor="#BBB", linewidth=0.5, zorder=0)
    # Chromosome label
    ax_gene.text(win_start + win_width * 0.005, chr_y, host_chr,
                 fontsize=7, va="center", ha="left", color="#666", style="italic")

    # Draw each gene
    for gi, gene in enumerate(genes):
        y = gene_y_positions[gi]
        color = GENE_COLORS[gi % len(GENE_COLORS)]
        gene_start = gene["start"]
        gene_end = gene["end"]
        gene_strand = gene["strand"]

        # Intron line
        draw_start = max(gene_start, win_start)
        draw_end = min(gene_end, win_end)
        ax_gene.plot([draw_start, draw_end], [y, y],
                     color=INTRON_COLOR, linewidth=1.2, zorder=1)

        # Strand direction arrows
        arrow_step = max(500, (draw_end - draw_start) // 15)
        for apos in range(int(draw_start) + 200, int(draw_end) - 200, arrow_step):
            dx = 150 if gene_strand == "+" else -150
            ax_gene.annotate("", xy=(apos + dx, y), xytext=(apos, y),
                             arrowprops=dict(arrowstyle="->", color=INTRON_COLOR,
                                             lw=0.6), zorder=1)

        # Exons
        for es, ee in gene["exons"]:
            if ee < win_start or es > win_end:
                continue
            es_c = max(es, win_start)
            ee_c = min(ee, win_end)
            ax_gene.barh(y, ee_c - es_c, left=es_c, height=GENE_HEIGHT,
                         color=UTR_COLOR, edgecolor=color,
                         linewidth=0.8, zorder=3)

        # CDS (darker fill)
        for cs, ce, frame in gene["cds"]:
            if ce < win_start or cs > win_end:
                continue
            cs_c = max(cs, win_start)
            ce_c = min(ce, win_end)
            ax_gene.barh(y, ce_c - cs_c, left=cs_c, height=GENE_HEIGHT,
                         color=color, edgecolor="black",
                         linewidth=0.7, alpha=0.85, zorder=4)

        # UTRs
        for us, ue in gene["utrs"]:
            if ue < win_start or us > win_end:
                continue
            us_c = max(us, win_start)
            ue_c = min(ue, win_end)
            ax_gene.barh(y, ue_c - us_c, left=us_c, height=GENE_HEIGHT * 0.7,
                         color=UTR_COLOR, edgecolor=color,
                         linewidth=0.5, zorder=3)

        # Gene label
        location = _classify_junction_location(junction_pos, gene)
        label_text = f"{gene['name']} ({gene_strand})"
        if location != "intergenic":
            label_text += f" [{location}]"

        # Place label above gene, slightly offset
        label_x = max(win_start + win_width * 0.01,
                      min((gene_start + gene_end) / 2, win_end - win_width * 0.15))
        ax_gene.text(label_x, y + GENE_HEIGHT / 2 + 0.08, label_text,
                     fontsize=7.5, va="bottom", ha="center", color=color,
                     fontweight="bold",
                     bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                               edgecolor=color, alpha=0.85, linewidth=0.5))

    # Junction line through gene track
    ax_gene.axvline(x=junction_pos, color=JUNCTION_COLOR, linewidth=2.5,
                    linestyle="-", zorder=20, alpha=0.9)

    # Junction label at top
    ax_gene.annotate(
        f"Junction: {host_chr}:{junction_pos:,}\n[{junction_type}] MAPQ={host_mapq}",
        xy=(junction_pos, max_y + GENE_SPACING),
        xytext=(junction_pos, max_y + GENE_SPACING + 0.2),
        fontsize=8.5, fontweight="bold", color=JUNCTION_COLOR,
        ha="center", va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF3E0",
                  edgecolor=JUNCTION_COLOR, alpha=0.95),
    )

    # Scale bar
    scale_len = _pick_scale_bar(win_width)
    sb_x = win_end - scale_len - win_width * 0.02
    sb_y = -0.3
    ax_gene.plot([sb_x, sb_x + scale_len], [sb_y, sb_y],
                 color="black", linewidth=2, zorder=10)
    ax_gene.plot([sb_x, sb_x], [sb_y - 0.05, sb_y + 0.05],
                 color="black", linewidth=1.5, zorder=10)
    ax_gene.plot([sb_x + scale_len, sb_x + scale_len], [sb_y - 0.05, sb_y + 0.05],
                 color="black", linewidth=1.5, zorder=10)
    ax_gene.text(sb_x + scale_len / 2, sb_y - 0.12, _fmt_kb(scale_len),
                 fontsize=7, ha="center", va="top", color="black")

    ax_gene.set_title(
        f"Junction Gene Context — {sample_name}\n"
        f"Construct: {element_name} | Confidence: {confidence}",
        fontsize=12, fontweight="bold", pad=12,
    )
    ax_gene.set_yticks([])
    ax_gene.tick_params(labelbottom=False, bottom=False)
    for spine in ax_gene.spines.values():
        spine.set_visible(False)

    # ==== ROW 2: Funnel contig diagram ====
    ax_contig = fig.add_subplot(gs[1])

    # Contig uses its OWN scale, centered so junction_on_contig aligns
    # with junction_pos from the gene track above.
    # We map contig coordinates → plot x via a linear transform.
    contig_visual_width = win_width * 0.65  # contig spans ~65% of plot width
    contig_scale = contig_visual_width / contig_len  # plot-units per contig-bp

    def contig_to_x(cp: float) -> float:
        """Map contig position to plot x-coordinate (junction-aligned)."""
        return junction_pos + (cp - junction_on_contig) * contig_scale

    ax_contig.set_xlim(win_start, win_end)
    ax_contig.set_ylim(-1.8, 2.0)

    # -- Contig backbone (gray bar) --
    bar_y = 0.0
    bar_h = 0.55
    contig_x0 = contig_to_x(0)
    contig_x1 = contig_to_x(contig_len)
    ax_contig.barh(bar_y, contig_x1 - contig_x0, left=contig_x0,
                   height=bar_h, color="#F5F5F5", edgecolor="#999",
                   linewidth=0.8, zorder=1)

    # -- Host portion --
    hx0 = contig_to_x(host_contig_start)
    hx1 = contig_to_x(host_contig_end)
    ax_contig.barh(bar_y, hx1 - hx0, left=hx0, height=bar_h,
                   color=HOST_COLOR, alpha=0.75, edgecolor=HOST_COLOR,
                   linewidth=1.2, zorder=5)
    host_mid = (hx0 + hx1) / 2
    ax_contig.text(host_mid, bar_y + 0.02,
                   f"Host\n{host_chr}:{host_start:,}-{host_end:,}",
                   ha="center", va="center", fontsize=7.5, fontweight="bold",
                   color="white", zorder=6)
    # Host length label below
    ax_contig.text(host_mid, bar_y - bar_h / 2 - 0.12,
                   _fmt_kb(host_aln_len), ha="center", va="top",
                   fontsize=7, color=HOST_COLOR, fontweight="bold")

    # -- Construct portion --
    cx0 = contig_to_x(const_contig_start)
    cx1 = contig_to_x(const_contig_end)
    ax_contig.barh(bar_y, cx1 - cx0, left=cx0, height=bar_h,
                   color=CONSTRUCT_COLOR, alpha=0.75, edgecolor=CONSTRUCT_COLOR,
                   linewidth=1.2, zorder=5)
    const_mid = (cx0 + cx1) / 2
    ax_contig.text(const_mid, bar_y + 0.02,
                   f"Construct\n{element_name}",
                   ha="center", va="center", fontsize=7.5, fontweight="bold",
                   color="white", zorder=6)
    # Construct length label below with arrow
    ax_contig.annotate(
        "", xy=(cx0, bar_y - bar_h / 2 - 0.15),
        xytext=(cx1, bar_y - bar_h / 2 - 0.15),
        arrowprops=dict(arrowstyle="<->", color=CONSTRUCT_COLOR, lw=1.5),
        zorder=8)
    ax_contig.text(const_mid, bar_y - bar_h / 2 - 0.30,
                   f"{_fmt_kb(construct_aln_len)} inserted",
                   ha="center", va="top", fontsize=8, color=CONSTRUCT_COLOR,
                   fontweight="bold", zorder=8)

    # -- Funnel/wedge from gene track junction to contig junction --
    junc_x = junction_pos  # same in both tracks (aligned by design)
    funnel_top = 1.8
    funnel_bot = bar_y + bar_h / 2
    funnel_narrow = win_width * 0.003  # narrow end at top

    # The funnel widens from the junction line toward the construct side
    if junction_type in ("LB", "5prime"):
        # Construct is LEFT of junction
        funnel_verts = [
            (junc_x - funnel_narrow, funnel_top),
            (junc_x + funnel_narrow, funnel_top),
            (junc_x, funnel_bot),
            (cx0, funnel_bot),
            (junc_x - funnel_narrow, funnel_top),
        ]
    else:
        # Construct is RIGHT of junction
        funnel_verts = [
            (junc_x - funnel_narrow, funnel_top),
            (junc_x + funnel_narrow, funnel_top),
            (cx1, funnel_bot),
            (junc_x, funnel_bot),
            (junc_x - funnel_narrow, funnel_top),
        ]

    funnel_poly = Polygon(funnel_verts, closed=True,
                          facecolor=CONSTRUCT_COLOR, alpha=0.10,
                          edgecolor=CONSTRUCT_COLOR, linewidth=0.8,
                          linestyle="--", zorder=2)
    ax_contig.add_patch(funnel_poly)

    # Junction line (vertically aligned with gene track)
    ax_contig.axvline(x=junction_pos, color=JUNCTION_COLOR, linewidth=2.5,
                      linestyle="-", zorder=20, alpha=0.9)
    ax_contig.plot(junction_pos, funnel_top, marker="v", markersize=8,
                   color=JUNCTION_COLOR, zorder=25, clip_on=False)

    # -- Sequence around junction --
    if contig_seq:
        jx = min(max(0, junction_on_contig), len(contig_seq) - 1)
        seq_s = max(0, jx - 15)
        seq_e = min(len(contig_seq), jx + 15)
        host_part = contig_seq[seq_s:jx]
        const_part = contig_seq[jx:seq_e]
        ax_contig.text(junction_pos, bar_y - bar_h / 2 - 0.65,
                       f"...{host_part}|{const_part}...",
                       ha="center", va="top", fontsize=7.5,
                       fontfamily="monospace",
                       bbox=dict(facecolor="white", edgecolor="#999",
                                 pad=3, boxstyle="round,pad=0.3"),
                       zorder=10)

    # Contig info
    ax_contig.text(win_start + win_width * 0.01, -1.55,
                   f"Contig: {contig_name}  ({_fmt_kb(contig_len)} total)",
                   fontsize=7.5, va="top", color="#555", style="italic")

    # Contig scale bar
    ctg_sb_len = _pick_scale_bar(int(contig_len))
    ctg_sb_x0 = contig_to_x(contig_len - ctg_sb_len - 5)
    ctg_sb_x1 = contig_to_x(contig_len - 5)
    sb_y2 = -1.3
    ax_contig.plot([ctg_sb_x0, ctg_sb_x1], [sb_y2, sb_y2],
                   color="black", linewidth=2, zorder=10)
    ax_contig.plot([ctg_sb_x0, ctg_sb_x0], [sb_y2 - 0.05, sb_y2 + 0.05],
                   color="black", linewidth=1.5, zorder=10)
    ax_contig.plot([ctg_sb_x1, ctg_sb_x1], [sb_y2 - 0.05, sb_y2 + 0.05],
                   color="black", linewidth=1.5, zorder=10)
    ax_contig.text((ctg_sb_x0 + ctg_sb_x1) / 2, sb_y2 - 0.12,
                   f"{_fmt_kb(ctg_sb_len)} (contig scale)",
                   fontsize=6.5, ha="center", va="top", color="#555")

    ax_contig.set_yticks([])
    ax_contig.tick_params(labelbottom=False, bottom=False)
    for spine in ax_contig.spines.values():
        spine.set_visible(False)

    # ==== ROW 3: Legend ====
    ax_leg = fig.add_subplot(gs[2])
    ax_leg.axis("off")

    legend_patches = [
        mpatches.Patch(color=HOST_COLOR, alpha=0.7, label="Host genome"),
        mpatches.Patch(color=CONSTRUCT_COLOR, alpha=0.7, label="Construct/T-DNA"),
        mpatches.Patch(color=CDS_COLOR, alpha=0.85, label="CDS"),
        mpatches.Patch(color=UTR_COLOR, label="Exon/UTR"),
        plt.Line2D([0], [0], color=INTRON_COLOR, linewidth=1.5, label="Intron"),
        plt.Line2D([0], [0], color=JUNCTION_COLOR, linewidth=2.5, label="Junction site"),
    ]
    ax_leg.legend(handles=legend_patches, loc="center", ncol=6,
                  fontsize=8.5, framealpha=0.9, edgecolor="#ccc")

    plt.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    log(f"Saved: {output_path}")


def _pick_scale_bar(win_width: int) -> int:
    """Pick a nice round scale bar length (~10-25% of window)."""
    candidates = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000]
    for c in candidates:
        if c >= win_width * 0.08 and c <= win_width * 0.30:
            return c
    # Fallback: pick closest to 15% of width
    target = int(win_width * 0.15)
    return min(candidates, key=lambda c: abs(c - target))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

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
    parser.add_argument("--confidence-filter", type=str, default=None,
                        help="Only plot junctions with this confidence or higher "
                             "(Low/Medium/High). Default: plot all.")
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

    # Filter by confidence if requested
    conf_order = {"Low": 0, "Medium": 1, "High": 2}
    if args.confidence_filter:
        min_conf = conf_order.get(args.confidence_filter, 0)
        junctions = [j for j in junctions
                     if conf_order.get(j.get("confidence", "Low"), 0) >= min_conf]

    # Read contigs
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

        gene_info = parse_gff3_region(
            args.gff, host_chr, junction_pos, junction_pos,
            flank=args.flank,
        )

        n_genes = len(gene_info.get("genes", []))
        log(f"    Found {n_genes} gene(s) in +/-{args.flank}bp")

        contig_seq = contig_seqs.get(contig_name)

        for ext in ("png", "pdf"):
            out_path = out_dir / f"junction_gene_{i+1}_{host_chr}_{junction_pos}.{ext}"
            plot_junction_with_gene(
                junction=junc,
                gene_info=gene_info,
                contig_seq=contig_seq,
                sample_name=args.sample_name,
                output_path=out_path,
            )

    log("All junction gene context plots generated.")


if __name__ == "__main__":
    main()
