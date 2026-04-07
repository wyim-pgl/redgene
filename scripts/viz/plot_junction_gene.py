#!/usr/bin/env python3
"""Junction visualization with gene/exon/CDS annotation + element map + depth.

Layout:
  Row 1: Gene track — genes in the region, junction centered
  Row 2: Full contig — all PAF element annotations shown
  Row 3: Junction sequence (color-coded host|construct)
  Row 4: Read depth profile on the contig
  Row 5: Legend

Usage:
  python plot_junction_gene.py \
    --junctions results/{sample}/s06_junction/junctions.tsv \
    --gff db/genome.gff3 \
    --contigs results/{sample}/s04_assembly/contigs.fasta \
    --construct-paf results/{sample}/s05_contig_map/{sample}_contigs_to_construct.paf \
    --host-paf results/{sample}/s05_contig_map/{sample}_contigs_to_host.paf \
    --reads-r1 results/{sample}/s03_extract/{sample}_construct_R1.fq.gz \
    --reads-r2 results/{sample}/s03_extract/{sample}_construct_R2.fq.gz \
    --sample-name {sample} --outdir results
"""

import argparse
import csv
import gzip
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon
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
    """Parse GFF3 for genes overlapping a region."""
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
# FASTA / PAF parsing
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


def parse_paf_for_contig(paf_path: Path, contig_name: str) -> list[dict]:
    """Parse PAF and return all primary alignments for a given contig."""
    hits = []
    with open(paf_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 12:
                continue
            if fields[0] != contig_name:
                continue
            # Only keep primary alignments (tp:A:P)
            tp = ""
            mapq = int(fields[11])
            for tag in fields[12:]:
                if tag.startswith("tp:A:"):
                    tp = tag[5:]
            if tp != "P":
                continue
            hits.append({
                "contig_start": int(fields[2]),
                "contig_end": int(fields[3]),
                "strand": fields[4],
                "target": fields[5],
                "target_len": int(fields[6]),
                "target_start": int(fields[7]),
                "target_end": int(fields[8]),
                "match_bp": int(fields[9]),
                "aln_len": int(fields[10]),
                "mapq": mapq,
            })
    return hits


def _element_short_name(target: str) -> str:
    """Extract short display name from element_db target string."""
    parts = target.split("|")
    if len(parts) >= 2:
        return parts[1]
    return target[:25]


def _element_category(target: str) -> str:
    """Extract category from element_db target string."""
    parts = target.split("|")
    if parts:
        return parts[0]
    return "unknown"


# ---------------------------------------------------------------------------
# Read depth via minimap2
# ---------------------------------------------------------------------------

def compute_contig_depth(
    contig_name: str, contig_seq: str,
    reads_r1: Path | None, reads_r2: Path | None,
) -> np.ndarray | None:
    """Map extracted reads to a single contig, return per-base depth array."""
    if not reads_r1 or not reads_r1.exists():
        return None

    try:
        import pysam
    except ImportError:
        log("pysam not available, skipping depth track")
        return None

    with tempfile.TemporaryDirectory(prefix="junc_depth_") as tmpdir:
        # Write single contig FASTA
        ref_fa = os.path.join(tmpdir, "ref.fa")
        with open(ref_fa, "w") as f:
            f.write(f">{contig_name}\n{contig_seq}\n")

        # Map reads with minimap2
        bam_path = os.path.join(tmpdir, "aln.bam")
        reads_args = [str(reads_r1)]
        if reads_r2 and reads_r2.exists():
            reads_args.append(str(reads_r2))

        cmd = (
            f"minimap2 -a -x sr -t 2 {ref_fa} {' '.join(reads_args)} 2>/dev/null"
            f" | samtools sort -@ 2 -o {bam_path} && samtools index {bam_path}"
        )
        result = subprocess.run(cmd, shell=True, capture_output=True, timeout=120)
        if result.returncode != 0:
            log(f"minimap2/samtools failed: {result.stderr.decode()[:200]}")
            return None

        # Compute depth
        depth = np.zeros(len(contig_seq), dtype=np.int32)
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for col in bam.pileup(contig_name, min_mapping_quality=0,
                                   min_base_quality=0):
                if 0 <= col.pos < len(depth):
                    depth[col.pos] = col.nsegments

    return depth


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

# Distinct colors for different construct element categories
ELEMENT_COLORS = {
    "vector": "#E53935",
    "promoter": "#FF9800",
    "terminator": "#9C27B0",
    "CDS": "#2196F3",
    "element-specific": "#F44336",
    "event-specific": "#E91E63",
    "taxon-specific": "#795548",
    "euginius": "#FF5722",
    "regulatory": "#00BCD4",
}


def _get_element_color(target: str) -> str:
    cat = _element_category(target)
    return ELEMENT_COLORS.get(cat, "#E53935")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _classify_junction_location(junction_pos: int, gene: dict) -> str:
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
    if abs(bp) >= 1000:
        return f"{bp/1000:.1f} kb"
    return f"{bp} bp"


def _pick_scale_bar(win_width: int) -> int:
    candidates = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000]
    for c in candidates:
        if c >= win_width * 0.08 and c <= win_width * 0.30:
            return c
    target = int(win_width * 0.15)
    return min(candidates, key=lambda c: abs(c - target))


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------

def plot_junction_with_gene(
    junction: dict,
    gene_info: dict,
    contig_seq: str | None,
    sample_name: str,
    output_path: Path,
    construct_hits: list[dict] | None = None,
    host_hits: list[dict] | None = None,
    depth: np.ndarray | None = None,
) -> None:
    """Create junction diagram: gene track + annotated contig + depth."""

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
    genes = sorted(genes, key=lambda g: g["start"])
    if len(genes) > 5:
        genes = sorted(genes, key=lambda g: abs((g["start"]+g["end"])/2 - junction_pos))[:5]
        genes = sorted(genes, key=lambda g: g["start"])

    # --- Window: center on junction ---
    flank = 10000
    all_positions = [junction_pos - flank, junction_pos + flank]
    for g in genes:
        all_positions.extend([g["start"], g["end"]])
    win_start = min(all_positions) - 500
    win_end = max(all_positions) + 500
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

    # --- Decide row heights ---
    has_elements = construct_hits and len(construct_hits) > 0
    has_depth = depth is not None and len(depth) > 0

    n_rows = 3  # gene + contig + legend
    ratios = [2.5, 2.5, 0.3]
    if has_depth:
        n_rows = 4
        ratios = [2.5, 2.5, 1.2, 0.3]

    fig = plt.figure(figsize=(16, 3 + sum(ratios) * 1.1))
    gs = gridspec.GridSpec(n_rows, 1, height_ratios=ratios, hspace=0.08)

    # ==== ROW 1: Gene track ====
    ax_gene = fig.add_subplot(gs[0])
    ax_gene.set_xlim(win_start, win_end)

    gene_y_positions = []
    GENE_HEIGHT = 0.35
    GENE_SPACING = 0.85
    for gi, gene in enumerate(genes):
        y = 0
        for prev_gi, prev_gene in enumerate(genes[:gi]):
            prev_y = gene_y_positions[prev_gi]
            if gene["start"] <= prev_gene["end"] + 200 and gene["end"] >= prev_gene["start"] - 200:
                y = max(y, prev_y + GENE_SPACING)
        gene_y_positions.append(y)

    max_y = max(gene_y_positions) if gene_y_positions else 0
    ax_gene.set_ylim(-0.8, max_y + GENE_SPACING + 0.3)

    chr_y = -0.5
    ax_gene.barh(chr_y, win_width, left=win_start, height=0.12,
                 color="#E0E0E0", edgecolor="#BBB", linewidth=0.5, zorder=0)
    ax_gene.text(win_start + win_width * 0.005, chr_y, host_chr,
                 fontsize=7, va="center", ha="left", color="#666", style="italic")

    for gi, gene in enumerate(genes):
        y = gene_y_positions[gi]
        color = GENE_COLORS[gi % len(GENE_COLORS)]
        gene_start = gene["start"]
        gene_end = gene["end"]
        gene_strand = gene["strand"]
        draw_start = max(gene_start, win_start)
        draw_end = min(gene_end, win_end)
        ax_gene.plot([draw_start, draw_end], [y, y],
                     color=INTRON_COLOR, linewidth=1.2, zorder=1)
        arrow_step = max(500, (draw_end - draw_start) // 15)
        for apos in range(int(draw_start) + 200, int(draw_end) - 200, arrow_step):
            dx = 150 if gene_strand == "+" else -150
            ax_gene.annotate("", xy=(apos + dx, y), xytext=(apos, y),
                             arrowprops=dict(arrowstyle="->", color=INTRON_COLOR, lw=0.6), zorder=1)
        for es, ee in gene["exons"]:
            if ee < win_start or es > win_end:
                continue
            ax_gene.barh(y, min(ee, win_end) - max(es, win_start),
                         left=max(es, win_start), height=GENE_HEIGHT,
                         color=UTR_COLOR, edgecolor=color, linewidth=0.8, zorder=3)
        for cs, ce, frame in gene["cds"]:
            if ce < win_start or cs > win_end:
                continue
            ax_gene.barh(y, min(ce, win_end) - max(cs, win_start),
                         left=max(cs, win_start), height=GENE_HEIGHT,
                         color=color, edgecolor="black", linewidth=0.7, alpha=0.85, zorder=4)
        for us, ue in gene["utrs"]:
            if ue < win_start or us > win_end:
                continue
            ax_gene.barh(y, min(ue, win_end) - max(us, win_start),
                         left=max(us, win_start), height=GENE_HEIGHT * 0.7,
                         color=UTR_COLOR, edgecolor=color, linewidth=0.5, zorder=3)
        location = _classify_junction_location(junction_pos, gene)
        label_text = f"{gene['name']} ({gene_strand})"
        if location != "intergenic":
            label_text += f" [{location}]"
        label_x = max(win_start + win_width * 0.01,
                      min((gene_start + gene_end) / 2, win_end - win_width * 0.15))
        ax_gene.text(label_x, y + GENE_HEIGHT / 2 + 0.08, label_text,
                     fontsize=7.5, va="bottom", ha="center", color=color,
                     fontweight="bold",
                     bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                               edgecolor=color, alpha=0.85, linewidth=0.5))

    ax_gene.axvline(x=junction_pos, color=JUNCTION_COLOR, linewidth=2.5,
                    linestyle="-", zorder=20, alpha=0.9)
    ax_gene.annotate(
        f"Junction: {host_chr}:{junction_pos:,}\n[{junction_type}] MAPQ={host_mapq}",
        xy=(junction_pos, max_y + GENE_SPACING),
        xytext=(junction_pos, max_y + GENE_SPACING + 0.2),
        fontsize=8.5, fontweight="bold", color=JUNCTION_COLOR,
        ha="center", va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF3E0",
                  edgecolor=JUNCTION_COLOR, alpha=0.95))

    scale_len = _pick_scale_bar(win_width)
    sb_x = win_end - scale_len - win_width * 0.02
    sb_y = -0.3
    ax_gene.plot([sb_x, sb_x + scale_len], [sb_y, sb_y], color="black", linewidth=2, zorder=10)
    ax_gene.plot([sb_x, sb_x], [sb_y - 0.05, sb_y + 0.05], color="black", linewidth=1.5, zorder=10)
    ax_gene.plot([sb_x + scale_len, sb_x + scale_len], [sb_y - 0.05, sb_y + 0.05], color="black", linewidth=1.5, zorder=10)
    ax_gene.text(sb_x + scale_len / 2, sb_y - 0.12, _fmt_kb(scale_len),
                 fontsize=7, ha="center", va="top", color="black")

    ax_gene.set_title(
        f"Junction Gene Context — {sample_name}\n"
        f"Construct: {element_name} | Confidence: {confidence}",
        fontsize=12, fontweight="bold", pad=12)
    ax_gene.set_yticks([])
    ax_gene.tick_params(labelbottom=False, bottom=False)
    for spine in ax_gene.spines.values():
        spine.set_visible(False)

    # ==== ROW 2: Full annotated contig ====
    ax_contig = fig.add_subplot(gs[1])
    contig_visual_width = win_width * 0.75
    contig_scale = contig_visual_width / contig_len

    def contig_to_x(cp: float) -> float:
        return junction_pos + (cp - junction_on_contig) * contig_scale

    ax_contig.set_xlim(win_start, win_end)

    # -- Element annotation tracks --
    # Collect all PAF hits for this contig
    all_elements = []
    if construct_hits:
        for hit in construct_hits:
            all_elements.append({
                "start": hit["contig_start"],
                "end": hit["contig_end"],
                "name": _element_short_name(hit["target"]),
                "category": _element_category(hit["target"]),
                "color": _get_element_color(hit["target"]),
                "type": "construct",
                "target_coords": f"{hit['target_start']:,}-{hit['target_end']:,}",
                "mapq": hit["mapq"],
                "aln_len": hit["aln_len"],
            })
    if host_hits:
        for hit in host_hits:
            all_elements.append({
                "start": hit["contig_start"],
                "end": hit["contig_end"],
                "name": hit["target"],
                "category": "host",
                "color": HOST_COLOR,
                "type": "host",
                "target_coords": f"{hit['target_start']:,}-{hit['target_end']:,}",
                "mapq": hit["mapq"],
                "aln_len": hit["aln_len"],
            })

    # Stack elements into rows to avoid overlap
    element_rows = []
    for elem in sorted(all_elements, key=lambda e: e["start"]):
        placed = False
        for row in element_rows:
            if all(elem["start"] >= r["end"] + 5 or elem["end"] <= r["start"] - 5 for r in row):
                row.append(elem)
                placed = True
                break
        if not placed:
            element_rows.append([elem])

    n_elem_rows = max(len(element_rows), 1)
    ELEM_HEIGHT = 0.4
    ELEM_SPACING = 0.55
    contig_bar_y = 0.0
    bar_h = 0.3

    # y-limits: elements above, sequence below
    y_top = contig_bar_y + (n_elem_rows) * ELEM_SPACING + 1.5
    y_bot = contig_bar_y - 2.0
    ax_contig.set_ylim(y_bot, y_top)

    # -- Contig backbone --
    contig_x0 = contig_to_x(0)
    contig_x1 = contig_to_x(contig_len)
    ax_contig.barh(contig_bar_y, contig_x1 - contig_x0, left=contig_x0,
                   height=bar_h, color="#F0F0F0", edgecolor="#999",
                   linewidth=0.8, zorder=1)

    # -- Draw element annotations above contig --
    used_categories = set()
    for ri, row in enumerate(element_rows):
        elem_y = contig_bar_y + bar_h / 2 + 0.2 + ri * ELEM_SPACING
        for elem in row:
            ex0 = contig_to_x(elem["start"])
            ex1 = contig_to_x(elem["end"])
            color = elem["color"]
            used_categories.add((elem["category"], color))

            # Draw the element bar
            ax_contig.barh(elem_y, ex1 - ex0, left=ex0, height=ELEM_HEIGHT,
                           color=color, alpha=0.7, edgecolor="black",
                           linewidth=0.6, zorder=5)

            # Label
            mid_x = (ex0 + ex1) / 2
            label = elem["name"]
            if elem["type"] == "host":
                label = f"{elem['name']}:{elem['target_coords']}"
                if len(label) > 30:
                    label = f"{elem['name'][:15]}..."

            # Check if element bar is wide enough for text
            bar_width_pts = abs(ex1 - ex0) / win_width * 16 * 72  # approx pts
            if bar_width_pts > 30:
                ax_contig.text(mid_x, elem_y, label,
                               ha="center", va="center", fontsize=5.5,
                               fontweight="bold", color="white",
                               clip_on=True, zorder=6)
            else:
                # Label above
                ax_contig.text(mid_x, elem_y + ELEM_HEIGHT / 2 + 0.03, label,
                               ha="center", va="bottom", fontsize=5,
                               color=color, fontweight="bold",
                               rotation=45, rotation_mode="anchor",
                               zorder=6)

    # -- Funnel from gene track to contig --
    junc_x = junction_pos
    funnel_top = y_top - 0.2
    funnel_bot = contig_bar_y + bar_h / 2
    funnel_narrow = win_width * 0.003

    if junction_type in ("LB", "5prime"):
        cx0_f = contig_to_x(0)
        funnel_verts = [
            (junc_x - funnel_narrow, funnel_top),
            (junc_x + funnel_narrow, funnel_top),
            (junc_x, funnel_bot),
            (cx0_f, funnel_bot),
            (junc_x - funnel_narrow, funnel_top),
        ]
    else:
        cx1_f = contig_to_x(contig_len)
        funnel_verts = [
            (junc_x - funnel_narrow, funnel_top),
            (junc_x + funnel_narrow, funnel_top),
            (cx1_f, funnel_bot),
            (junc_x, funnel_bot),
            (junc_x - funnel_narrow, funnel_top),
        ]
    funnel_poly = Polygon(funnel_verts, closed=True,
                          facecolor=CONSTRUCT_COLOR, alpha=0.06,
                          edgecolor=CONSTRUCT_COLOR, linewidth=0.6,
                          linestyle="--", zorder=2)
    ax_contig.add_patch(funnel_poly)

    # Junction line
    ax_contig.axvline(x=junction_pos, color=JUNCTION_COLOR, linewidth=2.5,
                      linestyle="-", zorder=20, alpha=0.9)
    ax_contig.plot(junction_pos, funnel_top, marker="v", markersize=8,
                   color=JUNCTION_COLOR, zorder=25, clip_on=False)

    # -- Sequence around junction --
    if contig_seq:
        jx = min(max(0, junction_on_contig), len(contig_seq) - 1)
        seq_flank = 50
        seq_s = max(0, jx - seq_flank)
        seq_e = min(len(contig_seq), jx + seq_flank)
        host_part = contig_seq[seq_s:jx]
        const_part = contig_seq[jx:seq_e]
        prefix = "..." if seq_s > 0 else ""
        suffix = "..." if seq_e < len(contig_seq) else ""

        seq_y = contig_bar_y - bar_h / 2 - 0.35
        ax_contig.text(junction_pos, seq_y,
                       f"{prefix}{host_part}",
                       ha="right", va="top", fontsize=6,
                       fontfamily="monospace", color=HOST_COLOR,
                       fontweight="bold", zorder=10)
        ax_contig.text(junction_pos, seq_y, "|",
                       ha="center", va="top", fontsize=6,
                       fontfamily="monospace", color=JUNCTION_COLOR,
                       fontweight="bold", zorder=10)
        ax_contig.text(junction_pos, seq_y,
                       f"{const_part}{suffix}",
                       ha="left", va="top", fontsize=6,
                       fontfamily="monospace", color=CONSTRUCT_COLOR,
                       fontweight="bold", zorder=10)
        full_seq = f"{prefix}{host_part}|{const_part}{suffix}"
        ax_contig.text(junction_pos, seq_y, full_seq,
                       ha="center", va="top", fontsize=6,
                       fontfamily="monospace", color="none",
                       bbox=dict(facecolor="white", edgecolor="#999",
                                 pad=3, boxstyle="round,pad=0.3"),
                       zorder=9)

    # Contig info & scale
    ax_contig.text(contig_x0, contig_bar_y - bar_h / 2 - 1.2,
                   f"Contig: {contig_name}  ({_fmt_kb(contig_len)})",
                   fontsize=7, va="top", color="#555", style="italic")

    ctg_sb_len = _pick_scale_bar(int(contig_len))
    ctg_sb_x0 = contig_to_x(contig_len - ctg_sb_len - 5)
    ctg_sb_x1 = contig_to_x(contig_len - 5)
    sb_y2 = contig_bar_y - bar_h / 2 - 1.3
    ax_contig.plot([ctg_sb_x0, ctg_sb_x1], [sb_y2, sb_y2], color="black", linewidth=2, zorder=10)
    ax_contig.plot([ctg_sb_x0, ctg_sb_x0], [sb_y2 - 0.05, sb_y2 + 0.05], color="black", linewidth=1.5, zorder=10)
    ax_contig.plot([ctg_sb_x1, ctg_sb_x1], [sb_y2 - 0.05, sb_y2 + 0.05], color="black", linewidth=1.5, zorder=10)
    ax_contig.text((ctg_sb_x0 + ctg_sb_x1) / 2, sb_y2 - 0.12,
                   f"{_fmt_kb(ctg_sb_len)} (contig)", fontsize=6, ha="center", va="top", color="#555")

    ax_contig.set_yticks([])
    ax_contig.tick_params(labelbottom=False, bottom=False)
    for spine in ax_contig.spines.values():
        spine.set_visible(False)

    # ==== ROW 3 (optional): Read depth ====
    if has_depth:
        ax_depth = fig.add_subplot(gs[2])
        ax_depth.set_xlim(win_start, win_end)

        # Convert contig positions to plot x-coordinates
        x_positions = np.array([contig_to_x(i) for i in range(len(depth))])
        mask = (x_positions >= win_start) & (x_positions <= win_end)

        if mask.any():
            x_vis = x_positions[mask]
            d_vis = depth[mask]

            # Color by host/construct region
            contig_positions = np.arange(len(depth))
            cp_vis = contig_positions[mask]

            # Build color arrays
            host_mask = np.zeros(len(cp_vis), dtype=bool)
            const_mask = np.zeros(len(cp_vis), dtype=bool)
            for hit in (host_hits or []):
                host_mask |= (cp_vis >= hit["contig_start"]) & (cp_vis < hit["contig_end"])
            for hit in (construct_hits or []):
                const_mask |= (cp_vis >= hit["contig_start"]) & (cp_vis < hit["contig_end"])

            # Fill areas
            ax_depth.fill_between(x_vis, 0, d_vis, where=host_mask,
                                  color=HOST_COLOR, alpha=0.4, step="mid", zorder=3)
            ax_depth.fill_between(x_vis, 0, d_vis, where=const_mask,
                                  color=CONSTRUCT_COLOR, alpha=0.4, step="mid", zorder=3)
            ax_depth.fill_between(x_vis, 0, d_vis, where=~(host_mask | const_mask),
                                  color="#999", alpha=0.3, step="mid", zorder=3)

            # Line on top
            ax_depth.plot(x_vis, d_vis, color="#333", linewidth=0.6, zorder=5)

            # Junction line
            ax_depth.axvline(x=junction_pos, color=JUNCTION_COLOR, linewidth=2,
                             linestyle="-", zorder=20, alpha=0.8)

            max_d = max(d_vis.max(), 1)
            ax_depth.set_ylim(0, max_d * 1.15)
            ax_depth.set_ylabel("Depth", fontsize=8)
            ax_depth.tick_params(labelsize=7)

            # Mean depth annotation
            mean_d = np.mean(d_vis[d_vis > 0]) if np.any(d_vis > 0) else 0
            ax_depth.text(win_end - win_width * 0.01, max_d * 1.05,
                          f"mean={mean_d:.0f}x  max={max_d}x",
                          ha="right", va="top", fontsize=7, color="#555")

        ax_depth.tick_params(labelbottom=False, bottom=False)
        for spine in ["top", "right", "bottom"]:
            ax_depth.spines[spine].set_visible(False)

    # ==== Last ROW: Legend ====
    ax_leg = fig.add_subplot(gs[-1])
    ax_leg.axis("off")

    legend_patches = [
        mpatches.Patch(color=HOST_COLOR, alpha=0.7, label="Host genome"),
        mpatches.Patch(color=CDS_COLOR, alpha=0.85, label="CDS"),
        mpatches.Patch(color=UTR_COLOR, label="Exon/UTR"),
        plt.Line2D([0], [0], color=INTRON_COLOR, linewidth=1.5, label="Intron"),
        plt.Line2D([0], [0], color=JUNCTION_COLOR, linewidth=2.5, label="Junction site"),
    ]
    # Add element category colors used
    for cat, color in sorted(used_categories):
        if cat != "host":
            legend_patches.append(mpatches.Patch(color=color, alpha=0.7, label=cat))

    ax_leg.legend(handles=legend_patches, loc="center",
                  ncol=min(len(legend_patches), 7),
                  fontsize=8, framealpha=0.9, edgecolor="#ccc")

    plt.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    log(f"Saved: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Junction visualization with gene/exon/CDS + element map + depth",
    )
    parser.add_argument("--junctions", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--contigs", type=Path, default=None)
    parser.add_argument("--construct-paf", type=Path, default=None,
                        help="PAF: contigs aligned to construct/element DB")
    parser.add_argument("--host-paf", type=Path, default=None,
                        help="PAF: contigs aligned to host genome")
    parser.add_argument("--reads-r1", type=Path, default=None,
                        help="Extracted reads R1 for depth calculation")
    parser.add_argument("--reads-r2", type=Path, default=None,
                        help="Extracted reads R2 for depth calculation")
    parser.add_argument("--sample-name", type=str, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=10000)
    parser.add_argument("--confidence-filter", type=str, default=None)
    args = parser.parse_args()

    out_dir = args.outdir / args.sample_name / "s06_junction"
    out_dir.mkdir(parents=True, exist_ok=True)

    junctions = []
    with open(args.junctions) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            junctions.append(row)
    if not junctions:
        log("No junctions found")
        return

    conf_order = {"Low": 0, "Medium": 1, "High": 2}
    if args.confidence_filter:
        min_conf = conf_order.get(args.confidence_filter, 0)
        junctions = [j for j in junctions
                     if conf_order.get(j.get("confidence", "Low"), 0) >= min_conf]

    contig_seqs = {}
    if args.contigs and args.contigs.exists():
        contig_seqs = read_contigs(args.contigs)
        log(f"Loaded {len(contig_seqs)} contigs")

    # Collect unique junction contigs for depth computation
    junction_contigs = set(j["contig_name"] for j in junctions)
    contig_depths = {}

    if args.reads_r1 and args.reads_r1.exists():
        for cname in junction_contigs:
            if cname in contig_seqs:
                log(f"  Computing depth for {cname}...")
                d = compute_contig_depth(cname, contig_seqs[cname],
                                         args.reads_r1, args.reads_r2)
                if d is not None:
                    contig_depths[cname] = d
                    log(f"    mean={np.mean(d):.0f}x, max={np.max(d)}x")

    log(f"Processing {len(junctions)} junctions")

    for i, junc in enumerate(junctions):
        host_chr = junc["host_chr"]
        junction_pos = int(junc["junction_pos_host"])
        contig_name = junc["contig_name"]

        log(f"  Junction {i+1}: {host_chr}:{junction_pos}")

        gene_info = parse_gff3_region(
            args.gff, host_chr, junction_pos, junction_pos, flank=args.flank)
        log(f"    Found {len(gene_info.get('genes', []))} gene(s)")

        contig_seq = contig_seqs.get(contig_name)

        construct_hits = None
        if args.construct_paf and args.construct_paf.exists():
            construct_hits = parse_paf_for_contig(args.construct_paf, contig_name)

        host_hits = None
        if args.host_paf and args.host_paf.exists():
            host_hits = parse_paf_for_contig(args.host_paf, contig_name)

        depth = contig_depths.get(contig_name)

        for ext in ("png", "pdf"):
            out_path = out_dir / f"junction_gene_{i+1}_{host_chr}_{junction_pos}.{ext}"
            plot_junction_with_gene(
                junction=junc,
                gene_info=gene_info,
                contig_seq=contig_seq,
                sample_name=args.sample_name,
                output_path=out_path,
                construct_hits=construct_hits,
                host_hits=host_hits,
                depth=depth,
            )

    log("All junction gene context plots generated.")


if __name__ == "__main__":
    main()
