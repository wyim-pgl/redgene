#!/usr/bin/env python3
"""Junction insertion site multi-track visualization.

Generates a genome browser-style multi-track figure for each detected
T-DNA insertion junction, inspired by T-DNAreader (Genome Biol 2025)
but implemented with pure matplotlib (no pyGenomeTracks dependency).

Tracks (top to bottom):
  1. Genomic coordinate axis with scale bar
  2. Gene model track (exons, CDS, UTR from GFF3)
  3. Host read depth track (from BAM, if available)
  4. T-DNA insertion site marker with orientation arrow
  5. Chimeric contig structure (host region vs construct region)
  6. Junction reads / soft-clipped reads (from BAM, if available)
  7. Amplicon primer locations (if amplicon TSV available)

Output:
  {outdir}/junction_track_{sample}_{contig}.png — one per junction contig
"""

from __future__ import annotations

import argparse
import gzip
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import os
os.environ["MPLBACKEND"] = "Agg"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import matplotlib.gridspec as gridspec
import numpy as np


# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------
COLORS = {
    "host": "#3B82F6",        # blue
    "construct": "#EF4444",   # red
    "junction": "#F59E0B",    # amber/gold
    "gene": "#10B981",        # emerald
    "exon": "#059669",        # darker emerald
    "cds": "#047857",         # even darker
    "utr": "#6EE7B7",        # light green
    "depth_fill": "#93C5FD",  # light blue
    "depth_line": "#2563EB",  # darker blue
    "fwd_primer": "#16A34A",  # green
    "rev_primer": "#EA580C",  # orange
    "read_host": "#60A5FA",   # lighter blue for host-aligned reads
    "read_clip": "#F87171",   # lighter red for clipped portion
    "contig_bg": "#F3F4F6",   # light gray
    "scale_bar": "#374151",   # dark gray
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def parse_junctions(tsv_path: Path) -> list[dict]:
    """Parse junctions.tsv."""
    junctions = []
    with open(tsv_path) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) >= len(header):
                junctions.append(dict(zip(header, fields)))
    return junctions


def load_contigs(fasta_path: Path) -> dict[str, str]:
    """Load contig sequences from FASTA."""
    contigs: dict[str, str] = {}
    name = ""
    parts: list[str] = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    contigs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name:
        contigs[name] = "".join(parts)
    return contigs


def load_paf(paf_path: Path) -> dict[str, list[dict]]:
    """Load PAF file into structured records by query name."""
    records: dict[str, list[dict]] = defaultdict(list)
    with open(paf_path) as fh:
        for line in fh:
            f = line.strip().split("\t")
            if len(f) < 12:
                continue
            rec = {
                "qname": f[0], "qlen": int(f[1]),
                "qstart": int(f[2]), "qend": int(f[3]),
                "strand": f[4],
                "tname": f[5], "tlen": int(f[6]),
                "tstart": int(f[7]), "tend": int(f[8]),
                "matches": int(f[9]), "alnlen": int(f[10]),
                "mapq": int(f[11]),
            }
            records[f[0]].append(rec)
    return records


def parse_gff3(gff_path: Path, region_chr: str, region_start: int, region_end: int) -> list[dict]:
    """Parse GFF3 for genes/mRNA/exon/CDS overlapping a region."""
    features = []
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.strip().split("\t")
            if len(f) < 9:
                continue
            chrom, source, ftype, start, end = f[0], f[1], f[2], int(f[3]), int(f[4])
            strand = f[6]
            if chrom != region_chr:
                continue
            if end < region_start or start > region_end:
                continue
            if ftype not in ("gene", "mRNA", "exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
                continue
            attrs = {}
            for kv in f[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v
            features.append({
                "chrom": chrom, "type": ftype, "start": start, "end": end,
                "strand": strand, "attrs": attrs,
            })
    return features


def load_depth_from_bam(bam_path: Path, region_chr: str, region_start: int,
                        region_end: int, bin_size: int = 100) -> tuple[np.ndarray, np.ndarray]:
    """Get read depth from BAM using samtools depth."""
    import subprocess
    region = f"{region_chr}:{region_start}-{region_end}"
    cmd = ["samtools", "depth", "-a", "-r", region, str(bam_path)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            return np.array([]), np.array([])
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return np.array([]), np.array([])

    positions = []
    depths = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            positions.append(int(parts[1]))
            depths.append(int(parts[2]))

    if not positions:
        return np.array([]), np.array([])

    pos = np.array(positions)
    dep = np.array(depths)

    # Bin for smoother display
    if bin_size > 1 and len(pos) > bin_size * 2:
        n_bins = max(1, (region_end - region_start) // bin_size)
        bin_edges = np.linspace(region_start, region_end, n_bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_depths = np.zeros(n_bins)
        indices = np.digitize(pos, bin_edges) - 1
        for i in range(n_bins):
            mask = indices == i
            if mask.any():
                bin_depths[i] = dep[mask].mean()
        return bin_centers, bin_depths

    return pos, dep.astype(float)


def load_amplicons(tsv_path: Path, sample: str) -> list[dict]:
    """Load amplicon primer data."""
    amplicons = []
    with open(tsv_path) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) >= len(header):
                amplicons.append(dict(zip(header, fields)))
    return amplicons


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def format_position(pos: int) -> str:
    """Format genomic position for display."""
    if pos >= 1_000_000:
        return f"{pos/1_000_000:.2f} Mb"
    elif pos >= 1_000:
        return f"{pos/1_000:.1f} kb"
    return str(pos)


def draw_scale_bar(ax, start: int, end: int):
    """Draw genomic coordinate scale bar."""
    span = end - start
    # Choose appropriate tick spacing
    if span > 100_000:
        tick_spacing = 50_000
    elif span > 10_000:
        tick_spacing = 5_000
    elif span > 1_000:
        tick_spacing = 500
    else:
        tick_spacing = 100

    first_tick = ((start // tick_spacing) + 1) * tick_spacing
    ticks = list(range(first_tick, end, tick_spacing))

    ax.set_xlim(start, end)
    ax.set_xticks(ticks)
    ax.set_xticklabels([format_position(t) for t in ticks], fontsize=7, rotation=30, ha="right")
    ax.tick_params(axis="x", length=3, width=0.5, pad=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_yticks([])

    # Scale bar
    bar_len = tick_spacing
    bar_x = start + span * 0.02
    ax.plot([bar_x, bar_x + bar_len], [0.5, 0.5], color=COLORS["scale_bar"], lw=2)
    ax.text(bar_x + bar_len / 2, 0.7, format_position(bar_len),
            ha="center", va="bottom", fontsize=7, color=COLORS["scale_bar"])
    ax.set_ylim(0, 1)


def draw_gene_track(ax, features: list[dict], region_start: int, region_end: int):
    """Draw gene models with exons and CDS."""
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-0.5, 2.5)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    genes = [f for f in features if f["type"] == "gene"]
    if not genes:
        ax.text((region_start + region_end) / 2, 1, "No annotated genes in region",
                ha="center", va="center", fontsize=8, color="#9CA3AF", style="italic")
        return

    for gene in genes:
        gs = max(gene["start"], region_start)
        ge = min(gene["end"], region_end)
        y = 1.0

        # Gene body line
        ax.plot([gs, ge], [y, y], color=COLORS["gene"], lw=1.5, solid_capstyle="butt")

        # Gene name
        name = gene["attrs"].get("Name", gene["attrs"].get("ID", ""))
        if name:
            ax.text((gs + ge) / 2, y + 0.6, name, ha="center", va="bottom",
                    fontsize=7, color=COLORS["gene"], fontweight="bold")

        # Strand arrow
        if gene["strand"] == "+":
            ax.annotate("", xy=(ge, y), xytext=(ge - (ge - gs) * 0.1, y),
                        arrowprops=dict(arrowstyle="->", color=COLORS["gene"], lw=1.5))
        elif gene["strand"] == "-":
            ax.annotate("", xy=(gs, y), xytext=(gs + (ge - gs) * 0.1, y),
                        arrowprops=dict(arrowstyle="->", color=COLORS["gene"], lw=1.5))

    # Draw exons
    exons = [f for f in features if f["type"] == "exon"]
    for exon in exons:
        es = max(exon["start"], region_start)
        ee = min(exon["end"], region_end)
        rect = Rectangle((es, 0.7), ee - es, 0.6, facecolor=COLORS["exon"],
                          edgecolor="none", alpha=0.7)
        ax.add_patch(rect)

    # Draw CDS (darker, thicker)
    cdss = [f for f in features if f["type"] == "CDS"]
    for cds in cdss:
        cs = max(cds["start"], region_start)
        ce = min(cds["end"], region_end)
        rect = Rectangle((cs, 0.6), ce - cs, 0.8, facecolor=COLORS["cds"],
                          edgecolor="black", linewidth=0.3, alpha=0.85)
        ax.add_patch(rect)

    ax.set_ylabel("Genes", fontsize=8, rotation=0, ha="right", va="center")


def draw_depth_track(ax, positions: np.ndarray, depths: np.ndarray,
                     region_start: int, region_end: int, junction_pos: int):
    """Draw read depth coverage track."""
    ax.set_xlim(region_start, region_end)

    if len(positions) == 0:
        ax.text((region_start + region_end) / 2, 0.5, "No BAM available",
                ha="center", va="center", fontsize=8, color="#9CA3AF", style="italic",
                transform=ax.get_xaxis_transform())
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        return

    ax.fill_between(positions, depths, alpha=0.4, color=COLORS["depth_fill"])
    ax.plot(positions, depths, color=COLORS["depth_line"], lw=0.8)

    # Junction line
    ax.axvline(x=junction_pos, color=COLORS["junction"], lw=1.5, ls="--", alpha=0.8)

    max_depth = max(depths) if len(depths) > 0 else 1
    ax.set_ylim(0, max_depth * 1.15)
    ax.set_ylabel("Depth", fontsize=8, rotation=0, ha="right", va="center")
    ax.tick_params(axis="y", labelsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_junction_marker(ax, junction_pos: int, junction_type: str, confidence: str,
                         host_chr: str, region_start: int, region_end: int):
    """Draw T-DNA insertion site marker."""
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Junction marker — diamond
    ax.plot(junction_pos, 0.5, marker="D", markersize=12, color=COLORS["junction"],
            markeredgecolor="black", markeredgewidth=0.8, zorder=5)

    # Label
    label = f"{junction_type} junction\n{host_chr}:{junction_pos:,}\n[{confidence}]"
    ax.text(junction_pos, 0.5, f"  {label}", va="center", ha="left",
            fontsize=7, fontweight="bold", color="#333333",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#FEF3C7", edgecolor=COLORS["junction"], alpha=0.9))

    ax.set_ylabel("T-DNA\nsite", fontsize=7, rotation=0, ha="right", va="center")


def draw_contig_track(ax, contig_name: str, contig_seq: str,
                      host_alns: list[dict], construct_alns: list[dict],
                      region_start: int, region_end: int, host_chr: str):
    """Draw chimeric contig structure showing host vs construct regions."""
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-0.5, 1.5)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Find the host alignment for this junction
    target_alns = [a for a in host_alns if a["tname"] == host_chr]
    if not target_alns:
        return

    best_h = max(target_alns, key=lambda a: a["qend"] - a["qstart"])
    best_c = max(construct_alns, key=lambda a: a["qend"] - a["qstart"]) if construct_alns else None

    if best_c is None:
        return

    # Map contig coordinates to genomic coordinates
    # The contig spans from host_start to host_end in genomic space
    h_qstart, h_qend = best_h["qstart"], best_h["qend"]
    h_tstart, h_tend = best_h["tstart"], best_h["tend"]
    c_qstart, c_qend = best_c["qstart"], best_c["qend"]

    qlen = best_h["qlen"]

    # Scale: genomic bp per contig bp
    if best_h["strand"] == "+":
        bp_per_q = (h_tend - h_tstart) / max(1, h_qend - h_qstart)
        contig_genome_start = h_tstart - h_qstart * bp_per_q
    else:
        bp_per_q = (h_tend - h_tstart) / max(1, h_qend - h_qstart)
        contig_genome_start = h_tend + h_qstart * bp_per_q - qlen * bp_per_q

    # Draw contig background
    cg_start = contig_genome_start
    cg_end = contig_genome_start + qlen * bp_per_q

    # Clamp to visible region
    vis_start = max(region_start, cg_start)
    vis_end = min(region_end, cg_end)

    if vis_end <= vis_start:
        return

    # Background
    rect = Rectangle((vis_start, 0.1), vis_end - vis_start, 0.8,
                      facecolor=COLORS["contig_bg"], edgecolor="#D1D5DB", lw=0.5)
    ax.add_patch(rect)

    # Host region (blue)
    h_gstart = contig_genome_start + h_qstart * bp_per_q
    h_gend = contig_genome_start + h_qend * bp_per_q
    h_vis_start = max(region_start, min(h_gstart, h_gend))
    h_vis_end = min(region_end, max(h_gstart, h_gend))
    rect_h = Rectangle((h_vis_start, 0.15), h_vis_end - h_vis_start, 0.7,
                         facecolor=COLORS["host"], edgecolor="none", alpha=0.7)
    ax.add_patch(rect_h)

    # Construct region (red)
    c_gstart = contig_genome_start + c_qstart * bp_per_q
    c_gend = contig_genome_start + c_qend * bp_per_q
    c_vis_start = max(region_start, min(c_gstart, c_gend))
    c_vis_end = min(region_end, max(c_gstart, c_gend))
    rect_c = Rectangle((c_vis_start, 0.15), c_vis_end - c_vis_start, 0.7,
                         facecolor=COLORS["construct"], edgecolor="none", alpha=0.7)
    ax.add_patch(rect_c)

    # Contig label
    short_name = contig_name.split("_length")[0] if "_length" in contig_name else contig_name
    ax.text(vis_start + 5, 1.2, f"{short_name} ({qlen} bp)", fontsize=6.5,
            color="#4B5563", fontweight="bold")

    # Element name
    if construct_alns:
        elem = best_c["tname"]
        short_elem = elem.split("|")[1] if "|" in elem else elem[:30]
        ax.text((c_vis_start + c_vis_end) / 2, 0.5, short_elem,
                ha="center", va="center", fontsize=6, color="white", fontweight="bold")

    ax.set_ylabel("Contig", fontsize=8, rotation=0, ha="right", va="center")


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------

def plot_junction_track(
    junction: dict,
    contig_seq: str,
    host_alns: list[dict],
    construct_alns: list[dict],
    gff_features: list[dict],
    depth_pos: np.ndarray,
    depth_vals: np.ndarray,
    amplicons: list[dict],
    sample_name: str,
    outdir: Path,
    flank: int = 5000,
) -> Path:
    """Create multi-track figure for one junction."""

    junction_pos = int(junction["junction_pos_host"])
    host_chr = junction["host_chr"]
    host_start = int(junction["host_start"])
    host_end = int(junction["host_end"])
    contig_name = junction["contig_name"]
    junction_type = junction.get("junction_type", "?")
    confidence = junction.get("confidence", "?")

    # Region to display
    region_start = max(0, junction_pos - flank)
    region_end = junction_pos + flank

    # Determine number of tracks
    has_depth = len(depth_pos) > 0
    has_genes = len(gff_features) > 0
    has_amplicon = len(amplicons) > 0

    n_tracks = 4  # scale, junction marker, contig, amplicon/info
    if has_genes:
        n_tracks += 1
    if has_depth:
        n_tracks += 1

    # Height ratios
    heights = []
    track_labels = []

    heights.append(0.8)  # scale bar
    track_labels.append("scale")
    if has_genes:
        heights.append(1.5)
        track_labels.append("genes")
    if has_depth:
        heights.append(2.0)
        track_labels.append("depth")
    heights.append(0.8)  # junction marker
    track_labels.append("junction")
    heights.append(1.5)  # contig structure
    track_labels.append("contig")
    if has_amplicon:
        heights.append(0.8)
        track_labels.append("amplicon")

    total_height = sum(heights) + 1.5
    fig = plt.figure(figsize=(14, total_height))
    gs = gridspec.GridSpec(len(heights), 1, height_ratios=heights, hspace=0.15)

    # Title
    construct_elem = junction.get("construct_element", "unknown")
    short_elem = construct_elem.split("|")[1] if "|" in construct_elem else construct_elem[:40]
    fig.suptitle(
        f"{sample_name}  —  {host_chr}:{junction_pos:,}  ({junction_type}, {confidence})\n"
        f"Construct: {short_elem}",
        fontsize=12, fontweight="bold", y=0.98,
    )

    idx = 0

    # Track 1: Scale bar
    ax_scale = fig.add_subplot(gs[idx])
    draw_scale_bar(ax_scale, region_start, region_end)
    idx += 1

    # Track 2: Gene models
    if has_genes:
        ax_genes = fig.add_subplot(gs[idx])
        draw_gene_track(ax_genes, gff_features, region_start, region_end)
        ax_genes.set_xlim(region_start, region_end)
        ax_genes.set_xticklabels([])
        idx += 1

    # Track 3: Read depth
    if has_depth:
        ax_depth = fig.add_subplot(gs[idx])
        draw_depth_track(ax_depth, depth_pos, depth_vals, region_start, region_end, junction_pos)
        ax_depth.set_xlim(region_start, region_end)
        ax_depth.set_xticklabels([])
        idx += 1

    # Track 4: Junction marker
    ax_junc = fig.add_subplot(gs[idx])
    draw_junction_marker(ax_junc, junction_pos, junction_type, confidence,
                         host_chr, region_start, region_end)
    idx += 1

    # Track 5: Contig structure
    ax_contig = fig.add_subplot(gs[idx])
    draw_contig_track(ax_contig, contig_name, contig_seq,
                      host_alns, construct_alns, region_start, region_end, host_chr)
    idx += 1

    # Track 6: Amplicon primers
    if has_amplicon:
        ax_amp = fig.add_subplot(gs[idx])
        ax_amp.set_xlim(region_start, region_end)
        ax_amp.set_ylim(0, 1)
        ax_amp.set_yticks([])
        for spine in ax_amp.spines.values():
            spine.set_visible(False)
        for amp in amplicons:
            # Show primer arrows at junction location
            ax_amp.annotate("", xy=(junction_pos - 50, 0.5),
                           xytext=(junction_pos - 200, 0.5),
                           arrowprops=dict(arrowstyle="->", color=COLORS["fwd_primer"], lw=2))
            ax_amp.annotate("", xy=(junction_pos + 50, 0.5),
                           xytext=(junction_pos + 200, 0.5),
                           arrowprops=dict(arrowstyle="->", color=COLORS["rev_primer"], lw=2))
            ax_amp.text(junction_pos, 0.15,
                       f"Amplicon: {amp.get('amplicon_len', '?')} bp",
                       ha="center", fontsize=7, color="#374151")
        ax_amp.set_ylabel("PCR", fontsize=7, rotation=0, ha="right", va="center")
        idx += 1

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLORS["host"], alpha=0.7, label="Host genome"),
        mpatches.Patch(facecolor=COLORS["construct"], alpha=0.7, label="Construct/T-DNA"),
        plt.Line2D([0], [0], color=COLORS["junction"], lw=1.5, ls="--", label="Junction breakpoint"),
    ]
    if has_genes:
        legend_elements.append(mpatches.Patch(facecolor=COLORS["cds"], label="CDS"))
    if has_depth:
        legend_elements.append(mpatches.Patch(facecolor=COLORS["depth_fill"], alpha=0.4, label="Read depth"))
    if has_amplicon:
        legend_elements.extend([
            plt.Line2D([0], [0], color=COLORS["fwd_primer"], lw=2, label="Fwd primer"),
            plt.Line2D([0], [0], color=COLORS["rev_primer"], lw=2, label="Rev primer"),
        ])

    fig.legend(handles=legend_elements, loc="lower center", ncol=len(legend_elements),
               fontsize=7, frameon=True, fancybox=True)

    fig.tight_layout(rect=[0.05, 0.04, 1, 0.94])

    # Save
    short_contig = contig_name.split("_length")[0].replace("_", "") if "_length" in contig_name else contig_name[:20]
    outpath = outdir / f"junction_track_{sample_name}_{short_contig}.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[track] Saved: {outpath}", file=sys.stderr)
    return outpath


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Multi-track junction visualization (genome browser style)")
    parser.add_argument("--junctions", required=True, type=Path)
    parser.add_argument("--contigs", required=True, type=Path)
    parser.add_argument("--host-paf", required=True, type=Path)
    parser.add_argument("--construct-paf", required=True, type=Path)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--gff", type=Path, default=None, help="GFF3 for gene models")
    parser.add_argument("--host-bam", type=Path, default=None, help="Host-aligned BAM for depth")
    parser.add_argument("--amplicons", type=Path, default=None, help="Amplicon TSV")
    parser.add_argument("--flank", type=int, default=5000, help="Flanking region in bp (default: 5000)")
    parser.add_argument("--min-confidence", default="Medium")
    args = parser.parse_args()

    confidence_order = {"Low": 0, "Medium": 1, "High": 2}
    min_conf = confidence_order.get(args.min_confidence, 1)

    junctions = parse_junctions(args.junctions)
    if not junctions:
        print("[track] No junctions found", file=sys.stderr)
        return

    contigs = load_contigs(args.contigs)
    host_paf = load_paf(args.host_paf)
    construct_paf = load_paf(args.construct_paf)

    # Load amplicons if available
    amplicon_data = []
    if args.amplicons and args.amplicons.exists():
        amplicon_data = load_amplicons(args.amplicons, args.sample_name)

    # Deduplicate: unique (contig, host_chr, junction_pos)
    seen = set()
    unique_junctions = []
    for j in junctions:
        key = (j["contig_name"], j["host_chr"], j["junction_pos_host"])
        conf = confidence_order.get(j.get("confidence", "Medium"), 1)
        if conf < min_conf:
            continue
        if key not in seen:
            seen.add(key)
            unique_junctions.append(j)

    print(f"[track] Generating tracks for {len(unique_junctions)} junctions", file=sys.stderr)

    args.outdir.mkdir(parents=True, exist_ok=True)

    for j in unique_junctions:
        cname = j["contig_name"]
        if cname not in contigs:
            continue

        junction_pos = int(j["junction_pos_host"])
        host_chr = j["host_chr"]
        region_start = max(0, junction_pos - args.flank)
        region_end = junction_pos + args.flank

        # Load GFF features
        gff_features = []
        if args.gff and args.gff.exists():
            gff_features = parse_gff3(args.gff, host_chr, region_start, region_end)

        # Load depth
        depth_pos, depth_vals = np.array([]), np.array([])
        if args.host_bam and args.host_bam.exists():
            depth_pos, depth_vals = load_depth_from_bam(
                args.host_bam, host_chr, region_start, region_end)

        # Filter amplicons for this junction
        matching_amps = [a for a in amplicon_data
                        if a.get("host_chr") == host_chr
                        and abs(int(a.get("junction_pos", 0)) - junction_pos) < 100]

        plot_junction_track(
            junction=j,
            contig_seq=contigs[cname],
            host_alns=host_paf.get(cname, []),
            construct_alns=construct_paf.get(cname, []),
            gff_features=gff_features,
            depth_pos=depth_pos,
            depth_vals=depth_vals,
            amplicons=matching_amps,
            sample_name=args.sample_name,
            outdir=args.outdir,
            flank=args.flank,
        )

    print(f"[track] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
