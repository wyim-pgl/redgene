#!/usr/bin/env python3
"""Generate a multi-panel sample summary figure for RedGene results.

Produces a publication-quality, single-page figure consolidating key results
from one sample's s09 targeted assembly run.  Panels:

  A  Chromosome ideogram with insertion sites marked
  B  Linear element map for each CANDIDATE insert
  C  Tier classification bar chart (transgene-positive vs negative, verdicts)
  D  Assembly convergence summary (rounds and final length per site)
  E  Element composition donut chart
  F  Coverage depth profile placeholder (or insert length waterfall)

Usage:
  python scripts/viz/plot_sample_summary.py \\
    --sample-name rice_G281 \\
    --insert-dir results/rice_G281/s05_insert_assembly \\
    --outdir results \\
    [--host-ref db/Osativa_323_v7.0.fa] \\
    [--gff db/Osativa_323_v7.0.gene_exons.gff3]
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import numpy as np


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def log(msg: str) -> None:
    print(f"[plot_sample_summary] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Colorblind-safe palette (Wong 2011, Nature Methods)
# ---------------------------------------------------------------------------

WONG = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "black": "#000000",
}

ELEMENT_COLORS: dict[str, str] = {
    "promoter": "#009E73",
    "cds": "#0072B2",
    "terminator": "#E69F00",
    "marker": "#CC79A7",
    "border": "#D55E00",
    "vector": "#56B4E9",
    "element-specific": "#F0E442",
    "event-specific": "#999999",
    "construct-specific": "#8B6914",  # dark goldenrod, visible on gray backbone
    "host-endogenous": "#CCCCCC",
    "unknown": "#FFFFFF",
}

VERDICT_COLORS: dict[str, str] = {
    "CANDIDATE": WONG["green"],
    "FALSE_POSITIVE": WONG["vermillion"],
    "TRANSGENE-NEGATIVE": "#BBBBBB",
}


# ---------------------------------------------------------------------------
# Global matplotlib style
# ---------------------------------------------------------------------------

def apply_style() -> None:
    """Set publication-quality defaults."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 8,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "legend.fontsize": 7,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 300,
        "pdf.fonttype": 42,  # TrueType in PDF (editable text)
        "ps.fonttype": 42,
    })


# ---------------------------------------------------------------------------
# Element classification (shared with plot_insert_structure.py)
# ---------------------------------------------------------------------------

def classify_element(element_name: str) -> str:
    """Classify a BLAST element hit into a display category."""
    nl = element_name.lower()
    if any(k in nl for k in ["promoter", "p-35s", "p-e35s", "p-fmv", "p-nos",
                              "p-ocs", "p-ubi", "p-act1", "p-ta29", "p-csvmv",
                              "p-ssuara", "p-cdpk", "p-rice_actin", "p-fmv34s"]):
        return "promoter"
    if any(k in nl for k in ["terminator", "t-nos", "t-35s", "t-e9", "t-ocs",
                              "t-pcambia", "t-pinii", "nos_ter", "polya"]):
        return "terminator"
    if any(k in nl for k in ["lb-tdna", "rb-tdna", "left_border",
                              "right_border"]):
        return "border"
    if nl.endswith("-lb") or nl.endswith("-rb"):
        return "border"
    if any(k in nl for k in ["nptii", "bar", "pat", "hpt", "aada",
                              "cp4-epsps", "epsps", "cry1", "cry2", "cry3",
                              "gus", "barnase", "barstar", "thaumatin",
                              "gat", "dsred", "mepsps", "cat-f", "rerio"]):
        return "cds"
    if any(k in nl for k in ["intron", "i-hsp", "ivs", "enhancer", "aint",
                              "adh1"]):
        return "marker"
    if "ql-con" in nl or "qt-con" in nl or "construct-specific" in nl:
        return "construct-specific"
    if "ql-ele" in nl or "qt-ele" in nl or "element-specific" in nl:
        return "element-specific"
    if "ql-eve" in nl or "qt-eve" in nl or "event-specific" in nl:
        return "event-specific"
    if "univec" in nl or "vector" in nl:
        return "vector"
    # If element name has no pipes (simple name, likely the construct ref)
    # and is short, treat as construct-specific
    if "|" not in element_name and len(element_name) < 60:
        return "construct-specific"
    return "unknown"


def short_element_name(element_name: str) -> str:
    """Extract a short display name from a full element identifier."""
    parts = element_name.split("|")
    if len(parts) >= 3:
        return parts[1]
    if "_" in element_name:
        return element_name.split("_")[0]
    return element_name


# ---------------------------------------------------------------------------
# Data parsing
# ---------------------------------------------------------------------------

def parse_stats(path: Path) -> dict[str, str]:
    """Parse s05_stats.txt key-value pairs."""
    stats: dict[str, str] = {}
    if not path.exists():
        log(f"Stats file not found: {path}")
        return stats
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2:
                stats[parts[0]] = parts[1]
    return stats


def extract_site_data(stats: dict[str, str]) -> list[dict[str, Any]]:
    """Extract per-site data from s05_stats.txt into a list of dicts."""
    sites: list[dict[str, Any]] = []
    # Discover site IDs from keys like insertion_XXXXX_insert_length
    seen_ids: set[str] = set()
    for key in stats:
        m = re.match(r"(insertion_\d+)_insert_length", key)
        if m:
            seen_ids.add(m.group(1))

    for site_id in sorted(seen_ids):
        sites.append({
            "site_id": site_id,
            "insert_length": int(stats.get(f"{site_id}_insert_length", "0")),
            "remaining_ns": int(stats.get(f"{site_id}_remaining_ns", "0")),
            "assembly_rounds": int(stats.get(f"{site_id}_assembly_rounds", "0")),
            "status": stats.get(f"{site_id}_status", "unknown"),
            "verdict": stats.get(f"{site_id}_verdict", "unknown"),
        })
    return sites


def parse_tier_classification(path: Path) -> list[dict[str, Any]]:
    """Parse site_tier_classification.tsv."""
    rows: list[dict[str, Any]] = []
    if not path.exists():
        log(f"Tier classification file not found: {path}")
        return rows
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append({
                "site_id": row["site_id"],
                "chrom": row["chrom"],
                "pos": int(row["pos"]),
                "transgene_positive": row["transgene_positive"] == "True",
                "hit_5p": row.get("hit_5p", ""),
                "hit_3p": row.get("hit_3p", ""),
                "hit_5p_source": row.get("hit_5p_source", ""),
                "hit_3p_source": row.get("hit_3p_source", ""),
            })
    return rows


def parse_element_annotation(path: Path) -> list[dict[str, Any]]:
    """Parse element_annotation.tsv from s09."""
    hits: list[dict[str, Any]] = []
    if not path.exists():
        log(f"Element annotation file not found: {path}")
        return hits
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            site_id = row["query"].replace("_assembled_insert", "")
            hits.append({
                "site_id": site_id,
                "element": row["element"],
                "identity": float(row["identity"]),
                "length": int(row["length"]),
                "q_start": int(row["q_start"]),
                "q_end": int(row["q_end"]),
                "category": classify_element(row["element"]),
                "short_name": short_element_name(row["element"]),
                "source": row.get("source", ""),
            })
    return hits


def parse_report(path: Path) -> dict[str, Any]:
    """Parse an insertion_XXXXX_report.txt for structured info."""
    data: dict[str, Any] = {
        "elements": [],
        "insertion_site": "",
        "structure": "",
        "linear_map": "",
        "foreign_elements": [],
        "host_endogenous_elements": [],
        "borders_found": 0,
    }
    if not path.exists():
        return data
    with open(path) as fh:
        in_elements = False
        in_foreign = False
        in_host_endo = False
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("Insertion site:"):
                data["insertion_site"] = line.split(":", 1)[1].strip()
            elif line.startswith("Structure:"):
                data["structure"] = line.split(":", 1)[1].strip()
            elif "Linear Map" in line:
                pass
            elif "host --[" in line or "host --" in line:
                data["linear_map"] = line.strip()
            elif line.startswith("T-DNA borders found:"):
                try:
                    data["borders_found"] = int(line.split(":")[1].strip())
                except ValueError:
                    pass
            elif line.startswith("Foreign elements"):
                in_foreign = True
                in_host_endo = False
                in_elements = False
            elif line.startswith("Host-endogenous"):
                in_host_endo = True
                in_foreign = False
                in_elements = False
            elif line.startswith("Detailed element positions"):
                in_elements = True
                in_foreign = False
                in_host_endo = False
            elif in_foreign and line.startswith("  +"):
                data["foreign_elements"].append(line.strip().lstrip("+ "))
            elif in_host_endo and line.startswith("  -"):
                data["host_endogenous_elements"].append(line.strip().lstrip("- "))
            elif in_elements and re.match(r"\s+\d+", line):
                cols = line.split()
                if len(cols) >= 5:
                    data["elements"].append({
                        "start": int(cols[0]),
                        "end": int(cols[1]),
                        "direction": cols[2],
                        "source": cols[3],
                        "element": " ".join(cols[4:]),
                    })
    return data


def read_fai(ref_path: Path) -> dict[str, int]:
    """Read chromosome sizes from .fai index."""
    fai = ref_path.parent / (ref_path.name + ".fai")
    sizes: dict[str, int] = {}
    if not fai.exists():
        log(f"FAI index not found: {fai}")
        return sizes
    with open(fai) as fh:
        for line in fh:
            cols = line.rstrip().split("\t")
            if len(cols) >= 2:
                sizes[cols[0]] = int(cols[1])
    return sizes


def read_fasta_lengths(fasta_path: Path) -> int:
    """Read a single-entry FASTA and return sequence length."""
    length = 0
    if not fasta_path.exists():
        return 0
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                length += len(line.rstrip())
    return length


def find_n_gaps(fasta_path: Path, min_gap: int = 10) -> list[tuple[int, int]]:
    """Find N-gap regions in a FASTA file."""
    gaps: list[tuple[int, int]] = []
    if not fasta_path.exists():
        return gaps
    seq_parts: list[str] = []
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq_parts.append(line.rstrip())
    seq = "".join(seq_parts)
    i = 0
    while i < len(seq):
        if seq[i].upper() == "N":
            start = i
            while i < len(seq) and seq[i].upper() == "N":
                i += 1
            if i - start >= min_gap:
                gaps.append((start, i))
        else:
            i += 1
    return gaps


# ---------------------------------------------------------------------------
# Helper: filter overlapping BLAST hits
# ---------------------------------------------------------------------------

def filter_overlapping_hits(hits: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Remove redundant overlapping BLAST hits, keeping longest first."""
    if not hits:
        return []
    sorted_hits = sorted(hits, key=lambda h: -(h["q_end"] - h["q_start"]))
    kept: list[dict[str, Any]] = []
    used: list[tuple[int, int]] = []
    for hit in sorted_hits:
        overlap = False
        hit_len = hit["q_end"] - hit["q_start"]
        if hit_len <= 0:
            continue
        for us, ue in used:
            ov_start = max(hit["q_start"], us)
            ov_end = min(hit["q_end"], ue)
            if ov_end > ov_start:
                frac = (ov_end - ov_start) / hit_len
                if frac > 0.5:
                    overlap = True
                    break
        if not overlap:
            kept.append(hit)
            used.append((hit["q_start"], hit["q_end"]))
    kept.sort(key=lambda h: h["q_start"])
    return kept


# ---------------------------------------------------------------------------
# Panel drawing functions
# ---------------------------------------------------------------------------

def _label_panel(ax: plt.Axes, letter: str) -> None:
    """Add a bold panel label at the top-left corner."""
    ax.text(
        -0.08, 1.12, letter,
        transform=ax.transAxes, fontsize=12, fontweight="bold",
        va="top", ha="left",
    )


def draw_panel_a_ideogram(
    ax: plt.Axes,
    tier_data: list[dict[str, Any]],
    site_data: list[dict[str, Any]],
    chrom_sizes: dict[str, int],
    sample_name: str,
) -> None:
    """Panel A: Chromosome ideogram with insertion sites marked."""
    _label_panel(ax, "A")
    ax.set_title("Insertion site overview", fontsize=10, fontweight="bold",
                 loc="left", pad=8)

    # Build verdict lookup from site_data
    verdict_map: dict[str, str] = {}
    for s in site_data:
        verdict_map[s["site_id"]] = s["verdict"]

    # Determine chromosomes to show
    if chrom_sizes:
        # Filter to main chromosomes (exclude ChrUn, ChrSy, scaffolds)
        chroms = [c for c in chrom_sizes
                  if re.match(r"^(Chr|chr)\d+$", c)]
        # Natural sort
        chroms.sort(key=lambda c: int(re.sub(r"\D", "", c) or "0"))
        sizes = {c: chrom_sizes[c] for c in chroms}
    else:
        # Infer from tier data
        chrom_set: set[str] = set()
        for row in tier_data:
            if re.match(r"^(Chr|chr)\d+$", row["chrom"]):
                chrom_set.add(row["chrom"])
        chroms = sorted(chrom_set, key=lambda c: int(re.sub(r"\D", "", c) or "0"))
        # Estimate sizes
        sizes = {c: 40_000_000 for c in chroms}

    if not chroms:
        ax.text(0.5, 0.5, "No chromosome data available",
                ha="center", va="center", fontsize=9, color="#999999",
                transform=ax.transAxes)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        return

    max_size = max(sizes.values())
    n_chr = len(chroms)
    bar_height = 0.6
    y_spacing = 1.2

    ax.set_xlim(-max_size * 0.12, max_size * 1.05)
    ax.set_ylim(-0.5, n_chr * y_spacing + 0.5)
    ax.invert_yaxis()

    # Build site positions indexed by chrom.
    # Only show transgene-positive sites to avoid overcrowding (there can be
    # thousands of negatives).
    sites_by_chrom: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in tier_data:
        if not row["transgene_positive"]:
            continue  # skip negatives on ideogram
        verdict = verdict_map.get(row["site_id"], "")
        if not verdict or verdict == "unknown":
            verdict = "CANDIDATE"
        sites_by_chrom[row["chrom"]].append({
            "pos": row["pos"],
            "verdict": verdict,
            "site_id": row["site_id"],
        })

    for i, chrom in enumerate(chroms):
        y = i * y_spacing
        size = sizes[chrom]

        # Chromosome bar
        ax.barh(y, size, height=bar_height, left=0,
                color="#E8E8E8", edgecolor="#AAAAAA", linewidth=0.5,
                zorder=1)

        # Chromosome label
        ax.text(-max_size * 0.01, y, chrom.replace("Chr", "Chr "),
                ha="right", va="center", fontsize=7, fontweight="bold")

        # Size label
        ax.text(size + max_size * 0.01, y,
                f"{size / 1e6:.1f} Mb", ha="left", va="center",
                fontsize=6, color="#888888")

        # Mark insertion sites (transgene-positive only)
        for site in sites_by_chrom.get(chrom, []):
            color = VERDICT_COLORS.get(site["verdict"], "#BBBBBB")
            marker_y = y - bar_height / 2 - 0.15
            ax.plot(site["pos"], marker_y, marker="v", color=color,
                    markersize=7, markeredgecolor="black",
                    markeredgewidth=0.4, zorder=3, clip_on=False)

    # Legend
    n_neg = sum(1 for r in tier_data if not r["transgene_positive"])
    legend_elements = [
        plt.Line2D([0], [0], marker="v", color="w", markerfacecolor=VERDICT_COLORS["CANDIDATE"],
                   markeredgecolor="black", markeredgewidth=0.4, markersize=6,
                   label="Candidate"),
        plt.Line2D([0], [0], marker="v", color="w", markerfacecolor=VERDICT_COLORS["FALSE_POSITIVE"],
                   markeredgecolor="black", markeredgewidth=0.4, markersize=6,
                   label="False positive"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=6,
              frameon=True, edgecolor="#CCCCCC", fancybox=False,
              handletextpad=0.3, borderpad=0.3,
              title=f"Transgene-positive only\n({n_neg:,} negatives omitted)",
              title_fontsize=5.5)

    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)


def draw_panel_b_insert_maps(
    ax: plt.Axes,
    site_data: list[dict[str, Any]],
    element_hits: list[dict[str, Any]],
    tier_data: list[dict[str, Any]],
    insert_dir: Path,
) -> None:
    """Panel B: Linear element maps for CANDIDATE inserts."""
    _label_panel(ax, "B")
    ax.set_title("Assembled insert element maps (candidates only)",
                 fontsize=10, fontweight="bold", loc="left", pad=8)

    # Filter to CANDIDATE sites
    candidates = [s for s in site_data if s["verdict"] == "CANDIDATE"]

    if not candidates:
        ax.text(0.5, 0.5, "No candidate insertions detected",
                ha="center", va="center", fontsize=9, color="#999999",
                transform=ax.transAxes)
        ax.axis("off")
        return

    # Build element hits grouped by site_id
    hits_by_site: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for h in element_hits:
        hits_by_site[h["site_id"]].append(h)

    # Build chrom/pos lookup from tier data
    loc_map: dict[str, str] = {}
    for row in tier_data:
        loc_map[row["site_id"]] = f"{row['chrom']}:{row['pos']:,}"

    n_sites = len(candidates)
    track_height = 0.6
    y_spacing = 1.4
    max_len = max(c["insert_length"] for c in candidates) if candidates else 1

    ax.set_xlim(-max_len * 0.03, max_len * 1.08)
    ax.set_ylim(-0.5, n_sites * y_spacing + 0.3)
    ax.invert_yaxis()

    for i, site in enumerate(candidates):
        sid = site["site_id"]
        ins_len = site["insert_length"]
        y = i * y_spacing

        # Backbone bar
        ax.barh(y, ins_len, height=track_height, left=0,
                color="#E8E8E8", edgecolor="#999999", linewidth=0.5,
                zorder=1)

        # N-gaps
        fasta_path = insert_dir / f"{sid}_insert.fasta"
        for gap_start, gap_end in find_n_gaps(fasta_path):
            ax.barh(y, gap_end - gap_start, height=track_height,
                    left=gap_start, color="#FADBD8", edgecolor="#D55E00",
                    linewidth=0.3, hatch="///", zorder=2)

        # Element annotations
        site_hits = filter_overlapping_hits(hits_by_site.get(sid, []))
        for hit in site_hits:
            cat = hit["category"]
            color = ELEMENT_COLORS.get(cat, ELEMENT_COLORS["unknown"])
            edge = "black" if color != "#FFFFFF" else "#AAAAAA"
            width = hit["q_end"] - hit["q_start"]
            rect = FancyBboxPatch(
                (hit["q_start"], y - track_height / 2),
                width, track_height,
                boxstyle="round,pad=0.01",
                facecolor=color, edgecolor=edge, linewidth=0.4,
                alpha=0.85, zorder=3,
            )
            ax.add_patch(rect)

            # Label
            label = hit["short_name"]
            if len(label) > 20:
                label = label[:19] + "..."
            if width > max_len * 0.06:
                ax.text(hit["q_start"] + width / 2, y, label,
                        ha="center", va="center", fontsize=5,
                        fontweight="bold", color="white", zorder=4,
                        clip_on=True)
            elif width > max_len * 0.02:
                ax.text(hit["q_start"] + width / 2,
                        y - track_height / 2 - 0.08, label,
                        ha="center", va="top", fontsize=4,
                        rotation=45, zorder=4, clip_on=True)

        # Site label
        loc = loc_map.get(sid, sid)
        label_text = f"{loc}  ({ins_len:,} bp)"
        ax.text(-max_len * 0.02, y, label_text, ha="right", va="center",
                fontsize=6, fontweight="bold")

    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Position in assembled insert (bp)", fontsize=7)

    # Scale bar
    scale = _nice_scale(max_len)
    bar_y = n_sites * y_spacing - 0.2
    ax.plot([0, scale], [bar_y, bar_y], color="black", linewidth=1.5,
            clip_on=False, zorder=5)
    ax.text(scale / 2, bar_y + 0.2, f"{scale:,} bp", ha="center",
            va="top", fontsize=6, zorder=5)


def draw_panel_c_tier_classification(
    ax: plt.Axes,
    tier_data: list[dict[str, Any]],
    site_data: list[dict[str, Any]],
) -> None:
    """Panel C: Tier classification summary bar chart."""
    _label_panel(ax, "C")
    ax.set_title("Site classification", fontsize=10, fontweight="bold",
                 loc="left", pad=8)

    if not tier_data:
        ax.text(0.5, 0.5, "No tier classification data",
                ha="center", va="center", fontsize=9, color="#999999",
                transform=ax.transAxes)
        ax.axis("off")
        return

    # Count categories
    n_total = len(tier_data)
    n_positive = sum(1 for r in tier_data if r["transgene_positive"])
    n_negative = n_total - n_positive

    # Among positive, count verdicts from site_data
    verdict_map: dict[str, str] = {}
    for s in site_data:
        verdict_map[s["site_id"]] = s["verdict"]

    n_candidate = 0
    n_false_positive = 0
    n_positive_no_assembly = 0
    for row in tier_data:
        if row["transgene_positive"]:
            v = verdict_map.get(row["site_id"], "")
            if v == "CANDIDATE":
                n_candidate += 1
            elif v == "FALSE_POSITIVE":
                n_false_positive += 1
            else:
                n_positive_no_assembly += 1

    # Draw horizontal stacked bar
    categories = ["Transgene-\nnegative", "Transgene-\npositive"]
    y_pos = [0, 1]
    bar_height = 0.55

    # Negative bar
    ax.barh(0, n_negative, height=bar_height, color="#BBBBBB",
            edgecolor="white", linewidth=0.5)
    if n_negative > 0:
        ax.text(n_negative / 2, 0, str(n_negative), ha="center", va="center",
                fontsize=8, fontweight="bold", color="white")

    # Positive: stacked
    left = 0
    for count, color, label in [
        (n_candidate, VERDICT_COLORS["CANDIDATE"], "Candidate"),
        (n_false_positive, VERDICT_COLORS["FALSE_POSITIVE"], "False positive"),
        (n_positive_no_assembly, "#F0E442", "Positive (no asm)"),
    ]:
        if count > 0:
            ax.barh(1, count, height=bar_height, left=left,
                    color=color, edgecolor="white", linewidth=0.5)
            if count >= 2 or n_positive <= 5:
                ax.text(left + count / 2, 1, str(count), ha="center",
                        va="center", fontsize=8, fontweight="bold",
                        color="white" if color != "#F0E442" else "black")
            left += count

    ax.set_yticks(y_pos)
    ax.set_yticklabels(categories, fontsize=7)
    ax.set_xlabel("Number of sites", fontsize=7)
    ax.set_xlim(0, max(n_total, n_positive, 1) * 1.1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Compact legend
    legend_patches = [
        mpatches.Patch(color=VERDICT_COLORS["CANDIDATE"], label="Candidate"),
        mpatches.Patch(color=VERDICT_COLORS["FALSE_POSITIVE"],
                       label="False positive"),
        mpatches.Patch(color="#BBBBBB", label="Transgene-negative"),
    ]
    if n_positive_no_assembly > 0:
        legend_patches.append(
            mpatches.Patch(color="#F0E442", label="Positive (no assembly)")
        )
    ax.legend(handles=legend_patches, loc="lower right", fontsize=6,
              frameon=True, edgecolor="#CCCCCC", fancybox=False,
              ncol=2, handletextpad=0.3, borderpad=0.3)


def draw_panel_d_convergence(
    ax: plt.Axes,
    site_data: list[dict[str, Any]],
    tier_data: list[dict[str, Any]],
) -> None:
    """Panel D: Assembly convergence summary."""
    _label_panel(ax, "D")
    ax.set_title("Assembly convergence", fontsize=10, fontweight="bold",
                 loc="left", pad=8)

    if not site_data:
        ax.text(0.5, 0.5, "No assembly data available",
                ha="center", va="center", fontsize=9, color="#999999",
                transform=ax.transAxes)
        ax.axis("off")
        return

    # Build chrom:pos lookup
    loc_map: dict[str, str] = {}
    for row in tier_data:
        loc_map[row["site_id"]] = f"{row['chrom']}:{row['pos']:,}"

    # Sort by insert length descending
    sorted_sites = sorted(site_data, key=lambda s: -s["insert_length"])

    labels: list[str] = []
    lengths: list[int] = []
    rounds: list[int] = []
    colors: list[str] = []

    for site in sorted_sites:
        sid = site["site_id"]
        loc = loc_map.get(sid, sid.replace("insertion_", ""))
        labels.append(loc)
        lengths.append(site["insert_length"])
        rounds.append(site["assembly_rounds"])
        colors.append(VERDICT_COLORS.get(site["verdict"], "#BBBBBB"))

    y_pos = np.arange(len(labels))
    bar_height = 0.55

    bars = ax.barh(y_pos, lengths, height=bar_height, color=colors,
                   edgecolor="white", linewidth=0.5)

    # Annotate with rounds
    for i, (ln, rd) in enumerate(zip(lengths, rounds)):
        ax.text(ln + max(lengths) * 0.01, i,
                f"{ln:,} bp  (r{rd})",
                ha="left", va="center", fontsize=6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlabel("Insert length (bp)", fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, max(lengths) * 1.35 if lengths else 1)
    ax.invert_yaxis()


def draw_panel_e_composition(
    ax: plt.Axes,
    element_hits: list[dict[str, Any]],
    site_data: list[dict[str, Any]],
) -> None:
    """Panel E: Element composition donut chart for CANDIDATE inserts."""
    _label_panel(ax, "E")
    ax.set_title("Element composition (candidates)",
                 fontsize=10, fontweight="bold", loc="left", pad=8)

    # Filter hits to CANDIDATE sites only
    candidate_ids = {s["site_id"] for s in site_data if s["verdict"] == "CANDIDATE"}
    candidate_hits = [h for h in element_hits if h["site_id"] in candidate_ids]

    if not candidate_hits:
        ax.text(0.5, 0.5, "No element annotations",
                ha="center", va="center", fontsize=9, color="#999999",
                transform=ax.transAxes)
        ax.axis("off")
        return

    # Sum aligned bases per category
    cat_bases: Counter[str] = Counter()
    for h in candidate_hits:
        cat = h["category"]
        cat_bases[cat] += h["length"]

    # Compute unmatched bases
    total_insert_bp = sum(s["insert_length"] for s in site_data
                         if s["site_id"] in candidate_ids)
    total_annotated = sum(cat_bases.values())
    unmatched = max(0, total_insert_bp - total_annotated)
    if unmatched > 0:
        cat_bases["unannotated"] = unmatched

    # Sort for consistent display
    order = ["promoter", "cds", "terminator", "marker", "border", "vector",
             "element-specific", "event-specific", "construct-specific",
             "host-endogenous", "unknown", "unannotated"]
    labels: list[str] = []
    sizes: list[int] = []
    colors: list[str] = []

    unannotated_color = "#F5F5F5"

    for cat in order:
        if cat in cat_bases and cat_bases[cat] > 0:
            labels.append(cat.replace("-", " ").capitalize())
            sizes.append(cat_bases[cat])
            if cat == "unannotated":
                colors.append(unannotated_color)
            else:
                colors.append(ELEMENT_COLORS.get(cat, "#DDDDDD"))

    # Remaining categories not in order
    for cat, val in cat_bases.items():
        if cat not in order and val > 0:
            labels.append(cat.replace("-", " ").capitalize())
            sizes.append(val)
            colors.append(ELEMENT_COLORS.get(cat, "#DDDDDD"))

    wedges, texts, autotexts = ax.pie(
        sizes, labels=None, colors=colors,
        autopct=lambda pct: f"{pct:.0f}%" if pct >= 5 else "",
        startangle=90, pctdistance=0.78,
        wedgeprops={"width": 0.45, "edgecolor": "white", "linewidth": 0.8},
        textprops={"fontsize": 6},
    )

    for at in autotexts:
        at.set_fontsize(6)
        at.set_fontweight("bold")

    # Center text
    total_kb = total_insert_bp / 1000
    ax.text(0, 0, f"{total_kb:.1f}\nkb",
            ha="center", va="center", fontsize=9, fontweight="bold")

    # Legend below
    legend_patches = [mpatches.Patch(color=c, label=f"{l} ({s:,} bp)")
                      for l, s, c in zip(labels, sizes, colors)]
    ax.legend(handles=legend_patches, loc="center left",
              bbox_to_anchor=(1.0, 0.5), fontsize=5.5,
              frameon=False, handletextpad=0.3, borderpad=0.2)

    ax.set_aspect("equal")


def draw_panel_f_insert_waterfall(
    ax: plt.Axes,
    site_data: list[dict[str, Any]],
    element_hits: list[dict[str, Any]],
    tier_data: list[dict[str, Any]],
) -> None:
    """Panel F: Insert coverage by annotated elements (waterfall/stacked)."""
    _label_panel(ax, "F")
    ax.set_title("Insert annotation coverage per site",
                 fontsize=10, fontweight="bold", loc="left", pad=8)

    candidates = [s for s in site_data if s["verdict"] == "CANDIDATE"]
    if not candidates:
        ax.text(0.5, 0.5, "No candidate insertions",
                ha="center", va="center", fontsize=9, color="#999999",
                transform=ax.transAxes)
        ax.axis("off")
        return

    # Build per-site annotation coverage
    hits_by_site: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for h in element_hits:
        hits_by_site[h["site_id"]].append(h)

    loc_map: dict[str, str] = {}
    for row in tier_data:
        loc_map[row["site_id"]] = f"{row['chrom']}:{row['pos']:,}"

    sorted_cands = sorted(candidates, key=lambda s: -s["insert_length"])
    labels: list[str] = []
    bar_data: list[dict[str, int]] = []

    all_cats_used: set[str] = set()
    for site in sorted_cands:
        sid = site["site_id"]
        loc = loc_map.get(sid, sid.replace("insertion_", ""))
        labels.append(loc)

        site_hits = filter_overlapping_hits(hits_by_site.get(sid, []))
        cat_bp: dict[str, int] = defaultdict(int)
        total_ann = 0
        for h in site_hits:
            bp = h["q_end"] - h["q_start"]
            cat_bp[h["category"]] += bp
            total_ann += bp
        unannotated = max(0, site["insert_length"] - total_ann)
        cat_bp["unannotated"] = unannotated
        bar_data.append(dict(cat_bp))
        all_cats_used.update(cat_bp.keys())

    # Stacked horizontal bars
    cat_order = ["promoter", "cds", "terminator", "marker", "border", "vector",
                 "element-specific", "event-specific", "construct-specific",
                 "host-endogenous", "unknown", "unannotated"]
    unannotated_color = "#F5F5F5"
    y_pos = np.arange(len(labels))
    bar_height = 0.55

    for cat in cat_order:
        if cat not in all_cats_used:
            continue
        lefts = np.zeros(len(labels))
        widths = np.array([bd.get(cat, 0) for bd in bar_data], dtype=float)

        # Compute left positions
        for prev_cat in cat_order:
            if prev_cat == cat:
                break
            if prev_cat in all_cats_used:
                lefts += np.array([bd.get(prev_cat, 0) for bd in bar_data],
                                 dtype=float)

        color = unannotated_color if cat == "unannotated" else \
            ELEMENT_COLORS.get(cat, "#DDDDDD")
        edge = "#AAAAAA" if color in ("#FFFFFF", "#F5F5F5") else "white"

        ax.barh(y_pos, widths, height=bar_height, left=lefts,
                color=color, edgecolor=edge, linewidth=0.3,
                label=cat.replace("-", " ").capitalize())

    # Total length labels
    for i, site in enumerate(sorted_cands):
        ax.text(site["insert_length"] + max(s["insert_length"]
                for s in sorted_cands) * 0.01,
                i, f"{site['insert_length']:,} bp",
                ha="left", va="center", fontsize=6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlabel("Position in insert (bp)", fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.invert_yaxis()
    max_len = max(s["insert_length"] for s in sorted_cands)
    ax.set_xlim(0, max_len * 1.18)

    # Legend
    used_cats = [c for c in cat_order if c in all_cats_used]
    legend_patches = []
    for cat in used_cats:
        color = unannotated_color if cat == "unannotated" else \
            ELEMENT_COLORS.get(cat, "#DDDDDD")
        legend_patches.append(
            mpatches.Patch(color=color,
                           label=cat.replace("-", " ").capitalize()))
    ax.legend(handles=legend_patches, loc="lower right", fontsize=5.5,
              frameon=True, edgecolor="#CCCCCC", fancybox=False, ncol=2,
              handletextpad=0.3, borderpad=0.3)


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _nice_scale(max_val: int) -> int:
    """Pick a round scale bar value."""
    if max_val > 10000:
        return 5000
    if max_val > 5000:
        return 2000
    if max_val > 2000:
        return 1000
    if max_val > 1000:
        return 500
    return 100


# ---------------------------------------------------------------------------
# Main figure assembly
# ---------------------------------------------------------------------------

def make_summary_figure(
    sample_name: str,
    insert_dir: Path,
    outdir: Path,
    host_ref: Path | None = None,
    gff: Path | None = None,
) -> None:
    """Build and save the multi-panel summary figure."""
    apply_style()

    # ---- Load data ----
    stats = parse_stats(insert_dir / "s05_stats.txt")
    site_data = extract_site_data(stats)
    tier_data = parse_tier_classification(insert_dir / "site_tier_classification.tsv")
    element_hits = parse_element_annotation(insert_dir / "element_annotation.tsv")

    log(f"Sample: {sample_name}")
    log(f"Sites in s05_stats: {len(site_data)}")
    log(f"Sites in tier classification: {len(tier_data)}")
    log(f"Element annotation hits: {len(element_hits)}")

    n_candidates = sum(1 for s in site_data if s["verdict"] == "CANDIDATE")
    n_fp = sum(1 for s in site_data if s["verdict"] == "FALSE_POSITIVE")
    log(f"Candidates: {n_candidates}, False positives: {n_fp}")

    # Chromosome sizes
    chrom_sizes: dict[str, int] = {}
    if host_ref and host_ref.exists():
        chrom_sizes = read_fai(host_ref)
        log(f"Loaded {len(chrom_sizes)} chromosome sizes from FAI")

    # ---- Layout ----
    # 3 rows x 2 columns:
    #   Row 0: Panel A (ideogram, spans full width)
    #   Row 1: Panel B (insert maps, spans full width)
    #   Row 2: Panel C (left) + Panel D (right)
    #   Row 3: Panel E (left) + Panel F (right)

    fig = plt.figure(figsize=(14, 18))
    fig.suptitle(f"RedGene Sample Summary: {sample_name}",
                 fontsize=14, fontweight="bold", y=0.98)

    # Compute dynamic height ratios
    n_chr = len([c for c in chrom_sizes if re.match(r"^(Chr|chr)\d+$", c)]) if chrom_sizes else 12
    n_cand_sites = max(n_candidates, 1)
    n_asm_sites = max(len(site_data), 1)

    # Normalized heights
    h_a = max(3.0, n_chr * 0.28)
    h_b = max(2.0, n_cand_sites * 0.7)
    h_cd = max(2.0, n_asm_sites * 0.35)
    h_ef = 3.0

    gs = gridspec.GridSpec(
        4, 2,
        figure=fig,
        height_ratios=[h_a, h_b, h_cd, h_ef],
        width_ratios=[1, 1],
        hspace=0.38,
        wspace=0.40,
        left=0.10,
        right=0.92,
        top=0.95,
        bottom=0.03,
    )

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :])
    ax_c = fig.add_subplot(gs[2, 0])
    ax_d = fig.add_subplot(gs[2, 1])
    ax_e = fig.add_subplot(gs[3, 0])
    ax_f = fig.add_subplot(gs[3, 1])

    # ---- Draw panels ----
    draw_panel_a_ideogram(ax_a, tier_data, site_data, chrom_sizes, sample_name)
    draw_panel_b_insert_maps(ax_b, site_data, element_hits, tier_data, insert_dir)
    draw_panel_c_tier_classification(ax_c, tier_data, site_data)
    draw_panel_d_convergence(ax_d, site_data, tier_data)
    draw_panel_e_composition(ax_e, element_hits, site_data)
    draw_panel_f_insert_waterfall(ax_f, site_data, element_hits, tier_data)

    # ---- Save ----
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ["png", "pdf"]:
        out_path = outdir / f"{sample_name}_summary.{ext}"
        fig.savefig(
            out_path,
            dpi=300 if ext == "png" else None,
            bbox_inches="tight",
        )
        log(f"Saved: {out_path}")

    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a multi-panel sample summary figure for RedGene.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python scripts/viz/plot_sample_summary.py \\
    --sample-name rice_G281 \\
    --insert-dir results/rice_G281/s05_insert_assembly \\
    --outdir results \\
    --host-ref db/Osativa_323_v7.0.fa
        """,
    )
    parser.add_argument("--sample-name", required=True,
                        help="Sample name (used for titles and output filenames)")
    parser.add_argument("--insert-dir", required=True,
                        help="Path to s05_insert_assembly output directory")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for summary figure")
    parser.add_argument("--host-ref", default=None,
                        help="Host reference FASTA (reads .fai for chrom sizes)")
    parser.add_argument("--gff", default=None,
                        help="GFF3 annotation (reserved for future gene context)")
    args = parser.parse_args()

    insert_dir = Path(args.insert_dir)
    if not insert_dir.exists():
        log(f"ERROR: s09 directory not found: {insert_dir}")
        sys.exit(1)

    host_ref = Path(args.host_ref) if args.host_ref else None
    gff = Path(args.gff) if args.gff else None

    make_summary_figure(
        sample_name=args.sample_name,
        insert_dir=insert_dir,
        outdir=Path(args.outdir),
        host_ref=host_ref,
        gff=gff,
    )
    log("Done.")


if __name__ == "__main__":
    main()
