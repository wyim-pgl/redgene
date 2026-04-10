#!/usr/bin/env python3
"""Visualize assembled insert structure from Step 5 annotation.

Draws a linear map of the reconstructed transgene insert showing:
  - Element annotations (promoters, CDS, terminators, borders) from BLAST
  - N-gap regions (unfilled portions)
  - T-DNA border motif positions
  - Multi-construct cassette boundaries (detected automatically)
  - CRL amplicon matches for event fingerprinting

Layout:
  Row 1: Insert backbone (gray bar with N-gaps marked)
  Row 2: Element annotations (colored blocks per category)
  Row 3: T-DNA border positions
  Row 4: Construct cassette boundaries (if multi-construct)
  Row 5: Legend

Usage:
  python plot_insert_structure.py \
    --insert results/{sample}/s05_insert_assembly/insert_only.fasta \
    --annotation results/{sample}/s05_insert_assembly/element_annotation.tsv \
    --borders results/{sample}/s05_insert_assembly/border_hits.tsv \
    --stats results/{sample}/s05_insert_assembly/s05_stats.txt \
    --sample-name {sample} --outdir results
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np


def log(msg: str) -> None:
    print(f"[plot_insert_structure] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Color scheme for element categories
# ---------------------------------------------------------------------------

CATEGORY_COLORS = {
    "promoter": "#E74C3C",       # red
    "terminator": "#3498DB",     # blue
    "coding": "#27AE60",         # green
    "border": "#E67E22",         # dark orange (T-DNA LB/RB)
    "intron": "#F1C40F",         # yellow (introns, enhancers)
    "insulator": "#8E44AD",      # dark purple (insulators, MAR)
    "transit": "#16A085",        # dark teal (transit peptides, CTP)
    "regulatory": "#F39C12",     # orange (other regulatory)
    "construct": "#9B59B6",      # purple (construct-specific)
    "element": "#1ABC9C",        # teal (element-specific amplicons)
    "event": "#E67E22",          # dark orange (event-specific)
    "taxon": "#95A5A6",          # gray (taxon-specific)
    "unknown": "#BDC3C7",        # light gray
}

CATEGORY_ORDER = [
    "promoter", "coding", "terminator", "border", "intron",
    "insulator", "transit", "regulatory",
    "construct", "element", "event", "taxon", "unknown",
]


def classify_element(element_name: str) -> str:
    """Classify an element BLAST hit into a display category."""
    name_lower = element_name.lower()

    # Promoters
    if any(k in name_lower for k in ["promoter", "p-35s", "p-e35s", "p-fmv",
                                      "p-nos", "p-ocs", "p-ubi", "p-act1",
                                      "p-ta29", "p-csvmv", "p-ssuara", "p-cdpk",
                                      "p-rice_actin", "p-fmv34s"]):
        return "promoter"
    # Terminators
    if any(k in name_lower for k in ["terminator", "t-nos", "t-35s", "t-e9",
                                      "t-ocs", "t-pcambia", "t-pinii",
                                      "nos_ter", "polya"]):
        return "terminator"
    # T-DNA borders
    if any(k in name_lower for k in ["lb-tdna", "rb-tdna", "left_border",
                                      "right_border"]):
        return "border"
    if name_lower.endswith("-lb") or name_lower.endswith("-rb"):
        return "border"
    # Introns / enhancers
    if any(k in name_lower for k in ["intron", "i-hsp", "ivs", "enhancer",
                                      "aint", "adh1"]):
        return "intron"
    # Insulators / MAR / scaffold attachment
    if any(k in name_lower for k in ["insulator", "mar", "sar",
                                      "scaffold_attachment", "rb7",
                                      "matrix_attachment"]):
        return "insulator"
    # Transit peptides / signal peptides
    if any(k in name_lower for k in ["ctp", "transit", "signal_peptide",
                                      "chloroplast_transit"]):
        return "transit"
    # Coding sequences
    if any(k in name_lower for k in ["nptii", "bar", "pat", "hpt", "aada",
                                      "cp4-epsps", "epsps", "cry1", "cry2",
                                      "cry3", "gus", "barnase", "barstar",
                                      "thaumatin", "gat", "dsred", "mepsps",
                                      "cat-f", "rerio"]):
        return "coding"
    # Other regulatory
    if any(k in name_lower for k in ["border", "ori", "backbone"]):
        return "regulatory"
    # qPCR amplicon categories
    if "ql-con" in name_lower or "qt-con" in name_lower or "construct-specific" in name_lower:
        return "construct"
    if "ql-ele" in name_lower or "qt-ele" in name_lower or "element-specific" in name_lower:
        return "element"
    if "ql-eve" in name_lower or "qt-eve" in name_lower or "event-specific" in name_lower:
        return "event"
    if "ql-tax" in name_lower or "qt-tax" in name_lower or "taxon-specific" in name_lower:
        return "taxon"
    return "unknown"


def short_name(element_name: str) -> str:
    """Extract a short display name from a full element identifier."""
    parts = element_name.split("|")
    if len(parts) >= 3:
        # EUginius format: type|ID|target|accession|desc|source
        return parts[1]  # ID
    if len(parts) == 1:
        # Simple name
        return element_name.split("_")[0] if "_" in element_name else element_name
    return parts[0]


# ---------------------------------------------------------------------------
# Data parsing
# ---------------------------------------------------------------------------

def read_fasta_single(path: Path) -> tuple[str, str]:
    """Read FASTA with a single sequence."""
    name = ""
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                parts.append(line)
    return name, "".join(parts)


def parse_annotation(path: Path) -> list[dict]:
    """Parse element_annotation.tsv from s09."""
    hits = []
    if not path.exists() or path.stat().st_size == 0:
        return hits
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            hits.append({
                "element": row["element"],
                "identity": float(row["identity"]),
                "length": int(row["length"]),
                "q_start": int(row["q_start"]),
                "q_end": int(row["q_end"]),
                "category": classify_element(row["element"]),
                "short_name": short_name(row["element"]),
            })
    # Sort by position
    hits.sort(key=lambda h: h["q_start"])
    return hits


def parse_borders(path: Path) -> list[dict]:
    """Parse border_hits.tsv from s09."""
    borders = []
    if not path.exists() or path.stat().st_size == 0:
        return borders
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip().split("\t")
            if len(cols) < 8:
                continue
            borders.append({
                "type": cols[0],  # RB_consensus or LB_consensus
                "identity": float(cols[2]),
                "length": int(cols[3]),
                "s_start": int(cols[6]),
                "s_end": int(cols[7]),
            })
    return borders


def parse_stats(path: Path) -> dict:
    """Parse s05_stats.txt."""
    stats = {}
    if not path.exists():
        return stats
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip().split("\t", 1)
            if len(parts) == 2:
                stats[parts[0]] = parts[1]
    return stats


def find_n_gaps(seq: str, min_gap: int = 10) -> list[tuple[int, int]]:
    """Find N-gap regions in the sequence."""
    gaps = []
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
# Multi-construct detection
# ---------------------------------------------------------------------------

def detect_cassettes(hits: list[dict], insert_len: int) -> list[dict]:
    """Detect construct cassette boundaries from element order.

    A cassette is promoter → CDS → terminator. If we see this pattern
    repeated, it indicates multi-construct insertion (tandem T-DNA).
    """
    # Find promoter-terminator pairs
    promoters = [h for h in hits if h["category"] == "promoter"]
    terminators = [h for h in hits if h["category"] == "terminator"]
    cds_hits = [h for h in hits if h["category"] == "coding"]

    if not promoters or not terminators:
        return []

    cassettes = []
    used_terms = set()

    for p in promoters:
        # Find the nearest downstream terminator
        best_term = None
        best_dist = float("inf")
        for t_idx, t in enumerate(terminators):
            if t_idx in used_terms:
                continue
            if t["q_start"] > p["q_end"]:
                dist = t["q_start"] - p["q_end"]
                if dist < best_dist:
                    best_dist = dist
                    best_term = (t_idx, t)

        if best_term is not None:
            t_idx, t = best_term
            used_terms.add(t_idx)

            # Find CDS between this promoter and terminator
            cassette_cds = [c for c in cds_hits
                            if c["q_start"] >= p["q_start"]
                            and c["q_end"] <= t["q_end"]]

            cassettes.append({
                "start": p["q_start"],
                "end": t["q_end"],
                "promoter": p["short_name"],
                "terminator": t["short_name"],
                "cds": [c["short_name"] for c in cassette_cds],
                "label": f"{p['short_name']}→{'→'.join(c['short_name'] for c in cassette_cds) if cassette_cds else '?'}→{t['short_name']}",
            })

    return cassettes


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_insert(
    insert_name: str,
    insert_seq: str,
    hits: list[dict],
    borders: list[dict],
    stats: dict,
    cassettes: list[dict],
    sample_name: str,
    outdir: Path,
) -> None:
    """Generate the insert structure visualization."""
    insert_len = len(insert_seq)
    n_gaps = find_n_gaps(insert_seq)

    # Determine number of tracks
    n_cassettes = len(cassettes)
    has_borders = len(borders) > 0
    n_rows = 3 + (1 if has_borders else 0) + (1 if n_cassettes > 1 else 0)

    fig_width = max(12, insert_len / 500)
    fig_width = min(fig_width, 24)
    fig_height = 2 + n_rows * 1.2

    fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, fig_height),
                              gridspec_kw={"height_ratios":
                                           [1.5] + [2.0] +
                                           ([0.8] if has_borders else []) +
                                           ([1.2] if n_cassettes > 1 else []) +
                                           [0.8]})
    if n_rows == 1:
        axes = [axes]

    for ax in axes:
        ax.set_xlim(-insert_len * 0.02, insert_len * 1.02)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax_idx = 0

    # ---- Track 1: Insert backbone ----
    ax = axes[ax_idx]
    ax_idx += 1

    # Draw backbone
    backbone_y = 0.5
    ax.barh(backbone_y, insert_len, height=0.3, left=0,
            color="#D5D8DC", edgecolor="#7F8C8D", linewidth=0.5)

    # Mark N-gaps
    for gap_start, gap_end in n_gaps:
        ax.barh(backbone_y, gap_end - gap_start, height=0.3, left=gap_start,
                color="#FADBD8", edgecolor="#E74C3C", linewidth=0.5,
                hatch="///")

    # Scale bar
    scale = 1000 if insert_len > 5000 else 500 if insert_len > 2000 else 100
    ax.plot([0, scale], [-0.1, -0.1], "k-", linewidth=2)
    ax.text(scale / 2, -0.35, f"{scale} bp", ha="center", fontsize=8)

    # Insert length label
    n_count = sum(1 for c in insert_seq if c.upper() == "N")
    status = stats.get("insertion_1_status", "unknown")
    ax.set_title(f"{sample_name} — Assembled Insert ({insert_len:,} bp, "
                 f"{n_count:,} Ns, {status})",
                 fontsize=12, fontweight="bold", pad=10)
    ax.set_ylim(-0.6, 1.0)
    ax.set_yticks([])
    ax.set_ylabel("Insert", fontsize=9)

    # ---- Track 2: Element annotations ----
    ax = axes[ax_idx]
    ax_idx += 1

    # Deduplicate overlapping hits: keep best identity per region
    filtered = _filter_overlapping_hits(hits)

    # Assign y-levels to avoid overlaps
    levels = _assign_levels(filtered)
    max_level = max(levels) if levels else 0

    for hit, level in zip(filtered, levels):
        cat = hit["category"]
        color = CATEGORY_COLORS.get(cat, CATEGORY_COLORS["unknown"])
        y = level * 0.7

        # Draw element block
        width = hit["q_end"] - hit["q_start"]
        rect = FancyBboxPatch(
            (hit["q_start"], y - 0.25), width, 0.5,
            boxstyle="round,pad=0.02",
            facecolor=color, edgecolor="black", linewidth=0.5,
            alpha=0.85,
        )
        ax.add_patch(rect)

        # Label (rotate if narrow)
        label = hit["short_name"]
        if len(label) > 15:
            label = label[:14] + "…"
        center_x = (hit["q_start"] + hit["q_end"]) / 2

        if width > insert_len * 0.05:
            ax.text(center_x, y, label, ha="center", va="center",
                    fontsize=6, fontweight="bold", color="white",
                    clip_on=True)
        else:
            ax.text(center_x, y + 0.4, label, ha="center", va="bottom",
                    fontsize=5, rotation=45, clip_on=True)

    ax.set_ylim(-0.5, (max_level + 1) * 0.7 + 0.5)
    ax.set_yticks([])
    ax.set_ylabel("Elements", fontsize=9)

    # ---- Track 3 (optional): T-DNA borders ----
    if has_borders:
        ax = axes[ax_idx]
        ax_idx += 1

        for b in borders:
            pos = (b["s_start"] + b["s_end"]) / 2
            color = "#E74C3C" if "RB" in b["type"] else "#2ECC71"
            label = "RB" if "RB" in b["type"] else "LB"
            ax.axvline(pos, color=color, linewidth=2, alpha=0.8)
            ax.text(pos, 0.7, label, ha="center", va="bottom",
                    fontsize=8, fontweight="bold", color=color)

        ax.set_ylim(0, 1.0)
        ax.set_yticks([])
        ax.set_ylabel("Borders", fontsize=9)

    # ---- Track 4 (optional): Cassette boundaries ----
    if n_cassettes > 1:
        ax = axes[ax_idx]
        ax_idx += 1

        cassette_colors = ["#AED6F1", "#ABEBC6", "#F9E79F", "#F5CBA7",
                           "#D7BDE2", "#A3E4D7"]

        for i, cas in enumerate(cassettes):
            color = cassette_colors[i % len(cassette_colors)]
            width = cas["end"] - cas["start"]
            rect = FancyBboxPatch(
                (cas["start"], 0.1), width, 0.6,
                boxstyle="round,pad=0.02",
                facecolor=color, edgecolor="black", linewidth=1,
                alpha=0.6,
            )
            ax.add_patch(rect)
            center = (cas["start"] + cas["end"]) / 2
            ax.text(center, 0.4, f"Cassette {i + 1}\n{cas['label']}",
                    ha="center", va="center", fontsize=6, fontweight="bold")

        ax.set_ylim(-0.1, 1.0)
        ax.set_yticks([])
        ax.set_ylabel(f"Cassettes\n({n_cassettes})", fontsize=9)
        ax.set_title(f"Multi-construct: {n_cassettes} cassettes detected",
                     fontsize=9, fontstyle="italic", loc="left")

    # ---- Legend track ----
    ax = axes[ax_idx]
    ax.set_visible(True)
    ax.axis("off")

    legend_patches = []
    used_cats = set(h["category"] for h in filtered)
    for cat in CATEGORY_ORDER:
        if cat in used_cats:
            legend_patches.append(
                mpatches.Patch(color=CATEGORY_COLORS[cat], label=cat.capitalize())
            )
    # Add N-gap legend
    legend_patches.append(
        mpatches.Patch(facecolor="#FADBD8", edgecolor="#E74C3C",
                       hatch="///", label="N-gap (unfilled)")
    )
    if has_borders:
        legend_patches.append(
            mpatches.Patch(color="#E74C3C", label="RB border")
        )
        legend_patches.append(
            mpatches.Patch(color="#2ECC71", label="LB border")
        )

    ax.legend(handles=legend_patches, loc="center", ncol=min(6, len(legend_patches)),
              fontsize=8, frameon=False)

    plt.tight_layout()

    # Save
    out_base = outdir / sample_name / "s05_insert_assembly"
    out_base.mkdir(parents=True, exist_ok=True)
    for ext in ["png", "pdf"]:
        out_path = out_base / f"insert_structure.{ext}"
        fig.savefig(out_path, dpi=200 if ext == "png" else None,
                    bbox_inches="tight")
        log(f"Saved: {out_path}")
    plt.close(fig)


def _filter_overlapping_hits(hits: list[dict]) -> list[dict]:
    """Remove redundant overlapping BLAST hits, keeping best identity."""
    if not hits:
        return []

    # Sort by length descending (prefer longer hits)
    sorted_hits = sorted(hits, key=lambda h: -(h["q_end"] - h["q_start"]))
    kept = []
    used_ranges: list[tuple[int, int]] = []

    for hit in sorted_hits:
        # Check overlap with already-kept hits
        overlap = False
        for ur_start, ur_end in used_ranges:
            overlap_start = max(hit["q_start"], ur_start)
            overlap_end = min(hit["q_end"], ur_end)
            if overlap_end > overlap_start:
                overlap_frac = (overlap_end - overlap_start) / (hit["q_end"] - hit["q_start"])
                if overlap_frac > 0.5:
                    overlap = True
                    break
        if not overlap:
            kept.append(hit)
            used_ranges.append((hit["q_start"], hit["q_end"]))

    kept.sort(key=lambda h: h["q_start"])
    return kept


def _assign_levels(hits: list[dict]) -> list[int]:
    """Assign y-levels to avoid visual overlap."""
    if not hits:
        return []

    levels = []
    level_ends: list[int] = []  # track rightmost x for each level

    for hit in hits:
        placed = False
        for lev_idx, lev_end in enumerate(level_ends):
            if hit["q_start"] > lev_end + 50:  # 50bp gap
                levels.append(lev_idx)
                level_ends[lev_idx] = hit["q_end"]
                placed = True
                break
        if not placed:
            levels.append(len(level_ends))
            level_ends.append(hit["q_end"])

    return levels


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualize assembled insert structure from Step 5")
    parser.add_argument("--insert", required=True,
                        help="Insert FASTA from s09 (insert_only.fasta)")
    parser.add_argument("--annotation", required=True,
                        help="element_annotation.tsv from s09")
    parser.add_argument("--borders", default=None,
                        help="border_hits.tsv from s09")
    parser.add_argument("--stats", default=None,
                        help="s05_stats.txt")
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    insert_path = Path(args.insert)
    if not insert_path.exists():
        log(f"Insert FASTA not found: {insert_path}")
        sys.exit(1)

    insert_name, insert_seq = read_fasta_single(insert_path)
    log(f"Insert: {insert_name}, {len(insert_seq):,} bp")

    # Parse annotation
    hits = parse_annotation(Path(args.annotation))
    log(f"Element hits: {len(hits)}")

    # Parse borders
    borders = []
    if args.borders:
        borders = parse_borders(Path(args.borders))
        log(f"Border hits: {len(borders)}")

    # Parse stats
    stats = {}
    if args.stats:
        stats = parse_stats(Path(args.stats))

    # Detect multi-construct cassettes
    cassettes = detect_cassettes(hits, len(insert_seq))
    if len(cassettes) > 1:
        log(f"Multi-construct detected: {len(cassettes)} cassettes")
        for i, cas in enumerate(cassettes):
            log(f"  Cassette {i + 1}: {cas['label']} "
                f"({cas['start']}-{cas['end']}, {cas['end'] - cas['start']}bp)")
    elif len(cassettes) == 1:
        log(f"Single cassette: {cassettes[0]['label']}")

    # Generate plot
    plot_insert(
        insert_name, insert_seq,
        hits, borders, stats, cassettes,
        args.sample_name, Path(args.outdir),
    )

    log("Done.")


if __name__ == "__main__":
    main()
