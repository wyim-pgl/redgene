#!/usr/bin/env python3
"""Plot 2: Junction structure diagram.

Linear diagram of chimeric contigs showing host/transgene annotation
at detected junctions.
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Plot junction structure diagram")
    parser.add_argument("--junctions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample-name", type=str, default="sample")
    args = parser.parse_args()

    # Read junctions
    junctions = []
    with open(args.junctions) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            junctions.append(row)

    if not junctions:
        print("No junctions to plot", file=sys.stderr)
        sys.exit(0)

    # Group by contig
    from collections import defaultdict
    contig_junctions = defaultdict(list)
    for j in junctions:
        contig_junctions[j["contig_name"]].append(j)

    n_contigs = len(contig_junctions)
    fig, ax = plt.subplots(1, 1, figsize=(14, max(2 + n_contigs * 1.5, 4)))
    ax.set_title(f"Junction Structure: {args.sample_name}",
                 fontsize=14, fontweight="bold")

    y_pos = 0
    yticks = []
    yticklabels = []
    host_color = "#4CAF50"
    construct_color = "#F44336"

    for contig_name, juncs in sorted(contig_junctions.items()):
        contig_len = int(juncs[0]["contig_len"])

        # Draw contig backbone
        ax.barh(y_pos, contig_len, height=0.6, color="#E0E0E0",
                edgecolor="black", linewidth=0.5)

        # Draw host regions
        for j in juncs:
            host_start = int(j["host_start"])
            host_end = int(j["host_end"])
            junction_type = j["junction_type"]
            host_chr = j["host_chr"]
            junction_pos = j["junction_pos_host"]

            # Approximate contig position for host alignment
            # (host part takes up part of the contig)
            host_region_len = abs(host_end - host_start)

            # Label
            label = f"{host_chr}:{junction_pos} [{junction_type}]"
            ax.annotate(label,
                       xy=(contig_len + 10, y_pos),
                       fontsize=8, va="center",
                       color=host_color if junction_type != "RB" else construct_color)

        # Contig name label
        short_name = contig_name.split("_cov_")[0] if "_cov_" in contig_name else contig_name
        yticks.append(y_pos)
        yticklabels.append(short_name)

        y_pos -= 1.5

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=8)
    ax.set_xlabel("Contig position (bp)")
    ax.set_xlim(-20, max(int(j["contig_len"]) for j in junctions) + 200)

    # Legend
    legend_patches = [
        mpatches.Patch(color=host_color, label="Host genome"),
        mpatches.Patch(color=construct_color, label="Construct/T-DNA"),
        mpatches.Patch(color="#E0E0E0", label="Contig backbone"),
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=9)

    # Summary table below
    summary_text = "Junction Summary:\n"
    for j in junctions:
        summary_text += (f"  {j['contig_name'][:30]}... "
                        f"{j['host_chr']}:{j['junction_pos_host']} "
                        f"[{j['junction_type']}] "
                        f"conf={j['confidence']}\n")
    fig.text(0.02, -0.02, summary_text, fontsize=7, family="monospace",
             va="top", transform=ax.transAxes)

    plt.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(args.output), dpi=150, bbox_inches="tight")
    print(f"[viz] Saved: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
