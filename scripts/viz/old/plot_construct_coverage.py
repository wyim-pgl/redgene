#!/usr/bin/env python3
"""Plot 1: Construct coverage profile.

Generates a depth-along-construct plot, showing which construct elements
are covered by reads, with color-coded regions by element type.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pysam


ELEMENT_COLORS = {
    "promoter": "#4CAF50",
    "terminator": "#FF9800",
    "cds": "#2196F3",
    "vector": "#9E9E9E",
    "element-specific": "#9C27B0",
    "construct-specific": "#F44336",
    "event-specific": "#00BCD4",
    "taxon-specific": "#795548",
}


def get_element_type(ref_name: str) -> str:
    """Extract element type from reference name."""
    parts = ref_name.split("|")
    if parts:
        etype = parts[0].strip()
        if etype in ELEMENT_COLORS:
            return etype
    return "vector"


def get_short_name(ref_name: str) -> str:
    """Get a short display name from reference header."""
    parts = ref_name.split("|")
    if len(parts) >= 2:
        return parts[1]
    return ref_name[:30]


def main():
    parser = argparse.ArgumentParser(description="Plot construct coverage profile")
    parser.add_argument("--bam", type=Path, required=True,
                        help="Construct-mapped BAM")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output plot file (PNG/PDF)")
    parser.add_argument("--sample-name", type=str, default="sample")
    parser.add_argument("--min-depth", type=int, default=1,
                        help="Min depth to include a reference")
    args = parser.parse_args()

    if not args.bam.exists():
        print(f"ERROR: BAM not found: {args.bam}", file=sys.stderr)
        sys.exit(1)

    # Collect depth using samtools depth (more reliable than pileup)
    import subprocess
    result = subprocess.run(
        ["samtools", "depth", "-a", str(args.bam)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"ERROR: samtools depth failed", file=sys.stderr)
        sys.exit(1)

    # Parse depth output: ref\tpos\tdepth
    from collections import defaultdict
    ref_depth_data = defaultdict(list)
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        ref_name, pos, depth = parts[0], int(parts[1]), int(parts[2])
        ref_depth_data[ref_name].append((pos, depth))

    ref_depths = {}
    ref_means = {}
    for ref_name, data in ref_depth_data.items():
        if not data:
            continue
        max_pos = max(p for p, _ in data)
        full_depth = np.zeros(max_pos + 1)
        for pos, depth in data:
            full_depth[pos - 1] = depth  # 1-based to 0-based
        mean_depth = np.mean(full_depth)
        if mean_depth >= args.min_depth:
            ref_depths[ref_name] = full_depth
            ref_means[ref_name] = mean_depth

    if not ref_depths:
        print("No references with sufficient depth", file=sys.stderr)
        sys.exit(0)

    # Sort by mean depth descending
    sorted_refs = sorted(ref_means.keys(), key=lambda r: ref_means[r], reverse=True)

    # Plot: one subplot per reference with depth >= min_depth, max 15
    n_refs = min(len(sorted_refs), 15)
    fig, axes = plt.subplots(n_refs, 1, figsize=(14, max(3 * n_refs, 6)),
                             squeeze=False)
    fig.suptitle(f"Construct Coverage Profile: {args.sample_name}",
                 fontsize=14, fontweight="bold", y=1.02)

    for i, ref_name in enumerate(sorted_refs[:n_refs]):
        ax = axes[i, 0]
        depths = ref_depths[ref_name]
        positions = np.arange(len(depths))
        etype = get_element_type(ref_name)
        color = ELEMENT_COLORS.get(etype, "#607D8B")
        short_name = get_short_name(ref_name)

        ax.fill_between(positions, depths, alpha=0.6, color=color)
        ax.plot(positions, depths, color=color, linewidth=0.5)
        ax.set_ylabel("Depth")
        ax.set_title(f"{short_name} ({etype}, {len(depths)}bp, "
                     f"mean={ref_means[ref_name]:.1f}x)",
                     fontsize=10, loc="left")
        ax.set_xlim(0, len(depths))
        ax.axhline(y=ref_means[ref_name], color="red", linestyle="--",
                   linewidth=0.8, alpha=0.6)

    axes[-1, 0].set_xlabel("Position (bp)")

    # Legend
    legend_patches = [mpatches.Patch(color=c, label=t)
                      for t, c in ELEMENT_COLORS.items()]
    fig.legend(handles=legend_patches, loc="upper right", fontsize=8,
              title="Element type")

    plt.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(args.output), dpi=150, bbox_inches="tight")
    print(f"[viz] Saved: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
