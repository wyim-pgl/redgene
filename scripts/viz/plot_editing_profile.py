#!/usr/bin/env python3
"""CRISPR editing profile visualization.

Generates a multi-panel figure showing per-position nucleotide composition,
read depth, and indel frequencies around gRNA target sites. Inspired by
CRISPResso2's nucleotide quilt approach, extended with integrated depth
tracks and modification highlighting.

Panels (top to bottom):
  1. Read depth bar chart with coverage gradient
  2. Nucleotide composition heatmap (A/T/G/C/del/ins frequencies)
  3. Modification frequency track (% reads with non-reference base)
  4. Reference sequence with gRNA, PAM, and cut site annotation

Usage:
  python plot_editing_profile.py \
    --treatment-bam results/{sample}/s04_host_map/{sample}_host.bam \
    --wt-bam results/WT/s04_host_map/WT_host.bam \
    --host-ref db/host.fa \
    --grna-targets results/{sample}/s06_indel/grna_targets.tsv \
    --editing-sites results/{sample}/s06_indel/editing_sites.tsv \
    --sample-name {sample} \
    --outdir results/{sample}/s06_indel
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import numpy as np


# ---------------------------------------------------------------------------
# CRISPResso2-inspired color scheme for nucleotides
# ---------------------------------------------------------------------------

NUC_COLORS = {
    "A": (127 / 255, 201 / 255, 127 / 255),   # green
    "T": (190 / 255, 174 / 255, 212 / 255),   # purple
    "G": (255 / 255, 255 / 255, 153 / 255),   # yellow
    "C": (253 / 255, 192 / 255, 134 / 255),   # orange
    "N": (200 / 255, 200 / 255, 200 / 255),   # gray
    "DEL": (0.12, 0.12, 0.12),                 # dark gray
    "INS": (193 / 255, 129 / 255, 114 / 255), # brown
}

# Modification highlight colors
MOD_COLOR = "#E53935"       # red for modifications
REF_MATCH_COLOR = "#E0E0E0" # light gray for reference-matching bases
GRNA_COLOR = "#2196F3"      # blue for gRNA
PAM_COLOR = "#FF5722"       # deep orange for PAM
CUT_COLOR = "#D32F2F"       # red for cut site


def log(msg: str) -> None:
    print(f"[plot_editing] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Pileup parsing
# ---------------------------------------------------------------------------

def parse_pileup_at_region(
    bam: Path, ref: Path, chrom: str, start: int, end: int,
) -> dict:
    """Parse samtools mpileup to get per-position base counts and indels.

    Returns dict keyed by position with:
      {pos: {"A": n, "T": n, "G": n, "C": n, "DEL": n, "INS": n,
             "dp": n, "ref_base": str, "ins_seqs": {seq: count},
             "del_seqs": {seq: count}}}
    """
    region = f"{chrom}:{start}-{end}"
    cmd = [
        "samtools", "mpileup", "-f", str(ref),
        "-r", region, "-d", "10000", "-q", "10", "-Q", "0",
        str(bam),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    data = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue

        pos = int(fields[1])
        ref_base = fields[2].upper()
        dp = int(fields[3])
        pileup_str = fields[4]

        counts = {"A": 0, "T": 0, "G": 0, "C": 0, "DEL": 0, "INS": 0}
        ins_seqs: dict[str, int] = {}
        del_seqs: dict[str, int] = {}

        i = 0
        while i < len(pileup_str):
            c = pileup_str[i]
            if c == "^":
                i += 2  # skip ^ + mapq char
                continue
            elif c == "$":
                i += 1
                continue
            elif c in ".,":
                counts[ref_base] += 1
                i += 1
            elif c.upper() in "ATGCN":
                counts[c.upper()] += 1
                i += 1
            elif c == "*":
                counts["DEL"] += 1
                i += 1
            elif c in "+-":
                # Indel notation without preceding base (shouldn't happen
                # in valid pileup, but skip the indel length + sequence)
                i += 1
                num_str = ""
                while i < len(pileup_str) and pileup_str[i].isdigit():
                    num_str += pileup_str[i]
                    i += 1
                if num_str:
                    i += int(num_str)
                continue
            else:
                i += 1
                continue

            # Check for indel after base
            if i < len(pileup_str) and pileup_str[i] in "+-":
                c2 = pileup_str[i]
                is_ins = c2 == "+"
                i += 1
                num_str = ""
                while i < len(pileup_str) and pileup_str[i].isdigit():
                    num_str += pileup_str[i]
                    i += 1
                if num_str:
                    n = int(num_str)
                    seq = pileup_str[i:i + n].upper()
                    i += n
                    if is_ins:
                        counts["INS"] += 1
                        ins_seqs[seq] = ins_seqs.get(seq, 0) + 1
                    else:
                        del_seqs[seq] = del_seqs.get(seq, 0) + 1

        data[pos] = {
            **counts,
            "dp": dp,
            "ref_base": ref_base,
            "ins_seqs": ins_seqs,
            "del_seqs": del_seqs,
        }

    return data


def get_ref_seq(ref: Path, chrom: str, start: int, end: int) -> str:
    """Get reference sequence for a region."""
    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        ["samtools", "faidx", str(ref), region],
        capture_output=True, text=True,
    )
    return "".join(result.stdout.strip().split("\n")[1:]).upper()


# ---------------------------------------------------------------------------
# Main plotting function
# ---------------------------------------------------------------------------

def plot_editing_profile(
    treatment_data: dict,
    wt_data: dict,
    ref_seq: str,
    chrom: str,
    start: int,
    end: int,
    grna_start: int | None,
    grna_end: int | None,
    grna_strand: str | None,
    pam_start: int | None,
    pam_end: int | None,
    cut_pos: int | None,
    editing_sites: list[dict],
    sample_name: str,
    output_path: Path,
    window: int = 30,
) -> None:
    """Generate the multi-panel editing profile figure."""

    # Determine plot range centered on cut site
    if cut_pos:
        plot_start = max(start, cut_pos - window)
        plot_end = min(end, cut_pos + window)
    else:
        plot_start = start
        plot_end = end

    positions = list(range(plot_start, plot_end + 1))
    n_pos = len(positions)

    if n_pos == 0:
        log("WARNING: No positions to plot")
        return

    # Extract data arrays
    t_depths = []
    wt_depths = []
    t_freqs = {b: [] for b in ["A", "T", "G", "C", "DEL", "INS"]}
    wt_freqs = {b: [] for b in ["A", "T", "G", "C", "DEL", "INS"]}
    mod_freqs = []  # modification frequency (treatment)
    ref_bases = []

    for pos in positions:
        t_d = treatment_data.get(pos, {"A": 0, "T": 0, "G": 0, "C": 0,
                                        "DEL": 0, "INS": 0, "dp": 0,
                                        "ref_base": "N"})
        w_d = wt_data.get(pos, {"A": 0, "T": 0, "G": 0, "C": 0,
                                 "DEL": 0, "INS": 0, "dp": 0,
                                 "ref_base": "N"})

        t_dp = max(t_d["dp"], 1)
        w_dp = max(w_d["dp"], 1)
        t_depths.append(t_d["dp"])
        wt_depths.append(w_d["dp"])

        rb = t_d["ref_base"]
        ref_bases.append(rb)

        for b in ["A", "T", "G", "C", "DEL", "INS"]:
            t_freqs[b].append(t_d.get(b, 0) / t_dp)
            wt_freqs[b].append(w_d.get(b, 0) / w_dp)

        # Modification frequency = fraction of non-reference reads
        ref_count = t_d.get(rb, 0) if rb in "ATGC" else 0
        mod_frac = 1.0 - (ref_count / t_dp) if t_dp > 0 else 0
        mod_freqs.append(mod_frac)

    t_depths = np.array(t_depths)
    wt_depths = np.array(wt_depths)
    mod_freqs = np.array(mod_freqs)

    # -----------------------------------------------------------------------
    # Create figure with GridSpec
    # -----------------------------------------------------------------------
    fig = plt.figure(figsize=(max(14, n_pos * 0.25), 12))
    gs = gridspec.GridSpec(
        5, 1,
        height_ratios=[1.5, 2.5, 2.5, 1.0, 0.8],
        hspace=0.08,
    )

    x = np.arange(n_pos)
    x_labels = [str(p) for p in positions]

    # Panel 1: Read depth (treatment + WT overlay)
    ax_depth = fig.add_subplot(gs[0])
    ax_depth.bar(x, t_depths, width=0.8, color="#1565C0", alpha=0.8,
                 label="Treatment", zorder=3)
    ax_depth.step(x, wt_depths, where="mid", color="#E65100", alpha=0.7,
                  linewidth=1.5, label="WT", zorder=4)
    ax_depth.set_ylabel("Read\nDepth", fontsize=9, rotation=0,
                        ha="right", va="center")
    ax_depth.set_xlim(-0.5, n_pos - 0.5)
    ax_depth.legend(loc="upper right", fontsize=7, framealpha=0.8)
    ax_depth.set_title(
        f"CRISPR Editing Profile — {sample_name}\n"
        f"{chrom}:{plot_start}-{plot_end}",
        fontsize=11, fontweight="bold", pad=10,
    )
    ax_depth.tick_params(labelbottom=False)
    ax_depth.spines["bottom"].set_visible(False)

    # Panel 2: Nucleotide composition — Treatment (stacked bars)
    ax_nuc_t = fig.add_subplot(gs[1])
    _draw_nucleotide_quilt(ax_nuc_t, x, n_pos, t_freqs, ref_bases,
                           positions, grna_start, grna_end, pam_start,
                           pam_end, cut_pos, "Treatment")
    ax_nuc_t.tick_params(labelbottom=False)

    # Panel 3: Nucleotide composition — WT (stacked bars)
    ax_nuc_w = fig.add_subplot(gs[2])
    _draw_nucleotide_quilt(ax_nuc_w, x, n_pos, wt_freqs, ref_bases,
                           positions, grna_start, grna_end, pam_start,
                           pam_end, cut_pos, "WT")
    ax_nuc_w.tick_params(labelbottom=False)

    # Panel 4: Modification frequency (treatment-specific)
    ax_mod = fig.add_subplot(gs[3])

    # Compute WT mod freq for subtraction
    wt_mod = []
    for i, pos in enumerate(positions):
        w_d = wt_data.get(pos, {"A": 0, "T": 0, "G": 0, "C": 0,
                                 "DEL": 0, "INS": 0, "dp": 0,
                                 "ref_base": "N"})
        w_dp = max(w_d["dp"], 1)
        rb = ref_bases[i]
        ref_count = w_d.get(rb, 0) if rb in "ATGC" else 0
        wt_mod.append(1.0 - (ref_count / w_dp))
    wt_mod = np.array(wt_mod)

    # Net modification = treatment - WT (clamp to 0)
    net_mod = np.maximum(mod_freqs - wt_mod, 0)

    # Color bars by significance
    bar_colors = []
    for i in range(n_pos):
        if net_mod[i] > 0.05 and t_depths[i] >= 3:
            bar_colors.append(MOD_COLOR)
        else:
            bar_colors.append("#BDBDBD")

    ax_mod.bar(x, net_mod * 100, width=0.8, color=bar_colors, zorder=3)
    ax_mod.axhline(y=5, color="#999", linestyle="--", linewidth=0.5, zorder=2)
    ax_mod.set_ylabel("Mod\n%", fontsize=9, rotation=0, ha="right", va="center")
    ax_mod.set_ylim(0, max(net_mod * 100) * 1.2 + 1 if max(net_mod) > 0 else 10)
    ax_mod.set_xlim(-0.5, n_pos - 0.5)
    ax_mod.tick_params(labelbottom=False)

    # Mark editing sites with arrows
    for site in editing_sites:
        site_pos = site.get("pos", 0)
        if plot_start <= site_pos <= plot_end:
            idx = site_pos - plot_start
            site_type = site.get("type", "")
            site_size = site.get("size", 0)
            indel_seq = site.get("indel_seq", "")
            label = f"{site_type[:3]}{site_size}bp"
            if indel_seq:
                label += f"\n{indel_seq}"
            ax_mod.annotate(
                label,
                xy=(idx, net_mod[idx] * 100),
                xytext=(idx, net_mod[idx] * 100 + 8),
                fontsize=6, fontweight="bold", color=MOD_COLOR,
                ha="center", va="bottom",
                arrowprops=dict(arrowstyle="->", color=MOD_COLOR, lw=1),
            )

    # Panel 5: Reference sequence annotation
    ax_ref = fig.add_subplot(gs[4])
    ax_ref.set_xlim(-0.5, n_pos - 0.5)
    ax_ref.set_ylim(0, 1)
    ax_ref.axis("off")

    # Draw reference bases as colored boxes
    for i, base in enumerate(ref_bases):
        pos = positions[i]
        # Background: gRNA or PAM region
        bg_color = "white"
        if grna_start and grna_end and grna_start <= pos <= grna_end:
            bg_color = "#BBDEFB"  # light blue
        if pam_start and pam_end and pam_start <= pos <= pam_end:
            bg_color = "#FFCCBC"  # light orange

        rect = mpatches.FancyBboxPatch(
            (i - 0.4, 0.15), 0.8, 0.7,
            boxstyle="round,pad=0.05",
            facecolor=bg_color, edgecolor="#ccc", linewidth=0.5,
        )
        ax_ref.add_patch(rect)

        # Base letter with nucleotide color
        color = NUC_COLORS.get(base, NUC_COLORS["N"])
        ax_ref.text(i, 0.5, base, ha="center", va="center",
                    fontsize=7, fontweight="bold",
                    fontfamily="monospace",
                    color="black",
                    bbox=dict(facecolor=(*color, 0.6), edgecolor="none",
                              pad=0.8, boxstyle="round"))

        # Cut site marker
        if cut_pos and pos == cut_pos:
            ax_ref.axvline(x=i, color=CUT_COLOR, linewidth=2,
                           linestyle="-", zorder=10)
            ax_ref.text(i, 1.05, "✂", ha="center", va="bottom",
                        fontsize=10, color=CUT_COLOR)

    # Position labels (every 5th)
    tick_positions = []
    tick_labels = []
    for i, pos in enumerate(positions):
        if pos % 5 == 0 or (cut_pos and pos == cut_pos):
            tick_positions.append(i)
            tick_labels.append(str(pos))

    ax_mod.set_xticks(tick_positions)
    ax_mod.set_xticklabels(tick_labels, fontsize=6, rotation=45)
    ax_mod.tick_params(labelbottom=True)

    # Legend annotations
    legend_elements = [
        mpatches.Patch(facecolor="#BBDEFB", edgecolor="#999",
                       label="gRNA target"),
        mpatches.Patch(facecolor="#FFCCBC", edgecolor="#999",
                       label="PAM"),
        plt.Line2D([0], [0], color=CUT_COLOR, linewidth=2,
                   label="Cut site"),
        mpatches.Patch(facecolor=MOD_COLOR, label="Significant mod (>5%)"),
    ]
    ax_ref.legend(handles=legend_elements, loc="upper center",
                  ncol=4, fontsize=7, framealpha=0.8,
                  bbox_to_anchor=(0.5, -0.1))

    plt.savefig(output_path, dpi=200, bbox_inches="tight",
                facecolor="white")
    plt.close()
    log(f"Saved: {output_path}")


def _draw_nucleotide_quilt(
    ax, x, n_pos, freqs, ref_bases, positions,
    grna_start, grna_end, pam_start, pam_end, cut_pos,
    label: str,
) -> None:
    """Draw stacked nucleotide composition bars (CRISPResso2 quilt style)."""

    base_order = ["A", "T", "G", "C", "DEL", "INS"]
    bottoms = np.zeros(n_pos)

    for base in base_order:
        vals = np.array(freqs[base])
        color = NUC_COLORS[base]

        # Dim reference-matching bases (CRISPResso2 shade_unchanged style)
        bar_colors = []
        bar_alphas = []
        for i in range(n_pos):
            if base == ref_bases[i]:
                bar_colors.append((*color, 0.3))  # dim reference match
            else:
                bar_colors.append((*color, 1.0))  # full color for variants
            bar_alphas.append(1.0)

        # Draw bars one by one for per-bar color control
        for i in range(n_pos):
            if vals[i] > 0.001:
                ax.bar(x[i], vals[i], width=0.85, bottom=bottoms[i],
                       color=bar_colors[i], edgecolor="none",
                       zorder=3)
        bottoms += vals

    # gRNA and PAM shading
    for i, pos in enumerate(positions):
        if grna_start and grna_end and grna_start <= pos <= grna_end:
            ax.axvspan(i - 0.45, i + 0.45, alpha=0.08,
                       color=GRNA_COLOR, zorder=1)
        if pam_start and pam_end and pam_start <= pos <= pam_end:
            ax.axvspan(i - 0.45, i + 0.45, alpha=0.08,
                       color=PAM_COLOR, zorder=1)

    # Cut site line
    if cut_pos:
        for i, pos in enumerate(positions):
            if pos == cut_pos:
                ax.axvline(x=i, color=CUT_COLOR, linewidth=1.5,
                           linestyle="--", alpha=0.6, zorder=5)

    ax.set_ylabel(f"{label}\nBase\nFreq", fontsize=8, rotation=0,
                  ha="right", va="center")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(-0.5, n_pos - 0.5)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=7)

    # Nucleotide legend
    nuc_legend = [
        mpatches.Patch(facecolor=(*NUC_COLORS["A"], 1), label="A"),
        mpatches.Patch(facecolor=(*NUC_COLORS["T"], 1), label="T"),
        mpatches.Patch(facecolor=(*NUC_COLORS["G"], 1), label="G"),
        mpatches.Patch(facecolor=(*NUC_COLORS["C"], 1), label="C"),
        mpatches.Patch(facecolor=(*NUC_COLORS["DEL"], 1), label="Del"),
        mpatches.Patch(facecolor=(*NUC_COLORS["INS"], 1), label="Ins"),
    ]
    ax.legend(handles=nuc_legend, loc="upper right", ncol=6,
              fontsize=6, framealpha=0.8)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="CRISPR editing profile visualization",
    )
    parser.add_argument("--treatment-bam", required=True, type=Path)
    parser.add_argument("--wt-bam", required=True, type=Path)
    parser.add_argument("--host-ref", required=True, type=Path)
    parser.add_argument("--grna-targets", required=True, type=Path,
                        help="grna_targets.tsv from step 8")
    parser.add_argument("--editing-sites", type=Path, default=None,
                        help="editing_sites.tsv from step 8")
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--window", type=int, default=30,
                        help="Window size around cut site (default: 30bp)")
    args = parser.parse_args()

    outdir = args.outdir / args.sample_name / "s06_indel"
    outdir.mkdir(parents=True, exist_ok=True)

    # Load gRNA targets (on-target sites only)
    targets = []
    with open(args.grna_targets) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("site_type") == "on-target":
                targets.append(row)

    if not targets:
        log("No on-target sites found in grna_targets.tsv")
        sys.exit(0)

    # Load editing sites
    editing_sites = []
    if args.editing_sites and args.editing_sites.exists():
        with open(args.editing_sites) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                try:
                    row["pos"] = int(row["pos"])
                    row["size"] = int(row.get("size", 0))
                except (ValueError, KeyError):
                    pass
                editing_sites.append(row)

    log(f"Generating editing profiles for {len(targets)} on-target site(s)")

    for target in targets:
        chrom = target["chrom"]
        grna_start = int(target["start"]) + 1  # 0-based to 1-based
        grna_end = int(target["end"])
        strand = target["strand"]
        cut_pos = int(target["cut_pos"])
        grna_idx = target.get("grna_idx", "1")
        grna_seq = target.get("grna_seq", "")

        # PAM position
        if strand == "+":
            pam_start = grna_end + 1
            pam_end = grna_end + 3
        else:
            pam_start = grna_start - 3
            pam_end = grna_start - 1

        # Region to plot
        region_start = max(1, cut_pos - args.window)
        region_end = cut_pos + args.window

        log(f"  gRNA {grna_idx} ({grna_seq}): {chrom}:{region_start}-{region_end}")

        # Parse pileups
        t_data = parse_pileup_at_region(
            args.treatment_bam, args.host_ref, chrom, region_start, region_end,
        )
        wt_data = parse_pileup_at_region(
            args.wt_bam, args.host_ref, chrom, region_start, region_end,
        )

        # Reference sequence
        ref_seq = get_ref_seq(args.host_ref, chrom, region_start, region_end)

        # Filter editing sites for this target
        target_edits = [
            s for s in editing_sites
            if s.get("grna_idx") == grna_idx
            and s.get("chrom") == chrom
        ]

        # Output path
        out_png = outdir / f"editing_profile_gRNA{grna_idx}.png"
        out_pdf = outdir / f"editing_profile_gRNA{grna_idx}.pdf"

        plot_editing_profile(
            treatment_data=t_data,
            wt_data=wt_data,
            ref_seq=ref_seq,
            chrom=chrom,
            start=region_start,
            end=region_end,
            grna_start=grna_start,
            grna_end=grna_end,
            grna_strand=strand,
            pam_start=pam_start,
            pam_end=pam_end,
            cut_pos=cut_pos,
            editing_sites=target_edits,
            sample_name=f"{args.sample_name} (gRNA{grna_idx}: {grna_seq})",
            output_path=out_png,
            window=args.window,
        )

        # Also save PDF
        plot_editing_profile(
            treatment_data=t_data,
            wt_data=wt_data,
            ref_seq=ref_seq,
            chrom=chrom,
            start=region_start,
            end=region_end,
            grna_start=grna_start,
            grna_end=grna_end,
            grna_strand=strand,
            pam_start=pam_start,
            pam_end=pam_end,
            cut_pos=cut_pos,
            editing_sites=target_edits,
            sample_name=f"{args.sample_name} (gRNA{grna_idx}: {grna_seq})",
            output_path=out_pdf,
            window=args.window,
        )

    log("All editing profiles generated.")


if __name__ == "__main__":
    main()
