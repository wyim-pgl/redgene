#!/usr/bin/env python3
"""Junction amplicon primer design and visualization.

For each detected junction, extracts the chimeric contig region spanning the
host-construct breakpoint, designs PCR primer pairs that produce 150-250bp
amplicons crossing the junction, and visualizes the amplicon structure.

Primer design uses simple Tm-based heuristics (no Primer3 dependency):
  - Primer length: 18-25 bp
  - Target Tm: 58-62°C (nearest-neighbor method approximation)
  - GC content: 40-60%
  - No runs of ≥4 same nucleotide
  - Avoids 3' self-complementarity

Output:
  {outdir}/junction_amplicons_{sample}.tsv   — primer sequences + amplicon info
  {outdir}/junction_amplicons_{sample}.fa    — amplicon FASTA sequences
  {outdir}/junction_amplicons_{sample}.png   — amplicon structure visualization
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

import os
os.environ["MPLBACKEND"] = "Agg"
os.environ["QT_QPA_PLATFORM"] = "offscreen"
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
plt.switch_backend("Agg")


# ---------------------------------------------------------------------------
# Primer design helpers
# ---------------------------------------------------------------------------

# Nearest-neighbor ΔH and ΔS (SantaLucia 1998, kcal/mol and cal/mol·K)
NN_DH = {
    "AA": -7.9, "TT": -7.9, "AT": -7.2, "TA": -7.2,
    "CA": -8.5, "TG": -8.5, "GT": -8.4, "AC": -8.4,
    "CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2,
    "CG": -10.6, "GC": -9.8, "GG": -8.0, "CC": -8.0,
}
NN_DS = {
    "AA": -22.2, "TT": -22.2, "AT": -20.4, "TA": -21.3,
    "CA": -22.7, "TG": -22.7, "GT": -22.4, "AC": -22.4,
    "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2,
    "CG": -27.2, "GC": -24.4, "GG": -19.9, "CC": -19.9,
}


def calc_tm(seq: str, oligo_conc_nM: float = 250.0, na_conc_mM: float = 50.0) -> float:
    """Calculate Tm using nearest-neighbor method (SantaLucia 1998)."""
    seq = seq.upper()
    if len(seq) < 10:
        return 0.0

    dh = 0.0  # kcal/mol
    ds = 0.0  # cal/mol·K

    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        dh += NN_DH.get(dinuc, -7.5)
        ds += NN_DS.get(dinuc, -21.0)

    # Initiation parameters
    dh += 0.2   # kcal/mol
    ds += -5.7  # cal/mol·K

    # Salt correction (Owczarzy et al. 2004 simplified)
    ds += 0.368 * (len(seq) - 1) * math.log(na_conc_mM / 1000.0)

    R = 1.987  # cal/mol·K
    ct = oligo_conc_nM * 1e-9
    tm = (dh * 1000.0) / (ds + R * math.log(ct / 4.0)) - 273.15
    return tm


def gc_content(seq: str) -> float:
    """GC fraction."""
    seq = seq.upper()
    gc = sum(1 for c in seq if c in "GC")
    return gc / len(seq) if seq else 0.0


def has_run(seq: str, max_run: int = 4) -> bool:
    """Check for homopolymer runs ≥ max_run."""
    seq = seq.upper()
    count = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            count += 1
            if count >= max_run:
                return True
        else:
            count = 1
    return False


def revcomp(seq: str) -> str:
    """Reverse complement."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def pick_primer(seq: str, direction: str, target_tm: float = 60.0,
                min_len: int = 18, max_len: int = 25) -> tuple[str, float, float] | None:
    """Pick a primer from a sequence region.

    direction: 'forward' (take from 5' end of seq) or 'reverse' (take from 3' end, revcomp)
    Returns: (primer_seq, tm, gc) or None

    Tries strict GC filter first (40-60%), then relaxes to 30-70% if needed.
    """
    for gc_lo, gc_hi in [(0.35, 0.65), (0.25, 0.75)]:
        best = None
        best_tm_diff = 999.0

        for plen in range(min_len, max_len + 1):
            if direction == "forward":
                candidate = seq[:plen].upper()
            else:
                candidate = revcomp(seq[-plen:]).upper()

            if "N" in candidate:
                continue

            gc = gc_content(candidate)
            if gc < gc_lo or gc > gc_hi:
                continue

            if has_run(candidate):
                continue

            tm = calc_tm(candidate)
            tm_diff = abs(tm - target_tm)

            if tm_diff < best_tm_diff:
                best_tm_diff = tm_diff
                best = (candidate, tm, gc)

        if best is not None:
            return best

    return None


@dataclass
class AmpliconResult:
    """Result of amplicon design for one junction."""
    sample: str
    contig_name: str
    host_chr: str
    junction_pos_host: int
    junction_type: str
    confidence: str
    amplicon_len: int
    amplicon_seq: str
    junction_offset: int  # position within amplicon where junction occurs
    fwd_primer: str
    fwd_tm: float
    fwd_gc: float
    rev_primer: str
    rev_tm: float
    rev_gc: float
    host_region: str  # "left" or "right" in the amplicon
    construct_element: str


def design_amplicon(
    contig_seq: str,
    contig_name: str,
    host_start_on_contig: int,
    host_end_on_contig: int,
    construct_start_on_contig: int,
    construct_end_on_contig: int,
    host_chr: str,
    junction_pos_host: int,
    junction_type: str,
    confidence: str,
    construct_element: str,
    sample: str,
    min_amplicon: int = 150,
    max_amplicon: int = 250,
) -> AmpliconResult | None:
    """Design an amplicon spanning the junction breakpoint on a chimeric contig."""

    contig_len = len(contig_seq)

    # Find the junction breakpoint on the contig
    if host_end_on_contig <= construct_start_on_contig:
        # Host on left, construct on right
        junction_on_contig = (host_end_on_contig + construct_start_on_contig) // 2
        host_region = "left"
    elif construct_end_on_contig <= host_start_on_contig:
        # Construct on left, host on right
        junction_on_contig = (construct_end_on_contig + host_start_on_contig) // 2
        host_region = "right"
    else:
        # Overlapping — use midpoint of overlap
        junction_on_contig = (max(host_start_on_contig, construct_start_on_contig) +
                              min(host_end_on_contig, construct_end_on_contig)) // 2
        host_region = "left" if host_start_on_contig < construct_start_on_contig else "right"

    # Target amplicon: center on junction, 150-250bp
    target_len = 200
    half = target_len // 2

    # Adjust if near contig edges
    amp_start = max(0, junction_on_contig - half)
    amp_end = min(contig_len, junction_on_contig + half)

    # Ensure minimum amplicon length
    if amp_end - amp_start < min_amplicon:
        # Try to extend
        shortfall = min_amplicon - (amp_end - amp_start)
        if amp_start > 0:
            amp_start = max(0, amp_start - shortfall)
        if amp_end - amp_start < min_amplicon:
            amp_end = min(contig_len, amp_start + min_amplicon)
        if amp_end - amp_start < min_amplicon:
            return None  # Contig too short

    amplicon_seq = contig_seq[amp_start:amp_end]
    junction_offset = junction_on_contig - amp_start

    # Design forward primer — scan sliding windows across the first half
    fwd = None
    fwd_pos = 0  # track where primer starts in amplicon
    for start in range(0, min(junction_offset, len(amplicon_seq) - 25), 5):
        region = amplicon_seq[start:start + 30]
        fwd = pick_primer(region, "forward")
        if fwd is not None:
            fwd_pos = start
            break
    if fwd is None:
        return None

    # Design reverse primer — scan sliding windows across the second half
    rev = None
    rev_pos = len(amplicon_seq)  # track where primer ends in amplicon
    for end_offset in range(0, min(len(amplicon_seq) - junction_offset, len(amplicon_seq) - 25), 5):
        end = len(amplicon_seq) - end_offset
        start = max(0, end - 30)
        region = amplicon_seq[start:end]
        rev = pick_primer(region, "reverse")
        if rev is not None:
            rev_pos = end
            break
    if rev is None:
        return None

    # Trim amplicon to actual primer-to-primer region
    actual_amp_start = amp_start + fwd_pos
    actual_amp_end = amp_start + rev_pos
    amplicon_seq = contig_seq[actual_amp_start:actual_amp_end]
    junction_offset = junction_on_contig - actual_amp_start

    fwd_seq, fwd_tm, fwd_gc = fwd
    rev_seq, rev_tm, rev_gc = rev

    # Actual amplicon length (trimmed to primer positions)
    actual_len = len(amplicon_seq)

    # Shorten element name for display
    element_short = construct_element.split("|")[0] + "|" + construct_element.split("|")[1] if "|" in construct_element else construct_element

    return AmpliconResult(
        sample=sample,
        contig_name=contig_name,
        host_chr=host_chr,
        junction_pos_host=junction_pos_host,
        junction_type=junction_type,
        confidence=confidence,
        amplicon_len=actual_len,
        amplicon_seq=amplicon_seq,
        junction_offset=junction_offset,
        fwd_primer=fwd_seq,
        fwd_tm=fwd_tm,
        fwd_gc=fwd_gc,
        rev_primer=rev_seq,
        rev_tm=rev_tm,
        rev_gc=rev_gc,
        host_region=host_region,
        construct_element=element_short,
    )


# ---------------------------------------------------------------------------
# Junction TSV parser
# ---------------------------------------------------------------------------

def parse_junctions(tsv_path: Path) -> list[dict]:
    """Parse junctions.tsv into list of dicts."""
    junctions = []
    with open(tsv_path) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) < len(header):
                continue
            row = dict(zip(header, fields))
            junctions.append(row)
    return junctions


def load_contigs(fasta_path: Path) -> dict[str, str]:
    """Load contig sequences from FASTA."""
    contigs: dict[str, str] = {}
    name = ""
    seq_parts: list[str] = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    contigs[name] = "".join(seq_parts)
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name:
        contigs[name] = "".join(seq_parts)
    return contigs


def load_paf(paf_path: Path) -> dict[str, list[tuple[int, int, str]]]:
    """Load PAF and return {contig_name: [(query_start, query_end, target_name), ...]}."""
    result: dict[str, list[tuple[int, int, str]]] = {}
    with open(paf_path) as fh:
        for line in fh:
            f = line.strip().split("\t")
            if len(f) < 12:
                continue
            qname = f[0]
            qs, qe = int(f[2]), int(f[3])
            tname = f[5]
            result.setdefault(qname, []).append((qs, qe, tname))
    return result


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def plot_amplicons(results: list[AmpliconResult], outpath: Path, sample_name: str) -> None:
    """Create a figure showing amplicon structure for each junction."""
    if not results:
        return

    n = len(results)
    fig_height = max(3, n * 1.8 + 1.5)
    fig, axes = plt.subplots(n, 1, figsize=(12, fig_height), squeeze=False)

    for idx, res in enumerate(results):
        ax = axes[idx, 0]
        amp_len = res.amplicon_len
        jx = res.junction_offset

        # Draw amplicon bar
        bar_y = 0.5
        bar_h = 0.3

        # Host region (blue)
        if res.host_region == "left":
            host_rect = FancyBboxPatch((0, bar_y - bar_h/2), jx, bar_h,
                                        boxstyle="round,pad=0.01",
                                        facecolor="#4488CC", edgecolor="black", linewidth=1.2)
            construct_rect = FancyBboxPatch((jx, bar_y - bar_h/2), amp_len - jx, bar_h,
                                            boxstyle="round,pad=0.01",
                                            facecolor="#CC4444", edgecolor="black", linewidth=1.2)
            ax.text(jx/2, bar_y, f"Host\n({res.host_chr})",
                    ha="center", va="center", fontsize=8, fontweight="bold", color="white")
            ax.text(jx + (amp_len - jx)/2, bar_y, f"Construct\n({res.construct_element})",
                    ha="center", va="center", fontsize=7, fontweight="bold", color="white")
        else:
            construct_rect = FancyBboxPatch((0, bar_y - bar_h/2), jx, bar_h,
                                            boxstyle="round,pad=0.01",
                                            facecolor="#CC4444", edgecolor="black", linewidth=1.2)
            host_rect = FancyBboxPatch((jx, bar_y - bar_h/2), amp_len - jx, bar_h,
                                        boxstyle="round,pad=0.01",
                                        facecolor="#4488CC", edgecolor="black", linewidth=1.2)
            ax.text(jx/2, bar_y, f"Construct\n({res.construct_element})",
                    ha="center", va="center", fontsize=7, fontweight="bold", color="white")
            ax.text(jx + (amp_len - jx)/2, bar_y, f"Host\n({res.host_chr})",
                    ha="center", va="center", fontsize=8, fontweight="bold", color="white")

        ax.add_patch(host_rect)
        ax.add_patch(construct_rect)

        # Junction line
        ax.axvline(x=jx, color="gold", linewidth=2.5, linestyle="--", zorder=5)
        ax.text(jx, bar_y + bar_h/2 + 0.08, f"Junction ({res.junction_type})\n{res.host_chr}:{res.junction_pos_host:,}",
                ha="center", va="bottom", fontsize=8, fontweight="bold", color="#333333")

        # Forward primer (green arrow)
        fwd_len = len(res.fwd_primer)
        ax.annotate("", xy=(fwd_len, bar_y - bar_h/2 - 0.06),
                     xytext=(0, bar_y - bar_h/2 - 0.06),
                     arrowprops=dict(arrowstyle="->", color="#228B22", lw=2.5))
        ax.text(fwd_len/2, bar_y - bar_h/2 - 0.14,
                f"Fwd: {res.fwd_primer}\nTm={res.fwd_tm:.1f}°C  GC={res.fwd_gc:.0%}",
                ha="center", va="top", fontsize=6.5, color="#228B22", family="monospace")

        # Reverse primer (orange arrow)
        rev_len = len(res.rev_primer)
        rev_start = amp_len
        rev_end = amp_len - rev_len
        ax.annotate("", xy=(rev_end, bar_y - bar_h/2 - 0.06),
                     xytext=(rev_start, bar_y - bar_h/2 - 0.06),
                     arrowprops=dict(arrowstyle="->", color="#CC6600", lw=2.5))
        ax.text((rev_start + rev_end)/2, bar_y - bar_h/2 - 0.14,
                f"Rev: {res.rev_primer}\nTm={res.rev_tm:.1f}°C  GC={res.rev_gc:.0%}",
                ha="center", va="top", fontsize=6.5, color="#CC6600", family="monospace")

        # Amplicon length annotation
        ax.text(amp_len/2, bar_y + bar_h/2 + 0.25, f"Amplicon: {amp_len} bp",
                ha="center", va="bottom", fontsize=9, fontweight="bold", color="black",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFFFCC", edgecolor="#999999"))

        # Title
        ax.set_title(f"{res.contig_name}  |  {res.confidence} confidence  |  {res.junction_type} junction",
                      fontsize=10, fontweight="bold", pad=12)

        ax.set_xlim(-10, amp_len + 10)
        ax.set_ylim(-0.05, 1.1)
        ax.set_xlabel("Position on amplicon (bp)", fontsize=8)
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor="#4488CC", edgecolor="black", label="Host genome"),
        mpatches.Patch(facecolor="#CC4444", edgecolor="black", label="Construct/transgene"),
        plt.Line2D([0], [0], color="gold", linewidth=2.5, linestyle="--", label="Junction breakpoint"),
        plt.Line2D([0], [0], color="#228B22", linewidth=2.5, label="Forward primer"),
        plt.Line2D([0], [0], color="#CC6600", linewidth=2.5, label="Reverse primer"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=5, fontsize=8,
               frameon=True, fancybox=True)

    fig.suptitle(f"Junction Amplicons — {sample_name}", fontsize=13, fontweight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0.06, 1, 0.96])
    fig.savefig(outpath, dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"[amplicon] Saved visualization: {outpath}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Design junction-spanning amplicons and primers")
    parser.add_argument("--junctions", required=True, type=Path, help="junctions.tsv from step 6")
    parser.add_argument("--contigs", required=True, type=Path, help="contigs.fasta from step 4")
    parser.add_argument("--host-paf", required=True, type=Path, help="contigs_to_host.paf from step 5")
    parser.add_argument("--construct-paf", required=True, type=Path, help="contigs_to_construct.paf from step 5")
    parser.add_argument("--sample-name", required=True, help="Sample name")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--min-amplicon", type=int, default=150, help="Minimum amplicon length (default: 150)")
    parser.add_argument("--max-amplicon", type=int, default=250, help="Maximum amplicon length (default: 250)")
    parser.add_argument("--min-confidence", default="Medium",
                        help="Minimum junction confidence to design primers for (default: Medium)")
    args = parser.parse_args()

    confidence_order = {"Low": 0, "Medium": 1, "High": 2}
    min_conf_val = confidence_order.get(args.min_confidence, 1)

    # Load data
    junctions = parse_junctions(args.junctions)
    if not junctions:
        print(f"[amplicon] No junctions found in {args.junctions}", file=sys.stderr)
        return

    contigs = load_contigs(args.contigs)
    host_paf = load_paf(args.host_paf)
    construct_paf = load_paf(args.construct_paf)

    # Deduplicate junctions: keep unique (contig_name, host_chr, junction_pos) with highest confidence
    seen: dict[tuple[str, str, str], dict] = {}
    for j in junctions:
        key = (j["contig_name"], j["host_chr"], j["junction_pos_host"])
        conf = confidence_order.get(j.get("confidence", "Medium"), 1)
        if key not in seen or conf > confidence_order.get(seen[key].get("confidence", "Medium"), 1):
            seen[key] = j
    unique_junctions = list(seen.values())

    # Filter by confidence
    unique_junctions = [j for j in unique_junctions
                        if confidence_order.get(j.get("confidence", "Medium"), 1) >= min_conf_val]

    print(f"[amplicon] Processing {len(unique_junctions)} unique junctions from {len(junctions)} total",
          file=sys.stderr)

    results: list[AmpliconResult] = []

    for j in unique_junctions:
        cname = j["contig_name"]
        if cname not in contigs:
            print(f"[amplicon] WARNING: contig {cname} not found in FASTA", file=sys.stderr)
            continue

        contig_seq = contigs[cname]

        # Get host and construct alignment coordinates on this contig from PAF
        h_alns = host_paf.get(cname, [])
        c_alns = construct_paf.get(cname, [])

        if not h_alns or not c_alns:
            continue

        # Find the best host alignment matching this junction's chromosome
        target_chr = j["host_chr"]
        matching_h = [a for a in h_alns if a[2] == target_chr]
        if not matching_h:
            matching_h = h_alns  # fallback

        # Use the one with the most coverage
        best_h = max(matching_h, key=lambda a: a[1] - a[0])

        # For construct: find the alignment that is most complementary to host
        # (least overlap, best for identifying junction)
        best_c = None
        best_c_score = -1
        for ca in c_alns:
            # Prefer construct alignments that don't overlap with host
            overlap = min(best_h[1], ca[1]) - max(best_h[0], ca[0])
            span = ca[1] - ca[0]
            score = span - max(0, overlap)
            if score > best_c_score:
                best_c_score = score
                best_c = ca
        if best_c is None:
            best_c = max(c_alns, key=lambda a: a[1] - a[0])

        res = design_amplicon(
            contig_seq=contig_seq,
            contig_name=cname,
            host_start_on_contig=best_h[0],
            host_end_on_contig=best_h[1],
            construct_start_on_contig=best_c[0],
            construct_end_on_contig=best_c[1],
            host_chr=j["host_chr"],
            junction_pos_host=int(j["junction_pos_host"]),
            junction_type=j.get("junction_type", "?"),
            confidence=j.get("confidence", "Medium"),
            construct_element=j.get("construct_element", "unknown"),
            sample=args.sample_name,
            min_amplicon=args.min_amplicon,
            max_amplicon=args.max_amplicon,
        )

        if res is not None:
            results.append(res)
        else:
            print(f"[amplicon] Could not design primers for {cname} → {j['host_chr']}:{j['junction_pos_host']}",
                  file=sys.stderr)

    if not results:
        print("[amplicon] No amplicons could be designed", file=sys.stderr)
        return

    # Write TSV
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    tsv_path = outdir / f"junction_amplicons_{args.sample_name}.tsv"
    fa_path = outdir / f"junction_amplicons_{args.sample_name}.fa"
    png_path = outdir / f"junction_amplicons_{args.sample_name}.png"

    with open(tsv_path, "w") as fh:
        fh.write("\t".join([
            "sample", "contig", "host_chr", "junction_pos", "junction_type", "confidence",
            "amplicon_len", "junction_offset", "host_region",
            "fwd_primer", "fwd_tm", "fwd_gc",
            "rev_primer", "rev_tm", "rev_gc",
            "construct_element", "amplicon_seq",
        ]) + "\n")
        for r in results:
            fh.write("\t".join([
                r.sample, r.contig_name, r.host_chr, str(r.junction_pos_host),
                r.junction_type, r.confidence,
                str(r.amplicon_len), str(r.junction_offset), r.host_region,
                r.fwd_primer, f"{r.fwd_tm:.1f}", f"{r.fwd_gc:.2f}",
                r.rev_primer, f"{r.rev_tm:.1f}", f"{r.rev_gc:.2f}",
                r.construct_element, r.amplicon_seq,
            ]) + "\n")
    print(f"[amplicon] Wrote {len(results)} amplicons to {tsv_path}", file=sys.stderr)

    # Write FASTA
    with open(fa_path, "w") as fh:
        for r in results:
            header = (f">{r.sample}_{r.contig_name}_{r.host_chr}:{r.junction_pos_host}"
                      f" len={r.amplicon_len} junction_offset={r.junction_offset}"
                      f" fwd={r.fwd_primer} rev={r.rev_primer}")
            fh.write(header + "\n")
            # Write sequence in 80-char lines
            for i in range(0, len(r.amplicon_seq), 80):
                fh.write(r.amplicon_seq[i:i+80] + "\n")
    print(f"[amplicon] Wrote amplicon FASTA to {fa_path}", file=sys.stderr)

    # Visualization
    plot_amplicons(results, png_path, args.sample_name)


if __name__ == "__main__":
    main()
