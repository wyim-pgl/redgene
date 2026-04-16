#!/usr/bin/env python3
"""Aggregate s05 assembly round statistics across all completed samples.

Emits two complementary views that inform the T3 ``max_rounds`` default:

* **Per-round growth contribution (from stderr logs)** — parses
  ``Round N growth: kmer=A, mm2=B, pilon=C, ssake=D`` lines written by
  ``scripts/s05_insert_assembly.py`` and sums the combined growth
  (kmer+mm2+pilon+ssake) at each round across all insertion sites.
  The decision signal is ``fraction_of_all`` for rounds > 3.

* **Per-site final round (from s05_stats.txt)** — reads the
  ``<site>_assembly_rounds`` entries to histogram how many rounds each site
  ultimately consumed (0 = merged during initial extension; 8 = hit
  ``max_rounds`` ceiling).

Both views are printed as TSV blocks on stdout so the call-site can pipe to
a file or a Markdown table. Re-run whenever new s05 logs/stats are available.

Usage:
    python scripts/measure_assembly_rounds.py [--results-dir results]
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

# Per-round growth log line emitted by s05_insert_assembly.py ~line 2573.
# Example: "[s05_insert_assembly]     Round 1 growth: kmer=249, mm2=1, pilon=0, ssake=0"
ROUND_RE = re.compile(
    r"Round\s+(\d+)\s+growth:\s*kmer=(\d+),\s*mm2=(\d+),\s*pilon=(\d+),\s*ssake=(\d+)",
    re.IGNORECASE,
)
SITE_RE = re.compile(r"=== Processing (insertion_\S+):")

# s05_stats.txt entry: "<site_id>_assembly_rounds\t<N>"
STATS_ROUND_RE = re.compile(r"^(insertion_\S+)_assembly_rounds\t(\d+)\s*$")


def parse_err_log(path: Path) -> list[tuple[str, int, int]]:
    """Return list of (sample_site_key, round_idx, combined_growth) tuples.

    The sample is inferred from the filename prefix is not reliable, so we
    just tag each row with the log file stem + current site header.
    """
    rows: list[tuple[str, int, int]] = []
    current_site = "unknown"
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return rows
    for line in text.splitlines():
        ms = SITE_RE.search(line)
        if ms:
            current_site = ms.group(1)
            continue
        mr = ROUND_RE.search(line)
        if mr:
            rnd = int(mr.group(1))
            combined = int(mr.group(2)) + int(mr.group(3)) + int(mr.group(4)) + int(mr.group(5))
            rows.append((f"{path.name}::{current_site}", rnd, combined))
    return rows


def parse_stats_file(path: Path) -> list[tuple[str, int]]:
    """Return [(site_id, final_round), ...] from an s05_stats.txt file."""
    out: list[tuple[str, int]] = []
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return out
    for line in text.splitlines():
        m = STATS_ROUND_RE.match(line)
        if m:
            out.append((m.group(1), int(m.group(2))))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results"),
        help="Directory containing per-sample s05 outputs + *.err logs",
    )
    args = ap.parse_args()

    results = args.results_dir
    if not results.is_dir():
        print(f"ERROR: results dir not found: {results}", file=sys.stderr)
        sys.exit(2)

    # ---- View A: per-round growth from stderr logs ---------------------
    err_logs = sorted(results.glob("*.err"))
    per_round_growth: dict[int, int] = defaultdict(int)
    per_round_nsites: dict[int, int] = defaultdict(int)
    sites_seen: set[str] = set()
    for log_path in err_logs:
        for site_key, rnd, combined in parse_err_log(log_path):
            per_round_growth[rnd] += combined
            if site_key not in sites_seen:
                sites_seen.add(site_key)
            # We count every (site, round) observation for per-site averages.
            per_round_nsites[rnd] += 1
    grand_total = sum(per_round_growth.values()) or 1

    print(f"# s05 per-round growth aggregate  (source: {len(err_logs)} *.err logs)")
    print(f"# distinct site-rounds parsed: {sum(per_round_nsites.values())}")
    print(f"# distinct insertion sites seen: {len(sites_seen)}")
    print("round\tn_site_rounds\tgrowth_bp_total\tfraction_of_all\tcumulative_fraction")
    cum = 0
    for rnd in sorted(per_round_growth):
        cum += per_round_growth[rnd]
        frac = per_round_growth[rnd] / grand_total
        cfrac = cum / grand_total
        print(
            f"{rnd}\t{per_round_nsites[rnd]}\t{per_round_growth[rnd]}\t"
            f"{frac:.4f}\t{cfrac:.4f}"
        )

    # ---- View B: per-site final round from s05_stats.txt ---------------
    stat_files = sorted(results.glob("*/s05_insert_assembly/s05_stats.txt"))
    final_round_hist: dict[int, int] = defaultdict(int)
    per_sample_counts: list[tuple[str, int, int]] = []  # (sample, n_sites, max_round)
    total_sites = 0
    for sf in stat_files:
        sample = sf.parent.parent.name
        rounds = [r for _, r in parse_stats_file(sf)]
        if not rounds:
            continue
        for r in rounds:
            final_round_hist[r] += 1
        total_sites += len(rounds)
        per_sample_counts.append((sample, len(rounds), max(rounds)))

    print()
    print(f"# per-site final assembly_rounds  (source: {len(stat_files)} s05_stats.txt)")
    print(f"# total insertion sites: {total_sites}")
    total = sum(final_round_hist.values()) or 1
    print("final_round\tn_sites\tfraction\tcumulative_fraction")
    cum_sites = 0
    for rnd in sorted(final_round_hist):
        cum_sites += final_round_hist[rnd]
        print(
            f"{rnd}\t{final_round_hist[rnd]}\t"
            f"{final_round_hist[rnd]/total:.4f}\t{cum_sites/total:.4f}"
        )

    print()
    print("# per-sample site counts (sample\tn_sites\tmax_final_round)")
    for sample, n, mr in per_sample_counts:
        print(f"{sample}\t{n}\t{mr}")


if __name__ == "__main__":
    main()
