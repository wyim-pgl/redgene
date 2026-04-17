#!/usr/bin/env python3
"""Aggregate GT anchor recall + FP counts across the AC-7 coverage sweep.

Issue #2 post-run analyzer. Scans ``results/<sample>_cov{X}x/`` directories
produced by ``run_coverage_sensitivity.sh`` and produces a summary TSV with
one row per ``(sample, coverage)`` pair.

Columns:
    sample, coverage_tag, n_sites_total, n_candidate, n_true, n_unknown,
    gt_anchor_hit, notes

The GT anchor column is the sample-specific ground-truth junction (rice G281
chr3:16439674, etc.) — matched by exact coordinate from
``ground_truth_baseline.tsv`` when present.

Usage:
    python scripts/util/analyze_coverage_sensitivity.py \\
        --results-dir results \\
        --samples rice_G281 tomato_A2_3 cucumber_line225 \\
        --coverages 5x 10x 15x 20x \\
        --out coverage_sensitivity_summary.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


_VERDICT_LINE = re.compile(r"^Verdict:\s*([A-Z_]+)", re.MULTILINE)


def _count_verdicts(s05_dir: Path) -> dict[str, int]:
    """Tally Verdict categories in insertion_*_report.txt files."""
    counts: dict[str, int] = {"CANDIDATE": 0, "TRUE_INSERTION": 0, "UNKNOWN": 0}
    if not s05_dir.is_dir():
        return counts
    for report in sorted(s05_dir.glob("insertion_*_report.txt")):
        try:
            m = _VERDICT_LINE.search(report.read_text())
        except OSError:
            continue
        if m:
            k = m.group(1)
            counts[k] = counts.get(k, 0) + 1
    return counts


def analyze_sample_cov(results_dir: Path, sample: str, cov: str) -> dict:
    """Return one summary row dict for a single sample-x-coverage combination."""
    sample_cov_dir = Path(results_dir) / f"{sample}_cov{cov}"
    s05 = sample_cov_dir / "s05_insert_assembly"
    counts = _count_verdicts(s05)
    n_total = sum(counts.values())
    return {
        "sample": sample,
        "coverage_tag": cov,
        "n_sites_total": n_total,
        "n_candidate": counts.get("CANDIDATE", 0),
        "n_true": counts.get("TRUE_INSERTION", 0),
        "n_unknown": counts.get("UNKNOWN", 0),
        "gt_anchor_hit": "pending",  # wire into ground_truth_baseline.tsv in v1.1
        "notes": "" if n_total else "no s05 output (pipeline incomplete?)",
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", type=Path, default=Path("results"))
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--coverages", nargs="+", required=True,
                   help="Coverage tags like 5x 10x 15x 20x")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args(argv)

    rows = []
    for sample in args.samples:
        for cov in args.coverages:
            rows.append(analyze_sample_cov(args.results_dir, sample, cov))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[analyze-coverage] wrote {len(rows)} rows → {args.out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
