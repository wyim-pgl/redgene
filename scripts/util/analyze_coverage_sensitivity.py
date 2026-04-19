#!/usr/bin/env python3
"""Aggregate GT anchor recall + FP counts across the AC-7 coverage sweep.

Issue #2 post-run analyzer. Scans ``results/<sample>_cov{X}x/`` directories
produced by ``run_coverage_sensitivity.sh`` and produces a summary TSV with
one row per ``(sample, coverage)`` pair.

Columns:
    sample, coverage_tag, n_sites_total, n_candidate, n_true, n_unknown,
    gt_anchor_hit, notes

``gt_anchor_hit`` is resolved against ``ground_truth_baseline.tsv`` when
supplied via ``--ground-truth``:
    HIT:<verdict>   matching insertion_<chrom>_<pos>_report.txt found
    MISS            GT exists for sample, no matching report file
    no_gt           GT row for sample not present in baseline TSV

Usage:
    python scripts/util/analyze_coverage_sensitivity.py \\
        --results-dir results \\
        --samples rice_G281 tomato_Cas9_A2_3 cucumber_line225 \\
        --coverages 5x 10x 15x 20x \\
        --ground-truth ground_truth_baseline.tsv \\
        --out coverage_sensitivity_summary.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


_VERDICT_LINE = re.compile(r"^Verdict:\s*([A-Z_]+)", re.MULTILINE)
_REPORT_STEM = re.compile(r"^insertion_(.*)_(\d+)_report\.txt$")


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


def _load_ground_truth(path: Path | None) -> dict[str, list[tuple[str, int]]]:
    """Return {sample: [(chrom, pos), ...]} from baseline TSV."""
    if path is None or not path.is_file():
        return {}
    gt: dict[str, list[tuple[str, int]]] = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = row.get("sample", "").strip()
            chrom = row.get("chrom", "").strip()
            pos_raw = row.get("position", "").strip()
            if not sample or not chrom or not pos_raw:
                continue
            try:
                pos = int(pos_raw)
            except ValueError:
                continue
            gt.setdefault(sample, []).append((chrom, pos))
    return gt


def _chrom_equal(a: str, b: str) -> bool:
    """Match chrom names ignoring NCBI version suffix (LKUO03001451 == LKUO03001451.1)."""
    strip = lambda s: re.sub(r"\.\d+$", "", s)
    return strip(a) == strip(b)


def _match_gt_report(
    s05_dir: Path,
    gt_entries: list[tuple[str, int]],
    pos_tolerance: int = 0,
) -> tuple[str, Path | None]:
    """Return (status, report_path) for the first GT hit.

    status is 'HIT:<verdict>' or 'MISS'. pos_tolerance allows fuzzy matching
    on junction position if needed (0 = exact).
    """
    if not gt_entries or not s05_dir.is_dir():
        return ("MISS", None)
    candidates: list[tuple[str, int, Path]] = []
    for report in s05_dir.glob("insertion_*_report.txt"):
        m = _REPORT_STEM.match(report.name)
        if not m:
            continue
        candidates.append((m.group(1), int(m.group(2)), report))
    for gt_chrom, gt_pos in gt_entries:
        for rep_chrom, rep_pos, path in candidates:
            if _chrom_equal(rep_chrom, gt_chrom) and abs(rep_pos - gt_pos) <= pos_tolerance:
                verdict = "UNKNOWN"
                try:
                    m = _VERDICT_LINE.search(path.read_text())
                    if m:
                        verdict = m.group(1)
                except OSError:
                    pass
                return (f"HIT:{verdict}", path)
    return ("MISS", None)


def analyze_sample_cov(
    results_dir: Path,
    sample: str,
    cov: str,
    ground_truth: dict[str, list[tuple[str, int]]] | None = None,
) -> dict:
    """Return one summary row dict for a single sample-x-coverage combination."""
    sample_cov_dir = Path(results_dir) / f"{sample}_cov{cov}"
    s05 = sample_cov_dir / "s05_insert_assembly"
    counts = _count_verdicts(s05)
    n_total = sum(counts.values())

    gt_entries = (ground_truth or {}).get(sample, [])
    if ground_truth is None:
        gt_status = "pending"
    elif not gt_entries:
        gt_status = "no_gt"
    else:
        gt_status, _ = _match_gt_report(s05, gt_entries)

    if n_total > 0:
        note = ""
    elif s05.is_dir():
        note = "pipeline complete, 0 insertion reports"
    else:
        note = "no s05 output (pipeline incomplete?)"

    return {
        "sample": sample,
        "coverage_tag": cov,
        "n_sites_total": n_total,
        "n_candidate": counts.get("CANDIDATE", 0),
        "n_true": counts.get("TRUE_INSERTION", 0),
        "n_unknown": counts.get("UNKNOWN", 0),
        "gt_anchor_hit": gt_status,
        "notes": note,
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", type=Path, default=Path("results"))
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--coverages", nargs="+", required=True,
                   help="Coverage tags like 5x 10x 15x 20x")
    p.add_argument("--ground-truth", type=Path, default=None,
                   help="ground_truth_baseline.tsv for gt_anchor_hit resolution")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args(argv)

    gt = _load_ground_truth(args.ground_truth) if args.ground_truth else None

    rows = []
    for sample in args.samples:
        for cov in args.coverages:
            rows.append(analyze_sample_cov(args.results_dir, sample, cov, gt))

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
