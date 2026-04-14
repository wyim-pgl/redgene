#!/usr/bin/env python3
"""Retroactively apply the UNKNOWN -> CANDIDATE_LOW_CONF heuristic to an
s05 results directory without re-running the pipeline.

Mirrors the logic now in scripts/s05_insert_assembly.py (commit f78912b):
  if verdict == UNKNOWN and n_borders >= 4 and largest_gap >= 1000:
      verdict = CANDIDATE_LOW_CONF

Updates both the per-site insertion_*_report.txt header line and the
aggregate s05_stats.txt verdict field. Original files are backed up with
a `.preheur` suffix on first run.

Usage:
    python scripts/postprocess_unknown_reclass.py results/soybean_AtYUCCA6/s05_insert_assembly
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

UNKNOWN_MIN_BORDERS = 4
UNKNOWN_MIN_FOREIGN_GAP = 1000

VERDICT_RE = re.compile(r"^Verdict: ([A-Z_]+)(?: — (.*))?$", re.MULTILINE)
BORDERS_RE = re.compile(r"^T-DNA borders found: (\d+)", re.MULTILINE)
GAP_RE = re.compile(r"largest non-host gap: ([\d,]+)bp")


def _backup(path: Path) -> None:
    bak = path.with_suffix(path.suffix + ".preheur")
    if not bak.exists():
        shutil.copy2(path, bak)


def reclass_report(report_path: Path) -> tuple[str, str] | None:
    """Return (old_verdict, new_verdict) if the report was rewritten."""
    text = report_path.read_text()
    m = VERDICT_RE.search(text)
    if not m:
        return None
    old_verdict = m.group(1)
    if old_verdict != "UNKNOWN":
        return None

    bm = BORDERS_RE.search(text)
    gm = GAP_RE.search(text)
    n_borders = int(bm.group(1)) if bm else 0
    largest_gap = int(gm.group(1).replace(",", "")) if gm else 0

    if n_borders < UNKNOWN_MIN_BORDERS or largest_gap < UNKNOWN_MIN_FOREIGN_GAP:
        return None

    new_reason = (
        f"no element annotations but {n_borders} T-DNA border motifs + "
        f"{largest_gap:,}bp non-host gap — likely sample-specific transgene "
        f"not in element_db"
    )
    new_line = f"Verdict: CANDIDATE_LOW_CONF — {new_reason}"
    new_text = VERDICT_RE.sub(new_line, text, count=1)

    _backup(report_path)
    report_path.write_text(new_text)
    return old_verdict, "CANDIDATE_LOW_CONF"


def update_stats(stats_path: Path, promoted_sites: set[str]) -> int:
    """Rewrite insertion_<site_id>_verdict rows in s05_stats.txt. Returns count."""
    if not stats_path.exists():
        return 0
    lines = stats_path.read_text().splitlines()
    changed = 0
    for i, line in enumerate(lines):
        if not line.startswith("insertion_") or "_verdict\t" not in line:
            continue
        key, _, value = line.partition("\t")
        site_id = key[: -len("_verdict")]
        if value == "UNKNOWN" and site_id in promoted_sites:
            lines[i] = f"{key}\tCANDIDATE_LOW_CONF"
            changed += 1
    if changed:
        _backup(stats_path)
        stats_path.write_text("\n".join(lines) + "\n")
    return changed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("s05_dir", type=Path,
                    help="Path to results/<sample>/s05_insert_assembly/")
    args = ap.parse_args()

    s05_dir: Path = args.s05_dir
    if not s05_dir.is_dir():
        print(f"ERROR: not a directory: {s05_dir}", file=sys.stderr)
        return 2

    reports = sorted(s05_dir.glob("insertion_*_report.txt"))
    if not reports:
        print(f"No insertion_*_report.txt under {s05_dir}", file=sys.stderr)
        return 1

    promoted: list[tuple[str, str]] = []
    for rpt in reports:
        result = reclass_report(rpt)
        if result is None:
            continue
        site_id = rpt.name[: -len("_report.txt")]
        promoted.append((site_id, rpt.name))
        print(f"  PROMOTED: {site_id}  (UNKNOWN -> CANDIDATE_LOW_CONF)")

    stats_changed = update_stats(
        s05_dir / "s05_stats.txt",
        {sid for sid, _ in promoted},
    )

    print()
    print(f"Reports scanned:       {len(reports)}")
    print(f"Reports promoted:      {len(promoted)}")
    print(f"s05_stats.txt updates: {stats_changed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
