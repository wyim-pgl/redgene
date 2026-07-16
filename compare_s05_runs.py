#!/usr/bin/env python3
"""Diff s05 verdicts between the April 2026 baselines and a re-run root.

Reads `s05_stats.txt`, never the per-site report files: s05 does not clean its
output directory, so a directory re-used across runs accumulates orphan
`insertion_<site>_report.txt` files that make a glob-and-count overstate the
site count. `write_stats` truncates `s05_stats.txt` every run, so it alone
describes exactly one run. (rice_G281's April directory carries 47 such orphans
across three runs — the bug that produced the "68 sites" figure this branch's
docs used to quote.)

Usage:
    python compare_s05_runs.py [NEW_ROOT]        # default results/hardened_20260710
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent
BASELINE_ROOT = REPO / "results"
DEFAULT_NEW_ROOT = REPO / "results" / "hardened_20260710"

SAMPLES = [
    "rice_G281",
    "tomato_WT",
    "tomato_Cas9_A2_1",
    "tomato_Cas9_A2_2",
    "tomato_Cas9_A2_3",
    "corn_ND207",
    "cucumber_line212",
    "cucumber_line224",
    "cucumber_line225",
    "soybean_AtYUCCA6",
    "soybean_UGT72E3",
]

# Documented true insertions, in the coordinate system of the reference in db/.
# Sources: benchmark.md "Test Datasets". Only bp-resolution entries appear here;
# the cucumber lines and soybean_AtYUCCA6 are documented at gene/chromosome
# granularity only, and soybean_AtYUCCA6's Glyma.* gene IDs are absent from the
# RefSeq GFF, so neither can be checked by coordinate.
#
# soybean_UGT72E3's coordinate was verified against db/Gmax_v4.0.gff3.gz:
# LOC102669442 (= Glyma.18g226800) spans NC_038254.2:51,789,501-51,906,722, and
# 51,882,860 falls inside it — so benchmark.md's coordinates are on the same
# assembly this pipeline maps against, not the source paper's.
#
# These are NOT pass/fail assertions. Three of the four are already absent or
# non-CANDIDATE in the April baselines; this script reports ground-truth status
# on both sides and fails only on a *regression* the hardening introduced.
GROUND_TRUTH: dict[str, list[tuple[str, int, str]]] = {
    "rice_G281": [("Chr3", 16_439_674, "2 copies, head-to-head")],
    "tomato_Cas9_A2_3": [("SLM_r2.0ch08", 65_107_378, "T-DNA junction per benchmark.md")],
    "corn_ND207": [("NC_050098.1", 181_367_276, "chr3, Bt construct")],
    "soybean_UGT72E3": [("NC_038254.2", 51_882_860, "Gm18, inside LOC102669442 CDS")],
}
# A discovered site counts as the ground-truth site if it lands within this many
# bases of the documented coordinate. Junctions are called at the clip boundary,
# so a true site can sit a few hundred bp off a paper's reported position.
GT_TOLERANCE_BP = 1_000

VERDICTS = ("CANDIDATE", "FALSE_POSITIVE", "UNKNOWN")
_VERDICT_KEY = re.compile(r"^insertion_(.+)_verdict$")


def nearest(verdicts: dict[str, str], chrom: str, pos: int) -> tuple[int, str, str] | None:
    """(distance, site_id, verdict) of the closest site on `chrom`, or None."""
    best = None
    for site_id, verdict in verdicts.items():
        c, _, p = site_id.rpartition("_")  # chrom names contain '_' and '.'
        if c != chrom or not p.isdigit():
            continue
        d = abs(int(p) - pos)
        if best is None or d < best[0]:
            best = (d, site_id, verdict)
    return best


def gt_status(verdicts: dict[str, str], chrom: str, pos: int) -> str:
    hit = nearest(verdicts, chrom, pos)
    if hit is None or hit[0] > GT_TOLERANCE_BP:
        return "not discovered"
    return hit[2]


def read_verdicts(stats_path: Path) -> dict[str, str] | None:
    """site_id -> verdict, or None if the run never produced stats."""
    if not stats_path.is_file():
        return None
    out: dict[str, str] = {}
    for line in stats_path.read_text().splitlines():
        key, _, value = line.partition("\t")
        m = _VERDICT_KEY.match(key)
        if m:
            out[m.group(1)] = value
    return out


def tally(verdicts: dict[str, str]) -> dict[str, int]:
    counts = {v: 0 for v in VERDICTS}
    for v in verdicts.values():
        counts[v] = counts.get(v, 0) + 1
    return counts


def main() -> int:
    new_root = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_NEW_ROOT

    # "reported" counts sites carrying a verdict, i.e. CAND + FP + UNK. This is
    # NOT the `insertion_sites` stat: that key has meant different things across
    # versions (Phase 1 junctions in the April 2026-04-12/13 runs, transgene-
    # positive sites later), so tomato_WT's April stats claim 6,875 sites while
    # recording zero verdicts. Verdict counts are comparable across versions.
    header = f"{'SAMPLE':<20} {'reported':>11}  {'CAND':>7}  {'FP':>7}  {'UNK':>7}"
    print(header)
    print("-" * len(header))

    details: list[str] = []
    truth: list[str] = []
    regressions: list[str] = []

    for sample in SAMPLES:
        old = read_verdicts(BASELINE_ROOT / sample / "s05_insert_assembly" / "s05_stats.txt")
        new = read_verdicts(new_root / sample / "s05_insert_assembly" / "s05_stats.txt")

        if new is None:
            print(f"{sample:<20} {'(re-run pending)':>11}")
            continue
        if old is None:
            print(f"{sample:<20} {'(no baseline)':>11}")
            continue

        o, n = tally(old), tally(new)

        def cell(a: int, b: int) -> str:
            return f"{a}→{b}" if a != b else str(a)

        print(
            f"{sample:<20} {cell(len(old), len(new)):>11}  "
            f"{cell(o['CANDIDATE'], n['CANDIDATE']):>7}  "
            f"{cell(o['FALSE_POSITIVE'], n['FALSE_POSITIVE']):>7}  "
            f"{cell(o['UNKNOWN'], n['UNKNOWN']):>7}"
        )

        changed = sorted(s for s in old.keys() & new.keys() if old[s] != new[s])
        gone = sorted(old.keys() - new.keys())
        added = sorted(new.keys() - old.keys())

        if changed or gone or added:
            details.append(f"\n{sample}:")
            for s in changed:
                details.append(f"    {s:<26} {old[s]} → {new[s]}")
            for s in gone:
                details.append(f"    {s:<26} {old[s]} → (no longer discovered)")
            for s in added:
                details.append(f"    {s:<26} (newly discovered) → {new[s]}")

        # Ground truth: report the status on both sides. Only a CANDIDATE that
        # stopped being CANDIDATE is a regression this branch caused; a site the
        # baseline already missed is a pre-existing recall gap, not a regression.
        for chrom, pos, note in GROUND_TRUTH.get(sample, []):
            was, now = gt_status(old, chrom, pos), gt_status(new, chrom, pos)
            arrow = f"{was} → {now}" if was != now else was
            truth.append(f"    {sample:<20} {chrom}:{pos:,}  {arrow}   ({note})")
            if was == "CANDIDATE" and now != "CANDIDATE":
                regressions.append(
                    f"    {sample} {chrom}:{pos:,} was CANDIDATE, is now {now}"
                )

    if details:
        print("\nPer-site changes:")
        print("\n".join(details))

    if truth:
        print("\nGround truth (documented in benchmark.md, checked within "
              f"±{GT_TOLERANCE_BP:,} bp):")
        print("\n".join(truth))

    print("\nVerdict:")
    if regressions:
        print("  REGRESSION — the hardening lost a true insertion:")
        print("\n".join(regressions))
        return 1
    print("  no ground-truth regression")
    return 0


if __name__ == "__main__":
    sys.exit(main())
