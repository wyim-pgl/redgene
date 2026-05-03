#!/usr/bin/env python3
"""Assign KCGP tentative_v0 IDs to sites in a sample's s05 output.

Reads results/{sample}/s05_insert_assembly/insertion_*_report.txt, sorts
CANDIDATE sites first then by host_chr (natural) and pos_5p, and writes
results/{sample}/kcgp_mapping.tsv with one row per discovered report.
Only CANDIDATE rows receive a 4-digit ordinal; FALSE_POSITIVE / UNKNOWN
rows carry kcgp_id = "-".
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from datetime import datetime
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kcgp_id import HOST_CODES, build_kcgp_id, _normalize_host_path  # noqa: E402

REPORT_RE = re.compile(r"^insertion_(.+)_(\d+)_report\.txt$")
VERDICT_PRIORITY = {
    "CANDIDATE": 0,
    "FALSE_POSITIVE": 1,
    "UNKNOWN": 2,
}
_NATURAL_RE = re.compile(r"(\d+|\D+)")


def _natural_key(s: str):
    parts = []
    for tok in _NATURAL_RE.findall(s):
        parts.append((0, int(tok)) if tok.isdigit() else (1, tok))
    return parts


def log(msg: str) -> None:
    print(f"[assign_kcgp_ids] {msg}", file=sys.stderr, flush=True)


def parse_report(path: Path) -> tuple[str, int, str] | None:
    m = REPORT_RE.match(path.name)
    if not m:
        log(f"skip {path.name}: filename does not match insertion_*_report.txt")
        return None
    host_chr, pos_str = m.group(1), m.group(2)
    pos = int(pos_str)
    verdict = "UNKNOWN"
    try:
        for line in path.read_text().splitlines():
            if line.startswith("Verdict:"):
                token = line[len("Verdict:") :].strip().split()[0]
                verdict = token.rstrip(".,—-")
                break
        else:
            log(f"warn {path.name}: no Verdict line, treating as UNKNOWN")
    except OSError as exc:
        log(f"warn {path.name}: read failed ({exc}), treating as UNKNOWN")
    return host_chr, pos, verdict


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--sample-dir", required=True, type=Path)
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--year", type=int, default=datetime.now().year)
    args = p.parse_args(argv)

    sample = args.sample_dir.name
    try:
        cfg = yaml.safe_load(args.config.read_text()) or {}
    except (OSError, yaml.YAMLError) as exc:
        log(f"FATAL: cannot read config: {exc}")
        return 1

    samples = cfg.get("samples", {}) or {}
    if sample not in samples:
        log(f"FATAL: sample {sample!r} not in {args.config}")
        return 1

    entry = samples[sample] or {}
    event = entry.get("event")
    if not event:
        log(
            f"FATAL: sample {sample!r} has no 'event:' field in {args.config}; "
            f"add 'event:' under samples.{sample} (e.g., event: G281)"
        )
        return 1

    host_ref = entry.get("host_reference")
    if not host_ref:
        log(f"FATAL: sample {sample!r} has no 'host_reference:' field")
        return 1

    norm = _normalize_host_path(host_ref)
    if norm not in HOST_CODES:
        log(
            f"FATAL: host '{norm}' (from {host_ref}) not in HOST_CODES; "
            f"known: {sorted(HOST_CODES)}"
        )
        return 1

    s05_dir = args.sample_dir / "s05_insert_assembly"
    reports: list[tuple[str, int, str, Path]] = []
    if s05_dir.exists():
        for report_path in sorted(s05_dir.glob("insertion_*_report.txt")):
            parsed = parse_report(report_path)
            if parsed is None:
                continue
            host_chr, pos, verdict = parsed
            reports.append((host_chr, pos, verdict, report_path))
    else:
        log(f"warn: {s05_dir} does not exist; emitting header-only mapping")

    reports.sort(key=lambda r: (
        VERDICT_PRIORITY.get(r[2], 99),
        _natural_key(r[0]),
        r[1],
    ))

    n_cand = sum(1 for r in reports if r[2] == "CANDIDATE")
    if n_cand > 9999:
        log(f"FATAL: {n_cand} CANDIDATE sites exceeds 4-digit ordinal range")
        return 1

    out_path = args.sample_dir / "kcgp_mapping.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ordinal = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "kcgp_id", "internal_site", "verdict", "host_chr", "pos_5p", "spec_version",
        ])
        for host_chr, pos, verdict, report_path in reports:
            if verdict == "CANDIDATE":
                ordinal += 1
                kid = build_kcgp_id(host_ref, event, args.year, ordinal)
            else:
                kid = "-"
            internal = f"{sample}_{host_chr}_{pos}"
            writer.writerow([kid, internal, verdict, host_chr, pos, "tentative_v0"])

    log(f"wrote {out_path} ({len(reports)} rows, {n_cand} CANDIDATE)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
