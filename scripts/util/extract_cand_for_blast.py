#!/usr/bin/env python3
"""Extract CANDIDATE insertion contigs for remote BLAST verification.

Issue #1 [AC-2] helper: scans ``results/<sample>/s05_insert_assembly/`` for
``insertion_*_report.txt`` files whose ``Verdict:`` line is ``CANDIDATE`` and
merges the matching ``insertion_*_insert.fasta`` records into a single
per-sample query FASTA suitable for ``blastn -remote -db nt`` classification.

This helper *only prepares* the query FASTA — it never calls ``blastn`` and
never submits to SLURM. ``scripts/util/prepare_remote_blast.sh`` emits the
accompanying SLURM job template. Submission is a manual, operator-driven step.

Usage:
    python scripts/util/extract_cand_for_blast.py \\
        --sample-dir results/cucumber_line225 \\
        --sample-name cucumber_line225 \\
        --out cucumber_line225_cand_inserts.fa
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable, List


_VERDICT_RE = re.compile(r"^Verdict:\s*CANDIDATE\b", re.MULTILINE)
_REPORT_RE = re.compile(r"^insertion_(?P<stem>.+)_report\.txt$")


def find_candidate_reports(s05_dir: Path) -> List[Path]:
    """Return report.txt paths whose Verdict is CANDIDATE.

    Reports with ``Verdict: CANDIDATE*`` (e.g. ``CANDIDATE — ...``) also match.
    """
    s05_dir = Path(s05_dir)
    if not s05_dir.is_dir():
        return []
    hits: List[Path] = []
    for report in sorted(s05_dir.glob("insertion_*_report.txt")):
        try:
            text = report.read_text()
        except OSError:
            continue
        if _VERDICT_RE.search(text):
            hits.append(report)
    return hits


def report_to_insert_fasta(report: Path) -> Path:
    """Map ``insertion_<stem>_report.txt`` → ``insertion_<stem>_insert.fasta``."""
    m = _REPORT_RE.match(report.name)
    if not m:
        raise ValueError(f"Not an insertion report filename: {report.name}")
    stem = m.group("stem")
    return report.parent / f"insertion_{stem}_insert.fasta"


def _iter_fasta_records(fasta: Path) -> Iterable[tuple[str, str]]:
    header: str | None = None
    buf: list[str] = []
    with open(fasta) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf)
                header = line[1:]
                buf = []
            elif line:
                buf.append(line)
    if header is not None:
        yield header, "".join(buf)


def merge_cand_inserts(
    sample_dir: Path,
    out_fasta: Path,
    sample_name: str,
) -> int:
    """Concatenate CAND insert FASTAs under ``sample_dir`` into ``out_fasta``.

    Each record header is re-tagged with ``{sample_name}|`` so BLAST hit TSVs
    remain traceable back to the source sample. Returns the number of CAND
    sites merged. If zero, writes an empty file (still a valid no-op query).
    """
    sample_dir = Path(sample_dir)
    out_fasta = Path(out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    s05 = sample_dir / "s05_insert_assembly"
    reports = find_candidate_reports(s05)

    if not reports:
        out_fasta.write_text("")
        return 0

    count = 0
    with open(out_fasta, "w") as outfh:
        for report in reports:
            fasta = report_to_insert_fasta(report)
            if not fasta.exists():
                print(
                    f"[extract_cand] WARN: missing insert FASTA for {report.name}",
                    file=sys.stderr,
                )
                continue
            stem = _REPORT_RE.match(report.name).group("stem")  # type: ignore[union-attr]
            for idx, (hdr, seq) in enumerate(_iter_fasta_records(fasta)):
                new_hdr = f"{sample_name}|{stem}|{idx}|{hdr}"
                outfh.write(f">{new_hdr}\n{seq}\n")
            count += 1
    return count


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample-dir", type=Path, required=True,
                   help="results/<sample>/ directory (contains s05_insert_assembly/)")
    p.add_argument("--sample-name", type=str, required=True,
                   help="Sample key, used as record-header prefix for BLAST traceability")
    p.add_argument("--out", type=Path, required=True,
                   help="Output merged FASTA (overwritten if present)")
    args = p.parse_args(argv)

    n = merge_cand_inserts(args.sample_dir, args.out, args.sample_name)
    print(f"[extract_cand] {n} CANDIDATE sites merged → {args.out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
