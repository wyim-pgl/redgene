#!/usr/bin/env python3
"""v1.1 PDF insertion report scaffold (Issue #6).

Generates a 7-section PDF summarizing a single sample's RedGene pipeline
output. This is a *scaffold*: several sections (junction diagram, CRISPR
editing) are placeholders that will be fleshed out in v1.1.

Sections:
    (1) Cover page        — sample name, date, audit_header summary
    (2) Sample summary    — site count, verdict distribution
    (3) Insertion sites   — tabular view parsed from insertion_*_report.txt
    (4) Junction diagram  — placeholder ("diagram pending v1.1")
    (5) Copy number       — optional, shown if results/<sample>/s07_copynumber/ exists
    (6) CRISPR editing    — placeholder ("pending v1.1")
    (7) Appendix          — pipeline commit, DB md5s, software versions

Dependency: matplotlib (already in the `redgene` micromamba env). NO reportlab.

Usage:
    python scripts/reports/insertion_pdf.py \\
        --sample rice_G281 --sample-dir results/rice_G281 \\
        --out results/rice_G281/rice_G281_insertion_report.pdf
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


_SITE_RE = re.compile(r"^Insertion site:\s*(?P<site>.+)$", re.MULTILINE)
_VERDICT_RE = re.compile(r"^Verdict:\s*(?P<v>[A-Z_]+)", re.MULTILINE)
_LEN_RE = re.compile(r"^Insert length:\s*(?P<len>[\d,]+)\s*bp", re.MULTILINE)


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------


def load_audit_header(sample_dir: Path) -> dict[str, Any]:
    """Load ``audit_header.json`` from the sample results dir, or return {}."""
    hdr_path = Path(sample_dir) / "audit_header.json"
    if not hdr_path.exists():
        return {}
    try:
        return json.loads(hdr_path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def parse_insertion_reports(s05_dir: Path) -> list[dict[str, Any]]:
    """Parse ``insertion_*_report.txt`` into structured rows."""
    s05_dir = Path(s05_dir)
    rows: list[dict[str, Any]] = []
    if not s05_dir.is_dir():
        return rows
    for report in sorted(s05_dir.glob("insertion_*_report.txt")):
        try:
            text = report.read_text()
        except OSError:
            continue
        site = (m := _SITE_RE.search(text)) and m.group("site").strip() or report.stem
        verdict = (m := _VERDICT_RE.search(text)) and m.group("v") or "UNSPECIFIED"
        length = (m := _LEN_RE.search(text)) and m.group("len") or "—"
        rows.append({
            "site": site,
            "verdict": verdict,
            "insert_length": length,
            "source": report.name,
        })
    return rows


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate verdict counts for the summary page."""
    counter = Counter(r["verdict"] for r in rows)
    return {
        "n_sites": len(rows),
        "by_verdict": dict(counter),
    }


# ---------------------------------------------------------------------------
# PDF rendering
# ---------------------------------------------------------------------------


def _page_text(pdf: PdfPages, title: str, body_lines: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis("off")
    ax.text(0.05, 0.95, title, fontsize=18, fontweight="bold",
            transform=ax.transAxes, verticalalignment="top")
    y = 0.88
    for line in body_lines:
        ax.text(0.05, y, line, fontsize=11, family="monospace",
                transform=ax.transAxes, verticalalignment="top")
        y -= 0.025
        if y < 0.05:
            break
    pdf.savefig(fig)
    plt.close(fig)


def _page_cover(pdf: PdfPages, sample_name: str, audit: dict[str, Any]) -> None:
    lines = [
        f"Sample:        {sample_name}",
        f"Generated:     {audit.get('generated_at', '—')}",
        f"Commit:        {audit.get('pipeline_commit', '—')}  (dirty={audit.get('pipeline_dirty', '—')})",
        "",
        "Input SHA-256:",
    ]
    for k, v in (audit.get("input_sha256") or {}).items():
        lines.append(f"  {k}: {v}")
    lines += ["", "DB manifest:"]
    for db in audit.get("db_manifest", []) or []:
        lines.append(
            f"  - {db.get('name')} md5={db.get('md5')} built={db.get('build_date')} seqs={db.get('seq_count')}"
        )
    if not audit:
        lines.append("  (audit_header.json not present — placeholder page)")
    _page_text(pdf, "RedGene Insertion Report", lines)


def _page_summary(pdf: PdfPages, summary: dict[str, Any]) -> None:
    lines = [
        f"Total sites parsed: {summary['n_sites']}",
        "",
        "Verdict distribution:",
    ]
    for v, n in sorted(summary.get("by_verdict", {}).items()):
        lines.append(f"  {v:<20s}  {n}")
    _page_text(pdf, "Sample Summary", lines)


def _page_sites_table(pdf: PdfPages, rows: list[dict[str, Any]]) -> None:
    body = [f"{'site':<44s}  {'verdict':<18s}  {'len_bp':>8s}  source"]
    body.append("-" * 100)
    for r in rows[:40]:  # first 40 fit one page; v1.1 adds pagination
        body.append(
            f"{r['site']:<44s}  {r['verdict']:<18s}  {r['insert_length']:>8s}  {r['source']}"
        )
    if not rows:
        body.append("(no insertion reports found)")
    _page_text(pdf, "Insertion Sites", body)


def _page_placeholder(pdf: PdfPages, title: str, note: str) -> None:
    _page_text(pdf, title, [note, "", "(scaffold — to be wired in v1.1)"])


def _page_copy_number(pdf: PdfPages, sample_dir: Path) -> None:
    s07 = Path(sample_dir) / "s07_copynumber"
    lines: list[str] = []
    if s07.is_dir():
        tsvs = sorted(s07.glob("*.tsv"))
        lines.append(f"s07_copynumber/ TSVs: {len(tsvs)}")
        for t in tsvs[:5]:
            lines.append(f"  - {t.name}")
    else:
        lines.append("(s07_copynumber/ not present — placeholder)")
    _page_text(pdf, "Copy Number", lines)


def _page_appendix(pdf: PdfPages, audit: dict[str, Any]) -> None:
    lines = [
        f"Pipeline commit: {audit.get('pipeline_commit', '—')}",
        f"Dirty: {audit.get('pipeline_dirty', '—')}",
        "",
        "Software versions:",
    ]
    for k, v in (audit.get("software_versions") or {}).items():
        lines.append(f"  {k}: {v}")
    lines += ["", "DB md5:"]
    for db in audit.get("db_manifest", []) or []:
        lines.append(f"  {db.get('name')}: {db.get('md5')}")
    if not audit:
        lines.append("(audit_header.json absent — appendix placeholder)")
    _page_text(pdf, "Appendix", lines)


def generate_pdf(
    *,
    sample_dir: Path,
    sample_name: str,
    out_pdf: Path,
) -> int:
    """Render the 7-section PDF. Returns the number of pages written."""
    sample_dir = Path(sample_dir)
    out_pdf = Path(out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    audit = load_audit_header(sample_dir)
    rows = parse_insertion_reports(sample_dir / "s05_insert_assembly")
    summary = summarize(rows)

    pages = 0
    with PdfPages(out_pdf) as pdf:
        _page_cover(pdf, sample_name, audit); pages += 1              # (1)
        _page_summary(pdf, summary); pages += 1                        # (2)
        _page_sites_table(pdf, rows); pages += 1                       # (3)
        _page_placeholder(pdf, "Junction Diagram",
                          "diagram pending v1.1"); pages += 1          # (4)
        _page_copy_number(pdf, sample_dir); pages += 1                 # (5)
        _page_placeholder(pdf, "CRISPR Editing",
                          "editing panel pending v1.1"); pages += 1    # (6)
        _page_appendix(pdf, audit); pages += 1                         # (7)
    return pages


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True, help="Sample name (for page titles)")
    p.add_argument("--sample-dir", type=Path, required=True,
                   help="results/<sample>/ directory")
    p.add_argument("--out", type=Path, required=True,
                   help="Output PDF path (parent directory auto-created)")
    args = p.parse_args(argv)

    n = generate_pdf(
        sample_dir=args.sample_dir,
        sample_name=args.sample,
        out_pdf=args.out,
    )
    print(f"[insertion-pdf] wrote {n} pages → {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
