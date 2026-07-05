#!/usr/bin/env python3
"""v1.1 PDF insertion report (Issue #6 / R-5).

Generates a multi-section PDF summarizing a single sample's RedGene pipeline
output.  Uses matplotlib's PdfPages backend (``reportlab`` is not available in
the ``redgene`` micromamba env; we intentionally stay pure-matplotlib).

Sections rendered (one page each unless noted):
    (1) Cover page          — sample name, date, audit_header summary
    (2) Sample information  — sample_id, host, construct, date
    (3) Pipeline version    — git commit, dirty flag, software versions
    (4) Sample summary      — site count, verdict distribution
    (5) Insertion sites     — paginated table (parsed from
        insertion_*_report.txt) with verdict / length / structure / elements
    (6) Junction diagrams   — per-candidate soft-clip anchor markers (one
        page; sites beyond first N are summarised)
    (7) Copy number         — matplotlib bar chart if s07_copynumber/
        copynumber.tsv exists, text block otherwise
    (8) CRISPR editing      — placeholder ("pending v1.1")
    (9) Appendix — audit    — pipeline commit, DB md5s, software versions
    (10) Appendix — CoC     — chain_of_custody.jsonl digest (only rendered
        when the file is present on disk)

The module exposes small pure-Python helpers (``load_audit_header``,
``parse_insertion_reports``, ``summarize``, ``load_copynumber``,
``load_coc_log``, ``extract_site_anchors``, ``build_sample_info``,
``build_pipeline_version``) so unit tests can assert each data block
independently without rendering a PDF.

Usage:
    python scripts/reports/insertion_pdf.py \\
        --sample rice_G281 --sample-dir results/rice_G281 \\
        --out results/rice_G281/rice_G281_insertion_report.pdf

Optional:
    --config config.yaml       (to look up host / construct references)
    --host HOST_LABEL          (overrides config)
    --construct CONSTRUCT_LBL  (overrides config)
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages

# Import the canonical "interesting" verdict set from the verdict module so the
# PDF renderer stays in sync when compute_verdict's vocabulary evolves.
try:
    from scripts.s05.verdict import REPORT_INTERESTING_VERDICTS
except ImportError:  # pragma: no cover - allow standalone invocation
    _scripts_dir = Path(__file__).resolve().parents[1]
    if str(_scripts_dir) not in sys.path:
        sys.path.insert(0, str(_scripts_dir))
    from s05.verdict import REPORT_INTERESTING_VERDICTS


_SITE_RE = re.compile(r"^Insertion site:\s*(?P<site>.+)$", re.MULTILINE)
_VERDICT_RE = re.compile(r"^Verdict:\s*(?P<v>[A-Z_]+)", re.MULTILINE)
_LEN_RE = re.compile(r"^Insert length:\s*(?P<len>[\d,]+)\s*bp", re.MULTILINE)
_STRUCT_RE = re.compile(r"^Structure:\s*(?P<s>.+)$", re.MULTILINE)
_FOREIGN_BLOCK_RE = re.compile(
    r"Foreign elements \(not in host genome\):\s*\d+\s*\n"
    r"(?P<body>(?:\s*[+\-*]\s.+\n?)+)",
    re.MULTILINE,
)
_SITE_ANCHOR_RE = re.compile(
    r"^Insertion site:\s*(?P<chrom>[^\s:]+):\s*(?P<pos>[\d,]+)",
    re.MULTILINE,
)


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


def _parse_foreign_elements(text: str) -> list[str]:
    """Return the list of foreign-element names from a report body."""
    m = _FOREIGN_BLOCK_RE.search(text)
    if not m:
        return []
    out: list[str] = []
    for line in m.group("body").splitlines():
        line = line.strip()
        if not line:
            continue
        # lines look like "  + NODE_2_length_..."
        parts = line.split(None, 1)
        if len(parts) == 2:
            out.append(parts[1].strip())
    return out


def parse_insertion_reports(s05_dir: Path) -> list[dict[str, Any]]:
    """Parse ``insertion_*_report.txt`` into structured rows.

    Each row has: ``site``, ``verdict``, ``insert_length`` (string, digits+commas
    as in the source, or ``"—"``), ``structure``, ``elements`` (list[str]),
    and ``source`` (original filename).
    """
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
        structure = (m := _STRUCT_RE.search(text)) and m.group("s").strip() or ""
        elements = _parse_foreign_elements(text)
        rows.append({
            "site": site,
            "verdict": verdict,
            "insert_length": length,
            "structure": structure,
            "elements": elements,
            "source": report.name,
        })
    return rows


def extract_site_anchors(rows: list[dict[str, Any]]) -> list[tuple[str, int]]:
    """Extract (chrom, pos) anchors from parsed site strings for diagramming.

    The site field looks like ``"Chr3:16,439,674-16,439,709 (35bp deletion)"``
    or simply ``"Chr10:13500285"``. Rows with unparseable sites are skipped.
    """
    anchors: list[tuple[str, int]] = []
    for r in rows:
        site = str(r.get("site", ""))
        m = re.match(r"^([^\s:]+):\s*([\d,]+)", site)
        if not m:
            continue
        chrom = m.group(1)
        try:
            pos = int(m.group(2).replace(",", ""))
        except ValueError:
            continue
        anchors.append((chrom, pos))
    return anchors


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate verdict counts for the summary page."""
    counter = Counter(r["verdict"] for r in rows)
    return {
        "n_sites": len(rows),
        "by_verdict": dict(counter),
    }


def load_copynumber(sample_dir: Path) -> dict[str, Any] | None:
    """Parse ``s07_copynumber/copynumber.tsv`` into a dict (first data row).

    Returns ``None`` if the file is missing or empty.
    """
    path = Path(sample_dir) / "s07_copynumber" / "copynumber.tsv"
    if not path.exists():
        return None
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                return dict(row)
    except (OSError, csv.Error):
        return None
    return None


def load_coc_log(sample_dir: Path) -> list[dict[str, Any]]:
    """Parse ``chain_of_custody.jsonl`` into a list of JSON entries.

    Missing / unreadable file → empty list.  Malformed individual lines are
    silently skipped.
    """
    path = Path(sample_dir) / "chain_of_custody.jsonl"
    if not path.exists():
        return []
    out: list[dict[str, Any]] = []
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    out.append(json.loads(line))
                except json.JSONDecodeError:
                    continue
    except OSError:
        return []
    return out


def load_editing_data(
    sample_dir: Path,
) -> tuple[str, list[dict[str, Any]] | None, list[dict[str, Any]] | None]:
    """Load step-6 CRISPR outputs and return a 3-state tuple.

    State machine:
      "missing"  → s06_indel/ absent, grna_targets.tsv absent, no on-target rows,
                   or any TSV parse failure (treat as not-applicable)
      "no_edits" → on-target gRNAs present, editing_sites.tsv missing or empty
      "ok"       → on-target gRNAs present and at least one edit row
    """
    s06 = Path(sample_dir) / "s06_indel"
    targets_path = s06 / "grna_targets.tsv"
    sites_path = s06 / "editing_sites.tsv"

    if not targets_path.exists():
        return ("missing", None, None)

    try:
        with open(targets_path, newline="") as fh:
            targets = list(csv.DictReader(fh, delimiter="\t"))
    except (OSError, csv.Error):
        return ("missing", None, None)

    on_targets = [t for t in targets if t.get("site_type") == "on-target"]
    if not on_targets:
        return ("missing", None, None)

    sites: list[dict[str, Any]] = []
    if sites_path.exists():
        try:
            with open(sites_path, newline="") as fh:
                sites = list(csv.DictReader(fh, delimiter="\t"))
        except (OSError, csv.Error):
            sites = []

    if not sites:
        return ("no_edits", on_targets, [])
    return ("ok", on_targets, sites)


def load_kcgp_mapping(sample_dir: Path) -> dict[str, Any] | None:
    """Return {'kcgp_id', 'spec_version'} for the lowest-ordinal CANDIDATE row, or None.

    Soft-fails on missing file, IO error, or malformed CSV — the PDF still renders.
    """
    path = Path(sample_dir) / "kcgp_mapping.tsv"
    if not path.exists():
        return None
    try:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("verdict") == "CANDIDATE" and row.get("kcgp_id", "-") != "-":
                    return {
                        "kcgp_id": row["kcgp_id"],
                        "spec_version": row.get("spec_version", "tentative_v0"),
                    }
    except (OSError, csv.Error):
        return None
    return None


def build_sample_info(
    sample_dir: Path,
    *,
    sample_name: str,
    host: str | None = None,
    construct: str | None = None,
) -> dict[str, Any]:
    """Assemble the sample-information block (section 2)."""
    audit = load_audit_header(sample_dir)
    return {
        "sample": sample_name,
        "host": host or "(unspecified)",
        "construct": construct or "(unspecified)",
        "date": audit.get("generated_at", "—"),
        "input_sha256": audit.get("input_sha256", {}) or {},
    }


def build_pipeline_version(audit: dict[str, Any]) -> dict[str, Any]:
    """Assemble the pipeline-version block (section 3)."""
    return {
        "pipeline_commit": audit.get("pipeline_commit", "—"),
        "pipeline_dirty": audit.get("pipeline_dirty", "—"),
        "software_versions": audit.get("software_versions", {}) or {},
        "db_manifest": audit.get("db_manifest", []) or [],
    }


# ---------------------------------------------------------------------------
# PDF rendering primitives
# ---------------------------------------------------------------------------


def _page_text(pdf: PdfPages, title: str, body_lines: list[str],
               *, line_fontsize: int = 10, family: str = "monospace") -> None:
    """Render a single-column text page. Lines that overflow are truncated."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis("off")
    ax.text(0.05, 0.95, title, fontsize=18, fontweight="bold",
            transform=ax.transAxes, verticalalignment="top")
    y = 0.88
    step = 0.022 if line_fontsize <= 9 else 0.025
    for line in body_lines:
        ax.text(0.05, y, line, fontsize=line_fontsize, family=family,
                transform=ax.transAxes, verticalalignment="top")
        y -= step
        if y < 0.04:
            print(
                f"[insertion_pdf] WARN: page '{title}' truncated at line {body_lines.index(line) + 1}/{len(body_lines)}",
                file=sys.stderr,
            )
            break
    pdf.savefig(fig)
    plt.close(fig)


def _paginate(lines: list[str], *, page_lines: int = 38) -> list[list[str]]:
    """Split body lines into pages of ``page_lines`` each."""
    pages: list[list[str]] = []
    if not lines:
        return [[]]
    for i in range(0, len(lines), page_lines):
        pages.append(lines[i:i + page_lines])
    return pages


# ---------------------------------------------------------------------------
# Section renderers
# ---------------------------------------------------------------------------


def _page_cover(
    pdf: PdfPages,
    sample_name: str,
    audit: dict[str, Any],
    sample_dir: Path | None = None,
) -> None:
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
            f"  - {db.get('name')} md5={db.get('md5')} built={db.get('build_date')} "
            f"seqs={db.get('seq_count')}"
        )
    if not audit:
        lines.append("  (audit_header.json not present — placeholder page)")

    if sample_dir is not None:
        kcgp = load_kcgp_mapping(sample_dir)
        if kcgp:
            lines = [
                f"⚠ KCGP ID ({kcgp['spec_version']} — pending official spec):",
                f"    {kcgp['kcgp_id']}",
                "",
            ] + lines

    _page_text(pdf, "RedGene Insertion Report", lines)


def _page_sample_info(pdf: PdfPages, info: dict[str, Any]) -> None:
    lines = [
        f"Sample ID:   {info['sample']}",
        f"Host:        {info['host']}",
        f"Construct:   {info['construct']}",
        f"Generated:   {info['date']}",
        "",
        "Input SHA-256:",
    ]
    if info.get("input_sha256"):
        for k, v in info["input_sha256"].items():
            lines.append(f"  {k}: {v}")
    else:
        lines.append("  (no audit_header.json)")
    _page_text(pdf, "Sample Information", lines)


def _page_pipeline_version(pdf: PdfPages, ver: dict[str, Any]) -> None:
    lines = [
        f"Pipeline commit: {ver['pipeline_commit']}",
        f"Working tree dirty: {ver['pipeline_dirty']}",
        "",
        "Software versions:",
    ]
    sv = ver.get("software_versions") or {}
    if sv:
        for k, v in sv.items():
            # Trim very long version strings to keep layout tidy.
            v_str = str(v).splitlines()[0][:100]
            lines.append(f"  {k:<10s}  {v_str}")
    else:
        lines.append("  (audit_header.json not present)")
    lines += ["", "Database manifest:"]
    for db in ver.get("db_manifest") or []:
        lines.append(
            f"  - {db.get('name')} md5={db.get('md5')} "
            f"built={db.get('build_date')} seqs={db.get('seq_count')}"
        )
    if not ver.get("db_manifest"):
        lines.append("  (none)")
    _page_text(pdf, "Pipeline Version", lines)


def _page_summary(pdf: PdfPages, summary: dict[str, Any]) -> None:
    lines = [
        f"Total sites parsed: {summary['n_sites']}",
        "",
        "Verdict distribution:",
    ]
    for v, n in sorted(summary.get("by_verdict", {}).items()):
        lines.append(f"  {v:<20s}  {n}")
    _page_text(pdf, "Sample Summary", lines)


def _sites_table_rows(rows: list[dict[str, Any]]) -> list[str]:
    header = (f"{'site':<38s}  {'verdict':<16s}  {'len_bp':>8s}  "
              f"{'structure':<24s}  elements")
    out = [header, "-" * 110]
    for r in rows:
        site = str(r["site"])[:38]
        verdict = str(r["verdict"])[:16]
        length = str(r["insert_length"])[:8]
        structure = str(r.get("structure", ""))[:24]
        elts = ",".join(r.get("elements") or [])[:28]
        out.append(f"{site:<38s}  {verdict:<16s}  {length:>8s}  "
                   f"{structure:<24s}  {elts}")
    if not rows:
        out.append("(no insertion reports found)")
    return out


def _page_sites_table(pdf: PdfPages, rows: list[dict[str, Any]]) -> int:
    """Render the insertion-sites table, paginating if >40 entries. Returns pages."""
    lines = _sites_table_rows(rows)
    pages = _paginate(lines, page_lines=42)
    for i, chunk in enumerate(pages, start=1):
        title = "Insertion Sites" if len(pages) == 1 \
            else f"Insertion Sites ({i}/{len(pages)})"
        _page_text(pdf, title, chunk, line_fontsize=8)
    return len(pages)


def _page_junction_diagram(pdf: PdfPages, rows: list[dict[str, Any]]) -> None:
    """Render a simple junction overview: scatter of anchor positions per chrom.

    This is NOT a per-read soft-clip pileup (that needs BAM access); it's a
    genome-scale overview of candidate anchors so reviewers can spot clustering
    and check against known ground-truth (e.g., rice Chr3:16,439,674).
    """
    # Filter to interesting verdicts (imported from scripts.s05.verdict so the
    # set stays in sync when compute_verdict adds / renames verdict labels).
    interesting_rows = [r for r in rows if r["verdict"] in REPORT_INTERESTING_VERDICTS]
    anchors = extract_site_anchors(interesting_rows or rows)

    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.set_title("Junction Diagram — candidate anchor positions",
                 fontsize=16, fontweight="bold", loc="left")

    if not anchors:
        ax.text(0.5, 0.5, "No candidate junctions parsed from s05 reports.",
                ha="center", va="center", fontsize=12, transform=ax.transAxes)
        ax.axis("off")
        pdf.savefig(fig)
        plt.close(fig)
        return

    chroms = sorted({c for c, _ in anchors})
    chrom_to_y = {c: i for i, c in enumerate(chroms)}

    for chrom, pos in anchors:
        y = chrom_to_y[chrom]
        ax.plot([pos], [y], marker="|", markersize=18,
                color="tab:red", linestyle="none")
    ax.set_yticks(list(chrom_to_y.values()))
    ax.set_yticklabels(chroms, fontsize=9)
    ax.set_xlabel("Genomic position (bp)")
    ax.set_ylabel("Chromosome / contig")
    ax.set_ylim(-1, len(chroms))
    ax.grid(True, axis="x", linestyle=":", alpha=0.4)
    ax.ticklabel_format(axis="x", style="plain")

    # Annotate the first few anchors with their coordinate.
    for chrom, pos in anchors[:8]:
        y = chrom_to_y[chrom]
        ax.annotate(f"{pos:,}", xy=(pos, y), xytext=(4, 6),
                    textcoords="offset points", fontsize=7)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _page_copy_number(pdf: PdfPages, sample_dir: Path) -> None:
    """Render a bar chart from copynumber.tsv, or a placeholder."""
    cn = load_copynumber(sample_dir)
    if not cn:
        _page_text(pdf, "Copy Number", ["(s07_copynumber/copynumber.tsv not present)"])
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 11),
                                   gridspec_kw={"height_ratios": [1, 2]})
    ax1.axis("off")
    sample = cn.get("sample", "—")
    desc = cn.get("copy_description", "—")
    conf = cn.get("confidence", "—")
    val = cn.get("site_validation", "—")
    ax1.text(0.0, 0.95, "Copy Number", fontsize=18, fontweight="bold",
             transform=ax1.transAxes, verticalalignment="top")
    lines = [
        f"Sample:            {sample}",
        f"Description:       {desc}",
        f"Estimated copies:  {cn.get('estimated_copies', '—')}",
        f"Depth ratio:       {cn.get('depth_ratio', '—')}",
        f"Candidate sites:   {cn.get('candidate_sites', '—')}",
        f"Site validation:   {val}",
        f"Confidence:        {conf}",
    ]
    y = 0.80
    for line in lines:
        ax1.text(0.02, y, line, fontsize=11, family="monospace",
                 transform=ax1.transAxes, verticalalignment="top")
        y -= 0.08

    # Bar chart: construct vs host median depth on the primary axis, depth
    # ratio on a secondary axis. Depth values (20-50 reads) and ratio (0.5-3.0)
    # are orders of magnitude apart; without twinx() the ratio bar vanishes.
    try:
        cd = float(cn.get("construct_median_depth", 0) or 0)
        hd = float(cn.get("host_median_depth", 0) or 0)
        ratio = float(cn.get("depth_ratio", 0) or 0)
    except ValueError:
        cd, hd, ratio = 0.0, 0.0, 0.0

    depth_labels = ["Construct median", "Host median"]
    depth_values = [cd, hd]
    x_depth = list(range(len(depth_labels)))
    bars_depth = ax2.bar(x_depth, depth_values,
                         color=["tab:blue", "tab:gray"], width=0.6)
    ax2.set_xticks(x_depth + [len(x_depth)])
    ax2.set_xticklabels(depth_labels + ["Depth ratio"])
    ax2.set_ylabel("Depth (reads)")
    ax2.set_title("Construct vs Host Coverage", fontsize=12)
    for bar, v in zip(bars_depth, depth_values):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                 f"{v:.2f}", ha="center", va="bottom", fontsize=10)
    ax2.grid(True, axis="y", linestyle=":", alpha=0.4)

    ax2r = ax2.twinx()
    bar_ratio = ax2r.bar([len(x_depth)], [ratio], color="tab:red", width=0.6)
    ax2r.set_ylabel("Depth ratio (construct / host)")
    for bar in bar_ratio:
        ax2r.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                  f"{ratio:.2f}", ha="center", va="bottom", fontsize=10)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def _page_placeholder(pdf: PdfPages, title: str, note: str) -> None:
    _page_text(pdf, title, [note, "", "(scaffold — to be wired in v1.1)"])


def _page_crispr_summary(
    pdf: PdfPages,
    on_targets: list[dict[str, Any]],
    sites: list[dict[str, Any]],
) -> None:
    """Section 8 summary page: gRNA inventory + per-gRNA edit highlights."""
    n_grna = len(on_targets)
    n_edits = len(sites)
    lines = [
        f"gRNAs (on-target): {n_grna}",
        f"Editing events:    {n_edits}",
        "",
        f"  {'gRNA':<5} {'Locus':<32} {'Edits':>6} {'Top freq':>10} {'Top type':<12}",
        f"  {'-'*5} {'-'*32} {'-'*6} {'-'*10} {'-'*12}",
    ]
    for t in on_targets:
        idx = t.get("grna_idx", "?")
        chrom = t.get("chrom", "?")
        cut = t.get("cut_pos", "?")
        strand = t.get("strand", "?")
        locus = f"{chrom}:{cut} ({strand})"
        if len(locus) > 32:
            locus = locus[:29] + "..."
        edits = [s for s in sites if s.get("grna_idx") == idx]
        if edits:
            top = max(edits, key=lambda r: float(r.get("freq", "0") or 0))
            try:
                freq_str = f"{float(top.get('freq', '0') or 0):.3f}"
            except ValueError:
                freq_str = "—"
            top_type = top.get("type", "—")
        else:
            freq_str = "—"
            top_type = "—"
        lines.append(
            f"  {str(idx):<5} {locus:<32} {len(edits):>6d} {freq_str:>10} {top_type:<12}"
        )
    _page_text(pdf, "CRISPR Editing — Summary", lines)


def _page_crispr_profile(
    pdf: PdfPages,
    sample_dir: Path,
    target: dict[str, Any],
) -> None:
    """Render one PDF page embedding editing_profile_gRNA{idx}.png.

    Falls back to a text page if the PNG is missing or unreadable.
    """
    idx = str(target["grna_idx"])
    seq = target["grna_seq"]
    chrom = target["chrom"]
    cut = target["cut_pos"]
    title = f"CRISPR Editing — gRNA {idx}"
    header = f"gRNA {idx} — {chrom}:{cut} ({seq})"

    png_path = Path(sample_dir) / "s06_indel" / f"editing_profile_gRNA{idx}.png"
    if not png_path.exists():
        _page_text(pdf, title, [
            header, "",
            f"(editing_profile_gRNA{idx}.png missing — "
            "scripts/viz/plot_editing_profile.py needs to run first)",
        ])
        return

    try:
        img = mpimg.imread(str(png_path))
    except (OSError, ValueError, SyntaxError):
        _page_text(pdf, title, [
            header, "",
            f"(editing_profile_gRNA{idx}.png unreadable — file may be corrupt)",
        ])
        return

    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis("off")
    ax.text(0.05, 0.96, title, fontsize=18, fontweight="bold",
            transform=ax.transAxes, verticalalignment="top")
    ax.text(0.05, 0.92, header, fontsize=10, family="monospace",
            transform=ax.transAxes, verticalalignment="top")
    img_ax = fig.add_axes((0.05, 0.05, 0.90, 0.83))
    img_ax.imshow(img)
    img_ax.axis("off")
    pdf.savefig(fig)
    plt.close(fig)


def _page_crispr_editing(pdf: PdfPages, sample_dir: Path) -> int:
    """Section 8 dispatcher. Returns the number of pages written (>=1)."""
    state, on_targets, sites = load_editing_data(sample_dir)

    if state == "missing":
        _page_text(pdf, "CRISPR Editing", [
            "(s06_indel CRISPR analysis was not run for this sample,",
            " or no on-target gRNA was defined.)",
        ])
        return 1

    if state == "no_edits":
        n = len(on_targets) if on_targets else 0
        _page_text(pdf, "CRISPR Editing", [
            f"On-target gRNAs analyzed: {n}",
            "Editing events detected:  0",
            "",
            "(no editing events passed step-6 thresholds)",
        ])
        return 1

    # state == "ok"
    pages = 0
    _page_crispr_summary(pdf, on_targets or [], sites or [])
    pages += 1
    for target in on_targets or []:
        _page_crispr_profile(pdf, sample_dir, target)
        pages += 1
    return pages


def _page_appendix_audit(pdf: PdfPages, audit: dict[str, Any]) -> None:
    lines = [
        f"Pipeline commit: {audit.get('pipeline_commit', '—')}",
        f"Dirty: {audit.get('pipeline_dirty', '—')}",
        "",
        "Software versions:",
    ]
    for k, v in (audit.get("software_versions") or {}).items():
        v_str = str(v).splitlines()[0][:100]
        lines.append(f"  {k:<10s}  {v_str}")
    lines += ["", "DB manifest:"]
    for db in audit.get("db_manifest", []) or []:
        lines.append(
            f"  - {db.get('name')} md5={db.get('md5')} "
            f"built={db.get('build_date')} seqs={db.get('seq_count')}"
        )
    if not audit:
        lines.append("(audit_header.json absent — appendix placeholder)")
    _page_text(pdf, "Appendix A — Audit Header", lines)


def _page_appendix_coc(pdf: PdfPages, entries: list[dict[str, Any]]) -> int:
    """Render the chain-of-custody log digest.

    Each entry is summarised to one line; long logs paginate. Returns number of
    pages emitted (0 if the log is empty / absent).
    """
    if not entries:
        return 0
    entries = sorted(entries, key=lambda e: e.get("ts", ""))
    body = ["event  step    ts                                  extra"]
    body.append("-" * 110)
    for e in entries:
        event = str(e.get("event", "—"))[:6]
        step = str(e.get("step", "—"))[:6]
        ts = str(e.get("ts", "—"))[:30]
        extra_bits = []
        if "exit_code" in e:
            extra_bits.append(f"exit={e['exit_code']}")
        if "wall_time_sec" in e:
            extra_bits.append(f"wall={e['wall_time_sec']}s")
        if "cmd" in e:
            extra_bits.append("cmd=" + str(e["cmd"])[:60])
        body.append(f"{event:<6s} {step:<6s} {ts:<32s}  " + " ".join(extra_bits))
    pages = _paginate(body, page_lines=44)
    for i, chunk in enumerate(pages, start=1):
        title = "Appendix B — Chain of Custody" if len(pages) == 1 \
            else f"Appendix B — Chain of Custody ({i}/{len(pages)})"
        _page_text(pdf, title, chunk, line_fontsize=7)
    return len(pages)


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


def _load_construct_inventory(s05b_dir: Path) -> dict[str, Any] | None:
    """Parse s05b outputs into {summary: dict, links: list[dict]} or None.

    Returns None when the step-5b output directory / summary is absent, so the
    PDF stays valid for samples that never ran step 5b.
    """
    s05b_dir = Path(s05b_dir)
    summary_p = s05b_dir / "construct_summary.txt"
    links_p = s05b_dir / "site_construct_links.tsv"
    if not summary_p.is_file():
        return None
    summary: dict[str, str] = {}
    for line in summary_p.read_text().splitlines():
        if line.startswith("#") or "\t" not in line:
            continue
        k, v = line.split("\t", 1)
        summary[k] = v
    links: list[dict[str, str]] = []
    if links_p.is_file():
        rows = links_p.read_text().splitlines()
        if rows:
            header = rows[0].split("\t")
            for r in rows[1:]:
                cells = r.split("\t")
                if len(cells) == len(header):
                    links.append(dict(zip(header, cells)))
    return {"summary": summary, "links": links}


def _page_construct_inventory(
    pdf: PdfPages, inv: dict[str, Any], sample_name: str
) -> None:
    """Render the Construct Inventory section (step 5b) as a text page."""
    summary = inv.get("summary", {})
    links = inv.get("links", [])
    body: list[str] = []
    body.append("Reconstructed transgene construct (step 5b)")
    body.append("")
    body.append(f"  Total construct bp     : {summary.get('total_construct_bp', 'NA')}")
    body.append(f"  Distinct elements      : {summary.get('distinct_elements', 'NA')}")
    body.append(f"  Element inventory      : {summary.get('element_inventory', 'NA')}")
    body.append(f"  Novel payload detected : {summary.get('novel_payload_detected', 'NA')}")
    body.append("")
    body.append("  Contigs by class:")
    for cls in ("construct", "chimeric", "foreign_unknown", "host"):
        body.append(f"    {cls:<16}: {summary.get(f'contigs_{cls}', '0')}")
    body.append("")
    body.append(f"  Sites linked to construct contigs: {summary.get('sites_linked', '0')}")
    body.append("")
    if links:
        body.append(f"  {'site_id':<22} {'contig':<14} {'side':<11} {'conf':<7}")
        for lk in links:
            body.append(
                f"  {lk.get('site_id', ''):<22} {lk.get('contig_id', ''):<14} "
                f"{lk.get('junction_side', ''):<11} {lk.get('link_confidence', ''):<7}"
            )
    else:
        body.append("  (no construct contigs linked to CANDIDATE sites)")
    _page_text(pdf, f"Construct Inventory — {sample_name}", body)


def generate_pdf(
    *,
    sample_dir: Path,
    sample_name: str,
    out_pdf: Path,
    host: str | None = None,
    construct: str | None = None,
) -> int:
    """Render the full PDF. Returns the number of pages written."""
    sample_dir = Path(sample_dir)
    out_pdf = Path(out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    audit = load_audit_header(sample_dir)
    rows = parse_insertion_reports(sample_dir / "s05_insert_assembly")
    summary = summarize(rows)
    sample_info = build_sample_info(sample_dir, sample_name=sample_name,
                                    host=host, construct=construct)
    pipeline_ver = build_pipeline_version(audit)
    coc_entries = load_coc_log(sample_dir)

    pages = 0
    with PdfPages(out_pdf) as pdf:
        _page_cover(pdf, sample_name, audit, sample_dir=sample_dir); pages += 1  # (1)
        _page_sample_info(pdf, sample_info); pages += 1                 # (2)
        _page_pipeline_version(pdf, pipeline_ver); pages += 1           # (3)
        _page_summary(pdf, summary); pages += 1                         # (4)
        pages += _page_sites_table(pdf, rows)                           # (5)
        _page_junction_diagram(pdf, rows); pages += 1                   # (6)
        construct_inv = _load_construct_inventory(
            sample_dir / "s05b_construct_assembly")
        if construct_inv is not None:
            _page_construct_inventory(pdf, construct_inv, sample_name); pages += 1  # (6b)
        _page_copy_number(pdf, sample_dir); pages += 1                  # (7)
        pages += _page_crispr_editing(pdf, sample_dir)                  # (8)
        _page_appendix_audit(pdf, audit); pages += 1                    # (9)
        pages += _page_appendix_coc(pdf, coc_entries)                   # (10, opt)
    return pages


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _resolve_host_construct_from_config(
    config_path: Path | None, sample_key: str
) -> tuple[str | None, str | None]:
    """Look up host_reference / construct_reference names from config.yaml."""
    if config_path is None or not Path(config_path).exists():
        return None, None
    try:
        import yaml  # local import; optional
    except Exception:  # pragma: no cover
        return None, None
    try:
        cfg = yaml.safe_load(Path(config_path).read_text()) or {}
    except Exception:  # pragma: no cover
        return None, None
    samples = cfg.get("samples", {}) or {}
    entry = samples.get(sample_key) or {}
    host = entry.get("host_reference")
    construct = entry.get("construct_reference")
    # Use just the filename stem for readability.
    if host:
        host = Path(host).name
    if construct:
        construct = Path(construct).name
    return host, construct


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=False, default=None,
                   help="Sample name (for page titles). Defaults to --sample-dir basename.")
    p.add_argument("--sample-dir", type=Path, required=True,
                   help="results/<sample>/ directory")
    p.add_argument("--out", type=Path, required=True,
                   help="Output PDF path (parent directory auto-created)")
    p.add_argument("--config", type=Path, default=None,
                   help="Optional path to config.yaml for host/construct lookup")
    p.add_argument("--host", type=str, default=None,
                   help="Override the host label in the sample-info section")
    p.add_argument("--construct", type=str, default=None,
                   help="Override the construct label in the sample-info section")
    args = p.parse_args(argv)

    sample_name = args.sample or args.sample_dir.name
    host = args.host
    construct = args.construct
    if (host is None or construct is None) and args.config is not None:
        cfg_host, cfg_construct = _resolve_host_construct_from_config(
            args.config, sample_name
        )
        host = host or cfg_host
        construct = construct or cfg_construct

    n = generate_pdf(
        sample_dir=args.sample_dir,
        sample_name=sample_name,
        out_pdf=args.out,
        host=host,
        construct=construct,
    )
    print(f"[insertion-pdf] wrote {n} pages → {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
