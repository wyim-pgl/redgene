#!/usr/bin/env python3
"""Step 11: MultiQC report generation.

Aggregates QC metrics from all pipeline steps into a single HTML report.
Collects outputs from:
  - fastp (step 1): adapter trimming, quality metrics
  - samtools flagstat (steps 2, 7): mapping statistics
  - samtools stats (step 7): insert size, GC content, coverage
  - SPAdes (step 4): assembly statistics (custom section)
  - Junction detection (step 6): summary table (custom section)
  - CRISPR editing (step 8): editing sites summary (custom section)
  - Copy number (step 10): estimates (custom section)

Usage:
  python s11_multiqc.py --outdir results --sample-name rice_G281
  python s11_multiqc.py --outdir results  # all samples
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def log(msg: str) -> None:
    print(f"[s11_multiqc] {msg}", file=sys.stderr, flush=True)


def generate_custom_data(outdir: Path, samples: list[str]) -> Path:
    """Generate MultiQC custom content files for pipeline-specific results."""

    mqc_dir = outdir / "multiqc_custom_data"
    mqc_dir.mkdir(parents=True, exist_ok=True)

    # ---- Insert assembly summary (from s05_stats.txt) ----
    site_rows = []
    for sample in samples:
        stats_file = outdir / sample / "s05_insert_assembly" / "s05_stats.txt"
        if not stats_file.exists():
            continue
        stats = {}
        with open(stats_file) as f:
            for line in f:
                parts = line.strip().split("\t", 1)
                if len(parts) == 2:
                    stats[parts[0]] = parts[1]

        n_sites = int(stats.get("insertion_sites", 0))
        n_candidates = sum(1 for k, v in stats.items()
                           if k.endswith("_verdict") and v == "CANDIDATE")
        n_fp = sum(1 for k, v in stats.items()
                   if k.endswith("_verdict") and v == "FALSE_POSITIVE")
        site_rows.append({
            "Sample": sample,
            "Sites": n_sites,
            "Candidates": n_candidates,
            "False positives": n_fp,
        })

    if site_rows:
        site_mqc = mqc_dir / "insert_assembly_mqc.tsv"
        with open(site_mqc, "w") as f:
            f.write("# id: 'insert_assembly'\n")
            f.write("# section_name: 'Insert Assembly'\n")
            f.write("# description: 'Targeted insert assembly — site detection and filtering'\n")
            f.write("# plot_type: 'table'\n")
            f.write("# pconfig:\n")
            f.write("#   id: 'insert_assembly_table'\n")
            headers = list(site_rows[0].keys())
            f.write("\t".join(headers) + "\n")
            for row in site_rows:
                f.write("\t".join(str(row[h]) for h in headers) + "\n")

    # ---- CRISPR editing summary ----
    edit_rows = []
    for sample in samples:
        edit_file = outdir / sample / "s06_indel" / "editing_sites.tsv"
        if not edit_file.exists():
            continue
        with open(edit_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            edits = list(reader)

        if not edits:
            continue

        for e in edits:
            edit_rows.append({
                "Sample": sample,
                "Chrom": e.get("chrom", ""),
                "Position": e.get("pos", ""),
                "Type": e.get("type", ""),
                "Size": e.get("size", ""),
                "Indel seq": e.get("indel_seq", ""),
                "Frequency": f"{float(e.get('freq', 0)):.1%}",
                "Depth": e.get("dp", ""),
                "Count": e.get("count", ""),
                "Zygosity": e.get("zygosity", ""),
                "gRNA": e.get("grna_idx", ""),
                "Site type": e.get("site_type", ""),
            })

    if edit_rows:
        edit_mqc = mqc_dir / "editing_sites_mqc.tsv"
        with open(edit_mqc, "w") as f:
            f.write("# id: 'editing_sites'\n")
            f.write("# section_name: 'CRISPR Editing Sites'\n")
            f.write("# description: 'Treatment-specific indels at gRNA target sites'\n")
            f.write("# plot_type: 'table'\n")
            f.write("# pconfig:\n")
            f.write("#   id: 'editing_sites_table'\n")
            headers = list(edit_rows[0].keys())
            f.write("\t".join(headers) + "\n")
            for row in edit_rows:
                f.write("\t".join(str(row[h]) for h in headers) + "\n")

    # ---- Copy number summary ----
    cn_rows = []
    for sample in samples:
        cn_file = outdir / sample / "s07_copynumber" / "copynumber.tsv"
        if not cn_file.exists():
            continue
        with open(cn_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                cn_rows.append({
                    "Sample": sample,
                    "Copy number": row.get("copy_number", ""),
                    "Construct depth": row.get("construct_depth", ""),
                    "Host depth": row.get("host_depth", ""),
                    "Marker": row.get("best_marker", ""),
                    "Confidence": row.get("confidence", ""),
                })

    if cn_rows:
        cn_mqc = mqc_dir / "copynumber_mqc.tsv"
        with open(cn_mqc, "w") as f:
            f.write("# id: 'copy_number'\n")
            f.write("# section_name: 'Copy Number Estimation'\n")
            f.write("# description: 'Depth-ratio copy number estimates'\n")
            f.write("# plot_type: 'table'\n")
            f.write("# pconfig:\n")
            f.write("#   id: 'copy_number_table'\n")
            headers = list(cn_rows[0].keys())
            f.write("\t".join(headers) + "\n")
            for row in cn_rows:
                f.write("\t".join(str(row[h]) for h in headers) + "\n")

    # ---- Variant effects summary ----
    eff_rows = []
    for sample in samples:
        eff_file = outdir / sample / "s06_indel" / "editing_effects.tsv"
        if not eff_file.exists():
            continue
        with open(eff_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                eff_rows.append({
                    "Sample": sample,
                    "Position": f"{row.get('chrom', '')}:{row.get('pos', '')}",
                    "Variant": f"{row.get('type', '')} {row.get('size', '')}bp",
                    "Effect": row.get("effect", ""),
                    "Gene": row.get("gene", ""),
                    "AA change": row.get("aa_change", ""),
                })

    if eff_rows:
        eff_mqc = mqc_dir / "variant_effects_mqc.tsv"
        with open(eff_mqc, "w") as f:
            f.write("# id: 'variant_effects'\n")
            f.write("# section_name: 'Variant Effects'\n")
            f.write("# description: 'Functional impact of CRISPR editing'\n")
            f.write("# plot_type: 'table'\n")
            f.write("# pconfig:\n")
            f.write("#   id: 'variant_effects_table'\n")
            headers = list(eff_rows[0].keys())
            f.write("\t".join(headers) + "\n")
            for row in eff_rows:
                f.write("\t".join(str(row[h]) for h in headers) + "\n")

    log(f"Custom data written to: {mqc_dir}")
    return mqc_dir


def run_multiqc(
    outdir: Path,
    samples: list[str],
    report_title: str = "RedGene Pipeline Report",
) -> Path:
    """Run MultiQC on pipeline outputs."""

    # Collect search directories (only main samples, not subsampled)
    search_dirs = []
    for sample in samples:
        sample_dir = outdir / sample
        if sample_dir.exists():
            search_dirs.append(str(sample_dir))

    if not search_dirs:
        log("No sample directories found")
        sys.exit(1)

    # Generate custom data
    custom_dir = generate_custom_data(outdir, samples)

    # Add custom data directory
    search_dirs.append(str(custom_dir))

    # MultiQC output
    mqc_outdir = outdir / "multiqc"
    mqc_outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "multiqc",
        "--force",
        "--outdir", str(mqc_outdir),
        "--title", report_title,
        "--filename", "redgene_report",
        "--no-data-dir",
    ] + search_dirs

    log(f"Running MultiQC on {len(search_dirs)} directories...")
    log(f"  CMD: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        log(f"MultiQC stderr: {result.stderr}")
        sys.exit(f"MultiQC failed with code {result.returncode}")

    # Print MultiQC stdout (module summary)
    for line in result.stderr.strip().split("\n"):
        if line.strip():
            log(f"  {line.strip()}")

    report_path = mqc_outdir / "redgene_report.html"
    log(f"Report: {report_path}")

    return report_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 11: MultiQC report generation",
    )
    parser.add_argument(
        "--outdir", type=Path, required=True,
        help="Base output directory (e.g., results)",
    )
    parser.add_argument(
        "--sample-name", type=str, default=None,
        help="Specific sample to report (default: auto-detect all)",
    )
    parser.add_argument(
        "--title", type=str, default="RedGene Pipeline Report",
        help="Report title",
    )
    args = parser.parse_args()

    # Discover samples
    if args.sample_name:
        samples = [args.sample_name]
    else:
        # Auto-detect: directories in outdir that have s01_qc
        samples = []
        for d in sorted(args.outdir.iterdir()):
            if d.is_dir() and (d / "s01_qc").exists():
                # Skip subsampled coverage test samples
                name = d.name
                if any(x in name for x in ("_15x", "_10x", "_5x", "_3x")):
                    continue
                samples.append(name)

    if not samples:
        log("No samples found")
        sys.exit(1)

    log(f"Samples: {', '.join(samples)}")

    report = run_multiqc(args.outdir, samples, args.title)
    log(f"Done. Report at: {report}")


if __name__ == "__main__":
    main()
