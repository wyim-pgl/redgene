#!/usr/bin/env python3
"""Step 1: QC and adapter trimming with fastp.

Runs fastp on paired-end Illumina reads to perform quality control,
adapter trimming, and quality filtering.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run fastp QC on paired-end reads"
    )
    parser.add_argument("--r1", type=Path, required=True,
                        help="Input forward reads (FASTQ, gzipped)")
    parser.add_argument("--r2", type=Path, required=True,
                        help="Input reverse reads (FASTQ, gzipped)")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="Base output directory")
    parser.add_argument("--threads", type=int, default=8,
                        help="Number of threads (default: 8)")
    parser.add_argument("--sample-name", type=str, required=True,
                        help="Sample name for output file naming")
    return parser.parse_args()


def run_fastp(r1: Path, r2: Path, outdir: Path, threads: int,
              sample_name: str) -> None:
    """Run fastp and report summary statistics."""
    # Create output directory
    step_dir = outdir / sample_name / "s01_qc"
    step_dir.mkdir(parents=True, exist_ok=True)

    out_r1 = step_dir / f"{sample_name}_R1.fq.gz"
    out_r2 = step_dir / f"{sample_name}_R2.fq.gz"
    html_report = step_dir / "fastp.html"
    json_report = step_dir / "fastp.json"

    # Validate inputs
    if not r1.exists():
        print(f"ERROR: R1 file not found: {r1}", file=sys.stderr)
        sys.exit(1)
    if not r2.exists():
        print(f"ERROR: R2 file not found: {r2}", file=sys.stderr)
        sys.exit(1)

    cmd = [
        "fastp",
        "--in1", str(r1),
        "--in2", str(r2),
        "--out1", str(out_r1),
        "--out2", str(out_r2),
        "--html", str(html_report),
        "--json", str(json_report),
        "--thread", str(threads),
        "--detect_adapter_for_pe",
        "--qualified_quality_phred", "20",
        "--length_required", "50",
    ]

    print(f"[s01_qc] Running fastp on sample: {sample_name}", file=sys.stderr)
    print(f"[s01_qc] Input R1: {r1}", file=sys.stderr)
    print(f"[s01_qc] Input R2: {r2}", file=sys.stderr)
    print(f"[s01_qc] Output dir: {step_dir}", file=sys.stderr)

    subprocess.run(cmd, check=True)

    # Parse and report stats from JSON
    if json_report.exists():
        with open(json_report) as fh:
            stats = json.load(fh)

        summary = stats.get("summary", {})
        before = summary.get("before_filtering", {})
        after = summary.get("after_filtering", {})

        print(f"[s01_qc] === QC Summary for {sample_name} ===", file=sys.stderr)
        print(f"[s01_qc] Before filtering:", file=sys.stderr)
        print(f"[s01_qc]   Total reads: {before.get('total_reads', 'N/A'):,}",
              file=sys.stderr)
        print(f"[s01_qc]   Total bases: {before.get('total_bases', 'N/A'):,}",
              file=sys.stderr)
        print(f"[s01_qc]   Q20 rate: {before.get('q20_rate', 'N/A')}",
              file=sys.stderr)
        print(f"[s01_qc]   Q30 rate: {before.get('q30_rate', 'N/A')}",
              file=sys.stderr)
        print(f"[s01_qc] After filtering:", file=sys.stderr)
        print(f"[s01_qc]   Total reads: {after.get('total_reads', 'N/A'):,}",
              file=sys.stderr)
        print(f"[s01_qc]   Total bases: {after.get('total_bases', 'N/A'):,}",
              file=sys.stderr)
        print(f"[s01_qc]   Q20 rate: {after.get('q20_rate', 'N/A')}",
              file=sys.stderr)
        print(f"[s01_qc]   Q30 rate: {after.get('q30_rate', 'N/A')}",
              file=sys.stderr)

        # Report filtering stats
        filtering = stats.get("filtering_result", {})
        passed = filtering.get("passed_filter_reads", 0)
        low_quality = filtering.get("low_quality_reads", 0)
        too_short = filtering.get("too_short_reads", 0)
        adapter_trimmed = stats.get("adapter_cutting", {}).get(
            "adapter_trimmed_reads", 0
        )
        print(f"[s01_qc]   Passed filter: {passed:,}", file=sys.stderr)
        print(f"[s01_qc]   Low quality: {low_quality:,}", file=sys.stderr)
        print(f"[s01_qc]   Too short: {too_short:,}", file=sys.stderr)
        print(f"[s01_qc]   Adapter trimmed: {adapter_trimmed:,}",
              file=sys.stderr)

    # Verify outputs exist
    for f in [out_r1, out_r2, html_report, json_report]:
        if not f.exists():
            print(f"ERROR: Expected output not created: {f}", file=sys.stderr)
            sys.exit(1)

    print(f"[s01_qc] Done. Outputs in {step_dir}", file=sys.stderr)


def main() -> None:
    args = parse_args()
    run_fastp(args.r1, args.r2, args.outdir, args.threads, args.sample_name)


if __name__ == "__main__":
    main()
