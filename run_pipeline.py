#!/usr/bin/env python3
"""GMO Positive Control Characterization Pipeline - Main Entry Point.

Runs the pipeline steps sequentially for each sample defined in config.yaml.
Each step invokes the corresponding script in scripts/ via subprocess.
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STEP_SCRIPTS: dict[str, str] = {
    "1": "scripts/s01_qc.py",
    "2": "scripts/s02_construct_map.py",
    "3": "scripts/s03_extract_reads.py",
    "4": "scripts/s04_assembly.py",
    "5": "scripts/s05_contig_map.py",
    "6": "scripts/s06_junction.py",
    "7": "scripts/s07_host_map.py",
    "8": "scripts/s08_indel_detection.py",
    "10": "scripts/s10_copynumber.py",
    "11": "scripts/s11_multiqc.py",
}

STEP_NAMES: dict[str, str] = {
    "1": "QC + trim (fastp)",
    "2": "Map reads to construct (bwa mem)",
    "3": "Extract construct-hitting reads + mates",
    "4": "Local assembly (SPAdes)",
    "5": "Map contigs to host (minimap2)",
    "6": "Chimeric contig / junction detection",
    "7": "Map all reads to host (bwa mem)",
    "8": "CRISPR indel detection (treatment vs WT)",
    "10": "Copy number estimation",
    "11": "MultiQC report generation",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stderr,
)
log = logging.getLogger("pipeline")


def parse_steps(step_spec: str) -> list[str]:
    """Parse a step specification like '1-6', '1,2,5', '1-6,10' into a sorted
    list of step keys.  Only steps present in STEP_SCRIPTS are accepted."""
    steps: set[str] = set()
    for part in step_spec.split(","):
        part = part.strip()
        if "-" in part:
            lo, hi = part.split("-", 1)
            lo_int, hi_int = int(lo), int(hi)
            for i in range(lo_int, hi_int + 1):
                steps.add(str(i))
        else:
            steps.add(part)

    # Validate
    unknown = steps - set(STEP_SCRIPTS)
    if unknown:
        sys.exit(f"ERROR: Unknown step(s): {', '.join(sorted(unknown, key=int))}. "
                 f"Available: {', '.join(sorted(STEP_SCRIPTS, key=int))}")

    # Return in execution order
    return sorted(steps, key=int)


def load_config(config_path: Path) -> dict[str, Any]:
    """Load and return the YAML configuration."""
    if not config_path.exists():
        sys.exit(f"ERROR: Config file not found: {config_path}")
    with open(config_path) as fh:
        cfg: dict[str, Any] = yaml.safe_load(fh)
    if "samples" not in cfg or not cfg["samples"]:
        sys.exit("ERROR: No samples defined in config.")
    return cfg


def resolve_samples(cfg: dict[str, Any], sample_arg: str | None) -> list[str]:
    """Return list of sample keys to process."""
    all_samples = list(cfg["samples"].keys())
    if sample_arg is None or sample_arg.lower() == "all":
        return all_samples
    if sample_arg not in cfg["samples"]:
        sys.exit(f"ERROR: Sample '{sample_arg}' not found in config. "
                 f"Available: {', '.join(all_samples)}")
    return [sample_arg]


def validate_sample_inputs(cfg: dict[str, Any], sample_key: str, base_dir: Path) -> None:
    """Check that essential input files for a sample exist."""
    sample = cfg["samples"][sample_key]
    reads = sample.get("reads", {})
    for tag in ("r1", "r2"):
        fq = reads.get(tag)
        if fq is None:
            sys.exit(f"ERROR: Sample '{sample_key}' missing reads.{tag} in config.")
        fq_path = base_dir / fq
        if not fq_path.exists():
            log.warning("Read file not found (may be created later): %s", fq_path)

    construct = sample.get("construct_reference")
    if construct is None:
        sys.exit(f"ERROR: Sample '{sample_key}' missing construct_reference in config.")

    host = sample.get("host_reference")
    if host is None:
        sys.exit(f"ERROR: Sample '{sample_key}' missing host_reference in config.")


def build_step_cmd(
    step: str,
    sample_key: str,
    sample_cfg: dict[str, Any],
    outdir: Path,
    threads: int,
    base_dir: Path,
) -> list[str]:
    """Build the command-line arguments for a specific step."""
    script = str(base_dir / STEP_SCRIPTS[step])
    sname = sample_key

    # Helper to resolve paths relative to base_dir
    def rp(p: str) -> str:
        return str(base_dir / p) if not Path(p).is_absolute() else p

    # Output path helpers for inter-step dependencies
    s01 = outdir / sname / "s01_qc"
    s02 = outdir / sname / "s02_construct_map"
    s03 = outdir / sname / "s03_extract"
    s04 = outdir / sname / "s04_assembly"
    s05 = outdir / sname / "s05_contig_map"
    s06 = outdir / sname / "s06_junction"
    s07 = outdir / sname / "s07_host_map"

    reads = sample_cfg.get("reads", {})
    construct_ref = rp(sample_cfg["construct_reference"])
    host_ref = rp(sample_cfg["host_reference"])

    if step == "1":
        return [sys.executable, script,
                "--r1", rp(reads["r1"]),
                "--r2", rp(reads["r2"]),
                "--outdir", str(outdir),
                "--threads", str(threads),
                "--sample-name", sname]
    elif step == "2":
        return [sys.executable, script,
                "--r1", str(s01 / f"{sname}_R1.fq.gz"),
                "--r2", str(s01 / f"{sname}_R2.fq.gz"),
                "--construct-ref", construct_ref,
                "--outdir", str(outdir),
                "--threads", str(threads),
                "--sample-name", sname]
    elif step == "3":
        return [sys.executable, script,
                "--bam", str(s02 / f"{sname}_construct.bam"),
                "--outdir", str(outdir),
                "--sample-name", sname,
                "--threads", str(threads)]
    elif step == "4":
        return [sys.executable, script,
                "--r1", str(s03 / f"{sname}_construct_R1.fq.gz"),
                "--r2", str(s03 / f"{sname}_construct_R2.fq.gz"),
                "--outdir", str(outdir),
                "--threads", str(threads),
                "--sample-name", sname]
    elif step == "5":
        return [sys.executable, script,
                "--contigs", str(s04 / "contigs.fasta"),
                "--host-ref", host_ref,
                "--construct-ref", construct_ref,
                "--outdir", str(outdir),
                "--threads", str(threads),
                "--sample-name", sname]
    elif step == "6":
        cmd = [sys.executable, script,
                "--host-paf", str(s05 / f"{sname}_contigs_to_host.paf"),
                "--construct-paf", str(s05 / f"{sname}_contigs_to_construct.paf"),
                "--contigs", str(s04 / "contigs.fasta"),
                "--outdir", str(outdir),
                "--sample-name", sname]
        # Lower identity threshold for element_db (fragmented references)
        if "element_db" in construct_ref or "_combined_db" in construct_ref:
            cmd.extend(["--min-identity", "0.70"])
        return cmd
    elif step == "7":
        return [sys.executable, script,
                "--r1", str(s01 / f"{sname}_R1.fq.gz"),
                "--r2", str(s01 / f"{sname}_R2.fq.gz"),
                "--host-ref", host_ref,
                "--outdir", str(outdir),
                "--threads", str(threads),
                "--sample-name", sname]
    elif step == "8":
        # Step 8 needs a WT BAM for comparison - look for it in config
        wt_sample = sample_cfg.get("wt_control")
        if wt_sample is None:
            # Try to auto-detect WT sample from config
            log.warning("No wt_control specified for %s; step 8 may fail", sname)
            wt_bam = Path("/dev/null")  # placeholder
        else:
            wt_bam = outdir / wt_sample / "s07_host_map" / f"{wt_sample}_host.bam"

        cmd = [sys.executable, script,
               "--treatment-bam", str(s07 / f"{sname}_host.bam"),
               "--wt-bam", str(wt_bam),
               "--host-ref", host_ref,
               "--outdir", str(outdir),
               "--sample-name", sname]

        # Add gRNA if specified in config
        grna = sample_cfg.get("grna")
        if grna:
            cmd.extend(["--grna", grna])

        # Add junctions if available
        junctions = s06 / "junctions.tsv"
        if junctions.exists():
            cmd.extend(["--junctions", str(junctions)])

        return cmd
    elif step == "10":
        return [sys.executable, script,
                "--construct-bam", str(s02 / f"{sname}_construct.bam"),
                "--host-bam", str(s07 / f"{sname}_host.bam"),
                "--construct-ref", construct_ref,
                "--host-ref", host_ref,
                "--outdir", str(outdir),
                "--sample-name", sname,
                "--junctions", str(s06 / "junctions.tsv")]
    elif step == "11":
        return [sys.executable, script,
                "--outdir", str(outdir),
                "--sample-name", sname,
                "--title", f"RedGene Report — {sname}"]
    else:
        sys.exit(f"ERROR: No command builder for step {step}")


def run_step(
    step: str,
    sample_key: str,
    sample_cfg: dict[str, Any],
    outdir: Path,
    threads: int,
    base_dir: Path,
    dry_run: bool,
) -> None:
    """Execute a single pipeline step for one sample."""
    script = base_dir / STEP_SCRIPTS[step]
    label = STEP_NAMES.get(step, step)

    if not script.exists() and not dry_run:
        sys.exit(f"ERROR: Script not found: {script}")

    cmd = build_step_cmd(step, sample_key, sample_cfg, outdir, threads, base_dir)

    log.info("Step %s: %s  [%s]", step, label, sample_key)
    if dry_run:
        log.info("  DRY-RUN: %s", " ".join(cmd))
        return

    log.info("  CMD: %s", " ".join(cmd))
    result = subprocess.run(cmd, cwd=str(base_dir))
    if result.returncode != 0:
        sys.exit(
            f"FAILED: Step {step} ({label}) for sample '{sample_key}' "
            f"exited with code {result.returncode}"
        )
    log.info("  Step %s completed successfully.", step)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="GMO Positive Control Characterization Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python run_pipeline.py                          # all samples, steps 1-6\n"
            "  python run_pipeline.py --sample rice_G281       # one sample\n"
            "  python run_pipeline.py --steps 1-6,10           # core + copy number\n"
            "  python run_pipeline.py --dry-run                # preview commands\n"
        ),
    )
    parser.add_argument(
        "--config", type=Path, default=Path("config.yaml"),
        help="Path to pipeline config YAML (default: config.yaml)",
    )
    parser.add_argument(
        "--sample", type=str, default=None,
        help="Sample key to process (default: all samples in config)",
    )
    parser.add_argument(
        "--steps", type=str, default="1-6",
        help="Steps to run, e.g. '1-6', '1-6,7,10' (default: 1-6)",
    )
    parser.add_argument(
        "--threads", type=int, default=None,
        help="Number of threads (default: from config or 8)",
    )
    parser.add_argument(
        "--outdir", type=Path, default=None,
        help="Output directory (default: from config or 'results')",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print commands without executing them",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    # Resolve base directory (where run_pipeline.py lives)
    base_dir = Path(__file__).resolve().parent

    # Resolve config path (relative paths are relative to base_dir)
    config_path = args.config if args.config.is_absolute() else base_dir / args.config
    cfg = load_config(config_path)

    # Determine output directory
    outdir = args.outdir
    if outdir is None:
        outdir = Path(cfg.get("output_dir", "results"))
    if not outdir.is_absolute():
        outdir = base_dir / outdir

    # Determine threads
    threads = args.threads
    if threads is None:
        threads = cfg.get("pipeline", {}).get("threads", 8)

    # Parse steps and samples
    steps = parse_steps(args.steps)
    samples = resolve_samples(cfg, args.sample)

    log.info("Pipeline configuration:")
    log.info("  Config:  %s", config_path)
    log.info("  Outdir:  %s", outdir)
    log.info("  Threads: %d", threads)
    log.info("  Steps:   %s", ", ".join(steps))
    log.info("  Samples: %s", ", ".join(samples))
    if args.dry_run:
        log.info("  Mode:    DRY-RUN")
    log.info("")

    # Validate inputs
    for sample_key in samples:
        validate_sample_inputs(cfg, sample_key, base_dir)

    # Run
    for sample_key in samples:
        log.info("=== Processing sample: %s ===", sample_key)
        if not args.dry_run:
            (outdir / sample_key).mkdir(parents=True, exist_ok=True)

        sample_cfg = cfg["samples"][sample_key]
        for step in steps:
            run_step(
                step=step,
                sample_key=sample_key,
                sample_cfg=sample_cfg,
                outdir=outdir,
                threads=threads,
                base_dir=base_dir,
                dry_run=args.dry_run,
            )

        log.info("=== Sample %s: all steps completed ===\n", sample_key)

    log.info("Pipeline finished successfully.")


if __name__ == "__main__":
    main()
