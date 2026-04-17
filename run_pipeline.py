#!/usr/bin/env python3
"""GMO Positive Control Characterization Pipeline - Main Entry Point.

Runs the pipeline steps sequentially for each sample defined in config.yaml.
Each step invokes the corresponding script in scripts/ via subprocess.

Pipeline (7 steps + optional step 4b):
  1. QC + trim (fastp)
  2. Map reads to construct + UniVec (bwa mem)
  3. Extract construct-hitting reads + mates
  4. Map all reads to host genome (bwa mem)  [bottleneck: 5-7h]
  4b. De novo assemble construct-hitting reads (s03 output) with SPAdes,
      producing a per-sample FASTA passed to step 5 via --extra-element-db.
  5. Targeted insert assembly + FP filtering  [CORE STEP]
  6. CRISPR indel detection (optional, needs WT control)
  7. Copy number estimation
"""

from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import yaml

from scripts._write_audit_header import write_audit_header

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STEP_SCRIPTS: dict[str, str] = {
    "1": "scripts/s01_qc.py",
    "2": "scripts/s02_construct_map.py",
    "3": "scripts/s03_extract_reads.py",
    "4": "scripts/s04_host_map.py",
    "4b": "scripts/s04b_construct_assembly.py",
    "5": "scripts/s05_insert_assembly.py",
    "6": "scripts/s06_indel.py",
    "7": "scripts/s07_copynumber.py",
}

STEP_NAMES: dict[str, str] = {
    "1": "QC + trim (fastp)",
    "2": "Map reads to construct + UniVec (bwa mem)",
    "3": "Extract construct-hitting reads + mates",
    "4": "Map all reads to host (bwa mem)",
    "4b": "De novo construct assembly (SPAdes, per-sample DB)",
    "5": "Targeted insert assembly + FP filtering",
    "6": "CRISPR indel detection (treatment vs WT)",
    "7": "Copy number estimation",
}

# Canonical execution order for steps. Used to expand ranges and to sort
# parsed steps. `4b` slots between `4` and `5` so `--steps 1-5` includes it.
STEP_ORDER: list[str] = ["1", "2", "3", "4", "4b", "5", "6", "7"]
STEP_INDEX: dict[str, int] = {s: i for i, s in enumerate(STEP_ORDER)}

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
    """Parse a step specification like '1-5', '1,2,5', '1-5,7', '4b' into a
    list of step keys in canonical execution order.

    Ranges are interpreted against STEP_ORDER, so '1-5' expands to
    ['1', '2', '3', '4', '4b', '5']. Only steps present in STEP_SCRIPTS are
    accepted. Individual tokens may be '4b' (or any other non-numeric step
    key) as well as plain integers.
    """
    steps: set[str] = set()
    for part in step_spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-", 1)
            lo = lo.strip()
            hi = hi.strip()
            if lo not in STEP_INDEX:
                sys.exit(f"ERROR: Unknown range start '{lo}' in --steps")
            if hi not in STEP_INDEX:
                sys.exit(f"ERROR: Unknown range end '{hi}' in --steps")
            lo_idx, hi_idx = STEP_INDEX[lo], STEP_INDEX[hi]
            if lo_idx > hi_idx:
                sys.exit(f"ERROR: Invalid range '{part}' (start > end)")
            for i in range(lo_idx, hi_idx + 1):
                steps.add(STEP_ORDER[i])
        else:
            steps.add(part)

    # Validate
    unknown = steps - set(STEP_SCRIPTS)
    if unknown:
        sys.exit(f"ERROR: Unknown step(s): {', '.join(sorted(unknown))}. "
                 f"Available: {', '.join(STEP_ORDER)}")

    # Return in canonical execution order
    return sorted(steps, key=lambda s: STEP_INDEX[s])


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
    no_remote_blast: bool = False,
    cfg: dict[str, Any] | None = None,
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
    s04 = outdir / sname / "s04_host_map"
    s04b = outdir / sname / "s04b_construct_asm"
    s05 = outdir / sname / "s05_insert_assembly"

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
                "--r1", str(s01 / f"{sname}_R1.fq.gz"),
                "--r2", str(s01 / f"{sname}_R2.fq.gz"),
                "--host-ref", host_ref,
                "--outdir", str(outdir),
                "--threads", str(threads),
                "--sample-name", sname]
    elif step == "4b":
        # De novo assemble construct-hitting reads (from s03) to build a
        # per-sample construct reference that step 5 can pass as
        # --extra-element-db. SPAdes gets ~2 GB/thread (half the SLURM
        # allocation of 4 GB/thread), leaving the rest as headroom for
        # SPAdes' graph-building memory overhead.
        memory_gb = max(8, threads * 4 // 2)
        return [sys.executable, script,
                "--r1", str(s03 / f"{sname}_construct_R1.fq.gz"),
                "--r2", str(s03 / f"{sname}_construct_R2.fq.gz"),
                "--outdir", str(outdir),
                "--sample-name", sname,
                "--threads", str(threads),
                "--memory-gb", str(memory_gb)]
    elif step == "5":
        # Check for WT-filtered reads (s03b output), fall back to s03
        s03_r1 = s03 / f"{sname}_filtered_R1.fq.gz"
        s03_r2 = s03 / f"{sname}_filtered_R2.fq.gz"
        if not s03_r1.exists():
            s03_r1 = s03 / f"{sname}_construct_R1.fq.gz"
            s03_r2 = s03 / f"{sname}_construct_R2.fq.gz"
        cmd = [sys.executable, script,
               "--host-bam", str(s04 / f"{sname}_host.bam"),
               "--host-ref", host_ref,
               "--element-db", construct_ref,
               "--construct-ref", construct_ref,
               "--s03-r1", str(s03_r1),
               "--s03-r2", str(s03_r2),
               "--outdir", str(outdir),
               "--sample-name", sname,
               "--threads", str(threads)]
        if no_remote_blast:
            cmd.append("--no-remote-blast")
        # If s04b produced a non-empty per-sample construct assembly, pass
        # it through as an extra element DB. Optional: samples that never
        # ran 4b simply skip this flag and run exactly as before.
        extra_db = s04b / "contigs.fasta"
        if extra_db.exists() and extra_db.stat().st_size > 0:
            cmd.extend(["--extra-element-db", str(extra_db)])
        # Always-on shared transgene reference DB (T5: gmo_combined_db_v2.fa -
        # cd-hit-est @ 0.95 dedup of common_payload + element_db + cas9_sgrna +
        # euginius_missing + payload_cds, with 4-way |src= tags for tier-based
        # merge in s05 classify_site_tiers). Path comes from
        # cfg.pipeline.common_payload_db, defaulting to the v2 tagged DB.
        cpd_rel = (cfg or {}).get("pipeline", {}).get(
            "common_payload_db", "element_db/gmo_combined_db_v2.fa"
        )
        cpd_path = base_dir / cpd_rel if not Path(cpd_rel).is_absolute() else Path(cpd_rel)
        if cpd_path.exists() and cpd_path.stat().st_size > 0:
            cmd.extend(["--common-payload-db", str(cpd_path)])
        # T10: auto-inject --mask-bed (host-endogenous ortholog regions from T9)
        # matched by host_reference filename stem. Missing BED → no-op.
        host_key = Path(host_ref).stem.lower()
        bed_map = {
            "osativa_323_v7.0": "rice_osativa_v7.bed",
            "slm_r2.0.pmol": "tomato_slm_r2.bed",
            "cucsat_b10v3": "cucumber_b10v3.bed",
            "zm_b73_v5": "corn_zm_b73_v5.bed",
            "gmax_v4.0": "soybean_gmax_v4.bed",
        }
        bed_name = next((v for k, v in bed_map.items() if k in host_key), None)
        if bed_name:
            bed_path = base_dir / "docs" / "host_masks" / bed_name
            if bed_path.exists():
                cmd.extend(["--mask-bed", str(bed_path)])
        return cmd
    elif step == "6":
        # Step 6 needs a WT BAM for comparison
        wt_sample = sample_cfg.get("wt_control")
        if wt_sample is None:
            log.warning("No wt_control specified for %s; step 6 may fail", sname)
            wt_bam = Path("/dev/null")
        else:
            wt_bam = outdir / wt_sample / "s04_host_map" / f"{wt_sample}_host.bam"

        cmd = [sys.executable, script,
               "--treatment-bam", str(s04 / f"{sname}_host.bam"),
               "--wt-bam", str(wt_bam),
               "--host-ref", host_ref,
               "--outdir", str(outdir),
               "--sample-name", sname]

        # Add gRNA if specified in config
        grna = sample_cfg.get("grna")
        if grna:
            cmd.extend(["--grna", grna])

        return cmd
    elif step == "7":
        # Copy number: junctions optional (use s05 site count if available)
        cmd = [sys.executable, script,
                "--construct-bam", str(s02 / f"{sname}_construct.bam"),
                "--host-bam", str(s04 / f"{sname}_host.bam"),
                "--construct-ref", construct_ref,
                "--host-ref", host_ref,
                "--outdir", str(outdir),
                "--sample-name", sname]
        # Pass s05 stats for cross-validation if available
        s05_stats = s05 / "s05_stats.txt"
        if s05_stats.exists():
            cmd.extend(["--site-stats", str(s05_stats)])
        return cmd
    else:
        sys.exit(f"ERROR: No command builder for step {step}")


def _fanout_step5(
    cmd: list[str],
    sample_key: str,
    outdir: Path,
    threads: int,
    base_dir: Path,
    dry_run: bool,
) -> None:
    """T8: run s05 Phase 1+1.5 inline, then submit SLURM array (Phase 2_3)
    and Phase 4 afterok job via scripts/submit_s05_array.sh.

    `cmd` is the full s05 command built by build_step_cmd("5", ...). We
    reuse it for the inline Phase 1_1.5 call (adding `--phase 1_1.5`) and
    parse out the path args we need to re-emit as env vars for the array
    wrapper (which builds its own s05 invocations for Phase 2_3 / 4).
    """
    # Parse out key=value-ish pairs from the s05 CLI so we can re-thread
    # them through to the array wrapper as env vars.  cmd layout is:
    #   [python, scripts/s05_insert_assembly.py, --host-bam, X, --host-ref, Y, ...]
    flag_map = {
        "--host-bam": "S05_HOST_BAM",
        "--host-ref": "S05_HOST_REF",
        "--element-db": "S05_ELEMENT_DB",
        "--construct-ref": "S05_CONSTRUCT_REF",
        "--extra-element-db": "S05_EXTRA_ELEMENT_DB",
        "--common-payload-db": "S05_COMMON_PAYLOAD_DB",
        "--s03-r1": "S05_S03_R1",
        "--s03-r2": "S05_S03_R2",
    }
    env_overrides: dict[str, str] = {}
    i = 2  # skip [python, script.py]
    has_no_remote_blast = False
    while i < len(cmd):
        tok = cmd[i]
        if tok == "--no-remote-blast":
            has_no_remote_blast = True
            i += 1
            continue
        if tok in flag_map and i + 1 < len(cmd):
            env_overrides[flag_map[tok]] = cmd[i + 1]
            i += 2
            continue
        i += 1
    if has_no_remote_blast:
        env_overrides["S05_NO_REMOTE_BLAST"] = "1"

    # Phase 1_1.5 inline: reuse exact same cmd, just append --phase.
    phase_cmd = list(cmd) + ["--phase", "1_1.5"]
    log.info("  [T8] Phase 1+1.5 inline: %s", " ".join(phase_cmd))
    if not dry_run:
        result = subprocess.run(phase_cmd, cwd=str(base_dir))
        if result.returncode != 0:
            sys.exit(
                f"FAILED: T8 Phase 1+1.5 for '{sample_key}' "
                f"exited with code {result.returncode}"
            )

    step_dir = outdir / sample_key / "s05_insert_assembly"
    array_cmd = [
        "bash",
        str(base_dir / "scripts" / "submit_s05_array.sh"),
        sample_key,
        str(step_dir),
        str(threads),
        str(outdir),
    ]
    env_str = " ".join(f"{k}={v}" for k, v in env_overrides.items())
    log.info("  [T8] Array submit: %s %s", env_str, " ".join(array_cmd))
    if dry_run:
        log.info("  [T8] (dry-run: sbatch calls suppressed)")
        return

    env = dict(os.environ)
    env.update(env_overrides)
    result = subprocess.run(array_cmd, cwd=str(base_dir), env=env)
    if result.returncode != 0:
        sys.exit(
            f"FAILED: T8 submit_s05_array.sh for '{sample_key}' "
            f"exited with code {result.returncode}"
        )
    log.info("  [T8] sbatch submissions complete; exiting without waiting.")


def run_step(
    step: str,
    sample_key: str,
    sample_cfg: dict[str, Any],
    outdir: Path,
    threads: int,
    base_dir: Path,
    dry_run: bool,
    no_remote_blast: bool = False,
    cfg: dict[str, Any] | None = None,
    fanout: bool = False,
) -> None:
    """Execute a single pipeline step for one sample."""
    script = base_dir / STEP_SCRIPTS[step]
    label = STEP_NAMES.get(step, step)

    if not script.exists() and not dry_run:
        sys.exit(f"ERROR: Script not found: {script}")

    cmd = build_step_cmd(step, sample_key, sample_cfg, outdir, threads, base_dir,
                         no_remote_blast=no_remote_blast, cfg=cfg)

    log.info("Step %s: %s  [%s]", step, label, sample_key)

    # T8 임시안: step 5 fan-out path.  Only applies to step 5.  v1.1 replaces
    # this with the full module split.
    if step == "5" and fanout:
        _fanout_step5(cmd, sample_key, outdir, threads, base_dir, dry_run)
        return

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
            "  python run_pipeline.py                          # all samples, steps 1-5 (includes 4b)\n"
            "  python run_pipeline.py --sample rice_G281       # one sample\n"
            "  python run_pipeline.py --steps 1-5,7            # core + copy number\n"
            "  python run_pipeline.py --steps 5                # re-run insert assembly\n"
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
        "--steps", type=str, default="1-5",
        help=(
            "Steps to run (1, 2, 3, 4, 4b, 5, 6, 7). "
            "Ranges expand across step 4b, so '1-5' includes 4b "
            "(de novo construct assembly). (default: 1-5)"
        ),
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
    parser.add_argument(
        "--no-remote-blast", action="store_true",
        help="Skip remote NCBI nt BLAST in step 5 (use local element_db only)",
    )
    parser.add_argument(
        "--fanout", action="store_true",
        help=(
            "T8 임시안 (v1.0 MVP): for step 5 only, run Phase 1+1.5 inline, "
            "then submit a SLURM array (one task per positive site) for "
            "Phase 2+3, and a Phase 4 afterok-dependency job. Exits after "
            "sbatch submission (does not wait for the array to finish). "
            "UGT72E3 path: ~48 h sequential → ~4 h fan-out."
        ),
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

        # AC-6 regulatory audit header (R-1..R-4): input SHA-256,
        # pipeline git commit + dirty flag, element-DB manifest, and
        # software version manifest. Written once per sample run.
        if not args.dry_run:
            try:
                reads_cfg = sample_cfg.get("reads", {})
                r1_path = Path(reads_cfg["r1"])
                r2_path = Path(reads_cfg["r2"])
                if not r1_path.is_absolute():
                    r1_path = base_dir / r1_path
                if not r2_path.is_absolute():
                    r2_path = base_dir / r2_path
                write_audit_header(
                    sample=sample_key,
                    reads_r1=r1_path,
                    reads_r2=r2_path,
                    db_manifest=base_dir / "element_db" / "gmo_combined_db_manifest.tsv",
                    out_path=outdir / sample_key / "audit_header.json",
                )
            except Exception as exc:  # pragma: no cover - best-effort audit
                log.warning("audit_header write failed for %s: %s",
                            sample_key, exc)

        for step in steps:
            run_step(
                step=step,
                sample_key=sample_key,
                sample_cfg=sample_cfg,
                outdir=outdir,
                threads=threads,
                base_dir=base_dir,
                dry_run=args.dry_run,
                no_remote_blast=args.no_remote_blast,
                cfg=cfg,
                fanout=args.fanout,
            )

        log.info("=== Sample %s: all steps completed ===\n", sample_key)

    log.info("Pipeline finished successfully.")


if __name__ == "__main__":
    main()
