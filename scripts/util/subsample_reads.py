#!/usr/bin/env python3
"""Fraction-first paired-end FASTQ subsampler for the AC-7 coverage sweep.

This is a thin wrapper around ``seqkit sample`` designed for the Issue #2
coverage-sensitivity batch. Compared to ``scripts/subsample_reads.py`` (the
original coverage-based subsampler, kept in place for backwards compat), this
helper:

* accepts ``--fraction`` directly (preferred, no genome-size math required)
* can *also* compute a fraction from ``--target-coverage`` + genome-size +
  total-pair-count when the operator wants a coverage-labeled output
* exposes a small Python API (``compute_fraction``, ``subsample``,
  ``layout_coverage_dir``) so tests and the sweep runner can import it

Usage (fraction mode):
    python scripts/util/subsample_reads.py \\
        --r1 in_R1.fq.gz --r2 in_R2.fq.gz \\
        --fraction 0.2 --seed 42 --output-prefix out_20pct

Usage (coverage-target mode):
    python scripts/util/subsample_reads.py \\
        --r1 in_R1.fq.gz --r2 in_R2.fq.gz \\
        --target-coverage 10 --genome-size 374000000 \\
        --total-pairs 50000000 --read-length 150 \\
        --output-prefix out_10x
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def compute_fraction(
    *,
    fraction: float | None,
    target_coverage: float | None,
    total_pairs: int | None,
    read_length: int,
    genome_size: int | None,
) -> float:
    """Return the seqkit ``--proportion`` value, clamped to ``(0, 1]``.

    If ``fraction`` is given, it is returned verbatim (after clamping).
    Otherwise ``(target_coverage, total_pairs, genome_size)`` must all be
    provided and the fraction is computed as
    ``target_coverage * genome_size / (2 * read_length * total_pairs)``.
    """
    if fraction is not None:
        if fraction <= 0:
            raise ValueError(f"fraction must be > 0, got {fraction}")
        return min(1.0, float(fraction))
    if target_coverage is None or total_pairs is None or genome_size is None:
        raise ValueError(
            "compute_fraction requires either --fraction or "
            "(--target-coverage + --total-pairs + --genome-size)"
        )
    pairs_needed = target_coverage * genome_size / (2 * read_length)
    return min(1.0, pairs_needed / total_pairs)


def subsample(
    *,
    r1: Path,
    r2: Path,
    fraction: float,
    seed: int,
    output_prefix: Path,
) -> int:
    """Run seqkit sample on both mates with a shared seed.

    Returns the number of pairs in the first output (best-effort — skips
    counting if seqkit stats fails, returns -1).
    """
    out_r1 = Path(f"{output_prefix}_R1.fq.gz")
    out_r2 = Path(f"{output_prefix}_R2.fq.gz")
    out_r1.parent.mkdir(parents=True, exist_ok=True)

    for in_path, out_path in ((r1, out_r1), (r2, out_r2)):
        subprocess.run(
            [
                "seqkit", "sample",
                "-p", str(fraction),
                "-s", str(seed),
                "-o", str(out_path),
                str(in_path),
            ],
            check=True,
        )

    # Count pairs in output (via seqkit stats TSV) — tolerant to format drift.
    try:
        result = subprocess.run(
            ["seqkit", "stats", "-T", str(out_r1)],
            capture_output=True, text=True, check=True,
        )
        lines = result.stdout.strip().splitlines()
        if len(lines) >= 2:
            return int(lines[1].split("\t")[3].replace(",", ""))
    except (subprocess.CalledProcessError, ValueError, IndexError):
        pass
    return -1


def layout_coverage_dir(
    *,
    base_outdir: Path,
    sample: str,
    coverage_tag: str,
    r1: Path,
    r2: Path,
) -> Path:
    """Create ``<base_outdir>/<sample>_cov<tag>/`` and symlink subsampled reads.

    The resulting directory is consumed by ``run_pipeline.py --sample
    <sample>_cov<tag>`` via its input-path override. We use symlinks (not
    copies) so disk usage stays minimal across the 3 hosts x 4 coverage grid.
    """
    base_outdir = Path(base_outdir)
    cov_dir = base_outdir / f"{sample}_cov{coverage_tag}"
    cov_dir.mkdir(parents=True, exist_ok=True)

    for in_path, tag in ((r1, "R1"), (r2, "R2")):
        link = cov_dir / f"{sample}_cov{coverage_tag}_{tag}.fq.gz"
        if link.exists() or link.is_symlink():
            link.unlink()
        try:
            link.symlink_to(Path(in_path).resolve())
        except OSError:
            # Fallback to copy when symlinks are disallowed (rare on NFS/GPFS)
            import shutil
            shutil.copy2(in_path, link)
    return cov_dir


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--r1", type=Path, required=True)
    p.add_argument("--r2", type=Path, required=True)
    p.add_argument("--fraction", type=float, default=None,
                   help="Sampling fraction (0-1). Overrides --target-coverage.")
    p.add_argument("--target-coverage", type=float, default=None,
                   help="Optional coverage-target mode. Requires --genome-size, "
                        "--total-pairs, --read-length.")
    p.add_argument("--genome-size", type=int, default=None)
    p.add_argument("--total-pairs", type=int, default=None)
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output-prefix", type=Path, required=True,
                   help="Writes {prefix}_R1.fq.gz, {prefix}_R2.fq.gz")
    args = p.parse_args(argv)

    frac = compute_fraction(
        fraction=args.fraction,
        target_coverage=args.target_coverage,
        total_pairs=args.total_pairs,
        read_length=args.read_length,
        genome_size=args.genome_size,
    )
    print(f"[subsample] fraction={frac:.4f} seed={args.seed}", file=sys.stderr)
    n = subsample(
        r1=args.r1, r2=args.r2,
        fraction=frac, seed=args.seed,
        output_prefix=args.output_prefix,
    )
    print(f"[subsample] wrote {n} pairs → {args.output_prefix}_R[12].fq.gz",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
