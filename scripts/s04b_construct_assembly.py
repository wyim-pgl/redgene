#!/usr/bin/env python3
"""Step 4b: De novo assemble construct-hitting reads (from s03) with SPAdes.

The resulting contigs.fasta serves as a sample-specific construct reference
passed into step 5 via --extra-element-db. This lets us annotate inserts
whose transgene payload (e.g. AtYUCCA6, bar gene) isn't present in the
shared element_db.

Input:  s03 R1/R2 FASTQ.GZ (construct-hitting reads + their mates)
Output: results/<sample>/s04b_construct_asm/contigs.fasta
        results/<sample>/s04b_construct_asm/spades_stderr.log

Returns exit code 0 in all paths: empty input, SPAdes success (contigs
copied), SPAdes failure on small input (empty contigs.fasta written).
This means downstream steps can always find `contigs.fasta` at the
expected path, even if assembly was infeasible.
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path


def _is_empty_fastq(path: Path) -> bool:
    """Empty-gzip or zero-read FASTQ detection. Missing/unreadable = empty."""
    try:
        with gzip.open(path, "rt") as fh:
            return fh.read(1) == ""
    except OSError:
        try:
            return path.stat().st_size == 0
        except OSError:
            return True  # missing or permission-denied → treat as empty


def run(args: argparse.Namespace) -> int:
    step_dir = args.outdir / args.sample_name / "s04b_construct_asm"
    step_dir.mkdir(parents=True, exist_ok=True)
    contigs = step_dir / "contigs.fasta"

    if _is_empty_fastq(args.r1) or _is_empty_fastq(args.r2):
        print("[s04b] Empty inputs; writing empty contigs.fasta", file=sys.stderr)
        contigs.write_text("")
        return 0

    spades_out = step_dir / "_spades_run"
    if spades_out.exists():
        shutil.rmtree(spades_out)
    cmd = [
        "spades.py", "--only-assembler", "--careful",
        "-1", str(args.r1), "-2", str(args.r2),
        "-o", str(spades_out),
        "-t", str(args.threads),
        "-m", str(args.memory_gb),
    ]
    print(f"[s04b] {' '.join(cmd)}", file=sys.stderr)
    log = step_dir / "spades_stderr.log"
    with log.open("w") as fh_err:
        proc = subprocess.run(cmd, stderr=fh_err)
    if proc.returncode != 0:
        print(
            f"[s04b] SPAdes exited {proc.returncode} (likely too few reads for "
            f"coverage); writing empty contigs.fasta. See {log}",
            file=sys.stderr,
        )
        contigs.write_text("")
        shutil.rmtree(spades_out, ignore_errors=True)
        return 0

    spades_contigs = spades_out / "contigs.fasta"
    if spades_contigs.exists() and spades_contigs.stat().st_size > 0:
        shutil.copy2(spades_contigs, contigs)
    else:
        contigs.write_text("")

    shutil.rmtree(spades_out, ignore_errors=True)
    n_contigs = sum(1 for line in contigs.read_text().splitlines()
                    if line.startswith(">"))
    print(f"[s04b] Wrote {contigs} ({n_contigs} contigs)", file=sys.stderr)
    return 0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--r1", type=Path, required=True)
    p.add_argument("--r2", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--memory-gb", type=int, default=16)
    return p.parse_args()


if __name__ == "__main__":
    sys.exit(run(parse_args()))
