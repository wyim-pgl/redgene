#!/usr/bin/env python3
"""Step 4b: De novo assemble construct-hitting reads (from s03) with SPAdes.

The resulting contigs.fasta serves as a sample-specific construct reference
passed into step 5 via --extra-element-db. This lets us annotate inserts
whose transgene payload (e.g. AtYUCCA6, bar gene) isn't present in the
shared element_db.

To curb cross-reactivity (Task 8 showed passing the raw 1,345-contig SPAdes
output inflated Phase 1.5 transgene-positive sites 25x), contigs are
filtered against curated marker DBs (common_payload.fa +
gmo_combined_db.fa): only contigs with at least one >=90% identity /
>=200 bp blastn hit survive. The unfiltered assembly is preserved as
``contigs_all.fasta`` for debugging; ``--no-filter`` restores legacy
behaviour.

Input:  s03 R1/R2 FASTQ.GZ (construct-hitting reads + their mates)
Output: results/<sample>/s04b_construct_asm/contigs.fasta       (filtered)
        results/<sample>/s04b_construct_asm/contigs_all.fasta   (raw SPAdes, when filtering is on)
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
from typing import Iterable, List, Sequence, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_MARKER_DBS = [
    REPO_ROOT / "element_db" / "common_payload.fa",
    REPO_ROOT / "element_db" / "gmo_combined_db.fa",
]


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


def _iter_fasta_records(path: Path) -> Iterable[Tuple[str, List[str]]]:
    """Yield (header_line, [sequence_lines]) tuples in input order."""
    header: str | None = None
    seq_lines: List[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_lines
                header = line.rstrip("\n")
                seq_lines = []
            else:
                seq_lines.append(line.rstrip("\n"))
    if header is not None:
        yield header, seq_lines


def _fasta_id(header: str) -> str:
    """Strip leading '>' and trailing comment; blastn qseqid semantics."""
    body = header[1:] if header.startswith(">") else header
    return body.split()[0] if body else body


def _is_nonempty_file(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def _filter_contigs_by_markers(
    contigs_fa: Path,
    marker_dbs: Sequence[Path],
    out_fa: Path,
    min_identity: float = 90.0,
    min_aln_len: int = 200,
) -> Tuple[int, int]:
    """Keep contigs with >=1 blastn hit passing (pident>=min_identity AND
    length>=min_aln_len) against ANY of ``marker_dbs``.

    - ``contigs_fa`` may be empty; returns (0, 0).
    - Missing/empty marker DBs are skipped (not an error).
    - Writes passing contigs to ``out_fa`` in input order.

    Returns (kept_count, total_count).
    """
    total = 0
    records: List[Tuple[str, List[str]]] = []
    for header, seq_lines in _iter_fasta_records(contigs_fa):
        records.append((header, seq_lines))
        total += 1

    if total == 0:
        out_fa.write_text("")
        return 0, 0

    passing: set[str] = set()
    for db in marker_dbs:
        if not _is_nonempty_file(Path(db)):
            print(f"[s04b] filter: skipping missing/empty marker DB {db}",
                  file=sys.stderr)
            continue
        cmd = [
            "blastn",
            "-query", str(contigs_fa),
            "-subject", str(db),
            "-outfmt", "6 qseqid sseqid pident length",
            "-evalue", "1e-10",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print(
                f"[s04b] filter: blastn vs {db} failed (rc={proc.returncode}): "
                f"{proc.stderr.strip()}",
                file=sys.stderr,
            )
            continue
        for line in proc.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            qseqid, _sseqid, pident_s, length_s = parts[0], parts[1], parts[2], parts[3]
            try:
                pident = float(pident_s)
                length = int(length_s)
            except ValueError:
                continue
            if pident >= min_identity and length >= min_aln_len:
                passing.add(qseqid)

    kept = 0
    with out_fa.open("w") as out:
        for header, seq_lines in records:
            if _fasta_id(header) in passing:
                out.write(header + "\n")
                for line in seq_lines:
                    out.write(line + "\n")
                kept += 1
    return kept, total


def _resolve_marker_dbs(raw: Sequence[Path] | None) -> List[Path]:
    if raw:
        return [Path(p) for p in raw]
    return [Path(p) for p in DEFAULT_MARKER_DBS]


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
    print(f"[s04b] SPAdes produced {n_contigs} contigs", file=sys.stderr)

    if args.no_filter:
        print(f"[s04b] --no-filter: keeping raw contigs.fasta ({n_contigs} contigs)",
              file=sys.stderr)
        return 0

    if n_contigs == 0:
        # Nothing to filter; leave empty contigs.fasta in place.
        print("[s04b] filter: 0 contigs to filter; leaving empty contigs.fasta",
              file=sys.stderr)
        return 0

    marker_dbs = _resolve_marker_dbs(args.marker_db)
    contigs_all = step_dir / "contigs_all.fasta"
    # Move raw -> contigs_all.fasta, then write filtered -> contigs.fasta
    shutil.move(str(contigs), str(contigs_all))
    kept, total = _filter_contigs_by_markers(
        contigs_all, marker_dbs, contigs,
        min_identity=args.min_identity, min_aln_len=args.min_aln_len,
    )
    db_names = ", ".join(str(p) for p in marker_dbs)
    print(
        f"[s04b] filter kept {kept}/{total} contigs via marker DBs "
        f"({db_names}; pident>={args.min_identity}, len>={args.min_aln_len})",
        file=sys.stderr,
    )
    if kept == 0:
        print("[s04b] WARNING: 0 contigs survived marker-DB filtering; "
              "downstream extra_db will be empty", file=sys.stderr)
    print(f"[s04b] Wrote {contigs} ({kept} contigs)", file=sys.stderr)
    return 0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--r1", type=Path, required=True)
    p.add_argument("--r2", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--memory-gb", type=int, default=16)
    p.add_argument(
        "--marker-db", type=Path, action="append", default=None,
        help="FASTA to filter contigs against (repeatable). Default: "
             "element_db/common_payload.fa + element_db/gmo_combined_db.fa.",
    )
    p.add_argument(
        "--min-identity", type=float, default=90.0,
        help="Minimum blastn percent identity for a marker hit (default 90.0).",
    )
    p.add_argument(
        "--min-aln-len", type=int, default=200,
        help="Minimum blastn alignment length (bp) for a marker hit (default 200).",
    )
    p.add_argument(
        "--no-filter", action="store_true",
        help="Skip marker-DB filtering; keep raw SPAdes contigs (legacy).",
    )
    return p.parse_args()


if __name__ == "__main__":
    sys.exit(run(parse_args()))
