"""Filter A — host-fraction false-positive check (Issue #4, Session 2).

Extracted verbatim from ``scripts/s05_insert_assembly.py`` so the BLAST-driven
host-fraction measurement can be unit-tested without the full monolith.

Filter A rationale
------------------
A genuine T-DNA insert has a contiguous foreign (non-host) stretch that
contains the construct sequence. False positives caused by host-derived
construct elements (e.g., rice Ubi1 promoter copies) lack this gap and
align almost end-to-end against the host. Filter A measures:

    host_fraction         = host-aligned bp / (insert_len - N_count)
    largest_foreign_gap   = longest contiguous non-host stretch (bp)

The monolith and ``scripts/s05/verdict.py`` Rule 1/6 consume the returned
tuple via ``FilterEvidence``. The N-count subtraction protects against
gap-fill placeholders inflating the denominator.

This module is intentionally side-effect free aside from the BLAST
subprocess call. The on-disk BLAST output filename pattern
(``_<insert_stem>_vs_host_chrom.tsv``) is shared with
``filter_c_chimeric._check_chimeric_assembly`` so the second filter can
reuse the first filter's output without re-running megablast.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from .primitives import read_fasta


# ---------------------------------------------------------------------------
# Thresholds (mirrored from the monolith so behaviour is bit-identical).
# ---------------------------------------------------------------------------
INSERT_HOST_MIN_PIDENT = 90.0   # min identity for host alignment to count


# ---------------------------------------------------------------------------
# BLAST-driven host-fraction measurement
# ---------------------------------------------------------------------------

def _blast_insert_vs_host(
    insert_fasta: Path,
    host_ref: Path,
    workdir: Path,
    threads: int = 4,
) -> tuple[float, int, int, int]:
    """BLAST assembled insert against host genome to measure host-fraction.

    Returns (host_fraction, host_covered_bp, insert_length, largest_foreign_gap).
    host_fraction = fraction of insert positions covered by host alignments.
    largest_foreign_gap = longest contiguous stretch NOT covered by host.
    """
    seqs = read_fasta(insert_fasta)
    if not seqs:
        return 0.0, 0, 0, 0
    insert_seq = list(seqs.values())[0]
    insert_len = len(insert_seq)
    if insert_len == 0:
        return 0.0, 0, 0, 0

    # Use _chrom variant filename so _check_chimeric_assembly can reuse it
    blast_out = workdir / f"_{insert_fasta.stem}_vs_host_chrom.tsv"
    if not blast_out.exists():
        result = subprocess.run(
            ["blastn", "-task", "megablast",
             "-query", str(insert_fasta), "-db", str(host_ref),
             "-outfmt", "6 qseqid qstart qend sseqid pident length",
             "-evalue", "1e-10", "-max_target_seqs", "10",
             "-num_threads", str(threads),
             "-out", str(blast_out)],
            stderr=subprocess.DEVNULL,
        )
    if not blast_out.exists():
        return 0.0, 0, insert_len, insert_len

    # Merge overlapping host-aligned intervals to compute non-redundant coverage
    intervals: list[tuple[int, int]] = []
    with open(blast_out) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue
            pident = float(cols[4])
            if pident < INSERT_HOST_MIN_PIDENT:
                continue
            q_start = int(cols[1])
            q_end = int(cols[2])
            lo, hi = min(q_start, q_end), max(q_start, q_end)
            intervals.append((lo, hi))

    if not intervals:
        return 0.0, 0, insert_len, insert_len

    # Merge overlapping intervals
    intervals.sort()
    merged: list[tuple[int, int]] = [intervals[0]]
    for lo, hi in intervals[1:]:
        if lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))

    host_bp = sum(hi - lo + 1 for lo, hi in merged)

    # Compute largest gap between host-aligned regions (= foreign/T-DNA region)
    # Include gaps at start and end of insert
    gaps: list[int] = []
    gaps.append(merged[0][0] - 1)                    # gap before first host hit
    for i in range(1, len(merged)):
        gaps.append(merged[i][0] - merged[i - 1][1] - 1)
    gaps.append(insert_len - merged[-1][1])           # gap after last host hit
    largest_gap = max(gaps) if gaps else 0

    # Exclude N-runs from denominator (Ns are gap-fill placeholders, not real sequence)
    n_count = insert_seq.upper().count("N")
    effective_len = insert_len - n_count
    if effective_len <= 0:
        return 0.0, host_bp, insert_len, largest_gap

    return host_bp / effective_len, host_bp, insert_len, largest_gap
