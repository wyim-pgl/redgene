"""Filter D — alternative-locus / construct-host-coverage check (Issue #4, Session 2).

Extracted verbatim from ``scripts/s05_insert_assembly.py`` so the BLAST-driven
construct + host combined-coverage analysis can be unit-tested without the
full monolith.

Filter D rationale
------------------
Real T-DNA inserts have low construct coverage (~10%, only border regions
align to the construct because the bulk of the insert is foreign cargo).
False positives created by construct-element homology — e.g., a host
chromosomal copy of CaMV 35S promoter mis-assembled together with reads
from the actual construct — show high construct coverage (~50%, the
"insert" essentially **is** a construct fragment) AND high combined
construct+host coverage (~85%+). Filter D rejects sites that meet both
fractions simultaneously.

The "altlocus" name comes from the design spec: this filter catches
inserts that are really host genomic DNA at an alternative (non-T-DNA)
locus, mis-assembled because the host genome contains a paralog of a
construct element. See ``scripts/s05/verdict.py`` Rule 4.

This module is intentionally side-effect free aside from the BLAST
subprocess call (`-subject` mode, no on-disk DB). The returned
``(is_fp, construct_frac, host_frac, combined_frac)`` tuple is consumed
by the assemble_and_evaluate_site() path and ``compute_verdict``.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from .primitives import log


# ---------------------------------------------------------------------------
# Thresholds (mirrored from the monolith so behaviour is bit-identical).
# ---------------------------------------------------------------------------
CONSTRUCT_HOST_MIN_COMBINED = 0.85  # min combined coverage (construct + host)
CONSTRUCT_MIN_FRACTION = 0.25       # min construct coverage to suspect FP
CONSTRUCT_HOST_MIN_PIDENT = 80.0    # identity threshold for construct BLAST


# ---------------------------------------------------------------------------
# BLAST-driven construct+host combined coverage check
# ---------------------------------------------------------------------------

def _check_construct_host_coverage(
    insert_fasta: Path,
    element_db: Path,
    host_fraction: float,
    host_bp: int,
    insert_len: int,
    n_count: int,
    workdir: Path,
    threads: int = 4,
) -> tuple[bool, float, float, float]:
    """Check if insert is fully explained by construct + host coverage.

    Real T-DNA inserts have low construct coverage (~10%, only border regions).
    False positives from construct-element homology have high construct
    coverage (~50%, the insert IS the construct fragment).

    Returns (is_fp, construct_frac, host_frac, combined_frac).

    Failure mode
    ------------
    ``is_fp=False`` means "this site is NOT a construct+host false positive", so
    a crashed blastn that returned ``False`` silently promoted the site toward
    CANDIDATE.  A failure is now logged with its stderr, and the host evidence
    already measured by Filter A is passed through untouched instead of being
    zeroed — the site is reported as *unmeasured by Filter D*, not as *clean*.
    """
    eff_len = insert_len - n_count
    if eff_len <= 0:
        return False, 0.0, 0.0, 0.0

    # BLAST insert vs construct/element_db.
    # `threads` is accepted for signature parity with the other filters but is
    # deliberately NOT forwarded: this is `-subject` mode, and blastn 2.16.0
    # answers `-num_threads` there with "'num_threads' is currently ignored when
    # 'subject' is specified" — the flag would only add a warning per site.
    result = subprocess.run(
        ["blastn", "-task", "megablast",
         "-query", str(insert_fasta), "-subject", str(element_db),
         "-outfmt", "6 qstart qend pident",
         "-evalue", "1e-5"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if result.returncode != 0:
        stderr_tail = (result.stderr or "")[-400:]
        log(f"  Filter D: construct blastn failed (rc={result.returncode}); "
            f"construct coverage unmeasured. stderr: {stderr_tail}")
        return False, 0.0, host_fraction, host_fraction

    # Merge overlapping intervals at >= threshold identity
    intervals: list[tuple[int, int]] = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        cols = line.split("\t")
        qs, qe, pid = int(cols[0]), int(cols[1]), float(cols[2])
        if pid >= CONSTRUCT_HOST_MIN_PIDENT:
            intervals.append((min(qs, qe), max(qs, qe)))

    if not intervals:
        return False, 0.0, host_fraction, host_fraction

    intervals.sort()
    merged: list[tuple[int, int]] = []
    for s, e in intervals:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))

    construct_bp = sum(e - s + 1 for s, e in merged)
    construct_frac = construct_bp / eff_len
    combined_frac = min(1.0, (construct_bp + host_bp) / eff_len)

    is_fp = (construct_frac >= CONSTRUCT_MIN_FRACTION
             and combined_frac >= CONSTRUCT_HOST_MIN_COMBINED)

    return is_fp, construct_frac, host_fraction, combined_frac
