"""Filter B — construct-flanking false-positive check (Issue #4, Session 1).

Extracted verbatim from ``scripts/s05_insert_assembly.py`` so the BLAST-driven
flanking discovery and the per-site overlap test can be unit-tested without
needing the full monolith.

Filter B rationale
------------------
Some construct references include host genomic flanking DNA at their ends
(e.g., rice_G281 has 201 bp of Chr11 at positions 8758-8958 appended to the
T-DNA cassette). When reads map to the host at those coordinates, they look
like a genuine soft-clip junction, but the "insert" sequence is identical to
the construct's own flanking region — a false detection, not a T-DNA
insertion. See ``docs/superpowers/specs/2026-04-19-s05-module-split-design.md``
§2.3 and ``scripts/s05/verdict.py`` Rule 2.

This module is intentionally side-effect free aside from the BLAST subprocess
call inside ``_find_construct_flanking_regions``. The "pure filter-B check"
(``filter_b_flanking_hit``) derives the ``FilterEvidence.flanking_hit`` tuple
consumed by ``compute_verdict`` and contains no I/O.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from .primitives import log


# ---------------------------------------------------------------------------
# Thresholds (mirrored from the monolith so behaviour is bit-identical).
# ---------------------------------------------------------------------------
# If the construct reference contains host genomic DNA at its ends (common
# when constructs are cloned with flanking), sites at those host coordinates
# are false detections.
CONSTRUCT_FLANK_PIDENT = 95.0   # min identity for construct→host flanking hit
CONSTRUCT_FLANK_MIN_LEN = 50    # min alignment length
CONSTRUCT_FLANK_SLOP = 500      # bp slop when checking site overlap


# ---------------------------------------------------------------------------
# BLAST-driven flanking discovery
# ---------------------------------------------------------------------------

def _find_construct_flanking_regions(
    construct_ref: Path,
    host_ref: Path,
    workdir: Path,
    threads: int = 4,
) -> list[tuple[str, int, int]]:
    """BLAST construct reference vs host to find host-flanking regions.

    Some construct references include host genomic flanking DNA at their ends
    (e.g., rice_G281 has 201bp of Chr11 at positions 8758-8958).  Soft-clip
    sites at these host coordinates are false detections.

    Returns list of (host_chr, start, end) tuples for flanking regions.
    """
    if not construct_ref.exists():
        return []

    blast_out = workdir / "_construct_vs_host_flanking.tsv"
    result = subprocess.run(
        ["blastn", "-task", "megablast",
         "-query", str(construct_ref), "-db", str(host_ref),
         "-outfmt", "6 qseqid qlen qstart qend sseqid sstart send pident length",
         "-evalue", "1e-20", "-max_target_seqs", "5",
         "-num_threads", str(threads),
         "-out", str(blast_out)],
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0 or not blast_out.exists():
        stderr_tail = (result.stderr or b"").decode("utf-8", errors="replace")[-400:]
        log(f"[filter_b] blastn rc={result.returncode}, skipping flanking scan. "
            f"stderr: {stderr_tail}")
        return []

    flanking: list[tuple[str, int, int]] = []
    with open(blast_out) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            qlen = int(cols[1])
            q_start = int(cols[2])
            q_end = int(cols[3])
            s_chr = cols[4]
            s_start = int(cols[5])
            s_end = int(cols[6])
            pident = float(cols[7])
            aln_len = int(cols[8])

            if pident < CONSTRUCT_FLANK_PIDENT or aln_len < CONSTRUCT_FLANK_MIN_LEN:
                continue

            # Only count hits near the ends of a construct entry (flanking),
            # not internal matches (e.g., host-derived promoters handled
            # separately by _filter_host_endogenous).
            near_start = q_start <= 100
            near_end = q_end >= qlen - 100
            if not (near_start or near_end):
                continue

            lo = min(s_start, s_end)
            hi = max(s_start, s_end)
            flanking.append((s_chr, lo, hi))

    if flanking:
        log(f"  Construct-flanking regions found: {len(flanking)}")
        for chrom, lo, hi in flanking:
            log(f"    {chrom}:{lo:,}-{hi:,} ({hi - lo + 1}bp)")

    return flanking


# ---------------------------------------------------------------------------
# Per-site overlap test
# ---------------------------------------------------------------------------

def _site_overlaps_flanking(
    site_chr: str,
    site_pos: int,
    flanking_regions: list[tuple[str, int, int]],
    slop: int = CONSTRUCT_FLANK_SLOP,
) -> tuple[bool, str]:
    """Check if a detection site overlaps a construct-flanking host region."""
    for chrom, lo, hi in flanking_regions:
        if chrom == site_chr and lo - slop <= site_pos <= hi + slop:
            return True, f"{chrom}:{lo:,}-{hi:,}"
    return False, ""


# ---------------------------------------------------------------------------
# Pure filter-B check (consumed by verdict.compute_verdict)
# ---------------------------------------------------------------------------

def filter_b_flanking_hit(
    site_chr: str,
    site_pos: int,
    flanking_regions: list[tuple[str, int, int]],
    slop: int = CONSTRUCT_FLANK_SLOP,
) -> tuple[str, int, int] | None:
    """Return the slop-expanded flanking hit that contains the site, or None.

    Pure function — no BLAST, no file I/O. Used by the caller that builds a
    :class:`~scripts.s05.verdict.FilterEvidence` to populate
    ``FilterEvidence.flanking_hit`` with coordinates whose simple
    ``lo <= site_pos <= hi`` check inside ``compute_verdict`` is equivalent to
    the slop-aware test in :func:`_site_overlaps_flanking`.
    """
    for chrom, lo, hi in flanking_regions:
        if chrom == site_chr and lo - slop <= site_pos <= hi + slop:
            return (chrom, lo - slop, hi + slop)
    return None
