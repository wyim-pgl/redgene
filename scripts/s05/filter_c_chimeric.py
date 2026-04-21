"""Filter C — multi-locus chimeric assembly check (Issue #4, Session 1).

Extracted verbatim from ``scripts/s05_insert_assembly.py`` so the BLAST-driven
off-target discovery can be unit-tested without the full monolith.

Filter C rationale
------------------
If an assembled insert's host-aligned portions map to >= 2 different
chromosomes (besides the site's own), the assembler merged reads from
unrelated loci via shared element homology (e.g., CaMV 35S promoter paralogs
scattered across the genome).  Strict per-BLAST-hit identity (>= 98%) is used
to distinguish actual chimeric DNA pieces from low-level element homologies
(~80-90%).  See ``scripts/s05/verdict.py`` Rule 3.

This module is intentionally side-effect free aside from the BLAST subprocess
call inside ``_check_chimeric_assembly``. The "pure filter-C check"
(``filter_c_offtarget_chrs``) builds the ``FilterEvidence.off_target_chrs``
list consumed by ``compute_verdict`` without any I/O of its own.
"""
from __future__ import annotations

import subprocess
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Thresholds (mirrored from the monolith so behaviour is bit-identical).
# ---------------------------------------------------------------------------
# Uses strict identity (>= 98%) to distinguish actual chimeric DNA pieces from
# low-level element homologies (e.g., 35S promoter paralogs at 80-90%).
CHIMERIC_MIN_PIDENT = 98.0      # strict identity for chimeric detection
CHIMERIC_MIN_OFFTARGET_BP = 150  # min bp on off-target chromosome to count


# ---------------------------------------------------------------------------
# BLAST-driven chimera detection
# ---------------------------------------------------------------------------

def _check_chimeric_assembly(
    insert_fasta: Path,
    host_ref: Path,
    site_chr: str,
    workdir: Path,
    threads: int = 4,
) -> tuple[bool, list[tuple[str, int]]]:
    """Check if assembled insert contains DNA from multiple host chromosomes.

    Returns (is_chimeric, off_target_hits) where off_target_hits is a list of
    (chromosome, aligned_bp) for chromosomes other than site_chr.
    """
    blast_out = workdir / f"_{insert_fasta.stem}_vs_host_chrom.tsv"
    # Reuse existing BLAST output if available (from _blast_insert_vs_host)
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
        if result.returncode != 0 or not blast_out.exists():
            return False, []

    # Accumulate aligned bp per chromosome (strict identity to avoid
    # counting element-level homologies as chimeric evidence)
    chr_bp: dict[str, int] = defaultdict(int)
    with open(blast_out) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue
            s_chr = cols[3]
            pident = float(cols[4])
            aln_len = int(cols[5])
            if pident >= CHIMERIC_MIN_PIDENT:
                chr_bp[s_chr] += aln_len

    # Find off-target chromosomes with significant coverage
    off_target: list[tuple[str, int]] = []
    for chrom, bp in sorted(chr_bp.items(), key=lambda x: -x[1]):
        if chrom != site_chr and bp >= CHIMERIC_MIN_OFFTARGET_BP:
            off_target.append((chrom, bp))

    is_chimeric = len(off_target) >= 2
    return is_chimeric, off_target


# ---------------------------------------------------------------------------
# Pure filter-C check (consumed by verdict.compute_verdict)
# ---------------------------------------------------------------------------

def filter_c_offtarget_chrs(
    is_chimeric: bool,
    off_target_hits: list[tuple[str, int]],
) -> list[tuple[str, int]]:
    """Return the off-target chromosome list for FilterEvidence.

    Pure function — no BLAST, no file I/O. Thin pass-through today; kept as a
    named helper so that future refactors (Issue #11 M-3 deprecation of
    ``FilterEvidence.is_chimeric``) can adjust the contract in one place.
    """
    # Rule 3 in compute_verdict is driven by ``len(off_target_chrs)``; the
    # ``is_chimeric`` boolean is retained only for backward compatibility.
    del is_chimeric  # suppress unused-argument lint; kept for call-site parity
    return list(off_target_hits)
