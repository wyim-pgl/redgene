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
call inside ``_check_chimeric_assembly``. The returned ``off_target`` list is
assigned directly to ``FilterEvidence.off_target_chrs`` by the caller and
consumed by ``compute_verdict`` Rule 3.
"""
from __future__ import annotations

import subprocess
from collections import defaultdict
from pathlib import Path

from .primitives import log


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
            stderr=subprocess.PIPE,
        )
        if result.returncode != 0 or not blast_out.exists():
            stderr_tail = (result.stderr or b"").decode("utf-8", errors="replace")[-400:]
            log(f"[filter_c] blastn rc={result.returncode}, skipping chimera check "
                f"for {insert_fasta.name}. stderr: {stderr_tail}")
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
