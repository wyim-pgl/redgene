"""Primitive helpers and dataclasses for the s05 package (Issue #4, Session 1).

Extracted verbatim from ``scripts/s05_insert_assembly.py`` (the monolith) so
every downstream submodule can import the same symbols without cycling back
through the entrypoint. Keeping them here means:

  * ``log`` — single stderr logger shared by every step-5 module.
  * ``revcomp`` — reverse-complement helper used by site discovery, assembly,
    and annotation modules.
  * ``read_fasta`` / ``write_fasta`` / ``_read_fq_seqs`` — tiny FASTA / FASTQ
    I/O routines used across read extraction, assembly, and annotation.
  * ``JunctionCluster``, ``InsertionSite``, ``LegacyJunction``, ``TierResult``
    — dataclasses consumed by site_discovery / classify / read_extraction /
    assembly.

No external side effects; no imports from the monolith or any later-stage
module (see ``tests/test_s05_import_dag.py`` for the enforced DAG).
"""
from __future__ import annotations

import gzip
import sys
from dataclasses import dataclass, field
from pathlib import Path


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

STEP = "s05_insert_assembly"


def log(msg: str) -> None:
    print(f"[{STEP}] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Reverse complement
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class JunctionCluster:
    """A cluster of soft-clipped reads at a genomic position."""
    host_chr: str
    position: int           # median position of cluster
    clip_direction: str     # 'right' (fwd read 3' clip) or 'left' (rev read 5' clip)
    clipped_seqs: list[str] = field(default_factory=list)
    consensus_clip: str = ""
    n_reads: int = 0
    maps_to_host: bool = False
    element_hits: list[str] = field(default_factory=list)


@dataclass
class InsertionSite:
    """An insertion site with paired 5'/3' junction clusters."""
    site_id: str
    host_chr: str
    junction_5p: JunctionCluster | None = None  # forward reads clipped on right
    junction_3p: JunctionCluster | None = None  # reverse reads clipped on left
    confidence: str = "low"
    seed_5p: str = ""       # consensus clip from 5' junction (insert start)
    seed_3p: str = ""       # consensus clip from 3' junction (insert end)
    is_validated: bool = False

    # Validation details
    clips_are_different: bool = False
    clips_not_in_host: bool = False
    has_element_hits: bool = False

    # Positions for read extraction
    pos_5p: int = 0
    pos_3p: int = 0


@dataclass
class LegacyJunction:
    """Junction from step 6 junctions.tsv (fallback mode)."""
    contig_name: str
    contig_len: int
    host_chr: str
    host_start: int
    host_end: int
    host_strand: str
    construct_element: str
    construct_start: int
    construct_end: int
    junction_pos_host: int
    junction_type: str
    confidence: str
    host_mapq: int


@dataclass
class TierResult:
    """Classification result for one insertion site."""
    site_id: str
    chrom: str = ""
    pos: int = 0
    transgene_positive: bool = False
    clip_5p_len: int = 0
    clip_3p_len: int = 0
    hit_5p: str = ""            # best transgene_db hit for 5p clip
    hit_5p_identity: float = 0
    hit_5p_aln_len: int = 0
    hit_5p_source: str = ""     # "element_db" or "univec"
    hit_3p: str = ""
    hit_3p_identity: float = 0
    hit_3p_aln_len: int = 0
    hit_3p_source: str = ""


# ---------------------------------------------------------------------------
# FASTA / FASTQ I/O
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = ""
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name:
        seqs[name] = "".join(parts)
    return seqs


def write_fasta(path: Path, name: str, seq: str, wrap: int = 80) -> None:
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), wrap):
            fh.write(seq[i:i + wrap] + "\n")


def _read_fq_seqs(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs: list[str] = []
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip().upper())
    return seqs


def _parse_src_tag(sseqid: str, default: str) -> tuple[str, str]:
    """Split a BLAST sseqid of the form 'accession|src=tag' into (accession, tag).

    v2 DB headers carry a '|src=<tag>' suffix to encode one of the 4-way source
    tags (element_db / payload / sample_contig / univec). Legacy headers lack
    the suffix - in that case the caller-supplied `default` is used.

    Returns (accession, src_tag). If no '|src=' separator is present,
    (sseqid, default) is returned unchanged.
    """
    accession, sep, src_tag = sseqid.partition("|src=")
    if not sep:
        return sseqid, default
    return accession, src_tag
