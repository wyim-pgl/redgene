#!/usr/bin/env python3
"""Step 9 — Targeted seed extension assembly for transgene insert reconstruction.

For each insertion site detected by step 6:
  1. Extract junction read seeds (soft-clipped T-DNA portions from host BAM)
  2. Extract candidate reads from junction region (reads + mates within flank)
  3. Build full k-mer index from candidate reads
  4. Extend 5' and 3' seeds toward each other via majority vote consensus
  5. Merge overlapping extensions into final insert
  6. Annotate via BLAST against element database

Algorithm: Pure-Python k-mer indexed seed extension with majority vote
consensus. Uses ALL k-mers (not prefix-only) for read recruitment, which
is feasible because the candidate read pool is small (~5K-50K reads from
the junction region, not 300K+ unmapped reads).

Inputs:
  - junctions.tsv from step 6
  - Host BAM from step 7
  - Element database for annotation

Outputs:
  - insert_only.fasta   — assembled insert sequence(s)
  - element_annotation.tsv — BLAST hits along insert
  - border_hits.tsv      — T-DNA border motif locations
  - s09_stats.txt        — convergence and assembly statistics
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pysam

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

STEP = "s09_targeted_assembly"


def log(msg: str) -> None:
    print(f"[{STEP}] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class Junction:
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
    junction_type: str  # LB or RB
    confidence: str
    host_mapq: int


@dataclass
class InsertionSite:
    """A paired (or single) insertion site grouping LB/RB junctions."""
    site_id: str
    junctions: list[Junction] = field(default_factory=list)
    host_chr: str = ""
    pos_5p: int = 0  # upstream junction position
    pos_3p: int = 0  # downstream junction position
    contig_5p: str = ""  # contig name for 5' junction
    contig_3p: str = ""  # contig name for 3' junction


# ---------------------------------------------------------------------------
# Junction parsing
# ---------------------------------------------------------------------------

def parse_junctions(tsv_path: Path) -> list[Junction]:
    """Parse junctions.tsv from step 6.

    Deduplicates by (contig_name, host_chr, junction_pos_host), keeping
    the entry with the highest MAPQ.  Each contig can have multiple element
    matches producing duplicate rows that would otherwise create spurious
    insertion sites.
    """
    raw: list[Junction] = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            raw.append(Junction(
                contig_name=row["contig_name"],
                contig_len=int(row["contig_len"]),
                host_chr=row["host_chr"],
                host_start=int(row["host_start"]),
                host_end=int(row["host_end"]),
                host_strand=row["host_strand"],
                construct_element=row["construct_element"],
                construct_start=int(row["construct_start"]),
                construct_end=int(row["construct_end"]),
                junction_pos_host=int(row["junction_pos_host"]),
                junction_type=row["junction_type"],
                confidence=row["confidence"],
                host_mapq=int(row["host_mapq"]),
            ))

    # Deduplicate: keep highest-MAPQ entry per (contig, chr, position)
    best: dict[tuple[str, str, int], Junction] = {}
    for j in raw:
        key = (j.contig_name, j.host_chr, j.junction_pos_host)
        if key not in best or j.host_mapq > best[key].host_mapq:
            best[key] = j

    return list(best.values())


def read_fasta(path: Path) -> dict[str, str]:
    """Read FASTA file into {name: sequence} dict."""
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
    """Write a single sequence to FASTA."""
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), wrap):
            fh.write(seq[i:i + wrap] + "\n")


# ---------------------------------------------------------------------------
# PAF parsing
# ---------------------------------------------------------------------------

@dataclass
class PafHit:
    query_name: str
    query_len: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_len: int
    target_start: int
    target_end: int
    mapq: int


def parse_paf(path: Path, contig_names: set[str] | None = None) -> list[PafHit]:
    """Parse PAF, optionally filtering to specific contigs."""
    hits = []
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip().split("\t")
            if len(cols) < 12:
                continue
            qname = cols[0]
            if contig_names and qname not in contig_names:
                continue
            hits.append(PafHit(
                query_name=qname,
                query_len=int(cols[1]),
                query_start=int(cols[2]),
                query_end=int(cols[3]),
                strand=cols[4],
                target_name=cols[5],
                target_len=int(cols[6]),
                target_start=int(cols[7]),
                target_end=int(cols[8]),
                mapq=int(cols[11]),
            ))
    return hits


def get_construct_portion(
    contig_name: str,
    contig_seq: str,
    host_hits: list[PafHit],
    construct_hits: list[PafHit],
) -> str:
    """Extract the construct-derived portion of a chimeric contig.

    Identifies host-aligned regions and returns the non-host portion,
    keeping a small overlap anchor (50bp).
    """
    contig_len = len(contig_seq)

    # Find host-covered ranges on the contig
    host_ranges = [(h.query_start, h.query_end) for h in host_hits
                   if h.query_name == contig_name]
    # Find construct-covered ranges
    construct_ranges = [(h.query_start, h.query_end) for h in construct_hits
                        if h.query_name == contig_name]

    if not construct_ranges:
        return contig_seq  # no construct alignment, return whole thing

    # Merge construct ranges to find the construct span
    construct_ranges.sort()
    merged_start = construct_ranges[0][0]
    merged_end = construct_ranges[-1][1]

    # If host is at the ends, trim with 50bp anchor overlap
    anchor = 50
    start = max(0, merged_start - anchor)
    end = min(contig_len, merged_end + anchor)

    return contig_seq[start:end]


# ---------------------------------------------------------------------------
# Insertion site grouping
# ---------------------------------------------------------------------------

def group_insertion_sites(junctions: list[Junction]) -> list[InsertionSite]:
    """Group junctions into insertion sites.

    Rules:
    - Same chr, positions within 50kb → pair as same insertion
    - Take highest-MAPQ junction per side
    - Different chromosomes → separate sites
    """
    if not junctions:
        return []

    # Use all junctions (MAPQ filtering already done via deduplication —
    # keeping highest-MAPQ per contig+chr+pos).  Low-MAPQ junctions may
    # still be valid insertion sites (e.g. rice G281 Chr3:16439719 MAPQ=10).
    # Group by chromosome
    by_chr: dict[str, list[Junction]] = {}
    for j in junctions:
        by_chr.setdefault(j.host_chr, []).append(j)

    sites: list[InsertionSite] = []
    site_idx = 0

    for chrom, juncs in by_chr.items():
        # Sort by position
        juncs.sort(key=lambda j: j.junction_pos_host)

        # Try to pair junctions within 50kb
        used = set()
        for i, j1 in enumerate(juncs):
            if i in used:
                continue
            paired = None
            for k, j2 in enumerate(juncs):
                if k in used or k == i:
                    continue
                if abs(j2.junction_pos_host - j1.junction_pos_host) <= 50000:
                    paired = k
                    break

            site_idx += 1
            site = InsertionSite(
                site_id=f"insertion_{site_idx}",
                host_chr=chrom,
            )

            if paired is not None:
                used.add(i)
                used.add(paired)
                p1, p2 = j1.junction_pos_host, juncs[paired].junction_pos_host
                site.pos_5p = min(p1, p2)
                site.pos_3p = max(p1, p2)
                if p1 <= p2:
                    site.contig_5p = j1.contig_name
                    site.contig_3p = juncs[paired].contig_name
                    site.junctions = [j1, juncs[paired]]
                else:
                    site.contig_5p = juncs[paired].contig_name
                    site.contig_3p = j1.contig_name
                    site.junctions = [juncs[paired], j1]
            else:
                used.add(i)
                site.pos_5p = j1.junction_pos_host
                site.pos_3p = 0
                site.contig_5p = j1.contig_name
                site.junctions = [j1]

            sites.append(site)

    return sites


# ---------------------------------------------------------------------------
# Candidate read extraction
# ---------------------------------------------------------------------------

def extract_candidate_reads(
    host_bam: Path,
    sites: list[InsertionSite],
    out_r1: Path,
    out_r2: Path,
    flank: int = 5000,
    threads: int = 4,
) -> int:
    """Extract candidate reads from host BAM for Pilon gap filling.

    Collects: junction-flanking reads + unmapped + mate-unmapped.
    Returns read pair count.
    """
    tmp_bams = []
    tmp_dir = out_r1.parent

    # Junction-flanking reads
    for site in sites:
        for junc in site.junctions:
            region = (f"{junc.host_chr}:"
                      f"{max(1, junc.junction_pos_host - flank)}-"
                      f"{junc.junction_pos_host + flank}")
            bam_tmp = tmp_dir / f"_junc_{junc.host_chr}_{junc.junction_pos_host}.bam"
            subprocess.run(
                ["samtools", "view", "-b", "-@", str(threads),
                 str(host_bam), region],
                stdout=open(bam_tmp, "wb"),
                stderr=subprocess.DEVNULL,
                check=True,
            )
            tmp_bams.append(bam_tmp)

    # Unmapped reads (construct-internal)
    unmapped_bam = tmp_dir / "_unmapped.bam"
    subprocess.run(
        ["samtools", "view", "-b", "-f", "4", "-@", str(threads), str(host_bam)],
        stdout=open(unmapped_bam, "wb"),
        stderr=subprocess.DEVNULL,
        check=True,
    )
    tmp_bams.append(unmapped_bam)

    # Mate-unmapped (one read maps near junction, mate in insert)
    mate_unmapped_bam = tmp_dir / "_mate_unmapped.bam"
    subprocess.run(
        ["samtools", "view", "-b", "-f", "8", "-@", str(threads), str(host_bam)],
        stdout=open(mate_unmapped_bam, "wb"),
        stderr=subprocess.DEVNULL,
        check=True,
    )
    tmp_bams.append(mate_unmapped_bam)

    # Merge all
    merged_bam = tmp_dir / "_candidate_merged.bam"
    if len(tmp_bams) > 1:
        subprocess.run(
            ["samtools", "merge", "-f", "-@", str(threads), str(merged_bam)]
            + [str(b) for b in tmp_bams],
            stderr=subprocess.DEVNULL,
            check=True,
        )
    else:
        shutil.copy2(tmp_bams[0], merged_bam)

    # Name-sort and convert to FASTQ
    nsort_bam = tmp_dir / "_candidate_nsort.bam"
    subprocess.run(
        ["samtools", "sort", "-n", "-@", str(threads),
         str(merged_bam), "-o", str(nsort_bam)],
        stderr=subprocess.DEVNULL,
        check=True,
    )
    subprocess.run(
        ["samtools", "fastq",
         "-1", str(out_r1), "-2", str(out_r2),
         "-s", "/dev/null", "-0", "/dev/null",
         str(nsort_bam)],
        stderr=subprocess.DEVNULL,
        check=True,
    )

    # Count reads
    count = 0
    with pysam.FastxFile(str(out_r1)) as fh:
        for _ in fh:
            count += 1

    # Cleanup temp BAMs
    for b in tmp_bams:
        b.unlink(missing_ok=True)
    merged_bam.unlink(missing_ok=True)
    nsort_bam.unlink(missing_ok=True)

    return count


# ---------------------------------------------------------------------------
# Pseudo-reference construction
# ---------------------------------------------------------------------------


def build_pseudo_reference(
    site: InsertionSite,
    contig_seqs: dict[str, str],
    host_hits: list[PafHit],
    construct_hits: list[PafHit],
    host_ref: Path,
    flank_size: int = 20000,
    gap_size: int = 30000,
) -> str:
    """Build pseudo-reference for one insertion site.

    Structure: [flank_5p] + [insert_5p] + [N-gap] + [insert_3p] + [flank_3p]

    For single-junction sites, BOTH genomic flanks around the junction
    position are used so Pilon can fill from both directions.

    Returns the assembled pseudo-reference sequence.
    """
    genome = pysam.FastaFile(str(host_ref))
    chr_len = genome.get_reference_length(site.host_chr)

    # Always extract BOTH flanks around the junction position(s)
    # 5' flank: upstream of the leftmost junction
    flank_5p_start = max(0, site.pos_5p - flank_size)
    flank_5p = genome.fetch(site.host_chr, flank_5p_start, site.pos_5p)

    # 3' flank: downstream of the rightmost junction (or same position for single)
    anchor_3p = site.pos_3p if site.pos_3p > 0 else site.pos_5p
    flank_3p_end = min(chr_len, anchor_3p + flank_size)
    flank_3p = genome.fetch(site.host_chr, anchor_3p, flank_3p_end)

    genome.close()

    # Extract construct portions from junction contigs
    insert_5p = ""
    if site.contig_5p and site.contig_5p in contig_seqs:
        insert_5p = get_construct_portion(
            site.contig_5p, contig_seqs[site.contig_5p],
            host_hits, construct_hits,
        )

    insert_3p = ""
    if site.contig_3p and site.contig_3p in contig_seqs and site.contig_3p != site.contig_5p:
        insert_3p = get_construct_portion(
            site.contig_3p, contig_seqs[site.contig_3p],
            host_hits, construct_hits,
        )

    # Assemble: flank_5p + insert_5p + N-gap + insert_3p + flank_3p
    pseudo = flank_5p + insert_5p + ("N" * gap_size) + insert_3p + flank_3p

    log(f"  Pseudo-ref: {len(flank_5p)}bp flank_5p + {len(insert_5p)}bp insert_5p"
        f" + {gap_size}N + {len(insert_3p)}bp insert_3p + {len(flank_3p)}bp flank_3p"
        f" = {len(pseudo)}bp total")

    return pseudo


# ---------------------------------------------------------------------------
# Seed extension assembler (k-mer indexed, TASR-like)
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


class SeedExtender:
    """Strand-aware k-mer extension with majority vote consensus.

    Key PE rules:
    - R1 forward: used as-is for extension
    - R2: reverse-complemented before use (restoring fragment strand)
    - Raw R2 (without RC) is NEVER indexed or used for extension
      → prevents chimeric assembly from head-to-head T-DNA constructs
    """

    def __init__(
        self,
        k: int = 15,
        min_overlap: int = 20,
        min_depth: int = 2,
        min_ratio: float = 0.7,
    ):
        self.k = k
        self.min_overlap = min_overlap
        self.min_depth = min_depth
        self.min_ratio = min_ratio
        self.seqs: list[str] = []  # R1 forward + R2 RC (fragment strand)
        # kmer → [(seq_idx, kmer_position), ...]
        self.kmer_index: dict[str, list[tuple[int, int]]] = defaultdict(list)

    def load_paired_reads(self, r1_path: Path, r2_path: Path) -> int:
        """Load paired FASTQ. R1=forward, R2=reverse-complement (fragment strand)."""
        r1_seqs = self._read_fastq_seqs(r1_path)
        r2_seqs = self._read_fastq_seqs(r2_path)
        k = self.k
        n = 0
        for r1, r2 in zip(r1_seqs, r2_seqs):
            if len(r1) < self.min_overlap or len(r2) < self.min_overlap:
                continue
            if "N" in r1 or "N" in r2:
                continue
            r2_rc = revcomp(r2)  # RC to restore fragment strand

            for seq in (r1, r2_rc):
                idx = len(self.seqs)
                self.seqs.append(seq)
                for p in range(len(seq) - k + 1):
                    self.kmer_index[seq[p:p + k]].append((idx, p))
            n += 1
        return n

    def add_seqs(self, seqs: list[str]) -> int:
        """Add pre-processed sequences to pool and index."""
        existing = set(self.seqs)
        k = self.k
        n_new = 0
        for seq in seqs:
            seq = seq.upper()
            if seq in existing or len(seq) < self.min_overlap or "N" in seq:
                continue
            existing.add(seq)
            idx = len(self.seqs)
            self.seqs.append(seq)
            for p in range(len(seq) - k + 1):
                self.kmer_index[seq[p:p + k]].append((idx, p))
            n_new += 1
        return n_new

    def _extend_right(self, seed: str, used: set[int]) -> str:
        """Extend seed to the right with used-read tracking.

        The `used` set is shared across iterations in run() to prevent
        cyclic re-recruitment of the same reads (which causes infinite
        looping through tandem repeats).
        """
        k = self.k
        min_ovl = self.min_overlap
        current = seed

        while True:
            if len(current) < k:
                break
            tail_kmer = current[-k:]
            hits = self.kmer_index.get(tail_kmer)
            if not hits:
                break

            votes: dict[int, Counter] = defaultdict(Counter)
            for seq_idx, kpos in hits:
                if seq_idx in used:
                    continue
                overlap = kpos + k
                if overlap < min_ovl:
                    continue
                read_seq = self.seqs[seq_idx]
                if current[-overlap:] == read_seq[:overlap]:
                    for j in range(overlap, len(read_seq)):
                        votes[j - overlap][read_seq[j]] += 1
                    used.add(seq_idx)

            if not votes:
                break

            ext = []
            for p in range(max(votes.keys()) + 1):
                if p not in votes:
                    break
                v = votes[p]
                total = sum(v.values())
                if total < self.min_depth:
                    break
                best_base, best_count = v.most_common(1)[0]
                if best_count / total < self.min_ratio:
                    break
                ext.append(best_base)

            if not ext:
                break
            current = current + "".join(ext)

        return current

    def run(self, seed: str, max_iterations: int = 100) -> str:
        """Iteratively extend seed until convergence.

        Tracks used reads across all iterations and both directions to
        prevent cyclic looping through repeated sequences.
        """
        current = seed
        used: set[int] = set()
        for i in range(max_iterations):
            # Right extension
            extended = self._extend_right(current, used)
            # Left extension (revcomp → right → revcomp)
            rc_extended = self._extend_right(revcomp(extended), used)
            extended = revcomp(rc_extended)

            growth = len(extended) - len(current)
            if growth == 0:
                log(f"      Iter {i + 1}: {len(current):,}bp → converged")
                break
            log(f"      Iter {i + 1}: {len(current):,} → {len(extended):,}bp "
                f"(+{growth})")
            current = extended
        return current

    @staticmethod
    def _read_fastq_seqs(path: Path) -> list[str]:
        opener = gzip.open if str(path).endswith(".gz") else open
        seqs: list[str] = []
        with opener(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    seqs.append(line.strip().upper())
        return seqs


# ---------------------------------------------------------------------------
# Candidate read extraction from junction region
# ---------------------------------------------------------------------------

def extract_junction_candidates(
    host_bam: Path,
    junctions: list[Junction],
    out_r1: Path,
    out_r2: Path,
    flank: int = 1000,
    threads: int = 4,
) -> int:
    """Extract candidate reads from host BAM near junction positions.

    Extracts reads mapping within ±flank of each junction position.
    Uses samtools view with region queries — this naturally includes reads
    whose mates are unmapped (T-DNA content) because the mapped end is in
    the junction region. Name-sorting and fastq conversion recovers both
    mates as paired reads.

    Returns read pair count.
    """
    tmp_dir = out_r1.parent

    # Build regions for junction-flanking reads
    regions = []
    seen = set()
    for j in junctions:
        start = max(1, j.junction_pos_host - flank)
        end = j.junction_pos_host + flank
        key = (j.host_chr, start, end)
        if key not in seen:
            regions.append(f"{j.host_chr}:{start}-{end}")
            seen.add(key)

    # Extract reads from junction regions into a single BAM
    region_bam = tmp_dir / "_junc_regions.bam"
    subprocess.run(
        ["samtools", "view", "-b", "-@", str(threads),
         str(host_bam)] + regions,
        stdout=open(region_bam, "wb"),
        stderr=subprocess.DEVNULL, check=True,
    )

    # Collect read names from junction regions, then extract both mates
    # This is needed because samtools view only returns reads physically
    # in the region — the mate at a distant locus (or unmapped) is not included.
    # We use samtools view -N (name list) on the full BAM to get both mates.
    namelist = tmp_dir / "_junc_readnames.txt"
    subprocess.run(
        ["samtools", "view", str(region_bam)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True,
    )
    # Extract read names (field 1) from region BAM
    result = subprocess.run(
        ["samtools", "view", str(region_bam)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True, text=True,
    )
    names = set()
    for line in result.stdout.splitlines():
        names.add(line.split("\t", 1)[0])
    with open(namelist, "w") as fh:
        for n in sorted(names):
            fh.write(n + "\n")
    log(f"    Junction region read names: {len(names):,}")

    # Extract both mates of all junction-region reads from the full BAM
    both_mates_bam = tmp_dir / "_junc_both_mates.bam"
    subprocess.run(
        ["samtools", "view", "-b", "-N", str(namelist),
         "-@", str(threads), str(host_bam)],
        stdout=open(both_mates_bam, "wb"),
        stderr=subprocess.DEVNULL, check=True,
    )

    # Name-sort and convert to paired FASTQ
    nsort_bam = tmp_dir / "_candidate_nsort.bam"
    subprocess.run(
        ["samtools", "sort", "-n", "-@", str(threads),
         str(both_mates_bam), "-o", str(nsort_bam)],
        stderr=subprocess.DEVNULL, check=True,
    )

    subprocess.run(
        ["samtools", "fastq", "-@", str(threads),
         "-1", str(out_r1), "-2", str(out_r2),
         "-s", "/dev/null", "-0", "/dev/null", str(nsort_bam)],
        stderr=subprocess.DEVNULL, check=True,
    )

    # Count
    count = 0
    with pysam.FastxFile(str(out_r1)) as fh:
        for _ in fh:
            count += 1

    # Cleanup
    for p in [region_bam, namelist, both_mates_bam, nsort_bam]:
        p.unlink(missing_ok=True)

    return count


# ---------------------------------------------------------------------------
# Junction read seed extraction
# ---------------------------------------------------------------------------

def extract_junction_seeds(
    host_bam: Path,
    junctions: list[Junction],
    host_ref: Path | None = None,
    host_flank_size: int = 1000,
    min_clip: int = 15,
    window: int = 10,
) -> list[tuple[str, str]]:
    """Extract T-DNA seed sequences from soft-clipped junction reads.

    For each junction position, collects soft-clipped reads from the host BAM,
    builds a majority-vote consensus from overlapping clips, and returns
    seed sequences with host genome flank prepended/appended.

    Host flank provides a known anchor for Pilon gap fill:
    - R-clip (T-DNA extends right): [host_1kb] + [T-DNA seed]
    - L-clip (T-DNA extends left):  [T-DNA seed] + [host_1kb]

    Returns list of (seed_name, seed_sequence) tuples.
    """
    seeds: list[tuple[str, str]] = []
    seen_positions: set[tuple[str, int]] = set()

    # Load host reference for flanking sequence
    genome = None
    if host_ref is not None:
        genome = pysam.FastaFile(str(host_ref))

    bam = pysam.AlignmentFile(str(host_bam), "rb")

    for j in junctions:
        chrom = j.host_chr
        jpos = j.junction_pos_host
        key = (chrom, jpos)
        if key in seen_positions:
            continue
        seen_positions.add(key)

        right_clips: list[tuple[int, str]] = []
        left_clips: list[tuple[int, str]] = []

        for read in bam.fetch(chrom, max(0, jpos - 200), jpos + 200):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            cigar = read.cigartuples
            if cigar is None:
                continue
            seq = read.query_sequence

            if cigar[-1][0] == 4 and cigar[-1][1] >= min_clip:
                sclip_len = cigar[-1][1]
                if abs(read.reference_end - jpos) <= window:
                    offset = read.reference_end - jpos
                    right_clips.append((offset, seq[-sclip_len:]))

            if cigar[0][0] == 4 and cigar[0][1] >= min_clip:
                sclip_len = cigar[0][1]
                if abs(read.reference_start - jpos) <= window:
                    offset = read.reference_start - jpos
                    left_clips.append((offset, seq[:sclip_len]))

        for direction, clips in [("R", right_clips), ("L", left_clips)]:
            if len(clips) < 2:
                continue
            votes: dict[int, Counter] = defaultdict(Counter)
            for offset, tdna in clips:
                for i, base in enumerate(tdna):
                    if direction == "R":
                        p = offset + i
                    else:
                        p = offset - len(tdna) + i
                    votes[p][base.upper()] += 1

            if direction == "R":
                consensus = []
                for p in range(min(votes), max(votes) + 1):
                    if p not in votes:
                        break
                    v = votes[p]
                    total = sum(v.values())
                    best, cnt = v.most_common(1)[0]
                    if total >= 2 and cnt / total >= 0.51:
                        consensus.append(best)
                    else:
                        break
            else:
                consensus = []
                for p in range(max(votes), min(votes) - 1, -1):
                    if p not in votes:
                        break
                    v = votes[p]
                    total = sum(v.values())
                    best, cnt = v.most_common(1)[0]
                    if total >= 2 and cnt / total >= 0.51:
                        consensus.append(best)
                    else:
                        break
                consensus.reverse()

            if len(consensus) >= 20:
                seed_name = f"junction_{chrom}_{jpos}_{direction}"
                tdna_seq = "".join(consensus)
                host_flank = ""
                flank_side = ""

                # Extract host genome flank at junction side
                if genome is not None:
                    chr_len = genome.get_reference_length(chrom)
                    if direction == "R":
                        # T-DNA extends right → host is LEFT of junction
                        flank_start = max(0, jpos - host_flank_size)
                        host_flank = genome.fetch(chrom, flank_start, jpos).upper()
                        flank_side = "prefix"
                        log(f"    Seed: {seed_name} = {len(host_flank)}bp host + "
                            f"{len(tdna_seq)}bp T-DNA (from {len(clips)} reads)")
                    else:
                        # T-DNA extends left → host is RIGHT of junction
                        flank_end = min(chr_len, jpos + host_flank_size)
                        host_flank = genome.fetch(chrom, jpos, flank_end).upper()
                        flank_side = "suffix"
                        log(f"    Seed: {seed_name} = {len(tdna_seq)}bp T-DNA + "
                            f"{len(host_flank)}bp host (from {len(clips)} reads)")
                else:
                    log(f"    Seed: {seed_name} = {len(tdna_seq)}bp "
                        f"(from {len(clips)} reads)")

                seeds.append((seed_name, tdna_seq, host_flank, flank_side))

    bam.close()
    if genome is not None:
        genome.close()
    return seeds


# ---------------------------------------------------------------------------
# K-mer recruitment from unmapped reads (replaces BWA-based recruitment)
# ---------------------------------------------------------------------------


def recruit_by_kmer(
    contig_seq: str,
    unmapped_r1: Path,
    unmapped_r2: Path,
    k: int = 25,
) -> tuple[list[str], list[str]]:
    """K-mer based paired recruitment from unmapped reads.

    Builds k-mer set from contig (both strands), scans unmapped R1 and R2_rc
    for matches. Returns (r1_seqs, r2_rc_seqs) — R2 is already RC'd.
    Raw R2 is NOT used for matching (head-to-head safety).
    """
    # Build contig k-mer set (both strands)
    contig_kmers: set[str] = set()
    for seq in (contig_seq, revcomp(contig_seq)):
        for i in range(len(seq) - k + 1):
            contig_kmers.add(seq[i:i + k])

    opener = gzip.open
    r1_seqs = _read_fq_seqs(unmapped_r1, opener)
    r2_seqs = _read_fq_seqs(unmapped_r2, opener)

    recruited_r1: list[str] = []
    recruited_r2rc: list[str] = []
    stride = max(1, k // 2)

    for r1, r2 in zip(r1_seqs, r2_seqs):
        if "N" in r1 or "N" in r2:
            continue
        r2_rc = revcomp(r2)

        hit = False
        for seq in (r1, r2_rc):  # check R1 and R2_rc only, NOT raw R2
            for i in range(0, len(seq) - k + 1, stride):
                if seq[i:i + k] in contig_kmers:
                    hit = True
                    break
            if hit:
                break

        if hit:
            recruited_r1.append(r1)
            recruited_r2rc.append(r2_rc)

    return recruited_r1, recruited_r2rc


def _read_fq_seqs(path: Path, opener=None) -> list[str]:
    """Read FASTQ sequences only."""
    if opener is None:
        opener = gzip.open if str(path).endswith(".gz") else open
    seqs: list[str] = []
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip().upper())
    return seqs


def extract_unmapped_paired(
    host_bam: Path,
    workdir: Path,
    threads: int = 4,
) -> tuple[Path, Path]:
    """Extract unmapped reads from host BAM as paired FASTQ (cached)."""
    r1 = workdir / "_unmapped_R1.fq.gz"
    r2 = workdir / "_unmapped_R2.fq.gz"
    if r1.exists() and r2.exists():
        return r1, r2
    log(f"  Extracting unmapped reads from host BAM...")
    subprocess.run(
        f"samtools view -f 4 -@ {threads} -b {host_bam} "
        f"| samtools sort -n -@ {threads} - "
        f"| samtools fastq -1 {r1} -2 {r2} -s /dev/null -0 /dev/null - "
        f"2>/dev/null",
        shell=True, check=True,
    )
    n = 0
    with gzip.open(r1, "rt") as fh:
        for _ in fh:
            n += 1
    log(f"  Unmapped reads: {n // 4:,} pairs")
    return r1, r2


# ---------------------------------------------------------------------------
# Pilon gap fill
# ---------------------------------------------------------------------------


def _write_pool_fastq(extender: SeedExtender, r1_out: Path, r2_out: Path) -> None:
    """Write extender's read pool as pseudo-paired FASTQ for Pilon mapping.

    Extender stores R1 (even idx) and R2_rc (odd idx) alternating.
    For Pilon, R2 needs to be in original sequencing orientation (RC of R2_rc).
    """
    with open(r1_out, "w") as f1, open(r2_out, "w") as f2:
        for i in range(0, len(extender.seqs) - 1, 2):
            r1 = extender.seqs[i]
            r2_rc = extender.seqs[i + 1]
            r2_raw = revcomp(r2_rc)
            name = f"read_{i // 2}"
            f1.write(f"@{name}/1\n{r1}\n+\n{'I' * len(r1)}\n")
            f2.write(f"@{name}/2\n{r2_raw}\n+\n{'I' * len(r2_raw)}\n")


def _find_end_anchors(
    contig: str,
    extender: SeedExtender,
    search_len: int = 50,
) -> tuple[str, str]:
    """Find PE mate anchor sequences at both contig ends.

    Right anchor: R1 reads (even idx) at contig right end → their R2_rc mates
                  are ~insert_size further right → use as right anchor.
    Left anchor:  R2_rc reads (odd idx) at contig left end → their R1 mates
                  are ~insert_size further left → use as left anchor.

    Returns (left_anchor, right_anchor). Empty string if no anchor found.
    """
    k = extender.k

    def _collect_mates(region: str, read_parity: int) -> list[str]:
        """Collect mate sequences for reads of given parity matching region k-mers."""
        kmers = {region[i:i + k] for i in range(len(region) - k + 1)}
        mates: list[str] = []
        seen: set[int] = set()
        for kmer in kmers:
            for seq_idx, kpos in extender.kmer_index.get(kmer, []):
                if seq_idx in seen or seq_idx % 2 != read_parity:
                    continue
                seen.add(seq_idx)
                mate_idx = seq_idx ^ 1  # even↔odd
                if mate_idx < len(extender.seqs):
                    mates.append(extender.seqs[mate_idx])
        return mates

    # Right anchor: R1 (even=0) at right end → R2_rc mates point further right
    right_mates = _collect_mates(contig[-search_len:], read_parity=0)
    right_anchor = max(right_mates, key=len) if right_mates else ""

    # Left anchor: R2_rc (odd=1) at left end → R1 mates point further left
    left_mates = _collect_mates(contig[:search_len], read_parity=1)
    left_anchor = max(left_mates, key=len) if left_mates else ""

    return left_anchor, right_anchor


def pilon_fill(
    contig: str,
    extender: SeedExtender,
    r1_fq: Path,
    r2_fq: Path,
    workdir: Path,
    insert_size: int = 500,
    threads: int = 4,
) -> str:
    """Build PE-anchored scaffold and run Pilon to extend contig ends.

    Finds PE mates of reads at contig ends, places them as anchors
    ~insert_size away, creating internal gaps that Pilon can fill.

    Scaffold: [left_anchor] + [insert_size N] + [contig] + [insert_size N] + [right_anchor]
                              ↑ internal gap ↑            ↑ internal gap ↑
    Both gaps have known sequence on both sides → Pilon fills them.

    Returns extended contig sequence.
    """
    workdir.mkdir(parents=True, exist_ok=True)

    left_anchor, right_anchor = _find_end_anchors(contig, extender)

    if not left_anchor and not right_anchor:
        log(f"    Pilon: no PE anchors found, skipping")
        return contig

    anchor_info = []
    if left_anchor:
        anchor_info.append(f"L={len(left_anchor)}bp")
    if right_anchor:
        anchor_info.append(f"R={len(right_anchor)}bp")
    log(f"    Pilon scaffold: {' '.join(anchor_info)} anchors, "
        f"{insert_size}N gaps, contig={len(contig):,}bp")

    # Build scaffold: [left_anchor] + NNN + [contig] + NNN + [right_anchor]
    gap = "N" * insert_size
    parts: list[str] = []
    if left_anchor:
        parts.extend([left_anchor, gap])
    parts.append(contig)
    if right_anchor:
        parts.extend([gap, right_anchor])
    scaffold = "".join(parts)

    ref_fa = workdir / "pilon_ref.fasta"
    write_fasta(ref_fa, "scaffold", scaffold)

    # minimap2 short-read mapping
    bam = workdir / "pilon_mapped.bam"
    subprocess.run(
        f"minimap2 -ax sr -t {threads} --secondary=no "
        f"{ref_fa} {r1_fq} {r2_fq} 2>/dev/null "
        f"| samtools sort -@ {threads} -o {bam}",
        shell=True, check=True,
    )
    subprocess.run(["samtools", "index", str(bam)],
                   stderr=subprocess.DEVNULL, check=True)

    # Count mapped reads
    r = subprocess.run(
        ["samtools", "view", "-c", "-F", "4", str(bam)],
        stdout=subprocess.PIPE, text=True, stderr=subprocess.DEVNULL,
    )
    n_mapped = int(r.stdout.strip()) if r.stdout.strip() else 0

    if n_mapped < 5:
        log(f"    Pilon: only {n_mapped} mapped reads, skipping")
        return contig

    # Pilon
    pilon_prefix = workdir / "pilon_out"
    pilon_log = workdir / "pilon.log"
    with open(pilon_log, "w") as logfh:
        subprocess.run(
            ["pilon", "--genome", str(ref_fa), "--frags", str(bam),
             "--output", str(pilon_prefix), "--fix", "gaps",
             "--mindepth", "2", "--gapmargin", "100000"],
            stdout=logfh, stderr=subprocess.STDOUT,
        )
    if pilon_log.exists():
        with open(pilon_log) as fh:
            for line in fh:
                line = line.strip()
                if any(k in line for k in ["Gap", "fix", "Total", "Confirmed"]):
                    log(f"    [Pilon] {line}")

    pilon_fa = Path(f"{pilon_prefix}.fasta")
    if not pilon_fa.exists():
        return contig

    seqs = read_fasta(pilon_fa)
    if not seqs:
        return contig

    filled = list(seqs.values())[0]

    # Strip any remaining edge N's (unfilled anchor gaps at edges)
    filled = filled.strip("N")

    if len(filled) < len(contig) * 0.5:
        log(f"    Pilon: result too short ({len(filled)} vs {len(contig)}), "
            f"keeping original")
        return contig

    return filled


# ---------------------------------------------------------------------------
# Seed extension assembly (alternating K-mer ext + Pilon)
# ---------------------------------------------------------------------------


def extract_junction_host_flanks(
    host_ref: Path,
    junctions: list[Junction],
    flank_len: int = 200,
) -> list[tuple[str, str, str]]:
    """Extract host genomic flanking sequences at each junction position.

    For each junction, extracts flank_len bp on BOTH sides of the junction:
    - left_flank: host sequence upstream of junction (host side for R-clip seeds)
    - right_flank: host sequence downstream of junction (host side for L-clip seeds)

    Returns list of (junction_id, left_flank_seq, right_flank_seq).
    """
    genome = pysam.FastaFile(str(host_ref))
    flanks: list[tuple[str, str, str]] = []
    seen: set[tuple[str, int]] = set()

    for j in junctions:
        key = (j.host_chr, j.junction_pos_host)
        if key in seen:
            continue
        seen.add(key)

        chr_len = genome.get_reference_length(j.host_chr)
        pos = j.junction_pos_host

        left_start = max(0, pos - flank_len)
        left_flank = genome.fetch(j.host_chr, left_start, pos).upper()

        right_end = min(chr_len, pos + flank_len)
        right_flank = genome.fetch(j.host_chr, pos, right_end).upper()

        jid = f"{j.host_chr}:{pos}"
        flanks.append((jid, left_flank, right_flank))
        log(f"    Host flank at {jid}: L={len(left_flank)}bp, R={len(right_flank)}bp")

    genome.close()
    return flanks


def _check_host_reached(
    contig: str,
    host_ref: Path | None,
    junctions: list,
    workdir: Path,
    check_len: int = 500,
    local_flank: int = 50000,
) -> tuple[bool, bool]:
    """Check if contig ends match the local host genome near known junctions.

    Extracts a local region (±local_flank bp around each junction) from the
    host genome and maps the contig's first/last check_len bp against it
    using minimap2.  This avoids whole-genome false positives from host-derived
    promoters (Ubi1, Act1) because only the immediate junction neighbourhood
    is queried.

    Returns (left_reached, right_reached).
    """
    if host_ref is None or len(contig) < check_len * 3:
        return False, False

    import pysam as _pysam

    # Build local reference from junction neighbourhood
    local_fa = workdir / "_host_local.fa"
    genome = _pysam.FastaFile(str(host_ref))
    seen: set[tuple[str, int]] = set()
    with open(local_fa, "w") as fout:
        for j in junctions:
            key = (j.host_chr, j.junction_pos_host)
            if key in seen:
                continue
            seen.add(key)
            chr_len = genome.get_reference_length(j.host_chr)
            start = max(0, j.junction_pos_host - local_flank)
            end = min(chr_len, j.junction_pos_host + local_flank)
            seq = genome.fetch(j.host_chr, start, end).upper()
            fout.write(f">{j.host_chr}:{start}-{end}\n{seq}\n")
    genome.close()

    # Write contig ends as query
    query_fa = workdir / "_host_check.fa"
    left_end = contig[:check_len].upper()
    right_end = contig[-check_len:].upper()
    with open(query_fa, "w") as fout:
        fout.write(f">left_end\n{left_end}\n")
        fout.write(f">right_end\n{right_end}\n")

    # Map with minimap2 (asm5 preset for high-identity genome matching)
    paf = workdir / "_host_check.paf"
    subprocess.run(
        ["minimap2", "-c", "--secondary=no", "-t", "2",
         str(local_fa), str(query_fa)],
        stdout=open(paf, "w"), stderr=subprocess.DEVNULL,
    )

    # Parse PAF — require ≥200bp match at ≥0.95 identity
    left_hit = False
    right_hit = False
    if paf.exists():
        with open(paf) as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 12:
                    continue
                qname = cols[0]
                match_bp = int(cols[9])
                block_len = int(cols[10])
                if match_bp < 200 or block_len < 200:
                    continue
                identity = match_bp / block_len
                if identity < 0.95:
                    continue
                target = cols[5]
                if qname == "left_end":
                    left_hit = True
                elif qname == "right_end":
                    right_hit = True

    # Cleanup
    for f in (local_fa, query_fa, paf):
        f.unlink(missing_ok=True)

    return left_hit, right_hit


def _detect_tandem_repeat(seq: str, min_period: int = 500, max_period: int = 5000,
                          ) -> tuple[str, int, int]:
    """Detect and trim tandem repeats in a sequence.

    Tests candidate periods from min_period to max_period. When a period with
    ≥95% identity is found, counts how many repeat copies exist at the end
    of the sequence and trims all but one.

    Returns (trimmed_seq, repeat_period, n_copies_removed).
    """
    if len(seq) < min_period * 3:
        return seq, 0, 0

    # Test candidate periods (step=1 for exact detection; ≤5K iterations is fast)
    best_period = 0
    for period in range(min_period, min(max_period + 1, len(seq) // 3)):
        # Check identity at this period over a diagnostic window
        window = min(period, len(seq) - period * 2)
        if window < 200:
            continue
        # Sample from the latter half of the sequence
        start = len(seq) - period * 2
        matches = 0
        total = 0
        for i in range(start, start + window):
            a = seq[i]
            b = seq[i + period]
            if a != "N" and b != "N":
                total += 1
                if a == b:
                    matches += 1
        if total > 200 and matches / total >= 0.95:
            best_period = period
            break  # shortest matching period

    if best_period == 0:
        return seq, 0, 0

    # Count how many complete copies exist at the end
    n_copies = 0
    pos = len(seq)
    while pos - best_period * 2 >= 0:
        chunk1 = seq[pos - best_period : pos]
        chunk2 = seq[pos - best_period * 2 : pos - best_period]
        matches = sum(1 for a, b in zip(chunk1, chunk2) if a == b and a != "N")
        total = sum(1 for a, b in zip(chunk1, chunk2) if a != "N" and b != "N")
        if total > 100 and matches / total >= 0.90:
            n_copies += 1
            pos -= best_period
        else:
            break

    if n_copies == 0:
        return seq, best_period, 0

    trim_len = n_copies * best_period
    trimmed = seq[: len(seq) - trim_len]
    return trimmed, best_period, n_copies


def seed_extension_assembly(
    seed_seqs: list[tuple[str, str, str, str]],  # (name, tdna_seq, host_flank, flank_side)
    candidate_r1: Path,
    candidate_r2: Path,
    element_db: Path,
    workdir: Path,
    site_id: str,
    host_bam: Path | None = None,
    host_ref: Path | None = None,
    all_junctions: list | None = None,
    threads: int = 4,
    k: int = 15,
    min_overlap: int = 20,
    min_depth: int = 2,
    min_ratio: float = 0.7,
    max_iterations: int = 100,
    max_rounds: int = 50,
) -> tuple[Path | None, int, str]:
    """Alternating K-mer extension + PE-anchored Pilon gap fill.

    For each round:
      1. K-mer recruit from unmapped reads
      2. K-mer extension (strand-aware, PE-aware)
      3. PE-anchored Pilon: find R2 mates at contig ends → build
         [anchor] + 500N + [contig] + 500N + [anchor] → Pilon fills gaps
      4. Check if contig ends reached host genome → stop if both ends hit

    Returns (contig_fasta_path, n_rounds, status).
    """
    workdir.mkdir(parents=True, exist_ok=True)

    if not seed_seqs:
        return None, 0, "no_seeds"

    # Initialize strand-aware extender with junction-region reads
    extender = SeedExtender(k=k, min_overlap=min_overlap,
                            min_depth=min_depth, min_ratio=min_ratio)
    log(f"  Loading candidate reads (strand-aware PE)...")
    n_pairs = extender.load_paired_reads(candidate_r1, candidate_r2)
    log(f"  Loaded {n_pairs:,} pairs ({len(extender.seqs):,} seqs)")
    log(f"  K-mer index: {len(extender.kmer_index):,} unique {k}-mers")

    # Extract unmapped paired reads (cached)
    unmapped_r1: Path | None = None
    unmapped_r2: Path | None = None
    if host_bam is not None:
        unmapped_r1, unmapped_r2 = extract_unmapped_paired(
            host_bam, workdir, threads=threads,
        )

    # Phase 1: K-mer extension of each seed (T-DNA only, no host flank)
    contigs: list[tuple[str, str]] = []
    best_host_flank = ""
    best_flank_side = ""
    for seed_name, tdna_seq, host_flank, flank_side in seed_seqs:
        log(f"  Extending: {seed_name} ({len(tdna_seq)}bp T-DNA)")
        extended = extender.run(tdna_seq, max_iterations=max_iterations)
        contigs.append((seed_name, extended))

    # Pick best seed (longest with element_db hits)
    best_name, best_seq = _pick_best_contig(contigs, element_db, workdir, site_id)
    if not best_seq:
        return None, 0, "no_assembly"

    # Retrieve host flank for the best seed
    for seed_name, tdna_seq, host_flank, flank_side in seed_seqs:
        if seed_name == best_name:
            best_host_flank = host_flank
            best_flank_side = flank_side
            break

    log(f"  Phase 1 best: {best_name} → {len(best_seq):,}bp T-DNA"
        f" + {len(best_host_flank)}bp host ({best_flank_side})")

    # Phase 2: Alternating recruitment + extension + PE-anchored Pilon
    # 'current' tracks T-DNA portion only; host flank is added for Pilon
    current = best_seq
    total_rounds = 0
    growth_history: list[int] = []  # track per-round growth for cycle detection

    for rnd in range(1, max_rounds + 1):
        prev_len = len(current)
        log(f"\n  --- Round {rnd} ---")

        # Step 1: K-mer recruit from unmapped pool (T-DNA k-mers only)
        if unmapped_r1 is not None and unmapped_r2 is not None:
            r1_new, r2rc_new = recruit_by_kmer(
                current, unmapped_r1, unmapped_r2, k=25,
            )
            if r1_new:
                n_added = extender.add_seqs(r1_new + r2rc_new)
                log(f"    Recruit: {len(r1_new)} pairs → {n_added} new seqs "
                    f"({len(extender.seqs):,} total)")
            else:
                log(f"    Recruit: 0 pairs")

        # Step 2: K-mer extension (T-DNA only)
        current = extender.run(current, max_iterations=max_iterations)
        ext_growth = len(current) - prev_len
        log(f"    K-mer ext: {prev_len:,} → {len(current):,}bp (+{ext_growth})")

        # Step 3: PE-anchored Pilon gap fill
        # Build full scaffold: host_flank + T-DNA (or T-DNA + host_flank)
        if best_flank_side == "prefix":
            pilon_contig = best_host_flank + current
        elif best_flank_side == "suffix":
            pilon_contig = current + best_host_flank
        else:
            pilon_contig = current

        pilon_dir = workdir / f"_pilon_r{rnd}"
        pilon_dir.mkdir(parents=True, exist_ok=True)
        pool_r1 = pilon_dir / "pool_R1.fq"
        pool_r2 = pilon_dir / "pool_R2.fq"
        _write_pool_fastq(extender, pool_r1, pool_r2)

        prev_pilon_full = len(pilon_contig)
        pilon_result = pilon_fill(
            pilon_contig, extender,
            pool_r1, pool_r2,
            pilon_dir, insert_size=500, threads=threads,
        )

        # Strip host flank back out to get T-DNA portion
        host_len = len(best_host_flank)
        if best_flank_side == "prefix" and host_len > 0:
            current = pilon_result[host_len:]
        elif best_flank_side == "suffix" and host_len > 0:
            current = pilon_result[:-host_len] if host_len > 0 else pilon_result
        else:
            current = pilon_result

        pilon_growth = len(current) - prev_len
        log(f"    Pilon: {prev_len:,} → {len(current):,}bp (+{pilon_growth})")

        shutil.rmtree(pilon_dir, ignore_errors=True)

        total_rounds = rnd
        total_growth = ext_growth + pilon_growth
        if total_growth == 0:
            log(f"    No growth → converged")
            break

        # Cycle detection (two levels):
        # 1) Fast: 2 consecutive rounds with identical growth → immediate stop
        # 2) Pattern: 3-round growth pattern repeats → stop
        growth_history.append(total_growth)
        cycle_detected = False
        if (len(growth_history) >= 2
                and growth_history[-1] == growth_history[-2]
                and growth_history[-1] > 0):
            log(f"    Constant growth ({growth_history[-1]}) for 2 rounds → cycling")
            cycle_detected = True
        elif len(growth_history) >= 6:
            recent = growth_history[-3:]
            earlier = growth_history[-6:-3]
            if recent == earlier:
                log(f"    Cyclic growth detected "
                    f"(pattern {recent} repeated) → stopping")
                cycle_detected = True

        if cycle_detected:
            # Trim tandem repeats from the assembled T-DNA
            trimmed, period, n_removed = _detect_tandem_repeat(current)
            if n_removed > 0:
                log(f"    Tandem repeat: period={period}bp, "
                    f"removed {n_removed} copies "
                    f"({len(current):,} → {len(trimmed):,}bp)")
                current = trimmed
            break

        # Step 4: Check if T-DNA far end reached host genome
        # Uses minimap2 against local host region (±50kb around junction)
        # Only check the FAR end (near end always has host from junction leak)
        if all_junctions and host_ref and len(current) >= 1000:
            left_hit, right_hit = _check_host_reached(
                current, host_ref, all_junctions, workdir,
                check_len=500,
            )
            if best_flank_side == "prefix":
                far_end_hit = right_hit
            elif best_flank_side == "suffix":
                far_end_hit = left_hit
            else:
                far_end_hit = left_hit and right_hit
            if far_end_hit:
                log(f"    >>> Far end reached host genome — insert complete!")
                break

    # Write final insert (T-DNA + host flank)
    if best_flank_side == "prefix":
        final_seq = best_host_flank + current
    elif best_flank_side == "suffix":
        final_seq = current + best_host_flank
    else:
        final_seq = current

    final_fa = workdir / f"{site_id}_insert.fasta"
    write_fasta(final_fa, f"{site_id}_assembled_insert", final_seq)
    log(f"  Final: {len(current):,}bp T-DNA + {len(best_host_flank)}bp host "
        f"= {len(final_seq):,}bp after {total_rounds} rounds")

    status = "complete" if len(current) > max(len(s) for _, s, *_ in seed_seqs) else "partial"
    return final_fa, total_rounds, status


def _pick_best_contig(
    contigs: list[tuple[str, str]],
    element_db: Path,
    workdir: Path,
    site_id: str,
) -> tuple[str, str]:
    """BLAST contigs vs element_db, pick longest with hits."""
    contigs_fa = workdir / f"{site_id}_contigs.fasta"
    with open(contigs_fa, "w") as fh:
        for name, seq in contigs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

    blast_out = workdir / f"{site_id}_blast.tsv"
    subprocess.run(
        ["blastn", "-query", str(contigs_fa), "-subject", str(element_db),
         "-outfmt", "6 qseqid sseqid pident length qstart qend",
         "-evalue", "1e-5", "-max_target_seqs", "10",
         "-out", str(blast_out)],
        stderr=subprocess.DEVNULL, check=True,
    )

    hit_names: set[str] = set()
    if blast_out.exists():
        with open(blast_out) as fh:
            for line in fh:
                hit_names.add(line.split("\t")[0])

    best_name, best_seq = "", ""
    for name, seq in contigs:
        if name in hit_names and len(seq) > len(best_seq):
            best_name, best_seq = name, seq
    if not best_seq:
        for name, seq in contigs:
            if len(seq) > len(best_seq):
                best_name, best_seq = name, seq

    return best_name, best_seq


# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------

def annotate_insert(
    insert_fasta: Path,
    element_db: Path,
    output_dir: Path,
    sample_name: str,
) -> tuple[Path, Path]:
    """Annotate insert with BLAST vs element_db + border motif search."""
    annotation_tsv = output_dir / "element_annotation.tsv"
    border_tsv = output_dir / "border_hits.tsv"

    # BLAST vs element_db
    # First make a blast DB if needed
    db_prefix = output_dir / "element_blastdb"
    subprocess.run(
        ["makeblastdb", "-in", str(element_db), "-dbtype", "nucl",
         "-out", str(db_prefix)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        check=True,
    )

    blast_out = output_dir / "blast_raw.tsv"
    subprocess.run(
        ["blastn", "-query", str(insert_fasta), "-db", str(db_prefix),
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
         "-evalue", "1e-5", "-max_target_seqs", "50",
         "-out", str(blast_out)],
        stderr=subprocess.DEVNULL,
        check=True,
    )

    # Parse and write annotation
    with open(blast_out) as fin, open(annotation_tsv, "w") as fout:
        fout.write("query\telement\tidentity\tlength\tq_start\tq_end\ts_start\ts_end\tevalue\n")
        hits = []
        for line in fin:
            cols = line.rstrip().split("\t")
            if len(cols) >= 10:
                hits.append(cols)
                fout.write(f"{cols[0]}\t{cols[1]}\t{cols[2]}\t{cols[3]}\t"
                           f"{cols[4]}\t{cols[5]}\t{cols[6]}\t{cols[7]}\t{cols[8]}\n")

    if hits:
        log(f"  BLAST hits: {len(hits)} elements found")
        # Sort by query position and show top hits
        hits.sort(key=lambda x: int(x[4]))
        for h in hits[:10]:
            log(f"    {h[4]}-{h[5]}: {h[1]} ({h[2]}% identity, {h[3]}bp)")
    else:
        log("  No BLAST hits found against element database")

    # T-DNA border motif search
    # Canonical border: TGGCAGGATATATTGTGGTGTAAAC (25bp)
    # Allow mismatches via short word BLAST
    border_fa = output_dir / "_borders.fa"
    with open(border_fa, "w") as fh:
        fh.write(">RB_consensus\nTGGCAGGATATATTGTGGTGTAAAC\n")
        fh.write(">LB_consensus\nTGGCAGGATATATTGTGGTGTAAAC\n")

    subprocess.run(
        ["blastn", "-query", str(border_fa), "-subject", str(insert_fasta),
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send",
         "-evalue", "1", "-word_size", "7",
         "-out", str(border_tsv)],
        stderr=subprocess.DEVNULL,
    )

    if border_tsv.stat().st_size > 0:
        log("  T-DNA border motifs found")
    else:
        log("  No border motifs found (may need manual inspection)")

    # Cleanup
    border_fa.unlink(missing_ok=True)
    blast_out.unlink(missing_ok=True)
    for ext in [".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto"]:
        p = Path(f"{db_prefix}{ext}")
        p.unlink(missing_ok=True)

    return annotation_tsv, border_tsv


# ---------------------------------------------------------------------------
# Stats output
# ---------------------------------------------------------------------------

def write_stats(
    stats_path: Path,
    sample_name: str,
    num_sites: int,
    candidate_reads: int,
    results: list[dict],
) -> None:
    """Write assembly statistics."""
    with open(stats_path, "w") as fh:
        fh.write(f"sample\t{sample_name}\n")
        fh.write(f"insertion_sites\t{num_sites}\n")
        fh.write(f"candidate_read_pairs\t{candidate_reads}\n")
        for r in results:
            prefix = r["site_id"]
            fh.write(f"{prefix}_insert_length\t{r['insert_length']}\n")
            fh.write(f"{prefix}_remaining_ns\t{r['remaining_ns']}\n")
            fh.write(f"{prefix}_assembly_rounds\t{r['rounds']}\n")
            fh.write(f"{prefix}_status\t{r['status']}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 9: Targeted seed extension insert assembly")
    parser.add_argument("--junctions", required=True,
                        help="junctions.tsv from step 6")
    parser.add_argument("--junction-contigs", default=None,
                        help="(unused, kept for CLI compat)")
    parser.add_argument("--host-paf", default=None,
                        help="(unused, kept for CLI compat)")
    parser.add_argument("--construct-paf", default=None,
                        help="(unused, kept for CLI compat)")
    parser.add_argument("--host-bam", required=True,
                        help="Host-mapped BAM from step 7")
    parser.add_argument("--host-ref", default=None,
                        help="Host reference FASTA (for junction-stop detection)")
    parser.add_argument("--element-db", required=True,
                        help="Element database FASTA for annotation")
    parser.add_argument("--construct-ref", default=None,
                        help="(unused, kept for CLI compat)")
    parser.add_argument("--s03-r1", default=None,
                        help="(unused, kept for CLI compat)")
    parser.add_argument("--s03-r2", default=None,
                        help="(unused, kept for CLI compat)")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--max-rounds", type=int, default=50,
                        help="Max alternating extension+Pilon rounds")
    parser.add_argument("--seed-k", type=int, default=15,
                        help="k-mer size for seed extension")
    parser.add_argument("--min-overlap", type=int, default=20)
    parser.add_argument("--min-depth", type=int, default=2)
    parser.add_argument("--min-ratio", type=float, default=0.7)
    parser.add_argument("--flank", type=int, default=1000,
                        help="Flank size for candidate read extraction (bp)")
    args = parser.parse_args()

    step_dir = Path(args.outdir) / args.sample_name / STEP
    step_dir.mkdir(parents=True, exist_ok=True)

    log(f"=== Step 9: Targeted Seed Extension Assembly for {args.sample_name} ===")

    # ---- Parse junctions ----
    junctions_path = Path(args.junctions)
    if not junctions_path.exists() or junctions_path.stat().st_size == 0:
        log("No junctions found (step 6 empty). Nothing to assemble.")
        write_stats(step_dir / "s09_stats.txt", args.sample_name, 0, 0, [])
        return

    junctions = parse_junctions(junctions_path)
    if not junctions:
        log("No junctions parsed. Exiting.")
        write_stats(step_dir / "s09_stats.txt", args.sample_name, 0, 0, [])
        return

    log(f"Parsed {len(junctions)} junctions")

    # ---- Group into insertion sites ----
    sites = group_insertion_sites(junctions)
    log(f"Grouped into {len(sites)} insertion site(s)")

    for site in sites:
        n_junc = len(site.junctions)
        paired = "paired" if site.pos_3p > 0 else "single"
        log(f"  {site.site_id}: {site.host_chr}:{site.pos_5p}"
            f"{f'-{site.pos_3p}' if site.pos_3p else ''}"
            f" ({n_junc} junction(s), {paired})")

    # ---- Detect cross-chromosome false positives ----
    contig_chrs: dict[str, set[str]] = {}
    for site in sites:
        for j in site.junctions:
            contig_chrs.setdefault(j.contig_name, set()).add(site.host_chr)

    cross_chr_contigs = {c for c, chrs in contig_chrs.items() if len(chrs) > 1}

    skip_sites: set[str] = set()
    if cross_chr_contigs:
        chr_scores: dict[str, int] = {}
        for site in sites:
            score = len(site.junctions) + (2 if site.pos_3p > 0 else 0)
            chr_scores[site.host_chr] = chr_scores.get(site.host_chr, 0) + score
        primary_chr = max(chr_scores, key=chr_scores.get)

        for site in sites:
            site_contigs = {j.contig_name for j in site.junctions}
            if (site_contigs & cross_chr_contigs
                    and site.host_chr != primary_chr):
                skip_sites.add(site.site_id)
                log(f"  Skipping {site.site_id} ({site.host_chr}:{site.pos_5p})"
                    f" — cross-chr false positive (contig also maps to {primary_chr})")

    # ---- Collect all junctions for host-reached detection ----
    all_junctions_for_host: list[Junction] = []
    if args.host_ref:
        all_junctions_for_host = [
            j for s in sites for j in s.junctions
            if s.site_id not in skip_sites
        ]
        if all_junctions_for_host:
            log("Extracting host flanks at junction positions...")
            extract_junction_host_flanks(
                Path(args.host_ref), all_junctions_for_host, flank_len=200,
            )

    # ---- Process each insertion site ----
    all_results = []
    for site in sites:
        log(f"\n=== Processing {site.site_id} ===")

        if site.site_id in skip_sites:
            all_results.append({
                "site_id": site.site_id,
                "insert_length": 0,
                "remaining_ns": 0,
                "rounds": 0,
                "status": "skipped_cross_chr",
            })
            continue

        # Step 1: Extract junction seeds (T-DNA soft-clips + 1kb host flank)
        seed_seqs = extract_junction_seeds(
            Path(args.host_bam), site.junctions,
            host_ref=Path(args.host_ref) if args.host_ref else None,
        )

        # Step 2: Extract candidate reads from junction region
        cand_r1 = step_dir / f"{site.site_id}_candidate_R1.fastq.gz"
        cand_r2 = step_dir / f"{site.site_id}_candidate_R2.fastq.gz"
        if not cand_r1.exists():
            n_cand = extract_junction_candidates(
                Path(args.host_bam), site.junctions,
                cand_r1, cand_r2,
                flank=args.flank, threads=args.threads,
            )
            log(f"  Candidate reads: {n_cand:,} pairs")
        else:
            n_cand = 0
            with pysam.FastxFile(str(cand_r1)) as fh:
                for _ in fh:
                    n_cand += 1
            log(f"  Candidate reads (cached): {n_cand:,} pairs")

        # Step 3: Alternating K-mer extension + PE-anchored Pilon
        assembly_fa, rounds, status = seed_extension_assembly(
            seed_seqs=seed_seqs,
            candidate_r1=cand_r1,
            candidate_r2=cand_r2,
            element_db=Path(args.element_db),
            workdir=step_dir,
            site_id=site.site_id,
            host_bam=Path(args.host_bam),
            host_ref=Path(args.host_ref) if args.host_ref else None,
            all_junctions=all_junctions_for_host or None,
            threads=args.threads,
            k=args.seed_k,
            min_overlap=args.min_overlap,
            min_depth=args.min_depth,
            min_ratio=args.min_ratio,
            max_iterations=100,
            max_rounds=args.max_rounds,
        )

        # Record result
        if assembly_fa and assembly_fa.exists():
            contigs = read_fasta(assembly_fa)
            if contigs:
                longest_name = max(contigs, key=lambda k: len(contigs[k]))
                longest_seq = contigs[longest_name]
                insert_len = len(longest_seq)
                insert_ns = longest_seq.upper().count("N")
                log(f"  Insert: {insert_len:,}bp → {status}")
            else:
                insert_len, insert_ns = 0, 0
                status = "no_assembly"
        else:
            insert_len, insert_ns = 0, 0
            status = "no_assembly"

        all_results.append({
            "site_id": site.site_id,
            "insert_length": insert_len,
            "remaining_ns": insert_ns,
            "rounds": rounds,
            "status": status,
        })

    # ---- Annotation ----
    combined_insert = step_dir / "insert_only.fasta"
    with open(combined_insert, "w") as fout:
        for site in sites:
            site_fa = step_dir / f"{site.site_id}_insert.fasta"
            if site_fa.exists():
                fout.write(open(site_fa).read())

    if combined_insert.exists() and combined_insert.stat().st_size > 0:
        log("\n=== Annotating insert ===")
        annotate_insert(
            combined_insert, Path(args.element_db),
            step_dir, args.sample_name,
        )

    # ---- Write stats ----
    write_stats(
        step_dir / "s09_stats.txt", args.sample_name,
        len(sites), 0, all_results,
    )

    log(f"\n=== Step 9 complete ===")
    log(f"Output: {step_dir}")
    for r in all_results:
        log(f"  {r['site_id']}: {r['insert_length']}bp, "
            f"{r['remaining_ns']} Ns, {r['rounds']} rounds → {r['status']}")


if __name__ == "__main__":
    main()
