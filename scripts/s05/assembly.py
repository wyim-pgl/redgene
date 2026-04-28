"""Phase 3 — Assembly (Issue #4, Session 4).

Extracted verbatim from ``scripts/s05_insert_assembly.py``. Owns the
k-mer extension + Pilon gap fill + minimap2/SSAKE refinement pipeline
that turns soft-clip seeds into a contiguous insert assembly.

Pipeline overview
-----------------
``assemble_insert`` is the public entry point invoked once per
``InsertionSite`` after Phase 1.5 classify. It orchestrates:

1. ``StrandAwareSeedExtender`` — k-mer index over candidate reads
2. ``recruit_by_kmer`` — pull reads matching seed k-mers
3. ``_minimap2_extend`` + ``_vote_extension`` — extend contigs by tail voting
4. ``_ssake_extend`` — fallback assembler for stalled extensions
5. ``pilon_fill`` — gap-fill via Pilon polishing
6. ``_check_merge`` — try to merge 5'/3' contigs into a single insert
7. ``check_host_termination`` — confirm host flanking reached
8. ``extract_foreign_reads`` + ``refine_with_foreign_reads`` — recruit
   the read pool relevant to the assembled insert and re-polish

This module is the highest-risk extraction in the s05 split: it owns
the core assembly logic that produces the FASTA consumed by Phase 4
annotation and the four post-assembly FP filters. Verbatim move +
``tests/fixtures/verdict_snapshots/`` are the primary regression
signal. The 3-host 15x e2e gate from the design spec is deferred
(hours of HPC compute); snapshot fixtures cover the verdict-level
surface that downstream consumers care about.

Deviation from verbatim policy (documented):
    ``assemble_insert`` calls ``extract_unmapped_paired`` which is still
    defined in the monolith (Phase 2 territory, Session 5). A module-top
    import would create a circular import during initialization (monolith
    imports assembly at top; assembly tries to import monolith before
    extract_unmapped_paired is defined). The fix: lazy import inside
    ``assemble_insert`` body — identical to Session 3's ``_parse_src_tag``
    cycle resolution. The call site is byte-identical; only the import
    location changes.
"""
from __future__ import annotations

import gzip
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import pysam

from .primitives import (
    log,
    revcomp,
    read_fasta,
    write_fasta,
    _read_fq_seqs,
    InsertionSite,
)


# ---------------------------------------------------------------------------
# Phase 3: K-mer Extension + Pilon Gap Fill
# ---------------------------------------------------------------------------

class StrandAwareSeedExtender:
    """K-mer based extension engine with PE strand awareness.

    INDEX RULES:
    - R1 sequences: indexed as-is (forward strand)
    - R2 sequences: reverse-complemented then indexed (fragment strand)
    - Raw R2 (as sequenced): NEVER indexed, NEVER used for extension

    EXTENSION RULES:
    - Track used reads to prevent infinite loops
    - Stop at branch points (base ratio < min_ratio)
    - Majority vote consensus for each new base
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
        self.seqs: list[str] = []
        self.kmer_index: dict[str, list[tuple[int, int]]] = defaultdict(list)

    def load_paired_reads(self, r1_path: Path, r2_path: Path) -> int:
        """Load R1 as forward, R2 as reverse-complement. Return pair count."""
        r1_seqs = _read_fq_seqs(r1_path)
        r2_seqs = _read_fq_seqs(r2_path)
        k = self.k
        n = 0
        for r1, r2 in zip(r1_seqs, r2_seqs):
            if len(r1) < self.min_overlap or len(r2) < self.min_overlap:
                continue
            if "N" in r1 or "N" in r2:
                continue
            r2_rc = revcomp(r2)

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
        """Extend seed to the right with used-read tracking."""
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

    def extend(self, seed: str, max_iterations: int = 100) -> str:
        """Extend seed bidirectionally. Return extended sequence."""
        current = seed
        used: set[int] = set()
        for i in range(max_iterations):
            # Right extension
            extended = self._extend_right(current, used)
            # Left extension (revcomp -> extend right -> revcomp)
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


# ---------------------------------------------------------------------------
# K-mer recruitment
# ---------------------------------------------------------------------------

def recruit_by_kmer(
    contig_seq: str,
    unmapped_r1: Path,
    unmapped_r2: Path,
    k: int = 25,
) -> tuple[list[str], list[str]]:
    """K-mer based paired recruitment from unmapped reads.

    Builds k-mer set from contig (both strands), scans R1 forward + R2 RC.
    Raw R2 is NOT used (head-to-head safety).
    Returns (r1_seqs, r2_rc_seqs).
    """
    contig_kmers: set[str] = set()
    for seq in (contig_seq, revcomp(contig_seq)):
        for i in range(len(seq) - k + 1):
            contig_kmers.add(seq[i:i + k])

    r1_seqs = _read_fq_seqs(unmapped_r1)
    r2_seqs = _read_fq_seqs(unmapped_r2)

    recruited_r1: list[str] = []
    recruited_r2rc: list[str] = []
    stride = max(1, k // 2)

    for r1, r2 in zip(r1_seqs, r2_seqs):
        if "N" in r1 or "N" in r2:
            continue
        r2_rc = revcomp(r2)

        hit = False
        for seq in (r1, r2_rc):
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


# ---------------------------------------------------------------------------
# Pilon gap fill (CRITICAL: [contig_5p] + NNN + [contig_3p] as ONE scaffold)
# ---------------------------------------------------------------------------

def pilon_fill(
    contig_5p: str,
    contig_3p: str,
    r1_fq: Path,
    r2_fq: Path,
    workdir: Path,
    gap_size: int = 1000,
    threads: int = 4,
) -> tuple[str, str, bool]:
    """Build scaffold [contig_5p]+NNN+[contig_3p] and run Pilon to fill gap.

    KEY INSIGHT: Pilon only fills INTERNAL gaps between two known sequences.
    Both contigs must be in ONE scaffold for Pilon to work.

    Returns (updated_5p, updated_3p, gap_filled).
    gap_filled=True means no N remains (complete assembly).
    """
    workdir.mkdir(parents=True, exist_ok=True)

    if not contig_5p or not contig_3p:
        return contig_5p, contig_3p, False

    # Build scaffold: [contig_5p] + NNN + [contig_3p]
    scaffold = contig_5p + ("N" * gap_size) + contig_3p
    ref_fa = workdir / "pilon_ref.fasta"
    write_fasta(ref_fa, "scaffold", scaffold)

    # Map reads with minimap2
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
        return contig_5p, contig_3p, False

    log(f"    Pilon: {n_mapped:,} reads mapped to scaffold "
        f"({len(contig_5p)}+{gap_size}N+{len(contig_3p)}={len(scaffold)}bp)")

    # Run Pilon
    pilon_prefix = workdir / "pilon_out"
    pilon_log_file = workdir / "pilon.log"
    with open(pilon_log_file, "w") as logfh:
        subprocess.run(
            ["pilon", "--genome", str(ref_fa), "--frags", str(bam),
             "--output", str(pilon_prefix), "--fix", "all",
             "--mindepth", "2", "--gapmargin", "100000"],
            stdout=logfh, stderr=subprocess.STDOUT,
        )

    if pilon_log_file.exists():
        with open(pilon_log_file) as fh:
            for line in fh:
                line = line.strip()
                if any(k in line for k in ["Gap", "fix", "Total", "Confirmed"]):
                    log(f"    [Pilon] {line}")

    pilon_fa = Path(f"{pilon_prefix}.fasta")
    if not pilon_fa.exists():
        return contig_5p, contig_3p, False

    seqs = read_fasta(pilon_fa)
    if not seqs:
        return contig_5p, contig_3p, False

    filled = list(seqs.values())[0]

    # Check if gap is fully filled (no N remaining)
    n_match = re.search(r"N{10,}", filled)
    if n_match is None:
        # Gap fully filled!
        log(f"    Pilon: gap completely filled! ({len(filled):,}bp)")
        return filled, "", True

    # Split at longest remaining N stretch
    n_runs = [(m.start(), m.end()) for m in re.finditer(r"N+", filled)]
    if not n_runs:
        return filled, "", True

    # Find longest N run
    longest = max(n_runs, key=lambda x: x[1] - x[0])
    new_5p = filled[:longest[0]]
    new_3p = filled[longest[1]:]

    growth_5p = len(new_5p) - len(contig_5p)
    growth_3p = len(new_3p) - len(contig_3p)
    remaining_gap = longest[1] - longest[0]
    log(f"    Pilon: 5p {len(contig_5p):,}→{len(new_5p):,}bp (+{growth_5p}), "
        f"3p {len(contig_3p):,}→{len(new_3p):,}bp (+{growth_3p}), "
        f"gap {gap_size}→{remaining_gap}N")

    return new_5p, new_3p, False


# ---------------------------------------------------------------------------
# Host genome termination check
# ---------------------------------------------------------------------------

def check_host_termination(
    contig_5p: str,
    contig_3p: str,
    host_ref: Path,
    workdir: Path,
    check_len: int = 100,
    min_identity: float = 0.95,
    min_match: int = 50,
) -> tuple[bool, bool]:
    """Check if contig ends have reached back to host genome.

    5' contig growing end = right end (extending into insert then through to host)
    3' contig growing end = left end (extending into insert then through to host)

    Returns (reached_5p, reached_3p).
    """
    if len(contig_5p) < check_len * 2 and len(contig_3p) < check_len * 2:
        return False, False

    query_fa = workdir / "_host_term.fa"
    with open(query_fa, "w") as fh:
        if len(contig_5p) >= check_len * 2:
            fh.write(f">growing_5p\n{contig_5p[-check_len:]}\n")
        if len(contig_3p) >= check_len * 2:
            fh.write(f">growing_3p\n{contig_3p[:check_len]}\n")

    paf = workdir / "_host_term.paf"
    subprocess.run(
        ["minimap2", "-c", "--secondary=no", "-t", "1",
         str(host_ref), str(query_fa)],
        stdout=open(paf, "w"), stderr=subprocess.DEVNULL,
    )

    reached_5p = False
    reached_3p = False
    if paf.exists():
        with open(paf) as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 12:
                    continue
                qname = cols[0]
                match_bp = int(cols[9])
                block_len = int(cols[10])
                if match_bp < min_match or block_len < min_match:
                    continue
                identity = match_bp / block_len
                if identity < min_identity:
                    continue
                if qname == "growing_5p":
                    reached_5p = True
                elif qname == "growing_3p":
                    reached_3p = True

    for f in (query_fa, paf):
        f.unlink(missing_ok=True)

    return reached_5p, reached_3p


# ---------------------------------------------------------------------------
# Foreign read refinement (minimap2 chaining + Pilon)
# ---------------------------------------------------------------------------

def extract_foreign_reads(
    s03_r1: Path,
    s03_r2: Path,
    host_ref: Path,
    workdir: Path,
    threads: int = 4,
) -> tuple[Path, Path]:
    """Extract s03 reads that DON'T map to host genome (construct-only reads).

    These 'foreign' reads contain transgene sequences not present in the host.
    By filtering out host-mapping reads, we enrich for construct-internal
    elements (hLF1, G6 EPSPS, Pepc, etc.) that enable Pilon to extend
    assemblies through regions the k-mer assembler couldn't resolve.
    """
    foreign_r1 = workdir / "_foreign_R1.fq.gz"
    foreign_r2 = workdir / "_foreign_R2.fq.gz"
    if foreign_r1.exists() and foreign_r2.exists():
        n = 0
        with gzip.open(foreign_r1, "rt") as fh:
            for _ in fh:
                n += 1
        log(f"  Foreign reads (cached): {n // 4:,} pairs")
        return foreign_r1, foreign_r2

    log(f"  Extracting foreign reads (s03 reads unmapped to host)...")

    # Map s03 reads to host
    host_bam = workdir / "_s03_to_host.bam"
    subprocess.run(
        f"bwa mem -t {threads} {host_ref} {s03_r1} {s03_r2} 2>/dev/null "
        f"| samtools sort -@ 4 -o {host_bam}",
        shell=True, check=True,
    )
    subprocess.run(["samtools", "index", str(host_bam)],
                   stderr=subprocess.DEVNULL, check=True)

    # Get properly mapped read names (MAPQ >= 20, both reads mapped)
    host_mapped = set()
    all_reads = set()
    with pysam.AlignmentFile(str(host_bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            all_reads.add(read.query_name)
            if (not read.is_unmapped and not read.mate_is_unmapped
                    and read.mapping_quality >= 20):
                host_mapped.add(read.query_name)

    foreign_names = all_reads - host_mapped
    log(f"  Total: {len(all_reads):,} pairs, "
        f"host: {len(host_mapped):,}, foreign: {len(foreign_names):,}")

    # Extract foreign reads as FASTQ
    readname_file = workdir / "_foreign_names.txt"
    with open(readname_file, "w") as fh:
        for name in foreign_names:
            fh.write(f"{name}\n")

    subprocess.run(
        f"samtools view -h -N {readname_file} {host_bam} "
        f"| samtools sort -n -@ 4 - "
        f"| samtools fastq -1 {foreign_r1} -2 {foreign_r2} "
        f"-s /dev/null -0 /dev/null - 2>/dev/null",
        shell=True, check=True,
    )

    # Cleanup
    host_bam.unlink(missing_ok=True)
    Path(f"{host_bam}.bai").unlink(missing_ok=True)
    readname_file.unlink(missing_ok=True)

    n = 0
    with gzip.open(foreign_r1, "rt") as fh:
        for _ in fh:
            n += 1
    log(f"  Foreign reads extracted: {n // 4:,} pairs")
    return foreign_r1, foreign_r2


def refine_with_foreign_reads(
    insert_fasta: Path,
    s03_r1: Path,
    s03_r2: Path,
    host_ref: Path,
    workdir: Path,
    threads: int = 4,
    max_rounds: int = 10,
) -> Path:
    """Refine assembled insert using minimap2 chaining + Pilon with foreign reads.

    Approach:
      1. Extract foreign reads (s03 reads unmapped to host)
      2. Iterative: minimap2 map → Pilon --fix all → check convergence
      3. Pilon's local reassembly can extend construct regions using
         foreign reads that the k-mer assembler couldn't resolve
         (especially palindromic/inverted repeat regions in head-to-head T-DNA)

    Returns path to refined FASTA.
    """
    refine_dir = workdir / "_foreign_refine"
    refine_dir.mkdir(parents=True, exist_ok=True)

    # Extract foreign reads
    foreign_r1, foreign_r2 = extract_foreign_reads(
        s03_r1, s03_r2, host_ref, refine_dir, threads=threads,
    )

    # Check if we have enough foreign reads
    n_foreign = 0
    with gzip.open(foreign_r1, "rt") as fh:
        for _ in fh:
            n_foreign += 1
    n_foreign //= 4
    if n_foreign < 10:
        log(f"  Too few foreign reads ({n_foreign}), skipping refinement")
        return insert_fasta

    # Start with current assembly
    seqs = read_fasta(insert_fasta)
    if not seqs:
        return insert_fasta
    seq_name = list(seqs.keys())[0]
    current_seq = list(seqs.values())[0]
    prev_len = len(current_seq)

    log(f"  Foreign read refinement: {prev_len:,}bp scaffold, "
        f"{n_foreign:,} foreign read pairs")

    for rnd in range(1, max_rounds + 1):
        scaffold_fa = refine_dir / f"scaffold_r{rnd}.fa"
        write_fasta(scaffold_fa, "insert_scaffold", current_seq)

        # Map foreign reads with minimap2 (k-mer chaining handles repeats)
        bam = refine_dir / f"r{rnd}.bam"
        subprocess.run(
            f"minimap2 -ax sr -t {threads} {scaffold_fa} "
            f"{foreign_r1} {foreign_r2} 2>/dev/null "
            f"| samtools sort -@ 4 -o {bam}",
            shell=True, check=True,
        )
        subprocess.run(["samtools", "index", str(bam)],
                       stderr=subprocess.DEVNULL, check=True)

        # Count mapped
        r = subprocess.run(
            ["samtools", "view", "-c", "-F", "4", str(bam)],
            stdout=subprocess.PIPE, text=True, stderr=subprocess.DEVNULL,
        )
        n_mapped = int(r.stdout.strip()) if r.stdout.strip() else 0

        if n_mapped < 5:
            log(f"    Round {rnd}: only {n_mapped} mapped, stopping")
            break

        # Run Pilon with --fix all
        pilon_prefix = refine_dir / f"pilon_r{rnd}"
        pilon_log = refine_dir / f"pilon_r{rnd}.log"
        with open(pilon_log, "w") as logfh:
            subprocess.run(
                ["pilon", "--genome", str(scaffold_fa),
                 "--frags", str(bam),
                 "--output", str(pilon_prefix),
                 "--fix", "all", "--mindepth", "1"],
                stdout=logfh, stderr=subprocess.STDOUT,
            )

        pilon_fa = Path(f"{pilon_prefix}.fasta")
        if not pilon_fa.exists():
            log(f"    Round {rnd}: Pilon failed")
            break

        new_seqs = read_fasta(pilon_fa)
        if not new_seqs:
            break

        new_seq = list(new_seqs.values())[0]
        new_len = len(new_seq)
        n_ns = new_seq.upper().count("N")

        log(f"    Round {rnd}: {prev_len:,}→{new_len:,}bp "
            f"(Ns={n_ns}, mapped={n_mapped})")

        # Check convergence
        if new_len == prev_len and new_seq == current_seq:
            log(f"    Converged at round {rnd}")
            break

        current_seq = new_seq
        prev_len = new_len

        # Cleanup intermediate files
        scaffold_fa.unlink(missing_ok=True)
        bam.unlink(missing_ok=True)
        Path(f"{bam}.bai").unlink(missing_ok=True)
        pilon_fa.unlink(missing_ok=True)
        pilon_log.unlink(missing_ok=True)

    # Write refined result
    refined_fa = workdir / insert_fasta.name
    write_fasta(refined_fa, seq_name, current_seq)

    orig_seqs = read_fasta(insert_fasta)
    orig_len = len(list(orig_seqs.values())[0])
    final_len = len(current_seq)
    final_ns = current_seq.upper().count("N")
    log(f"  Refinement: {orig_len:,}→{final_len:,}bp, "
        f"Ns: {list(orig_seqs.values())[0].upper().count('N')}→{final_ns}")

    # Cleanup refine directory
    shutil.rmtree(refine_dir, ignore_errors=True)

    # Update the original file
    write_fasta(insert_fasta, seq_name, current_seq)
    return insert_fasta


# ---------------------------------------------------------------------------
# Overlap / merge detection
# ---------------------------------------------------------------------------

def _minimap2_extend(
    contig: str,
    r1_fq: Path,
    r2_fq: Path,
    workdir: Path,
    threads: int = 4,
    min_clip: int = 20,
    min_depth: int = 2,
    min_ratio: float = 0.6,
) -> str:
    """Extend contig using minimap2 soft-clip consensus.

    Maps reads to contig, collects soft-clipped tails beyond contig ends,
    and extends by majority vote. Handles repeats that k-mer extension
    cannot resolve because minimap2 uses minimizer chaining with gap
    penalties for placement, not exact k-mer lookup.

    Returns extended contig (may be unchanged if no extension possible).
    """
    if len(contig) < 100:
        return contig

    ref_fa = workdir / "_mm2ext_ref.fa"
    write_fasta(ref_fa, "contig", contig)

    bam_path = workdir / "_mm2ext.bam"
    subprocess.run(
        f"minimap2 -ax sr -t {threads} --secondary=no "
        f"{ref_fa} {r1_fq} {r2_fq} 2>/dev/null "
        f"| samtools sort -@ 2 -o {bam_path}",
        shell=True, check=True,
    )
    subprocess.run(["samtools", "index", str(bam_path)],
                   stderr=subprocess.DEVNULL, check=True)

    contig_len = len(contig)

    # Collect right-side soft-clip extensions (reads extending past contig end)
    right_tails: list[str] = []
    # Collect left-side soft-clip extensions (reads extending before contig start)
    left_tails: list[str] = []

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            cigar = read.cigartuples
            if not cigar:
                continue
            seq = read.query_sequence
            if not seq:
                continue

            # Right extension: read maps near contig end, has right soft-clip
            # cigar[-1] == (4, clip_len) means soft-clip at right end
            if (cigar[-1][0] == 4 and cigar[-1][1] >= min_clip
                    and read.reference_end is not None
                    and read.reference_end >= contig_len - 10):
                clip_len = cigar[-1][1]
                tail = seq[-clip_len:]
                right_tails.append(tail)

            # Left extension: read maps near contig start, has left soft-clip
            # cigar[0] == (4, clip_len) means soft-clip at left end
            if (cigar[0][0] == 4 and cigar[0][1] >= min_clip
                    and read.reference_start is not None
                    and read.reference_start <= 10):
                clip_len = cigar[0][1]
                tail = seq[:clip_len]
                left_tails.append(tail)

    # Cleanup
    for f in [ref_fa, bam_path, Path(f"{bam_path}.bai")]:
        f.unlink(missing_ok=True)

    extended = contig

    # Right extension via majority vote on soft-clip tails
    if len(right_tails) >= min_depth:
        ext = _vote_extension(right_tails, min_depth, min_ratio)
        if ext:
            extended = extended + ext

    # Left extension via majority vote on soft-clip tails (reverse direction)
    if len(left_tails) >= min_depth:
        # Reverse tails so we vote from the clip boundary outward
        rev_tails = [t[::-1] for t in left_tails]
        ext = _vote_extension(rev_tails, min_depth, min_ratio)
        if ext:
            extended = ext[::-1] + extended

    return extended


def _vote_extension(
    tails: list[str], min_depth: int = 2, min_ratio: float = 0.6,
) -> str:
    """Majority-vote consensus from a list of tail sequences.

    Each tail starts at the extension boundary and grows outward.
    Returns the consensus extension string (may be empty).
    """
    ext_bases: list[str] = []
    pos = 0
    while True:
        votes: Counter = Counter()
        for t in tails:
            if pos < len(t):
                votes[t[pos]] += 1
        total = sum(votes.values())
        if total < min_depth:
            break
        best_base, best_count = votes.most_common(1)[0]
        if best_count / total < min_ratio:
            break
        ext_bases.append(best_base)
        pos += 1
    return "".join(ext_bases)


def _ssake_extend(
    contig: str,
    r1_fq: Path,
    r2_fq: Path,
    workdir: Path,
    min_overlap: int = 20,
) -> str:
    """Extend contig using SSAKE's overlap-layout-consensus.

    SSAKE uses a prefix tree for rapid read overlap detection,
    complementing k-mer extension (greedy, exact-overlap) and
    minimap2 extension (chain-based, soft-clip). SSAKE can resolve
    moderately repetitive regions where k-mer extension stalls due
    to its all-at-once prefix search over the full read set.

    Runs via 'micromamba run -n ssake SSAKE' (separate Python 2.7 env).
    Returns extended contig (unchanged if SSAKE fails or no growth).
    """
    if len(contig) < 50:
        return contig

    ssake_dir = workdir / "_ssake"
    ssake_dir.mkdir(parents=True, exist_ok=True)

    # SSAKE input: seed FASTA (-s) + reads FASTA (-f)
    seed_fa = ssake_dir / "seed.fa"
    with open(seed_fa, "w") as fh:
        fh.write(f">seed\n{contig}\n")

    # Convert paired FASTQ to interleaved FASTA for SSAKE
    reads_fa = ssake_dir / "reads.fa"
    n_reads = 0
    with open(reads_fa, "w") as out:
        for fq_path in (r1_fq, r2_fq):
            if not fq_path.exists():
                continue
            with open(fq_path) as fq:
                while True:
                    header = fq.readline()
                    if not header:
                        break
                    seq = fq.readline().strip()
                    fq.readline()  # +
                    fq.readline()  # qual
                    if seq and len(seq) >= min_overlap:
                        out.write(f">r{n_reads}\n{seq}\n")
                        n_reads += 1

    if n_reads < 10:
        shutil.rmtree(ssake_dir, ignore_errors=True)
        return contig

    # Run SSAKE: -s seed, -f reads, -m min_overlap, -o 1 (1 pass), -w 1
    ssake_prefix = ssake_dir / "ssake_out"
    try:
        subprocess.run(
            ["micromamba", "run", "-n", "ssake", "SSAKE",
             "-f", str(reads_fa),
             "-s", str(seed_fa),
             "-b", str(ssake_prefix),
             "-m", str(min_overlap),
             "-o", "1",        # single pass
             "-w", "1",        # min reads for base call
             ],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            timeout=120,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError):
        shutil.rmtree(ssake_dir, ignore_errors=True)
        return contig

    # Read SSAKE contigs — pick longest that contains original seed
    ssake_contigs = Path(f"{ssake_prefix}.contigs")
    best = contig
    if ssake_contigs.exists():
        seqs = read_fasta(ssake_contigs)
        for name, seq in seqs.items():
            if len(seq) > len(best):
                # Check the SSAKE contig extends the seed (contains it or overlaps)
                if contig[:50] in seq or contig[-50:] in seq:
                    best = seq

    shutil.rmtree(ssake_dir, ignore_errors=True)
    return best


def _check_merge(contig_5p: str, contig_3p: str, min_overlap: int = 30) -> str | None:
    """Check if contig_5p and contig_3p overlap. Return merged seq or None."""
    if not contig_5p or not contig_3p:
        return None

    # Check if 3' end of 5p overlaps with 5' end of 3p
    max_check = min(len(contig_5p), len(contig_3p), 500)
    for ovl in range(max_check, min_overlap - 1, -1):
        if contig_5p[-ovl:] == contig_3p[:ovl]:
            merged = contig_5p + contig_3p[ovl:]
            log(f"    MERGE: {ovl}bp overlap detected! "
                f"→ {len(merged):,}bp merged contig")
            return merged
    return None


# ---------------------------------------------------------------------------
# Write pool FASTQ for Pilon
# ---------------------------------------------------------------------------

def _write_pool_fastq(extender: StrandAwareSeedExtender,
                      r1_out: Path, r2_out: Path) -> None:
    """Write extender's read pool as pseudo-paired FASTQ for Pilon mapping."""
    with open(r1_out, "w") as f1, open(r2_out, "w") as f2:
        for i in range(0, len(extender.seqs) - 1, 2):
            r1 = extender.seqs[i]
            r2_rc = extender.seqs[i + 1]
            r2_raw = revcomp(r2_rc)
            name = f"read_{i // 2}"
            f1.write(f"@{name}/1\n{r1}\n+\n{'I' * len(r1)}\n")
            f2.write(f"@{name}/2\n{r2_raw}\n+\n{'I' * len(r2_raw)}\n")


# ---------------------------------------------------------------------------
# Main assembly loop
# ---------------------------------------------------------------------------

def assemble_insert(
    site: InsertionSite,
    candidate_r1: Path,
    candidate_r2: Path,
    host_bam: Path,
    host_ref: Path,
    element_db: Path,
    workdir: Path,
    threads: int = 4,
    max_rounds: int = 5,  # T3: reduced from 8 per T2 measurement (docs/measurements/max_rounds_study.md)
    ext_k: int = 15,
    recruit_k: int = 25,
    gap_size: int = 1000,
    s03_r1: Path | None = None,
    s03_r2: Path | None = None,
) -> tuple[Path | None, int, str]:
    """Main loop alternating k-mer extension and Pilon gap fill.

    Returns (fasta_path, n_rounds, status).
    Status: 'complete', 'partial', 'no_seeds', 'no_assembly'.
    """
    # Lazy import to avoid circular import: extract_unmapped_paired is still
    # defined in the monolith (Phase 2 territory, Session 5). A module-top
    # import would fail during initialization because the monolith imports
    # assembly.py before extract_unmapped_paired is defined in its own body.
    # This local import is safe: by call time, the monolith is fully loaded.
    # Pattern mirrors conftest.py (repo root → sys.path) and per_site.py shim.
    try:
        from scripts.s05_insert_assembly import extract_unmapped_paired
    except ImportError:
        import sys as _sys
        from pathlib import Path as _Path
        _repo = _Path(__file__).resolve().parents[2]
        if str(_repo) not in _sys.path:
            _sys.path.insert(0, str(_repo))
        from scripts.s05_insert_assembly import extract_unmapped_paired

    workdir.mkdir(parents=True, exist_ok=True)

    seed_5p = site.seed_5p
    seed_3p = site.seed_3p

    if not seed_5p and not seed_3p:
        log(f"  No seeds for {site.site_id}")
        return None, 0, "no_seeds"

    log(f"  Seeds: 5p={len(seed_5p)}bp, 3p={len(seed_3p)}bp")

    # Initialize extender
    extender = StrandAwareSeedExtender(
        k=ext_k, min_overlap=20, min_depth=2, min_ratio=0.7,
    )

    # Load construct-extracted reads (s03) as PRIMARY pool — these are the
    # reads that actually hit the construct, much more targeted than host BAM
    if s03_r1 and s03_r2 and s03_r1.exists() and s03_r2.exists():
        log(f"  Loading construct-extracted reads (s03)...")
        n_s03 = extender.load_paired_reads(s03_r1, s03_r2)
        log(f"  S03 reads: {n_s03:,} pairs")
    else:
        n_s03 = 0

    # Load junction-region candidate reads (supplementary)
    log(f"  Loading junction-region reads...")
    n_pairs = extender.load_paired_reads(candidate_r1, candidate_r2)
    log(f"  Total pool: {len(extender.seqs):,} seqs, "
        f"{len(extender.kmer_index):,} unique {ext_k}-mers")

    # Extract unmapped reads (cached)
    unmapped_r1, unmapped_r2 = extract_unmapped_paired(
        host_bam, workdir, threads=threads,
    )

    # Initial k-mer extension of both seeds
    contig_5p = seed_5p
    contig_3p = seed_3p

    if contig_5p:
        log(f"  Initial extension of 5p seed ({len(contig_5p)}bp)...")
        contig_5p = extender.extend(contig_5p, max_iterations=100)
        log(f"  5p: {len(seed_5p)} → {len(contig_5p):,}bp")

    if contig_3p:
        log(f"  Initial extension of 3p seed ({len(contig_3p)}bp)...")
        contig_3p = extender.extend(contig_3p, max_iterations=100)
        log(f"  3p: {len(seed_3p)} → {len(contig_3p):,}bp")

    # Check immediate merge
    merged = _check_merge(contig_5p, contig_3p)
    if merged:
        final_fa = workdir / f"{site.site_id}_insert.fasta"
        write_fasta(final_fa, f"{site.site_id}_assembled_insert", merged)
        log(f"  COMPLETE after initial extension: {len(merged):,}bp")
        return final_fa, 0, "complete"

    # Main iterative loop
    total_rounds = 0
    growth_history: list[int] = []
    status = "partial"

    for rnd in range(1, max_rounds + 1):
        prev_5p_len = len(contig_5p)
        prev_3p_len = len(contig_3p)
        log(f"\n  --- Round {rnd} ---")

        # Step A: K-mer recruit from unmapped pool
        recruit_contig = contig_5p + contig_3p
        r1_new, r2rc_new = recruit_by_kmer(
            recruit_contig, unmapped_r1, unmapped_r2, k=recruit_k,
        )
        if r1_new:
            n_added = extender.add_seqs(r1_new + r2rc_new)
            log(f"    Recruit: {len(r1_new)} pairs → {n_added} new seqs")
        else:
            log(f"    Recruit: 0 pairs")

        # Step B: K-mer extension (strand-aware)
        if contig_5p:
            contig_5p = extender.extend(contig_5p, max_iterations=100)
        if contig_3p:
            contig_3p = extender.extend(contig_3p, max_iterations=100)

        kmer_growth = (len(contig_5p) - prev_5p_len) + (len(contig_3p) - prev_3p_len)
        log(f"    K-mer ext: 5p={prev_5p_len:,}→{len(contig_5p):,}, "
            f"3p={prev_3p_len:,}→{len(contig_3p):,} (+{kmer_growth})")

        # Step B2: minimap2 soft-clip extension (handles repeats k-mer can't)
        mm2_dir = workdir / f"_mm2_r{rnd}"
        mm2_dir.mkdir(parents=True, exist_ok=True)
        pool_r1_mm2 = mm2_dir / "pool_R1.fq"
        pool_r2_mm2 = mm2_dir / "pool_R2.fq"
        _write_pool_fastq(extender, pool_r1_mm2, pool_r2_mm2)

        pre_mm2_5p = len(contig_5p)
        pre_mm2_3p = len(contig_3p)
        if contig_5p:
            contig_5p = _minimap2_extend(
                contig_5p, pool_r1_mm2, pool_r2_mm2, mm2_dir,
                threads=threads,
            )
        if contig_3p:
            contig_3p = _minimap2_extend(
                contig_3p, pool_r1_mm2, pool_r2_mm2, mm2_dir,
                threads=threads,
            )
        mm2_growth = (len(contig_5p) - pre_mm2_5p) + (len(contig_3p) - pre_mm2_3p)
        if mm2_growth > 0:
            log(f"    mm2 ext:  5p={pre_mm2_5p:,}→{len(contig_5p):,}, "
                f"3p={pre_mm2_3p:,}→{len(contig_3p):,} (+{mm2_growth})")
            # Recruit new reads matching the mm2-extended region
            r1_mm2, r2rc_mm2 = recruit_by_kmer(
                contig_5p + contig_3p, unmapped_r1, unmapped_r2, k=recruit_k,
            )
            if r1_mm2:
                n_mm2 = extender.add_seqs(r1_mm2 + r2rc_mm2)
                if n_mm2:
                    log(f"    mm2 recruit: {n_mm2} new seqs")
        shutil.rmtree(mm2_dir, ignore_errors=True)

        # Step C: Check merge
        merged = _check_merge(contig_5p, contig_3p)
        if merged:
            final_fa = workdir / f"{site.site_id}_insert.fasta"
            write_fasta(final_fa, f"{site.site_id}_assembled_insert", merged)
            log(f"  COMPLETE: contigs merged → {len(merged):,}bp (round {rnd})")
            return final_fa, rnd, "complete"

        # Step D: Pilon gap fill
        pre_pilon_5p = len(contig_5p)
        pre_pilon_3p = len(contig_3p)
        pilon_growth = 0
        if contig_5p and contig_3p:
            pilon_dir = workdir / f"_pilon_r{rnd}"
            pilon_dir.mkdir(parents=True, exist_ok=True)

            pool_r1 = pilon_dir / "pool_R1.fq"
            pool_r2 = pilon_dir / "pool_R2.fq"
            _write_pool_fastq(extender, pool_r1, pool_r2)

            contig_5p, contig_3p, gap_filled = pilon_fill(
                contig_5p, contig_3p,
                pool_r1, pool_r2,
                pilon_dir,
                gap_size=gap_size,
                threads=threads,
            )

            shutil.rmtree(pilon_dir, ignore_errors=True)
            pilon_growth = (len(contig_5p) - pre_pilon_5p) + (len(contig_3p) - pre_pilon_3p)

            if gap_filled:
                final_seq = contig_5p
                final_fa = workdir / f"{site.site_id}_insert.fasta"
                write_fasta(final_fa, f"{site.site_id}_assembled_insert", final_seq)
                log(f"  COMPLETE: Pilon filled gap → {len(final_seq):,}bp (round {rnd})")
                return final_fa, rnd, "complete"

        # Step D2: SSAKE overlap-layout-consensus extension
        ssake_dir = workdir / f"_ssake_r{rnd}"
        ssake_dir.mkdir(parents=True, exist_ok=True)
        ssake_r1 = ssake_dir / "pool_R1.fq"
        ssake_r2 = ssake_dir / "pool_R2.fq"
        _write_pool_fastq(extender, ssake_r1, ssake_r2)

        pre_ssake_5p = len(contig_5p)
        pre_ssake_3p = len(contig_3p)
        if contig_5p:
            contig_5p = _ssake_extend(contig_5p, ssake_r1, ssake_r2, ssake_dir)
        if contig_3p:
            contig_3p = _ssake_extend(contig_3p, ssake_r1, ssake_r2, ssake_dir)
        ssake_growth = (len(contig_5p) - pre_ssake_5p) + (len(contig_3p) - pre_ssake_3p)
        if ssake_growth > 0:
            log(f"    SSAKE:   5p={pre_ssake_5p:,}→{len(contig_5p):,}, "
                f"3p={pre_ssake_3p:,}→{len(contig_3p):,} (+{ssake_growth})")
            # Check merge after SSAKE
            merged = _check_merge(contig_5p, contig_3p)
            if merged:
                final_fa = workdir / f"{site.site_id}_insert.fasta"
                write_fasta(final_fa, f"{site.site_id}_assembled_insert", merged)
                log(f"  COMPLETE: SSAKE+merge → {len(merged):,}bp (round {rnd})")
                shutil.rmtree(ssake_dir, ignore_errors=True)
                return final_fa, rnd, "complete"
        shutil.rmtree(ssake_dir, ignore_errors=True)

        # Step E: Check termination — only converge when ALL 4 steps show zero growth
        total_rounds = rnd
        log(f"    Round {rnd} growth: kmer={kmer_growth}, mm2={mm2_growth}, "
            f"pilon={pilon_growth}, ssake={ssake_growth}")

        if kmer_growth == 0 and mm2_growth == 0 and pilon_growth == 0 and ssake_growth == 0:
            log(f"    All 4 assemblers show zero growth → converged")
            status = "converged"
            break

        # Check if contig ends reached host genome
        reached_5p, reached_3p = check_host_termination(
            contig_5p, contig_3p, host_ref, workdir,
        )
        if reached_5p:
            log(f"    5' contig reached host genome!")
        if reached_3p:
            log(f"    3' contig reached host genome!")
        if reached_5p and reached_3p:
            log(f"    >>> Both ends at host genome — insert fully traversed!")
            status = "complete"
            break

        # Cycle detection — per-step growth pattern over 3 rounds
        total_growth = kmer_growth + mm2_growth + pilon_growth + ssake_growth
        growth_history.append(total_growth)
        if len(growth_history) >= 3:
            last3 = growth_history[-3:]
            if last3[0] == last3[1] == last3[2] and last3[0] > 0:
                log(f"    Same total growth ({last3[0]}) for 3 rounds → cycling")
                status = "converged"
                break
        if len(growth_history) >= 6:
            recent = growth_history[-3:]
            earlier = growth_history[-6:-3]
            if recent == earlier:
                log(f"    Cyclic growth pattern detected → stopping")
                status = "converged"
                break

    # Write final result
    if contig_5p and contig_3p:
        # Output both contigs separately if not merged
        final_seq = contig_5p + ("N" * 100) + contig_3p
        final_fa = workdir / f"{site.site_id}_insert.fasta"
        write_fasta(final_fa, f"{site.site_id}_assembled_insert", final_seq)
        log(f"  Final: 5p={len(contig_5p):,}bp + 3p={len(contig_3p):,}bp "
            f"= {len(final_seq):,}bp after {total_rounds} rounds → {status}")
    elif contig_5p:
        final_fa = workdir / f"{site.site_id}_insert.fasta"
        write_fasta(final_fa, f"{site.site_id}_assembled_insert", contig_5p)
        log(f"  Final: {len(contig_5p):,}bp (5p only) after {total_rounds} rounds")
    elif contig_3p:
        final_fa = workdir / f"{site.site_id}_insert.fasta"
        write_fasta(final_fa, f"{site.site_id}_assembled_insert", contig_3p)
        log(f"  Final: {len(contig_3p):,}bp (3p only) after {total_rounds} rounds")
    else:
        return None, total_rounds, "no_assembly"

    # Phase 3b: Foreign read refinement (minimap2 chaining + Pilon)
    # After k-mer assembly converges, use foreign reads (s03 unmapped to host)
    # with iterative minimap2+Pilon to extend construct regions.
    # This resolves palindromic/inverted repeat regions that k-mer extension
    # cannot handle (e.g., head-to-head T-DNA structures).
    if s03_r1 and s03_r2 and s03_r1.exists() and s03_r2.exists():
        log(f"\n  === Foreign read refinement (minimap2 + Pilon) ===")
        final_fa = refine_with_foreign_reads(
            insert_fasta=final_fa,
            s03_r1=s03_r1,
            s03_r2=s03_r2,
            host_ref=host_ref,
            workdir=workdir,
            threads=threads,
            max_rounds=10,
        )

    return final_fa, total_rounds, status
