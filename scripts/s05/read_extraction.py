"""Phase 2 — Read Extraction (Issue #4, Session 5).

Extracted verbatim from ``scripts/s05_insert_assembly.py``. Owns the
candidate read extraction pipeline that pulls clip-anchored reads + their
mates from the host BAM for downstream Phase 3 assembly. Stage 2 in the
s05 DAG; imports from ``.primitives`` (stage 0) — no imports from
``.site_discovery`` (stage 1) are needed because ``_extract_seeds_at_positions``
is defined in site_discovery (it is solely called by
``legacy_junctions_to_sites``, a Phase 1 function).

Functions
---------
- ``extract_candidate_reads`` — public entry; extracts reads + mates near
  insertion sites, writes paired FASTQs for assembly
- ``extract_unmapped_paired`` — pulls unmapped read pairs (recovery for
  sites where one mate failed to map)

Replaces the per_site.py shim re-exports pending Session 7 shim retirement.
"""
from __future__ import annotations

import gzip
import os
import subprocess
from pathlib import Path

import pysam

from .primitives import (
    log,
    InsertionSite,
)


# ---------------------------------------------------------------------------
# Phase 2: Candidate Read Extraction
# ---------------------------------------------------------------------------

def extract_candidate_reads(
    host_bam: Path,
    site: InsertionSite,
    out_r1: Path,
    out_r2: Path,
    flank: int = 5000,
    threads: int = 4,
    min_clip: int = 20,
    junction_window: int = 20,
) -> int:
    """Extract reads for assembly from junction regions, phased by allele.

    For hemizygous insertions, mixing WT-allele and insertion-allele reads
    causes assembly ambiguity. We phase reads using soft-clip diagnostics:

    1. Scan junction position(s) ±junction_window
    2. INSERTION-allele: reads with soft-clip at junction position
       (these span the host→T-DNA boundary)
    3. WT-allele: reads that span the junction position WITHOUT a clip
       at that boundary (they read straight through → no insertion)
    4. AMBIGUOUS: distal reads (>junction_window from junction).
       These are kept ONLY IF their mate is not classified as WT.
       A distal read paired with a WT-classified mate is excluded
       (the WT evidence from the mate overrides the distal's ambiguity).
    5. Mate propagation: if any read in a pair is insertion-allele,
       the whole pair is insertion-allele (insertion wins ties).
    6. Output: insertion-allele + ambiguous-without-WT-mate pairs,
       excluding pairs where any mate has WT evidence and no mate
       has insertion evidence (exclude = wt_reads - insertion_reads).
    """
    tmp_dir = out_r1.parent

    # Build regions
    regions = []
    positions = []
    if site.pos_5p > 0:
        positions.append(site.pos_5p)
    if site.pos_3p > 0:
        positions.append(site.pos_3p)
    if not positions:
        return 0

    for pos in positions:
        start = max(1, pos - flank)
        end = pos + flank
        regions.append(f"{site.host_chr}:{start}-{end}")

    # Extract reads from junction regions
    region_bam = tmp_dir / f"_{site.site_id}_regions.bam"
    with open(region_bam, "wb") as bam_fh:
        subprocess.run(
            ["samtools", "view", "-b", "-@", str(threads),
             str(host_bam)] + regions,
            stdout=bam_fh,
            stderr=subprocess.DEVNULL, check=True,
        )

    # Phase reads by allele using pysam
    insertion_reads: set[str] = set()
    wt_reads: set[str] = set()
    distal_reads: set[str] = set()

    with pysam.AlignmentFile(str(region_bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            if qname is None:
                continue
            cigar = read.cigartuples
            if not cigar:
                distal_reads.add(qname)
                continue

            ref_start = read.reference_start  # 0-based
            ref_end = read.reference_end       # exclusive
            if ref_start is None or ref_end is None:
                distal_reads.add(qname)
                continue

            classified = False
            for pos in positions:
                pos0 = pos - 1  # convert to 0-based
                # Distance check: is read NEAR this junction at all?
                if ref_end < pos0 - junction_window or ref_start > pos0 + junction_window:
                    continue

                # INSERTION diagnostic: soft-clip end touches junction
                left_clip = cigar[0][1] if cigar[0][0] == 4 else 0
                right_clip = cigar[-1][1] if cigar[-1][0] == 4 else 0

                # Right soft-clip → break point is at ref_end → matches 5p junction
                is_insertion = False
                if right_clip >= min_clip and abs(ref_end - pos0) <= junction_window:
                    is_insertion = True
                # Left soft-clip → break point is at ref_start → matches 3p junction
                if left_clip >= min_clip and abs(ref_start - pos0) <= junction_window:
                    is_insertion = True

                if is_insertion:
                    insertion_reads.add(qname)
                    classified = True
                    break

                # WT diagnostic: spans junction (start before & end after) WITHOUT clip
                # at that boundary
                if ref_start <= pos0 - junction_window and ref_end >= pos0 + junction_window:
                    # No clip at the junction boundary
                    if right_clip < min_clip and left_clip < min_clip:
                        wt_reads.add(qname)
                        classified = True
                        break

            if not classified:
                distal_reads.add(qname)

    # Mate propagation: insertion wins over WT.
    # samtools view -N fetches BOTH mates by name, so we only need a name set.
    # A pair is excluded if ANY mate has WT evidence (spans junction w/o clip)
    # AND no mate is insertion-classified. Distal does NOT rescue WT: the
    # 3' mate of a pair is ~300bp away from the junction and gets classified
    # as distal, but that doesn't contradict the 5' mate's WT observation.
    all_names = insertion_reads | wt_reads | distal_reads
    exclude_names = wt_reads - insertion_reads
    keep_names = all_names - exclude_names

    log(f"    Phasing: {len(insertion_reads):,} insertion-allele, "
        f"{len(wt_reads):,} WT-allele, {len(distal_reads):,} distal/ambiguous")
    log(f"    Excluded {len(exclude_names):,} pure-WT pairs")

    namelist = tmp_dir / f"_{site.site_id}_names.txt"
    with open(namelist, "w") as fh:
        for n in sorted(keep_names):
            fh.write(n + "\n")
    log(f"    Junction region reads (phased): {len(keep_names):,} read names")

    # Extract both mates from full BAM
    both_mates_bam = tmp_dir / f"_{site.site_id}_mates.bam"
    with open(both_mates_bam, "wb") as bam_fh:
        subprocess.run(
            ["samtools", "view", "-b", "-N", str(namelist),
             "-@", str(threads), str(host_bam)],
            stdout=bam_fh,
            stderr=subprocess.DEVNULL, check=True,
        )

    # Also get all-unmapped pairs (both mates unmapped, flag 12)
    unmapped_bam = tmp_dir / f"_{site.site_id}_unmapped.bam"
    with open(unmapped_bam, "wb") as bam_fh:
        subprocess.run(
            ["samtools", "view", "-b", "-f", "12", "-@", str(threads), str(host_bam)],
            stdout=bam_fh,
            stderr=subprocess.DEVNULL, check=True,
        )

    # Merge
    merged_bam = tmp_dir / f"_{site.site_id}_merged.bam"
    subprocess.run(
        ["samtools", "merge", "-f", "-@", str(threads), str(merged_bam),
         str(both_mates_bam), str(unmapped_bam)],
        stderr=subprocess.DEVNULL, check=True,
    )

    # Name-sort and convert to paired FASTQ
    nsort_bam = tmp_dir / f"_{site.site_id}_nsort.bam"
    subprocess.run(
        ["samtools", "sort", "-n", "-@", str(threads),
         str(merged_bam), "-o", str(nsort_bam)],
        stderr=subprocess.DEVNULL, check=True,
    )
    subprocess.run(
        ["samtools", "fastq", "-@", str(threads),
         "-1", str(out_r1), "-2", str(out_r2),
         "-s", "/dev/null", "-0", "/dev/null", str(nsort_bam)],
        stderr=subprocess.DEVNULL, check=True,
    )

    # Count reads
    count = 0
    opener = gzip.open if str(out_r1).endswith(".gz") else open
    with opener(str(out_r1), "rt") as fh:
        for _ in fh:
            count += 1
    count = count // 4  # FASTQ: 4 lines per read

    # Cleanup temp files
    for p in [region_bam, namelist, both_mates_bam, unmapped_bam,
              merged_bam, nsort_bam]:
        p.unlink(missing_ok=True)

    return count


def extract_unmapped_paired(
    host_bam: Path,
    workdir: Path,
    threads: int = 4,
) -> tuple[Path, Path]:
    """Extract unmapped reads from host BAM as paired FASTQ (cached)."""
    r1 = workdir / "_unmapped_R1.fq.gz"
    r2 = workdir / "_unmapped_R2.fq.gz"
    if r1.exists() and r2.exists():
        n = 0
        with gzip.open(r1, "rt") as fh:
            for _ in fh:
                n += 1
        log(f"  Unmapped reads (cached): {n // 4:,} pairs")
        return r1, r2

    log(f"  Extracting unmapped reads from host BAM...")
    # The cache is site-independent and shared on `workdir` across the fan-out
    # (issue #18). Two concurrent sites that both miss the cache would otherwise
    # write `r1`/`r2` at the same time and interleave into a corrupt gzip. Write
    # to PID-unique temp files and atomically rename into place; a late writer
    # simply replaces an identical file.
    tmp_r1 = r1.with_suffix(r1.suffix + f".tmp.{os.getpid()}")
    tmp_r2 = r2.with_suffix(r2.suffix + f".tmp.{os.getpid()}")
    subprocess.run(
        f"samtools view -f 4 -@ {threads} -b {host_bam} "
        f"| samtools sort -n -@ {threads} - "
        f"| samtools fastq -1 {tmp_r1} -2 {tmp_r2} -s /dev/null -0 /dev/null - "
        f"2>/dev/null",
        shell=True, check=True,
    )
    os.replace(tmp_r1, r1)
    os.replace(tmp_r2, r2)
    n = 0
    with gzip.open(r1, "rt") as fh:
        for _ in fh:
            n += 1
    log(f"  Unmapped reads: {n // 4:,} pairs")
    return r1, r2
