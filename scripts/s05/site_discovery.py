"""Phase 1 — Site Discovery (Issue #4, Session 5).

Extracted verbatim from ``scripts/s05_insert_assembly.py``. Owns the
soft-clip junction detection pipeline that turns a host-mapped BAM
into a list of ``InsertionSite`` candidates. Stage 1 in the s05 DAG;
imports only from ``.primitives`` (stage 0) and ``.classify`` (stage 1,
same-stage import allowed).

Pipeline overview
-----------------
``find_softclip_junctions`` is the public entry point. It scans the
host BAM for soft-clipped reads, clusters their clips into per-position
``JunctionCluster`` objects, and builds ``InsertionSite`` records for
downstream Phase 1.5 classify and Phase 2/3 read-extraction +
assembly.

Helpers
-------
- ``_build_consensus`` — majority-vote consensus over a clip stack
- ``_batch_check_maps_to_host`` — BLAST clip seqs to host, mark host-mapping
- ``_extract_seeds_at_positions`` — fills InsertionSite seed_5p / seed_3p
  from BAM at known junction positions (called by legacy_junctions_to_sites)
- ``_apply_mask_bed`` — pre-filter sites by user-supplied BED mask (T10)
- ``parse_legacy_junctions`` / ``legacy_junctions_to_sites`` — fallback
  path that ingests step-6 junctions.tsv when running step 5 alone

Architectural note
------------------
``_extract_seeds_at_positions`` is defined here (not in read_extraction)
because it is solely called by ``legacy_junctions_to_sites`` (Phase 1).
It calls ``_build_consensus`` defined in this same module, so no import
from a later stage is needed. Moving it to read_extraction.py (stage 2)
would require site_discovery (stage 1) to import from a later stage —
a DAG violation. This deviation from the session-5 plan is documented
in the commit message.

Replaces the v1.0 shim that re-exported these names from the monolith.
"""
from __future__ import annotations

import csv
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import pysam

from .primitives import (
    log,
    InsertionSite,
    JunctionCluster,
    LegacyJunction,
)


# ---------------------------------------------------------------------------
# Phase 1: Soft-clip Junction Detection
# ---------------------------------------------------------------------------

def _build_consensus(seqs: list[str], direction: str) -> str:
    """Build majority-vote consensus from a list of clipped sequences.

    For 'right' clips: sequences are aligned from their left (start of insert).
    For 'left' clips: sequences are aligned from their right (end of insert).
    """
    if not seqs:
        return ""

    if direction == "right":
        # Align from left
        votes: dict[int, Counter] = defaultdict(Counter)
        for seq in seqs:
            for i, base in enumerate(seq):
                votes[i][base.upper()] += 1
        consensus = []
        for p in range(max(votes.keys()) + 1):
            if p not in votes:
                break
            v = votes[p]
            total = sum(v.values())
            best, cnt = v.most_common(1)[0]
            if total >= 2 and cnt / total >= 0.51:
                consensus.append(best)
            else:
                break
        return "".join(consensus)
    else:
        # Align from right
        votes = defaultdict(Counter)
        max_len = max(len(s) for s in seqs)
        for seq in seqs:
            offset = max_len - len(seq)
            for i, base in enumerate(seq):
                votes[offset + i][base.upper()] += 1
        consensus = []
        for p in range(max_len):
            if p not in votes:
                consensus.append("N")
                continue
            v = votes[p]
            total = sum(v.values())
            best, cnt = v.most_common(1)[0]
            if total >= 2 and cnt / total >= 0.51:
                consensus.append(best)
            else:
                consensus.append("N")
        # Trim leading/trailing N
        result = "".join(consensus).strip("N")
        return result


def _batch_check_maps_to_host(seqs: dict[str, str], host_ref: Path, workdir: Path,
                              min_identity: float = 0.90,
                              min_coverage: float = 0.80) -> set[str]:
    """Batch check which sequences map to host genome using one minimap2 call.

    Args:
        seqs: dict of {name: sequence} to check
        Returns: set of names that map to host
    """
    if not seqs:
        return set()

    query_fa = workdir / "_clip_check_batch.fa"
    with open(query_fa, "w") as fh:
        for name, seq in seqs.items():
            if len(seq) >= 20:
                fh.write(f">{name}\n{seq}\n")

    result = subprocess.run(
        ["minimap2", "-c", "--secondary=no", "-t", "4",
         str(host_ref), str(query_fa)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    query_fa.unlink(missing_ok=True)

    maps_to_host: set[str] = set()
    for line in result.stdout.splitlines():
        cols = line.strip().split("\t")
        if len(cols) < 12:
            continue
        qname = cols[0]
        q_len = int(cols[1])
        match_bp = int(cols[9])
        block_len = int(cols[10])
        if block_len == 0:
            continue
        identity = match_bp / block_len
        coverage = (int(cols[3]) - int(cols[2])) / q_len if q_len > 0 else 0
        if identity >= min_identity and coverage >= min_coverage:
            maps_to_host.add(qname)
    return maps_to_host


# ---------------------------------------------------------------------------
# _batch_check_element_hits moved to scripts/s05/classify.py
# (Issue #4 Session 4) and re-imported at the top of this module.
# ---------------------------------------------------------------------------
def find_softclip_junctions(
    host_bam: Path,
    host_ref: Path,
    element_db: Path | None,
    workdir: Path,
    min_clip: int = 20,
    cluster_window: int = 50,
    extra_dbs: list[Path] | None = None,
) -> list[InsertionSite]:
    """Scan host BAM for insertion sites using bidirectional soft-clip analysis.

    Three conditions for a validated insertion site:
    1. Bidirectional soft-clips at the same genomic position
    2. The two clip sequences are DIFFERENT (insert, not SV)
    3. Neither clip maps to host genome (foreign sequence)
    """
    from .classify import _batch_check_element_hits

    workdir.mkdir(parents=True, exist_ok=True)
    log("Phase 1: Scanning host BAM for soft-clip junctions...")

    # Step 1: Collect all soft-clipped reads
    right_clips: dict[str, list[tuple[int, str]]] = defaultdict(list)  # chr -> [(pos, seq)]
    left_clips: dict[str, list[tuple[int, str]]] = defaultdict(list)

    bam = pysam.AlignmentFile(str(host_bam), "rb")
    n_right = 0
    n_left = 0

    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        cigar = read.cigartuples
        if cigar is None:
            continue
        seq = read.query_sequence
        if seq is None:
            continue

        chrom = read.reference_name

        # Right clip (3' end clipped — forward read approaching junction)
        if cigar[-1][0] == 4 and cigar[-1][1] >= min_clip:
            clip_len = cigar[-1][1]
            clip_pos = read.reference_end  # position where clip starts
            clip_seq = seq[-clip_len:]
            right_clips[chrom].append((clip_pos, clip_seq))
            n_right += 1

        # Left clip (5' end clipped — reverse read approaching junction)
        if cigar[0][0] == 4 and cigar[0][1] >= min_clip:
            clip_len = cigar[0][1]
            clip_pos = read.reference_start  # position where clip ends
            clip_seq = seq[:clip_len]
            left_clips[chrom].append((clip_pos, clip_seq))
            n_left += 1

    bam.close()
    log(f"  Collected {n_right:,} right-clips, {n_left:,} left-clips")

    # Step 2: Cluster by position (min_depth=3 to filter noise)
    MIN_CLUSTER_DEPTH = 3

    def _cluster(clips: list[tuple[int, str]], window: int) -> list[tuple[int, list[str]]]:
        if not clips:
            return []
        clips.sort(key=lambda x: x[0])
        clusters = []
        current_seqs = [clips[0][1]]
        current_positions = [clips[0][0]]

        for pos, seq in clips[1:]:
            if pos - current_positions[0] <= window:
                current_seqs.append(seq)
                current_positions.append(pos)
            else:
                if len(current_seqs) >= MIN_CLUSTER_DEPTH:
                    median_pos = sorted(current_positions)[len(current_positions) // 2]
                    clusters.append((median_pos, current_seqs))
                current_seqs = [seq]
                current_positions = [pos]
        if len(current_seqs) >= MIN_CLUSTER_DEPTH:
            median_pos = sorted(current_positions)[len(current_positions) // 2]
            clusters.append((median_pos, current_seqs))

        return clusters

    # Step 3: Pair forward/reverse clusters at same position
    sites: list[InsertionSite] = []
    used_ids: set[str] = set()

    for chrom in set(list(right_clips.keys()) + list(left_clips.keys())):
        r_clusters = _cluster(right_clips.get(chrom, []), cluster_window)
        l_clusters = _cluster(left_clips.get(chrom, []), cluster_window)

        # Try to pair right and left clusters within window
        used_l = set()
        for r_pos, r_seqs in r_clusters:
            best_l = None
            best_dist = cluster_window + 1
            for li, (l_pos, l_seqs) in enumerate(l_clusters):
                if li in used_l:
                    continue
                dist = abs(r_pos - l_pos)
                if dist <= cluster_window and dist < best_dist:
                    best_l = li
                    best_dist = dist

            if best_l is not None:
                l_pos, l_seqs = l_clusters[best_l]
                used_l.add(best_l)

                r_consensus = _build_consensus(r_seqs, "right")
                l_consensus = _build_consensus(l_seqs, "left")

                if len(r_consensus) < min_clip or len(l_consensus) < min_clip:
                    continue

                anchor_pos = min(r_pos, l_pos)
                sid = f"insertion_{chrom}_{anchor_pos}"
                if sid in used_ids:
                    suffix = 2
                    while f"{sid}_{suffix}" in used_ids:
                        suffix += 1
                    sid = f"{sid}_{suffix}"
                used_ids.add(sid)
                jc_5p = JunctionCluster(
                    host_chr=chrom, position=r_pos,
                    clip_direction="right",
                    clipped_seqs=r_seqs, consensus_clip=r_consensus,
                    n_reads=len(r_seqs),
                )
                jc_3p = JunctionCluster(
                    host_chr=chrom, position=l_pos,
                    clip_direction="left",
                    clipped_seqs=l_seqs, consensus_clip=l_consensus,
                    n_reads=len(l_seqs),
                )

                site = InsertionSite(
                    site_id=sid,
                    host_chr=chrom,
                    junction_5p=jc_5p,
                    junction_3p=jc_3p,
                    seed_5p=r_consensus,
                    seed_3p=l_consensus,
                    pos_5p=min(r_pos, l_pos),
                    pos_3p=max(r_pos, l_pos),
                )
                sites.append(site)
            else:
                # Single-direction cluster: only keep high-depth (≥5 reads)
                if len(r_seqs) >= 5:
                    r_consensus = _build_consensus(r_seqs, "right")
                    if len(r_consensus) >= min_clip:
                        sid = f"insertion_{chrom}_{r_pos}"
                        if sid in used_ids:
                            suffix = 2
                            while f"{sid}_{suffix}" in used_ids:
                                suffix += 1
                            sid = f"{sid}_{suffix}"
                        used_ids.add(sid)
                        jc_5p = JunctionCluster(
                            host_chr=chrom, position=r_pos,
                            clip_direction="right",
                            clipped_seqs=r_seqs, consensus_clip=r_consensus,
                            n_reads=len(r_seqs),
                        )
                        site = InsertionSite(
                            site_id=sid,
                            host_chr=chrom,
                            junction_5p=jc_5p,
                            seed_5p=r_consensus,
                            pos_5p=r_pos,
                            confidence="low",
                        )
                        sites.append(site)

        # Unpaired left clusters (only high-depth)
        for li, (l_pos, l_seqs) in enumerate(l_clusters):
            if li in used_l:
                continue
            if len(l_seqs) < 5:
                continue
            l_consensus = _build_consensus(l_seqs, "left")
            if len(l_consensus) >= min_clip:
                sid = f"insertion_{chrom}_{l_pos}"
                if sid in used_ids:
                    suffix = 2
                    while f"{sid}_{suffix}" in used_ids:
                        suffix += 1
                    sid = f"{sid}_{suffix}"
                used_ids.add(sid)
                jc_3p = JunctionCluster(
                    host_chr=chrom, position=l_pos,
                    clip_direction="left",
                    clipped_seqs=l_seqs, consensus_clip=l_consensus,
                    n_reads=len(l_seqs),
                )
                site = InsertionSite(
                    site_id=sid,
                    host_chr=chrom,
                    junction_3p=jc_3p,
                    seed_3p=l_consensus,
                    pos_5p=l_pos,
                    confidence="low",
                )
                sites.append(site)

    log(f"  Found {len(sites)} candidate insertion sites")

    # Step 4: Validate — batched external tool calls for speed
    # Separate paired (bidirectional) vs single-direction sites
    paired_sites = [s for s in sites
                    if s.junction_5p is not None and s.junction_3p is not None]
    single_sites = [s for s in sites
                    if s.junction_5p is None or s.junction_3p is None]

    log(f"  Paired (bidirectional): {len(paired_sites)}, "
        f"Single-direction: {len(single_sites)}")

    # Condition 2 filter on paired sites: clips must be different
    cond2_passed: list[InsertionSite] = []
    for site in paired_sites:
        cmp_len = min(20, len(site.seed_5p), len(site.seed_3p))
        site.clips_are_different = site.seed_5p[:cmp_len] != site.seed_3p[:cmp_len]
        if not site.clips_are_different:
            continue
        cond2_passed.append(site)

    log(f"  After clip-difference filter: {len(cond2_passed)} paired sites")

    # Condition 3: batch host-mapping check (one minimap2 call)
    host_check_seqs: dict[str, str] = {}
    for site in cond2_passed:
        host_check_seqs[f"{site.site_id}_5p"] = site.seed_5p
        host_check_seqs[f"{site.site_id}_3p"] = site.seed_3p

    maps_to_host_set = _batch_check_maps_to_host(host_check_seqs, host_ref, workdir)
    log(f"  Host-mapping check: {len(maps_to_host_set)} clips map to host")

    validated_sites: list[InsertionSite] = []
    for site in cond2_passed:
        maps_5p = f"{site.site_id}_5p" in maps_to_host_set
        maps_3p = f"{site.site_id}_3p" in maps_to_host_set
        if site.junction_5p:
            site.junction_5p.maps_to_host = maps_5p
        if site.junction_3p:
            site.junction_3p.maps_to_host = maps_3p
        site.clips_not_in_host = not maps_5p and not maps_3p
        if not site.clips_not_in_host:
            continue
        site.is_validated = True
        site.confidence = "high"
        validated_sites.append(site)

    # Batch element-hit check for all validated + single-direction sites
    element_check_seqs: dict[str, str] = {}
    for site in validated_sites:
        element_check_seqs[f"{site.site_id}_5p"] = site.seed_5p
        element_check_seqs[f"{site.site_id}_3p"] = site.seed_3p
    for site in single_sites:
        clip_seq = site.seed_5p or site.seed_3p
        if clip_seq:
            element_check_seqs[f"{site.site_id}_clip"] = clip_seq

    if element_db and element_check_seqs:
        all_element_hits = _batch_check_element_hits(
            element_check_seqs, element_db, workdir, extra_dbs=extra_dbs)

        for site in validated_sites:
            hits_5p = all_element_hits.get(f"{site.site_id}_5p", [])
            hits_3p = all_element_hits.get(f"{site.site_id}_3p", [])
            if site.junction_5p:
                site.junction_5p.element_hits = hits_5p
            if site.junction_3p:
                site.junction_3p.element_hits = hits_3p
            site.has_element_hits = bool(hits_5p or hits_3p)

        for site in single_sites:
            clip_hits = all_element_hits.get(f"{site.site_id}_clip", [])
            if clip_hits:
                if site.junction_5p:
                    site.junction_5p.element_hits = clip_hits
                if site.junction_3p:
                    site.junction_3p.element_hits = clip_hits
                site.has_element_hits = True
                site.confidence = "medium"
                validated_sites.append(site)
                log(f"  {site.site_id}: single-direction but element hit "
                    f"({clip_hits[0]}) → medium confidence")

    for site in validated_sites:
        log(f"  {site.site_id}: {site.host_chr}:{site.pos_5p}-{site.pos_3p} "
            f"VALIDATED [{site.confidence}] "
            f"(5p={len(site.seed_5p)}bp, 3p={len(site.seed_3p)}bp, "
            f"element={'yes' if site.has_element_hits else 'no'})")

    # Sort: element-hit sites FIRST (primary), then paired > single, then read support
    # Element DB hits are the strongest signal for T-DNA insertion vs SV
    validated_sites.sort(key=lambda s: (
        0 if s.has_element_hits else 1,
        0 if s.confidence == "high" else (1 if s.confidence == "medium" else 2),
        -(s.junction_5p.n_reads if s.junction_5p else 0)
        - (s.junction_3p.n_reads if s.junction_3p else 0),
    ))

    log(f"  Validated: {len(validated_sites)} insertion sites "
        f"({sum(1 for s in validated_sites if s.has_element_hits)} with element hits)")

    return validated_sites


# ---------------------------------------------------------------------------
# T10: Pre-mask BED intersect (host-endogenous ortholog regions)
# ---------------------------------------------------------------------------

#: Tier value written onto sites that intersect a T9 host-endogenous mask
#: region.  Centralised here so downstream consumers (generate_report audit
#: appender, tests, future modules) can import the tag instead of duplicating
#: the string literal.  Changing this constant is a wire-format break — update
#: the audit TSVs and regression fixtures in lock-step.
MASKED_SOURCE_TAG = "FALSE_NEGATIVE_MASKED"


def _apply_mask_bed(sites: list, bed_path: Path) -> list:
    """Tag sites whose ``(chr, pos)`` falls inside a BED mask region.

    Sites are NOT dropped — their ``tier`` is changed to
    :data:`MASKED_SOURCE_TAG` (``FALSE_NEGATIVE_MASKED``) and ``mask_element``
    holds the element_name of the overlapping BED region.  This preserves the
    audit trail for regulatory review while preventing CANDIDATE promotion on
    host-endogenous ortholog regions (T9 E-01/E-02).

    Interval semantics
    ------------------
    Regions are treated as **half-open** ``[start, end)`` — i.e. ``start`` is
    inclusive and ``end`` is exclusive.  This matches BED convention and the
    ``start <= pos < end`` check below; a site at ``pos == end`` is NOT
    considered masked.

    BED format (from ``scripts/build_element_mask_bed.sh``):
    ``chr<TAB>start<TAB>end<TAB>name[<TAB>score]``.  An empty or missing BED
    file (including a BED with zero parseable regions) is a no-op.

    Parameters
    ----------
    sites
        A list of site records.  The list itself is returned (not copied); the
        records inside are mutated in place — see *Mutation contract* below.
    bed_path
        Path to the pre-computed mask BED.

    Mutation contract
    -----------------
    Accepts heterogeneous site records and mutates the *original* objects:

    * ``dict`` with ``chr``/``pos`` keys (test + pipeline JSON form) — sets
      ``s["tier"]`` and ``s["mask_element"]``.
    * Any object with ``host_chr`` + ``pos_5p`` attributes (e.g. the
      ``InsertionSite`` dataclass) — sets ``tier`` and ``mask_element`` via
      ``setattr``.  When ``pos_5p`` is falsy the helper falls back to
      ``pos_3p``.

    Callers that need the prior value should snapshot it before invoking this
    function; the helper does not save or restore previous ``tier`` values.

    Complexity
    ----------
    Regions are grouped by chromosome into a list sorted by ``start`` plus a
    parallel ``start`` array; per-site lookup uses :func:`bisect.bisect_right`
    and scans only intervals whose ``start <= pos``.  Runtime is therefore
    ``O(R log R + S log R)`` for ``R`` regions and ``S`` sites, instead of the
    ``O(S * R)`` worst case of a linear scan.
    """
    import bisect

    if not bed_path.exists() or bed_path.stat().st_size == 0:
        return sites
    raw_regions: dict[str, list[tuple[int, int, str]]] = {}
    for raw in bed_path.read_text().splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        cols = raw.split("\t")
        if len(cols) < 3:
            continue
        try:
            start = int(cols[1])
            end = int(cols[2])
        except ValueError:
            continue
        chrom = cols[0]
        name = cols[3] if len(cols) > 3 else "masked"
        raw_regions.setdefault(chrom, []).append((start, end, name))
    if not raw_regions:
        return sites
    # Sort regions by start and stash a parallel starts array for bisect.
    # NOTE: the intervals may still overlap each other; bisect narrows the
    # scan to candidates with start <= pos, then we check end > pos.
    index: dict[str, tuple[list[int], list[tuple[int, int, str]]]] = {}
    for chrom, ivs in raw_regions.items():
        ivs.sort(key=lambda t: t[0])
        index[chrom] = ([iv[0] for iv in ivs], ivs)
    for s in sites:
        if isinstance(s, dict):
            chrom = s.get("chr")
            pos = s.get("pos")
        else:
            chrom = getattr(s, "host_chr", None)
            pos = getattr(s, "pos_5p", None)
            if not pos:
                pos = getattr(s, "pos_3p", None)
        if chrom is None or pos is None:
            continue
        entry = index.get(chrom)
        if entry is None:
            continue
        starts, ivs = entry
        hi = bisect.bisect_right(starts, pos)
        # Walk backwards through candidates whose start <= pos; because the
        # regions in the mask BED are short (<= a few kb) and only a handful
        # overlap any given pos, this inner loop is effectively O(1).
        for idx in range(hi - 1, -1, -1):
            start, end, name = ivs[idx]
            if end <= pos:
                # Half-open: pos >= end means outside; keep scanning earlier
                # intervals because an earlier start may still extend past pos.
                continue
            # start <= pos < end -> match.
            if isinstance(s, dict):
                s["tier"] = MASKED_SOURCE_TAG
                s["mask_element"] = name
            else:
                setattr(s, "tier", MASKED_SOURCE_TAG)
                setattr(s, "mask_element", name)
            break
    return sites


# ---------------------------------------------------------------------------
# Transgene-positive identification
# ---------------------------------------------------------------------------
# Strategy: Instead of asking "is this clip in the host?" (fails with wrong
# cultivar reference), ask "is this clip from a transgene?"
# Uses blastn-short against transgene_db (element DB + UniVec vectors).
# Any transgene hit = assembly target. No hit = skip.
#
# NOTE: TierResult dataclass moved to scripts/s05/primitives.py
# (Issue #4 Session 1) and re-imported at the top of this module.


# ---------------------------------------------------------------------------
# Phase 1.5 (Classify) functions moved to scripts/s05/classify.py
# (Issue #4 Session 4) and re-imported at the top of this module for
# backward compatibility. _SRC_TIER / _TIER2_SRCS / _UNKNOWN_SRC_WARNED
# module state migrated with them. _TIER2_SRCS is derived from _SRC_TIER
# to stay in sync with gmo_combined_db_v2.fa cd-hit clustering; see
# scripts/s05/classify.py for the full audit trail comment.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Fallback: Parse junctions.tsv from step 6
# ---------------------------------------------------------------------------

def parse_legacy_junctions(tsv_path: Path) -> list[LegacyJunction]:
    """Parse junctions.tsv from step 6 as fallback."""
    raw: list[LegacyJunction] = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            raw.append(LegacyJunction(
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

    # Deduplicate
    best: dict[tuple[str, str, int], LegacyJunction] = {}
    for j in raw:
        key = (j.contig_name, j.host_chr, j.junction_pos_host)
        if key not in best or j.host_mapq > best[key].host_mapq:
            best[key] = j
    return list(best.values())


def _extract_seeds_at_positions(
    site: InsertionSite,
    host_bam: Path,
    positions: list[int],
    chrom: str,
    min_clip: int,
    window: int,
) -> None:
    """Extract soft-clip seeds at known junction positions."""
    bam = pysam.AlignmentFile(str(host_bam), "rb")

    for jpos in positions:
        right_clips: list[str] = []
        left_clips: list[str] = []

        for read in bam.fetch(chrom, max(0, jpos - 200), jpos + 200):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            cigar = read.cigartuples
            if cigar is None:
                continue
            seq = read.query_sequence
            if seq is None:
                continue

            if cigar[-1][0] == 4 and cigar[-1][1] >= min_clip:
                if abs(read.reference_end - jpos) <= window:
                    right_clips.append(seq[-cigar[-1][1]:])

            if cigar[0][0] == 4 and cigar[0][1] >= min_clip:
                if abs(read.reference_start - jpos) <= window:
                    left_clips.append(seq[:cigar[0][1]])

        if right_clips and len(right_clips) >= 2:
            consensus = _build_consensus(right_clips, "right")
            if len(consensus) >= min_clip:
                site.junction_5p = JunctionCluster(
                    host_chr=chrom, position=jpos, clip_direction="right",
                    clipped_seqs=right_clips, consensus_clip=consensus,
                    n_reads=len(right_clips),
                )
                if not site.seed_5p or len(consensus) > len(site.seed_5p):
                    site.seed_5p = consensus

        if left_clips and len(left_clips) >= 2:
            consensus = _build_consensus(left_clips, "left")
            if len(consensus) >= min_clip:
                site.junction_3p = JunctionCluster(
                    host_chr=chrom, position=jpos, clip_direction="left",
                    clipped_seqs=left_clips, consensus_clip=consensus,
                    n_reads=len(left_clips),
                )
                if not site.seed_3p or len(consensus) > len(site.seed_3p):
                    site.seed_3p = consensus

    bam.close()


def legacy_junctions_to_sites(
    junctions: list[LegacyJunction],
    host_bam: Path,
    min_clip: int = 15,
    window: int = 10,
) -> list[InsertionSite]:
    """Convert step 6 junctions to InsertionSites with soft-clip seeds."""
    if not junctions:
        return []

    # Group by chromosome
    by_chr: dict[str, list[LegacyJunction]] = {}
    for j in junctions:
        by_chr.setdefault(j.host_chr, []).append(j)

    sites: list[InsertionSite] = []
    used_ids: set[str] = set()

    for chrom, juncs in by_chr.items():
        juncs.sort(key=lambda j: j.junction_pos_host)

        # Pair within 50kb
        used = set()
        for i, j1 in enumerate(juncs):
            if i in used:
                continue
            paired_idx = None
            for k, j2 in enumerate(juncs):
                if k in used or k == i:
                    continue
                if abs(j2.junction_pos_host - j1.junction_pos_host) <= 50000:
                    paired_idx = k
                    break

            anchor_pos = j1.junction_pos_host
            sid = f"insertion_{chrom}_{anchor_pos}"
            if sid in used_ids:
                suffix = 2
                while f"{sid}_{suffix}" in used_ids:
                    suffix += 1
                sid = f"{sid}_{suffix}"
            used_ids.add(sid)
            site = InsertionSite(
                site_id=sid,
                host_chr=chrom,
            )

            positions = [j1.junction_pos_host]
            if paired_idx is not None:
                used.add(i)
                used.add(paired_idx)
                positions.append(juncs[paired_idx].junction_pos_host)
            else:
                used.add(i)

            site.pos_5p = min(positions)
            site.pos_3p = max(positions) if len(positions) > 1 else 0

            # Extract seeds from soft-clips at junction positions
            _extract_seeds_at_positions(site, host_bam, positions, chrom,
                                        min_clip, window)
            sites.append(site)

    return sites
