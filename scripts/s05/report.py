"""Phase 5 — Report + Stats (Issue #4, Session 6).

Extracted verbatim from ``scripts/s05_insert_assembly.py``. Owns the
per-site report generation (insertion_report.tsv with verdict labels)
and the sample-level stats TSV. Stage 8 in the s05 DAG.

Functions
---------
- ``generate_report`` — writes per-site insertion_report.tsv with
  embedded CANDIDATE/FALSE_POSITIVE/UNKNOWN verdict logic. NOTE: the
  verdict logic is duplicated with ``scripts/s05/verdict.py::compute_verdict``;
  Session 7 cleanup may consolidate.
- ``write_stats`` — sample-level summary TSV (insertion_sites count,
  candidate_read_pairs, per-site result rows).

Replaces the v1.0 ``annotate_report.py`` shim's monolith re-export of
these two symbols.
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

from .primitives import (
    log,
    read_fasta,
    InsertionSite,
)
from .verdict import (
    FilterEvidence,
    VerdictRules,
    compute_verdict,
    select_canonical_elements,
)
from .filter_a_host import (
    INSERT_HOST_MIN_PIDENT,
    _blast_insert_vs_host,
)
from .filter_b_flank import (
    CONSTRUCT_FLANK_SLOP,
    _site_overlaps_flanking,
)
from .filter_c_chimeric import (
    _check_chimeric_assembly,
)
from .filter_d_altlocus import (
    _check_construct_host_coverage,
)


def generate_report(
    insert_fasta: Path,
    annotation_tsv: Path,
    border_tsv: Path,
    site: InsertionSite,
    # NOTE: verdict logic (CANDIDATE/FALSE_POSITIVE/UNKNOWN) is embedded here.
    # If a verdict-only mode is needed later, extract into compute_verdict().
    n_rounds: int,
    status: str,
    output_dir: Path,
    host_endogenous_ids: set[str] | None = None,
    host_ref: Path | None = None,
    element_db: Path | None = None,
    threads: int = 4,
    construct_flanking: list[tuple[str, int, int]] | None = None,
    rules: VerdictRules | None = None,
) -> tuple[Path, str]:
    """Generate human-readable linear map report."""
    report_path = output_dir / f"{site.site_id}_report.txt"

    # Read insert sequence
    seqs = read_fasta(insert_fasta)
    if not seqs:
        with open(report_path, "w") as fh:
            fh.write("No insert assembled.\n")
        return report_path, "NO_ASSEMBLY"

    insert_name = list(seqs.keys())[0]
    insert_seq = seqs[insert_name]
    insert_len = len(insert_seq)
    n_count = insert_seq.upper().count("N")

    # Read annotation hits (with orientation from s_start/s_end)
    # elements: (q_start, q_end, element_name, s_strand, source)
    elements: list[tuple[int, int, str, str, str]] = []
    orientation_by_elem: dict[str, set[str]] = defaultdict(set)
    # Best BLAST identity (pident, 0-100) seen for each element. Feeds the
    # canonical-triplet identity gate below.  A legacy annotation TSV without an
    # `identity` column cannot be gated — do not silently veto every element.
    identity_by_elem: dict[str, float] = {}
    has_identity_col = True
    if annotation_tsv.exists():
        with open(annotation_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            has_identity_col = "identity" in (reader.fieldnames or [])
            for row in reader:
                try:
                    q_start = int(row["q_start"])
                    q_end = int(row["q_end"])
                    s_start = int(row["s_start"])
                    s_end = int(row["s_end"])
                    elem = row["element"]
                    source = row.get("source", "element_db")
                    # Filter to current insert only
                    if row["query"] != insert_name:
                        continue
                    strand = "+" if s_start < s_end else "-"
                    elements.append((q_start, q_end, elem, strand, source))
                    orientation_by_elem[elem].add(strand)
                except (ValueError, KeyError):
                    continue
                try:
                    identity = float(row.get("identity") or 0.0)
                except (TypeError, ValueError):
                    identity = 0.0
                identity_by_elem[elem] = max(
                    identity_by_elem.get(elem, 0.0), identity,
                )
    elements.sort(key=lambda x: x[0])

    # Deduplicate overlapping elements (keep longest per region)
    deduped: list[tuple[int, int, str, str, str]] = []
    for start, end, elem, strand, source in elements:
        overlap = False
        for i, (ds, de, dn, _ds, _src) in enumerate(deduped):
            if start < de and end > ds:
                if (end - start) > (de - ds):
                    deduped[i] = (start, end, elem, strand, source)
                overlap = True
                break
        if not overlap:
            deduped.append((start, end, elem, strand, source))
    elements = sorted(deduped, key=lambda x: x[0])

    # Count element occurrences for multi-construct detection
    elem_counts: Counter = Counter()
    for _, _, elem, _, _ in elements:
        elem_counts[elem] += 1

    # ---- Build element-level evidence (no I/O) ----
    host_endo = host_endogenous_ids or set()
    unique_elems = set(elem_counts)
    foreign_elems: set[str] = set()
    endo_elems: set[str] = set()
    for elem in unique_elems:
        if elem in host_endo:
            endo_elems.add(elem)
        else:
            foreign_elems.add(elem)
    # sources_by_element: last-seen source tag per element name.
    sources_by_element: dict[str, str] = {
        elem: source for _, _, elem, _, source in elements
    }

    # ---- Gather BLAST filter evidence (lazy: only run when useful) ----
    # Evidence is collected for CANDIDATE sites and UNKNOWN sites (host-only
    # reclassification). All four filter values default to 0/None so that
    # compute_verdict receives a complete FilterEvidence for every path.
    host_fraction = 0.0
    host_bp = 0
    largest_gap = 0
    flanking_hit_str = ""           # human-readable for report diagnostics
    flanking_hit_tup: tuple[str, int, int] | None = None  # structured for verdict
    is_chimeric = False
    off_target_chrs: list[tuple[str, int]] = []
    construct_frac = 0.0
    combined_frac = 0.0

    # Determine whether this site is worth running BLASTs on.
    # Sites that are all-host-endogenous skip BLAST, matching old behaviour;
    # canonical_triplet override (highest priority in compute_verdict) still
    # works because ev.host_fraction defaults to 0.0 < cand_host_fraction_max.
    _needs_blast = bool(unique_elems and foreign_elems) or not unique_elems
    # "not unique_elems" → UNKNOWN path, needs BLAST for host-only reclassification.
    # "foreign_elems" → CANDIDATE path, needs BLAST for all four FP filters.

    if _needs_blast and host_ref is not None:
        # Filter A evidence: host-fraction + largest non-host gap
        host_fraction, host_bp, _, largest_gap = _blast_insert_vs_host(
            insert_fasta, host_ref, output_dir, threads=threads,
        )
        if not unique_elems:
            log(f"  UNKNOWN reclassification: BLASTed {insert_fasta.name} vs host")

    if _needs_blast and construct_flanking:
        # Filter B evidence: find first overlapping flanking region
        site_pos = site.pos_5p
        _, flanking_hit_str = _site_overlaps_flanking(
            site.host_chr, site_pos, construct_flanking,
        )
        for chrom, lo, hi in construct_flanking:
            slop = CONSTRUCT_FLANK_SLOP
            if chrom == site.host_chr and lo - slop <= site_pos <= hi + slop:
                # Store slop-expanded coords so compute_verdict's simple
                # lo <= site_pos <= hi boundary check is equivalent.
                flanking_hit_tup = (chrom, lo - slop, hi + slop)
                break

    if _needs_blast and host_ref is not None:
        # Filter C evidence: off-target chromosome spans
        is_chimeric, off_target_chrs = _check_chimeric_assembly(
            insert_fasta, host_ref, site.host_chr, output_dir, threads=threads,
        )

    if _needs_blast and host_ref is not None and element_db is not None:
        # Filter D evidence: construct + host combined coverage
        _, construct_frac, _, combined_frac = _check_construct_host_coverage(
            insert_fasta, element_db,
            host_fraction, host_bp, insert_len, n_count,
            output_dir, threads=threads,
        )

    # ---- Delegate verdict to compute_verdict (Issue #3 full wire-in) ----
    # Build a FilterEvidence bundle from all gathered evidence and call the
    # pure function.  This replaces the old inline if/else verdict chain plus
    # the _apply_canonical_override post-hoc call.
    _verdict_rules = rules if rules is not None else VerdictRules()

    # Rule 1 outranks every FP filter, so only elements whose best BLAST hit
    # clears `canonical_triplet_min_identity` may complete a canonical triplet.
    if has_identity_col:
        matched_canonical = select_canonical_elements(
            {elem: identity_by_elem.get(elem, 0.0) for elem in unique_elems},
            _verdict_rules.canonical_triplet_min_identity,
        )
        for elem in sorted(unique_elems - matched_canonical):
            log(f"  canonical gate: {elem} best identity "
                f"{identity_by_elem.get(elem, 0.0):.1f}% below "
                f"{_verdict_rules.canonical_triplet_min_identity:.0%} — "
                f"not eligible for triplet promotion")
    else:
        if unique_elems:
            log("  canonical gate: annotation TSV has no 'identity' column — "
                "triplet promotion left ungated")
        matched_canonical = set(unique_elems)

    _ev = FilterEvidence(
        elements=list(unique_elems),
        host_bp=host_bp,
        host_fraction=host_fraction,
        largest_gap=largest_gap,
        flanking_hit=flanking_hit_tup,
        off_target_chrs=off_target_chrs,
        construct_frac=construct_frac,
        combined_frac=combined_frac,
        is_chimeric=is_chimeric,
        site_chr=site.host_chr,
        site_pos=site.pos_5p,
        matched_canonical=matched_canonical,
        sources_by_element=sources_by_element,
        host_endogenous_elements=endo_elems,
    )
    verdict, verdict_reason = compute_verdict(_ev, _verdict_rules)

    # For report diagnostics: restore the human-readable flanking string if
    # compute_verdict returned a flanking FP (flanking_hit_str may differ from
    # the structured tuple representation used internally).
    # The reporting block below uses flanking_hit_str for display; verdict_reason
    # already contains the full human-readable text from compute_verdict.

    # Build linear map with orientation arrows
    linear_parts = []
    for _, _, elem, strand, _ in elements:
        short_name = elem.split("|")[0] if "|" in elem else (
            elem.split("_")[0] if "_" in elem else elem)
        arrow = "→" if strand == "+" else "←"
        linear_parts.append(f"[{short_name}{arrow}]")

    # Detect head-to-head / tandem arrangement
    # Head-to-head: same element found in BOTH orientations
    structure = "single-copy"
    bidirectional_elems = {e for e, strands in orientation_by_elem.items()
                          if "+" in strands and "-" in strands
                          and len(e) > 30}  # skip short amplicons
    if bidirectional_elems:
        structure = "head-to-head 2-copy T-DNA"
    elif len(elements) > 3:
        if elem_counts.most_common(1)[0][1] >= 2:
            structure = f"multi-copy (≥{elem_counts.most_common(1)[0][1]} copies)"

    # Read border hits.
    # `border_hits.tsv` is written once per s05 run, for ALL assembled inserts —
    # its subject column (2) names the insert each hit belongs to. Counting every
    # line made each site's report claim the borders found on every other site
    # (rice_G281 Chr3:16,439,674 reported "T-DNA borders found: 10" while all ten
    # hits belonged to three other inserts).
    n_borders = 0
    if border_tsv.exists():
        with open(border_tsv) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 2 and cols[1] == insert_name:
                    n_borders += 1

    # Determine deletion size
    deletion_size = abs(site.pos_3p - site.pos_5p) if site.pos_3p > 0 else 0

    with open(report_path, "w") as fh:
        fh.write("=" * 70 + "\n")
        fh.write("RedGene Insert Assembly & Annotation Report\n")
        fh.write("=" * 70 + "\n")
        pos_str = f"{site.host_chr}:{site.pos_5p:,}"
        if site.pos_3p > 0 and site.pos_3p != site.pos_5p:
            pos_str += f"-{site.pos_3p:,}"
        fh.write(f"Insertion site: {pos_str}")
        if deletion_size > 0:
            fh.write(f" ({deletion_size}bp deletion)")
        fh.write("\n")
        fh.write(f"Insert length: {insert_len:,} bp")
        if n_count > 0:
            fh.write(f" ({n_count} unresolved N's)")
        fh.write("\n")
        fh.write(f"Assembly status: {status.upper()} (round {n_rounds})\n")
        fh.write(f"Structure: {structure}\n")
        fh.write(f"Verdict: {verdict}")
        if verdict_reason:
            fh.write(f" — {verdict_reason}")
        fh.write("\n")
        if n_borders > 0:
            fh.write(f"T-DNA borders found: {n_borders}\n")
        fh.write("\n")

        if linear_parts:
            fh.write("--- Linear Map ---\n")
            fh.write("5' host --")
            line = ""
            for part in linear_parts:
                if len(line) + len(part) > 60:
                    fh.write(line + "\n         ")
                    line = ""
                line += part + "--"
            fh.write(line + " 3' host\n")
        else:
            fh.write("--- No element annotations found ---\n")

        fh.write("=" * 70 + "\n")

        # Detailed element list
        if elements:
            fh.write("\nDetailed element positions:\n")
            fh.write(f"{'Start':>8}  {'End':>8}  {'Dir':>3}  {'Src':>5}  {'Element'}\n")
            fh.write("-" * 70 + "\n")
            for start, end, elem, strand, source in elements:
                src_tag = "local" if source == "element_db" else "NCBI"
                fh.write(f"{start:>8}  {end:>8}  {strand:>3}  {src_tag:>5}  {elem}\n")

        # Multi-copy detection
        multi = {e: c for e, c in elem_counts.items() if c >= 2}
        if multi:
            fh.write(f"\nMulti-copy elements detected:\n")
            for elem, count in multi.items():
                fh.write(f"  {elem}: {count} copies\n")

        # Head-to-head evidence
        if bidirectional_elems:
            fh.write(f"\nHead-to-head evidence (same element in both orientations):\n")
            for elem in sorted(bidirectional_elems):
                fh.write(f"  {elem}\n")

        # Foreign vs host-endogenous element breakdown
        if host_endo and unique_elems:
            fh.write(f"\nForeign elements (not in host genome): {len(foreign_elems)}\n")
            for elem in sorted(foreign_elems):
                fh.write(f"  + {elem}\n")
            if endo_elems:
                fh.write(f"Host-endogenous elements "
                         f"(excluded at DB-level BLAST): {len(endo_elems)}\n")
                for elem in sorted(endo_elems):
                    fh.write(f"  - {elem}\n")

        # Post-assembly filter diagnostics
        if host_fraction > 0:
            fh.write(f"\nHost-fraction analysis: {host_fraction:.1%} host "
                     f"({host_bp:,}/{insert_len - n_count:,}bp at "
                     f"≥{INSERT_HOST_MIN_PIDENT}% identity), "
                     f"largest non-host gap: {largest_gap:,}bp\n")
        if flanking_hit_str:
            fh.write(f"Construct-flanking overlap: site overlaps {flanking_hit_str}\n")
        if off_target_chrs:
            off_str = ", ".join(f"{c}:{bp}bp" for c, bp in off_target_chrs)
            tag = " → chimeric" if is_chimeric else ""
            fh.write(f"Off-target chromosomes in assembly: {off_str}{tag}\n")
        if construct_frac > 0:
            fh.write(f"Construct coverage: {construct_frac:.1%}, "
                     f"combined (construct+host): {combined_frac:.1%}\n")

    log(f"  Report written: {report_path} [{verdict}]")
    return report_path, verdict


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
            fh.write(f"{prefix}_verdict\t{r.get('verdict', 'NO_ASSEMBLY')}\n")
