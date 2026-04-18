"""Pure verdict logic extracted from s05_insert_assembly.generate_report().

T6 (v1.0 MVP) — isolates the CANDIDATE / FALSE_POSITIVE / UNKNOWN decision
tree into a side-effect-free function so it can be unit-tested without BLAST,
SPAdes, or file I/O.

Issue #3 (v1.1 full wire-in) completes the call-site replacement inside
generate_report and reconciles the filter priority between the old inline code
and this module.

Filter priority (canonical, first-match wins)
=============================================

  1.  canonical_triplet promotion (runs FIRST).
      If every element in a configured triplet is present in
      ev.matched_canonical AND ev.host_fraction < cand_host_fraction_max,
      return CANDIDATE immediately.  Rationale: a canonical transgene marker
      set (e.g. bar + P-CaMV35S + T-ocs) is strong evidence of a real
      insertion regardless of host fraction or other FP signals.

  1.5 host_endogenous (all elements are host-derived → FP).
      If every annotated element matched the host genome at Tier 1/2 BLAST
      threshold, the assembly is host-genomic noise.  Placed after Rule 1 so
      that a canonical triplet still promotes even if every element happens to
      have a host BLAST hit (unusual but possible with cultivar-derived
      promoters like P-Act1).

  2.  Filter B — construct-flanking overlap (FP).
      If the insertion site falls within a region where the construct
      reference contains host genomic DNA at its ends (detected by
      construct→host BLAST at the call site), the site is a false detection
      of the host–construct junction rather than a true insertion.

  3.  Filter C — multi-locus chimeric assembly (FP).
      If the assembled insert's host-aligned portions span ≥2 chromosomes
      other than the site chromosome, the assembler merged reads from
      unrelated loci via shared element homology.

  4.  Filter D — construct+host coverage (FP).
      If construct coverage ≥ fp_construct_frac_min AND combined (construct +
      host) coverage ≥ fp_combined_frac_min, the insert is host genomic DNA
      with construct-element homology, not a real T-DNA.

  5.  Filter A — host fraction + small gap (FP).
      If host_fraction ≥ fp_host_fraction_min AND largest_gap <
      fp_largest_gap_max, the insert is essentially a host genomic fragment
      with a negligible foreign stretch.
      IMPORTANT: A site with HIGH host_fraction but a LARGE gap (≥500 bp) is
      a genuine junction contig (host→T-DNA→host spans host DNA on both
      sides).  Filter A must come AFTER B/C/D so that chimeric assemblies
      with large gaps are still caught by Filter C.

  6.  Elements present (all FP filters survived) → CANDIDATE.
      If elements are annotated and no FP filter fired, the site is
      CANDIDATE.  host_fraction is NOT re-checked here; Rule 5 already
      handles the "high host_fraction + small gap" case.  The canonical
      counter-example is rice_G281 Chr3:16,439,674 (87.4% host, 1,024 bp
      gap): Filter A does not fire (gap ≥ 500 bp), so the site is CANDIDATE.

  7.  No elements + host-only insert → FP reclassification.
      If no elements are annotated and the insert is mostly host DNA with
      negligible construct match, reclassify UNKNOWN → FALSE_POSITIVE.

  8.  Fallthrough → UNKNOWN.

The companion `config_loader.load_verdict_rules()` materialises the thresholds
from config.yaml.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FilterEvidence:
    """Immutable bundle of post-assembly evidence feeding the verdict tree.

    The caller (generate_report in s05_insert_assembly.py) assembles this from
    BLAST results, flanking-region checks, chimera detection, and element
    annotation. Keeping it a plain dataclass means compute_verdict() can be
    exercised by hand in tests without stubbing any external tools.
    """

    # Annotated element names on the assembled insert (may be empty).
    elements: list[str] = field(default_factory=list)
    # Host BLAST coverage on the insert.
    host_bp: int = 0
    host_fraction: float = 0.0
    # Longest non-host stretch on the insert (bp).
    largest_gap: int = 0
    # ("chrom", lo, hi) of a construct-flanking hit overlapping the site, or None.
    flanking_hit: Optional[tuple[str, int, int]] = None
    # Off-target chromosomes touched by the assembled insert: [(chrom, bp), ...].
    off_target_chrs: list[tuple[str, int]] = field(default_factory=list)
    # Construct-element coverage on the insert (0..1).
    construct_frac: float = 0.0
    # host_fraction + construct_frac (non-overlapping merge).
    combined_frac: float = 0.0
    # DEPRECATED (Issue #11 M-3): `is_chimeric` is informational only; the
    # chimeric filter (rule 3 in compute_verdict) is driven entirely by
    # len(off_target_chrs) >= rules.fp_off_target_chrs_min. Kept on the
    # dataclass for backward-compatible dict-expansion in existing tests;
    # slated for removal in v1.1 once callers migrate off the kwarg.
    is_chimeric: bool = False
    # Insertion-site coordinates, used to evaluate flanking overlap messages.
    site_chr: str = ""
    site_pos: int = 0
    # Element names that intersect one of the canonical triplets in rules.
    matched_canonical: set[str] = field(default_factory=set)
    # Element → source-tier tag (sample_contig / element_db / common_payload / univec).
    sources_by_element: dict[str, str] = field(default_factory=dict)
    # Element names excluded at Tier 1/2 DB-level host BLAST (Phase 4 host filter).
    host_endogenous_elements: set[str] = field(default_factory=set)


@dataclass(frozen=True)
class VerdictRules:
    """Thresholds + canonical triplets loaded from config.yaml.

    Defaults mirror the compiled-in constants in s05_insert_assembly.py
    (INSERT_HOST_FRACTION, INSERT_MIN_FOREIGN_GAP, UNKNOWN_HOST_MIN_FRACTION,
    UNKNOWN_MAX_CONSTRUCT_FRAC, ...) so extracting the logic does not change
    verdicts for any existing sample.
    """

    cand_host_fraction_max: float = 0.80
    cand_largest_gap_min: int = 500
    fp_host_fraction_min: float = 0.80
    fp_largest_gap_max: int = 500
    fp_off_target_chrs_min: int = 2
    fp_combined_frac_min: float = 0.85
    fp_construct_frac_min: float = 0.25
    unknown_to_fp_host_fraction_min: float = 0.85
    unknown_to_fp_construct_frac_max: float = 0.05
    canonical_triplets: dict[str, set[str]] = field(default_factory=dict)
    canonical_triplet_min_identity: float = 0.90


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_matching_triplet(
    matched: set[str],
    triplets: dict[str, set[str]],
) -> Optional[str]:
    """Return the name of any canonical triplet fully covered by matched, else None."""
    for name, triplet in triplets.items():
        if triplet and triplet.issubset(matched):
            return name
    return None


def _flanking_overlaps_site(
    flanking_hit: Optional[tuple[str, int, int]],
    site_chr: str,
    site_pos: int,
) -> bool:
    """True iff flanking_hit is set and its interval covers the site."""
    if flanking_hit is None:
        return False
    chrom, lo, hi = flanking_hit
    if site_chr and chrom != site_chr:
        return False
    return lo <= site_pos <= hi


# ---------------------------------------------------------------------------
# Pure verdict function
# ---------------------------------------------------------------------------

def compute_verdict(
    ev: FilterEvidence,
    rules: VerdictRules,
) -> tuple[str, str]:
    """Classify an insertion site as CANDIDATE / FALSE_POSITIVE / UNKNOWN.

    Pure function — no I/O, no subprocess, no BLAST. Matches the decision
    tree currently embedded in s05_insert_assembly.generate_report()
    (lines 3197-3514). Rule priority (first match wins):

      1. canonical-triplet match + host_fraction < cand_host_fraction_max
         -> CANDIDATE.  A triplet matches when every element in the
         configured triplet (e.g. ``{bar, P-CaMV35S, T-ocs}``) is present
         in ``ev.matched_canonical``; see ``_find_matching_triplet``.
      1.5. all annotated elements are host-endogenous (Phase 4 host filter)
         -> FALSE_POSITIVE (host-endogenous)
      2. site overlaps a construct-flanking region
         -> FALSE_POSITIVE (flanking)
      3. off_target_chrs >= fp_off_target_chrs_min
         -> FALSE_POSITIVE (chimeric).  The boolean ``ev.is_chimeric`` is
         informational only; this rule is driven entirely by the
         ``off_target_chrs`` count (Issue #11 M-3).
      4. construct_frac >= fp_construct_frac_min AND combined_frac >= fp_combined_frac_min
         -> FALSE_POSITIVE (construct+host)
      5. host_fraction >= fp_host_fraction_min AND largest_gap < fp_largest_gap_max
         -> FALSE_POSITIVE (host+gap)
      6. elements present (survived all FP filters above)
         -> CANDIDATE.  host_fraction is NOT re-checked here because Rule 5
         already handles the high-host + small-gap case.  High host_fraction
         with a large gap is a genuine junction contig (e.g. rice_G281 Chr3
         87.4% host, 1,024 bp gap).
      7. no elements AND host_fraction >= unknown_to_fp_host_fraction_min
         AND construct_frac <= unknown_to_fp_construct_frac_max
         -> FALSE_POSITIVE (host-only)
      8. else -> UNKNOWN

    Reason strings are currently English-only. v1.1 will route them through
    a ``REASON_KEYS`` lookup (see Issue #11 M-5) so KO translations can be
    swapped in without touching rule logic.
    """
    # Rule 1: canonical triplet promotion (runs first so CRISPR/marker-positive
    # assemblies with some host contamination still get flagged CANDIDATE).
    triplet_name = _find_matching_triplet(
        ev.matched_canonical, rules.canonical_triplets,
    )
    if triplet_name is not None and ev.host_fraction < rules.cand_host_fraction_max:
        return ("CANDIDATE",
                f"canonical_triplet[{triplet_name}] matched "
                f"({sorted(ev.matched_canonical)}); host_fraction="
                f"{ev.host_fraction:.0%} below {rules.cand_host_fraction_max:.0%}")

    # Rule 1.5: Phase 4 host-endogenous filter.
    # If every unique annotated element was excluded at Tier 1/2 DB-level
    # host BLAST, the whole assembly is host-endogenous noise.
    if ev.elements and ev.host_endogenous_elements:
        unique_elems = set(ev.elements)
        if unique_elems and unique_elems <= ev.host_endogenous_elements:
            return ("FALSE_POSITIVE",
                    "all annotated elements are host-endogenous "
                    "(matched host genome at Tier 1/2 BLAST threshold)")

    # Rule 2: construct-flanking overlap (Filter B).
    if _flanking_overlaps_site(ev.flanking_hit, ev.site_chr, ev.site_pos):
        assert ev.flanking_hit is not None
        chrom, lo, hi = ev.flanking_hit
        return ("FALSE_POSITIVE",
                f"site {ev.site_chr}:{ev.site_pos:,} overlaps construct-flanking "
                f"region {chrom}:{lo:,}-{hi:,} "
                f"\u2014 host DNA in construct reference")

    # Rule 3: multi-locus chimeric assembly (Filter C).
    if len(ev.off_target_chrs) >= rules.fp_off_target_chrs_min:
        off_str = ", ".join(f"{c}:{bp}bp" for c, bp in ev.off_target_chrs)
        return ("FALSE_POSITIVE",
                f"chimeric assembly \u2014 host-aligned portions span "
                f"{len(ev.off_target_chrs)} off-target chromosomes ({off_str})")

    # Rule 4: construct + host coverage fully explains insert (Filter D).
    if (ev.construct_frac >= rules.fp_construct_frac_min
            and ev.combined_frac >= rules.fp_combined_frac_min):
        return ("FALSE_POSITIVE",
                f"insert fully explained by construct ({ev.construct_frac:.0%}) "
                f"+ host ({ev.host_fraction:.0%}) = {ev.combined_frac:.0%} combined "
                f"coverage \u2014 host genomic DNA with construct-element homology")

    # Rule 5: host fraction high + no meaningful foreign gap (Filter A).
    if (ev.host_fraction >= rules.fp_host_fraction_min
            and ev.largest_gap < rules.fp_largest_gap_max):
        return ("FALSE_POSITIVE",
                f"assembled insert is {ev.host_fraction:.0%} host genome "
                f"({ev.host_bp:,}bp) with only {ev.largest_gap}bp non-host gap "
                f"(need \u2265{rules.fp_largest_gap_max}bp for real T-DNA)")

    # Rule 6: elements annotated and all FP filters survived → CANDIDATE.
    # NOTE: host_fraction is intentionally NOT checked here.  Filter A (Rule 5)
    # already handles the "high host_fraction + small gap" FP case.  A site
    # with high host_fraction but a LARGE gap (≥fp_largest_gap_max) is a genuine
    # junction contig (host→T-DNA→host) and must be CANDIDATE.  The canonical
    # example is rice_G281 Chr3:16,439,674 (87.4% host, 1,024 bp gap): Filter A
    # does not fire (gap ≥ 500 bp), so the site must be CANDIDATE, not UNKNOWN.
    if ev.elements:
        return ("CANDIDATE",
                f"{len(ev.elements)} element(s) annotated; host_fraction="
                f"{ev.host_fraction:.0%} (all FP filters passed)")

    # Rule 7: no elements, insert looks host-only (UNKNOWN reclassification).
    if (not ev.elements
            and ev.host_fraction >= rules.unknown_to_fp_host_fraction_min
            and ev.construct_frac <= rules.unknown_to_fp_construct_frac_max):
        return ("FALSE_POSITIVE",
                f"no element annotations; insert is {ev.host_fraction:.0%} host "
                f"genome ({ev.host_bp:,}bp) with {ev.construct_frac:.0%} construct "
                f"\u2014 host genomic DNA")

    # Rule 8: fallthrough.
    if not ev.elements:
        return ("UNKNOWN", "no element annotations")
    return ("UNKNOWN",
            f"{len(ev.elements)} element(s) annotated but evidence inconclusive "
            f"(host_fraction={ev.host_fraction:.0%})")
