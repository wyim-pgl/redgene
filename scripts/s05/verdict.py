"""Pure verdict logic extracted from s05_insert_assembly.generate_report().

T6 (v1.0 MVP) — isolates the CANDIDATE / FALSE_POSITIVE / UNKNOWN decision
tree into a side-effect-free function so it can be unit-tested without BLAST,
SPAdes, or file I/O.

The companion `config_loader.load_verdict_rules()` materialises the thresholds
from config.yaml. Call-site replacement inside generate_report is deferred to
v1.1 (per team-consensus §2.2); v1.0 ships the pure function plus tests only.
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
    # True when the chimera check (Filter C) flagged the assembly.
    is_chimeric: bool = False
    # Insertion-site coordinates, used to evaluate flanking overlap messages.
    site_chr: str = ""
    site_pos: int = 0
    # Element names that intersect one of the canonical triplets in rules.
    matched_canonical: set[str] = field(default_factory=set)
    # Element → source-tier tag (sample_contig / element_db / common_payload / univec).
    sources_by_element: dict[str, str] = field(default_factory=dict)


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
         -> CANDIDATE
      2. site overlaps a construct-flanking region
         -> FALSE_POSITIVE (flanking)
      3. off_target_chrs >= fp_off_target_chrs_min
         -> FALSE_POSITIVE (chimeric)
      4. construct_frac >= fp_construct_frac_min AND combined_frac >= fp_combined_frac_min
         -> FALSE_POSITIVE (construct+host)
      5. host_fraction >= fp_host_fraction_min AND largest_gap < fp_largest_gap_max
         -> FALSE_POSITIVE (host+gap)
      6. elements present AND host_fraction < cand_host_fraction_max
         -> CANDIDATE
      7. no elements AND host_fraction >= unknown_to_fp_host_fraction_min
         AND construct_frac <= unknown_to_fp_construct_frac_max
         -> FALSE_POSITIVE (host-only)
      8. else -> UNKNOWN
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

    # Rule 6: elements annotated and host fraction acceptable.
    if ev.elements and ev.host_fraction < rules.cand_host_fraction_max:
        return ("CANDIDATE",
                f"{len(ev.elements)} element(s) annotated; host_fraction="
                f"{ev.host_fraction:.0%} below {rules.cand_host_fraction_max:.0%}")

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
