"""TDD tests for Issue #3 — compute_verdict filter priority reconciliation.

Priority order (canonical, first-match wins):
  1. canonical_triplet (promotion, runs first)
  1.5. host_endogenous (all-elements-are-host → FP)
  2. Filter B: construct-flanking overlap → FP
  3. Filter C: multi-locus chimeric → FP
  4. Filter D: construct+host coverage explains insert → FP
  5. Filter A: host_fraction+gap → FP
  6. elements present (any FP filter survived) → CANDIDATE
  7. UNKNOWN → FP reclassification (host-only, no elements)
  8. fallthrough → UNKNOWN

Key design decision: Rule 6 must NOT gate on host_fraction < cand_host_fraction_max.
A site that has elements and survives all FP filters (1-5) is CANDIDATE regardless
of host_fraction, because:
- Rule 5 (Filter A) already catches the "high host_fraction + small gap" FP case.
- When host_fraction is high BUT the gap is large (≥fp_largest_gap_max), the assembly
  spans genuine T-DNA (host→T-DNA→host junction contig). This is a true CANDIDATE.
- The rice_G281 Chr3:16,439,674 site (87.4% host, 1,024bp gap) is the canonical example:
  Filter A does not fire (gap ≥ 500), so the site MUST be CANDIDATE, not UNKNOWN.
"""
from __future__ import annotations

import pytest

from scripts.s05.verdict import compute_verdict, FilterEvidence, VerdictRules


_RULES = VerdictRules(
    cand_host_fraction_max=0.80,
    fp_host_fraction_min=0.80,
    fp_largest_gap_max=500,
    fp_off_target_chrs_min=2,
    fp_combined_frac_min=0.85,
    fp_construct_frac_min=0.25,
    unknown_to_fp_host_fraction_min=0.85,
    unknown_to_fp_construct_frac_max=0.05,
    canonical_triplets={"rice_G281": {"hLF1", "P-Gt1", "T-nos"},
                        "default": {"bar", "P-CaMV35S", "T-ocs"}},
)


# ---------------------------------------------------------------------------
# Priority rule P-1: Rule 6 must fire even when host_fraction >= threshold,
# as long as no FP filter (A–D) has fired first.  This is the rice_G281
# Chr3:16,439,674 regression case (87.4% host, 1,024 bp gap, elements present).
# ---------------------------------------------------------------------------

def test_p1_rule6_candidate_high_host_large_gap():
    """Elements + high host_fraction + large gap = CANDIDATE (rice_G281 case).

    Filter A does not fire (gap ≥ fp_largest_gap_max).  No other FP filter
    fires.  Rule 6 must return CANDIDATE.
    """
    ev = FilterEvidence(
        elements=["NODE_2_length_2383_cov_25"],
        host_bp=11_688, host_fraction=0.874, largest_gap=1_024,
        flanking_hit=None, off_target_chrs=[("Chr2", 746)],
        construct_frac=0.102, combined_frac=0.976,
        site_chr="Chr3", site_pos=16_439_674,
    )
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "CANDIDATE", (
        f"Expected CANDIDATE for rice_G281 canonical case; got {verdict!r} ({reason!r}). "
        f"Rule 6 must not gate on host_fraction < cand_host_fraction_max when "
        f"Filter A already handles the high-host + small-gap FP case."
    )


def test_p1_rule6_candidate_low_host():
    """Standard case: elements + low host_fraction = CANDIDATE."""
    ev = FilterEvidence(
        elements=["bar", "P-CaMV35S"],
        host_fraction=0.30, largest_gap=5_000,
        construct_frac=0.10, combined_frac=0.40,
    )
    verdict, _ = compute_verdict(ev, _RULES)
    assert verdict == "CANDIDATE"


def test_p1_rule6_candidate_exactly_at_threshold():
    """Edge: host_fraction == fp_host_fraction_min (0.80) with large gap = CANDIDATE."""
    ev = FilterEvidence(
        elements=["bar"],
        host_fraction=0.80, largest_gap=501,   # gap just barely ≥ fp_largest_gap_max
        construct_frac=0.05, combined_frac=0.85,
    )
    # combined_frac == 0.85 and construct_frac < 0.25 → Filter D does NOT fire.
    # Filter A: host_fraction==0.80 AND gap 501 ≥ 500 → does NOT fire.
    # Rule 6 must fire.
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "CANDIDATE", f"Got {verdict!r}: {reason!r}"


# ---------------------------------------------------------------------------
# Priority rule P-2: Filter A (Rule 5) must come AFTER Filters B, C, D
# in compute_verdict, and BEFORE Rule 6.  Verify through specific ordering
# scenarios where only one filter fires at a time.
# ---------------------------------------------------------------------------

def test_p2_filter_b_before_a():
    """Filter B (flanking) fires; Filter A data present but irrelevant."""
    ev = FilterEvidence(
        elements=["bar"],
        host_fraction=0.50, largest_gap=200,  # Filter A would fire if checked
        flanking_hit=("Chr1", 1_000, 2_000),
        site_chr="Chr1", site_pos=1_500,
    )
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"
    assert "flanking" in reason.lower() or "overlaps" in reason.lower()


def test_p2_filter_c_chimeric_fires():
    """Filter C (chimeric) fires with ≥2 off-target chromosomes."""
    ev = FilterEvidence(
        elements=["bar"],
        host_fraction=0.30, largest_gap=5_000,
        off_target_chrs=[("Chr5", 300), ("Chr7", 200)],
    )
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"
    assert "chimeric" in reason.lower()


def test_p2_filter_d_construct_host_fires():
    """Filter D fires when construct+host fully explains insert."""
    ev = FilterEvidence(
        elements=["bar"],
        host_fraction=0.60, largest_gap=5_000,
        construct_frac=0.30, combined_frac=0.90,
    )
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"
    assert "construct" in reason.lower() or "combined" in reason.lower()


def test_p2_filter_a_fires_small_gap():
    """Filter A fires: high host_fraction + small gap."""
    ev = FilterEvidence(
        elements=["bar"],
        host_fraction=0.85, largest_gap=300,  # gap < 500
        construct_frac=0.05, combined_frac=0.90,
    )
    # combined_frac=0.90 ≥ 0.85, but construct_frac=0.05 < 0.25 → Filter D does NOT fire.
    # Filter A: 0.85 ≥ 0.80 AND 300 < 500 → fires.
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"
    assert "host" in reason.lower()


# ---------------------------------------------------------------------------
# Priority rule P-3: canonical_triplet beats ALL FP filters (highest priority).
# ---------------------------------------------------------------------------

def test_p3_canonical_triplet_beats_filter_a():
    """Canonical triplet wins even when Filter A would fire."""
    ev = FilterEvidence(
        elements=["bar", "P-CaMV35S", "T-ocs"],
        matched_canonical={"bar", "P-CaMV35S", "T-ocs"},
        host_fraction=0.60,   # below cand_host_fraction_max
        largest_gap=200,      # Filter A would fire without canonical
        construct_frac=0.05, combined_frac=0.65,
    )
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "CANDIDATE"
    assert "canonical_triplet" in reason


def test_p3_canonical_triplet_does_not_beat_high_host():
    """Canonical triplet does NOT promote when host_fraction >= cand_host_fraction_max."""
    ev = FilterEvidence(
        elements=["bar", "P-CaMV35S", "T-ocs"],
        matched_canonical={"bar", "P-CaMV35S", "T-ocs"},
        host_fraction=0.90,   # above cand_host_fraction_max → triplet gate fails
        largest_gap=200,
        construct_frac=0.05, combined_frac=0.95,
    )
    verdict, _ = compute_verdict(ev, _RULES)
    # canonical_triplet rule (1) fails due to high host_fraction;
    # Filter A fires (0.90 ≥ 0.80, gap 200 < 500).
    assert verdict == "FALSE_POSITIVE"


# ---------------------------------------------------------------------------
# Priority rule P-4: host_endogenous (Rule 1.5) beats Filters B/C/D/A/6
# but NOT canonical_triplet (Rule 1).
# ---------------------------------------------------------------------------

def test_p4_host_endogenous_beats_rule6():
    """All elements endogenous → FP even though elements are present."""
    ev = FilterEvidence(
        elements=["OsActin1", "OsUbi1"],
        host_endogenous_elements={"OsActin1", "OsUbi1"},
        host_fraction=0.30, largest_gap=5_000,
    )
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"
    assert "host-endogenous" in reason


def test_p4_canonical_triplet_beats_host_endogenous():
    """Canonical triplet (Rule 1) must take priority over host_endogenous (Rule 1.5)."""
    ev = FilterEvidence(
        elements=["bar", "P-CaMV35S", "T-ocs"],
        matched_canonical={"bar", "P-CaMV35S", "T-ocs"},
        host_endogenous_elements={"bar", "P-CaMV35S", "T-ocs"},
        host_fraction=0.30, largest_gap=5_000,
    )
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "CANDIDATE"
    assert "canonical_triplet" in reason
