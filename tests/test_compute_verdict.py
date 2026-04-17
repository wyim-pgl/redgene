"""Tests for scripts.s05.verdict.compute_verdict (pure function).

T6 pytest scenarios — 8 cases covering the decision tree:
CANDIDATE / FALSE_POSITIVE (A/B/C/D) / UNKNOWN reclassification / canonical triplet.
"""
from scripts.s05.verdict import compute_verdict, FilterEvidence, VerdictRules


def _ev(**kw):
    defaults = dict(
        elements=[], host_bp=0, host_fraction=0.0, largest_gap=9999,
        flanking_hit=None, off_target_chrs=[], construct_frac=0.0,
        combined_frac=0.0, is_chimeric=False,
        site_chr="Chr1", site_pos=1000,
        matched_canonical=set(), sources_by_element={},
    )
    defaults.update(kw)
    return FilterEvidence(**defaults)


_RULES = VerdictRules(
    cand_host_fraction_max=0.80,
    cand_largest_gap_min=500,
    fp_host_fraction_min=0.80,
    fp_largest_gap_max=500,
    fp_off_target_chrs_min=2,
    fp_combined_frac_min=0.85,
    fp_construct_frac_min=0.25,
    unknown_to_fp_host_fraction_min=0.85,
    unknown_to_fp_construct_frac_max=0.05,
    canonical_triplets={"default": {"bar", "P-CaMV35S", "T-ocs"}},
    canonical_triplet_min_identity=0.90,
)


def test_candidate_basic():
    ev = _ev(elements=["bar", "P-CaMV35S"], host_fraction=0.3, largest_gap=5000)
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "CANDIDATE"


def test_fp_a_host_fraction_small_gap():
    ev = _ev(host_fraction=0.85, largest_gap=300)
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"
    assert "85" in reason or "0.85" in reason


def test_fp_b_flanking_overlap():
    ev = _ev(flanking_hit=("Chr11", 8758, 8958), site_chr="Chr11", site_pos=8800)
    verdict, _ = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"


def test_fp_c_multi_locus_chimeric():
    ev = _ev(off_target_chrs=[("Chr5", 200), ("Chr7", 150)], is_chimeric=True)
    verdict, _ = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"


def test_fp_d_construct_host_explain():
    ev = _ev(construct_frac=0.35, host_fraction=0.60, combined_frac=0.90)
    verdict, _ = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"


def test_unknown_to_fp_host_only():
    ev = _ev(elements=[], host_fraction=0.90, construct_frac=0.02)
    verdict, _ = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"


def test_unknown_preserved():
    ev = _ev(elements=[], host_fraction=0.50, construct_frac=0.10)
    verdict, _ = compute_verdict(ev, _RULES)
    assert verdict == "UNKNOWN"


def test_canonical_triplet_promotion_from_sample_contig():
    # elements 가 존재하고 canonical triplet 모두 매칭되면 CANDIDATE 승격
    ev = _ev(
        elements=["bar", "P-CaMV35S", "T-ocs"],
        sources_by_element={"bar": "sample_contig",
                            "P-CaMV35S": "sample_contig",
                            "T-ocs": "sample_contig"},
        matched_canonical={"bar", "P-CaMV35S", "T-ocs"},
        host_fraction=0.60,
    )
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "CANDIDATE"
    assert "canonical_triplet" in reason
