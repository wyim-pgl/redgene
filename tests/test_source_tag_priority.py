"""Test 4-way source tag + element_db-family > univec priority (Task T5).

These tests guard the tag-tiered merge in classify_site_tiers. All three
non-univec sources (element_db, payload, sample_contig) share tier 2 and
collectively beat tier-1 univec at any bitscore tie. Within tier 2, ties
fall back to strict '>' bitscore comparison.
"""
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "s05", REPO / "scripts/s05_insert_assembly.py"
)
s05 = importlib.util.module_from_spec(_spec)
sys.modules["s05"] = s05
_spec.loader.exec_module(s05)

_should_replace = s05._should_replace


def test_element_db_beats_univec_at_tied_bitscore():
    # BUG-3 regression guard - already passing under legacy semantic
    assert _should_replace(
        existing={"src": "univec", "bit": 100},
        new_src="element_db", new_bit=100,
    ) is True


def test_payload_beats_univec():
    assert _should_replace(
        existing={"src": "univec", "bit": 100},
        new_src="payload", new_bit=100,
    ) is True


def test_sample_contig_beats_univec():
    assert _should_replace(
        existing={"src": "univec", "bit": 100},
        new_src="sample_contig", new_bit=100,
    ) is True


def test_element_db_family_tie_keeps_incumbent():
    # element_db / payload / sample_contig share tier 2; tie resolves by strict '>'
    assert _should_replace(
        existing={"src": "element_db", "bit": 150},
        new_src="payload", new_bit=150,
    ) is False


def test_higher_bitscore_wins_within_tier():
    # Within the element_db-family tier, strict '>' bitscore replaces incumbent.
    # Cross-tier (univec -> element_db-family) is already covered above;
    # the legacy BUG-3 guard (test_extra_element_db) also locks down the
    # rule that a lower-tier univec NEVER supplants a tier-2 incumbent,
    # regardless of bitscore.
    assert _should_replace(
        existing={"src": "element_db", "bit": 100},
        new_src="payload", new_bit=200,
    ) is True


def test_site_with_payload_only_hits_is_transgene_positive():
    """Regression guard: after cd-hit clustering, CaMV35S et al may only appear
    under payload tag. classify_site_tiers' has_element_hit must accept tier-2
    sources (element_db/payload/sample_contig), not just literal 'element_db'.

    BUG-3 companion: cluster 0 of gmo_combined_db_v2.fa.clstr shows P-CaMV35S
    payload rep absorbed 3 element_db amplicons. Without this check expansion,
    rice_G281 Chr3:16,439,674 (CaMV35S-containing) would be dropped.
    """
    # Simulate the internal call path: a site whose only tier-2 hits are payload
    hit_5p = {"source": "payload", "bitscore": 150}
    hit_3p = {"source": "payload", "bitscore": 140}

    # Use the same logic as classify_site_tiers line ~1112:
    # has_element_hit = True iff at least one side has a tier-2 source
    TIER2 = {k for k, v in s05._SRC_TIER.items() if v >= 2}
    assert "payload" in TIER2
    assert "element_db" in TIER2
    assert "sample_contig" in TIER2
    assert "univec" not in TIER2

    has_element_hit = (hit_5p.get("source") in TIER2) or \
                      (hit_3p.get("source") in TIER2)
    assert has_element_hit is True

    # Also guard the module-level _TIER2_SRCS frozenset used by the real
    # classify_site_tiers call (single source of truth).
    assert s05._TIER2_SRCS == frozenset(TIER2)


def test_site_with_univec_only_hits_not_transgene_positive():
    """Same check, other direction: univec-only site is NOT transgene-positive."""
    hit_5p = {"source": "univec", "bitscore": 200}
    hit_3p = {"source": "univec", "bitscore": 200}

    TIER2 = {k for k, v in s05._SRC_TIER.items() if v >= 2}

    has_element_hit = (hit_5p.get("source") in TIER2) or \
                      (hit_3p.get("source") in TIER2)
    assert has_element_hit is False
