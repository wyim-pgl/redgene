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
