"""Tests for _apply_canonical_override helper in scripts.s05_insert_assembly.

Issue #3 v1.0 wire-in (scoped): canonical_triplet override after existing
verdict pipeline. Promotes FALSE_POSITIVE/UNKNOWN -> CANDIDATE when a canonical
transgene triplet (e.g., bar + P-CaMV35S + T-ocs) is fully present and
host_fraction is below the CAND threshold.
"""
import sys
from pathlib import Path

# scripts/s05_insert_assembly.py imports at module-scope (pysam, etc.) but
# canonical override helper is pure. Add scripts/ to sys.path and import the
# symbol without triggering BLAST/pysam code paths (only the helper is used).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from s05.verdict import VerdictRules  # noqa: E402
from s05_insert_assembly import _apply_canonical_override  # noqa: E402


_RULES = VerdictRules(
    cand_host_fraction_max=0.80,
    canonical_triplets={
        "default": {"bar", "P-CaMV35S", "T-ocs"},
        "rice_G281": {"hLF1", "P-Gt1", "T-nos"},
    },
)


def test_override_promotes_fp_when_triplet_matches_and_host_low():
    verdict, reason = _apply_canonical_override(
        verdict="FALSE_POSITIVE",
        reason="filter A: 85% host",
        unique_elems={"bar", "P-CaMV35S", "T-ocs", "nptII"},
        host_fraction=0.30,
        rules=_RULES,
    )
    assert verdict == "CANDIDATE"
    assert "canonical_triplet" in reason
    assert "override" in reason.lower()


def test_override_skips_when_host_fraction_too_high():
    verdict, reason = _apply_canonical_override(
        verdict="FALSE_POSITIVE",
        reason="filter A: 90% host",
        unique_elems={"bar", "P-CaMV35S", "T-ocs"},
        host_fraction=0.90,
        rules=_RULES,
    )
    assert verdict == "FALSE_POSITIVE"
    assert reason == "filter A: 90% host"


def test_override_skips_when_triplet_incomplete():
    verdict, reason = _apply_canonical_override(
        verdict="FALSE_POSITIVE",
        reason="filter A",
        unique_elems={"bar", "P-CaMV35S"},  # missing T-ocs
        host_fraction=0.30,
        rules=_RULES,
    )
    assert verdict == "FALSE_POSITIVE"


def test_override_skips_when_rules_is_none():
    verdict, reason = _apply_canonical_override(
        verdict="UNKNOWN",
        reason="",
        unique_elems={"bar", "P-CaMV35S", "T-ocs"},
        host_fraction=0.30,
        rules=None,
    )
    assert verdict == "UNKNOWN"


def test_override_skips_when_unique_elems_empty():
    verdict, reason = _apply_canonical_override(
        verdict="UNKNOWN",
        reason="no elements",
        unique_elems=set(),
        host_fraction=0.30,
        rules=_RULES,
    )
    assert verdict == "UNKNOWN"


def test_override_matches_rice_g281_triplet():
    # hLF1 + P-Gt1 + T-nos — different from default
    verdict, reason = _apply_canonical_override(
        verdict="FALSE_POSITIVE",
        reason="filter B",
        unique_elems={"hLF1", "P-Gt1", "T-nos", "other"},
        host_fraction=0.50,
        rules=_RULES,
    )
    assert verdict == "CANDIDATE"
    assert "canonical_triplet" in reason


def test_override_passthrough_candidate():
    # Already CANDIDATE, triplet also matches — keep CANDIDATE but re-annotate
    verdict, reason = _apply_canonical_override(
        verdict="CANDIDATE",
        reason="elements annotated",
        unique_elems={"bar", "P-CaMV35S", "T-ocs"},
        host_fraction=0.30,
        rules=_RULES,
    )
    assert verdict == "CANDIDATE"
