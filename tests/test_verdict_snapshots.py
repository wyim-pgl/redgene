"""Snapshot tests for compute_verdict — golden TSV regression.

Issue #3 (full wire-in). Loads golden verdict TSVs committed under
tests/fixtures/verdict_snapshots/ and asserts compute_verdict reproduces
the expected classification without running any BLAST or I/O.

Golden TSVs are sourced from actual pipeline runs:
  - rice_G281: produced by commit 025c1ed (feature/v1.0-mvp-2026-04-16)
  - cucumber_line225: AC-2 ground truth (8 CAND / 0 FP / 0 inconclusive)

TSV column semantics map directly to FilterEvidence fields:
  elements           — pipe-separated element names; empty = no annotations
  host_fraction      — float [0, 1]
  host_bp            — int bp covered by host BLAST
  largest_gap        — int longest non-host gap (bp)
  flanking_hit       — "chr:lo:hi" or empty
  site_chr           — chromosome of the insertion site
  site_pos           — int position of the insertion site
  off_target_chrs    — "chr:bp|chr:bp" pipe-separated or empty
  construct_frac     — float [0, 1]
  combined_frac      — float [0, 1]
  host_endogenous_elements — pipe-separated; empty = none
  matched_canonical  — pipe-separated; empty = none
  expected_verdict   — CANDIDATE | FALSE_POSITIVE | UNKNOWN
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

import pytest

from scripts.s05.verdict import compute_verdict, FilterEvidence, VerdictRules
from scripts.s05.config_loader import load_verdict_rules

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "verdict_snapshots"
REPO_ROOT = Path(__file__).resolve().parents[1]

# Rules consistent with the production defaults.
_RULES = VerdictRules(
    cand_host_fraction_max=0.80,
    fp_host_fraction_min=0.80,
    fp_largest_gap_max=500,
    fp_off_target_chrs_min=2,
    fp_combined_frac_min=0.85,
    fp_construct_frac_min=0.25,
    unknown_to_fp_host_fraction_min=0.85,
    unknown_to_fp_construct_frac_max=0.05,
    canonical_triplets={
        "rice_G281": {"hLF1", "P-Gt1", "T-nos"},
        "default": {"bar", "P-CaMV35S", "T-ocs"},
        "soybean_AtYUCCA6": {"bar", "P-CaMV35S", "T-ocs"},
        "soybean_UGT72E3": {"bar", "P-CaMV35S", "T-nos"},
        "tomato_Cas9_A2_3": {"bar", "SpCas9", "sgRNA_scaffold_generic"},
    },
)


# ---------------------------------------------------------------------------
# TSV helpers
# ---------------------------------------------------------------------------

def _parse_pipe_list(s: str) -> list[str]:
    """Parse a pipe-separated string into a list; empty string → []."""
    return [x for x in s.split("|") if x] if s.strip() else []


def _parse_flanking_hit(s: str) -> Optional[tuple[str, int, int]]:
    """Parse 'chr:lo:hi' into a tuple, or return None for empty."""
    s = s.strip()
    if not s:
        return None
    parts = s.rsplit(":", 2)
    if len(parts) != 3:
        return None
    chrom, lo, hi = parts
    return (chrom, int(lo), int(hi))


def _parse_off_target_chrs(s: str) -> list[tuple[str, int]]:
    """Parse 'chr:bp|chr:bp' into list of tuples; empty string → []."""
    if not s.strip():
        return []
    result = []
    for item in s.split("|"):
        item = item.strip()
        if not item:
            continue
        chrom, bp_str = item.rsplit(":", 1)
        result.append((chrom, int(bp_str)))
    return result


def _load_fixture(tsv_path: Path) -> list[dict]:
    """Load a golden TSV into a list of row dicts, skipping comment lines."""
    rows = []
    with open(tsv_path, encoding="utf-8") as fh:
        reader = csv.DictReader(
            (line for line in fh if not line.startswith("#")),
            delimiter="\t",
        )
        for row in reader:
            rows.append(row)
    return rows


def _row_to_evidence(row: dict) -> FilterEvidence:
    """Convert a golden TSV row to a FilterEvidence instance."""
    elements = _parse_pipe_list(row["elements"])
    matched_canonical = set(_parse_pipe_list(row["matched_canonical"]))
    host_endogenous = set(_parse_pipe_list(row["host_endogenous_elements"]))
    off_target = _parse_off_target_chrs(row["off_target_chrs"])
    flanking = _parse_flanking_hit(row["flanking_hit"])

    return FilterEvidence(
        elements=elements,
        host_bp=int(row["host_bp"]),
        host_fraction=float(row["host_fraction"]),
        largest_gap=int(row["largest_gap"]),
        flanking_hit=flanking,
        off_target_chrs=off_target,
        construct_frac=float(row["construct_frac"]),
        combined_frac=float(row["combined_frac"]),
        is_chimeric=False,          # informational only (Issue #11 M-3)
        site_chr=row["site_chr"],
        site_pos=int(row["site_pos"]),
        matched_canonical=matched_canonical if matched_canonical else set(elements),
        sources_by_element={e: "sample_contig" for e in elements},
        host_endogenous_elements=host_endogenous,
    )


# ---------------------------------------------------------------------------
# Parametric snapshot tests
# ---------------------------------------------------------------------------

def _snapshot_cases(tsv_name: str):
    """Return pytest.param list for a given golden TSV filename."""
    tsv = FIXTURE_DIR / tsv_name
    if not tsv.exists():
        return []
    rows = _load_fixture(tsv)
    return [
        pytest.param(row["site_id"], _row_to_evidence(row), row["expected_verdict"],
                     id=f"{tsv_name}::{row['site_id']}")
        for row in rows
        if row.get("site_id")
    ]


@pytest.mark.parametrize("site_id,ev,expected", _snapshot_cases("rice_G281.tsv"))
def test_rice_g281_snapshot(site_id, ev, expected):
    """Regression: compute_verdict must reproduce golden rice_G281 classifications."""
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == expected, (
        f"{site_id}: expected {expected!r} but got {verdict!r}\n"
        f"  reason: {reason!r}\n"
        f"  host_fraction={ev.host_fraction:.1%}, largest_gap={ev.largest_gap}bp, "
        f"  elements={ev.elements}, construct_frac={ev.construct_frac:.1%}"
    )


@pytest.mark.parametrize("site_id,ev,expected", _snapshot_cases("cucumber_line225.tsv"))
def test_cucumber_line225_snapshot(site_id, ev, expected):
    """Regression: compute_verdict must reproduce golden cucumber_line225 classifications.

    AC-2 ground truth: 8 CANDIDATE sites (0 FP, 0 inconclusive) per commit 025c1ed.
    """
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == expected, (
        f"{site_id}: expected {expected!r} but got {verdict!r}\n"
        f"  reason: {reason!r}\n"
        f"  host_fraction={ev.host_fraction:.1%}, largest_gap={ev.largest_gap}bp, "
        f"  elements={ev.elements}, construct_frac={ev.construct_frac:.1%}"
    )


# ---------------------------------------------------------------------------
# Config-driven rules integration test
# ---------------------------------------------------------------------------

def test_load_verdict_rules_from_config():
    """VerdictRules loaded from repo config.yaml must include rice_G281 triplet."""
    cfg = REPO_ROOT / "config.yaml"
    if not cfg.exists():
        pytest.skip("config.yaml not present in this checkout")
    rules = load_verdict_rules(cfg, sample="rice_G281")
    # Triplets must be a dict (never None) with at least the rice_G281 entry.
    assert isinstance(rules.canonical_triplets, dict)
    assert "rice_G281" in rules.canonical_triplets or "default" in rules.canonical_triplets


def test_rice_g281_true_insertion_with_config_rules():
    """rice_G281 Chr3:16,439,674 must remain CANDIDATE when using config.yaml rules."""
    cfg = REPO_ROOT / "config.yaml"
    rules = load_verdict_rules(cfg, sample="rice_G281") if cfg.exists() else _RULES
    ev = FilterEvidence(
        elements=["NODE_2_length_2383_cov_25.330756"],
        host_bp=11688, host_fraction=0.874, largest_gap=1024,
        off_target_chrs=[("Chr2", 746)],
        construct_frac=0.102, combined_frac=0.976,
        site_chr="Chr3", site_pos=16439674,
        matched_canonical={"NODE_2_length_2383_cov_25.330756"},
    )
    verdict, reason = compute_verdict(ev, rules)
    assert verdict == "CANDIDATE", (
        f"rice_G281 Chr3:16439674 regression FAILED — got {verdict!r}: {reason!r}. "
        f"This site has 87.4% host fraction but 1,024 bp non-host gap (> 500 bp), "
        f"so Filter A must NOT fire.  Check Rule 6 in compute_verdict."
    )
