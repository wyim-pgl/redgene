"""Unit tests for Phase 1 site discovery — the soft-clip junction detector.

Before this file `scripts/s05/site_discovery.py` had no direct coverage at all:
`_build_consensus`, the position clusterer, and the right/left cluster pairing
were exercised only through full-pipeline runs that need a host BAM, minimap2
and a reference genome.

Covers:
  * A4 — the left-clip consensus must not leak internal ``N`` into a seed.
  * B1 — duplicate / QC-fail / low-MAPQ reads must not donate soft clips.
  * B2 — cluster pairing must be deterministic and distance-optimal.
  * B3 — cluster depth thresholds must be callable parameters.
"""
from __future__ import annotations

import pysam
import pytest

from scripts.s05.site_discovery import (
    _build_consensus,
    _cluster_positions,
    _pair_clusters,
    _read_donates_clip,
)


# ---------------------------------------------------------------------------
# A4 — consensus must never hand an internal N to the assembler
# ---------------------------------------------------------------------------

def test_left_consensus_drops_internal_n():
    """Two equal-depth clips disagreeing at one internal base.

    The tie makes that column ambiguous.  ``seed_3p`` feeds
    StrandAwareSeedExtender, which skips any sequence containing ``N``, so the
    consensus must stop at the ambiguity rather than embed it.  Left clips abut
    the host on their RIGHT end, so the informative part is the suffix.
    """
    consensus = _build_consensus(["AAACAAA", "AAAGAAA"], "left")

    assert "N" not in consensus
    assert consensus == "AAA"


def test_left_consensus_keeps_full_agreement():
    consensus = _build_consensus(["ACGTACG", "ACGTACG", "ACGTACG"], "left")

    assert consensus == "ACGTACG"


def test_left_consensus_trims_shallow_leading_bases():
    """Positions supported by a single clip are not consensus — they are trimmed."""
    consensus = _build_consensus(["TTTTGGGG", "GGGG"], "left")

    assert consensus == "GGGG"


def test_right_consensus_stops_at_first_ambiguity():
    """The right branch already truncates; assert the contract explicitly."""
    consensus = _build_consensus(["ACGTAAA", "ACGTCCC"], "right")

    assert consensus == "ACGT"


def test_empty_clip_stack_yields_empty_consensus():
    assert _build_consensus([], "left") == ""
    assert _build_consensus([], "right") == ""


# ---------------------------------------------------------------------------
# B3 — clustering with an explicit depth threshold
# ---------------------------------------------------------------------------

def test_cluster_positions_drops_shallow_clusters():
    clips = [(100, "AAA"), (101, "AAA"), (500, "CCC")]

    clusters = _cluster_positions(clips, window=50, min_depth=2)

    assert len(clusters) == 1
    assert clusters[0][0] == 101  # median of [100, 101]


def test_cluster_positions_min_depth_one_keeps_singletons():
    clips = [(100, "AAA"), (500, "CCC")]

    clusters = _cluster_positions(clips, window=50, min_depth=1)

    assert [pos for pos, _ in clusters] == [100, 500]


def test_cluster_positions_splits_beyond_window():
    clips = [(100, "A"), (110, "A"), (400, "C"), (405, "C")]

    clusters = _cluster_positions(clips, window=50, min_depth=2)

    assert [pos for pos, _ in clusters] == [110, 405]


# ---------------------------------------------------------------------------
# B2 — pairing must be distance-optimal, not iteration-order-dependent
# ---------------------------------------------------------------------------

def test_pairing_prefers_the_closer_right_cluster():
    """Regression for the order-dependent greedy.

    Walking right clusters in positional order lets r=100 claim l=108 (dist 8)
    and strands r=110 (dist 2).  Global distance ordering must pair 110 <-> 108.
    """
    pairs = _pair_clusters(r_positions=[100, 110], l_positions=[108], window=50)

    assert pairs == [(1, 0)]


def test_pairing_is_order_independent():
    forward = _pair_clusters(r_positions=[100, 110], l_positions=[108], window=50)
    # Same clusters, right side enumerated the other way around.
    reverse = _pair_clusters(r_positions=[110, 100], l_positions=[108], window=50)

    assert forward == [(1, 0)]
    assert reverse == [(0, 0)]


def test_pairing_respects_the_window():
    assert _pair_clusters(r_positions=[100], l_positions=[400], window=50) == []


def test_pairing_matches_each_cluster_at_most_once():
    pairs = _pair_clusters(r_positions=[100, 101], l_positions=[100, 102], window=50)

    r_used = [r for r, _ in pairs]
    l_used = [l for _, l in pairs]
    assert len(pairs) == 2
    assert sorted(r_used) == [0, 1]
    assert sorted(l_used) == [0, 1]


# ---------------------------------------------------------------------------
# B1 — which reads are allowed to donate a soft clip
# ---------------------------------------------------------------------------

def _make_read(*, mapq: int = 60, duplicate: bool = False, qcfail: bool = False,
               secondary: bool = False, supplementary: bool = False,
               unmapped: bool = False) -> pysam.AlignedSegment:
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100_000}]}
    )
    read = pysam.AlignedSegment(header)
    read.query_name = "r1"
    read.query_sequence = "ACGT" * 15
    read.query_qualities = pysam.qualitystring_to_array("I" * 60)
    read.reference_id = 0
    read.reference_start = 1000
    read.cigarstring = "40M20S"
    read.mapping_quality = mapq
    read.is_unmapped = unmapped
    read.is_duplicate = duplicate
    read.is_qcfail = qcfail
    read.is_secondary = secondary
    read.is_supplementary = supplementary
    return read


def test_ordinary_read_donates_a_clip():
    assert _read_donates_clip(_make_read(), min_mapq=0) is True


def test_pcr_duplicate_is_rejected():
    assert _read_donates_clip(_make_read(duplicate=True), min_mapq=0) is False


def test_qcfail_read_is_rejected():
    assert _read_donates_clip(_make_read(qcfail=True), min_mapq=0) is False


def test_secondary_and_supplementary_are_rejected():
    assert _read_donates_clip(_make_read(secondary=True), min_mapq=0) is False
    assert _read_donates_clip(_make_read(supplementary=True), min_mapq=0) is False


def test_unmapped_read_is_rejected():
    assert _read_donates_clip(_make_read(unmapped=True), min_mapq=0) is False


@pytest.mark.parametrize("mapq,min_mapq,expected", [
    (0, 0, True),     # default threshold keeps today's behaviour
    (0, 20, False),
    (19, 20, False),
    (20, 20, True),
    (60, 20, True),
])
def test_mapq_threshold_is_enforced(mapq, min_mapq, expected):
    assert _read_donates_clip(_make_read(mapq=mapq), min_mapq=min_mapq) is expected
