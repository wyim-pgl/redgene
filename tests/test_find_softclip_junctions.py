"""End-to-end regression for `find_softclip_junctions` on a synthetic BAM.

The helper-level tests in `test_site_discovery.py` pin the consensus, the
clusterer and the pairing in isolation.  Nothing pinned the function that wires
them together: `compute_verdict`'s snapshot fixtures sit *downstream* of site
discovery, so a change to which sites Phase 1 emits produces no snapshot drift.

This builds a 6 kb host, plants a 400 bp foreign insert at one locus, and writes
a coordinate-sorted BAM whose reads carry the soft clips a real bidirectional
junction produces.  It then asserts the discovered site, its seeds, and the
read-quality and depth gates.
"""
from __future__ import annotations

import random
import shutil

import pysam
import pytest

from scripts.s05.site_discovery import find_softclip_junctions


HOST_LEN = 6000
JUNCTION = 3000
ANCHOR = 40          # bp of host alignment per read
CLIP = 30            # bp of insert carried as a soft clip


def _dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture(scope="module")
def genome():
    host = _dna(HOST_LEN, seed=11)
    insert = _dna(400, seed=99)
    return host, insert


def _write_bam(path, host: str, reads: list[tuple[int, str, str]]) -> None:
    """reads: [(reference_start, cigarstring, sequence, ...extra flags via dict)]"""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chrT", "LN": len(host)}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (start, cigar, seq, flags) in enumerate(reads):
            a = pysam.AlignedSegment(out.header)
            a.query_name = f"read{i}"
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.reference_id = 0
            a.reference_start = start
            a.cigarstring = cigar
            a.mapping_quality = flags.get("mapq", 60)
            a.is_duplicate = flags.get("duplicate", False)
            a.is_qcfail = flags.get("qcfail", False)
            out.write(a)
    pysam.index(str(path))


def _junction_reads(host: str, insert: str, depth: int, **flags):
    """`depth` right-clip reads and `depth` left-clip reads at JUNCTION."""
    reads = []
    for _ in range(depth):
        # 5' side: host ends at JUNCTION, then the insert starts -> right clip
        seq = host[JUNCTION - ANCHOR:JUNCTION] + insert[:CLIP]
        reads.append((JUNCTION - ANCHOR, f"{ANCHOR}M{CLIP}S", seq, flags))
    for _ in range(depth):
        # 3' side: the insert ends, then host resumes at JUNCTION -> left clip
        seq = insert[-CLIP:] + host[JUNCTION:JUNCTION + ANCHOR]
        reads.append((JUNCTION, f"{CLIP}S{ANCHOR}M", seq, flags))
    return sorted(reads, key=lambda r: r[0])


@pytest.fixture()
def workspace(tmp_path, genome):
    host, insert = genome
    host_ref = tmp_path / "host.fa"
    host_ref.write_text(">chrT\n" + host + "\n")
    return tmp_path, host_ref, host, insert


pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None, reason="minimap2 not on PATH"
)


def test_bidirectional_junction_is_discovered(workspace):
    tmp_path, host_ref, host, insert = workspace
    bam = tmp_path / "host.bam"
    _write_bam(bam, host, _junction_reads(host, insert, depth=3))

    sites = find_softclip_junctions(bam, host_ref, None, tmp_path / "work",
                                    min_clip=20)

    assert len(sites) == 1
    site = sites[0]
    assert site.host_chr == "chrT"
    assert site.pos_5p == JUNCTION
    assert site.pos_3p == JUNCTION
    assert site.is_validated is True
    assert site.confidence == "high"


def test_discovered_seeds_are_the_insert_ends(workspace):
    tmp_path, host_ref, host, insert = workspace
    bam = tmp_path / "host.bam"
    _write_bam(bam, host, _junction_reads(host, insert, depth=3))

    site = find_softclip_junctions(bam, host_ref, None, tmp_path / "work",
                                   min_clip=20)[0]

    assert site.seed_5p == insert[:CLIP]
    assert site.seed_3p == insert[-CLIP:]
    assert "N" not in site.seed_5p and "N" not in site.seed_3p


def test_cluster_below_min_depth_is_not_reported(workspace):
    tmp_path, host_ref, host, insert = workspace
    bam = tmp_path / "host.bam"
    _write_bam(bam, host, _junction_reads(host, insert, depth=2))

    sites = find_softclip_junctions(bam, host_ref, None, tmp_path / "work",
                                    min_clip=20, min_cluster_depth=3)

    assert sites == []


def test_min_cluster_depth_can_be_lowered(workspace):
    """The knob the coverage-sensitivity work needs."""
    tmp_path, host_ref, host, insert = workspace
    bam = tmp_path / "host.bam"
    _write_bam(bam, host, _junction_reads(host, insert, depth=2))

    sites = find_softclip_junctions(bam, host_ref, None, tmp_path / "work",
                                    min_clip=20, min_cluster_depth=2)

    assert len(sites) == 1
    assert sites[0].pos_5p == JUNCTION


def test_duplicate_reads_do_not_support_a_junction(workspace):
    """3 duplicate-flagged pairs must not clear min_cluster_depth=3."""
    tmp_path, host_ref, host, insert = workspace
    bam = tmp_path / "host.bam"
    _write_bam(bam, host, _junction_reads(host, insert, depth=3, duplicate=True))

    sites = find_softclip_junctions(bam, host_ref, None, tmp_path / "work",
                                    min_clip=20)

    assert sites == []


def test_qcfail_reads_do_not_support_a_junction(workspace):
    tmp_path, host_ref, host, insert = workspace
    bam = tmp_path / "host.bam"
    _write_bam(bam, host, _junction_reads(host, insert, depth=3, qcfail=True))

    sites = find_softclip_junctions(bam, host_ref, None, tmp_path / "work",
                                    min_clip=20)

    assert sites == []


def test_min_mapq_filters_low_quality_anchors(workspace):
    tmp_path, host_ref, host, insert = workspace
    bam = tmp_path / "host.bam"
    _write_bam(bam, host, _junction_reads(host, insert, depth=3, mapq=5))

    assert find_softclip_junctions(bam, host_ref, None, tmp_path / "work",
                                   min_clip=20, min_mapq=20) == []
    # ...and the same BAM still yields the site at the default floor.
    assert len(find_softclip_junctions(bam, host_ref, None, tmp_path / "work2",
                                       min_clip=20, min_mapq=0)) == 1


def _two_clip_reads(host, clip_5p: str, clip_3p: str, depth: int = 3):
    reads = []
    for _ in range(depth):
        reads.append((JUNCTION - ANCHOR, f"{ANCHOR}M{len(clip_5p)}S",
                      host[JUNCTION - ANCHOR:JUNCTION] + clip_5p, {}))
    for _ in range(depth):
        reads.append((JUNCTION, f"{len(clip_3p)}S{ANCHOR}M",
                      clip_3p + host[JUNCTION:JUNCTION + ANCHOR], {}))
    return sorted(reads, key=lambda r: r[0])


def test_identical_clips_are_rejected_as_not_an_insert(workspace):
    """Condition 2: both sides showing the same sequence is an SV, not an insert."""
    tmp_path, host_ref, host, insert = workspace
    bam = tmp_path / "host.bam"
    _write_bam(bam, host, _two_clip_reads(host, insert[:60], insert[:60]))

    assert find_softclip_junctions(bam, host_ref, None, tmp_path / "work",
                                   min_clip=20) == []


def test_clips_that_map_back_to_host_are_rejected(workspace):
    """Condition 3: host-derived clips are a rearrangement, not foreign DNA.

    The two clips are taken from *different* distant host loci so that condition
    2 (clips must differ) passes and the host-mapping check is what rejects the
    site.  The companion assertion swaps in foreign sequence to prove the same
    read layout is otherwise discovered — i.e. this test fails for the right
    reason.
    """
    tmp_path, host_ref, host, insert = workspace

    host_bam = tmp_path / "host_derived.bam"
    _write_bam(host_bam, host, _two_clip_reads(host, host[5000:5200], host[800:1000]))
    assert find_softclip_junctions(host_bam, host_ref, None, tmp_path / "w1",
                                   min_clip=20) == []

    foreign_bam = tmp_path / "foreign.bam"
    _write_bam(foreign_bam, host, _two_clip_reads(host, insert[:200], insert[200:400]))
    assert len(find_softclip_junctions(foreign_bam, host_ref, None, tmp_path / "w2",
                                       min_clip=20)) == 1
