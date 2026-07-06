"""Regression tests for pure, tool-independent s05 helpers that previously had
no direct coverage (annotation hit-merging / blast6 parse, primitives I/O)."""
from scripts.s05 import annotation, primitives


# --- primitives ---
def test_revcomp_acgt():
    assert primitives.revcomp("AAACG") == "CGTTT"
    assert primitives.revcomp("ACGT") == "ACGT"  # self-complement


def test_read_write_fasta_roundtrip(tmp_path):
    p = tmp_path / "x.fa"
    primitives.write_fasta(p, "seq1", "ACGTACGTAC", wrap=4)
    assert primitives.read_fasta(p) == {"seq1": "ACGTACGTAC"}


def test_read_fasta_first_token_is_id(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(">seq1 description here\nACGT\n>seq2\nTTTT\n")
    assert primitives.read_fasta(p) == {"seq1": "ACGT", "seq2": "TTTT"}


def test_parse_src_tag():
    assert primitives._parse_src_tag("acc123|src=univec", "def") == ("acc123", "univec")
    assert primitives._parse_src_tag("plain_acc", "element_db") == ("plain_acc", "element_db")


# --- annotation._parse_blast6 ---
def test_parse_blast6_strips_src_tag_and_filters_short(tmp_path):
    p = tmp_path / "b.tsv"
    p.write_text(
        "c1\tbar|src=element_db\t99.0\t120\t1\t120\t1\t120\t1e-50\t200\n"
        "c1\tshort\t99.0\t20\t1\t20\t1\t20\t1e-5\t40\n"          # len 20 < 30 → dropped
        "c1\ttrunc\t99.0\t50\t1\t50\t1\t50\n"                     # <10 cols → dropped
    )
    hits = annotation._parse_blast6(p)
    assert len(hits) == 1
    assert hits[0]["subject"] == "bar"          # |src= stripped
    assert hits[0]["src_tag"] == "element_db"
    assert hits[0]["length"] == 120


def test_parse_blast6_missing_file_returns_empty(tmp_path):
    assert annotation._parse_blast6(tmp_path / "nope.tsv") == []


# --- annotation._merge_annotations ---
def _hit(query, qs, qe, bitscore, subject):
    return {"query": query, "q_start": qs, "q_end": qe,
            "bitscore": bitscore, "subject": subject}


def test_merge_annotations_dominance_and_nonoverlap():
    local = [_hit("c1", 1, 100, 200.0, "barA"),
             _hit("c1", 10, 90, 150.0, "barB")]      # 81 bp, >80% within barA → dominated
    remote = [_hit("c1", 200, 300, 180.0, "ntX")]    # non-overlapping → kept
    merged = annotation._merge_annotations(local, remote)
    subjects = [h["subject"] for h in merged]
    assert "barA" in subjects
    assert "barB" not in subjects
    assert "ntX" in subjects
    assert all("source" in h for h in merged)


def test_merge_annotations_prefers_element_db_on_tie():
    local = [_hit("c1", 1, 100, 150.0, "elemHit")]
    remote = [_hit("c1", 1, 100, 150.0, "ntHit")]
    merged = annotation._merge_annotations(local, remote)
    assert len(merged) == 1
    assert merged[0]["subject"] == "elemHit"
    assert merged[0]["source"] == "element_db"


def test_merge_annotations_empty():
    assert annotation._merge_annotations([], []) == []
