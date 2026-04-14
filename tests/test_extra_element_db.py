"""Verify s05 _batch_check_element_hits can consult extra FASTAs."""
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "s05", REPO / "scripts/s05_insert_assembly.py"
)
s05 = importlib.util.module_from_spec(spec)
sys.modules["s05"] = s05
spec.loader.exec_module(s05)


# Three distinct high-complexity 50bp sequences — chosen so BLAST's default
# dust filter can't mask them. Kept at module scope so tests stay in sync.
PRIMARY_SEQ = "ATCGATCGATCGAAGCTTGGATCCAAGCTAGCTAGCTAGAACCGGTTAACC"
COMMON_SEQ  = "GTACCGGTACCGGTTAACGCTAGCATGCATGCATGCAAGCTTGGATCGATC"
SAMPLE_SEQ  = "TGCAAGTTCGATCGTACGTAGCTAGCATCGATCGTTAAGGCCTAGCTTGCA"
UNRELATED   = "CACACACAGTGTGTGTTTACACACACAGTGTGTGTTTACACACACAGTGTA"


def test_batch_check_element_hits_reads_extra_db(tmp_path):
    """Single extra DB via the new list-based signature."""
    primary = tmp_path / "primary.fa"
    primary.write_text(">elem_A\n" + PRIMARY_SEQ + "\n")
    extra = tmp_path / "extra.fa"
    extra.write_text(">bar|X17220.1\n" + COMMON_SEQ + "\n")
    # Query sequence matches extra, not primary
    seqs = {"clip_q1": COMMON_SEQ}
    hits = s05._batch_check_element_hits(
        seqs, primary, tmp_path, extra_dbs=[extra],
    )
    assert "clip_q1" in hits, f"no hits found: {hits}"
    hit_str = " ".join(hits["clip_q1"])
    assert "bar" in hit_str, f"expected 'bar' in hits, got: {hit_str}"


def test_batch_check_element_hits_without_extras_ignores_extra_seq(tmp_path):
    """When extra_dbs is None, a query matching only a would-be extra seq
    must NOT produce hits. Guards against accidental extra-DB consultation."""
    primary = tmp_path / "primary.fa"
    primary.write_text(">elem_A\n" + PRIMARY_SEQ + "\n")

    seqs = {"clip_only_extra": COMMON_SEQ}
    hits = s05._batch_check_element_hits(seqs, primary, tmp_path, extra_dbs=None)
    assert "clip_only_extra" not in hits or hits["clip_only_extra"] == []


def test_batch_check_element_hits_consults_multiple_dbs(tmp_path):
    """Two extra DBs (shared common_payload + per-sample contigs): hits
    from both must merge into the result."""
    primary = tmp_path / "primary.fa"
    primary.write_text(">elem_A\n" + PRIMARY_SEQ + "\n")
    common = tmp_path / "common.fa"
    common.write_text(">bar|X17220.1\n" + COMMON_SEQ + "\n")
    sample_specific = tmp_path / "sample.fa"
    sample_specific.write_text(">node_1\n" + SAMPLE_SEQ + "\n")

    seqs = {
        "hit_common":  COMMON_SEQ,
        "hit_sample":  SAMPLE_SEQ,
    }
    hits = s05._batch_check_element_hits(
        seqs, primary, tmp_path,
        extra_dbs=[common, sample_specific],
    )
    assert "hit_common" in hits, f"common-payload hit missing: {hits}"
    assert any("bar" in h for h in hits["hit_common"]), \
        f"expected 'bar' in common hits, got: {hits['hit_common']}"
    assert "hit_sample" in hits, f"sample-specific hit missing: {hits}"
    assert any("node_1" in h for h in hits["hit_sample"]), \
        f"expected 'node_1' in sample hits, got: {hits['hit_sample']}"
