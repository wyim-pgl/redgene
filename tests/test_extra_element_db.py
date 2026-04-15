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
        "hit_common":     COMMON_SEQ,
        "hit_sample":     SAMPLE_SEQ,
        "unrelated_clip": UNRELATED,
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
    # Negative: unrelated clip matches none of the three DBs
    assert "unrelated_clip" not in hits or hits["unrelated_clip"] == [], \
        f"unrelated_clip should not hit any DB, got: {hits.get('unrelated_clip')}"


def test_annotate_insert_uses_extra_dbs(tmp_path):
    """Phase 3 annotate_insert must consult extra_dbs too, else s04b
    contigs / common_payload payloads (bar, AtYUCCA6, etc.) silently
    drop out of element_annotation.tsv and every site lands on the
    'no element annotations' UNKNOWN branch."""
    # Minimal insert containing two distinct 50bp subsequences, chosen
    # to survive BLAST's default dust filter.
    insert_fa = tmp_path / "insert.fasta"
    primary_seq = "ATCGATCGATCGAAGCTTGGATCCAAGCTAGCTAGCTAGAACCGGTTAACC"
    extra_seq   = "TGCAAGTTCGATCGTACGTAGCTAGCATCGATCGTTAAGGCCTAGCTTGCA"
    insert_fa.write_text(">insert1\n" + primary_seq + extra_seq + "\n")

    primary_db = tmp_path / "primary.fa"
    primary_db.write_text(">elem_A\n" + primary_seq + "\n")

    extra_db = tmp_path / "extra.fa"
    extra_db.write_text(">bar|X17220.1\n" + extra_seq + "\n")

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    s05.annotate_insert(
        insert_fa, primary_db, out_dir,
        sample_name="t",
        no_remote_blast=True,
        extra_dbs=[extra_db],
    )
    ann = (out_dir / "element_annotation.tsv").read_text()
    assert "elem_A" in ann, f"primary-DB hit missing from annotation: {ann}"
    assert "bar" in ann, f"extra-DB hit missing from annotation: {ann}"
