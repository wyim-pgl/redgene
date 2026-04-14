"""Verify s05 _batch_check_element_hits can consult a second FASTA."""
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


def test_batch_check_element_hits_reads_extra_db(tmp_path):
    # Use varied sequences; BLAST filters low-complexity homopolymers
    primary_seq = (
        "ATCGATCGATCGAAGCTTGGATCCAAGCTAGCTAGCTAGAACCGGTTAACC"  # 50bp
    )
    extra_seq = (
        "GTACCGGTACCGGTTAACGCTAGCATGCATGCATGCAAGCTTGGATCGATC"  # 50bp, distinct
    )
    primary = tmp_path / "primary.fa"
    primary.write_text(">elem_A\n" + primary_seq + "\n")
    extra = tmp_path / "extra.fa"
    extra.write_text(">bar|X17220.1\n" + extra_seq + "\n")
    # Query sequence matches extra, not primary
    seqs = {"clip_q1": extra_seq}
    hits = s05._batch_check_element_hits(
        seqs, primary, tmp_path, extra_db=extra,
    )
    assert "clip_q1" in hits, f"no hits found: {hits}"
    hit_str = " ".join(hits["clip_q1"])
    assert "bar" in hit_str, f"expected 'bar' in hits, got: {hit_str}"
