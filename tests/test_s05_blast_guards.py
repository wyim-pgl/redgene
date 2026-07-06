"""Robustness guards: local/host BLAST failures must degrade gracefully,
not crash the step or trust a partial output file."""
import subprocess
from pathlib import Path

from scripts.s05 import annotation, filter_a_host


class _FakeProc:
    def __init__(self, returncode: int):
        self.returncode = returncode
        self.stdout = ""
        self.stderr = ""


def test_run_local_blast_makeblastdb_failure_returns_empty(tmp_path, monkeypatch):
    q = tmp_path / "insert.fasta"
    q.write_text(">c1\nACGTACGTACGT\n")
    db = tmp_path / "elem.fa"
    db.write_text(">e1\nACGTACGTACGT\n")

    monkeypatch.setattr(annotation.subprocess, "run",
                        lambda *a, **k: _FakeProc(1))
    # Must not raise CalledProcessError; returns [] gracefully.
    assert annotation._run_local_blast(q, db, tmp_path) == []


def test_run_local_blast_blastn_failure_returns_empty(tmp_path, monkeypatch):
    q = tmp_path / "insert.fasta"
    q.write_text(">c1\nACGTACGTACGT\n")
    db = tmp_path / "elem.fa"
    db.write_text(">e1\nACGTACGTACGT\n")

    calls = {"n": 0}

    def fake_run(cmd, *a, **k):
        calls["n"] += 1
        # makeblastdb succeeds (call 1), blastn fails (call 2)
        return _FakeProc(0 if calls["n"] == 1 else 1)

    monkeypatch.setattr(annotation.subprocess, "run", fake_run)
    assert annotation._run_local_blast(q, db, tmp_path) == []


def test_blast_insert_vs_host_failure_fails_open(tmp_path, monkeypatch):
    insert = tmp_path / "ins.fasta"
    insert.write_text(">c1\n" + "ACGT" * 50 + "\n")  # 200 bp
    host = tmp_path / "host.fa"
    host.write_text(">chr1\nACGTACGTAC\n")

    # blastn fails and writes no output file → must fail open to no-host-coverage.
    monkeypatch.setattr(filter_a_host.subprocess, "run",
                        lambda *a, **k: _FakeProc(1))
    frac, host_bp, insert_len, gap = filter_a_host._blast_insert_vs_host(
        insert, host, tmp_path)
    assert frac == 0.0
    assert host_bp == 0
    assert insert_len == 200
    assert gap == 200  # whole insert unmeasured → treated as one foreign gap


def test_blast_insert_vs_host_success_path_unchanged(tmp_path, monkeypatch):
    """On returncode 0 with a real host hit, behaviour is unchanged."""
    insert = tmp_path / "ins.fasta"
    insert.write_text(">c1\n" + "A" * 100 + "\n")
    host = tmp_path / "host.fa"
    host.write_text(">chr1\nAAAA\n")

    out_path = tmp_path / "_ins_vs_host_chrom.tsv"

    def fake_run(cmd, *a, **k):
        # emulate blastn writing a host hit covering positions 1..100
        out_path.write_text("c1\t1\t100\tchr1\t99.0\t100\n")
        return _FakeProc(0)

    monkeypatch.setattr(filter_a_host.subprocess, "run", fake_run)
    frac, host_bp, insert_len, gap = filter_a_host._blast_insert_vs_host(
        insert, host, tmp_path)
    assert host_bp == 100
    assert frac == 1.0
    assert gap == 0
