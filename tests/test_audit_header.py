"""Verify the AC-6 audit header writer produces the 4 regulatory fields.

Covers regulatory MUST R-1/R-2/R-3/R-4 per docs/team-review/team-consensus.md:
  R-1 input_sha256  — SHA-256 of R1 and R2
  R-2 pipeline_commit + pipeline_dirty — git commit hash and clean/dirty flag
  R-3 db_manifest   — parsed element-DB manifest TSV
  R-4 software_versions — versions of bwa/minimap2/samtools/etc.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest  # noqa: F401 - imported for pytest fixture discovery parity

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts._write_audit_header import write_audit_header  # noqa: E402


def test_audit_header_contains_4_fields(tmp_path):
    r1 = tmp_path / "reads_R1.fq.gz"
    r1.write_bytes(b"@r\nACGT\n+\n!!!!\n")
    r2 = tmp_path / "reads_R2.fq.gz"
    r2.write_bytes(b"@r\nACGT\n+\n!!!!\n")
    db_manifest = tmp_path / "element_db_manifest.tsv"
    db_manifest.write_text(
        "name\tmd5\tbuild_date\tseq_count\n"
        "element_db\tdeadbeef\t2026-04-16\t146\n"
    )

    out = tmp_path / "audit.json"
    write_audit_header(
        sample="rice_G281",
        reads_r1=r1,
        reads_r2=r2,
        db_manifest=db_manifest,
        out_path=out,
    )

    data = json.loads(out.read_text())

    # R-1 input_sha256 — 64-char hex digest
    assert "input_sha256" in data
    assert data["input_sha256"]["r1"].startswith(
        ("a", "b", "c", "d", "e", "f") + tuple("0123456789")
    )
    assert len(data["input_sha256"]["r1"]) == 64  # SHA-256 hex
    assert len(data["input_sha256"]["r2"]) == 64

    # R-2 pipeline_commit + pipeline_dirty
    assert "pipeline_commit" in data
    assert isinstance(data["pipeline_commit"], str)
    assert len(data["pipeline_commit"]) >= 7  # git short-hash lower bound
    assert "pipeline_dirty" in data
    assert isinstance(data["pipeline_dirty"], bool)

    # R-3 db_manifest
    assert "db_manifest" in data
    assert data["db_manifest"][0]["name"] == "element_db"
    assert data["db_manifest"][0]["md5"] == "deadbeef"

    # R-4 software_versions
    assert "software_versions" in data
    assert any(k.startswith("bwa") for k in data["software_versions"])


def test_audit_header_missing_manifest_returns_empty_list(tmp_path):
    """Graceful handling when element_db manifest TSV does not yet exist
    (T5 will create it; T1 must not hard-fail before then)."""
    r1 = tmp_path / "R1.fq.gz"
    r1.write_bytes(b"@r\nA\n+\n!\n")
    r2 = tmp_path / "R2.fq.gz"
    r2.write_bytes(b"@r\nA\n+\n!\n")
    missing_manifest = tmp_path / "does_not_exist.tsv"
    out = tmp_path / "audit.json"

    write_audit_header(
        sample="test",
        reads_r1=r1,
        reads_r2=r2,
        db_manifest=missing_manifest,
        out_path=out,
    )

    data = json.loads(out.read_text())
    assert data["db_manifest"] == []
