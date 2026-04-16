"""Audit header writer for RedGene v1.0 regulatory compliance.

Writes per-sample ``audit_header.json`` with four regulatory MUST fields
(team-consensus.md §2.1 item 1, §5):

    R-1  input_sha256        — SHA-256 of R1/R2 fastq inputs
    R-2  pipeline_commit/_dirty — git HEAD hash + clean/dirty flag
    R-3  db_manifest         — parsed element-DB manifest TSV (may be empty)
    R-4  software_versions   — first-line version strings of key tools

Called from ``run_pipeline.py`` at the start of each per-sample loop iteration
so every run is uniquely fingerprinted and reproducible for quarantine audits.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_commit_and_dirty() -> tuple[str, bool]:
    """Return (HEAD SHA, True if working tree is dirty).

    Falls back to ("unknown", False) if the repo is unavailable (e.g. tarball
    deploy), to avoid breaking pipeline startup in unusual environments.
    """
    repo = Path(__file__).resolve().parent.parent
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(repo), text=True, stderr=subprocess.DEVNULL,
        ).strip()
        status = subprocess.check_output(
            ["git", "status", "--porcelain"],
            cwd=str(repo), text=True, stderr=subprocess.DEVNULL,
        )
        return commit, bool(status.strip())
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("[audit] warning: git metadata unavailable", file=sys.stderr)
        return "unknown", False


def _software_versions() -> dict[str, str]:
    tools: dict[str, list[str]] = {
        "bwa": ["bwa"],
        "minimap2": ["minimap2", "--version"],
        "samtools": ["samtools", "--version"],
        "blastn": ["blastn", "-version"],
        "spades": ["spades.py", "--version"],
        "cd-hit-est": ["cd-hit-est", "-h"],
        "fastp": ["fastp", "--version"],
        "python": [sys.executable, "--version"],
    }
    out: dict[str, str] = {}
    for name, cmd in tools.items():
        try:
            res = subprocess.run(cmd, capture_output=True, text=True,
                                 timeout=5, check=False)
            lines = (res.stdout + res.stderr).strip().splitlines()
            out[name] = lines[0] if lines else "unknown"
        except (FileNotFoundError, subprocess.TimeoutExpired):
            out[name] = "not-found"
    return out


def _parse_manifest(path: Path) -> list[dict[str, str]]:
    """Parse a TSV manifest with a header row. Missing file -> empty list.

    The T5 task (element_db/gmo_combined_db_manifest.tsv) will generate this
    file; until then we must not raise.
    """
    if not path.exists():
        return []
    lines = path.read_text().strip().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    entries: list[dict[str, str]] = []
    for row in lines[1:]:
        cols = row.split("\t")
        if len(cols) == len(header):
            entries.append(dict(zip(header, cols)))
    return entries


def write_audit_header(
    *,
    sample: str,
    reads_r1: Path,
    reads_r2: Path,
    db_manifest: Path,
    out_path: Path,
) -> None:
    """Serialize the 4-field audit header to ``out_path`` as indented JSON."""
    commit, dirty = _git_commit_and_dirty()
    data = {
        "sample": sample,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "input_sha256": {
            "r1": _sha256_file(reads_r1),
            "r2": _sha256_file(reads_r2),
        },
        "pipeline_commit": commit,
        "pipeline_dirty": dirty,
        "db_manifest": _parse_manifest(db_manifest),
        "software_versions": _software_versions(),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(data, indent=2))
    print(f"[audit] wrote {out_path}", file=sys.stderr)
