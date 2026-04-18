"""Regression tests for Issue #5 — element_db/build_crl_amplicons.sh.

This script reproducibly converts `element_db/crl_amplicons_raw.tsv` into
`element_db/crl_amplicons.fa` using the T5 4-way source-tag convention so
the v1.1 DB rebuild can absorb the 82 CRL GMOMETHODS amplicons as a 5th
source tag alongside EUginius / corn borders / common payload / cas9.

v1.0 scope (this commit): build script + unit tests. The actual DB rebuild
and rice_G281 regression sweep are deferred to v1.1 (see Issue #5 comment).
"""
from __future__ import annotations

import re
import subprocess
import textwrap
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "element_db" / "build_crl_amplicons.sh"
RAW_TSV = REPO_ROOT / "element_db" / "crl_amplicons_raw.tsv"
TAG_SCRIPT = REPO_ROOT / "element_db" / "tag_and_concat.py"


def test_script_exists_and_executable():
    assert SCRIPT.exists(), SCRIPT
    # `chmod +x` so the build can be invoked as-is from cron / Makefile.
    assert SCRIPT.stat().st_mode & 0o111, f"{SCRIPT} must be executable"


def test_script_uses_set_euo_pipefail():
    text = SCRIPT.read_text()
    assert re.search(r"set\s+-euo\s+pipefail", text), (
        "bash build scripts must `set -euo pipefail` (T4 convention)"
    )


def test_script_atomic_move_pattern():
    """Must stage output via mktemp -p same-fs, then mv (T4 atomic contract)."""
    text = SCRIPT.read_text()
    assert re.search(r"mktemp.*-p\s+\"\$\(dirname", text), (
        "build_crl_amplicons.sh must stage via `mktemp -p \"$(dirname \"$OUT\")\"`"
    )
    assert re.search(r"^\s*mv\s+\"\$TMPOUT\"\s+\"\$OUT\"", text, re.MULTILINE), (
        "build_crl_amplicons.sh must finish with `mv \"$TMPOUT\" \"$OUT\"`"
    )


def test_script_converts_single_row(tmp_path):
    """Smoke-test: fabricate a 1-row TSV and run the script end-to-end."""
    # Copy the script into a tmp layout with the expected neighbour files.
    eldb = tmp_path / "element_db"
    eldb.mkdir()
    (eldb / "crl_amplicons_raw.tsv").write_text(
        "method_id\ttype\ttarget\taccession\tlength\tsequence\n"
        "QL-TEST-001\tconstruct-specific\tCaMV P-35S / TEST\tAX123\t12\tACGTACGTACGT\n"
    )
    # The real script locates inputs relative to its own directory via
    # `cd "$(dirname "$0")"`, so copy it into the tmp element_db/.
    (eldb / "build_crl_amplicons.sh").write_bytes(SCRIPT.read_bytes())
    (eldb / "build_crl_amplicons.sh").chmod(0o755)
    result = subprocess.run(
        ["bash", str(eldb / "build_crl_amplicons.sh")],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"stderr={result.stderr}"
    out = eldb / "crl_amplicons.fa"
    assert out.exists()
    headers = [ln for ln in out.read_text().splitlines() if ln.startswith(">")]
    assert len(headers) == 1, headers
    # Header must carry 5+ pipe-separated fields so downstream parsers keep
    # working, and must include the accession + amplicon length annotation.
    h = headers[0]
    assert h.startswith(">construct-specific|QL-TEST-001|"), h
    assert "AX123" in h
    assert "amplicon_12bp" in h or "amplicon_12" in h
    assert "crl_v1" in h or "crl" in h.lower()


def test_script_handles_missing_accession(tmp_path):
    """Rows with `accession == none` must use the literal `CRL-GMOMETHODS`
    fallback so the header never contains the bare string `none`."""
    eldb = tmp_path / "element_db"
    eldb.mkdir()
    (eldb / "crl_amplicons_raw.tsv").write_text(
        "method_id\ttype\ttarget\taccession\tlength\tsequence\n"
        "QL-TEST-002\tconstruct-specific\tT-nos / PG\tnone\t8\tACGTACGT\n"
    )
    (eldb / "build_crl_amplicons.sh").write_bytes(SCRIPT.read_bytes())
    (eldb / "build_crl_amplicons.sh").chmod(0o755)
    result = subprocess.run(
        ["bash", str(eldb / "build_crl_amplicons.sh")],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"stderr={result.stderr}"
    header = next(
        ln for ln in (eldb / "crl_amplicons.fa").read_text().splitlines()
        if ln.startswith(">")
    )
    assert "none" not in header, (
        "raw `none` accession must be replaced with CRL-GMOMETHODS fallback"
    )
    assert "CRL-GMOMETHODS" in header


def test_script_output_covers_all_accepted_rows(tmp_path):
    """Re-running the build against the committed raw TSV must produce a
    header set that is a superset of the currently-committed CRL FASTA.

    The committed fasta was hand-curated (v0) and drops 3 taxon/element
    duplicates (QT-ELE-00-001, QT-TAX-BV-001, QT-TAX-ZM-003). The v1.0
    build script pulls in every raw row so the reproducibility contract is
    "TSV is the single source of truth" — the legacy hand-drops get merged
    back in and will be removed at cd-hit-est time during the v1.1 DB
    rebuild. The test therefore asserts the rebuilt set covers every
    committed header, not strict equality.
    """
    committed_fa = REPO_ROOT / "element_db" / "crl_amplicons.fa"
    committed_raw = REPO_ROOT / "element_db" / "crl_amplicons_raw.tsv"
    if not committed_fa.exists() or not committed_raw.exists():
        import pytest
        pytest.skip("committed CRL files not present")

    eldb = tmp_path / "element_db"
    eldb.mkdir()
    (eldb / "crl_amplicons_raw.tsv").write_bytes(committed_raw.read_bytes())
    (eldb / "build_crl_amplicons.sh").write_bytes(SCRIPT.read_bytes())
    (eldb / "build_crl_amplicons.sh").chmod(0o755)
    result = subprocess.run(
        ["bash", str(eldb / "build_crl_amplicons.sh")],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"stderr={result.stderr}"
    rebuilt_headers = {
        ln for ln in (eldb / "crl_amplicons.fa").read_text().splitlines()
        if ln.startswith(">")
    }
    committed_headers = {
        ln for ln in committed_fa.read_text().splitlines()
        if ln.startswith(">")
    }
    # Rebuilt set must cover every committed header (method_id collision =
    # immediate drift warning).  Extras are acceptable — they get deduped
    # during v1.1 cd-hit-est.
    missing = committed_headers - rebuilt_headers
    assert not missing, (
        f"rebuilt CRL FASTA missing {len(missing)} committed headers: "
        f"{sorted(missing)[:3]}..."
    )


def test_tag_and_concat_srcmap_includes_crl():
    """scripts/tag_and_concat.py SRCMAP must recognise crl_amplicons.fa so
    future `make all` incorporates it with a `|src=crl` suffix."""
    text = TAG_SCRIPT.read_text()
    assert "crl_amplicons.fa" in text, (
        "tag_and_concat.py SRCMAP missing `crl_amplicons.fa` — Issue #5 5th tag"
    )
    # The tag value must be `crl` per the 5-way vocabulary documented in the
    # Issue #10 M-3 whitelist.
    assert re.search(r'"crl_amplicons\.fa"\s*:\s*"crl"', text), (
        "crl_amplicons.fa must map to `crl` source tag"
    )


# ---- New tests for Makefile + DB integration (Issue #5 v1.1 AC) ----------

MAKEFILE = REPO_ROOT / "element_db" / "Makefile"
DB_V2 = REPO_ROOT / "element_db" / "gmo_combined_db_v2.fa"


def test_makefile_srcs_includes_crl_amplicons():
    """element_db/Makefile SRCS list must include crl_amplicons.fa so
    `make all` passes it to tag_and_concat.py and into gmo_combined_db_v2.fa."""
    mk = MAKEFILE.read_text()
    m = re.search(r"^SRCS\s*:=\s*(.+)$", mk, re.MULTILINE)
    assert m, "Makefile missing SRCS := line"
    srcs = m.group(1)
    assert "crl_amplicons.fa" in srcs, (
        "Makefile SRCS must list crl_amplicons.fa so DB rebuild includes CRL sequences"
    )


@pytest.mark.skipif(not DB_V2.exists(), reason="gmo_combined_db_v2.fa not built yet")
def test_combined_db_contains_crl_src_tag():
    """After `make all`, gmo_combined_db_v2.fa must contain at least one
    header with `|src=crl` (the 5th source tag from EU CRL GMOMETHODS)."""
    headers = [
        ln for ln in DB_V2.read_text().splitlines()
        if ln.startswith(">")
    ]
    crl_headers = [h for h in headers if "|src=crl" in h]
    assert crl_headers, (
        f"gmo_combined_db_v2.fa has no |src=crl headers; "
        "run `make clean all` in element_db/ to rebuild with CRL amplicons"
    )


@pytest.mark.skipif(not DB_V2.exists(), reason="gmo_combined_db_v2.fa not built yet")
def test_combined_db_crl_count_within_range():
    """cd-hit-est dedup at 95% must retain between 60 and 85 CRL sequences
    (82 input after hand-curation; some share >=95% identity with EUginius
    seqs already in the DB but the majority should survive deduplication)."""
    headers = [
        ln for ln in DB_V2.read_text().splitlines()
        if ln.startswith(">")
    ]
    crl_count = sum(1 for h in headers if "|src=crl" in h)
    assert 60 <= crl_count <= 85, (
        f"Expected 60–85 CRL sequences after cd-hit-est dedup, got {crl_count}. "
        "Check element_db/gmo_combined_db_v2.fa.clstr for cluster details."
    )
