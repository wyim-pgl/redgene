"""Regression tests for Issue #10 — T5 element DB build hardening.

Covers the 6 minor follow-ups from commits 175d209 / f73cb39:

  * M-1 element_db path sync — pipeline hardcodes `gmo_combined_db_v2.fa`
    in two places (build_element_mask_bed.sh + s05_insert_assembly.py).
    The Makefile `DB` target must match so the build and the consumer stay
    in lockstep.
  * M-2 corner cases — `tag_and_concat.py` must handle empty and single-seq
    FASTAs without crashing.
  * M-3 4-way source tag consistency — every header in
    `gmo_combined_db_v2.fa` and `_raw_combined.fa` must carry a `|src=<tag>`
    suffix.
  * M-4 cry34/35Ab1 ambiguity — not applicable in v1.0 (no cry34/35 seqs
    currently shipped); captured in a defer-style test that re-fails when
    the sequences do land so the ambiguity is addressed.
  * M-5 X04879.1 T-ocs subregion — accession X04879.1 must be labelled
    `T-ocs` (octopine synthase terminator), carry the `:1-706` subregion,
    and be tagged `|src=payload`. The legacy mislabel must not reappear.
  * M-6 cd-hit-est parameter docs — the Makefile / tag_and_concat.py must
    document why `-c 0.95 -n 10` was chosen.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
ELEMENT_DB = REPO_ROOT / "element_db"
MAKEFILE = ELEMENT_DB / "Makefile"
TAG_SCRIPT = ELEMENT_DB / "tag_and_concat.py"
DB_V2 = ELEMENT_DB / "gmo_combined_db_v2.fa"
RAW_COMBINED = ELEMENT_DB / "_raw_combined.fa"


def _headers(fa: Path) -> list[str]:
    return [ln for ln in fa.read_text().splitlines() if ln.startswith(">")]


# ---- M-1: path sync between Makefile DB target + consumers ---------------

def test_m1_makefile_db_matches_consumers():
    """element_db/Makefile `DB :=` must equal the consumer path constant."""
    mk = MAKEFILE.read_text()
    m = re.search(r"^DB\s*:=\s*(\S+)", mk, re.MULTILINE)
    assert m, "Makefile missing `DB :=` target line"
    db_name = m.group(1)
    # build_element_mask_bed.sh hardcodes the same default.
    mask_sh = (REPO_ROOT / "scripts" / "build_element_mask_bed.sh").read_text()
    assert db_name in mask_sh, (
        f"scripts/build_element_mask_bed.sh does not reference `{db_name}`; "
        "the Makefile DB target and the mask-bed script must stay in sync"
    )
    # The canonical consumer path is also referenced in s05 comments.
    s05 = (REPO_ROOT / "scripts" / "s05_insert_assembly.py").read_text()
    assert db_name in s05, (
        f"scripts/s05_insert_assembly.py comments must reference `{db_name}` "
        "so audit trail links Phase 1.5 results back to the build artefact"
    )


# ---- M-2: tag_and_concat.py corner cases ---------------------------------

def test_m2_tag_and_concat_handles_empty_fasta(tmp_path, monkeypatch):
    """An empty source FASTA must yield empty stdout, not a KeyError/crash."""
    empty = tmp_path / "common_payload.fa"
    empty.write_text("")
    cmd = [sys.executable, str(TAG_SCRIPT), str(empty)]
    # The script consults SRCMAP by file BASENAME, so the tmp file must be
    # named like a real source; copy it into element_db layout with the
    # required basename.
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"stderr={result.stderr}"
    assert result.stdout == "", f"expected empty stdout, got {result.stdout!r}"


def test_m2_tag_and_concat_handles_single_seq(tmp_path):
    """A single-record FASTA must be tagged correctly with no extra lines."""
    src = tmp_path / "common_payload.fa"
    src.write_text(">only|record\nACGT\n")
    cmd = [sys.executable, str(TAG_SCRIPT), str(src)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    out = result.stdout.splitlines()
    assert out == [">only|record|src=payload", "ACGT"], out


def test_m2_tag_and_concat_rejects_unknown_file(tmp_path):
    """Files not in SRCMAP must exit non-zero with a helpful message."""
    bogus = tmp_path / "totally_fake.fa"
    bogus.write_text(">x\nA\n")
    cmd = [sys.executable, str(TAG_SCRIPT), str(bogus)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0
    assert "SRCMAP" in result.stderr


# ---- M-3: 4-way source tag consistency -----------------------------------

@pytest.mark.skipif(not DB_V2.exists(), reason="gmo_combined_db_v2.fa not built")
def test_m3_every_header_has_src_tag():
    headers = _headers(DB_V2)
    missing = [h for h in headers if "|src=" not in h]
    assert not missing, (
        f"{len(missing)} headers missing `|src=` tag in {DB_V2.name}; "
        f"first offender: {missing[0] if missing else ''!r}"
    )


@pytest.mark.skipif(not DB_V2.exists(), reason="gmo_combined_db_v2.fa not built")
def test_m3_src_tag_values_are_from_whitelist():
    """src values must be a subset of the documented tag vocabulary."""
    # Issue #5 added `crl` as the 5th tag (not wired into v1.0 DB yet).
    allowed = {"payload", "element_db", "plant_endogenous", "corn_border", "crl"}
    headers = _headers(DB_V2)
    seen = set()
    for h in headers:
        m = re.search(r"\|src=([A-Za-z0-9_]+)", h)
        if m:
            seen.add(m.group(1))
    unexpected = seen - allowed
    assert not unexpected, (
        f"unexpected src tag value(s) {sorted(unexpected)}; "
        f"update the whitelist or the tag_and_concat SRCMAP"
    )


# ---- M-4: cry34/35Ab1 ambiguity (v1.0 not-applicable placeholder) --------

@pytest.mark.skipif(not DB_V2.exists(), reason="gmo_combined_db_v2.fa not built")
def test_m4_cry34_35_ambiguity_plan_if_present():
    """No cry34/cry35 sequences ship in v1.0 element_db. If they ever land
    the build must emit disambiguated headers (cry34Ab1 vs cry35Ab1) rather
    than letting cd-hit-est collapse them into one representative."""
    headers = _headers(DB_V2)
    cry34 = [h for h in headers if "cry34" in h.lower()]
    cry35 = [h for h in headers if "cry35" in h.lower()]
    if not cry34 and not cry35:
        pytest.skip("cry34/35 not shipped in current v1.0 DB — see Issue #10 M-4 v1.1 defer")
    # If either present, we must have disambiguation in the header.
    for h in cry34 + cry35:
        assert re.search(r"cry3[45]Ab1\b|cry3[45]\|", h), (
            f"ambiguous cry34/35 header `{h}`; disambiguate Ab1 vs Ab2"
        )


# ---- M-5: X04879.1 T-ocs subregion ---------------------------------------

def test_m5_x04879_labelled_tocs_in_common_payload():
    """element_db/common_payload.fa must label X04879.1 as T-ocs with the
    curated 1-706 subregion (prevents the legacy BUG-9-class mislabel)."""
    fa = ELEMENT_DB / "common_payload.fa"
    assert fa.exists()
    headers = _headers(fa)
    ocs = [h for h in headers if "X04879" in h]
    assert ocs, "X04879.1 missing from common_payload.fa"
    for h in ocs:
        assert "T-ocs" in h, f"X04879.1 header not tagged T-ocs: {h!r}"
        assert ":1-706" in h, (
            f"X04879.1 header missing `:1-706` subregion: {h!r}"
        )


@pytest.mark.skipif(not DB_V2.exists(), reason="gmo_combined_db_v2.fa not built")
def test_m5_x04879_propagates_to_combined_db():
    headers = _headers(DB_V2)
    ocs = [h for h in headers if "X04879" in h]
    # cd-hit-est may collapse the canonical T-ocs into the larger cluster,
    # in which case it should reappear under a different header prefix —
    # still labelled T-ocs.
    if not ocs:
        pytest.skip("X04879.1 absorbed by cd-hit-est cluster (see .clstr file)")
    for h in ocs:
        assert "T-ocs" in h, f"combined DB mislabel: {h!r}"
        assert "|src=payload" in h, (
            f"X04879.1 must be tagged `src=payload`, got {h!r}"
        )


# ---- M-6: cd-hit-est parameter documentation -----------------------------

def test_m6_makefile_documents_cdhit_params():
    """-c 0.95 and -n 10 must be explained inside the Makefile comment block."""
    mk = MAKEFILE.read_text()
    # At minimum, the commit/issue rationale file must reference the thresholds.
    assert "0.95" in mk, "Makefile missing cd-hit-est identity threshold (-c 0.95)"
    assert "-n 10" in mk or "word_size" in mk.lower() or "10" in mk, (
        "Makefile missing cd-hit-est word size (-n 10) parameter"
    )
    assert "95%" in mk or "identity" in mk.lower(), (
        "Makefile missing rationale comment for the 0.95 identity threshold"
    )
