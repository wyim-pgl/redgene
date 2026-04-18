"""Regression tests for Issue #9 — element_db/build_common_payload.sh hardening.

Covers three minor hardening items that were deferred from the T4 merge:
  * I-1 unterminated-line guard: the bash `while IFS=$'\t' read` loop must
    process a final row that lacks a trailing newline. Without `|| [[ -n ... ]]`
    the last record would be silently skipped.
  * I-2 schema doc: `element_db/common_payload_schema.md` must exist and
    describe all five TSV columns so future maintainers do not have to reverse
    engineer the column order.
  * I-3 mktemp same-fs: the temp file must be created next to the final output
    (`mktemp -p "$(dirname "$OUT")"`) so the subsequent `mv` is an atomic
    rename instead of a cross-filesystem copy.
"""
from __future__ import annotations

import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "element_db" / "build_common_payload.sh"
SCHEMA = REPO_ROOT / "element_db" / "common_payload_schema.md"


def _script_text() -> str:
    return SCRIPT.read_text()


def test_build_common_payload_script_exists():
    assert SCRIPT.exists(), SCRIPT


def test_i1_read_loop_has_unterminated_line_guard():
    """I-1: the `while read` loop must tolerate a missing final newline.

    The canonical bash idiom is ``while ... read -r acc ... || [[ -n "$acc" ]]``.
    Without the trailing ``|| [[ -n "$acc" ]]`` clause ``read`` returns non-zero
    when the final line lacks a newline, and the loop body is skipped for that
    record (silent data loss).
    """
    text = _script_text()
    # Allow a little whitespace flexibility but require the guard on the same
    # line as the `while ... read` construct.
    assert re.search(r"while\s+IFS=.*read\s+-r\s+acc[^\n]*\|\|\s*\[\[\s*-n\s*\"\$acc\"", text), (
        "while-read loop missing `|| [[ -n \"$acc\" ]]` unterminated-line guard"
    )


def test_i2_schema_doc_exists_and_lists_all_columns():
    """I-2: schema doc must be present and enumerate the 5 TSV columns."""
    assert SCHEMA.exists(), f"missing {SCHEMA} (I-2 schema documentation)"
    text = SCHEMA.read_text()
    for col in ("accession", "purpose", "seq_start", "seq_stop", "notes"):
        assert col in text, f"schema doc missing column `{col}`"


def test_i3_mktemp_same_filesystem():
    """I-3: mktemp must target the same directory as the final output.

    Using the system-default ``/tmp`` can make the final ``mv`` a cross-device
    copy, which is slow and non-atomic (a crash mid-copy leaves a half-written
    file at the destination).  ``mktemp -p "$(dirname "$OUT")"`` keeps the
    staging file on the same filesystem so ``mv`` is an atomic rename.
    """
    text = _script_text()
    # Accept `mktemp -p "$(dirname "$OUT")"` or an equivalent with OUT quoted.
    assert re.search(r"mktemp\s+-p\s+\"\$\(dirname\s+\"\$OUT\"\)\"", text), (
        "mktemp call must use `-p \"$(dirname \"$OUT\")\"` for atomic rename"
    )
