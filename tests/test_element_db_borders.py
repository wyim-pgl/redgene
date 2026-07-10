"""The T-DNA border records in the tracked element DBs must hold real borders.

The `canonical_v1` entries were built by slicing AF485783.1 at the wrong
coordinates. AF485783.1 is the pBI121 binary vector (14,758 bp); its annotated
borders are on the complement strand:

    T-DNA right border  complement(2454..2478)  -> TGACAGGATATATTGGCGGGTAAAC
    T-DNA left border   complement(8621..8646)  -> TGGCAGGATATATTGTGGTGTAAACA

The shipped records instead carried bases 1-25 of the record
(`TGAGCGTCGCAAAGGCGCTCGGTCT`) as "LB" and a 25-bp window near 14,040
(`GGCCTCGGCCTGAGAGCCAAAACAC`) as "RB". Neither is a border repeat: both lack the
invariant `CAGGATATAT` core, and neither occurs in any construct in `db/`.

These entries feed `classify.py`'s blastn-short tier check (min alignment 20 bp),
so a 25-bp record is short enough to matter there.
"""
from __future__ import annotations

from pathlib import Path

import pytest


# From AF485783.1 (pBI121) FEATURES, complement strand — the accession the
# record headers already cite. The LB feature spans 26 bp; the 25-bp repeat is
# its first 25 bases. Independently: GenBank J01825 (LB) / J01826 (RB).
CANONICAL_RB = "TGACAGGATATATTGGCGGGTAAAC"
CANONICAL_LB = "TGGCAGGATATATTGTGGTGTAAAC"

#: The wrong sequences that shipped as canonical_v1.
BOGUS = {"TGAGCGTCGCAAAGGCGCTCGGTCT", "GGCCTCGGCCTGAGAGCCAAAACAC"}

#: Tracked FASTAs that carry the border records.
DB_FILES = [
    Path("element_db/gmo_combined_db_v2.fa"),
    Path("db/gmo_all_combined_db.fa"),
    Path("db/gmo_corn_combined_db.fa"),
]


def _border_records(path: Path) -> dict[str, str]:
    """Return {'LB-TDNA': seq, 'RB-TDNA': seq} for the canonical border entries."""
    out: dict[str, str] = {}
    name, buf = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            fields = line[1:].split("|")
            name = fields[1] if len(fields) > 1 else line[1:]
            name = name if name in ("LB-TDNA", "RB-TDNA") else None
            buf = []
        elif name:
            buf.append(line.strip().upper())
    if name:
        out[name] = "".join(buf)
    return out


@pytest.mark.parametrize("db", DB_FILES, ids=lambda p: p.name)
def test_border_records_hold_the_canonical_repeats(db):
    if not db.exists():
        pytest.skip(f"{db} not present")

    records = _border_records(db)

    assert records.get("RB-TDNA") == CANONICAL_RB, f"{db}: RB-TDNA is wrong"
    assert records.get("LB-TDNA") == CANONICAL_LB, f"{db}: LB-TDNA is wrong"


@pytest.mark.parametrize("db", DB_FILES, ids=lambda p: p.name)
def test_bogus_border_sequences_are_gone(db):
    if not db.exists():
        pytest.skip(f"{db} not present")

    text = db.read_text().upper()

    for bad in BOGUS:
        assert bad not in text, f"{db} still carries the mis-sliced AF485783.1 window {bad}"


@pytest.mark.parametrize("db", DB_FILES, ids=lambda p: p.name)
def test_border_records_carry_the_invariant_core(db):
    """Every T-DNA border repeat contains the invariant CAGGATATAT core."""
    if not db.exists():
        pytest.skip(f"{db} not present")

    for name, seq in _border_records(db).items():
        assert "CAGGATATAT" in seq, f"{db}: {name} lacks the border core"
        assert len(seq) == 25, f"{db}: {name} is {len(seq)}bp, expected 25"


def test_manifest_md5_matches_the_database():
    """The audit manifest must not advertise a stale checksum."""
    import hashlib

    db = Path("element_db/gmo_combined_db_v2.fa")
    manifest = Path("element_db/gmo_combined_db_manifest.tsv")
    if not (db.exists() and manifest.exists()):
        pytest.skip("element DB or manifest not present")

    rows = [ln.split("\t") for ln in manifest.read_text().splitlines() if ln.strip()]
    header, data = rows[0], rows[1]
    fields = dict(zip(header, data))

    actual_md5 = hashlib.md5(db.read_bytes()).hexdigest()
    actual_count = db.read_text().count(">")

    assert fields["md5"] == actual_md5
    assert int(fields["seq_count"]) == actual_count
