"""Unit tests for the T-DNA border motif scan (A1).

Historically `annotate_insert` wrote two FASTA records, `RB_consensus` and
`LB_consensus`, holding the *same* 25-mer.  Every genuine border therefore
produced two identical BLAST rows, the RB/LB labels meant nothing, and the
second border repeat that both of this project's constructs actually carry was
never queried.

Provenance of the two repeats asserted below — record-scoped scan of the first
FASTA record of each construct in `db/`:

    db/G281_construct.fa  >rice_G281   (8,958 bp)
        278  TGACAGGATATATTGGCGGGTAAAC
      6,556  TGGCAGGATATATTGTGGTGTAAAC
    db/Cas9_construct.fa  >tomato_A2_3 (16,419 bp)
          0  TGGCAGGATATATTGTGGTGTAAAC
     10,141  TGACAGGATATATTGGCGGGTAAAC
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from scripts.s05.annotation import (
    TDNA_BORDER_REPEATS,
    _dedupe_border_hits,
    _run_border_blast,
    _write_border_query,
)
from scripts.s05.primitives import read_fasta, revcomp


# ---------------------------------------------------------------------------
# The repeat table itself
# ---------------------------------------------------------------------------

def test_border_repeats_are_distinct():
    seqs = [seq for _, seq in TDNA_BORDER_REPEATS]

    assert len(seqs) >= 2
    assert len(set(seqs)) == len(seqs), "duplicate border query (the A1 bug)"


def test_border_repeats_have_distinct_labels():
    labels = [label for label, _ in TDNA_BORDER_REPEATS]

    assert len(set(labels)) == len(labels)


def test_border_repeats_are_25mers():
    for label, seq in TDNA_BORDER_REPEATS:
        assert len(seq) == 25, f"{label} is {len(seq)}bp, expected the 25-bp repeat"
        assert set(seq) <= set("ACGT"), f"{label} has non-ACGT bases"


def test_labels_match_the_primary_genbank_deposits():
    """RB/LB assignment is sourced, not guessed.

    GenBank J01825 — "Ti plasmid (from A.tumefaciens, nopaline strain T37),
    T-DNA 5' (left) border" — contains TGGCAGGATATATTGTGGTGTAAAC.
    GenBank J01826 — "...T-DNA 3' (right) border", COMMENT "right border at
    base 158" — contains TGACAGGATATATTGGCGGGTAAAC at bases 158-182.
    Independently reproduced by pCAMBIA-1300 (AF234296) and pBIN19 (U09365).
    """
    by_label = dict(TDNA_BORDER_REPEATS)

    assert by_label["TDNA_RB"] == "TGACAGGATATATTGGCGGGTAAAC"
    assert by_label["TDNA_LB"] == "TGGCAGGATATATTGTGGTGTAAAC"


def test_right_border_is_listed_first():
    """T-DNA transfer runs RB -> LB; keep the table in that order."""
    assert TDNA_BORDER_REPEATS[0][0] == "TDNA_RB"


def test_border_query_fasta_has_one_record_per_repeat(tmp_path):
    query = tmp_path / "_borders.fa"

    _write_border_query(query)

    records = read_fasta(query)
    assert len(records) == len(TDNA_BORDER_REPEATS)
    assert len(set(records.values())) == len(TDNA_BORDER_REPEATS)


# ---------------------------------------------------------------------------
# HSP de-duplication — the fix for the ~2x inflated border count
# ---------------------------------------------------------------------------

def _row(qseqid, pident, length, sstart, send, sseqid="insert_1"):
    return [qseqid, sseqid, str(pident), str(length), "1", str(length),
            str(sstart), str(send)]


def test_overlapping_hsps_at_one_locus_collapse_to_one_row():
    rows = [
        _row("TDNA_border_A", 100.0, 25, 300, 324),
        _row("TDNA_border_B", 72.0, 25, 300, 324),   # cross-hit, same locus
    ]

    deduped = _dedupe_border_hits(rows)

    assert len(deduped) == 1
    assert deduped[0][0] == "TDNA_border_A"  # best hit wins


def test_distinct_loci_are_both_kept():
    rows = [
        _row("TDNA_border_A", 100.0, 25, 300, 324),
        _row("TDNA_border_B", 100.0, 25, 9000, 9024),
    ]

    deduped = _dedupe_border_hits(rows)

    assert len(deduped) == 2


def test_dedupe_is_strand_agnostic():
    """A minus-strand hit reports sstart > send; it is the same locus."""
    rows = [
        _row("TDNA_border_A", 100.0, 25, 324, 300),
        _row("TDNA_border_A", 98.0, 25, 300, 324),
    ]

    deduped = _dedupe_border_hits(rows)

    assert len(deduped) == 1


def test_dedupe_keeps_the_higher_scoring_hit():
    rows = [
        _row("TDNA_border_B", 80.0, 20, 300, 319),
        _row("TDNA_border_A", 100.0, 25, 300, 324),
    ]

    deduped = _dedupe_border_hits(rows)

    assert len(deduped) == 1
    assert deduped[0][0] == "TDNA_border_A"


def test_dedupe_preserves_outfmt6_column_count():
    """scripts/viz/plot_insert_structure.py requires >= 8 columns."""
    deduped = _dedupe_border_hits([_row("TDNA_border_A", 100.0, 25, 300, 324)])

    assert len(deduped[0]) >= 8


def test_dedupe_of_nothing_is_nothing():
    assert _dedupe_border_hits([]) == []


def test_dedupe_tie_is_broken_independently_of_input_order():
    """Two equally scoring hits at one locus must resolve the same way.

    Otherwise the surviving label depends on the order blastn happened to emit
    its HSPs.
    """
    a = _row("TDNA_border_A", 100.0, 25, 300, 324)
    b = _row("TDNA_border_B", 100.0, 25, 300, 324)

    forward = _dedupe_border_hits([a, b])
    reverse = _dedupe_border_hits([b, a])

    assert forward == reverse
    assert len(forward) == 1


# ---------------------------------------------------------------------------
# End-to-end against real blastn
# ---------------------------------------------------------------------------

@pytest.mark.skipif(shutil.which("blastn") is None, reason="blastn not on PATH")
def test_scan_finds_each_border_exactly_once(tmp_path):
    """An insert carrying both repeats once must report exactly two loci.

    The pre-fix code reported four rows for this input (two identical queries
    x two loci), which is what inflated `T-DNA borders found:`.
    """
    (_, rep_a), (_, rep_b) = TDNA_BORDER_REPEATS[0], TDNA_BORDER_REPEATS[1]
    filler = "ACGTTGCA" * 60          # 480 bp of neutral spacer
    insert = rep_a + filler + revcomp(rep_b)

    insert_fa = tmp_path / "insert.fa"
    insert_fa.write_text(f">insert_1\n{insert}\n")
    border_tsv = tmp_path / "border_hits.tsv"

    _run_border_blast(insert_fa, tmp_path, border_tsv)

    rows = [ln for ln in border_tsv.read_text().splitlines() if ln.strip()]
    assert len(rows) == 2, f"expected 2 border loci, got {len(rows)}:\n" + "\n".join(rows)
    assert {r.split("\t")[0] for r in rows} == {
        TDNA_BORDER_REPEATS[0][0], TDNA_BORDER_REPEATS[1][0]
    }


@pytest.mark.skipif(shutil.which("blastn") is None, reason="blastn not on PATH")
def test_short_partial_hits_are_not_borders(tmp_path):
    """A 12-bp scrap of a repeat is not a border.

    Querying two similar 25-mers at `-word_size 7` surfaces short spurious HSPs
    that the single-query version never produced.  Anything much shorter than
    the 25-bp repeat must be discarded, or the border count regains the
    inflation the dedupe removed.
    """
    _, rep_a = TDNA_BORDER_REPEATS[0]
    insert = "ACGTTGCA" * 30 + rep_a[:12] + "TGCAACGT" * 30

    insert_fa = tmp_path / "insert.fa"
    insert_fa.write_text(f">insert_1\n{insert}\n")
    border_tsv = tmp_path / "border_hits.tsv"

    _run_border_blast(insert_fa, tmp_path, border_tsv)

    assert border_tsv.read_text().strip() == ""


@pytest.mark.skipif(shutil.which("blastn") is None, reason="blastn not on PATH")
@pytest.mark.skipif(not Path("db/G281_construct.fa").exists(),
                    reason="db/G281_construct.fa not present")
def test_real_rice_construct_reports_exactly_two_borders(tmp_path):
    """rice_G281's construct reference is byte-identical to pCAMBIA-1300 (AF234296).

    AF234296 annotates `misc_feature 298..323 /note="right border T-DNA repeat"`
    and `misc_feature 6557..6582 /note="left border repeat from C58 T-DNA"`.
    The 25-mers themselves sit at 279 (RB) and 6,557 (LB) — the depositor's RB
    feature is offset ~19 bp onto the 3' flank, and there is no other border
    repeat in that region.

    Regression for two defects at once: the old scan reported 4 rows (2 loci x
    2 identical queries), and a naive two-query fix reports 4 again (2 loci + 2
    short spurious HSPs at 4,321 and 7,747).
    """
    record: list[str] = []
    started = False
    for line in Path("db/G281_construct.fa").read_text().splitlines():
        if line.startswith(">"):
            if started:
                break
            started = True
            continue
        record.append(line.strip())

    insert_fa = tmp_path / "rice_G281.fa"
    insert_fa.write_text(">rice_G281\n" + "".join(record) + "\n")
    border_tsv = tmp_path / "border_hits.tsv"

    _run_border_blast(insert_fa, tmp_path, border_tsv)

    rows = [ln.split("\t") for ln in border_tsv.read_text().splitlines() if ln.strip()]
    at = {r[0]: min(int(r[6]), int(r[7])) for r in rows}

    assert len(rows) == 2, f"expected 2 border loci, got {rows}"
    assert at == {"TDNA_RB": 279, "TDNA_LB": 6557}, f"got {at}"


@pytest.mark.skipif(shutil.which("blastn") is None, reason="blastn not on PATH")
def test_scan_on_borderless_insert_writes_empty_file(tmp_path):
    insert_fa = tmp_path / "insert.fa"
    insert_fa.write_text(">insert_1\n" + "AC" * 300 + "\n")
    border_tsv = tmp_path / "border_hits.tsv"

    _run_border_blast(insert_fa, tmp_path, border_tsv)

    assert border_tsv.exists()
    assert border_tsv.read_text().strip() == ""
