"""A3 — `canonical_triplet_min_identity` must actually gate Rule 1.

The knob is loaded by `config_loader`, documented in `config.yaml`, and asserted
by `test_verdict_hardening.py` — but `generate_report` used to hand
`compute_verdict` *every* annotated element as `matched_canonical`, regardless of
BLAST identity.  Rule 1 is the highest-priority rule and overrides Filters B/C/D,
so a low-identity spurious triplet could force a CANDIDATE verdict.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from scripts.s05.report import generate_report
from scripts.s05.primitives import InsertionSite
from scripts.s05.verdict import VerdictRules, select_canonical_elements


TRIPLET = {"triplet": {"bar", "P-CaMV35S", "T-ocs"}}


# ---------------------------------------------------------------------------
# The pure helper
# ---------------------------------------------------------------------------

def test_elements_above_threshold_are_selected():
    selected = select_canonical_elements(
        {"bar": 99.5, "P-CaMV35S": 97.0}, min_identity=0.90,
    )

    assert selected == {"bar", "P-CaMV35S"}


def test_elements_below_threshold_are_dropped():
    selected = select_canonical_elements(
        {"bar": 99.5, "P-CaMV35S": 71.2}, min_identity=0.90,
    )

    assert selected == {"bar"}


def test_threshold_is_inclusive():
    assert select_canonical_elements({"bar": 90.0}, min_identity=0.90) == {"bar"}


def test_percent_scale_threshold_is_accepted():
    """0.90 and 90.0 must mean the same thing."""
    as_fraction = select_canonical_elements({"bar": 95.0}, min_identity=0.90)
    as_percent = select_canonical_elements({"bar": 95.0}, min_identity=90.0)

    assert as_fraction == as_percent == {"bar"}


def test_zero_threshold_selects_everything():
    selected = select_canonical_elements({"bar": 12.0}, min_identity=0.0)

    assert selected == {"bar"}


def test_missing_identity_is_not_selected():
    assert select_canonical_elements({}, min_identity=0.90) == set()


# ---------------------------------------------------------------------------
# Wired through generate_report
# ---------------------------------------------------------------------------

def _write_inputs(tmp_path: Path, identity: float) -> tuple[Path, Path, Path]:
    insert_fa = tmp_path / "insert.fa"
    insert_fa.write_text(">contig_1\n" + "ACGT" * 250 + "\n")

    annotation = tmp_path / "element_annotation.tsv"
    lines = ["query\telement\tidentity\tlength\tq_start\tq_end\ts_start\ts_end\tevalue\tsource"]
    for i, elem in enumerate(["bar", "P-CaMV35S", "T-ocs"]):
        start = 1 + i * 200
        lines.append(
            f"contig_1\t{elem}\t{identity}\t150\t{start}\t{start + 149}"
            f"\t1\t150\t1e-50\telement_db"
        )
    annotation.write_text("\n".join(lines) + "\n")

    border = tmp_path / "border_hits.tsv"
    border.write_text("")
    return insert_fa, annotation, border


def _run(tmp_path: Path, identity: float, min_identity: float) -> str:
    insert_fa, annotation, border = _write_inputs(tmp_path, identity)
    site = InsertionSite(site_id="insertion_Chr1_1000", host_chr="Chr1", pos_5p=1000)
    rules = VerdictRules(
        canonical_triplets=TRIPLET,
        canonical_triplet_min_identity=min_identity,
    )

    _, verdict = generate_report(
        insert_fa, annotation, border, site,
        n_rounds=1, status="complete", output_dir=tmp_path,
        host_ref=None, element_db=None, construct_flanking=None, rules=rules,
    )
    return (tmp_path / "insertion_Chr1_1000_report.txt").read_text()


def test_high_identity_triplet_promotes_via_rule_1(tmp_path):
    report = _run(tmp_path, identity=99.5, min_identity=0.90)

    assert "Verdict: CANDIDATE" in report
    assert "canonical_triplet[triplet]" in report


def test_low_identity_triplet_does_not_promote_via_rule_1(tmp_path):
    """All three elements present, but every hit is a 71% match.

    The site may still be CANDIDATE via Rule 6 (elements present, no FP filter
    fired) — what must NOT happen is the Rule 1 override, which would have
    beaten Filters B/C/D had they been evaluated.
    """
    report = _run(tmp_path, identity=71.0, min_identity=0.90)

    assert "canonical_triplet" not in report
    assert "element(s) annotated" in report


def test_default_rules_still_promote_a_perfect_triplet(tmp_path):
    """Guard against the gate accidentally rejecting everything."""
    report = _run(tmp_path, identity=100.0, min_identity=VerdictRules().canonical_triplet_min_identity)

    assert "canonical_triplet[triplet]" in report


def test_annotation_without_an_identity_column_is_not_gated(tmp_path):
    """A legacy element_annotation.tsv carries no `identity` column.

    Identity cannot be measured, so the gate must not silently disable canonical
    promotion for every element in the file.
    """
    insert_fa = tmp_path / "insert.fa"
    insert_fa.write_text(">contig_1\n" + "ACGT" * 250 + "\n")

    annotation = tmp_path / "element_annotation.tsv"
    lines = ["query\telement\tlength\tq_start\tq_end\ts_start\ts_end\tevalue\tsource"]
    for i, elem in enumerate(["bar", "P-CaMV35S", "T-ocs"]):
        start = 1 + i * 200
        lines.append(
            f"contig_1\t{elem}\t150\t{start}\t{start + 149}\t1\t150\t1e-50\telement_db"
        )
    annotation.write_text("\n".join(lines) + "\n")

    border = tmp_path / "border_hits.tsv"
    border.write_text("")
    site = InsertionSite(site_id="insertion_Chr1_1000", host_chr="Chr1", pos_5p=1000)
    rules = VerdictRules(canonical_triplets=TRIPLET, canonical_triplet_min_identity=0.90)

    generate_report(
        insert_fa, annotation, border, site,
        n_rounds=1, status="complete", output_dir=tmp_path,
        host_ref=None, element_db=None, construct_flanking=None, rules=rules,
    )
    report = (tmp_path / "insertion_Chr1_1000_report.txt").read_text()

    assert "canonical_triplet[triplet]" in report


def test_unparseable_identity_does_not_drop_the_element(tmp_path):
    """A garbage identity value must not remove the element from the linear map."""
    insert_fa = tmp_path / "insert.fa"
    insert_fa.write_text(">contig_1\n" + "ACGT" * 250 + "\n")

    annotation = tmp_path / "element_annotation.tsv"
    annotation.write_text(
        "query\telement\tidentity\tlength\tq_start\tq_end\ts_start\ts_end\tevalue\tsource\n"
        "contig_1\tbar\tN/A\t150\t1\t150\t1\t150\t1e-50\telement_db\n"
    )
    border = tmp_path / "border_hits.tsv"
    border.write_text("")
    site = InsertionSite(site_id="insertion_Chr1_1000", host_chr="Chr1", pos_5p=1000)

    generate_report(
        insert_fa, annotation, border, site,
        n_rounds=1, status="complete", output_dir=tmp_path,
        host_ref=None, element_db=None, construct_flanking=None,
        rules=VerdictRules(),
    )
    report = (tmp_path / "insertion_Chr1_1000_report.txt").read_text()

    assert "bar" in report
    assert "1 element(s) annotated" in report
