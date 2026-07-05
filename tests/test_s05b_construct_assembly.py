"""Unit tests for scripts/s05b_construct_assembly.py (Step 5b)."""
import hashlib
import importlib.util
import shutil
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "s05b",
    Path(__file__).resolve().parent.parent / "scripts" / "s05b_construct_assembly.py",
)
s05b = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(s05b)

requires_blastn = pytest.mark.skipif(
    shutil.which("blastn") is None, reason="blastn not on PATH"
)


def _complex_seq(n: int) -> str:
    """Deterministic, high-complexity DNA (avoids blastn dust masking)."""
    bases = "ACGT"
    out: list[str] = []
    h = hashlib.sha256(b"redgene-s05b").digest()
    while len(out) < n:
        h = hashlib.sha256(h).digest()
        for b in h:
            out.append(bases[b & 3])
            if len(out) >= n:
                break
    return "".join(out)


# --- Task 1: interval merge / coverage / blast6 parse ---
def test_merge_intervals_overlap_and_adjacent():
    assert s05b._merge_intervals([(1, 5), (4, 8), (10, 12), (12, 15)]) == [(1, 8), (10, 15)]


def test_merge_intervals_empty():
    assert s05b._merge_intervals([]) == []


def test_coverage_fraction_basic():
    assert s05b._coverage_fraction([(1, 5), (4, 8), (10, 15)], 20) == pytest.approx(0.7)


def test_coverage_fraction_zero_length():
    assert s05b._coverage_fraction([(1, 5)], 0) == 0.0


def test_parse_blast6_filters_by_pident():
    text = (
        "c1\thostX\t95.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t180\n"
        "c1\thostX\t70.0\t50\t0\t0\t200\t249\t1\t50\t1e-10\t90\n"
        "c2\thostX\t99.0\t30\t0\t0\t60\t31\t1\t30\t1e-20\t60\n"
    )
    hits = s05b._parse_blast6(text, min_pident=90.0)
    assert hits["c1"] == [(1, 100)]
    assert hits["c2"] == [(31, 60)]
    assert "hostX" not in hits


# --- Task 2: classification ---
def test_classify_chimeric_both_present():
    assert s05b._classify_contig(0.5, 0.5) == "chimeric"
    assert s05b._classify_contig(0.20, 0.20) == "chimeric"


def test_classify_construct_element_dominant():
    assert s05b._classify_contig(0.10, 0.80) == "construct"
    assert s05b._classify_contig(0.0, 0.50) == "construct"


def test_classify_host_dominant():
    assert s05b._classify_contig(0.95, 0.0) == "host"
    assert s05b._classify_contig(0.70, 0.10) == "host"


def test_classify_foreign_unknown():
    assert s05b._classify_contig(0.10, 0.10) == "foreign_unknown"
    assert s05b._classify_contig(0.0, 0.0) == "foreign_unknown"
    assert s05b._classify_contig(0.10, 0.30) == "foreign_unknown"


# --- Task 3: blastn coverage runners ---
@requires_blastn
def test_run_blastn_self_match(tmp_path):
    seq = _complex_seq(320)
    q = tmp_path / "q.fa"
    q.write_text(f">c1\n{seq}\n")
    subj = tmp_path / "s.fa"
    subj.write_text(f">ref1\n{seq}\n")
    hits = s05b._run_blastn(q, [subj], min_pident=90.0)
    assert "c1" in hits
    assert s05b._coverage_fraction(hits["c1"], len(seq)) > 0.9


def test_run_blastn_skips_missing_subject(tmp_path):
    q = tmp_path / "q.fa"
    q.write_text(">c1\nACGTACGTAC\n")
    hits = s05b._run_blastn(q, [tmp_path / "nope.fa"], min_pident=90.0)
    assert hits == {}


def test_contig_lengths(tmp_path):
    fa = tmp_path / "contigs.fa"
    fa.write_text(">c1\nACGTACGT\n>c2\nAAA\n")
    assert s05b._contig_lengths(fa) == {"c1": 8, "c2": 3}


# --- Task 4: s05 site parsing ---
def test_parse_s05_sites_filters_candidate(tmp_path):
    s05 = tmp_path / "s05_insert_assembly"
    s05.mkdir()
    (s05 / "insertion_Chr3_16439674_report.txt").write_text(
        "Insertion site: Chr3:16439674\nVerdict: CANDIDATE\n"
    )
    (s05 / "insertion_Chr3_31434949_report.txt").write_text(
        "Insertion site: Chr3:31434949\nVerdict: FALSE_POSITIVE\n"
    )
    (s05 / "insertion_SLM_r2.0ch08_65107378_report.txt").write_text(
        "Verdict: CANDIDATE\n"
    )
    sites = s05b._parse_s05_sites(s05)
    ids = {s.site_id for s in sites}
    assert ids == {"Chr3_16439674", "SLM_r2.0ch08_65107378"}
    by_id = {s.site_id: s for s in sites}
    assert by_id["Chr3_16439674"].chrom == "Chr3"
    assert by_id["Chr3_16439674"].pos == 16439674
    assert by_id["SLM_r2.0ch08_65107378"].chrom == "SLM_r2.0ch08"
    assert by_id["SLM_r2.0ch08_65107378"].pos == 65107378


def test_parse_s05_sites_missing_dir(tmp_path):
    assert s05b._parse_s05_sites(tmp_path / "absent") == []


# --- Task 5: PAF parse + linkage ---
def test_parse_paf_minimal():
    paf = "c1\t300\t0\t300\t+\tChr3\t40000000\t16439600\t16439900\t280\t300\t60\n"
    alns = s05b._parse_paf(paf)
    assert len(alns) == 1
    assert alns[0].contig_id == "c1"
    assert alns[0].chrom == "Chr3"
    assert alns[0].tstart == 16439600
    assert alns[0].tend == 16439900


def test_link_contigs_downstream_and_confidence():
    sites = [s05b.Site("Chr3_16439674", "Chr3", 16439674, "CANDIDATE")]
    alns = [s05b.Aln("c1", "Chr3", 16439700, 16450000)]
    links = s05b._link_contigs_to_sites(alns, sites, {"c1": "chimeric"})
    assert len(links) == 1
    assert links[0].site_id == "Chr3_16439674"
    assert links[0].contig_id == "c1"
    assert links[0].junction_side == "downstream"
    assert links[0].link_confidence == "high"


def test_link_contigs_upstream_medium():
    sites = [s05b.Site("Chr3_16439674", "Chr3", 16439674, "CANDIDATE")]
    alns = [s05b.Aln("c2", "Chr3", 16430000, 16439700)]
    links = s05b._link_contigs_to_sites(alns, sites, {"c2": "construct"})
    assert links[0].junction_side == "upstream"
    assert links[0].link_confidence == "medium"


def test_link_contigs_no_match_wrong_chrom_or_far():
    sites = [s05b.Site("Chr3_16439674", "Chr3", 16439674, "CANDIDATE")]
    alns = [
        s05b.Aln("c3", "Chr9", 16439700, 16450000),
        s05b.Aln("c4", "Chr3", 17000000, 17000300),
    ]
    assert s05b._link_contigs_to_sites(
        alns, sites, {"c3": "construct", "c4": "construct"}
    ) == []


# --- Task 6: output writers ---
def test_write_outputs_basic(tmp_path):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">c_host\nAAAA\n>c_con\nCCCC\n>c_chim\nGGGG\n")
    rows = [
        s05b.ContigRow("c_host", 4, "host", 0.95, 0.0, "", 0.0, 0, ""),
        s05b.ContigRow("c_con", 4, "construct", 0.0, 0.9, "bar", 99.0, 4, ""),
        s05b.ContigRow("c_chim", 4, "chimeric", 0.4, 0.4, "nptII", 88.0, 4, ""),
    ]
    links = [s05b.Link("Chr3_16439674", "CANDIDATE", "c_chim", "downstream", "high")]
    out = tmp_path / "s05b"
    out.mkdir()
    s05b.write_outputs(out, contigs, rows, links)

    fasta = (out / "construct_contigs.fasta").read_text()
    assert ">c_con" in fasta and ">c_chim" in fasta
    assert ">c_host" not in fasta
    assert "class=construct" in fasta

    cls_tsv = (out / "contig_classification.tsv").read_text().splitlines()
    assert cls_tsv[0].startswith("contig_id\t")
    assert any(line.startswith("c_host\t") for line in cls_tsv)

    links_tsv = (out / "site_construct_links.tsv").read_text().splitlines()
    assert links_tsv[0].startswith("site_id\t")
    assert any("c_chim" in line for line in links_tsv)

    summary = (out / "construct_summary.txt").read_text()
    assert "bar" in summary and "nptII" in summary


def test_write_outputs_empty(tmp_path):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text("")
    out = tmp_path / "s05b"
    out.mkdir()
    s05b.write_outputs(out, contigs, [], [])
    assert (out / "construct_contigs.fasta").read_text() == ""
    assert (out / "contig_classification.tsv").read_text().splitlines()[0].startswith("contig_id\t")
    assert (out / "site_construct_links.tsv").read_text().splitlines()[0].startswith("site_id\t")
    assert (out / "construct_summary.txt").exists()


# --- Task 7: run() orchestration exit-0 contract + CLI ---
def test_run_empty_inputs_exit0(tmp_path):
    args = s05b.parse_args_from([
        "--contigs-all", str(tmp_path / "absent_contigs.fasta"),
        "--s03-r1", str(tmp_path / "absent_R1.fq.gz"),
        "--s03-r2", str(tmp_path / "absent_R2.fq.gz"),
        "--s05-dir", str(tmp_path / "absent_s05"),
        "--host-ref", str(tmp_path / "absent_host.fa"),
        "--element-db", str(tmp_path / "absent_db.fa"),
        "--outdir", str(tmp_path / "results"),
        "--sample-name", "tsample",
        "--threads", "1",
    ])
    rc = s05b.run(args)
    assert rc == 0
    out = tmp_path / "results" / "tsample" / "s05b_construct_assembly"
    assert (out / "construct_contigs.fasta").read_text() == ""
    assert (out / "contig_classification.tsv").exists()
    assert (out / "site_construct_links.tsv").exists()
    assert (out / "construct_summary.txt").exists()


def test_help_smoke(capsys):
    with pytest.raises(SystemExit) as exc:
        s05b.parse_args_from(["--help"])
    assert exc.value.code == 0
    assert "construct" in capsys.readouterr().out.lower()
