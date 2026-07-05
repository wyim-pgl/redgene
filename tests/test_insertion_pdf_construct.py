"""Tests for the step-8 PDF Construct Inventory section (step 5b integration)."""
import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "insertion_pdf",
    Path(__file__).resolve().parent.parent / "scripts" / "reports" / "insertion_pdf.py",
)
ipdf = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(ipdf)


def _write_s05b(d: Path):
    d.mkdir(parents=True, exist_ok=True)
    (d / "construct_summary.txt").write_text(
        "# Step 5b construct inventory\n"
        "total_construct_bp\t12000\n"
        "contigs_construct\t3\ncontigs_chimeric\t1\n"
        "contigs_foreign_unknown\t1\ncontigs_host\t10\n"
        "distinct_elements\t2\nelement_inventory\tbar,nptII\n"
        "sites_linked\t1\nnovel_payload_detected\tyes\n"
    )
    (d / "site_construct_links.tsv").write_text(
        "site_id\tverdict\tcontig_id\tjunction_side\tlink_confidence\n"
        "Chr3_16439674\tCANDIDATE\tc_chim\tdownstream\thigh\n"
    )
    (d / "contig_classification.tsv").write_text(
        "contig_id\tlength\tclass\thost_frac\telement_frac\ttop_element\t"
        "element_pident\telement_alnlen\tremote_hit\n"
        "c_chim\t4000\tchimeric\t0.4\t0.4\tnptII\t88.0\t400\t\n"
    )


def test_load_construct_inventory_present(tmp_path):
    _write_s05b(tmp_path / "s05b_construct_assembly")
    inv = ipdf._load_construct_inventory(tmp_path / "s05b_construct_assembly")
    assert inv is not None
    assert inv["summary"]["element_inventory"] == "bar,nptII"
    assert inv["summary"]["novel_payload_detected"] == "yes"
    assert len(inv["links"]) == 1
    assert inv["links"][0]["contig_id"] == "c_chim"


def test_load_construct_inventory_absent(tmp_path):
    assert ipdf._load_construct_inventory(tmp_path / "absent") is None


def test_page_construct_inventory_renders(tmp_path):
    """The page renderer runs without error and adds a page when inventory present."""
    from matplotlib.backends.backend_pdf import PdfPages

    _write_s05b(tmp_path / "s05b_construct_assembly")
    inv = ipdf._load_construct_inventory(tmp_path / "s05b_construct_assembly")
    out_pdf = tmp_path / "out.pdf"
    with PdfPages(out_pdf) as pdf:
        ipdf._page_construct_inventory(pdf, inv, "tsample")
    assert out_pdf.stat().st_size > 0
