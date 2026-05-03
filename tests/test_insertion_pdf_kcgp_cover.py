"""Tests for KCGP cover-page integration in scripts/reports/insertion_pdf.py."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "reports" / "insertion_pdf.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("insertion_pdf", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["insertion_pdf"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def _write_mapping(path: Path, rows: list[tuple[str, str, str]]) -> None:
    """rows: list of (kcgp_id, verdict, internal_site). Other columns stubbed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["kcgp_id\tinternal_site\tverdict\thost_chr\tpos_5p\tspec_version"]
    for kid, verdict, internal in rows:
        lines.append(f"{kid}\t{internal}\t{verdict}\tChr1\t100\ttentative_v0")
    path.write_text("\n".join(lines) + "\n")


def test_load_kcgp_mapping_absent(tmp_path):
    mod = _load_module()
    assert mod.load_kcgp_mapping(tmp_path) is None


def test_load_kcgp_mapping_with_candidate(tmp_path):
    mod = _load_module()
    _write_mapping(tmp_path / "kcgp_mapping.tsv", [
        ("KCGP-OSA-G281-26-0001", "CANDIDATE", "rice_G281_Chr3_16439674"),
        ("-", "FALSE_POSITIVE", "rice_G281_Chr1_8923"),
    ])
    out = mod.load_kcgp_mapping(tmp_path)
    assert out is not None
    assert out["kcgp_id"] == "KCGP-OSA-G281-26-0001"
    assert out["spec_version"] == "tentative_v0"


def test_load_kcgp_mapping_corrupt(tmp_path):
    mod = _load_module()
    p = tmp_path / "kcgp_mapping.tsv"
    p.write_bytes(b"\x00\x01garbage")
    result = mod.load_kcgp_mapping(tmp_path)
    assert result is None or isinstance(result, dict)


def test_cover_with_mapping(tmp_path, monkeypatch):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    _write_mapping(tmp_path / "kcgp_mapping.tsv", [
        ("KCGP-OSA-G281-26-0001", "CANDIDATE", "rice_G281_Chr3_16439674"),
    ])

    captured = {}
    def fake_page_text(pdf, title, lines, **kwargs):
        captured["title"] = title
        captured["lines"] = lines

    monkeypatch.setattr(mod, "_page_text", fake_page_text)

    out_pdf = tmp_path / "cover.pdf"
    with PdfPages(out_pdf) as pdf:
        mod._page_cover(pdf, "rice_G281", {"sample": "rice_G281"}, sample_dir=tmp_path)

    body = "\n".join(captured["lines"])
    assert "KCGP ID" in body
    assert "tentative_v0" in body
    assert "KCGP-OSA-G281-26-0001" in body


def test_cover_skips_dash_only_mapping(tmp_path, monkeypatch):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    _write_mapping(tmp_path / "kcgp_mapping.tsv", [
        ("-", "FALSE_POSITIVE", "rice_G281_Chr1_8923"),
        ("-", "UNKNOWN", "rice_G281_Chr1_99"),
    ])

    captured = {}
    def fake_page_text(pdf, title, lines, **kwargs):
        captured["lines"] = lines

    monkeypatch.setattr(mod, "_page_text", fake_page_text)

    out_pdf = tmp_path / "cover.pdf"
    with PdfPages(out_pdf) as pdf:
        mod._page_cover(pdf, "rice_G281", {"sample": "rice_G281"}, sample_dir=tmp_path)

    body = "\n".join(captured["lines"])
    assert "KCGP ID" not in body
