"""Test _apply_mask_bed tags sites in host endogenous regions."""
from __future__ import annotations

from scripts.s05.site_discovery import _apply_mask_bed  # new in T10


def test_site_in_mask_gets_tagged(tmp_path):
    bed = tmp_path / "mask.bed"
    bed.write_text("Chr3\t100\t500\tP-Ubi1-maize\t99.5\n")
    sites = [
        {"chr": "Chr3", "pos": 300, "tier": "transgene-positive"},
        {"chr": "Chr3", "pos": 600, "tier": "transgene-positive"},
        {"chr": "Chr5", "pos": 300, "tier": "transgene-positive"},
    ]
    out = _apply_mask_bed(sites, bed)
    assert out[0]["tier"] == "FALSE_NEGATIVE_MASKED"
    assert out[0]["mask_element"] == "P-Ubi1-maize"
    assert out[1]["tier"] == "transgene-positive"  # 600 > 500, outside region
    assert out[2]["tier"] == "transgene-positive"  # different chr


def test_empty_bed_is_noop(tmp_path):
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    sites = [{"chr": "Chr3", "pos": 300, "tier": "transgene-positive"}]
    out = _apply_mask_bed(sites, bed)
    assert out == sites
