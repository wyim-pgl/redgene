"""Test _apply_mask_bed tags sites in host endogenous regions."""
from __future__ import annotations

from scripts.s05.site_discovery import _apply_mask_bed  # new in T10
from scripts.s05_insert_assembly import MASKED_SOURCE_TAG


def test_site_in_mask_gets_tagged(tmp_path):
    bed = tmp_path / "mask.bed"
    bed.write_text("Chr3\t100\t500\tP-Ubi1-maize\t99.5\n")
    sites = [
        {"chr": "Chr3", "pos": 300, "tier": "transgene-positive"},
        {"chr": "Chr3", "pos": 600, "tier": "transgene-positive"},
        {"chr": "Chr5", "pos": 300, "tier": "transgene-positive"},
    ]
    out = _apply_mask_bed(sites, bed)
    assert out[0]["tier"] == MASKED_SOURCE_TAG
    assert out[0]["mask_element"] == "P-Ubi1-maize"
    assert out[1]["tier"] == "transgene-positive"  # 600 > 500, outside region
    assert out[2]["tier"] == "transgene-positive"  # different chr


def test_empty_bed_is_noop(tmp_path):
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    sites = [{"chr": "Chr3", "pos": 300, "tier": "transgene-positive"}]
    out = _apply_mask_bed(sites, bed)
    assert out == sites


# Issue #13 — hardening additions below -------------------------------------

def test_masked_source_tag_constant_value():
    """Literal should stay in sync with the historical wire-format tag."""
    assert MASKED_SOURCE_TAG == "FALSE_NEGATIVE_MASKED"


def test_interval_is_half_open_start_inclusive_end_exclusive(tmp_path):
    """BED spec: [start, end) — start is inclusive, end is exclusive."""
    bed = tmp_path / "mask.bed"
    bed.write_text("Chr3\t100\t500\tregionA\n")
    sites = [
        {"chr": "Chr3", "pos": 100, "tier": "transgene-positive"},  # inclusive start
        {"chr": "Chr3", "pos": 499, "tier": "transgene-positive"},  # inside
        {"chr": "Chr3", "pos": 500, "tier": "transgene-positive"},  # exclusive end
    ]
    out = _apply_mask_bed(sites, bed)
    assert out[0]["tier"] == MASKED_SOURCE_TAG
    assert out[1]["tier"] == MASKED_SOURCE_TAG
    assert out[2]["tier"] == "transgene-positive"  # 500 is OUTSIDE [100, 500)


def test_dataclass_branch_uses_host_chr_and_pos_5p(tmp_path):
    """Non-dict site objects must be tagged via setattr on tier/mask_element."""
    class FakeSite:
        def __init__(self, host_chr: str, pos_5p: int):
            self.host_chr = host_chr
            self.pos_5p = pos_5p
            self.tier = "transgene-positive"

    bed = tmp_path / "mask.bed"
    bed.write_text("Chr3\t100\t500\tdataclass-hit\n")
    sites = [FakeSite("Chr3", 300), FakeSite("Chr3", 700)]
    out = _apply_mask_bed(sites, bed)
    assert out[0].tier == MASKED_SOURCE_TAG
    assert out[0].mask_element == "dataclass-hit"
    assert out[1].tier == "transgene-positive"
    assert not hasattr(out[1], "mask_element")


def test_dataclass_branch_falls_back_to_pos_3p(tmp_path):
    """When pos_5p is None/0 the helper should use pos_3p as the site pos."""
    class FakeSite:
        def __init__(self, host_chr: str, pos_5p, pos_3p: int):
            self.host_chr = host_chr
            self.pos_5p = pos_5p
            self.pos_3p = pos_3p
            self.tier = "transgene-positive"

    bed = tmp_path / "mask.bed"
    bed.write_text("Chr3\t100\t500\tfallback3p\n")
    # pos_5p is 0 (falsy) -> helper should consult pos_3p=300.
    sites = [FakeSite("Chr3", 0, 300)]
    out = _apply_mask_bed(sites, bed)
    assert out[0].tier == MASKED_SOURCE_TAG
    assert out[0].mask_element == "fallback3p"


def test_many_regions_lookup_is_fast(tmp_path):
    """Sanity check: large region set should not exhibit pathological slowness.

    This does NOT strictly prove O(log R) but it does lock in a generous upper
    bound so a future regression to O(R) per site would be caught.
    """
    import time

    bed = tmp_path / "mask.bed"
    lines = [f"Chr1\t{i * 1000}\t{i * 1000 + 500}\treg{i}\n" for i in range(20_000)]
    bed.write_text("".join(lines))
    sites = [{"chr": "Chr1", "pos": i * 1000 + 100, "tier": "t"} for i in range(2_000)]
    t0 = time.perf_counter()
    out = _apply_mask_bed(sites, bed)
    elapsed = time.perf_counter() - t0
    # Every site should hit a region; lookup + tagging together must stay
    # under a generous ceiling even on slow CI.
    assert all(s["tier"] == MASKED_SOURCE_TAG for s in out)
    assert elapsed < 2.0, f"_apply_mask_bed too slow: {elapsed:.2f}s"
