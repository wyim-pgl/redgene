# compute_verdict Filter Priority

**Issue #3 (v1.1 full wire-in) — 2026-04-18**

This document records the canonical rule-priority order for `compute_verdict`
in `scripts/s05/verdict.py` and the rationale for each decision.  It serves
as the authoritative reference for anyone modifying filter logic.

---

## Canonical priority (first-match wins)

| Rule | Name | When it fires | Returns |
|------|------|---------------|---------|
| 1 | canonical_triplet | every element in a configured triplet is in `matched_canonical` **and** `host_fraction < cand_host_fraction_max` | CANDIDATE |
| 1.5 | host_endogenous | every annotated element matched host genome at Tier 1/2 BLAST | FALSE_POSITIVE |
| 2 | Filter B: flanking | site falls inside a construct-flanking host region (slop-expanded) | FALSE_POSITIVE |
| 3 | Filter C: chimeric | host-aligned portions span ≥2 off-target chromosomes | FALSE_POSITIVE |
| 4 | Filter D: construct+host | construct_frac ≥ 0.25 AND combined_frac ≥ 0.85 | FALSE_POSITIVE |
| 5 | Filter A: host+gap | host_fraction ≥ 0.80 AND largest_gap < 500 bp | FALSE_POSITIVE |
| 6 | elements present | elements annotated; all FP filters passed | CANDIDATE |
| 7 | UNKNOWN→FP | no elements, host_fraction ≥ 0.85, construct_frac ≤ 0.05 | FALSE_POSITIVE |
| 8 | fallthrough | — | UNKNOWN |

---

## Rationale for non-obvious orderings

### Rule 1 (canonical_triplet) runs before host_endogenous (Rule 1.5)

A canonical transgene marker set (e.g. `{bar, P-CaMV35S, T-ocs}`) provides
extremely strong biological evidence of a real insertion.  Even if one of
those markers is also present in the host genome (e.g. `P-CaMV35S` as a
CpG-island homolog), the presence of ALL three together is not a host
artifact.  Placing Rule 1 first ensures such sites are not silently discarded
by the host_endogenous rule.

### Filter A (host+gap) runs AFTER Filters B, C, D (Rules 2–4)

The old inline code in `generate_report()` ran Filter A first (short-circuit
efficiency).  `compute_verdict` places it after B/C/D because:

1. Filters B and C detect qualitatively different FP mechanisms (flanking
   host DNA in the construct reference; multi-locus chimeric assembly) that
   should be named explicitly in the reason string rather than collapsing to
   the generic "host genome" message.
2. In practice, the final CANDIDATE/FP outcome is **identical** between the
   two orderings because the filters are mutually exclusive in the common
   cases.  Sites that fail Filter A (high host + small gap) never pass
   Filters B or C.

### Rule 6 does NOT re-check host_fraction

The original `compute_verdict` had `ev.host_fraction < cand_host_fraction_max`
in Rule 6 as an extra guard.  This was **incorrect** and caused a regression
on rice_G281 Chr3:16,439,674 (host_fraction=87.4%, gap=1,024 bp):

- Filter A does not fire because the gap exceeds 500 bp (genuine T-DNA span).
- Rule 6 would return UNKNOWN instead of CANDIDATE because 0.874 ≥ 0.80.

The correct semantics are: *if a site has annotated elements and survived all
FP filters (Rules 1–5), it is CANDIDATE regardless of host_fraction.*
Filter A already handles the "high host_fraction + tiny gap" FP case.  The
fix (Issue #3 TDD) removes the redundant check from Rule 6.

---

## Golden regression cases

| Sample | Site | Expected | Key evidence |
|--------|------|----------|-------------|
| rice_G281 | Chr3:16,439,674 | CANDIDATE | 87.4% host, 1,024 bp gap — Filter A must NOT fire |
| rice_G281 | Chr11:12,123,277 | FALSE_POSITIVE | Filter B flanking: Chr11:12,121,778–12,123,978 (±500 slop) |
| rice_G281 | Chr10:13,500,285 | FALSE_POSITIVE | Filter A: 96% host, 452 bp gap < 500 bp |
| rice_G281 | Chr1:29,159,516 | FALSE_POSITIVE | Filter D: 92% construct |
| cucumber_line225 | LKUO03001451.1:6,501 | CANDIDATE | 3.2% host, 5,403 bp gap, element annotated |
| cucumber_line225 | LKUO03002166.1:547,872 | FALSE_POSITIVE | Filter D: 87% construct |

---

## Implementation note: flanking_hit slop encoding

`_flanking_overlaps_site` in `verdict.py` uses a simple boundary check
(`lo <= site_pos <= hi`).  The caller (`generate_report`) must expand the
raw construct-flanking coordinates by `CONSTRUCT_FLANK_SLOP = 500 bp` before
storing them in `FilterEvidence.flanking_hit`.  This keeps verdict.py
side-effect-free (no constants import) while preserving the ±500 bp slop
matching semantics.
