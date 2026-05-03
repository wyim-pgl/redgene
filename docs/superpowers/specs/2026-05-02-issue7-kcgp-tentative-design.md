# Issue #7 — KCGP Nomenclature (tentative_v0) Design

**Date**: 2026-05-02
**Status**: Spec — pending implementation plan
**GitHub issue**: [wyim-pgl/redgene#7](https://github.com/wyim-pgl/redgene/issues/7)
**Implements**: R-6 (KCGP nomenclature, partial)
**Spec version**: `tentative_v0` — official KCGP nomenclature spec from team lead is still pending. This design produces ID strings under a clearly-marked placeholder format that can be regenerated under the official spec without schema migration.

## Goal

Provide a deterministic, reproducible site-ID scheme so insertion sites can be referenced in regulatory contexts. The official KCGP spec is undelivered, so this implementation:

1. Emits IDs under a placeholder format (`KCGP-{HOST3}-{EVENT}-{YY}-{SITE4}`).
2. Always carries a `spec_version = "tentative_v0"` tag in every output.
3. Surfaces the placeholder warning in the PDF report so downstream readers cannot mistake it for the official ID.
4. Lives in standalone files so the regenerate-under-official-spec path is one CLI invocation away.

## Non-Goals

- Implementing the official KCGP spec. That waits on the spec document from the team lead.
- Wiring KCGP IDs into Chain-of-Custody (`coc_logger.py`). CoC logs are step-level, not site-level — there is no natural integration point. Defer to v1 once the spec is known.
- Adding KCGP ID columns to the PDF sites table. Cover-page surfacing is enough for the tentative phase.
- Auto-running ID assignment from `run_pipeline.py`. Standalone post-hoc script keeps blast radius small under a tentative spec.

## Architecture

```
config.yaml (one new field per sample)
  └─ samples.{sample}.event: "G281"

scripts/util/kcgp_id.py (new, ~80 lines)
  ├─ HOST_CODES = {host_ref_basename → 3-letter code} (5 hosts)
  ├─ build_kcgp_id(host_ref, event, year, site_ord, spec="tentative_v0") → str
  ├─ parse_kcgp_id(kcgp_id: str) → dict
  └─ _normalize_host_path(host_ref_path) → str

scripts/util/assign_kcgp_ids.py (new, ~120 lines)
  └─ CLI: --sample-dir results/{sample} --config config.yaml [--year YYYY]
     Reads insertion_*_report.txt → sorts → emits kcgp_mapping.tsv

scripts/reports/insertion_pdf.py (one new loader, one cover-page edit)
  ├─ load_kcgp_mapping(sample_dir) → dict | None
  └─ _page_cover() → prepends KCGP ID lines when mapping exists

tests/test_kcgp_id.py (new)
tests/test_assign_kcgp_ids.py (new)
tests/test_insertion_pdf_kcgp_cover.py (new)
```

### Invariants

- `parse_kcgp_id(build_kcgp_id(...))` is the identity for round-trip.
- Same inputs always produce the same ID (deterministic, no random seed).
- `kcgp_mapping.tsv` absence is non-fatal everywhere — the PDF and the rest of the pipeline must continue working as before.
- The string `tentative_v0` appears in every mapping TSV row and on every PDF cover that surfaces a KCGP ID, so external readers cannot mistake the placeholder for the official spec.
- Only verdict `CANDIDATE` sites receive a 4-digit ordinal. `FALSE_POSITIVE` and `UNKNOWN` sites get `-` in the `kcgp_id` column to keep the audit trail without polluting the ordinal space. Verdict labels match exactly what the existing `insertion_*_report.txt` files write today.

## ID Format

`KCGP-{HOST3}-{EVENT}-{YY}-{SITE4}`

| Field | Width | Source | Validation |
|-------|-------|--------|------------|
| `HOST3` | 3 chars | `HOST_CODES[normalize(host_reference)]` | KeyError if host not registered |
| `EVENT` | variable | `config.yaml` `samples.{name}.event` | `^[A-Za-z0-9_]+$`; underscore allowed because the field separator is `-` |
| `YY` | 2 digits | `--year` arg or `datetime.now().year % 100` | `00`-`99` |
| `SITE4` | 4 digits | per-CANDIDATE ordinal after sorting | `0001`-`9999`; ValueError on overflow |

`HOST_CODES`:
```python
HOST_CODES = {
    "Osativa_323_v7.0":   "OSA",   # rice
    "SLM_r2.0.pmol":      "SLE",   # tomato
    "CucSat_B10v3":       "CSA",   # cucumber
    "Zm_B73_v5":          "ZMA",   # corn
    "Gmax_v4.0":          "GMA",   # soybean
}
```

`_normalize_host_path("db/Osativa_323_v7.0.fa")` strips the directory and a single trailing `.fa`/`.fasta` extension, returning the dict key.

Example IDs:
- `rice_G281` → `KCGP-OSA-G281-26-0001`
- `tomato_Cas9_A2_3` (event=`A2_3`) → `KCGP-SLE-A2_3-26-0001`
- `cucumber_line225` → `KCGP-CSA-line225-26-0001`

## `assign_kcgp_ids.py` Behavior

**Input discovery**: glob `{sample_dir}/s05_insert_assembly/insertion_*_report.txt`.

**Per-report parse**:
- Filename → `host_chr` and `pos_5p` (regex `insertion_(.+)_(\d+)_report\.txt`).
- Body → first `Verdict:` line. Whitespace-trim the value.

**Sort key**: `(verdict_priority, host_chr_natural, pos_5p)`.
- `verdict_priority`: `CANDIDATE` first (priority 0), then `FALSE_POSITIVE` (1), `UNKNOWN` (2), anything else (3).
- `host_chr_natural`: split into alphanumeric runs and compare numerically where applicable so `Chr2 < Chr10`.
- `pos_5p`: integer ascending.

**Ordinal assignment**: enumerate sorted `CANDIDATE` rows from 1, assign `SITE4 = f"{n:04d}"`. All other verdicts: `kcgp_id = "-"`.

**Output**: `{sample_dir}/kcgp_mapping.tsv` with header
`kcgp_id\tinternal_site\tverdict\thost_chr\tpos_5p\tspec_version`. One row per discovered report. Header always written, even when zero rows.

### Edge cases

| Case | Behavior |
|------|----------|
| `s05_insert_assembly/` absent | Exit 0; write header-only `kcgp_mapping.tsv`; log warning to stderr. |
| Zero `insertion_*_report.txt` files | Exit 0; write header-only `kcgp_mapping.tsv`. |
| Sample missing from `config.yaml` | Exit 1 with explanatory message. |
| Sample present but no `event` field | Exit 1 with `add 'event:' field under samples.{name}` hint. |
| `host_reference` not in `HOST_CODES` | Exit 1 listing the registered keys. |
| `> 9999` CANDIDATE sites | Exit 1; the 4-digit field cannot represent the result. (Not expected at any realistic biological scale.) |
| Filename does not match the regex | Log warning, skip the file, continue. |
| `Verdict:` line missing from a report | Treat as verdict `UNKNOWN`, log warning. |

### CLI

```bash
python scripts/util/assign_kcgp_ids.py \
  --sample-dir results/rice_G281 \
  --config config.yaml \
  [--year 2026]
```

`--year` defaults to `datetime.now().year`. Sample name is derived from `--sample-dir` basename. Output path is hard-coded to `{sample_dir}/kcgp_mapping.tsv`.

## PDF Cover Integration

`scripts/reports/insertion_pdf.py`:

**New loader** (alongside `load_audit_header`, `load_copynumber`, etc.):
```python
def load_kcgp_mapping(sample_dir: Path) -> dict[str, Any] | None:
    """Return {'kcgp_id', 'spec_version'} for the lowest-ordinal CANDIDATE row, or None."""
    path = Path(sample_dir) / "kcgp_mapping.tsv"
    if not path.exists():
        return None
    try:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("verdict") == "CANDIDATE" and row.get("kcgp_id", "-") != "-":
                    return {
                        "kcgp_id": row["kcgp_id"],
                        "spec_version": row["spec_version"],
                    }
    except (OSError, csv.Error):
        return None
    return None
```

**`_page_cover` change**: prepend two lines to the existing `lines` list when a mapping is found:
```python
kcgp = load_kcgp_mapping(sample_dir)
if kcgp:
    lines = [
        f"⚠ KCGP ID ({kcgp['spec_version']} — pending official spec):",
        f"    {kcgp['kcgp_id']}",
        "",
    ] + lines
```

When the mapping is absent the cover renders identically to today, protecting the existing scaffold tests.

## Error Handling

Match the soft-failure pattern of the other PDF loaders: `(OSError, csv.Error)` returns `None`, no traceback. The CLI scripts (`assign_kcgp_ids.py`) emit explicit error messages and non-zero exit codes for misconfiguration so the operator knows to fix `config.yaml` or `HOST_CODES`.

## Testing

The 255-test baseline must remain green. New test files:

### `tests/test_kcgp_id.py`

| Test | Assertion |
|------|-----------|
| `test_build_basic` | `build_kcgp_id("db/Osativa_323_v7.0.fa", "G281", 2026, 1)` returns `"KCGP-OSA-G281-26-0001"` |
| `test_build_all_5_hosts` | All 5 entries in `HOST_CODES` produce well-formed IDs |
| `test_build_unknown_host_raises` | `db/unknown.fa` raises `KeyError` |
| `test_round_trip` | `parse_kcgp_id(build_kcgp_id(...))` returns the original fields, for 10 inputs |
| `test_parse_invalid_format` | `KCGP-XX` (too few fields) raises `ValueError` |
| `test_event_with_underscore` | `event="A2_3"` produces `"KCGP-SLE-A2_3-26-0001"` and round-trips |
| `test_normalize_host_path` | strips `db/` prefix and `.fa`/`.fasta` suffix |
| `test_site_overflow` | `site_ord=10000` raises `ValueError` |

### `tests/test_assign_kcgp_ids.py`

| Test | Fixture | Assertion |
|------|---------|-----------|
| `test_no_s05_dir` | empty sample_dir | exit 0; mapping TSV is header-only |
| `test_no_reports` | empty `s05_insert_assembly/` | exit 0; header-only |
| `test_assigns_ordinal_to_candidate_only` | 1 CANDIDATE + 1 FALSE_POSITIVE + 1 UNKNOWN reports | CANDIDATE row has `KCGP-...-0001`; FALSE_POSITIVE/UNKNOWN rows have `-` |
| `test_sort_order_natural` | reports for `Chr2:100`, `Chr10:200`, `Chr1:50` (all CANDIDATE) | ordinals `0001`/`0002`/`0003` go to Chr1, Chr2, Chr10 in that order |
| `test_missing_event_in_config` | config without `event` field | exit 1, message mentions `event:` |
| `test_unknown_host_in_config` | config with `db/foo.fa` | exit 1, lists known hosts |

### `tests/test_insertion_pdf_kcgp_cover.py`

| Test | Fixture | Assertion |
|------|---------|-----------|
| `test_cover_without_mapping` | sample_dir with no `kcgp_mapping.tsv` | `_page_cover` body identical to current behavior (regression guard via monkeypatched `_page_text` capture) |
| `test_cover_with_mapping` | mapping TSV with 1 CANDIDATE row | captured `lines` list contains `"KCGP ID (tentative_v0"` and the actual ID string |
| `test_cover_skips_dash_rows` | mapping TSV where all rows have `kcgp_id == "-"` | cover body unchanged (no KCGP block) |
| `test_load_kcgp_mapping_absent` | no file | returns `None` |
| `test_load_kcgp_mapping_corrupt` | malformed TSV bytes | returns `None`, no exception |

## Acceptance Criteria

After implementation, run on the three validation samples already used for Issue #6:

1. Add `event` field to `config.yaml` for `rice_G281` (`G281`), `tomato_Cas9_A2_3` (`A2_3`), `cucumber_line225` (`line225`).
2. Run `python scripts/util/assign_kcgp_ids.py --sample-dir results/{sample} --config config.yaml` for each.
3. Verify each `results/{sample}/kcgp_mapping.tsv` exists and contains at least one `KCGP-...-0001` row (the lone CANDIDATE for each sample).
4. Regenerate PDFs: `python run_pipeline.py --sample {sample} --steps 8`.
5. Confirm the cover page on each PDF shows `⚠ KCGP ID (tentative_v0 — pending official spec): KCGP-XYZ-...-0001`.
6. Comment on Issue #7 with the 3 mapping TSV summaries and the migration path note ("regenerate via assign_kcgp_ids.py once spec arrives").

Issue #7 stays open after this lands. The official spec is still pending; closing requires the v1 implementation.

## Migration Path (when official spec arrives)

1. Replace `HOST_CODES` mapping and `build_kcgp_id` body to match the official spec.
2. Bump default `spec_version` from `"tentative_v0"` to `"v1"`.
3. Re-run `assign_kcgp_ids.py` against every sample. The `kcgp_mapping.tsv` schema is stable; only the `kcgp_id` column values and `spec_version` change.
4. Regenerate PDFs.
5. Drop the `⚠` warning prefix from `_page_cover`.
6. Wire into CoC if the official spec defines a site-level field there.

No schema migration is required because the mapping TSV is a flat per-sample artifact regenerated from the canonical `insertion_*_report.txt` files. The internal pipeline never reads `kcgp_id` back; only humans and the PDF do.

## Out of Scope

- KCGP IDs in CoC log entries (no natural site-level field; tracked under Issue #8 follow-up if needed).
- KCGP column in PDF sites table (cover-page surfacing is enough for tentative).
- KCGP ID in JSON audit header.
- Cross-sample uniqueness checks. The official spec may require a global registry; for tentative_v0 IDs are only unique within a sample.
- Step-8 PDF regeneration as a side effect of `assign_kcgp_ids.py`. Operators run the steps separately.

## Open Questions

None at spec time. Implementation may surface natural-sort edge cases for unusual chromosome names (contigs like `B10v3_unplaced_contig_1234`); the implementer can fall back to lexical sort with a stderr warning if natural sort produces unexpected ordering.
