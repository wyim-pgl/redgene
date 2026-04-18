# `common_payload_manifest.tsv` schema

The manifest drives `element_db/build_common_payload.sh`, which fetches a
small curated set of canonical transgene-payload sequences from NCBI Entrez
and writes `element_db/common_payload.fa`. The pipeline wires the resulting
FASTA in as an always-on `--common-payload-db` for s05 (see `run_pipeline.py`).

## Columns (tab-separated, 5 columns)

| # | Column      | Required? | Description                                                                                          |
|---|-------------|-----------|------------------------------------------------------------------------------------------------------|
| 1 | `accession` | yes       | NCBI nucleotide accession (with version), e.g. `X17220.1`.                                           |
| 2 | `purpose`   | yes       | Short canonical element name used in the FASTA header, e.g. `bar`, `nptII`, `P-CaMV35S`, `T-nos`.    |
| 3 | `seq_start` | optional  | 1-based inclusive fetch start. Set together with `seq_stop` to slice multi-element accessions.       |
| 4 | `seq_stop`  | optional  | 1-based inclusive fetch end. Required iff `seq_start` is set.                                        |
| 5 | `notes`     | optional  | Free-text rationale / source / BUG reference. Ignored by the build script; present for audit trail. |

- The first row must be the header `accession<TAB>purpose<TAB>seq_start<TAB>seq_stop<TAB>notes`.
- Blank rows and rows whose first column is `accession` are skipped.
- `seq_start` / `seq_stop` are mandatory for any accession that contains more
  than one payload element (e.g. `V00087.1` carries both P-nos and T-nos);
  omitting the range silently mis-annotates reads (BUG-9). Use a 1-based
  inclusive coordinate system.

## FASTA header format

The build script emits one of two header shapes depending on whether the
range columns are populated:

- Whole-accession fetch:
  `>${purpose}|${accession}`
- Ranged fetch:
  `>${purpose}|${accession}:${seq_start}-${seq_stop}`

Downstream s05 annotation relies on the `|` separator to split the `purpose`
tag from the provenance accession, so neither field may contain `|`.

## Atomic write contract

`build_common_payload.sh` stages all `efetch` output into a temp file located
in the same directory as `common_payload.fa` (`mktemp -p "$(dirname "$OUT")" ...`).
On `set -euo pipefail` failure the script aborts before the final `mv`,
leaving the previous `common_payload.fa` intact. On success the `mv` is an
atomic rename because both paths share a filesystem — downstream consumers
therefore never observe a half-written FASTA.

## Re-running

Only re-run the build when the manifest changes. Requires the NCBI Entrez
CLI (`efetch`) on `$PATH`. A 0.4 s inter-request sleep keeps the loop under
NCBI's per-IP rate limit.
