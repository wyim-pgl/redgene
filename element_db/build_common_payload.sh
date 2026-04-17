#!/usr/bin/env bash
# Fetch common transgene payload sequences from NCBI Entrez.
# Re-run only when element_db/common_payload_manifest.tsv changes.
#
# Manifest columns (tab-separated):
#   accession  purpose  seq_start  seq_stop  notes
#
# seq_start / seq_stop are optional. If both are provided, fetch is restricted
# to that subregion (1-based, inclusive). This is REQUIRED for multi-element
# accessions such as V00087.1 (nos operon containing both P-nos and T-nos) to
# avoid silent mis-annotation (BUG-9).
set -euo pipefail
cd "$(dirname "$0")"

MANIFEST=common_payload_manifest.tsv
OUT=common_payload.fa

TMPOUT=$(mktemp)
trap 'rm -f "$TMPOUT"' EXIT

while IFS=$'\t' read -r acc purpose seq_start seq_stop notes; do
    [[ "$acc" == "accession" ]] && continue
    [[ -z "$acc" ]] && continue
    if [[ -n "$seq_start" && -n "$seq_stop" ]]; then
        header=">${purpose}|${acc}:${seq_start}-${seq_stop}"
        echo "Fetching $acc:$seq_start-$seq_stop ($purpose)..." >&2
        efetch -db nuccore -id "$acc" -format fasta \
            -seq_start "$seq_start" -seq_stop "$seq_stop" \
            | awk -v h="$header" 'NR==1{print h; next}{print}' \
            >> "$TMPOUT"
    else
        header=">${purpose}|${acc}"
        echo "Fetching $acc ($purpose)..." >&2
        efetch -db nuccore -id "$acc" -format fasta \
            | awk -v h="$header" 'NR==1{print h; next}{print}' \
            >> "$TMPOUT"
    fi
    sleep 0.4
done < "$MANIFEST"

# Only publish on success. If efetch failed anywhere above, set -e aborted
# the script and this line never runs, leaving the prior $OUT intact.
mv "$TMPOUT" "$OUT"

echo "Wrote $OUT with $(grep -c '^>' "$OUT") sequences"
