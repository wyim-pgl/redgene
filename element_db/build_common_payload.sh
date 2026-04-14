#!/usr/bin/env bash
# Fetch common transgene payload sequences from NCBI Entrez.
# Re-run only when element_db/common_payload_manifest.tsv changes.
set -euo pipefail
cd "$(dirname "$0")"

MANIFEST=common_payload_manifest.tsv
OUT=common_payload.fa

TMPOUT=$(mktemp)
trap 'rm -f "$TMPOUT"' EXIT

while IFS=$'\t' read -r acc purpose notes; do
    [[ "$acc" == "accession" ]] && continue
    [[ -z "$acc" ]] && continue
    header="${purpose}|${acc}"
    echo "Fetching $acc ($purpose)..." >&2
    efetch -db nuccore -id "$acc" -format fasta \
        | awk -v h=">${header}" 'NR==1{print h; next}{print}' \
        >> "$TMPOUT"
    sleep 0.4
done < "$MANIFEST"

# Only publish on success. If efetch failed anywhere above, set -e aborted
# the script and this line never runs, leaving the prior $OUT intact.
mv "$TMPOUT" "$OUT"

echo "Wrote $OUT with $(grep -c '^>' "$OUT") sequences"
