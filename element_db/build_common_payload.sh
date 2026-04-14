#!/usr/bin/env bash
# Fetch common transgene payload sequences from NCBI Entrez.
# Re-run only when element_db/common_payload_manifest.tsv changes.
set -euo pipefail
cd "$(dirname "$0")"

MANIFEST=common_payload_manifest.tsv
OUT=common_payload.fa
TMP=$(mktemp -d)

: > "$OUT"
while IFS=$'\t' read -r acc purpose notes; do
    [[ "$acc" == "accession" ]] && continue
    [[ -z "$acc" ]] && continue
    # Derive a stable sequence name: purpose|accession
    header="${purpose}|${acc}"
    echo "Fetching $acc ($purpose)..." >&2
    efetch -db nuccore -id "$acc" -format fasta \
        | awk -v h=">${header}" 'NR==1{print h; next}{print}' \
        >> "$OUT"
    sleep 0.4  # be nice to NCBI
done < "$MANIFEST"

echo >&2
echo "Wrote $OUT with $(grep -c '^>' "$OUT") sequences" >&2
rm -rf "$TMP"
