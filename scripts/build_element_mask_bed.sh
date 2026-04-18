#!/bin/bash
# scripts/build_element_mask_bed.sh
# Generate BED of regions in host genome that match element_db >=98% identity >=200bp.
# Usage: build_element_mask_bed.sh <host_fasta> <out_bed> [min_id] [min_len]
set -euo pipefail

HOST="$1"
OUT_BED="$2"
MIN_ID="${3:-98}"
MIN_LEN="${4:-200}"

ELEMENT_DB="${ELEMENT_DB:-element_db/gmo_combined_db_v2.fa}"
THREADS="${THREADS:-8}"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# makeblastdb if index missing
if [ ! -f "${HOST}.nsq" ]; then
    echo "[build_element_mask_bed] makeblastdb on $HOST (one-time)"
    makeblastdb -in "$HOST" -dbtype nucl
fi

echo "[build_element_mask_bed] megablast $ELEMENT_DB vs $HOST (threads=$THREADS)"
blastn -task megablast -query "$ELEMENT_DB" -db "$HOST" \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
    -evalue 1e-10 -num_threads "$THREADS" -out "$WORK/hits.tsv"

echo "[build_element_mask_bed] filter >=${MIN_ID}% id >=${MIN_LEN}bp + sort"
awk -v mi="$MIN_ID" -v ml="$MIN_LEN" -F'\t' \
    '$3 >= mi && $4 >= ml {
        if ($7 < $8) print $2 "\t" ($7-1) "\t" $8 "\t" $1 "\t" $3;
        else         print $2 "\t" ($8-1) "\t" $7 "\t" $1 "\t" $3
    }' "$WORK/hits.tsv" \
    | sort -k1,1 -k2,2n \
    > "$WORK/raw.bed"

echo "[build_element_mask_bed] bedtools merge (use ';' separator for element names)"
if [ -s "$WORK/raw.bed" ]; then
    bedtools merge -i "$WORK/raw.bed" -c 4,5 -o distinct,max -delim ";" \
        > "$OUT_BED"
else
    : > "$OUT_BED"
fi

echo "[build_element_mask_bed] done. $(wc -l < "$OUT_BED") regions -> $OUT_BED"
