#!/usr/bin/env bash
# element_db/build_crl_amplicons.sh
#
# Issue #5 (v1.0 scope: script only; DB rebuild deferred to v1.1).
#
# Converts the committed element_db/crl_amplicons_raw.tsv (82+ EU CRL
# GMOMETHODS amplicon records — construct/element/event/taxon-specific)
# into element_db/crl_amplicons.fa using the T5 four-way source-tag
# convention.  The downstream `tag_and_concat.py` SRCMAP maps this file
# onto the `crl` source tag (the 5th tag alongside payload / element_db /
# plant_endogenous / corn_border).
#
# Input TSV schema (6 cols, tab-separated, header row):
#   method_id  type  target  accession  length  sequence
#
# Output FASTA header (6 pipe-separated fields, no spaces in any field):
#   >${type}|${method_id}|${target_slug}|${accession_or_fallback}|amplicon_${length}bp|crl_v1
# where:
#   * target_slug  = target with ' ' -> '_' (accession-safe slugging)
#   * accession_or_fallback = accession when non-"none" / non-empty, else
#                             the literal "CRL-GMOMETHODS" fallback string
#
# Atomic-write contract is identical to build_common_payload.sh:
#   1. stage output into $(mktemp -p "$(dirname "$OUT")")
#   2. mv "$TMPOUT" "$OUT" only on success so a partial write can never
#      clobber the previously-good FASTA.

set -euo pipefail
cd "$(dirname "$0")"

RAW=crl_amplicons_raw.tsv
OUT=crl_amplicons.fa

if [ ! -s "$RAW" ]; then
    echo "error: $RAW missing or empty" >&2
    exit 2
fi

TMPOUT=$(mktemp -p "$(dirname "$OUT")" crl_amplicons.XXXXXX.fa)
trap 'rm -f "$TMPOUT"' EXIT

# Column indices (1-based): 1=method_id 2=type 3=target 4=accession 5=length 6=sequence
awk -F '\t' '
    BEGIN { OFS=""; }
    NR==1 { next }                    # skip header
    NF < 6 { next }                   # malformed rows
    {
        method_id = $1
        type      = $2
        target    = $3
        acc       = $4
        len       = $5
        seq       = $6

        # Accession-safe slug: replace any whitespace with `_`.
        gsub(/[[:space:]]+/, "_", target)

        # Fallback for rows where the CRL lab did not publish an accession
        # (literal string "none" in the TSV).  Matches the convention used
        # in the hand-curated v0 element_db/crl_amplicons.fa.
        if (acc == "" || acc == "none") acc = "CRL-GMOMETHODS"

        print ">", type, "|", method_id, "|", target, "|", acc, "|amplicon_", len, "bp|crl_v1"
        # Hard-wrap the sequence at 80 cols to match the sibling build
        # scripts (efetch also emits 80-col FASTA by default).
        while (length(seq) > 0) {
            print substr(seq, 1, 80)
            seq = substr(seq, 81)
        }
    }
' "$RAW" > "$TMPOUT"

mv "$TMPOUT" "$OUT"
echo "Wrote $OUT with $(grep -c '^>' "$OUT") sequences"
