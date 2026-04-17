#!/usr/bin/env python3
"""Concat input FASTAs and append a |src=<tag> suffix to every header.

Used by element_db/Makefile to prepare the pre-dedup input for cd-hit-est.
The suffix is parsed downstream by scripts/s05_insert_assembly.py so
hits can be tiered by source (Task T5).
"""
import sys
from pathlib import Path

# Filename -> source tag.  Keep in sync with scripts/s05_insert_assembly.py
# `_SRC_TIER` keys (Issue #10 M-1).  The path sync points are:
#   * element_db/Makefile            — SRCS list is the input side
#   * this SRCMAP                    — filename -> tag mapping
#   * scripts/s05_insert_assembly.py — `_SRC_TIER` consumer + priority
#   * scripts/build_element_mask_bed.sh — reads the final DB via `$DB`
# Changing any of these requires updating the others.
SRCMAP = {
    "common_payload.fa":   "payload",
    "payload_cds.fa":      "payload",
    "gmo_combined_db.fa":  "element_db",
    "cas9_sgrna.fa":       "element_db",
    "euginius_missing.fa": "element_db",
    # Issue #5: 5th source tag — EU CRL GMOMETHODS amplicons (82 seqs
    # after hand-curation; v1.1 rebuild will pull all ~85 raw rows and
    # let cd-hit-est handle the 3 duplicates the legacy FASTA dropped).
    "crl_amplicons.fa":    "crl",
}


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: tag_and_concat.py <fa> [<fa>...]", file=sys.stderr)
        return 2
    for fa in argv[1:]:
        name = Path(fa).name
        tag = SRCMAP.get(name)
        if tag is None:
            print(f"error: unknown source file {name}; add to SRCMAP",
                  file=sys.stderr)
            return 2
        for line in Path(fa).read_text().splitlines():
            if line.startswith(">"):
                # Preserve full original header; append tag suffix.
                sys.stdout.write(f"{line}|src={tag}\n")
            else:
                sys.stdout.write(f"{line}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
