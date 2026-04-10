#!/usr/bin/env python3
"""Filter UniVec FASTA to keep only vectors, removing adapters/primers/linkers.

Reads NCBI UniVec (gnl|uv|ACCESSION|description format), removes sequences
whose headers match adapter/primer/linker keywords, and writes a filtered
FASTA with 'univec|' prefix on each sequence name.

Two modes:
  --mode all      Keep all vectors (remove adapters/primers/linkers only)
  --mode plant    Keep only plant expression/binary vectors (for GMO detection)

Usage:
    python scripts/filter_univec.py \
        --input db/UniVec.fa \
        --output db/univec_vectors.fa

    python scripts/filter_univec.py \
        --input db/UniVec.fa \
        --output db/univec_plant_vectors.fa \
        --mode plant
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


# Keywords indicating adapter/primer/linker sequences to REMOVE
REMOVE_KEYWORDS = [
    r"\badapter\b",
    r"\bprimer\b",
    r"\blinker\b",
    r"\bPCR\b",
    r"\bIllumina\b",
    r"\bIon\s*Torrent\b",
    r"\b454\b",
    r"\bSOLiD\b",
    r"\bbarcode\b",
    r"\bindex\b",
]

# Compiled pattern (case-insensitive)
_REMOVE_RE = re.compile("|".join(REMOVE_KEYWORDS), re.IGNORECASE)

# Plant expression/binary vector keywords (KEEP in plant mode)
PLANT_KEYWORDS = [
    r"\bbinary\b",
    r"\bpCAMBIA\b",
    r"\bpPZP\d",
    r"\bpGreen\b",
    r"\bpSoup\b",
    r"\bpBI\d",
    r"\bpBIN\d",
    r"\bBin\s*19\b",
    r"\bplant\s+transform",
    r"\bT-DNA\b",
    r"\bGUS\s+gene\s+fusion\b",
    r"\bpART\b",
    r"\bpORE\b",
    r"\bpMDC\b",
    r"\bpEarley\b",
    r"\bpGWB\b",
    r"\bpCLEAN\b",
    r"\bpMAA\b",
]

_PLANT_RE = re.compile("|".join(PLANT_KEYWORDS), re.IGNORECASE)


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Parse FASTA into list of (header, sequence) tuples."""
    records: list[tuple[str, str]] = []
    header = ""
    seq_parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header:
        records.append((header, "".join(seq_parts)))
    return records


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter UniVec to vectors only (remove adapters/primers/linkers)")
    parser.add_argument("--input", required=True, help="Input UniVec FASTA")
    parser.add_argument("--output", required=True, help="Output filtered FASTA")
    parser.add_argument("--mode", choices=["all", "plant"], default="all",
                        help="'all': keep all vectors; 'plant': keep only plant binary/expression vectors")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    records = parse_fasta(input_path)
    total = len(records)
    print(f"Total sequences in UniVec: {total}", file=sys.stderr)

    kept: list[tuple[str, str]] = []
    removed_categories: dict[str, int] = {}

    for header, seq in records:
        # Step 1: always remove adapters/primers/linkers
        match = _REMOVE_RE.search(header)
        if match:
            keyword = match.group(0).lower()
            removed_categories[keyword] = removed_categories.get(keyword, 0) + 1
            continue

        # Step 2: in plant mode, only keep plant binary/expression vectors
        if args.mode == "plant":
            if not _PLANT_RE.search(header):
                removed_categories["non-plant-vector"] = removed_categories.get("non-plant-vector", 0) + 1
                continue

        kept.append((header, seq))

    removed = total - len(kept)

    # Write filtered output with univec| prefix
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        for header, seq in kept:
            # Prefix the sequence name with univec|
            fh.write(f">univec|{header}\n")
            # Write sequence in 70-char lines
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")

    print(f"Mode:    {args.mode}", file=sys.stderr)
    print(f"Kept:    {len(kept)}", file=sys.stderr)
    print(f"Removed: {removed}", file=sys.stderr)
    print(f"Categories removed:", file=sys.stderr)
    for cat, count in sorted(removed_categories.items(), key=lambda x: -x[1]):
        print(f"  {cat}: {count}", file=sys.stderr)
    print(f"Output written to: {output_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
