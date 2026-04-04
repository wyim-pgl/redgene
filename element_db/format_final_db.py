#!/usr/bin/env python3
"""format_final_db.py

Merge gmo_elements_db.fa (curated elements), euginius_fullseq.fa (full source sequences),
and euginius_amplicons.fa (fallback amplicons) into a unified FASTA database.

Priority order:
  1. Curated elements (gmo_elements_db.fa)
  2. Full sequences from EUginius BLAST (euginius_fullseq.fa)
  3. Amplicon-only fallback (euginius_amplicons.fa) - used when full seq not available

Header format (matching gmo_elements_db.fa):
  >class|name|organism_or_target|accession|description|version

Output:
  - gmo_combined_db.fa    (merged FASTA)
  - gmo_combined_db.tsv   (summary table)
"""

import csv
import re
import sys
from pathlib import Path
from Bio import SeqIO


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--elements-fa", default="gmo_elements_db.fa")
    parser.add_argument("--fullseq-fa", default="euginius_fullseq.fa")
    parser.add_argument("--amplicons-fa", default="euginius_amplicons.fa")
    parser.add_argument("--primers-tsv", default="euginius_primers.tsv")
    parser.add_argument("--output-fa", default="gmo_combined_db.fa")
    parser.add_argument("--output-tsv", default="gmo_combined_db.tsv")
    args = parser.parse_args()

    records = []
    seen_seqs = set()
    seen_methods = set()  # track method codes to avoid amplicon+fullseq duplication

    # Load curated elements
    if Path(args.elements_fa).exists():
        for rec in SeqIO.parse(args.elements_fa, "fasta"):
            seq = str(rec.seq).upper()
            records.append({
                "header": rec.id,
                "seq": seq,
                "source": "curated_element",
                "length": len(seq),
            })
            seen_seqs.add(seq)
        print(f"[INFO] Loaded {len(records)} curated elements", file=sys.stderr)

    # Load primer data for metadata
    primer_data = {}
    if Path(args.primers_tsv).exists():
        with open(args.primers_tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                primer_data[row["method_code"]] = row

    # Load full sequences (priority over amplicons)
    fullseq_count = 0
    if Path(args.fullseq_fa).exists():
        for rec in SeqIO.parse(args.fullseq_fa, "fasta"):
            seq = str(rec.seq).upper()
            if seq in seen_seqs:
                continue
            seen_seqs.add(seq)
            # Extract method code from header to track
            parts = rec.id.split("|")
            if len(parts) > 1:
                seen_methods.add(parts[1])
            records.append({
                "header": rec.id,
                "seq": seq,
                "source": "euginius_fullseq",
                "length": len(seq),
            })
            fullseq_count += 1
        print(f"[INFO] Added {fullseq_count} unique EUginius full sequences",
              file=sys.stderr)

    # Load EUginius amplicons (fallback for methods without full sequence)
    amp_count = 0
    if Path(args.amplicons_fa).exists():
        for rec in SeqIO.parse(args.amplicons_fa, "fasta"):
            seq = str(rec.seq).upper()
            if seq in seen_seqs:
                continue
            # Skip if we already have full sequence for this method
            parts = rec.id.split("|")
            method_code = parts[1] if len(parts) > 1 else ""
            if method_code in seen_methods:
                continue
            seen_seqs.add(seq)
            records.append({
                "header": rec.id,
                "seq": seq,
                "source": "euginius_amplicon",
                "length": len(seq),
            })
            amp_count += 1
    print(f"[INFO] Added {amp_count} unique EUginius amplicons (fallback)",
          file=sys.stderr)

    # Write combined FASTA
    with open(args.output_fa, "w") as f:
        for rec in records:
            f.write(f">{rec['header']}\n")
            f.write(wrap_fasta(rec["seq"]) + "\n")
    print(f"[INFO] Combined FASTA: {args.output_fa} ({len(records)} sequences)",
          file=sys.stderr)

    # Write summary TSV
    fieldnames = ["header_id", "source", "length_bp", "sequence_type"]
    with open(args.output_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in records:
            parts = rec["header"].split("|")
            writer.writerow({
                "header_id": parts[1] if len(parts) > 1 else rec["header"],
                "source": rec["source"],
                "length_bp": rec["length"],
                "sequence_type": parts[0] if parts else "unknown",
            })
    print(f"[INFO] Summary TSV: {args.output_tsv}", file=sys.stderr)

    # Stats
    total_bp = sum(r["length"] for r in records)
    print(f"\n[INFO] SUMMARY", file=sys.stderr)
    print(f"  Curated elements: {len(records) - fullseq_count - amp_count}",
          file=sys.stderr)
    print(f"  EUginius full sequences: {fullseq_count}", file=sys.stderr)
    print(f"  EUginius amplicons (fallback): {amp_count}", file=sys.stderr)
    print(f"  Total sequences: {len(records)}", file=sys.stderr)
    print(f"  Total bp: {total_bp:,}", file=sys.stderr)


if __name__ == "__main__":
    main()
