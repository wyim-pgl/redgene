#!/usr/bin/env python3
"""fetch_amplicons.py

Use EUginius primer data to fetch real amplicon sequences from NCBI.

Strategy:
  1. Methods with clean amplicons (no N) → use directly
  2. Methods with NNN amplicons → BLAST forward primer via qblast(blastn-short)
     to find source accession, then fetch amplicon region

Output:
  - euginius_amplicons.fa  (FASTA, format matching gmo_elements_db.fa)
  - euginius_amplicons.tsv (summary table)

Usage:
  python fetch_amplicons.py \
      --input euginius_primers.tsv \
      --output euginius_amplicons.fa \
      --email wyim@unr.edu
"""

import argparse
import csv
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def clean_seq(seq: str) -> str:
    """Remove non-alpha chars, uppercase."""
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def has_degenerate(seq: str) -> bool:
    """Check if sequence has degenerate bases (not ACGT or N)."""
    return bool(re.search(r"[^ACGTN]", seq))


def resolve_degenerate(seq: str) -> str:
    """Replace degenerate IUPAC bases with first matching base for BLAST."""
    iupac = {
        "R": "A", "Y": "C", "S": "G", "W": "A", "K": "G", "M": "A",
        "B": "C", "D": "A", "H": "A", "V": "A",
    }
    return "".join(iupac.get(c, c) for c in seq)


def revcomp(seq: str) -> str:
    """Reverse complement."""
    comp = str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN")
    return seq.translate(comp)[::-1]


def find_primer_in_seq(primer: str, target: str) -> Optional[int]:
    """Find primer position in target allowing degenerate bases."""
    primer_resolved = resolve_degenerate(primer)
    # Try forward
    pos = target.find(primer_resolved)
    if pos >= 0:
        return pos
    # Try with mismatches (up to 2)
    plen = len(primer_resolved)
    for i in range(len(target) - plen + 1):
        mismatches = sum(1 for a, b in zip(primer_resolved, target[i:i+plen]) if a != b)
        if mismatches <= 2:
            return i
    return None


def blast_primer(primer_seq: str, email: str, database: str = "nt") -> Optional[Dict]:
    """BLAST a primer sequence against NCBI and return top hit info."""
    Entrez.email = email
    query = resolve_degenerate(primer_seq)

    try:
        result_handle = NCBIWWW.qblast(
            "blastn", database, query,
            word_size=7,
            expect=1000,
            megablast=False,
            hitlist_size=3,
        )
        blast_records = NCBIXML.parse(result_handle)
        record = next(blast_records)

        if not record.alignments:
            return None

        alignment = record.alignments[0]
        hsp = alignment.hsps[0]

        # Extract accession from title
        title = alignment.title
        acc_match = re.search(r"gi\|\d+\|[^|]*\|([^|]+)\|", title)
        if not acc_match:
            acc_match = re.search(r"([A-Z]{1,2}_?\d+\.\d+)", title)
        accession = acc_match.group(1).rstrip("|") if acc_match else title.split()[0]

        return {
            "accession": accession,
            "title": title,
            "sbjct_start": hsp.sbjct_start,
            "sbjct_end": hsp.sbjct_end,
            "identity": hsp.identities,
            "align_len": hsp.align_length,
            "strand": "+" if hsp.sbjct_start < hsp.sbjct_end else "-",
        }
    except Exception as e:
        print(f"  BLAST error: {e}", file=sys.stderr)
        return None


def fetch_region(accession: str, start: int, end: int, email: str) -> Optional[str]:
    """Fetch a sequence region from NCBI."""
    Entrez.email = email
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=accession,
            rettype="fasta", retmode="text",
            seq_start=min(start, end),
            seq_stop=max(start, end),
        )
        record = SeqIO.read(handle, "fasta")
        handle.close()
        seq = str(record.seq).upper()
        if start > end:  # reverse strand
            seq = revcomp(seq)
        return seq
    except Exception as e:
        print(f"  Fetch error for {accession}:{start}-{end}: {e}", file=sys.stderr)
        return None


def safe_header(text: str) -> str:
    """Clean text for FASTA header."""
    if not text:
        return "NA"
    text = text.replace("|", "_").replace("/", "_").replace("\\", "_")
    text = text.replace("'", "").replace('"', "")
    text = re.sub(r"\s+", "_", text.strip())
    return text[:80] or "NA"


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    parser = argparse.ArgumentParser(description="Fetch amplicon sequences from NCBI using EUginius primers")
    parser.add_argument("--input", default="euginius_primers.tsv")
    parser.add_argument("--output", default="euginius_amplicons.fa")
    parser.add_argument("--summary", default="euginius_amplicons.tsv")
    parser.add_argument("--email", required=True)
    parser.add_argument("--max-blast", type=int, default=0, help="Max BLAST queries (0=unlimited)")
    parser.add_argument("--skip-blast", action="store_true", help="Only use clean amplicons, skip BLAST")
    args = parser.parse_args()

    with open(args.input) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    print(f"[INFO] Loaded {len(rows)} methods", file=sys.stderr)

    records = []
    summary = []
    blast_count = 0
    blast_cache: Dict[str, Optional[Dict]] = {}  # primer_seq -> blast result

    for i, row in enumerate(rows):
        code = row["method_code"]
        scope = row["target_scope"]
        target = row["target_name"]
        fwd = clean_seq(row["forward_primer_seq"])
        rev = clean_seq(row["reverse_primer_seq"])
        amp_seq = clean_seq(row.get("amplicon_seq", ""))
        amp_size = row.get("amplicon_size", "")

        print(f"[{i+1}/{len(rows)}] {code} ({scope})", file=sys.stderr)

        source = ""
        final_seq = ""
        accession = "EUginius"
        status = ""

        # Strategy 1: clean amplicon available
        if amp_seq and "N" not in amp_seq and len(amp_seq) > 20:
            final_seq = amp_seq
            source = "euginius_amplicon"
            status = "OK_CLEAN"
            print(f"  -> Clean amplicon: {len(final_seq)} bp", file=sys.stderr)

        # Strategy 2: BLAST forward primer
        elif not args.skip_blast and fwd and len(fwd) >= 15:
            if args.max_blast > 0 and blast_count >= args.max_blast:
                status = "SKIPPED_BLAST_LIMIT"
                print(f"  -> Skipped (blast limit reached)", file=sys.stderr)
            else:
                # Check cache
                if fwd in blast_cache:
                    hit = blast_cache[fwd]
                    print(f"  -> Using cached BLAST result", file=sys.stderr)
                else:
                    print(f"  -> BLASTing fwd primer ({len(fwd)} bp)...", file=sys.stderr)
                    hit = blast_primer(fwd, args.email)
                    blast_cache[fwd] = hit
                    blast_count += 1
                    time.sleep(0.5)  # rate limit

                if hit:
                    accession = hit["accession"]
                    # Calculate amplicon region
                    try:
                        size = int(amp_size) if amp_size and amp_size.isdigit() else 500
                    except ValueError:
                        size = 500

                    if hit["strand"] == "+":
                        start = hit["sbjct_start"]
                        end = start + size - 1
                    else:
                        end = hit["sbjct_start"]
                        start = end - size + 1
                        if start < 1:
                            start = 1

                    fetched = fetch_region(accession, start, end, args.email)
                    if fetched:
                        final_seq = fetched
                        source = f"NCBI:{accession}:{start}-{end}"
                        status = "OK_BLAST"
                        print(f"  -> Fetched: {len(final_seq)} bp from {accession}", file=sys.stderr)
                        time.sleep(0.3)
                    else:
                        status = "FAILED_FETCH"
                        print(f"  -> Failed to fetch from {accession}", file=sys.stderr)
                else:
                    status = "NO_BLAST_HIT"
                    print(f"  -> No BLAST hit", file=sys.stderr)
        else:
            status = "SKIPPED"
            print(f"  -> Skipped", file=sys.stderr)

        if final_seq:
            header = (
                f">{safe_header(scope)}|{safe_header(code)}|"
                f"{safe_header(target)}|{accession}|"
                f"amplicon_{amp_size}bp|euginius_v1"
            )
            records.append((header, final_seq))

        summary.append({
            "method_code": code,
            "target_scope": scope,
            "target_name": target,
            "amplicon_size": amp_size,
            "actual_length": len(final_seq) if final_seq else 0,
            "source": source,
            "status": status,
            "accession": accession,
        })

    # Write FASTA
    with open(args.output, "w") as f:
        for header, seq in records:
            f.write(header + "\n")
            f.write(wrap_fasta(seq) + "\n")
    print(f"\n[INFO] FASTA: {args.output} ({len(records)} sequences)", file=sys.stderr)

    # Write summary
    with open(args.summary, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    print(f"[INFO] Summary: {args.summary}", file=sys.stderr)

    # Stats
    ok = sum(1 for s in summary if s["status"].startswith("OK"))
    fail = sum(1 for s in summary if "FAIL" in s["status"] or "NO_" in s["status"])
    skip = sum(1 for s in summary if "SKIP" in s["status"])
    print(f"\n[INFO] OK: {ok}, Failed: {fail}, Skipped: {skip}", file=sys.stderr)
    print(f"[INFO] BLAST queries made: {blast_count}", file=sys.stderr)


if __name__ == "__main__":
    main()
