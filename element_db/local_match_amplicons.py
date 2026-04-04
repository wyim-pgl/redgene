#!/usr/bin/env python3
"""local_match_amplicons.py

Match EUginius primers against known sequences to extract amplicons.
Uses local sequence matching (no BLAST needed) against:
  1. gmo_elements_db.fa (our curated DB)
  2. NCBI sequences fetched by accession from element.list

For remaining unmatched, fetches from NCBI using primer BLAST.

Output: appends to euginius_amplicons.fa
"""

import csv
import re
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Seq import Seq


def clean_seq(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def resolve_degenerate(seq: str) -> str:
    iupac = {
        "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
        "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
        "H": "[ACT]", "V": "[ACG]",
    }
    return "".join(iupac.get(c, c) for c in seq.upper())


def find_primer_regex(primer: str, target: str) -> int:
    """Find primer in target using regex for degenerate bases. Returns start position or -1."""
    pattern = resolve_degenerate(primer)
    try:
        m = re.search(pattern, target, re.IGNORECASE)
        if m:
            return m.start()
    except re.error:
        pass
    return -1


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def safe_header(text: str) -> str:
    if not text:
        return "NA"
    text = text.replace("|", "_").replace("/", "_").replace("\\", "_")
    text = text.replace("'", "").replace('"', "")
    text = re.sub(r"\s+", "_", text.strip())
    return text[:80] or "NA"


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def load_local_db(fasta_path: str) -> dict:
    """Load FASTA sequences into dict keyed by record id."""
    seqs = {}
    if not Path(fasta_path).exists():
        return seqs
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seqs[rec.id] = str(rec.seq).upper()
    return seqs


def try_local_match(fwd: str, rev: str, amp_size: int, db_seqs: dict):
    """Try to find amplicon in local database sequences."""
    rev_rc = revcomp(rev)

    for seq_id, seq in db_seqs.items():
        fwd_pos = find_primer_regex(fwd, seq)
        if fwd_pos < 0:
            continue

        rev_pos = find_primer_regex(rev_rc, seq)
        if rev_pos < 0:
            # Try finding reverse primer on reverse strand
            rev_pos_fwd = find_primer_regex(rev, seq)
            if rev_pos_fwd >= 0 and rev_pos_fwd > fwd_pos:
                amplicon = seq[fwd_pos:rev_pos_fwd + len(rev)]
                if amp_size == 0 or abs(len(amplicon) - amp_size) < amp_size * 0.3:
                    return amplicon, seq_id
            continue

        if rev_pos > fwd_pos:
            amplicon = seq[fwd_pos:rev_pos + len(rev_rc)]
            if amp_size == 0 or abs(len(amplicon) - amp_size) < amp_size * 0.3:
                return amplicon, seq_id

    return None, None


def fetch_accession_seq(accession: str, email: str) -> str:
    """Fetch full sequence from NCBI by accession."""
    Entrez.email = email
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession,
                               rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq).upper()
    except Exception as e:
        print(f"  Fetch error {accession}: {e}", file=sys.stderr)
        return ""


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="euginius_primers.tsv")
    parser.add_argument("--local-db", default="gmo_elements_db.fa")
    parser.add_argument("--element-list", default="element.list")
    parser.add_argument("--output", default="euginius_amplicons.fa")
    parser.add_argument("--summary", default="euginius_amplicons.tsv")
    parser.add_argument("--email", required=True)
    args = parser.parse_args()

    # Load primer data
    with open(args.input) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    # Filter to only methods that need BLAST (no clean amplicon)
    need_match = []
    for r in rows:
        amp = clean_seq(r.get("amplicon_seq", ""))
        if not amp or "N" in amp or len(amp) <= 20:
            need_match.append(r)

    print(f"[INFO] {len(need_match)} methods need amplicon sequences", file=sys.stderr)

    # Load local DB
    db_seqs = load_local_db(args.local_db)
    print(f"[INFO] Local DB: {len(db_seqs)} sequences", file=sys.stderr)

    # Load known accessions from element.list for additional sequences
    known_accessions = set()
    if Path(args.element_list).exists():
        with open(args.element_list) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 5 and parts[4] not in ("MANUAL", "accession"):
                    known_accessions.add(parts[4])
    print(f"[INFO] Known accessions: {len(known_accessions)}", file=sys.stderr)

    # Fetch and add known accession sequences to local DB
    acc_cache = {}
    for acc in sorted(known_accessions):
        if acc in acc_cache:
            continue
        print(f"  Fetching accession {acc}...", file=sys.stderr)
        seq = fetch_accession_seq(acc, args.email)
        if seq:
            acc_cache[acc] = seq
            db_seqs[f"ncbi|{acc}"] = seq
        time.sleep(0.4)
    print(f"[INFO] Total local+NCBI sequences: {len(db_seqs)}", file=sys.stderr)

    # Read existing amplicons file to avoid duplicates
    existing_codes = set()
    if Path(args.output).exists():
        for rec in SeqIO.parse(args.output, "fasta"):
            parts = rec.id.split("|")
            if len(parts) >= 2:
                existing_codes.add(parts[1])
    print(f"[INFO] Existing amplicons: {len(existing_codes)}", file=sys.stderr)

    # Match primers
    matched = []
    unmatched = []

    for i, row in enumerate(need_match):
        code = row["method_code"]
        code_safe = safe_header(code)

        if code_safe in existing_codes:
            continue

        fwd = clean_seq(row["forward_primer_seq"])
        rev = clean_seq(row["reverse_primer_seq"])
        amp_size_str = row.get("amplicon_size", "")
        amp_size = int(amp_size_str) if amp_size_str and amp_size_str.isdigit() else 0

        print(f"[{i+1}/{len(need_match)}] {code}", file=sys.stderr, end="")

        amplicon, source_id = try_local_match(fwd, rev, amp_size, db_seqs)

        if amplicon:
            matched.append((row, amplicon, source_id))
            print(f" -> MATCHED ({len(amplicon)} bp from {source_id[:40]})", file=sys.stderr)
        else:
            unmatched.append(row)
            print(f" -> no match", file=sys.stderr)

    print(f"\n[INFO] Matched: {len(matched)}, Unmatched: {len(unmatched)}", file=sys.stderr)

    # Append matched amplicons to FASTA
    with open(args.output, "a") as f:
        for row, amplicon, source_id in matched:
            scope = safe_header(row["target_scope"])
            code = safe_header(row["method_code"])
            target = safe_header(row["target_name"])
            amp_size = row.get("amplicon_size", "NA")

            # Extract accession from source_id
            acc = source_id.split("|")[1] if "|" in source_id else source_id.split("|")[0]

            header = f">{scope}|{code}|{target}|{acc}|amplicon_{amp_size}bp|euginius_v1"
            f.write(header + "\n")
            f.write(wrap_fasta(amplicon) + "\n")

    # Update summary
    new_summary = []
    if Path(args.summary).exists():
        with open(args.summary) as f:
            reader = csv.DictReader(f, delimiter="\t")
            new_summary = list(reader)

    for row, amplicon, source_id in matched:
        code = row["method_code"]
        # Update existing entry or add new
        found = False
        for s in new_summary:
            if s["method_code"] == code:
                s["actual_length"] = str(len(amplicon))
                s["source"] = f"local_match:{source_id}"
                s["status"] = "OK_LOCAL"
                found = True
                break
        if not found:
            new_summary.append({
                "method_code": code,
                "target_scope": row["target_scope"],
                "target_name": row["target_name"],
                "amplicon_size": row.get("amplicon_size", ""),
                "actual_length": str(len(amplicon)),
                "source": f"local_match:{source_id}",
                "status": "OK_LOCAL",
                "accession": source_id,
            })

    if new_summary:
        with open(args.summary, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=new_summary[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(new_summary)

    # Report unmatched
    if unmatched:
        print(f"\n[INFO] Unmatched methods ({len(unmatched)}):", file=sys.stderr)
        for r in unmatched[:20]:
            print(f"  {r['method_code']:40s} {r['target_scope']:20s} {r['target_name'][:50]}", file=sys.stderr)
        if len(unmatched) > 20:
            print(f"  ... and {len(unmatched)-20} more", file=sys.stderr)

    total = len(existing_codes) + len(matched)
    print(f"\n[INFO] Total amplicons in {args.output}: {total}", file=sys.stderr)


if __name__ == "__main__":
    main()
