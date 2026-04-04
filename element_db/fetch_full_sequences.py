#!/usr/bin/env python3
"""fetch_full_sequences.py

For each EUginius amplicon (OK_CLEAN), find the source accession via BLAST
and download the FULL accession sequence (not just the amplicon region).

Strategy:
  1. Try local match: check if amplicon is a subsequence of curated DB entries
  2. BLAST amplicon against NCBI nt to find source accession
  3. Fetch full sequence for that accession

Output:
  - euginius_fullseq.fa   (FASTA with full source sequences)
  - euginius_fullseq.tsv  (summary table)

Usage:
  python fetch_full_sequences.py \
      --amplicons euginius_amplicons.fa \
      --amplicons-tsv euginius_amplicons.tsv \
      --local-db gmo_elements_db.fa \
      --email wyim@unr.edu
"""

import argparse
import csv
import re
import sys
import time
from pathlib import Path
from typing import Dict, Optional, Tuple

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq


def clean_seq(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq).upper()


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


def load_fasta_dict(path: str) -> Dict[str, str]:
    """Load FASTA as {header_id: sequence}."""
    seqs = {}
    if not Path(path).exists():
        return seqs
    for rec in SeqIO.parse(path, "fasta"):
        seqs[rec.description] = str(rec.seq).upper()
    return seqs


def find_in_local_db(amplicon: str, db_seqs: Dict[str, str]) -> Optional[str]:
    """Check if amplicon is a subsequence of any local DB entry. Returns header."""
    amp_rc = revcomp(amplicon)
    for header, seq in db_seqs.items():
        if amplicon in seq or amp_rc in seq:
            return header
    return None


VECTOR_KEYWORDS = [
    "vector", "plasmid", "cloning", "pbi", "pcambia", "pbluescript",
    "puc", "pgem", "pet-", "pbr322", "pgreen", "pbin", "psoup",
    "binary vector", "shuttle vector", "expression vector", "t-dna",
    "transformation vector", "construct", "recombinant plasmid",
]

# Max length for a single gene/element (skip if larger = likely vector/genome)
MAX_ELEMENT_LENGTH = 8000


def is_vector(title: str, length: int) -> bool:
    """Check if an accession is likely a vector/plasmid rather than a gene."""
    title_lower = title.lower()
    for kw in VECTOR_KEYWORDS:
        if kw in title_lower:
            return True
    if length > MAX_ELEMENT_LENGTH:
        return True
    return False


def blast_sequence(seq: str, email: str, database: str = "nt",
                   program: str = "blastn") -> Optional[Dict]:
    """BLAST a sequence against NCBI and return top hit info.
    Skips vector/plasmid hits; prefers gene/element accessions."""
    Entrez.email = email
    try:
        result_handle = NCBIWWW.qblast(
            program, database, seq,
            megablast=True,
            hitlist_size=10,
            expect=1e-10,
        )
        blast_records = NCBIXML.parse(result_handle)
        record = next(blast_records)

        if not record.alignments:
            return None

        # Collect all good hits, preferring non-vector ones
        candidates = []
        for alignment in record.alignments[:10]:
            hsp = alignment.hsps[0]
            title = alignment.title

            # Extract accession
            acc_match = re.search(r"([A-Z]{1,2}_?\d+\.\d+)", title)
            accession = acc_match.group(1) if acc_match else title.split()[0]

            identity_pct = hsp.identities / hsp.align_length * 100
            if identity_pct < 95:
                continue

            info = {
                "accession": accession,
                "title": title,
                "identity_pct": identity_pct,
                "align_len": hsp.align_length,
                "sbjct_start": hsp.sbjct_start,
                "sbjct_end": hsp.sbjct_end,
                "sbjct_length": alignment.length,
                "is_vector": is_vector(title, alignment.length),
            }
            candidates.append(info)

        if not candidates:
            return None

        # Prefer non-vector hits
        non_vector = [c for c in candidates if not c["is_vector"]]
        if non_vector:
            # Among non-vectors, prefer shorter (more likely individual gene)
            non_vector.sort(key=lambda x: x["sbjct_length"])
            return non_vector[0]

        # All hits are vectors - return None (will fall back to amplicon)
        print(f"  -> All BLAST hits are vectors/plasmids, skipping",
              file=sys.stderr)
        return None

    except Exception as e:
        print(f"  BLAST error: {e}", file=sys.stderr)
        return None


def fetch_full_sequence(accession: str, email: str) -> Optional[Tuple[str, str]]:
    """Fetch full sequence from NCBI. Returns (sequence, description)."""
    Entrez.email = email
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=accession,
            rettype="fasta", retmode="text",
        )
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq).upper(), record.description
    except Exception as e:
        print(f"  Fetch error for {accession}: {e}", file=sys.stderr)
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Fetch full source sequences for EUginius amplicons")
    parser.add_argument("--amplicons", default="euginius_amplicons.fa")
    parser.add_argument("--amplicons-tsv", default="euginius_amplicons.tsv")
    parser.add_argument("--local-db", default="gmo_elements_db.fa")
    parser.add_argument("--output-fa", default="euginius_fullseq.fa")
    parser.add_argument("--output-tsv", default="euginius_fullseq.tsv")
    parser.add_argument("--email", required=True)
    parser.add_argument("--max-blast", type=int, default=0,
                        help="Max BLAST queries (0=unlimited)")
    parser.add_argument("--skip-local", action="store_true",
                        help="Skip already-local-matched entries")
    args = parser.parse_args()

    # Load amplicons
    amplicons = load_fasta_dict(args.amplicons)
    print(f"[INFO] Loaded {len(amplicons)} amplicons", file=sys.stderr)

    # Load amplicon summary to know status
    amp_status = {}
    if Path(args.amplicons_tsv).exists():
        with open(args.amplicons_tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                amp_status[row["method_code"]] = row

    # Load local DB
    local_db = load_fasta_dict(args.local_db)
    print(f"[INFO] Local DB: {len(local_db)} sequences", file=sys.stderr)

    results = []
    blast_count = 0
    seen_accessions = {}  # accession -> (seq, desc) cache

    for i, (header, amp_seq) in enumerate(amplicons.items()):
        parts = header.split("|")
        scope = parts[0] if len(parts) > 0 else "unknown"
        method_code = parts[1] if len(parts) > 1 else header
        target = parts[2] if len(parts) > 2 else "unknown"

        # Restore / from _ in method_code for lookup
        method_code_orig = method_code.replace("_", "/", method_code.count("_"))
        # Try exact match first
        status_row = amp_status.get(method_code_orig)
        if not status_row:
            # Try with safe_header version
            for k in amp_status:
                if safe_header(k) == method_code:
                    status_row = amp_status[k]
                    method_code_orig = k
                    break

        print(f"\n[{i+1}/{len(amplicons)}] {method_code} ({scope}, {len(amp_seq)}bp)",
              file=sys.stderr)

        source_type = ""
        full_seq = ""
        full_desc = ""
        accession = ""

        # Strategy 1: Check if amplicon is in local curated DB
        local_match = find_in_local_db(amp_seq, local_db)
        if local_match:
            full_seq = local_db[local_match]
            # Extract accession from header
            lparts = local_match.split("|")
            accession = lparts[3] if len(lparts) > 3 else "local"
            full_desc = local_match
            source_type = "local_fullseq"
            print(f"  -> Local match: {local_match[:60]}... ({len(full_seq)}bp)",
                  file=sys.stderr)

        # Strategy 2: BLAST against NCBI
        elif args.max_blast == 0 or blast_count < args.max_blast:
            print(f"  -> BLASTing amplicon ({len(amp_seq)}bp)...", file=sys.stderr)
            hit = blast_sequence(amp_seq, args.email)
            blast_count += 1
            time.sleep(1.0)  # rate limit

            if hit:
                accession = hit["accession"]
                print(f"  -> BLAST hit: {accession} "
                      f"({hit['identity_pct']:.1f}% identity, "
                      f"subject len={hit['sbjct_length']})",
                      file=sys.stderr)

                # Fetch full sequence (use cache if available)
                if accession in seen_accessions:
                    full_seq, full_desc = seen_accessions[accession]
                    print(f"  -> Using cached full sequence ({len(full_seq)}bp)",
                          file=sys.stderr)
                else:
                    result = fetch_full_sequence(accession, args.email)
                    if result:
                        full_seq, full_desc = result
                        seen_accessions[accession] = (full_seq, full_desc)
                        print(f"  -> Fetched full sequence: {len(full_seq)}bp",
                              file=sys.stderr)
                        time.sleep(0.5)
                    else:
                        print(f"  -> Failed to fetch full sequence", file=sys.stderr)
                source_type = "ncbi_blast_fullseq"
            else:
                print(f"  -> No BLAST hit", file=sys.stderr)
                source_type = "NO_HIT"
        else:
            source_type = "SKIPPED_BLAST_LIMIT"
            print(f"  -> Skipped (blast limit)", file=sys.stderr)

        results.append({
            "method_code": method_code_orig or method_code,
            "scope": scope,
            "target": target,
            "amplicon_len": len(amp_seq),
            "fullseq_len": len(full_seq) if full_seq else 0,
            "accession": accession,
            "source_type": source_type,
            "description": full_desc[:100] if full_desc else "",
            "full_seq": full_seq,
        })

    # Write full sequence FASTA (deduplicate by accession)
    written_accessions = set()
    with open(args.output_fa, "w") as f:
        for r in results:
            if not r["full_seq"]:
                continue
            # Deduplicate by accession+sequence
            dedup_key = r["accession"] + "|" + r["scope"]
            if dedup_key in written_accessions:
                continue
            written_accessions.add(dedup_key)

            header = (
                f">{r['scope']}|{safe_header(r['method_code'])}|"
                f"{safe_header(r['target'])}|{r['accession']}|"
                f"fullseq_{r['fullseq_len']}bp|euginius_v1"
            )
            f.write(header + "\n")
            f.write(wrap_fasta(r["full_seq"]) + "\n")

    n_written = len(written_accessions)
    print(f"\n[INFO] Full sequences FASTA: {args.output_fa} "
          f"({n_written} sequences)", file=sys.stderr)

    # Write summary TSV
    with open(args.output_tsv, "w", newline="") as f:
        fieldnames = ["method_code", "scope", "target", "amplicon_len",
                      "fullseq_len", "accession", "source_type", "description"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in results:
            row = {k: v for k, v in r.items() if k != "full_seq"}
            writer.writerow(row)
    print(f"[INFO] Summary: {args.output_tsv}", file=sys.stderr)

    # Stats
    local_ok = sum(1 for r in results if r["source_type"] == "local_fullseq")
    blast_ok = sum(1 for r in results if r["source_type"] == "ncbi_blast_fullseq")
    no_hit = sum(1 for r in results if r["source_type"] == "NO_HIT")
    skipped = sum(1 for r in results if "SKIP" in r["source_type"])
    print(f"\n[INFO] SUMMARY:", file=sys.stderr)
    print(f"  Local match (full seq in curated DB): {local_ok}", file=sys.stderr)
    print(f"  BLAST + NCBI fetch: {blast_ok}", file=sys.stderr)
    print(f"  No BLAST hit: {no_hit}", file=sys.stderr)
    print(f"  Skipped: {skipped}", file=sys.stderr)
    print(f"  Total BLAST queries: {blast_count}", file=sys.stderr)


if __name__ == "__main__":
    main()
