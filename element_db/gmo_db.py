#!/usr/bin/env python3
"""build_gmo_element_db.py

Build a curated GMO transgenic element reference FASTA database
by downloading sequences from NCBI GenBank and combining with
manually curated sequences.

Sources:
    - Won Yim curated element list (2025)
    - Debode et al. 2019 (Sci Rep 9:15595)
    - EUginius database element names
    - JRC GMOMETHODS / GMO-Amplicons element names
    - Published literature (accessions in master list)

Usage:
    python build_gmo_element_db.py \
        --master element_master_list.tsv \
        --manual manual_sequences.fa \
        --output gmo_elements_db.fa \
        --email your@email.com \
        --log build_log.tsv

Output:
    1. gmo_elements_db.fa       - Combined FASTA ready for bwa index
    2. build_log.tsv             - Download log with status for each element
    3. element_summary.tsv       - Summary table for supplementary material

Requirements:
    pip install biopython
"""

import argparse
import csv
import sys
import time
import os
from datetime import datetime

try:
    from Bio import Entrez, SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("ERROR: Biopython is required. Install with:")
    print("  pip install biopython")
    sys.exit(1)


def parse_master_list(filepath):
    """Parse the TSV master list of GMO elements."""
    elements = []
    with open(filepath, "r") as f:
        reader = csv.DictReader(
            (row for row in f if not row.startswith("#")),
            delimiter="\t"
        )
        for row in reader:
            elements.append(row)
    return elements


def parse_manual_fasta(filepath):
    """Parse manually curated FASTA sequences into a dict keyed by header_id."""
    manual = {}
    if not os.path.exists(filepath):
        return manual
    for record in SeqIO.parse(filepath, "fasta"):
        # Extract header_id from the FASTA header
        # Format: class|header_id|organism|...
        parts = record.id.split("|")
        if len(parts) >= 2:
            header_id = parts[1]
            manual[header_id] = record
    return manual


def fetch_from_ncbi(accession, start, end, email, retries=3):
    """
    Fetch a sequence region from NCBI GenBank.

    Parameters:
        accession: GenBank accession (e.g., V00141.1)
        start: 1-based start coordinate (or "FULL")
        end: 1-based end coordinate (or "FULL")
        email: email for NCBI Entrez
        retries: number of retry attempts

    Returns:
        SeqRecord or None
    """
    Entrez.email = email

    for attempt in range(retries):
        try:
            if start == "FULL" or end == "FULL":
                # Fetch entire sequence
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=accession,
                    rettype="fasta",
                    retmode="text"
                )
            else:
                # Fetch specific region (Entrez uses 1-based coordinates)
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=accession,
                    rettype="fasta",
                    retmode="text",
                    seq_start=int(start),
                    seq_stop=int(end)
                )

            record = SeqIO.read(handle, "fasta")
            handle.close()
            return record

        except Exception as e:
            print(f"  Attempt {attempt + 1}/{retries} failed for {accession}: {e}")
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))  # Exponential backoff

    return None


def build_header(element):
    """Build standardized FASTA header from element dict."""
    parts = [
        element["element_class"],
        element["header_id"],
        element["organism"],
        element["accession"],
        element["boundary_rule"].replace(" ", "_").replace(";", "")[:80],
        "canonical_v1"
    ]
    return "|".join(parts)


def main():
    parser = argparse.ArgumentParser(
        description="Build GMO transgenic element reference FASTA database"
    )
    parser.add_argument(
        "--master", required=True,
        help="Path to element_master_list.tsv"
    )
    parser.add_argument(
        "--manual", required=True,
        help="Path to manual_sequences.fa"
    )
    parser.add_argument(
        "--output", default="gmo_elements_db.fa",
        help="Output FASTA file (default: gmo_elements_db.fa)"
    )
    parser.add_argument(
        "--email", required=True,
        help="Email for NCBI Entrez (required by NCBI)"
    )
    parser.add_argument(
        "--log", default="build_log.tsv",
        help="Build log file (default: build_log.tsv)"
    )
    parser.add_argument(
        "--summary", default="element_summary.tsv",
        help="Summary table for supplementary material"
    )
    args = parser.parse_args()

    print(f"{'=' * 60}")
    print(f"GMO Transgenic Element Database Builder")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'=' * 60}")

    # Parse inputs
    elements = parse_master_list(args.master)
    manual_seqs = parse_manual_fasta(args.manual)
    print(f"\nLoaded {len(elements)} elements from master list")
    print(f"Loaded {len(manual_seqs)} manual sequences")

    # Process each element
    records = []
    log_entries = []
    summary_entries = []

    for i, elem in enumerate(elements):
        header_id = elem["header_id"]
        accession = elem["accession"]
        start = elem["seq_start"]
        end = elem["seq_end"]
        element_name = elem["element_name"]
        element_class = elem["element_class"]

        print(f"\n[{i+1}/{len(elements)}] {element_class}: {element_name} ({header_id})")

        if accession == "MANUAL":
            # Use manually curated sequence
            if header_id in manual_seqs:
                manual_rec = manual_seqs[header_id]
                # Rebuild with standardized header
                fasta_header = build_header(elem)
                new_record = SeqRecord(
                    manual_rec.seq,
                    id=fasta_header,
                    description=""
                )
                records.append(new_record)
                seq_len = len(manual_rec.seq)
                status = "OK_MANUAL"
                print(f"  -> Manual sequence: {seq_len} bp")
            else:
                status = "MISSING_MANUAL"
                seq_len = 0
                print(f"  -> WARNING: Manual sequence not found for {header_id}")

        else:
            # Fetch from NCBI
            print(f"  -> Fetching {accession} [{start}:{end}] from NCBI...")
            record = fetch_from_ncbi(accession, start, end, args.email)

            if record is not None:
                fasta_header = build_header(elem)
                new_record = SeqRecord(
                    record.seq,
                    id=fasta_header,
                    description=""
                )
                records.append(new_record)
                seq_len = len(record.seq)
                status = "OK_NCBI"
                print(f"  -> Downloaded: {seq_len} bp")

                # NCBI rate limit: max 3 requests per second
                time.sleep(0.4)
            else:
                status = "FAILED_NCBI"
                seq_len = 0
                print(f"  -> ERROR: Failed to download {accession}")

        # Log entry
        log_entries.append({
            "element_class": element_class,
            "element_name": element_name,
            "header_id": header_id,
            "accession": accession,
            "coordinates": f"{start}-{end}",
            "seq_length": seq_len,
            "status": status,
            "timestamp": datetime.now().isoformat()
        })

        # Summary entry
        summary_entries.append({
            "Element class": element_class,
            "Element name": element_name,
            "FASTA header ID": header_id,
            "Source organism": elem["organism"],
            "GenBank accession": accession if accession != "MANUAL" else "curated",
            "Coordinates": f"{start}-{end}" if accession != "MANUAL" else "full",
            "Length (bp)": seq_len,
            "Boundary rule": elem["boundary_rule"],
            "Source note": elem["source_note"]
        })

    # Write output FASTA
    print(f"\n{'=' * 60}")
    print(f"Writing output files...")

    with open(args.output, "w") as f:
        SeqIO.write(records, f, "fasta")
    print(f"  FASTA: {args.output} ({len(records)} sequences)")

    # Write build log
    if log_entries:
        with open(args.log, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=log_entries[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(log_entries)
        print(f"  Log: {args.log}")

    # Write summary table
    if summary_entries:
        with open(args.summary, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=summary_entries[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_entries)
        print(f"  Summary: {args.summary}")

    # Report
    ok_count = sum(1 for e in log_entries if e["status"].startswith("OK"))
    fail_count = sum(1 for e in log_entries if "FAILED" in e["status"] or "MISSING" in e["status"])
    total_bp = sum(e["seq_length"] for e in log_entries)

    print(f"\n{'=' * 60}")
    print(f"SUMMARY")
    print(f"  Total elements: {len(elements)}")
    print(f"  Successfully retrieved: {ok_count}")
    print(f"  Failed/missing: {fail_count}")
    print(f"  Total sequence: {total_bp:,} bp")
    print(f"{'=' * 60}")

    if fail_count > 0:
        print(f"\nWARNING: {fail_count} elements could not be retrieved.")
        print("Check build_log.tsv for details.")
        print("Add missing sequences to manual_sequences.fa and re-run.")

    # Print bwa index command
    print(f"\nNext steps:")
    print(f"  # Index for BWA")
    print(f"  bwa index {args.output}")
    print(f"")
    print(f"  # Or index for minimap2")
    print(f"  minimap2 -d {args.output.replace('.fa', '.mmi')} {args.output}")


if __name__ == "__main__":
    main()
