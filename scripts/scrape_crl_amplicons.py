#!/usr/bin/env python3
"""Scrape amplicon sequences from EU CRL GMOMETHODS database.

Downloads all method pages and extracts amplicon sequences into FASTA format.
"""

import re
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path


BASE_URL = "https://gmo-crl.jrc.ec.europa.eu/gmomethods"


def fetch_page(url: str, retries: int = 3) -> str:
    """Fetch URL content with retries."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "User-Agent": "Mozilla/5.0 (RedGene pipeline, academic research)"
            })
            with urllib.request.urlopen(req, timeout=30) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, TimeoutError) as e:
            print(f"  Retry {attempt+1}/{retries}: {e}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    return ""


def get_method_ids() -> list[str]:
    """Collect all method IDs from paginated listing."""
    method_ids = []
    page = 1
    while True:
        url = f"{BASE_URL}/?page={page}"
        print(f"Fetching listing page {page}...", file=sys.stderr)
        html = fetch_page(url)
        if not html:
            break

        # Match method IDs like QL-CON-00-001, QL-ELE-xxx, QT-TAX-xxx, QL-EVE-xxx
        ids = re.findall(r'href="/gmomethods/show/([A-Z0-9][A-Za-z0-9_-]+)/"', html)
        if not ids:
            break

        method_ids.extend(ids)
        print(f"  Found {len(ids)} methods on page {page}", file=sys.stderr)

        # Check if there's a next page
        if f"page={page + 1}" not in html:
            break
        page += 1
        time.sleep(0.5)

    return list(dict.fromkeys(method_ids))  # deduplicate preserving order


def extract_amplicon(method_id: str, html: str) -> dict | None:
    """Extract amplicon sequence and metadata from a method page."""
    # Look for the amplicon sequence - typically a block of ACGT characters
    # Try multiple patterns

    # Pattern 1: Large block of nucleotides (amplicon)
    seq_matches = re.findall(r'(?<![a-zA-Z])([ACGTacgtNnRYSWKMBDHV\s]{50,}?)(?![a-zA-Z])', html)

    # Filter to only valid DNA sequences
    amplicon_seq = None
    for match in seq_matches:
        cleaned = re.sub(r'\s+', '', match).upper()
        # Must be mostly ACGT
        acgt_count = sum(1 for c in cleaned if c in "ACGTN")
        if len(cleaned) >= 50 and acgt_count / len(cleaned) > 0.95:
            amplicon_seq = cleaned
            break

    if not amplicon_seq:
        return None

    # Extract GenBank accession
    genbank = ""
    gb_match = re.search(r'(?:GenBank|Accession|NCBI)[^>]*?([A-Z]{1,2}\d{5,8})', html)
    if gb_match:
        genbank = gb_match.group(1)

    # Extract target description
    desc = ""
    # Try to get the method name/description
    title_match = re.search(r'<h[12][^>]*>([^<]+)</h[12]>', html)
    if title_match:
        desc = title_match.group(1).strip()

    # Extract amplicon length
    length = len(amplicon_seq)

    return {
        "method_id": method_id,
        "sequence": amplicon_seq,
        "length": length,
        "genbank": genbank,
        "description": desc,
    }


def main():
    outdir = Path("element_db")
    outdir.mkdir(exist_ok=True)

    print("=== Collecting method IDs ===", file=sys.stderr)
    method_ids = get_method_ids()
    print(f"Found {len(method_ids)} methods total", file=sys.stderr)

    amplicons = []
    failed = []

    print("\n=== Fetching amplicon sequences ===", file=sys.stderr)
    for i, mid in enumerate(method_ids):
        url = f"{BASE_URL}/show/{mid}/"
        print(f"[{i+1}/{len(method_ids)}] {mid}...", file=sys.stderr, end=" ")

        html = fetch_page(url)
        if not html:
            print("FETCH FAILED", file=sys.stderr)
            failed.append(mid)
            continue

        result = extract_amplicon(mid, html)
        if result:
            amplicons.append(result)
            print(f"OK ({result['length']}bp)", file=sys.stderr)
        else:
            print("NO SEQUENCE FOUND", file=sys.stderr)
            failed.append(mid)

        time.sleep(0.3)  # be polite

    # Write FASTA
    outfile = outdir / "crl_amplicons.fa"
    with open(outfile, "w") as f:
        for amp in amplicons:
            header = f">{amp['method_id']}|{amp['genbank']}|{amp['length']}bp|crl_gmomethods"
            if amp['description']:
                header += f"|{amp['description']}"
            f.write(f"{header}\n")
            # Wrap sequence at 80 chars
            seq = amp['sequence']
            for j in range(0, len(seq), 80):
                f.write(seq[j:j+80] + "\n")

    print(f"\n=== Results ===", file=sys.stderr)
    print(f"Total methods: {len(method_ids)}", file=sys.stderr)
    print(f"Amplicons extracted: {len(amplicons)}", file=sys.stderr)
    print(f"Failed/no sequence: {len(failed)}", file=sys.stderr)
    print(f"Output: {outfile}", file=sys.stderr)

    if failed:
        print(f"\nFailed methods: {', '.join(failed)}", file=sys.stderr)


if __name__ == "__main__":
    main()
