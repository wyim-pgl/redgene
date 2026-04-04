"""EUginius method scraper - refactored

Input:
  - EUginius_method_code.txt

Outputs:
  - euginius_primers.tsv
  - euginius_primers.fasta
  - euginius_failures.tsv
  - optional HTML cache directory

Usage:
  python euginius_scraper_refactored.py \
      --input EUginius_method_code.txt \
      --out-tsv euginius_primers.tsv \
      --out-fasta euginius_primers.fasta \
      --out-fail euginius_failures.tsv \
      --cache-dir html_cache
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, List, Optional, Tuple
from urllib.parse import quote

import requests
from bs4 import BeautifulSoup
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


BASE_URL = "https://www.euginius.eu/euginius/pages/method_detailview.jsf?method={}"


USER_AGENT = (
    "Mozilla/5.0 (X11; Linux x86_64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/123.0 Safari/537.36"
)


@dataclass
class MethodRecord:
    source_db: str = "EUginius"
    method_code: str = ""
    method_name: str = ""
    method_type: str = ""
    target_scope: str = ""
    target_name: str = ""
    description: str = ""
    forward_primer_name: str = ""
    forward_primer_seq: str = ""
    reverse_primer_name: str = ""
    reverse_primer_seq: str = ""
    probe_name: str = ""
    probe_seq: str = ""
    amplicon_size: str = ""
    amplicon_seq: str = ""
    source_url: str = ""
    notes: str = ""


@dataclass
class FailureRecord:
    method_code: str
    source_url: str
    error: str


def build_session(total_retries: int = 3, backoff: float = 1.0) -> requests.Session:
    retry = Retry(
        total=total_retries,
        connect=total_retries,
        read=total_retries,
        status=total_retries,
        backoff_factor=backoff,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


def read_method_codes(path: Path) -> List[str]:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    seen = set()
    result = []

    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            code = raw.strip()
            if not code:
                continue
            if code.startswith("#"):
                continue
            if code not in seen:
                seen.add(code)
                result.append(code)

    return result


def clean_text(text: Optional[str]) -> str:
    if not text:
        return ""
    text = text.replace("\xa0", " ")
    text = re.sub(r"\r", "", text)
    text = re.sub(r"[ \t\f\v]+", " ", text)
    text = re.sub(r"\n\s*\n+", "\n", text)
    return text.strip()


def normalize_seq(seq: str) -> str:
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def safe_filename(text: str) -> str:
    text = text.strip()
    text = re.sub(r"[^\w.\-]+", "_", text)
    return text[:180] if text else "NA"


def safe_header_token(text: str) -> str:
    if not text:
        return "NA"
    text = text.strip()
    text = text.replace("|", "_")
    text = text.replace("/", "_")
    text = text.replace("\\", "_")
    text = text.replace("'", "")
    text = text.replace('"', "")
    text = re.sub(r"\s+", "_", text)
    return text or "NA"


def method_url(method_code: str) -> str:
    return BASE_URL.format(quote(method_code, safe=""))


def fetch_html(
    session: requests.Session,
    method_code: str,
    timeout: int = 40,
    sleep_sec: float = 0.5,
) -> Tuple[str, str]:
    url = method_url(method_code)
    resp = session.get(url, timeout=timeout)
    resp.raise_for_status()
    html = resp.text
    if sleep_sec > 0:
        time.sleep(sleep_sec)
    return html, url


def maybe_save_html(cache_dir: Optional[Path], method_code: str, html: str) -> None:
    if cache_dir is None:
        return
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / f"{safe_filename(method_code)}.html"
    out.write_text(html, encoding="utf-8")


def soup_text(html: str) -> str:
    soup = BeautifulSoup(html, "html.parser")
    return clean_text(soup.get_text("\n", strip=True))


def extract_block(text: str, start_label: str, end_labels: Iterable[str]) -> str:
    idx = text.find(start_label)
    if idx == -1:
        return ""
    sub = text[idx + len(start_label):]
    hits = []
    for label in end_labels:
        pos = sub.find(label)
        if pos != -1:
            hits.append(pos)
    end = min(hits) if hits else len(sub)
    return clean_text(sub[:end])


def first_nonempty(*values: str) -> str:
    for v in values:
        if v:
            return v
    return ""


def parse_name_seq_block(block: str) -> Tuple[str, str]:
    if not block:
        return "", ""

    name = ""
    seq = ""

    m_name = re.search(r"Name:\s*(.+?)(?:\n|Sequence:)", block, flags=re.I | re.S)
    if m_name:
        name = clean_text(m_name.group(1))

    m_seq = re.search(r"Sequence:\s*([ACGTUacgtuNnRYSWKMBDHVryswkmbdhv \n\r\t\-\.\(\)]+?)(?:\s*(?:Size:|Name:|$))", block, flags=re.I)
    if m_seq:
        seq = normalize_seq(m_seq.group(1))

    return name, seq


def parse_amplicon_block(block: str) -> Tuple[str, str]:
    if not block:
        return "", ""

    size = ""
    seq = ""

    m_size = re.search(r"Size:\s*([0-9]+)", block, flags=re.I)
    if m_size:
        size = m_size.group(1)

    m_seq = re.search(r"Sequence:\s*([ACGTUacgtuNnRYSWKMBDHVryswkmbdhv \n\r\t\-\.\(\)]+?)(?:\s*(?:Size:|Name:|$))", block, flags=re.I)
    if m_seq:
        seq = normalize_seq(m_seq.group(1))

    return size, seq


def infer_target_scope(method_code: str, method_type: str) -> str:
    mt = (method_type or "").lower()

    if "event" in mt:
        return "event-specific"
    if "construct" in mt:
        return "construct-specific"
    if "element" in mt:
        return "element-specific"
    if "taxon" in mt or "species" in mt:
        return "taxon-specific"
    if "plant" in mt:
        return "plant-specific"

    if "-EVE-" in method_code:
        return "event-specific"
    if "-CON-" in method_code:
        return "construct-specific"
    if "-ELE-" in method_code:
        return "element-specific"
    if "-TAX-" in method_code:
        return "taxon-specific"
    if "-PLN-" in method_code:
        return "plant-specific"
    if "-BAC-" in method_code:
        return "bacterial-control"

    return ""


def parse_method_page(html: str, method_code: str, url: str) -> MethodRecord:
    text = soup_text(html)

    record = MethodRecord(
        method_code=method_code,
        source_url=url,
    )

    record.method_name = extract_block(
        text,
        "Name:",
        ["Description:", "Comment:", "Validation:", "Standardisation:", "Type:"]
    )

    record.description = extract_block(
        text,
        "Description:",
        [
            "Comment:",
            "Validation:",
            "Standardisation:",
            "Type:",
            "Target GMO name:",
            "Target GMO names:",
            "Target Species name:",
            "Target DNA element:",
            "Oligonucleotides:",
            "Documents",
        ]
    )

    record.method_type = extract_block(
        text,
        "Type:",
        [
            "Target GMO name:",
            "Target GMO names:",
            "Target Species name:",
            "Target DNA element:",
            "Oligonucleotides:",
            "Documents",
        ]
    )

    record.target_scope = infer_target_scope(method_code, record.method_type)

    record.target_name = first_nonempty(
        extract_block(text, "Target GMO name:", ["Target GMO names:", "Target Species name:", "Target DNA element:", "Oligonucleotides:", "Documents"]),
        extract_block(text, "Target GMO names:", ["Target Species name:", "Target DNA element:", "Oligonucleotides:", "Documents"]),
        extract_block(text, "Target Species name:", ["Target DNA element:", "Oligonucleotides:", "Documents"]),
        extract_block(text, "Target DNA element:", ["Oligonucleotides:", "Documents"]),
    )

    fwd_block = extract_block(text, "Forward Primer", ["Reverse Primer", "Probe", "Amplicon:", "Documents"])
    rev_block = extract_block(text, "Reverse Primer", ["Probe", "Amplicon:", "Documents"])
    probe_block = extract_block(text, "Probe", ["Amplicon:", "Documents"])
    amp_block = extract_block(text, "Amplicon:", ["Documents"])

    record.forward_primer_name, record.forward_primer_seq = parse_name_seq_block(fwd_block)
    record.reverse_primer_name, record.reverse_primer_seq = parse_name_seq_block(rev_block)
    record.probe_name, record.probe_seq = parse_name_seq_block(probe_block)
    record.amplicon_size, record.amplicon_seq = parse_amplicon_block(amp_block)

    parsed_any = any([
        record.forward_primer_seq,
        record.reverse_primer_seq,
        record.probe_seq,
        record.amplicon_seq,
    ])

    if not parsed_any:
        record.notes = "No oligo or amplicon sequence parsed. Inspect page manually."

    return record


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_tsv(records: List[MethodRecord], out_path: Path) -> None:
    fieldnames = list(MethodRecord().__dict__.keys())
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in records:
            writer.writerow(asdict(rec))


def write_failures(failures: List[FailureRecord], out_path: Path) -> None:
    fieldnames = ["method_code", "source_url", "error"]
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in failures:
            writer.writerow(asdict(rec))


def write_fasta(records: List[MethodRecord], out_path: Path) -> None:
    with out_path.open("w", encoding="utf-8") as fh:
        for rec in records:
            entries = [
                ("forward_primer", rec.forward_primer_name, rec.forward_primer_seq),
                ("reverse_primer", rec.reverse_primer_name, rec.reverse_primer_seq),
                ("probe", rec.probe_name, rec.probe_seq),
                ("amplicon", "NA", rec.amplicon_seq),
            ]
            for seq_type, name, seq in entries:
                if not seq:
                    continue
                header = (
                    f">euginius|{safe_header_token(rec.method_code)}|"
                    f"{seq_type}|{safe_header_token(name)}"
                )
                fh.write(header + "\n")
                fh.write(wrap_fasta(seq) + "\n")


def run_scrape(
    method_codes: List[str],
    session: requests.Session,
    cache_dir: Optional[Path],
    timeout: int,
    sleep_sec: float,
) -> Tuple[List[MethodRecord], List[FailureRecord]]:
    records: List[MethodRecord] = []
    failures: List[FailureRecord] = []

    total = len(method_codes)

    for i, code in enumerate(method_codes, start=1):
        url = method_url(code)
        print(f"[{i}/{total}] {code}", file=sys.stderr)
        try:
            html, _ = fetch_html(session, code, timeout=timeout, sleep_sec=sleep_sec)
            maybe_save_html(cache_dir, code, html)
            record = parse_method_page(html, code, url)
            records.append(record)
        except Exception as exc:
            failures.append(FailureRecord(
                method_code=code,
                source_url=url,
                error=str(exc),
            ))

    return records, failures


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Scrape EUginius method detail pages from method code list.")
    p.add_argument("--input", default="EUginius_method_code.txt", help="Input text file with one method code per line")
    p.add_argument("--out-tsv", default="euginius_primers.tsv", help="Output TSV path")
    p.add_argument("--out-fasta", default="euginius_primers.fasta", help="Output FASTA path")
    p.add_argument("--out-fail", default="euginius_failures.tsv", help="Output failure TSV path")
    p.add_argument("--cache-dir", default="", help="Optional directory to save raw HTML pages")
    p.add_argument("--timeout", type=int, default=40, help="HTTP timeout in seconds")
    p.add_argument("--sleep", type=float, default=0.6, help="Sleep between requests in seconds")
    p.add_argument("--retries", type=int, default=3, help="HTTP retry count")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    input_path = Path(args.input)
    out_tsv = Path(args.out_tsv)
    out_fasta = Path(args.out_fasta)
    out_fail = Path(args.out_fail)
    cache_dir = Path(args.cache_dir) if args.cache_dir else None

    method_codes = read_method_codes(input_path)
    if not method_codes:
        raise ValueError(f"No method codes found in {input_path}")

    print(f"[INFO] Loaded {len(method_codes)} method codes", file=sys.stderr)

    session = build_session(total_retries=args.retries, backoff=1.0)

    records, failures = run_scrape(
        method_codes=method_codes,
        session=session,
        cache_dir=cache_dir,
        timeout=args.timeout,
        sleep_sec=args.sleep,
    )

    write_tsv(records, out_tsv)
    write_fasta(records, out_fasta)
    write_failures(failures, out_fail)

    print(f"[INFO] Parsed methods : {len(records)}", file=sys.stderr)
    print(f"[INFO] Failures      : {len(failures)}", file=sys.stderr)
    print(f"[INFO] TSV           : {out_tsv}", file=sys.stderr)
    print(f"[INFO] FASTA         : {out_fasta}", file=sys.stderr)
    print(f"[INFO] Failure TSV   : {out_fail}", file=sys.stderr)
    if cache_dir:
        print(f"[INFO] HTML cache    : {cache_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
