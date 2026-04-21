#!/usr/bin/env python3
"""Extract non-host gap segments from an assembled insert FASTA given
the host BLAST hit TSV. Writes gaps >= min_gap bp to a single FASTA.

Host TSV columns (qname, qstart, qend, sname, identity, aln_len).
"""
import argparse
import sys
from pathlib import Path


def read_fasta(path):
    name, seq = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line)
    if name:
        yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--insert-fasta", type=Path, required=True)
    ap.add_argument("--host-tsv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-gap", type=int, default=200)
    args = ap.parse_args()

    seqs = dict(read_fasta(args.insert_fasta))
    if len(seqs) != 1:
        print(f"Expected one sequence, got {len(seqs)}", file=sys.stderr)
        return 2
    (insert_name, insert_seq), = seqs.items()
    insert_len = len(insert_seq)

    # Collect host hit intervals (1-based -> 0-based half-open)
    intervals = []
    for line in args.host_tsv.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 3:
            continue
        try:
            qs, qe = int(f[1]), int(f[2])
        except ValueError:
            continue
        s, e = (qs - 1, qe) if qs <= qe else (qe - 1, qs)
        intervals.append((s, e))

    # Merge overlapping
    intervals.sort()
    merged = []
    for s, e in intervals:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))

    # Complement -> gaps
    gaps = []
    prev = 0
    for s, e in merged:
        if s > prev:
            gaps.append((prev, s))
        prev = e
    if prev < insert_len:
        gaps.append((prev, insert_len))

    with args.out.open("w") as fh:
        n = 0
        for s, e in gaps:
            if e - s < args.min_gap:
                continue
            n += 1
            header = f"{insert_name}_gap{n}_{s+1}-{e}"
            fh.write(f">{header}\n")
            sub = insert_seq[s:e]
            for i in range(0, len(sub), 80):
                fh.write(sub[i:i+80] + "\n")
            print(f"  gap{n}: {s+1}-{e} ({e-s}bp) -> {header}",
                  file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main() or 0)
