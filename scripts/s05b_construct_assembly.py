#!/usr/bin/env python3
"""Step 5b: Reconstruct and characterize the inserted transgene construct.

Reuses s04b's pre-filter SPAdes contigs (``contigs_all.fasta``), classifies each
contig by host/element BLAST coverage (assemble-then-classify), annotates
non-host contigs (local element BLAST always; opt-in remote nt for unknowns),
and links contigs back to s05 CANDIDATE insertion sites via minimap2.

Output: results/<sample>/s05b_construct_assembly/
  construct_contigs.fasta       non-host contigs (chimeric written whole)
  contig_classification.tsv     all contigs: class + coverage fractions + annotation
  site_construct_links.tsv      s05 site_id -> contig links
  construct_summary.txt         sample-level inventory

Exit code 0 in all input-shortfall paths (empty/missing inputs still produce the
expected output files), so downstream steps can always find them.

Self-contained (like s04b): a local ``_read_fasta`` is used instead of importing
``scripts/s05/primitives`` so the step does not pull the heavy s05 package
(pysam et al.) at import time.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple

# --- classification thresholds (see spec 2026-06-30-step5b-construct-assembly) ---
MIN_CONTIG_LEN = 200
HOST_PIDENT = 90.0
ELEMENT_PIDENT = 80.0
CHIMERIC_MIN = 0.20
ELEMENT_MIN = 0.50
HOST_MIN = 0.70
LINK_SLOP = 100


def log(msg: str) -> None:
    print(f"[s05b] {msg}", file=sys.stderr)


def _read_fasta(path: Path) -> dict[str, str]:
    """name -> sequence (first whitespace-delimited token of the header is the id)."""
    seqs: dict[str, str] = {}
    name = ""
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts)
                name = line[1:].split()[0] if len(line) > 1 else ""
                parts = []
            else:
                parts.append(line)
    if name:
        seqs[name] = "".join(parts)
    return seqs


def _is_nonempty_file(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------
def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping/adjacent 1-based inclusive intervals."""
    if not intervals:
        return []
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def _coverage_fraction(hits: list[tuple[int, int]], length: int) -> float:
    """Merged covered bp / length. 0.0 when length <= 0."""
    if length <= 0:
        return 0.0
    covered = sum(end - start + 1 for start, end in _merge_intervals(hits))
    return covered / length


def _parse_blast6(text: str, min_pident: float) -> dict[str, list[tuple[int, int]]]:
    """qseqid -> [(qstart, qend)] for rows passing pident >= min_pident.

    outfmt: 6 qseqid sseqid pident length mismatch gapopen qstart qend
            sstart send evalue bitscore
    """
    out: dict[str, list[tuple[int, int]]] = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        qseqid = parts[0]
        try:
            pident = float(parts[2])
            qstart = int(parts[6])
            qend = int(parts[7])
        except ValueError:
            continue
        if pident < min_pident:
            continue
        out.setdefault(qseqid, []).append((min(qstart, qend), max(qstart, qend)))
    return out


def _classify_contig(host_frac: float, element_frac: float) -> str:
    """Assign a contig class from host/element coverage fractions.

    Priority (first match wins):
      1. both >= CHIMERIC_MIN          -> chimeric (junction-spanning)
      2. element_frac >= ELEMENT_MIN   -> construct
      3. host_frac >= HOST_MIN         -> host (discarded from FASTA)
      4. otherwise                     -> foreign_unknown
    """
    if host_frac >= CHIMERIC_MIN and element_frac >= CHIMERIC_MIN:
        return "chimeric"
    if element_frac >= ELEMENT_MIN:
        return "construct"
    if host_frac >= HOST_MIN:
        return "host"
    return "foreign_unknown"


# ---------------------------------------------------------------------------
# BLAST / alignment runners
# ---------------------------------------------------------------------------
def _run_blastn(
    query: Path, subjects: list[Path], min_pident: float
) -> dict[str, list[tuple[int, int]]]:
    """Pool blastn hits of `query` vs each subject FASTA. qseqid -> [(qstart,qend)]."""
    pooled: dict[str, list[tuple[int, int]]] = {}
    if not _is_nonempty_file(query):
        return pooled
    for subj in subjects:
        if not _is_nonempty_file(Path(subj)):
            log(f"blastn: skipping missing/empty subject {subj}")
            continue
        cmd = [
            "blastn", "-query", str(query), "-subject", str(subj),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                       "qstart qend sstart send evalue bitscore",
            "-evalue", "1e-10",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            log(f"blastn vs {subj} failed (rc={proc.returncode}): {proc.stderr.strip()}")
            continue
        for qseqid, hits in _parse_blast6(proc.stdout, min_pident).items():
            pooled.setdefault(qseqid, []).extend(hits)
    return pooled


def _contig_lengths(contigs_fa: Path) -> dict[str, int]:
    if not _is_nonempty_file(contigs_fa):
        return {}
    return {cid: len(seq) for cid, seq in _read_fasta(contigs_fa).items()}


def _top_element_hits(
    contigs_fa: Path, elem_dbs: list[Path]
) -> dict[str, tuple[str, float, int]]:
    """contig_id -> (sseqid, pident, alnlen) for the best (highest pident) element hit."""
    best: dict[str, tuple[str, float, int]] = {}
    if not _is_nonempty_file(contigs_fa):
        return best
    for db in elem_dbs:
        if not _is_nonempty_file(Path(db)):
            continue
        cmd = ["blastn", "-query", str(contigs_fa), "-subject", str(db),
               "-outfmt", "6 qseqid sseqid pident length", "-evalue", "1e-10"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            continue
        for line in proc.stdout.splitlines():
            p = line.split("\t")
            if len(p) < 4:
                continue
            qid, sid = p[0], p[1]
            try:
                pid, aln = float(p[2]), int(p[3])
            except ValueError:
                continue
            cur = best.get(qid)
            if cur is None or pid > cur[1]:
                best[qid] = (sid, pid, aln)
    return best


# ---------------------------------------------------------------------------
# s05 site parsing
# ---------------------------------------------------------------------------
class Site(NamedTuple):
    site_id: str
    chrom: str
    pos: int
    verdict: str


def _parse_s05_sites(s05_dir: Path, verdicts: set[str] | None = None) -> list[Site]:
    """Parse insertion_*_report.txt into Site rows, filtered by verdict.

    Filename form: insertion_<chrom>_<pos>_report.txt ; <chrom> may contain
    underscores, so site_id = text between 'insertion_' and '_report', and
    chrom/pos split on the LAST underscore of site_id.
    """
    wanted = {"CANDIDATE"} if verdicts is None else verdicts
    sites: list[Site] = []
    if not s05_dir.is_dir():
        return sites
    for report in sorted(s05_dir.glob("insertion_*_report.txt")):
        name = report.name
        if not name.startswith("insertion_") or not name.endswith("_report.txt"):
            continue
        site_id = name[len("insertion_"):-len("_report.txt")]
        if "_" not in site_id:
            continue
        chrom, pos_s = site_id.rsplit("_", 1)
        try:
            pos = int(pos_s)
        except ValueError:
            continue
        verdict = ""
        for line in report.read_text().splitlines():
            if line.startswith("Verdict:"):
                verdict = line.split(":", 1)[1].strip()
                break
        if verdict in wanted:
            sites.append(Site(site_id, chrom, pos, verdict))
    return sites


# ---------------------------------------------------------------------------
# PAF parse + per-site linkage
# ---------------------------------------------------------------------------
class Aln(NamedTuple):
    contig_id: str
    chrom: str
    tstart: int
    tend: int


def _parse_paf(text: str) -> list[Aln]:
    alns: list[Aln] = []
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        try:
            alns.append(Aln(parts[0], parts[5], int(parts[7]), int(parts[8])))
        except ValueError:
            continue
    return alns


class Link(NamedTuple):
    site_id: str
    verdict: str
    contig_id: str
    junction_side: str
    link_confidence: str


def _link_contigs_to_sites(
    alns: list[Aln],
    sites: list[Site],
    contig_class: dict[str, str],
    slop: int = LINK_SLOP,
) -> list[Link]:
    """Link contigs to s05 sites where a host alignment abuts the site pos."""
    links: list[Link] = []
    seen: set[tuple[str, str]] = set()
    for site in sites:
        for aln in alns:
            if aln.chrom != site.chrom:
                continue
            near_start = abs(aln.tstart - site.pos) <= slop
            near_end = abs(aln.tend - site.pos) <= slop
            if not (near_start or near_end):
                continue
            key = (site.site_id, aln.contig_id)
            if key in seen:
                continue
            seen.add(key)
            side = "downstream" if near_start else "upstream"
            conf = "high" if contig_class.get(aln.contig_id) == "chimeric" else "medium"
            links.append(Link(site.site_id, site.verdict, aln.contig_id, side, conf))
    return links


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------
class ContigRow(NamedTuple):
    contig_id: str
    length: int
    cls: str
    host_frac: float
    element_frac: float
    top_element: str
    element_pident: float
    element_alnlen: int
    remote_hit: str


_CLS_TSV_HEADER = (
    "contig_id\tlength\tclass\thost_frac\telement_frac\t"
    "top_element\telement_pident\telement_alnlen\tremote_hit"
)
_LINKS_TSV_HEADER = "site_id\tverdict\tcontig_id\tjunction_side\tlink_confidence"


def write_outputs(
    step_dir: Path,
    contigs_fa: Path,
    rows: list[ContigRow],
    links: list[Link],
) -> None:
    step_dir.mkdir(parents=True, exist_ok=True)
    seqs = _read_fasta(contigs_fa) if _is_nonempty_file(contigs_fa) else {}

    # FASTA: non-host contigs, whole, annotated header
    with (step_dir / "construct_contigs.fasta").open("w") as fh:
        for r in rows:
            if r.cls == "host":
                continue
            seq = seqs.get(r.contig_id)
            if not seq:
                continue
            fh.write(
                f">{r.contig_id} class={r.cls} top_element={r.top_element or 'NA'} "
                f"host_frac={r.host_frac:.3f} element_frac={r.element_frac:.3f}\n"
            )
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    # classification TSV: ALL rows (host included, for transparency)
    with (step_dir / "contig_classification.tsv").open("w") as fh:
        fh.write(_CLS_TSV_HEADER + "\n")
        for r in rows:
            fh.write(
                f"{r.contig_id}\t{r.length}\t{r.cls}\t{r.host_frac:.3f}\t"
                f"{r.element_frac:.3f}\t{r.top_element}\t{r.element_pident:.1f}\t"
                f"{r.element_alnlen}\t{r.remote_hit}\n"
            )

    # links TSV
    with (step_dir / "site_construct_links.tsv").open("w") as fh:
        fh.write(_LINKS_TSV_HEADER + "\n")
        for lk in links:
            fh.write(
                f"{lk.site_id}\t{lk.verdict}\t{lk.contig_id}\t"
                f"{lk.junction_side}\t{lk.link_confidence}\n"
            )

    # summary
    non_host = [r for r in rows if r.cls != "host"]
    by_cls: dict[str, int] = {}
    for r in rows:
        by_cls[r.cls] = by_cls.get(r.cls, 0) + 1
    elements = sorted({r.top_element for r in non_host if r.top_element})
    construct_bp = sum(r.length for r in non_host)
    has_unknown = any(r.cls == "foreign_unknown" for r in rows)
    with (step_dir / "construct_summary.txt").open("w") as fh:
        fh.write("# Step 5b construct inventory\n")
        fh.write(f"total_construct_bp\t{construct_bp}\n")
        for cls in ("construct", "chimeric", "foreign_unknown", "host"):
            fh.write(f"contigs_{cls}\t{by_cls.get(cls, 0)}\n")
        fh.write(f"distinct_elements\t{len(elements)}\n")
        fh.write(f"element_inventory\t{','.join(elements) if elements else 'NA'}\n")
        fh.write(f"sites_linked\t{len({lk.site_id for lk in links})}\n")
        fh.write(f"novel_payload_detected\t{'yes' if has_unknown else 'no'}\n")


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def _run_minimap2(contigs_fa: Path, host_ref: Path, threads: int) -> list[Aln]:
    if not _is_nonempty_file(contigs_fa) or not _is_nonempty_file(host_ref):
        return []
    if shutil.which("minimap2") is None:
        log("minimap2 not on PATH; skipping per-site linkage")
        return []
    cmd = ["minimap2", "-x", "asm10", "-t", str(threads), str(host_ref), str(contigs_fa)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        log(f"minimap2 failed (rc={proc.returncode}): {proc.stderr.strip()[:200]}")
        return []
    return _parse_paf(proc.stdout)


def _obtain_contigs(args: argparse.Namespace, step_dir: Path) -> Path:
    """Reuse s04b contigs_all.fasta, else SPAdes-fallback on s03 reads, else empty."""
    if _is_nonempty_file(args.contigs_all):
        log(f"using s04b contigs: {args.contigs_all}")
        return args.contigs_all
    local = step_dir / "contigs.fasta"
    r1, r2 = args.s03_r1, args.s03_r2
    if not (_is_nonempty_file(r1) and _is_nonempty_file(r2)):
        log("no contigs_all and no s03 reads; empty contigs")
        local.write_text("")
        return local
    if shutil.which("spades.py") is None:
        log("spades.py not on PATH; empty contigs")
        local.write_text("")
        return local
    spades_out = step_dir / "_spades_run"
    if spades_out.exists():
        shutil.rmtree(spades_out)
    cmd = ["spades.py", "--only-assembler", "--careful",
           "-1", str(r1), "-2", str(r2), "-o", str(spades_out),
           "-t", str(args.threads), "-m", str(args.memory_gb)]
    log(" ".join(cmd))
    proc = subprocess.run(cmd, stderr=subprocess.DEVNULL)
    sp_contigs = spades_out / "contigs.fasta"
    if proc.returncode == 0 and _is_nonempty_file(sp_contigs):
        shutil.copy2(sp_contigs, local)
    else:
        log(f"SPAdes fallback produced no contigs (rc={proc.returncode}); empty")
        local.write_text("")
    shutil.rmtree(spades_out, ignore_errors=True)
    return local


def _length_filtered(contigs_fa: Path, step_dir: Path, min_len: int) -> Path:
    """Write a length-filtered copy; returns its path (may be empty)."""
    out = step_dir / "contigs_filtered.fasta"
    seqs = _read_fasta(contigs_fa) if _is_nonempty_file(contigs_fa) else {}
    with out.open("w") as fh:
        for cid, seq in seqs.items():
            if len(seq) >= min_len:
                fh.write(f">{cid}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
    return out


def _annotate_remote_unknowns(
    rows: list[ContigRow], contigs_fa: Path, dbg: Path
) -> None:
    """Run remote nt on foreign_unknown contigs; mutate rows in place. Never raises."""
    unknown_ids = [r.contig_id for r in rows if r.cls == "foreign_unknown"]
    if not unknown_ids:
        return
    seqs = _read_fasta(contigs_fa)
    q = dbg / "unknown.fasta"
    with q.open("w") as fh:
        for cid in unknown_ids:
            if seqs.get(cid):
                fh.write(f">{cid}\n{seqs[cid]}\n")
    if not _is_nonempty_file(q):
        return
    cmd = ["blastn", "-query", str(q), "-db", "nt", "-remote",
           "-evalue", "1e-10", "-max_target_seqs", "5",
           "-outfmt", "6 qseqid stitle pident length"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except (subprocess.TimeoutExpired, OSError) as exc:
        log(f"remote nt blast failed/timeout: {exc}; local-only annotation")
        return
    if proc.returncode != 0:
        log(f"remote nt blast rc={proc.returncode}; local-only annotation")
        return
    hit: dict[str, str] = {}
    for line in proc.stdout.splitlines():
        p = line.split("\t")
        if len(p) >= 2 and p[0] not in hit:
            hit[p[0]] = p[1][:120]
    for i, r in enumerate(rows):
        if r.contig_id in hit:
            rows[i] = r._replace(remote_hit=hit[r.contig_id])


def run(args: argparse.Namespace) -> int:
    step_dir = args.outdir / args.sample_name / "s05b_construct_assembly"
    step_dir.mkdir(parents=True, exist_ok=True)
    dbg = step_dir / "_classification_blast"
    dbg.mkdir(exist_ok=True)

    raw_contigs = _obtain_contigs(args, step_dir)
    contigs = _length_filtered(raw_contigs, step_dir, MIN_CONTIG_LEN)
    lengths = _contig_lengths(contigs)

    if not lengths:
        log("0 contigs after length filter; writing empty outputs")
        write_outputs(step_dir, contigs, [], [])
        return 0

    host_hits = _run_blastn(contigs, [args.host_ref], HOST_PIDENT) if args.host_ref else {}
    elem_dbs = list(args.element_db or [])
    elem_hits = _run_blastn(contigs, elem_dbs, ELEMENT_PIDENT)
    top_elem = _top_element_hits(contigs, elem_dbs)

    rows: list[ContigRow] = []
    contig_class: dict[str, str] = {}
    for cid, length in lengths.items():
        hf = _coverage_fraction(host_hits.get(cid, []), length)
        ef = _coverage_fraction(elem_hits.get(cid, []), length)
        cls = _classify_contig(hf, ef)
        contig_class[cid] = cls
        te = top_elem.get(cid, ("", 0.0, 0))
        rows.append(ContigRow(cid, length, cls, hf, ef, te[0], te[1], te[2], ""))

    if args.remote_blast:
        _annotate_remote_unknowns(rows, contigs, dbg)

    # per-site linkage: build a non-host FASTA, map to host, link
    tmp_nh = step_dir / "_nonhost_for_link.fasta"
    seqs = _read_fasta(contigs)
    with tmp_nh.open("w") as fh:
        for r in rows:
            if r.cls != "host" and seqs.get(r.contig_id):
                fh.write(f">{r.contig_id}\n{seqs[r.contig_id]}\n")
    sites = _parse_s05_sites(args.s05_dir) if args.s05_dir else []
    alns = (
        _run_minimap2(tmp_nh, args.host_ref, args.threads)
        if (sites and args.host_ref) else []
    )
    links = _link_contigs_to_sites(alns, sites, contig_class)
    tmp_nh.unlink(missing_ok=True)

    write_outputs(step_dir, contigs, rows, links)
    log(f"done: {len(rows)} contigs, {len(links)} site links")
    return 0


def parse_args_from(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--contigs-all", type=Path, required=True,
                   help="s04b contigs_all.fasta (pre-filter SPAdes). Reused if present.")
    p.add_argument("--s03-r1", type=Path, required=True, help="SPAdes-fallback R1")
    p.add_argument("--s03-r2", type=Path, required=True, help="SPAdes-fallback R2")
    p.add_argument("--s05-dir", type=Path, default=None,
                   help="s05_insert_assembly dir for per-site linkage")
    p.add_argument("--host-ref", type=Path, default=None)
    p.add_argument("--element-db", type=Path, action="append", default=None,
                   help="element FASTA (repeatable)")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--memory-gb", type=int, default=16)
    p.add_argument("--remote-blast", action="store_true",
                   help="Run NCBI nt remote BLAST on foreign_unknown contigs (opt-in).")
    return p.parse_args(argv)


def parse_args() -> argparse.Namespace:
    return parse_args_from(None)


if __name__ == "__main__":
    sys.exit(run(parse_args()))
