"""Phase 4 — Annotation (Issue #4, Session 3).

Extracted verbatim from ``scripts/s05_insert_assembly.py`` so the BLAST-driven
annotation pipeline (local element_db + extras + optional NCBI nt) can be
unit-tested and re-imported without pulling in the full monolith.

Pipeline overview
-----------------
``annotate_insert`` is the public entry point. It runs:

1. Local BLAST against the primary ``element_db`` and any ``extra_dbs``
   (e.g. ``element_db/common_payload.fa`` and per-sample s04b SPAdes
   contigs) via ``_run_local_blast``, writing parsed hits as a list of
   dicts produced by ``_parse_blast6``.
2. Optional remote BLAST against NCBI ``nt`` via ``_run_remote_blast``
   (skipped when ``no_remote_blast=True``; default True for unit tests).
3. Merge of the two streams via ``_merge_annotations`` using best-bitscore
   per query region with element_db preferred at ties (Task T7 priority
   rule, see ``tests/test_extra_element_db.py`` regression guards).
4. Emission of ``element_annotation.tsv`` and a T-DNA border-motif scan
   into ``border_hits.tsv`` (consumed by the Session 6 report module).

The ``|src=<tag>`` suffix that the v2 element DB bakes into every header
is stripped here via ``_parse_src_tag`` (re-exported by
``scripts/s05/classify.py``); the extracted tag is carried in the
``src_tag`` field of each hit dict for downstream consumers.

This module does no verdict-style classification — it only produces hits.
The monolith's ``classify_site_tiers`` consumes these hits later.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from .primitives import log, _parse_src_tag


# ---------------------------------------------------------------------------
# BLAST output parsing
# ---------------------------------------------------------------------------

def _parse_blast6(path: Path, min_len: int = 30) -> list[dict]:
    """Parse BLAST outfmt-6 with 10+ columns into list of dicts.

    T5: also strip the |src=<tag> suffix that the v2 DB bakes into every
    header so the ``subject`` field stays clean for element_annotation.tsv.
    The extracted tag is carried in the ``src_tag`` key for downstream
    consumers; callers that don't need it just ignore it.
    """
    hits: list[dict] = []
    if not path.exists():
        return hits
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip().split("\t")
            if len(cols) < 10:
                continue
            aln_len = int(cols[3])
            if aln_len < min_len:
                continue
            subject = cols[1]
            subject_clean, src_tag = _parse_src_tag(subject, default="")
            hits.append({
                "query": cols[0], "subject": subject_clean,
                "identity": float(cols[2]), "length": aln_len,
                "q_start": int(cols[4]), "q_end": int(cols[5]),
                "s_start": int(cols[6]), "s_end": int(cols[7]),
                "evalue": float(cols[8]), "bitscore": float(cols[9]),
                "src_tag": src_tag,
            })
    return hits


# ---------------------------------------------------------------------------
# Local + remote BLAST execution
# ---------------------------------------------------------------------------

def _run_local_blast(
    insert_fasta: Path, element_db: Path, output_dir: Path,
    tag: str | None = None,
) -> list[dict]:
    """Local BLAST vs element_db. Fast, covers known GMO elements.

    ``tag`` disambiguates intermediate filenames when the function is
    called multiple times in the same ``output_dir`` (e.g. once for the
    primary element_db and again for each extra DB). Defaults to the
    stem of ``element_db`` so callers can omit it.
    """
    suffix = tag if tag is not None else element_db.stem
    db_prefix = output_dir / f"_element_blastdb_{suffix}"

    def _cleanup_db() -> None:
        for ext in [".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto", ".njs"]:
            Path(f"{db_prefix}{ext}").unlink(missing_ok=True)

    # Local BLAST is best-effort annotation: a failure (disk full, makeblastdb
    # OOM, missing binary) must degrade to "no hits" rather than crash the whole
    # step — a site can still get a verdict from element tiers. Mirrors the
    # graceful handling in _run_remote_blast.
    mk = subprocess.run(
        ["makeblastdb", "-in", str(element_db), "-dbtype", "nucl",
         "-out", str(db_prefix)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    if mk.returncode != 0:
        log(f"  Local BLAST ({suffix}): makeblastdb failed "
            f"(rc={mk.returncode}); skipping (0 hits)")
        _cleanup_db()
        return []
    blast_out = output_dir / f"_local_blast_{suffix}.tsv"
    bl = subprocess.run(
        ["blastn", "-query", str(insert_fasta), "-db", str(db_prefix),
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
         "-evalue", "1e-5", "-max_target_seqs", "50",
         "-out", str(blast_out)],
        stderr=subprocess.DEVNULL,
    )
    if bl.returncode != 0:
        log(f"  Local BLAST ({suffix}): blastn failed "
            f"(rc={bl.returncode}); skipping (0 hits)")
        blast_out.unlink(missing_ok=True)
        _cleanup_db()
        return []
    hits = _parse_blast6(blast_out)
    # Cleanup
    blast_out.unlink(missing_ok=True)
    _cleanup_db()
    log(f"  Local BLAST ({suffix}): {len(hits)} hits")
    return hits


def _run_remote_blast(
    insert_fasta: Path, output_dir: Path,
    timeout: int = 600, max_retries: int = 2,
) -> list[dict]:
    """Remote BLAST vs NCBI nt. Slower, but annotates unknown regions."""
    blast_out = output_dir / "_remote_blast.tsv"
    for attempt in range(1, max_retries + 1):
        log(f"  Remote BLAST vs NCBI nt (attempt {attempt}/{max_retries})...")
        try:
            proc = subprocess.run(
                ["blastn", "-query", str(insert_fasta), "-db", "nt",
                 "-remote",
                 "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle",
                 "-evalue", "1e-10", "-max_target_seqs", "10",
                 "-out", str(blast_out)],
                stderr=subprocess.PIPE, timeout=timeout,
            )
            if proc.returncode == 0 and blast_out.exists():
                break
            log(f"    Remote BLAST returned code {proc.returncode}")
        except subprocess.TimeoutExpired:
            log(f"    Remote BLAST timed out ({timeout}s)")
        except FileNotFoundError:
            log("    blastn not found — skipping remote BLAST")
            return []
    else:
        log("  Remote BLAST failed after all retries — skipping")
        return []

    # Parse — outfmt has 11 columns (extra stitle), but _parse_blast6 needs 10+
    hits: list[dict] = []
    if blast_out.exists():
        with open(blast_out) as fh:
            for line in fh:
                cols = line.rstrip().split("\t")
                if len(cols) < 10:
                    continue
                aln_len = int(cols[3])
                if aln_len < 30:
                    continue
                # Use stitle (col 10) as a readable subject description
                stitle = cols[10] if len(cols) > 10 else cols[1]
                hits.append({
                    "query": cols[0],
                    "subject": f"{cols[1]}|{stitle}",
                    "identity": float(cols[2]), "length": aln_len,
                    "q_start": int(cols[4]), "q_end": int(cols[5]),
                    "s_start": int(cols[6]), "s_end": int(cols[7]),
                    "evalue": float(cols[8]), "bitscore": float(cols[9]),
                })
        blast_out.unlink(missing_ok=True)

    log(f"  Remote BLAST: {len(hits)} hits from NCBI nt")
    return hits


# ---------------------------------------------------------------------------
# Hit merging
# ---------------------------------------------------------------------------

def _merge_annotations(
    local_hits: list[dict], remote_hits: list[dict],
) -> list[dict]:
    """Merge local + remote hits, keeping best bitscore per query region.

    For each position along the insert, prefer the hit with highest bitscore.
    Local hits from element_db are preferred at equal score (more specific names).
    """
    # Tag source
    for h in local_hits:
        h["source"] = "element_db"
    for h in remote_hits:
        h["source"] = "ncbi_nt"

    all_hits = local_hits + remote_hits
    if not all_hits:
        return []

    # Group by query sequence
    from collections import defaultdict as _dd
    by_query: dict[str, list[dict]] = _dd(list)
    for h in all_hits:
        by_query[h["query"]].append(h)

    merged: list[dict] = []
    for qname, hits in by_query.items():
        # Sort by bitscore descending, prefer element_db on tie
        hits.sort(key=lambda h: (-h["bitscore"],
                                  0 if h["source"] == "element_db" else 1))

        # Greedy interval selection: pick best non-overlapping hits
        # (allow 80% reciprocal overlap to keep significant alternatives)
        selected: list[dict] = []
        covered: list[tuple[int, int]] = []

        for h in hits:
            qs, qe = min(h["q_start"], h["q_end"]), max(h["q_start"], h["q_end"])
            h_len = qe - qs + 1

            # Check overlap with already-selected hits
            dominated = False
            for cs, ce in covered:
                overlap = max(0, min(qe, ce) - max(qs, cs) + 1)
                if overlap > 0.80 * h_len:
                    dominated = True
                    break
            if not dominated:
                selected.append(h)
                covered.append((qs, qe))

        merged.extend(selected)

    merged.sort(key=lambda h: (h["query"], h["q_start"]))
    return merged


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def annotate_insert(
    insert_fasta: Path,
    element_db: Path,
    output_dir: Path,
    sample_name: str,
    no_remote_blast: bool = False,
    extra_dbs: list[Path] | None = None,
) -> tuple[Path, Path]:
    """Annotate insert with local element_db BLAST + remote NCBI nt BLAST.

    Runs both in sequence (local is fast, remote may take 1-5 min),
    then merges by best bitscore per region. Output format is unchanged
    for downstream report generation.

    ``extra_dbs`` mirrors the Phase 1.5 plumbing (common_payload.fa and/or
    per-sample s04b SPAdes contigs): each extra DB is BLASTed separately
    and its hits are concatenated into the local-hit stream before the
    best-bitscore merge, so sample-specific payloads (e.g. bar, AtYUCCA6,
    full T-DNA backbone contigs) surface in ``element_annotation.tsv``
    alongside the shared EUginius catalogue.
    """
    annotation_tsv = output_dir / "element_annotation.tsv"
    border_tsv = output_dir / "border_hits.tsv"

    # ---- Local BLAST (element_db + any extra DBs) ----
    local_hits = _run_local_blast(
        insert_fasta, element_db, output_dir, tag="primary",
    )
    for i, edb in enumerate(extra_dbs or []):
        if edb is None or not edb.exists() or edb.stat().st_size == 0:
            continue
        extra_hits = _run_local_blast(
            insert_fasta, edb, output_dir, tag=f"extra{i}_{edb.stem}",
        )
        local_hits.extend(extra_hits)

    # ---- Remote BLAST (NCBI nt) ----
    if no_remote_blast:
        log("  Remote BLAST skipped (--no-remote-blast)")
        remote_hits = []
    else:
        remote_hits = _run_remote_blast(insert_fasta, output_dir)

    # ---- Merge ----
    merged = _merge_annotations(local_hits, remote_hits)
    log(f"  Merged annotation: {len(merged)} regions "
        f"({sum(1 for h in merged if h['source'] == 'element_db')} local, "
        f"{sum(1 for h in merged if h['source'] == 'ncbi_nt')} remote)")

    # ---- Write annotation TSV ----
    with open(annotation_tsv, "w") as fout:
        fout.write("query\telement\tidentity\tlength\tq_start\tq_end\t"
                   "s_start\ts_end\tevalue\tsource\n")
        for h in merged:
            fout.write(f"{h['query']}\t{h['subject']}\t{h['identity']}\t"
                       f"{h['length']}\t{h['q_start']}\t{h['q_end']}\t"
                       f"{h['s_start']}\t{h['s_end']}\t{h['evalue']}\t"
                       f"{h['source']}\n")

    if merged:
        for h in merged[:15]:
            src = "L" if h["source"] == "element_db" else "R"
            log(f"    {h['q_start']:>6}-{h['q_end']:<6} [{src}] "
                f"{h['subject'][:60]} ({h['identity']:.1f}%, {h['length']}bp)")
    else:
        log("  No BLAST hits found")

    # ---- T-DNA border motif search ----
    _run_border_blast(insert_fasta, output_dir, border_tsv)

    return annotation_tsv, border_tsv


# ---------------------------------------------------------------------------
# T-DNA border motif scan
# ---------------------------------------------------------------------------

#: The two 25-bp imperfect direct repeats that delimit the nopaline
#: (pTiT37 / pTiC58) T-DNA.  Listed in transfer order: VirD2 nicks at the right
#: border first, and the T-strand is transferred RB -> LB.
#:
#: Assignment is sourced from the primary deposits, not inferred:
#:   * GenBank J01825 — DEFINITION "Ti plasmid (from A.tumefaciens, nopaline
#:     strain T37), T-DNA 5' (left) border" — contains the LB 25-mer.
#:   * GenBank J01826 — DEFINITION "...T-DNA 3' (right) border", COMMENT
#:     "right border at base 158" — the RB 25-mer sits at bases 158-182.
#:   (Yadav, Vanderleyden, Bennett, Barnes & Chilton 1982, PNAS 79:6322.)
#: Independently reproduced by two mainstream binary vectors that cite those
#: records: pBIN19 (U09365) and pCAMBIA-1300 (AF234296).  AF234296 annotates
#: `left border repeat from C58 T-DNA` at 6557..6582 and `right border T-DNA
#: repeat` at 298..323 — and ``db/G281_construct.fa``'s first record is
#: byte-identical to AF234296, with the RB 25-mer at 279 and the LB 25-mer at
#: 6,557 exactly.  ``db/Cas9_construct.fa`` carries LB at 1 and RB at 10,142.
#:
#: The octopine T-DNA (pTi15955, X00493) contains neither repeat; its four
#: border repeats are 24-mers with different variable arms.
#:
#: The element DBs cite AF485783.1 (pBI121) for their own RB/LB records; that
#: accession's annotated borders are complement(2454..2478) = RB and
#: complement(8621..8646) = LB, and they equal the two 25-mers below.  Those
#: records used to hold bases 1-25 and ~14,040-14,065 of AF485783.1 instead —
#: a coordinate slip — and were corrected in the same change that added this
#: table (see tests/test_element_db_borders.py).
TDNA_BORDER_REPEATS: tuple[tuple[str, str], ...] = (
    ("TDNA_RB", "TGACAGGATATATTGGCGGGTAAAC"),
    ("TDNA_LB", "TGGCAGGATATATTGTGGTGTAAAC"),
)

#: Minimum aligned length (bp) for an HSP to count as a border repeat.  RB and
#: LB share only the invariant ``CAGGATATAT`` core and a terminal ``TAAAC``, so
#: at ``-word_size 7`` blastn also emits ~9 bp scraps wherever that core occurs
#: (e.g. rice_G281 positions 4,321 and 7,747).  Those are not borders.  Genuine
#: borders — including the cross-hit of one repeat onto the other's locus at
#: ~84% identity — align across the full 25 bp.
BORDER_MIN_ALIGN_LEN = 20


def _write_border_query(path: Path) -> None:
    """Write the border-repeat query FASTA, one record per distinct repeat."""
    with open(path, "w") as fh:
        for label, seq in TDNA_BORDER_REPEATS:
            fh.write(f">{label}\n{seq}\n")


def _dedupe_border_hits(rows: list[list[str]]) -> list[list[str]]:
    """Collapse BLAST HSPs that describe the same border locus on the insert.

    A 25-bp query at ``-word_size 7`` produces several overlapping HSPs per
    locus, and the two repeats are similar enough to cross-hit each other's
    locus.  Left unmerged, ``generate_report`` — which counts raw lines — reports
    several times as many borders as the insert actually carries.

    Rows are outfmt-6 ``qseqid sseqid pident length qstart qend sstart send``.
    Hits are grouped by subject and overlapping subject interval (orientation is
    ignored: a minus-strand hit reports ``sstart > send``), and the highest
    scoring hit — ``pident x length`` — represents the locus.  Ties resolve on
    ``(qseqid, span)`` so the surviving label never depends on the order blastn
    happened to emit its HSPs.
    """
    def _score(row: list[str]) -> float:
        return float(row[2]) * int(row[3])

    def _span(row: list[str]) -> tuple[int, int]:
        s, e = int(row[6]), int(row[7])
        return (min(s, e), max(s, e))

    # Best-first, so the first row claiming a locus is the one that represents it.
    best_first = sorted(rows, key=lambda r: (-_score(r), r[0], _span(r)))

    kept: list[list[str]] = []
    claimed: dict[str, list[tuple[int, int]]] = {}
    for row in best_first:
        subject = row[1]
        lo, hi = _span(row)
        spans = claimed.setdefault(subject, [])
        if any(lo <= s_hi and hi >= s_lo for s_lo, s_hi in spans):
            continue
        spans.append((lo, hi))
        kept.append(row)

    return sorted(kept, key=lambda r: (r[1], _span(r)[0]))


def _run_border_blast(insert_fasta: Path, output_dir: Path,
                      border_tsv: Path) -> Path:
    """BLAST both T-DNA border repeats against the insert; write distinct loci.

    Writes ``border_tsv`` in outfmt-6 (>= 8 columns, as
    ``scripts/viz/plot_insert_structure.py`` requires), one row per border
    locus.  An empty file means no border motif was found.
    """
    border_fa = output_dir / "_borders.fa"
    _write_border_query(border_fa)
    raw_tsv = output_dir / "_border_hits_raw.tsv"

    try:
        result = subprocess.run(
            ["blastn", "-query", str(border_fa), "-subject", str(insert_fasta),
             "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send",
             "-evalue", "1", "-word_size", "7",
             "-out", str(raw_tsv)],
            stderr=subprocess.PIPE, text=True,
        )
        if result.returncode != 0:
            stderr_tail = (result.stderr or "")[-400:]
            log(f"  Border scan: blastn failed (rc={result.returncode}); "
                f"no borders reported. stderr: {stderr_tail}")
            border_tsv.write_text("")
            return border_tsv

        rows: list[list[str]] = []
        if raw_tsv.exists():
            for line in raw_tsv.read_text().splitlines():
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8:
                    continue
                try:
                    if int(cols[3]) < BORDER_MIN_ALIGN_LEN:
                        continue
                except ValueError:
                    continue
                rows.append(cols)

        deduped = _dedupe_border_hits(rows)
        border_tsv.write_text(
            "".join("\t".join(row) + "\n" for row in deduped)
        )

        if deduped:
            loci = ", ".join(f"{r[0]}@{min(int(r[6]), int(r[7]))}" for r in deduped)
            log(f"  T-DNA border motifs found: {len(deduped)} locus/loci ({loci})")
        else:
            log("  No border motifs found (may need manual inspection)")
    finally:
        border_fa.unlink(missing_ok=True)
        raw_tsv.unlink(missing_ok=True)

    return border_tsv
