# s05 Session 3 — Annotation Module Extraction (Issue #4)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract Phase 4 (Annotation) — 5 functions plus their `# Phase 4: Annotation` section header — from `scripts/s05_insert_assembly.py` into `scripts/s05/annotation.py`. Single commit per spec §4 Session 3.

**Architecture:** `annotation.py` lives at DAG stage 4 (already pre-registered in `tests/test_s05_import_dag.py`). It depends on:
- `scripts/s05/primitives.py` (stage 0) — `log`
- `scripts/s05/classify.py` (stage 1) — `_parse_src_tag` (used by `_parse_blast6` to strip `|src=<tag>` suffixes from BLAST sseqids)

To make `_parse_src_tag` accessible to `annotation.py`, the existing `classify.py` shim is extended to also re-export it from the monolith. This is consistent with the shim's role as a stage-1 surface; the symbol stays defined in the monolith until Session 4 fully extracts classify.

The legacy `annotate_report.py` shim at stage 4 is updated to source `annotate_insert` from the new `annotation` module while continuing to re-export `generate_report` and `write_stats` from the monolith (those belong to Session 6: report.py).

**Tech Stack:** Python 3.11, pytest 9, micromamba env `redgene`. BLAST+ subprocess. No new deps.

---

## Reference: Session 1/2 Pattern

Each filter module (`filter_a_host.py`, `filter_b_flank.py`, `filter_c_chimeric.py`, `filter_d_altlocus.py`):
1. Verbatim function copy (same algorithm, same constants, same comments)
2. Only allowed s05 sibling imports follow the DAG (stage X imports from stages 0..X)
3. Module docstring explains rationale + cites verdict.py rule when relevant
4. `from __future__ import annotations`
5. Threshold constants (if any) mirrored from monolith verbatim
6. Monolith re-imports symbols in BOTH `try` / `except ImportError` branches; original definition replaced with a 3-line tombstone comment

Session 3 follows the same pattern, with one extension: it also extends the existing `classify.py` shim's re-export list (mechanical, no behavior change).

---

## File Map

| File | Action | Purpose |
|------|--------|---------|
| `scripts/s05/annotation.py` | **Create** | 5 functions: `_parse_blast6`, `_run_local_blast`, `_run_remote_blast`, `_merge_annotations`, `annotate_insert` |
| `scripts/s05/classify.py` | **Modify** (shim) | Add `_parse_src_tag` to existing re-export list (1-line change) |
| `scripts/s05/annotate_report.py` | **Modify** (shim) | Source `annotate_insert` from `.annotation`; keep `generate_report`/`write_stats` from monolith |
| `scripts/s05/__init__.py` | **Modify** | Re-export the 5 new annotation symbols |
| `scripts/s05_insert_assembly.py` | **Modify** | Add new imports to BOTH ImportError branches; delete original 5 function definitions; add tombstone comment block |

---

## Task 1: Extract Annotation Module

**Files:** as listed above.

- [ ] **Step 1: Verify pre-state — pytest baseline + DAG test + git HEAD**

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"
git log --oneline -3
pytest tests/ -q
pytest tests/test_s05_import_dag.py -v
pytest tests/test_extra_element_db.py::test_annotate_insert_uses_extra_dbs -v
```

Expected:
- `git log` HEAD = `6658fb8 Issue #4 Session 2 [2/2]: extract filter_d_altlocus.py ...` (Session 2 close-out)
- pytest: **237 passed, 1 skipped**
- DAG: **12 PASSED**
- BLAST smoke: **1 PASSED** (this single test will be the Session 3 exit-gate "BLAST path smoke")

If any check disagrees, stop and report — Session 2 may not have landed cleanly.

- [ ] **Step 2: Capture original function bodies for byte-identity verification**

Save the pre-edit monolith into a tempfile so the implementer can `diff` after the move. Use the current commit (HEAD) as the source of truth:

```bash
git show HEAD:scripts/s05_insert_assembly.py > /tmp/s05_pre_session3.py
wc -l /tmp/s05_pre_session3.py
# Expected: 3879 lines
```

Then extract the four annotation function bodies (line ranges from the current HEAD):
- `_parse_blast6`: monolith lines 2797–2828
- `_run_local_blast`: lines 2829–2862 (look for blank lines as boundary)
- `_run_remote_blast`: lines 2864–2919
- `_merge_annotations`: lines 2920–2974
- `annotate_insert`: lines 2976–3083

Use `awk 'NR>=START && NR<=END' /tmp/s05_pre_session3.py > /tmp/<func>.py` to capture each block. The implementer should use these for byte-identity checks during Step 3 self-review.

- [ ] **Step 3: Create `scripts/s05/annotation.py` with verbatim functions**

Create `scripts/s05/annotation.py` with this EXACT structure. Each function body must be byte-identical to its source in the monolith — preserve every comment, every blank line, every whitespace pattern.

```python
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
from collections import defaultdict
from pathlib import Path

from .classify import _parse_src_tag
from .primitives import log


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
    subprocess.run(
        ["makeblastdb", "-in", str(element_db), "-dbtype", "nucl",
         "-out", str(db_prefix)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
    )
    blast_out = output_dir / f"_local_blast_{suffix}.tsv"
    subprocess.run(
        ["blastn", "-query", str(insert_fasta), "-db", str(db_prefix),
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
         "-evalue", "1e-5", "-max_target_seqs", "50",
         "-out", str(blast_out)],
        stderr=subprocess.DEVNULL, check=True,
    )
    hits = _parse_blast6(blast_out)
    # Cleanup
    blast_out.unlink(missing_ok=True)
    for ext in [".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto", ".njs"]:
        Path(f"{db_prefix}{ext}").unlink(missing_ok=True)
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
    border_fa = output_dir / "_borders.fa"
    with open(border_fa, "w") as fh:
        fh.write(">RB_consensus\nTGGCAGGATATATTGTGGTGTAAAC\n")
        fh.write(">LB_consensus\nTGGCAGGATATATTGTGGTGTAAAC\n")

    subprocess.run(
        ["blastn", "-query", str(border_fa), "-subject", str(insert_fasta),
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send",
         "-evalue", "1", "-word_size", "7",
         "-out", str(border_tsv)],
        stderr=subprocess.DEVNULL,
    )

    if border_tsv.exists() and border_tsv.stat().st_size > 0:
        log("  T-DNA border motifs found")
    else:
        log("  No border motifs found (may need manual inspection)")

    border_fa.unlink(missing_ok=True)
    return annotation_tsv, border_tsv
```

**Important byte-identity notes:**
- `_parse_blast6` was originally calling `_parse_src_tag` as a bare name (because both functions lived in the monolith). In the new module, `_parse_src_tag` is imported from `.classify` at module top — the call site in the function body remains the bare name `_parse_src_tag(subject, default="")`.
- `_merge_annotations` retains its pre-existing `from collections import defaultdict as _dd` inside the function body (a minor redundancy with the module-level `from collections import defaultdict` import). This is **deliberate verbatim preservation** — do NOT clean it up.
- The module-level `from collections import defaultdict` is harmless because no top-level code uses it; if it offends, leave it: matches Session 1/2 verbatim policy.

Also note: `defaultdict` at the top is technically unused at module level (only `_dd` inside `_merge_annotations` uses it). To keep the byte-identity rule strict, you may either include it (mirrors that the original monolith imports `defaultdict` at top — verify this with `grep "from collections import defaultdict" /tmp/s05_pre_session3.py`) or omit it (the original may not import it at module level — if so, omit it and let the function-local import stand alone). **Run the grep first; mirror what the monolith actually has.** Whichever path you take, the function bodies must remain byte-identical.

- [ ] **Step 4: Update `scripts/s05/classify.py` shim — add `_parse_src_tag`**

Open `scripts/s05/classify.py`. Current content:

```python
"""Phase 1.5 classify re-exports (shim, v1.0)."""
from scripts.s05_insert_assembly import (  # noqa: F401
    classify_site_tiers,
    _batch_check_element_hits,
    _should_replace,
    _filter_host_endogenous,
)
```

Add `_parse_src_tag` to the re-export list (alphabetize for cleanliness):

```python
"""Phase 1.5 classify re-exports (shim, v1.0).

Session 3 (2026-04-28): added ``_parse_src_tag`` to the re-export list so
``scripts/s05/annotation.py`` can consume it without crossing back into
the monolith. Full classify extraction is scheduled for Session 4 per
``docs/superpowers/specs/2026-04-19-s05-module-split-design.md`` §4.
"""
from scripts.s05_insert_assembly import (  # noqa: F401
    classify_site_tiers,
    _batch_check_element_hits,
    _filter_host_endogenous,
    _parse_src_tag,
    _should_replace,
)
```

The DAG test will continue to pass: `classify.py` does NOT import from any s05 sibling, only from the (out-of-package) monolith.

- [ ] **Step 5: Update `scripts/s05/annotate_report.py` shim**

Current:
```python
"""Phase 4 annotate + report re-exports (shim, v1.0)."""
from scripts.s05_insert_assembly import (  # noqa: F401
    annotate_insert,
    generate_report,
    write_stats,
)
```

Source `annotate_insert` from the new module; keep `generate_report` / `write_stats` from monolith (Session 6 territory):

```python
"""Phase 4 annotate + report re-exports (shim, v1.0).

Session 3 (2026-04-28): ``annotate_insert`` now lives in
``scripts/s05/annotation.py``; this shim sources it from there while
``generate_report`` and ``write_stats`` continue to ride on the monolith
until Session 6 carves out ``report.py``.
"""
from .annotation import annotate_insert  # noqa: F401
from scripts.s05_insert_assembly import (  # noqa: F401
    generate_report,
    write_stats,
)
```

The DAG: `annotate_report` is stage 4, `annotation` is stage 4 — same-stage import is allowed (test rule is strict `>`).

- [ ] **Step 6: Update `scripts/s05/__init__.py` — re-export the 5 annotation symbols**

After the existing `from .filter_d_altlocus import (...)` block, add:

```python
from .annotation import (  # noqa: F401
    _parse_blast6,
    _run_local_blast,
    _run_remote_blast,
    _merge_annotations,
    annotate_insert,
)
```

Update the docstring (currently mentions Sessions 1 + 2) to add a Session 3 line:

```
Session 3 (2026-04-28) added ``annotation`` (Phase 4 BLAST-driven
element annotation), retiring the ``annotate_report.py`` shim's
re-export of ``annotate_insert``.
```

- [ ] **Step 7: Wire monolith — add imports to BOTH ImportError branches**

In `scripts/s05_insert_assembly.py`, locate the existing `try: ... except ImportError:` block at the top (lines ~47–125 after Session 2). Add a new import group AFTER the `filter_d_altlocus` block in BOTH branches:

```python
    from s05.annotation import (
        _parse_blast6,
        _run_local_blast,
        _run_remote_blast,
        _merge_annotations,
        annotate_insert,
    )
```

Both branches must be modified identically (one inside `try:`, one inside `except ImportError:`).

- [ ] **Step 8: Delete original 5 functions from monolith + add tombstone**

Delete the following blocks from `scripts/s05_insert_assembly.py`:
- The `# Phase 4: Annotation` section header (3 lines, currently at ~2792–2794)
- `def _parse_blast6(...)` (currently lines ~2797–2828)
- `def _run_local_blast(...)` (currently lines ~2829–2862)
- `def _run_remote_blast(...)` (currently lines ~2864–2919)
- `def _merge_annotations(...)` (currently lines ~2920–2974)
- `def annotate_insert(...)` (currently lines ~2976–3083)

Replace all six deletions (one section header + 5 functions) with a single tombstone block placed where the section header used to be:

```python
# ---------------------------------------------------------------------------
# Phase 4 (Annotation) functions moved to scripts/s05/annotation.py
# (Issue #4 Session 3) and re-imported at the top of this module for
# backward compatibility with existing callers. The monolith's
# ``classify_site_tiers`` continues to consume the resulting hits.
# ---------------------------------------------------------------------------
```

- [ ] **Step 9: Run DAG test (must show 13 PASSED)**

```bash
pytest tests/test_s05_import_dag.py -v
```

Expected output: **13 PASSED**. The new `[annotation.py]` parametrize entry must appear and PASS. If `annotation.py` triggers a DAG violation, check that its only s05-internal imports are `from .primitives import log` and `from .classify import _parse_src_tag` — nothing else from the package.

If the test fails because `classify.py` now has a violation: rerun the test with verbose output and confirm the message points at `classify.py` (it shouldn't — classify only imports from the monolith, not from s05 siblings). If it does, you've accidentally introduced an `from .annotation import ...` in `classify.py` (which would be a cycle — annotation imports from classify too).

- [ ] **Step 10: Run BLAST path smoke (per spec exit gate)**

```bash
pytest tests/test_extra_element_db.py::test_annotate_insert_uses_extra_dbs -v
```

Expected: **1 passed**. This test creates real FASTA inputs, calls `s05.annotate_insert(...)` with `no_remote_blast=True`, and verifies the resulting `element_annotation.tsv` contains both `elem_A` (primary DB hit) and `bar` (extra DB hit). Failure here is a **real regression** — investigate immediately.

- [ ] **Step 11: Run full pytest + verdict snapshots**

```bash
pytest tests/ -q
pytest tests/test_verdict_snapshots.py -v
```

Expected:
- Full suite: **238 passed, 1 skipped** (237 baseline + 1 new DAG parametrize entry for `annotation.py`)
- Verdict snapshots: 26 PASSED with zero drift

If pytest goes higher than 238 (e.g., 239), check whether you accidentally added a new test file. If it goes lower than 238, snapshot drift or a real regression — investigate.

- [ ] **Step 12: Verify CLI smoke**

```bash
python scripts/s05_insert_assembly.py --help > /dev/null && echo OK
```

Expected: prints `OK` and exits 0.

- [ ] **Step 13: Commit (on `main` branch)**

```bash
git add scripts/s05/annotation.py scripts/s05/__init__.py \
        scripts/s05/classify.py scripts/s05/annotate_report.py \
        scripts/s05_insert_assembly.py
git commit -m "$(cat <<'EOF'
Issue #4 Session 3: extract annotation.py from s05 monolith

Phase 4 annotation pipeline (5 functions: _parse_blast6,
_run_local_blast, _run_remote_blast, _merge_annotations, annotate_insert)
moves verbatim into scripts/s05/annotation.py. Monolith re-imports all
five symbols so existing callers (classify_site_tiers consumes the
resulting hits, generate_report writes them to element_annotation.tsv)
remain bit-identical.

annotation.py imports `_parse_src_tag` from `.classify`; the existing
classify.py shim is extended to re-export it (mechanical, 1-symbol
addition). Full classify extraction is Session 4 per spec.

annotate_report.py shim now sources `annotate_insert` from
`.annotation`; `generate_report` / `write_stats` continue to ride on
the monolith pending Session 6 (report.py extraction).

DAG: annotation pre-registered at stage 4 in tests/test_s05_import_dag.py;
parametrized DAG test picks it up automatically.

Tests: 237 PASS + 1 skipped → 238 PASS + 1 skipped (DAG parametrize +1).
BLAST path smoke (test_annotate_insert_uses_extra_dbs) green. No
verdict logic touched; 26 verdict snapshots unchanged.

Per docs/superpowers/specs/2026-04-19-s05-module-split-design.md §4 Session 3.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 14: Final exit-gate verification**

```bash
pytest tests/ -q
pytest tests/test_s05_import_dag.py -v
pytest tests/test_extra_element_db.py -v
python scripts/s05_insert_assembly.py --help > /dev/null && echo "CLI OK"
git log --oneline -3
wc -l scripts/s05_insert_assembly.py scripts/s05/*.py
```

Expected:
- pytest: 238 PASS + 1 skipped
- DAG: 13 PASSED
- test_extra_element_db: all PASSED (BLAST smoke green)
- CLI: `CLI OK`
- monolith line count: drops by ~290 lines (5 functions + section header) → ~3589 lines
- new module: `scripts/s05/annotation.py` ~290 lines

---

## Risks and rollback

| Risk | Mitigation |
|------|------------|
| `_parse_src_tag` cyclic import (annotation → classify, classify → monolith → ?) | classify.py shim only does `from scripts.s05_insert_assembly import ...` — that's an out-of-package import, not subject to DAG. No cycle. |
| `importlib.util` test loading breaks (test_extra_element_db loads monolith as `s05`) | The pattern already works in Session 1/2 (filters use the same trick). Session 3 just adds another import group; no new mechanism. |
| Snapshot drift | All 5 functions are verbatim copies. If snapshots drift, the move was not truly verbatim — revert and diff. |
| BLAST smoke fails | The `_parse_src_tag` indirection through classify shim is the only behavioral change. If smoke fails, check that `classify.py` actually exposes `_parse_src_tag` (test: `python -c "from scripts.s05.classify import _parse_src_tag; print(_parse_src_tag('a\|src=b','default'))"` → `('a','b')`). |
| `defaultdict` import drift | Step 3 explicitly says to mirror what the monolith currently has. Run the `grep` before deciding whether to add the module-level import. |

Rollback: single commit. `git revert <sha>` or `git reset --hard HEAD~1` returns to baseline (237 + 1 skipped, DAG 12).
