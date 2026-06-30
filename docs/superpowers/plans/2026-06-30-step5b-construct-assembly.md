# Step 5b — Construct Assembly Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in pipeline sub-step `5b` that reconstructs and characterizes the inserted transgene construct from Illumina reads — a global classified/annotated contig inventory plus per-site links to s05 CANDIDATE insertion sites.

**Architecture:** A single standalone script `scripts/s05b_construct_assembly.py` reuses s04b's pre-filter SPAdes contigs (`contigs_all.fasta`), classifies each contig by host/element BLAST coverage (assemble-then-classify), annotates non-host contigs (local element BLAST always; opt-in remote nt for unknowns), and links contigs to s05 sites via minimap2. Wired into `run_pipeline.py` after step 5. The step-8 PDF gains a "Construct Inventory" section.

**Tech Stack:** Python 3.11, BLAST+ (`blastn`), minimap2, SPAdes (fallback), pysam-free; pytest for tests. Follows `scripts/s04b_construct_assembly.py` conventions (argparse CLI, `pathlib.Path`, stderr logging, exit-0 contract).

## Global Constraints

- Python 3.11, type hints preferred; all paths `pathlib.Path`; logging to stderr.
- **Exit-0 contract:** the step never raises on missing/empty inputs; it always writes the expected output files (empty FASTA + header-only TSVs) so downstream (step 8) can find them. Return code 0 in all input-shortfall paths.
- **No s05 behavior change:** the existing 26 verdict snapshot fixtures in `tests/fixtures/verdict_snapshots/` MUST remain byte-identical. 5b reads s05 output, never writes into `s05_insert_assembly/`.
- Concrete thresholds (module constants, exact values): `MIN_CONTIG_LEN = 200`, `HOST_PIDENT = 90.0`, `ELEMENT_PIDENT = 80.0`, `CHIMERIC_MIN = 0.20`, `ELEMENT_MIN = 0.50`, `HOST_MIN = 0.70`, `LINK_SLOP = 100`.
- Classification priority (first match wins): `chimeric` (both fracs ≥ CHIMERIC_MIN) → `construct` (element_frac ≥ ELEMENT_MIN) → `host` (host_frac ≥ HOST_MIN) → `foreign_unknown`.
- Output dir: `results/{sample}/s05b_construct_assembly/`.
- The pipeline step dict in `run_pipeline.py` is named `STEP_SCRIPTS` (not STEP_MAP); the order list is `STEP_ORDER`.
- Environment for running tests: `export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"` then `pytest`.
- Spec: `docs/superpowers/specs/2026-06-30-step5b-construct-assembly-design.md`.

---

## File Structure

| File | Responsibility |
|------|----------------|
| `scripts/s05b_construct_assembly.py` (create) | The 5b step: pure helpers (interval merge, coverage, classify, blast6/PAF parse, s05 site parse, linkage), subprocess wrappers (blastn, minimap2, SPAdes fallback), output writers, `run()`/argparse |
| `tests/test_s05b_construct_assembly.py` (create) | Unit tests for all pure helpers + writers + exit-0 contract; one blastn integration smoke guarded by `shutil.which` |
| `run_pipeline.py` (modify) | Add `"5b"` to `STEP_ORDER` + `STEP_SCRIPTS`, `build_step_cmd` branch, `--steps` help |
| `tests/test_run_pipeline_step5b.py` (create) | Dispatch test: `build_step_cmd("5b", …)` argv + `parse_steps` includes 5b |
| `scripts/reports/insertion_pdf.py` (modify) | "Construct Inventory" PDF section (present/absent paths) |
| `tests/test_insertion_pdf_construct.py` (create) | PDF section renders with 5b output; omitted when absent |
| `scripts/util/ground_truth_baseline.py` (move target) | Relocated benchmark tool (from `scripts/s09_ground_truth_baseline.py`) |
| `CLAUDE.md`, `docs/measurements/coverage_sensitivity.md` (modify) | Document 5b; fix moved-tool references |

---

## Task 1: Pure helpers — interval merge, coverage fraction, blast6 parse

**Files:**
- Create: `scripts/s05b_construct_assembly.py`
- Test: `tests/test_s05b_construct_assembly.py`

**Interfaces:**
- Produces:
  - `_merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]` — merges overlapping/adjacent 1-based inclusive intervals, sorted.
  - `_coverage_fraction(hits: list[tuple[int, int]], length: int) -> float` — merged covered bp / `length`; `0.0` if `length <= 0`.
  - `_parse_blast6(text: str, min_pident: float) -> dict[str, list[tuple[int, int]]]` — maps `qseqid` → list of `(qstart, qend)` (normalized so start ≤ end) for rows with `pident >= min_pident`. Expects outfmt `6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_s05b_construct_assembly.py
import importlib.util
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "s05b", Path(__file__).resolve().parent.parent / "scripts" / "s05b_construct_assembly.py"
)
s05b = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(s05b)


def test_merge_intervals_overlap_and_adjacent():
    assert s05b._merge_intervals([(1, 5), (4, 8), (10, 12), (12, 15)]) == [(1, 8), (10, 15)]


def test_merge_intervals_empty():
    assert s05b._merge_intervals([]) == []


def test_coverage_fraction_basic():
    # covered bp = (1..8)=8 + (10..15)=6 = 14 over length 20 → 0.7
    assert s05b._coverage_fraction([(1, 5), (4, 8), (10, 15)], 20) == pytest.approx(0.7)


def test_coverage_fraction_zero_length():
    assert s05b._coverage_fraction([(1, 5)], 0) == 0.0


def test_parse_blast6_filters_by_pident():
    text = (
        "c1\thostX\t95.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t180\n"
        "c1\thostX\t70.0\t50\t0\t0\t200\t249\t1\t50\t1e-10\t90\n"   # below 90 → dropped
        "c2\thostX\t99.0\t30\t0\t0\t60\t31\t1\t30\t1e-20\t60\n"      # qstart>qend → normalized
    )
    hits = s05b._parse_blast6(text, min_pident=90.0)
    assert hits["c1"] == [(1, 100)]
    assert hits["c2"] == [(31, 60)]
    assert "hostX" not in hits  # keyed by qseqid only
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_s05b_construct_assembly.py -q`
Expected: FAIL — `FileNotFoundError`/`ModuleNotFoundError` (script doesn't exist yet).

- [ ] **Step 3: Write minimal implementation**

Create `scripts/s05b_construct_assembly.py`:

```python
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
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_s05b_construct_assembly.py -q`
Expected: PASS (5 tests).

- [ ] **Step 5: Commit**

```bash
git add scripts/s05b_construct_assembly.py tests/test_s05b_construct_assembly.py
git commit -m "Step 5b [1/10]: pure helpers (interval merge, coverage, blast6 parse)"
```

---

## Task 2: Contig classification

**Files:**
- Modify: `scripts/s05b_construct_assembly.py`
- Test: `tests/test_s05b_construct_assembly.py`

**Interfaces:**
- Produces: `_classify_contig(host_frac: float, element_frac: float) -> str` returning one of `"chimeric"`, `"construct"`, `"host"`, `"foreign_unknown"` by the priority rule.

- [ ] **Step 1: Write the failing test**

```python
def test_classify_chimeric_both_present():
    assert s05b._classify_contig(0.5, 0.5) == "chimeric"
    # boundary: both exactly at CHIMERIC_MIN
    assert s05b._classify_contig(0.20, 0.20) == "chimeric"


def test_classify_construct_element_dominant():
    # element high, host below chimeric floor → construct
    assert s05b._classify_contig(0.10, 0.80) == "construct"
    assert s05b._classify_contig(0.0, 0.50) == "construct"  # boundary ELEMENT_MIN


def test_classify_host_dominant():
    assert s05b._classify_contig(0.95, 0.0) == "host"
    assert s05b._classify_contig(0.70, 0.10) == "host"  # boundary HOST_MIN, element below chimeric


def test_classify_foreign_unknown():
    assert s05b._classify_contig(0.10, 0.10) == "foreign_unknown"
    assert s05b._classify_contig(0.0, 0.0) == "foreign_unknown"
    # element present but below ELEMENT_MIN, host below HOST_MIN, not both >= chimeric
    assert s05b._classify_contig(0.10, 0.30) == "foreign_unknown"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_s05b_construct_assembly.py -k classify -q`
Expected: FAIL — `AttributeError: module 's05b' has no attribute '_classify_contig'`.

- [ ] **Step 3: Write minimal implementation**

Add to `scripts/s05b_construct_assembly.py` (after `_parse_blast6`):

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_s05b_construct_assembly.py -k classify -q`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add scripts/s05b_construct_assembly.py tests/test_s05b_construct_assembly.py
git commit -m "Step 5b [2/10]: contig classification with concrete thresholds"
```

---

## Task 3: BLAST coverage runners (host_frac / element_frac per contig)

**Files:**
- Modify: `scripts/s05b_construct_assembly.py`
- Test: `tests/test_s05b_construct_assembly.py`

**Interfaces:**
- Consumes: `_parse_blast6`, `_coverage_fraction` (Task 1); `read_fasta` from `scripts/s05/primitives.py`.
- Produces:
  - `_run_blastn(query: Path, subjects: list[Path], min_pident: float) -> dict[str, list[tuple[int, int]]]` — runs `blastn -subject` against each FASTA in `subjects`, pools and returns merged `qseqid -> [(qstart,qend)]`. Missing/empty subjects skipped. blastn failure on a subject → warn + skip (returns what succeeded).
  - `_contig_lengths(contigs_fa: Path) -> dict[str, int]` — contig_id → length via `read_fasta`.

- [ ] **Step 1: Write the failing test** (integration smoke, guarded by tool availability)

```python
import shutil

requires_blastn = pytest.mark.skipif(
    shutil.which("blastn") is None, reason="blastn not on PATH"
)


@requires_blastn
def test_run_blastn_self_match(tmp_path):
    seq = "ACGT" * 80  # 320 bp
    q = tmp_path / "q.fa"
    q.write_text(f">c1\n{seq}\n")
    subj = tmp_path / "s.fa"
    subj.write_text(f">ref1\n{seq}\n")
    hits = s05b._run_blastn(q, [subj], min_pident=90.0)
    assert "c1" in hits
    frac = s05b._coverage_fraction(hits["c1"], len(seq))
    assert frac > 0.9


def test_run_blastn_skips_missing_subject(tmp_path):
    q = tmp_path / "q.fa"
    q.write_text(">c1\nACGTACGTAC\n")
    missing = tmp_path / "nope.fa"
    hits = s05b._run_blastn(q, [missing], min_pident=90.0)
    assert hits == {}  # nothing to match, no crash


def test_contig_lengths(tmp_path):
    fa = tmp_path / "contigs.fa"
    fa.write_text(">c1\nACGTACGT\n>c2\nAAA\n")
    assert s05b._contig_lengths(fa) == {"c1": 8, "c2": 3}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_s05b_construct_assembly.py -k "blastn or contig_lengths" -q`
Expected: FAIL — `_run_blastn` / `_contig_lengths` not defined.

- [ ] **Step 3: Write minimal implementation**

Add imports and helpers to `scripts/s05b_construct_assembly.py`. At the top, after the stdlib imports:

```python
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))
from s05.primitives import read_fasta  # noqa: E402  (pure FASTA reader, no side effects)
```

Then add:

```python
def _is_nonempty_file(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


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
    return {cid: len(seq) for cid, seq in read_fasta(contigs_fa).items()}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_s05b_construct_assembly.py -k "blastn or contig_lengths" -q`
Expected: PASS (the self-match test runs when blastn is present, else skips).

- [ ] **Step 5: Commit**

```bash
git add scripts/s05b_construct_assembly.py tests/test_s05b_construct_assembly.py
git commit -m "Step 5b [3/10]: blastn coverage runners + contig lengths"
```

---

## Task 4: Parse s05 CANDIDATE sites

**Files:**
- Modify: `scripts/s05b_construct_assembly.py`
- Test: `tests/test_s05b_construct_assembly.py`

**Interfaces:**
- Produces:
  - `Site` — a `typing.NamedTuple` with fields `site_id: str`, `chrom: str`, `pos: int`, `verdict: str`.
  - `_parse_s05_sites(s05_dir: Path, verdicts: set[str] | None = None) -> list[Site]` — scans `insertion_*_report.txt`, derives `site_id` from the filename (`insertion_<chrom>_<pos>_report.txt` → site_id `<chrom>_<pos>`), reads the `Verdict:` line, filters to `verdicts` (default `{"CANDIDATE"}`). Missing dir → `[]`.

- [ ] **Step 1: Write the failing test**

```python
def test_parse_s05_sites_filters_candidate(tmp_path):
    s05 = tmp_path / "s05_insert_assembly"
    s05.mkdir()
    (s05 / "insertion_Chr3_16439674_report.txt").write_text(
        "Insertion site: Chr3:16439674\nVerdict: CANDIDATE\n"
    )
    (s05 / "insertion_Chr3_31434949_report.txt").write_text(
        "Insertion site: Chr3:31434949\nVerdict: FALSE_POSITIVE\n"
    )
    (s05 / "insertion_SLM_r2.0ch08_65107378_report.txt").write_text(
        "Verdict: CANDIDATE\n"
    )
    sites = s05b._parse_s05_sites(s05)
    ids = {s.site_id for s in sites}
    assert ids == {"Chr3_16439674", "SLM_r2.0ch08_65107378"}
    by_id = {s.site_id: s for s in sites}
    assert by_id["Chr3_16439674"].chrom == "Chr3"
    assert by_id["Chr3_16439674"].pos == 16439674
    assert by_id["SLM_r2.0ch08_65107378"].chrom == "SLM_r2.0ch08"
    assert by_id["SLM_r2.0ch08_65107378"].pos == 65107378


def test_parse_s05_sites_missing_dir(tmp_path):
    assert s05b._parse_s05_sites(tmp_path / "absent") == []
```

Note the second fixture proves chrom names may themselves contain underscores
(`SLM_r2.0ch08`), so `chrom`/`pos` split must be on the **last** underscore.

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_s05b_construct_assembly.py -k parse_s05 -q`
Expected: FAIL — `_parse_s05_sites` not defined.

- [ ] **Step 3: Write minimal implementation**

Add `from typing import NamedTuple` to the imports, then:

```python
class Site(NamedTuple):
    site_id: str
    chrom: str
    pos: int
    verdict: str


def _parse_s05_sites(s05_dir: Path, verdicts: set[str] | None = None) -> list[Site]:
    """Parse insertion_*_report.txt into Site rows, filtered by verdict.

    Filename form: insertion_<chrom>_<pos>_report.txt ; <chrom> may contain
    underscores, so site_id = stuff between 'insertion_' and '_report', and
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_s05b_construct_assembly.py -k parse_s05 -q`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add scripts/s05b_construct_assembly.py tests/test_s05b_construct_assembly.py
git commit -m "Step 5b [4/10]: parse s05 CANDIDATE sites (last-underscore split)"
```

---

## Task 5: PAF parse + per-site linkage

**Files:**
- Modify: `scripts/s05b_construct_assembly.py`
- Test: `tests/test_s05b_construct_assembly.py`

**Interfaces:**
- Consumes: `Site` (Task 4).
- Produces:
  - `Aln` — `NamedTuple(contig_id: str, chrom: str, tstart: int, tend: int)` (0-based PAF target coords, half-open as emitted by minimap2).
  - `_parse_paf(text: str) -> list[Aln]` — PAF cols: qname(0) … tname(5) tstart(7) tend(8).
  - `Link` — `NamedTuple(site_id: str, verdict: str, contig_id: str, junction_side: str, link_confidence: str)`.
  - `_link_contigs_to_sites(alns: list[Aln], sites: list[Site], contig_class: dict[str, str], slop: int = LINK_SLOP) -> list[Link]` — for each (aln, site) where `aln.chrom == site.chrom` and the aln start or end is within `slop` of `site.pos`, emit a Link. `junction_side = "downstream"` if `|tstart - pos| <= slop`, else `"upstream"` (aln end near pos). `link_confidence = "high"` if `contig_class[contig] == "chimeric"` else `"medium"`. Deterministic order: sites then contigs as given. De-duplicate identical (site_id, contig_id) keeping the first.

- [ ] **Step 1: Write the failing test**

```python
def test_parse_paf_minimal():
    paf = "c1\t300\t0\t300\t+\tChr3\t40000000\t16439600\t16439900\t280\t300\t60\n"
    alns = s05b._parse_paf(paf)
    assert len(alns) == 1
    assert alns[0].contig_id == "c1"
    assert alns[0].chrom == "Chr3"
    assert alns[0].tstart == 16439600
    assert alns[0].tend == 16439900


def test_link_contigs_downstream_and_confidence():
    sites = [s05b.Site("Chr3_16439674", "Chr3", 16439674, "CANDIDATE")]
    # aln tstart=16439700 within 100 of pos 16439674 -> downstream
    alns = [s05b.Aln("c1", "Chr3", 16439700, 16450000)]
    links = s05b._link_contigs_to_sites(alns, sites, {"c1": "chimeric"})
    assert len(links) == 1
    assert links[0].site_id == "Chr3_16439674"
    assert links[0].contig_id == "c1"
    assert links[0].junction_side == "downstream"
    assert links[0].link_confidence == "high"


def test_link_contigs_upstream_medium():
    sites = [s05b.Site("Chr3_16439674", "Chr3", 16439674, "CANDIDATE")]
    # aln tend=16439700 within 100 of pos -> upstream; construct class -> medium
    alns = [s05b.Aln("c2", "Chr3", 16430000, 16439700)]
    links = s05b._link_contigs_to_sites(alns, sites, {"c2": "construct"})
    assert links[0].junction_side == "upstream"
    assert links[0].link_confidence == "medium"


def test_link_contigs_no_match_wrong_chrom_or_far():
    sites = [s05b.Site("Chr3_16439674", "Chr3", 16439674, "CANDIDATE")]
    alns = [
        s05b.Aln("c3", "Chr9", 16439700, 16450000),       # wrong chrom
        s05b.Aln("c4", "Chr3", 17000000, 17000300),       # far away
    ]
    assert s05b._link_contigs_to_sites(alns, sites, {"c3": "construct", "c4": "construct"}) == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_s05b_construct_assembly.py -k "paf or link" -q`
Expected: FAIL — `_parse_paf`/`Aln`/`_link_contigs_to_sites` not defined.

- [ ] **Step 3: Write minimal implementation**

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_s05b_construct_assembly.py -k "paf or link" -q`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add scripts/s05b_construct_assembly.py tests/test_s05b_construct_assembly.py
git commit -m "Step 5b [5/10]: PAF parse + per-site contig linkage"
```

---

## Task 6: Output writers

**Files:**
- Modify: `scripts/s05b_construct_assembly.py`
- Test: `tests/test_s05b_construct_assembly.py`

**Interfaces:**
- Consumes: `Link` (Task 5), `read_fasta` (Task 3).
- Produces:
  - `ContigRow` — `NamedTuple(contig_id, length: int, cls: str, host_frac: float, element_frac: float, top_element: str, element_pident: float, element_alnlen: int, remote_hit: str)`.
  - `write_outputs(step_dir: Path, contigs_fa: Path, rows: list[ContigRow], links: list[Link]) -> None` — writes `construct_contigs.fasta` (non-host contigs, whole, annotated header), `contig_classification.tsv` (ALL rows incl. host, header row), `site_construct_links.tsv` (header row), `construct_summary.txt`. Empty `rows` still writes headers + empty FASTA.

- [ ] **Step 1: Write the failing test**

```python
def test_write_outputs_basic(tmp_path):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">c_host\nAAAA\n>c_con\nCCCC\n>c_chim\nGGGG\n")
    rows = [
        s05b.ContigRow("c_host", 4, "host", 0.95, 0.0, "", 0.0, 0, ""),
        s05b.ContigRow("c_con", 4, "construct", 0.0, 0.9, "bar", 99.0, 4, ""),
        s05b.ContigRow("c_chim", 4, "chimeric", 0.4, 0.4, "nptII", 88.0, 4, ""),
    ]
    links = [s05b.Link("Chr3_16439674", "CANDIDATE", "c_chim", "downstream", "high")]
    out = tmp_path / "s05b"
    out.mkdir()
    s05b.write_outputs(out, contigs, rows, links)

    fasta = (out / "construct_contigs.fasta").read_text()
    assert ">c_con" in fasta and ">c_chim" in fasta
    assert ">c_host" not in fasta            # host excluded from FASTA
    assert "class=construct" in fasta

    cls_tsv = (out / "contig_classification.tsv").read_text().splitlines()
    assert cls_tsv[0].startswith("contig_id\t")
    assert any(line.startswith("c_host\t") for line in cls_tsv)  # host kept in TSV

    links_tsv = (out / "site_construct_links.tsv").read_text().splitlines()
    assert links_tsv[0].startswith("site_id\t")
    assert any("c_chim" in line for line in links_tsv)

    summary = (out / "construct_summary.txt").read_text()
    assert "bar" in summary and "nptII" in summary


def test_write_outputs_empty(tmp_path):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text("")
    out = tmp_path / "s05b"
    out.mkdir()
    s05b.write_outputs(out, contigs, [], [])
    assert (out / "construct_contigs.fasta").read_text() == ""
    assert (out / "contig_classification.tsv").read_text().splitlines()[0].startswith("contig_id\t")
    assert (out / "site_construct_links.tsv").read_text().splitlines()[0].startswith("site_id\t")
    assert (out / "construct_summary.txt").exists()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_s05b_construct_assembly.py -k write_outputs -q`
Expected: FAIL — `ContigRow`/`write_outputs` not defined.

- [ ] **Step 3: Write minimal implementation**

```python
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
    seqs = read_fasta(contigs_fa) if _is_nonempty_file(contigs_fa) else {}
    by_id = {r.contig_id: r for r in rows}

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

    # classification TSV: ALL rows (host included)
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_s05b_construct_assembly.py -k write_outputs -q`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add scripts/s05b_construct_assembly.py tests/test_s05b_construct_assembly.py
git commit -m "Step 5b [6/10]: output writers (FASTA, classification/links TSV, summary)"
```

---

## Task 7: Orchestration `run()` + argparse + minimap2/SPAdes wrappers + exit-0 contract

**Files:**
- Modify: `scripts/s05b_construct_assembly.py`
- Test: `tests/test_s05b_construct_assembly.py`

**Interfaces:**
- Consumes: every helper above.
- Produces:
  - `_run_minimap2(contigs_fa: Path, host_ref: Path, threads: int) -> list[Aln]` — runs `minimap2 -x asm10 -t {threads}` host_ref→contigs; returns parsed `Aln`. Missing tool/ref or failure → `[]` + warn.
  - `_obtain_contigs(args, step_dir) -> Path` — returns `contigs_all.fasta` if present/non-empty, else runs SPAdes fallback on s03 reads, else writes empty contigs file. Always returns a path (possibly empty).
  - `_resolve_remote_hit(...)` is folded into `run`; remote nt is opt-in via `args.remote_blast` and wrapped so failure never propagates.
  - `run(args: argparse.Namespace) -> int` — orchestrates; returns 0 in all shortfall paths.
  - `parse_args() -> argparse.Namespace`, `if __name__ == "__main__": sys.exit(run(parse_args()))`.

- [ ] **Step 1: Write the failing test** (exit-0 contract + CLI smoke; no external tools needed)

```python
def test_run_empty_inputs_exit0(tmp_path):
    # No contigs-all, no s03 reads, no s05 dir, no host ref → must still exit 0
    args = s05b.parse_args_from([
        "--contigs-all", str(tmp_path / "absent_contigs.fasta"),
        "--s03-r1", str(tmp_path / "absent_R1.fq.gz"),
        "--s03-r2", str(tmp_path / "absent_R2.fq.gz"),
        "--s05-dir", str(tmp_path / "absent_s05"),
        "--host-ref", str(tmp_path / "absent_host.fa"),
        "--element-db", str(tmp_path / "absent_db.fa"),
        "--outdir", str(tmp_path / "results"),
        "--sample-name", "tsample",
        "--threads", "1",
    ])
    rc = s05b.run(args)
    assert rc == 0
    out = tmp_path / "results" / "tsample" / "s05b_construct_assembly"
    assert (out / "construct_contigs.fasta").read_text() == ""
    assert (out / "contig_classification.tsv").exists()
    assert (out / "site_construct_links.tsv").exists()
    assert (out / "construct_summary.txt").exists()


def test_help_smoke(capsys):
    with pytest.raises(SystemExit) as exc:
        s05b.parse_args_from(["--help"])
    assert exc.value.code == 0
    assert "construct" in capsys.readouterr().out.lower()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_s05b_construct_assembly.py -k "exit0 or help_smoke" -q`
Expected: FAIL — `parse_args_from`/`run` not defined.

- [ ] **Step 3: Write minimal implementation**

Add the minimap2 wrapper, SPAdes fallback, orchestration, and argparse. Use a
`parse_args_from(argv)` helper so tests can pass an explicit argv.

```python
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
    seqs = read_fasta(contigs_fa) if _is_nonempty_file(contigs_fa) else {}
    with out.open("w") as fh:
        for cid, seq in seqs.items():
            if len(seq) >= min_len:
                fh.write(f">{cid}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
    return out


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
    elem_dbs = [db for db in (args.element_db or [])]
    elem_hits = _run_blastn(contigs, elem_dbs, ELEMENT_PIDENT)

    # Top element annotation: re-run element blastn capturing sseqid/pident/length
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

    # optional remote nt on foreign_unknown
    if args.remote_blast:
        _annotate_remote_unknowns(rows, contigs, dbg)

    # per-site linkage (non-host contigs only)
    nonhost_fa = step_dir / "construct_contigs.fasta"
    # write_outputs writes the FASTA; minimap needs it, so write first then link,
    # then rewrite with links. Simpler: build a temp non-host FASTA here.
    tmp_nh = step_dir / "_nonhost_for_link.fasta"
    seqs = read_fasta(contigs)
    with tmp_nh.open("w") as fh:
        for r in rows:
            if r.cls != "host" and seqs.get(r.contig_id):
                fh.write(f">{r.contig_id}\n{seqs[r.contig_id]}\n")
    sites = _parse_s05_sites(args.s05_dir) if args.s05_dir else []
    alns = _run_minimap2(tmp_nh, args.host_ref, args.threads) if (sites and args.host_ref) else []
    links = _link_contigs_to_sites(alns, sites, contig_class)
    tmp_nh.unlink(missing_ok=True)

    write_outputs(step_dir, contigs, rows, links)
    log(f"done: {len(rows)} contigs, {len(links)} site links")
    return 0


def _top_element_hits(contigs_fa: Path, elem_dbs: list[Path]) -> dict[str, tuple[str, float, int]]:
    """contig_id -> (sseqid, pident, alnlen) for the best element hit."""
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


def _annotate_remote_unknowns(rows: list[ContigRow], contigs_fa: Path, dbg: Path) -> None:
    """Run remote nt on foreign_unknown contigs; mutate rows in place. Never raises."""
    unknown_ids = [r.contig_id for r in rows if r.cls == "foreign_unknown"]
    if not unknown_ids:
        return
    seqs = read_fasta(contigs_fa)
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


def parse_args_from(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--contigs-all", type=Path, required=True,
                   help="s04b contigs_all.fasta (pre-filter SPAdes). Reused if present.")
    p.add_argument("--s03-r1", type=Path, required=True, help="SPAdes-fallback R1")
    p.add_argument("--s03-r2", type=Path, required=True, help="SPAdes-fallback R2")
    p.add_argument("--s05-dir", type=Path, default=None, help="s05_insert_assembly dir for linkage")
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_s05b_construct_assembly.py -q`
Expected: PASS (all tasks 1-7 tests).

- [ ] **Step 5: Commit**

```bash
git add scripts/s05b_construct_assembly.py tests/test_s05b_construct_assembly.py
git commit -m "Step 5b [7/10]: run() orchestration, minimap2/SPAdes wrappers, exit-0 contract"
```

---

## Task 8: Wire `5b` into `run_pipeline.py` + dispatch test

**Files:**
- Modify: `run_pipeline.py` (`STEP_ORDER` line 65; `STEP_SCRIPTS` dict lines ~40-47; `build_step_cmd` add branch after the `"5"` branch ends ~line 300; `--steps` help text)
- Create: `tests/test_run_pipeline_step5b.py`

**Interfaces:**
- Consumes: `run_pipeline.build_step_cmd`, `run_pipeline.parse_steps`.
- Produces: `STEP_SCRIPTS["5b"] == "scripts/s05b_construct_assembly.py"`; a `build_step_cmd("5b", …)` argv containing the 5b script and the expected flags.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_run_pipeline_step5b.py
from pathlib import Path

import run_pipeline


def test_step_order_has_5b_after_5():
    order = run_pipeline.STEP_ORDER
    assert "5b" in order
    assert order.index("5b") == order.index("5") + 1


def test_parse_steps_excludes_5b_from_1_5_but_includes_in_1_5b():
    assert "5b" not in run_pipeline.parse_steps("1-5")
    assert "5b" in run_pipeline.parse_steps("1-5b")
    assert run_pipeline.parse_steps("5b") == {"5b"}


def test_build_step_cmd_5b_argv(tmp_path):
    sample_cfg = {
        "host_reference": "db/host.fa",
        "construct_reference": "db/construct.fa",
        "reads": {"r1": "x_R1.fq.gz", "r2": "x_R2.fq.gz"},
    }
    cmd = run_pipeline.build_step_cmd(
        "5b", "tsample", sample_cfg,
        outdir=tmp_path / "results", threads=4, base_dir=Path("/repo"),
        no_remote_blast=True,
    )
    joined = " ".join(cmd)
    assert "scripts/s05b_construct_assembly.py" in joined
    assert "--contigs-all" in cmd and "--s05-dir" in cmd
    assert "--sample-name" in cmd and "tsample" in cmd
    assert "--remote-blast" not in cmd      # no_remote_blast=True suppresses it


def test_build_step_cmd_5b_remote_when_enabled(tmp_path):
    sample_cfg = {
        "host_reference": "db/host.fa",
        "construct_reference": "db/construct.fa",
        "reads": {"r1": "x_R1.fq.gz", "r2": "x_R2.fq.gz"},
    }
    cmd = run_pipeline.build_step_cmd(
        "5b", "tsample", sample_cfg,
        outdir=tmp_path / "results", threads=4, base_dir=Path("/repo"),
        no_remote_blast=False,
    )
    assert "--remote-blast" in cmd
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_run_pipeline_step5b.py -q`
Expected: FAIL — `"5b"` not in `STEP_ORDER`; `build_step_cmd` raises `KeyError` on `STEP_SCRIPTS["5b"]`.

- [ ] **Step 3: Write minimal implementation**

In `run_pipeline.py`:

(a) `STEP_SCRIPTS` dict — add the entry after the `"5"` line:

```python
    "5": "scripts/s05_insert_assembly.py",
    "5b": "scripts/s05b_construct_assembly.py",
    "6": "scripts/s06_indel.py",
```

(b) `STEP_ORDER` (line 65) — insert `"5b"` after `"5"`:

```python
STEP_ORDER: list[str] = ["1", "2", "3", "4", "4b", "5", "5b", "6", "7", "8"]
```

(c) `build_step_cmd` — add a branch immediately after the `"5"` branch's
`return cmd` (around line 300), before `elif step == "6":`:

```python
    elif step == "5b":
        # Construct assembly: reconstruct + characterize the inserted construct.
        # Reuses s04b contigs_all.fasta; links to s05 CANDIDATE sites.
        memory_gb = max(8, threads * 4 // 2)
        cmd = [sys.executable, script,
               "--contigs-all", str(s04b / "contigs_all.fasta"),
               "--s03-r1", str(s03 / f"{sname}_construct_R1.fq.gz"),
               "--s03-r2", str(s03 / f"{sname}_construct_R2.fq.gz"),
               "--s05-dir", str(s05),
               "--host-ref", host_ref,
               "--element-db", rp("element_db/gmo_combined_db.fa"),
               "--element-db", rp("element_db/common_payload.fa"),
               "--outdir", str(outdir),
               "--sample-name", sname,
               "--threads", str(threads),
               "--memory-gb", str(memory_gb)]
        if not no_remote_blast:
            cmd.append("--remote-blast")
        return cmd
```

(d) `--steps` help text — update the argparse help string to mention `5b`
(find the `--steps` argument near line 660 and append "5b" to its listed steps).

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_run_pipeline_step5b.py -q`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add run_pipeline.py tests/test_run_pipeline_step5b.py
git commit -m "Step 5b [8/10]: wire 5b into run_pipeline (STEP_ORDER, dispatch, --steps)"
```

---

## Task 9: Step 8 PDF "Construct Inventory" section

**Files:**
- Modify: `scripts/reports/insertion_pdf.py`
- Create: `tests/test_insertion_pdf_construct.py`

**Interfaces:**
- Consumes: `s05b_construct_assembly/` outputs (`construct_summary.txt`, `contig_classification.tsv`, `site_construct_links.tsv`).
- Produces:
  - `_load_construct_inventory(s05b_dir: Path) -> dict | None` — returns parsed summary + links, or `None` when the dir/files are absent.
  - `_page_construct_inventory(pdf, inv: dict, ...) -> None` — renders one matplotlib page. Called from the existing `generate_pdf` dispatcher only when `_load_construct_inventory` is non-None.

First inspect the current dispatcher to match the established section pattern
(e.g. the Section 8 CRISPR panel added in Issue #6): read
`scripts/reports/insertion_pdf.py` and `tests/test_insertion_pdf_section8.py` to
mirror naming, the `PdfPages` usage, and the present/absent guard convention.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_insertion_pdf_construct.py
from pathlib import Path

import importlib.util

_SPEC = importlib.util.spec_from_file_location(
    "insertion_pdf",
    Path(__file__).resolve().parent.parent / "scripts" / "reports" / "insertion_pdf.py",
)
ipdf = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(ipdf)


def _write_s05b(d: Path):
    d.mkdir(parents=True, exist_ok=True)
    (d / "construct_summary.txt").write_text(
        "# Step 5b construct inventory\n"
        "total_construct_bp\t12000\n"
        "contigs_construct\t3\ncontigs_chimeric\t1\n"
        "contigs_foreign_unknown\t1\ncontigs_host\t10\n"
        "distinct_elements\t2\nelement_inventory\tbar,nptII\n"
        "sites_linked\t1\nnovel_payload_detected\tyes\n"
    )
    (d / "site_construct_links.tsv").write_text(
        "site_id\tverdict\tcontig_id\tjunction_side\tlink_confidence\n"
        "Chr3_16439674\tCANDIDATE\tc_chim\tdownstream\thigh\n"
    )
    (d / "contig_classification.tsv").write_text(
        "contig_id\tlength\tclass\thost_frac\telement_frac\ttop_element\t"
        "element_pident\telement_alnlen\tremote_hit\n"
        "c_chim\t4000\tchimeric\t0.4\t0.4\tnptII\t88.0\t400\t\n"
    )


def test_load_construct_inventory_present(tmp_path):
    _write_s05b(tmp_path / "s05b_construct_assembly")
    inv = ipdf._load_construct_inventory(tmp_path / "s05b_construct_assembly")
    assert inv is not None
    assert inv["summary"]["element_inventory"] == "bar,nptII"
    assert inv["summary"]["novel_payload_detected"] == "yes"
    assert len(inv["links"]) == 1
    assert inv["links"][0]["contig_id"] == "c_chim"


def test_load_construct_inventory_absent(tmp_path):
    assert ipdf._load_construct_inventory(tmp_path / "absent") is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_insertion_pdf_construct.py -q`
Expected: FAIL — `_load_construct_inventory` not defined.

- [ ] **Step 3: Write minimal implementation**

Add to `scripts/reports/insertion_pdf.py` (loader + page renderer). Loader:

```python
def _load_construct_inventory(s05b_dir: Path) -> dict | None:
    """Parse s05b outputs into {summary: dict, links: list[dict]} or None."""
    summary_p = s05b_dir / "construct_summary.txt"
    links_p = s05b_dir / "site_construct_links.tsv"
    if not summary_p.is_file():
        return None
    summary: dict[str, str] = {}
    for line in summary_p.read_text().splitlines():
        if line.startswith("#") or "\t" not in line:
            continue
        k, v = line.split("\t", 1)
        summary[k] = v
    links: list[dict] = []
    if links_p.is_file():
        rows = links_p.read_text().splitlines()
        if rows:
            header = rows[0].split("\t")
            for r in rows[1:]:
                cells = r.split("\t")
                if len(cells) == len(header):
                    links.append(dict(zip(header, cells)))
    return {"summary": summary, "links": links}
```

Then add `_page_construct_inventory(pdf, inv, sample_name)` rendering a matplotlib
table page (mirror the table-drawing helper used by the Section 8 page; reuse the
existing `_text_page`/table helper rather than inventing a new one), and call it
from `generate_pdf` guarded by:

```python
    inv = _load_construct_inventory(out_root / sample_name / "s05b_construct_assembly")
    if inv is not None:
        _page_construct_inventory(pdf, inv, sample_name)
```

(Use the actual `generate_pdf` parameter names found when reading the file in
Step 0; the path root variable may be named `outdir`/`out_root` — match it.)

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_insertion_pdf_construct.py -q`
Expected: PASS (2 tests). Also run the existing PDF tests to confirm no regression:
`pytest tests/test_insertion_pdf_section8.py tests/test_insertion_pdf_kcgp_cover.py -q` → PASS.

- [ ] **Step 5: Commit**

```bash
git add scripts/reports/insertion_pdf.py tests/test_insertion_pdf_construct.py
git commit -m "Step 5b [9/10]: step-8 PDF Construct Inventory section"
```

---

## Task 10: Move `s09_ground_truth_baseline.py` → `util/`, docs, full suite

**Files:**
- Move: `scripts/s09_ground_truth_baseline.py` → `scripts/util/ground_truth_baseline.py`
- Modify: any `tests/` and `docs/measurements/coverage_sensitivity.md` references to the old path
- Modify: `CLAUDE.md` (document step 5b in the step table + flow)

**Interfaces:** none (relocation + docs).

- [ ] **Step 1: Find references to the old path**

Run:
```bash
grep -rn "s09_ground_truth_baseline" --include="*.py" --include="*.md" --include="*.sh" \
  scripts tests docs CLAUDE.md 2>/dev/null | grep -v "docs/superpowers/plans/"
```
Record each hit. (Historical plan docs under `docs/superpowers/plans/` are
immutable — do NOT edit them.)

- [ ] **Step 2: Move the file with git**

```bash
git mv scripts/s09_ground_truth_baseline.py scripts/util/ground_truth_baseline.py
```

- [ ] **Step 3: Update live references**

For each non-plan reference found in Step 1, update the path
`scripts/s09_ground_truth_baseline.py` → `scripts/util/ground_truth_baseline.py`
(and bare `s09_ground_truth_baseline` invocations likewise). If a test imports it
by file path, fix the path string.

- [ ] **Step 4: Document step 5b in CLAUDE.md**

In `CLAUDE.md`, add a row to the Key scripts table:

```markdown
| `scripts/s05b_construct_assembly.py` | Step 5b: reconstruct + characterize inserted construct (global inventory + per-site links). Reuses s04b contigs_all.fasta; assemble-then-classify; opt-in remote nt. |
```

Update the `STEP_ORDER` line in CLAUDE.md to
`["1", "2", "3", "4", "4b", "5", "5b", "6", "7", "8"]` and add a one-line flow note:
`→ [5b] construct assembly (opt-in: global construct inventory + per-site links)`.

- [ ] **Step 5: Run the full suite + smoke the step CLI**

```bash
export PATH="/data/gpfs/assoc/pgl/bin/conda/conda_envs/redgene/bin:$PATH"
pytest tests/ -q
python scripts/s05b_construct_assembly.py --help >/dev/null && echo "CLI OK"
pytest tests/test_s05_import_dag.py tests/test_verdict_snapshots.py -q   # snapshots: 0 drift
```
Expected: full suite PASS (prior baseline + new tests), `CLI OK`, snapshot tests PASS with zero drift.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "Step 5b [10/10]: relocate ground_truth_baseline to util/, document 5b, full suite green"
```

---

## Self-Review

**Spec coverage:**
- §2 pipeline position (STEP_ORDER 5b, opt-in) → Task 8 ✓
- §3 inputs (contigs_all reuse, SPAdes fallback, host/element DB, s05 sites, remote opt-in) → Tasks 3,4,7,8 ✓
- §4 algorithm (obtain/length-filter, classify with thresholds, annotate, linkage) → Tasks 1-7 ✓
- §5 outputs (4 files + headers + whole chimeric) → Task 6 ✓
- §6 exit-0/edge handling → Tasks 6,7 (empty-input tests) ✓
- §7 run_pipeline wiring → Task 8 ✓
- §8 step-8 PDF section → Task 9 ✓
- §9 name-collision cleanup → Task 10 ✓
- §10 module structure (standalone, reuse primitives.read_fasta) → Tasks 1,3 ✓
- §11 testing (classification, helpers, site parse, linkage, edge/exit-0, --help, dispatch, PDF) → Tasks 1-9 ✓
- §12 out-of-scope (code-advancement track) → not in plan ✓

**Placeholder scan:** Task 9 intentionally defers the exact `_page_construct_inventory` body to a read-the-file step because it must mirror an existing matplotlib helper whose name/signature is in `insertion_pdf.py`; the loader (the testable logic) is fully specified. All other steps contain complete code.

**Type consistency:** `ContigRow`, `Site`, `Aln`, `Link` field names are used identically across Tasks 4-7. `_run_blastn` signature (query, subjects:list, min_pident) is consistent between Tasks 3 and 7. `STEP_SCRIPTS`/`STEP_ORDER` names verified against `run_pipeline.py:65,165`.
