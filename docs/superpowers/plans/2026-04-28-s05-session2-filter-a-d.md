# s05 Session 2 — Filter A + Filter D Extraction (Issue #4)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract `_blast_insert_vs_host` (Filter A) into `scripts/s05/filter_a_host.py` and `_check_construct_host_coverage` (Filter D) into `scripts/s05/filter_d_altlocus.py`, mirroring the Session 1 pattern (verbatim move + monolith re-import).

**Architecture:** Each filter module is a side-effect-free module (BLAST subprocess only) that depends only on `scripts/s05/primitives` (stage 0). Filters live at DAG stage 5. The monolith `scripts/s05_insert_assembly.py` imports symbols back into its namespace via the existing `try/except ImportError` block so all callers continue to work unchanged. The `tests/test_s05_import_dag.py` already pre-registers `filter_a_host` and `filter_d_altlocus` at stage 5, so no test edit is needed.

**Tech Stack:** Python 3.11, pytest 9, micromamba env `redgene`. BLAST+ subprocess. No new deps.

---

## Reference: Session 1 Pattern (filter_b_flank.py / filter_c_chimeric.py)

Each filter module:
1. Verbatim function copy (same algorithm, same constants)
2. `from .primitives import log` (only allowed s05 sibling import; primitives is stage 0)
3. Module-level threshold constants (mirrored from monolith)
4. Module docstring explains rationale + cites verdict.py rule
5. `from __future__ import annotations`

Monolith changes:
1. Add filter symbols to existing `try: from s05.X import ...` block (both branches)
2. Replace original function definition with a comment pointer
3. Constants stay defined at top of monolith (mirrored), so callers using bare `INSERT_HOST_FRACTION` etc. don't break

---

## File Map

| File | Action | Purpose |
|------|--------|---------|
| `scripts/s05/filter_a_host.py` | **Create** | `_blast_insert_vs_host` + `INSERT_HOST_MIN_PIDENT` |
| `scripts/s05/filter_d_altlocus.py` | **Create** | `_check_construct_host_coverage` + `CONSTRUCT_HOST_MIN_PIDENT`, `CONSTRUCT_HOST_MIN_COMBINED`, `CONSTRUCT_MIN_FRACTION` |
| `scripts/s05/__init__.py` | **Modify** | Re-export new symbols (keep public package API stable) |
| `scripts/s05_insert_assembly.py` | **Modify** | Two `try/except ImportError` branches gain new imports; original function bodies replaced with pointer comments |

---

## Task 1: Extract Filter A (`_blast_insert_vs_host`)

**Files:**
- Create: `scripts/s05/filter_a_host.py`
- Modify: `scripts/s05/__init__.py`
- Modify: `scripts/s05_insert_assembly.py:48-101` (import block) and `:3130-3209` (original function)

- [ ] **Step 1: Verify pre-state — pytest baseline + DAG test**

Run:
```bash
cd /data/gpfs/assoc/pgl/develop/redgene
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene
pytest tests/ -q
pytest tests/test_s05_import_dag.py -v
```
Expected:
- pytest: 235 PASS, 1 skipped (~20s)
- DAG test: 10 PASSED

Record the count exactly (`grep "passed"`); Step 9 reverify must match.

- [ ] **Step 2: Create `scripts/s05/filter_a_host.py` with verbatim function**

Create file with this exact content:

```python
"""Filter A — host-fraction false-positive check (Issue #4, Session 2).

Extracted verbatim from ``scripts/s05_insert_assembly.py`` so the BLAST-driven
host-fraction measurement can be unit-tested without the full monolith.

Filter A rationale
------------------
A genuine T-DNA insert has a contiguous foreign (non-host) stretch that
contains the construct sequence. False positives caused by host-derived
construct elements (e.g., rice Ubi1 promoter copies) lack this gap and
align almost end-to-end against the host. Filter A measures:

    host_fraction         = host-aligned bp / (insert_len - N_count)
    largest_foreign_gap   = longest contiguous non-host stretch (bp)

The monolith and ``scripts/s05/verdict.py`` Rule 1/6 consume the returned
tuple via ``FilterEvidence``. The N-count subtraction protects against
gap-fill placeholders inflating the denominator.

This module is intentionally side-effect free aside from the BLAST
subprocess call. The on-disk BLAST output filename pattern
(``_<insert_stem>_vs_host_chrom.tsv``) is shared with
``filter_c_chimeric._check_chimeric_assembly`` so the second filter can
reuse the first filter's output without re-running megablast.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from .primitives import read_fasta


# ---------------------------------------------------------------------------
# Thresholds (mirrored from the monolith so behaviour is bit-identical).
# ---------------------------------------------------------------------------
INSERT_HOST_MIN_PIDENT = 90.0   # min identity for host alignment to count


# ---------------------------------------------------------------------------
# BLAST-driven host-fraction measurement
# ---------------------------------------------------------------------------

def _blast_insert_vs_host(
    insert_fasta: Path,
    host_ref: Path,
    workdir: Path,
    threads: int = 4,
) -> tuple[float, int, int, int]:
    """BLAST assembled insert against host genome to measure host-fraction.

    Returns (host_fraction, host_covered_bp, insert_length, largest_foreign_gap).
    host_fraction = fraction of insert positions covered by host alignments.
    largest_foreign_gap = longest contiguous stretch NOT covered by host.
    """
    seqs = read_fasta(insert_fasta)
    if not seqs:
        return 0.0, 0, 0, 0
    insert_seq = list(seqs.values())[0]
    insert_len = len(insert_seq)
    if insert_len == 0:
        return 0.0, 0, 0, 0

    # Use _chrom variant filename so _check_chimeric_assembly can reuse it
    blast_out = workdir / f"_{insert_fasta.stem}_vs_host_chrom.tsv"
    if not blast_out.exists():
        result = subprocess.run(
            ["blastn", "-task", "megablast",
             "-query", str(insert_fasta), "-db", str(host_ref),
             "-outfmt", "6 qseqid qstart qend sseqid pident length",
             "-evalue", "1e-10", "-max_target_seqs", "10",
             "-num_threads", str(threads),
             "-out", str(blast_out)],
            stderr=subprocess.DEVNULL,
        )
    if not blast_out.exists():
        return 0.0, 0, insert_len, insert_len

    # Merge overlapping host-aligned intervals to compute non-redundant coverage
    intervals: list[tuple[int, int]] = []
    with open(blast_out) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue
            pident = float(cols[4])
            if pident < INSERT_HOST_MIN_PIDENT:
                continue
            q_start = int(cols[1])
            q_end = int(cols[2])
            lo, hi = min(q_start, q_end), max(q_start, q_end)
            intervals.append((lo, hi))

    if not intervals:
        return 0.0, 0, insert_len, insert_len

    # Merge overlapping intervals
    intervals.sort()
    merged: list[tuple[int, int]] = [intervals[0]]
    for lo, hi in intervals[1:]:
        if lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))

    host_bp = sum(hi - lo + 1 for lo, hi in merged)

    # Compute largest gap between host-aligned regions (= foreign/T-DNA region)
    # Include gaps at start and end of insert
    gaps: list[int] = []
    gaps.append(merged[0][0] - 1)                    # gap before first host hit
    for i in range(1, len(merged)):
        gaps.append(merged[i][0] - merged[i - 1][1] - 1)
    gaps.append(insert_len - merged[-1][1])           # gap after last host hit
    largest_gap = max(gaps) if gaps else 0

    # Exclude N-runs from denominator (Ns are gap-fill placeholders, not real sequence)
    n_count = insert_seq.upper().count("N")
    effective_len = insert_len - n_count
    if effective_len <= 0:
        return 0.0, host_bp, insert_len, largest_gap

    return host_bp / effective_len, host_bp, insert_len, largest_gap
```

**Important:** the function body must be byte-identical to the original (lines 3130–3209 of the monolith). The only change vs. the original is that `read_fasta` is now imported from `.primitives` instead of being a module-level helper. The function itself never calls `log()` — `subprocess.DEVNULL` swallows stderr and the function silently returns zeros on failure — so `log` is intentionally **not** imported.

- [ ] **Step 3: Add re-exports in `scripts/s05/__init__.py`**

After the `from .filter_c_chimeric import (...)` block (currently lines 30–34), add:

```python
from .filter_a_host import (  # noqa: F401
    INSERT_HOST_MIN_PIDENT,
    _blast_insert_vs_host,
)
```

Update the docstring's session-1 sentence (lines 4–9 of `__init__.py`) to mention Session 2:
```
Session 2 (2026-04-28) added ``filter_a_host`` and ``filter_d_altlocus``
filter modules.
```

- [ ] **Step 4: Wire monolith to import from new module (both ImportError branches)**

In `scripts/s05_insert_assembly.py`, locate the existing two `try/except ImportError` import blocks (currently spanning lines 47–101). Add a new import group **after** the `filter_c_chimeric` block in BOTH branches:

```python
    from s05.filter_a_host import (
        INSERT_HOST_MIN_PIDENT,
        _blast_insert_vs_host,
    )
```

(Use `from s05.filter_a_host import ...` exactly — matches the existing convention; the `s05` prefix already works because of the `sys.path.insert` in the fallback branch.)

Then delete the original module-level constants and function:
- Delete the line `INSERT_HOST_MIN_PIDENT = 90.0   # min identity for host alignment to count` near line 132 of the monolith
- Delete the entire `def _blast_insert_vs_host(...)` definition body (currently lines 3130–3209). Replace with this single comment block (matches Session 1 style):

```python
# _blast_insert_vs_host moved to scripts/s05/filter_a_host.py
# (Issue #4 Session 2) and re-imported at the top of this module for
# backward compatibility with existing callers.

```

- [ ] **Step 5: Run DAG test (must still PASS — `filter_a_host` already in STAGE)**

Run:
```bash
pytest tests/test_s05_import_dag.py -v
```
Expected: **11 PASSED** (was 10; the parametrized `test_module_imports_respect_dag` now picks up the new file).

Specifically expect a new line:
```
tests/test_s05_import_dag.py::test_module_imports_respect_dag[filter_a_host.py] PASSED
```

If this fails because `filter_a_host` imports from `verdict` or any stage > 5, fix the import (only `primitives` is allowed).

- [ ] **Step 6: Run full pytest — must match baseline (no behaviour change)**

Run:
```bash
pytest tests/ -q
```
Expected: same as Step 1 baseline (235 PASS + 1 skipped). The DAG suite gains 1 test, so the new total is 236 PASS + 1 skipped. **Snapshot tests in `tests/fixtures/verdict_snapshots/` MUST still pass — they are the primary regression signal.**

- [ ] **Step 7: Verify monolith CLI smoke**

Run:
```bash
python scripts/s05_insert_assembly.py --help > /dev/null && echo OK
```
Expected: prints `OK` and exits 0. (Confirms imports resolve in standalone mode.)

- [ ] **Step 8: Commit**

```bash
git add scripts/s05/filter_a_host.py scripts/s05/__init__.py scripts/s05_insert_assembly.py
git commit -m "$(cat <<'EOF'
Issue #4 Session 2 [1/2]: extract filter_a_host.py from s05 monolith

_blast_insert_vs_host (Filter A — host-fraction + foreign-gap measurement)
moves verbatim into scripts/s05/filter_a_host.py. INSERT_HOST_MIN_PIDENT
threshold migrates with it. Monolith re-imports both symbols so existing
callers and the assemble_and_evaluate_site() path remain bit-identical.

DAG: filter_a_host already pre-registered at stage 5 in
tests/test_s05_import_dag.py (Session 1 forward-looking entries). Test
suite picks up the new module via the parametrize glob.

Tests: 235 PASS + 1 skipped → 236 PASS + 1 skipped (DAG parametrize +1).
No snapshot drift; no verdict logic touched.

Per docs/superpowers/specs/2026-04-19-s05-module-split-design.md §4 Session 2.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Extract Filter D (`_check_construct_host_coverage`)

**Files:**
- Create: `scripts/s05/filter_d_altlocus.py`
- Modify: `scripts/s05/__init__.py`
- Modify: `scripts/s05_insert_assembly.py` (two ImportError branches + delete original definition + delete original constants)

- [ ] **Step 1: Create `scripts/s05/filter_d_altlocus.py` with verbatim function**

Create file:

```python
"""Filter D — alternative-locus / construct-host-coverage check (Issue #4, Session 2).

Extracted verbatim from ``scripts/s05_insert_assembly.py`` so the BLAST-driven
construct + host combined-coverage analysis can be unit-tested without the
full monolith.

Filter D rationale
------------------
Real T-DNA inserts have low construct coverage (~10%, only border regions
align to the construct because the bulk of the insert is foreign cargo).
False positives created by construct-element homology — e.g., a host
chromosomal copy of CaMV 35S promoter mis-assembled together with reads
from the actual construct — show high construct coverage (~50%, the
"insert" essentially **is** a construct fragment) AND high combined
construct+host coverage (~85%+). Filter D rejects sites that meet both
fractions simultaneously.

The "altlocus" name comes from the design spec: this filter catches
inserts that are really host genomic DNA at an alternative (non-T-DNA)
locus, mis-assembled because the host genome contains a paralog of a
construct element. See ``scripts/s05/verdict.py`` Rule 4.

This module is intentionally side-effect free aside from the BLAST
subprocess call (`-subject` mode, no on-disk DB). The returned
``(is_fp, construct_frac, host_frac, combined_frac)`` tuple is consumed
by the assemble_and_evaluate_site() path and ``compute_verdict``.
"""
from __future__ import annotations

import subprocess
from pathlib import Path


# ---------------------------------------------------------------------------
# Thresholds (mirrored from the monolith so behaviour is bit-identical).
# ---------------------------------------------------------------------------
CONSTRUCT_HOST_MIN_COMBINED = 0.85  # min combined coverage (construct + host)
CONSTRUCT_MIN_FRACTION = 0.25       # min construct coverage to suspect FP
CONSTRUCT_HOST_MIN_PIDENT = 80.0    # identity threshold for construct BLAST


# ---------------------------------------------------------------------------
# BLAST-driven construct+host combined coverage check
# ---------------------------------------------------------------------------

def _check_construct_host_coverage(
    insert_fasta: Path,
    element_db: Path,
    host_fraction: float,
    host_bp: int,
    insert_len: int,
    n_count: int,
    workdir: Path,
    threads: int = 4,
) -> tuple[bool, float, float, float]:
    """Check if insert is fully explained by construct + host coverage.

    Real T-DNA inserts have low construct coverage (~10%, only border regions).
    False positives from construct-element homology have high construct
    coverage (~50%, the insert IS the construct fragment).

    Returns (is_fp, construct_frac, host_frac, combined_frac).
    """
    eff_len = insert_len - n_count
    if eff_len <= 0:
        return False, 0.0, 0.0, 0.0

    # BLAST insert vs construct/element_db
    blast_out = workdir / f"_{insert_fasta.stem}_vs_construct.tsv"
    result = subprocess.run(
        ["blastn", "-task", "megablast",
         "-query", str(insert_fasta), "-subject", str(element_db),
         "-outfmt", "6 qstart qend pident",
         "-evalue", "1e-5"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    if result.returncode != 0:
        return False, 0.0, 0.0, 0.0

    # Merge overlapping intervals at >= threshold identity
    intervals: list[tuple[int, int]] = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        cols = line.split("\t")
        qs, qe, pid = int(cols[0]), int(cols[1]), float(cols[2])
        if pid >= CONSTRUCT_HOST_MIN_PIDENT:
            intervals.append((min(qs, qe), max(qs, qe)))

    if not intervals:
        return False, 0.0, host_fraction, host_fraction

    intervals.sort()
    merged: list[tuple[int, int]] = []
    for s, e in intervals:
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))

    construct_bp = sum(e - s + 1 for s, e in merged)
    construct_frac = construct_bp / eff_len
    combined_frac = min(1.0, (construct_bp + host_bp) / eff_len)

    is_fp = (construct_frac >= CONSTRUCT_MIN_FRACTION
             and combined_frac >= CONSTRUCT_HOST_MIN_COMBINED)

    return is_fp, construct_frac, host_fraction, combined_frac
```

Function body MUST be byte-identical to monolith lines 3065–3127. The only difference is `from .primitives import log` is **not** added — Filter D never calls log() (verified by inspection). `blast_out` variable is computed but never written to disk because the BLAST output goes to stdout via `result.stdout`; this is a latent dead local in the original — preserve it verbatim to keep the diff minimal.

- [ ] **Step 2: Re-export in `scripts/s05/__init__.py`**

After the `from .filter_a_host import (...)` block from Task 1, add:

```python
from .filter_d_altlocus import (  # noqa: F401
    CONSTRUCT_HOST_MIN_COMBINED,
    CONSTRUCT_MIN_FRACTION,
    CONSTRUCT_HOST_MIN_PIDENT,
    _check_construct_host_coverage,
)
```

- [ ] **Step 3: Wire monolith imports (both ImportError branches)**

In `scripts/s05_insert_assembly.py`, after the `filter_a_host` import block from Task 1, add to BOTH branches:

```python
    from s05.filter_d_altlocus import (
        CONSTRUCT_HOST_MIN_COMBINED,
        CONSTRUCT_MIN_FRACTION,
        CONSTRUCT_HOST_MIN_PIDENT,
        _check_construct_host_coverage,
    )
```

Delete the three module-level threshold definitions in the monolith (currently lines 150–152):
```
CONSTRUCT_HOST_MIN_COMBINED = 0.85  # min combined coverage (construct + host)
CONSTRUCT_MIN_FRACTION = 0.25       # min construct coverage to suspect FP
CONSTRUCT_HOST_MIN_PIDENT = 80.0    # identity threshold for construct BLAST
```

Delete the entire `def _check_construct_host_coverage(...)` definition (currently lines 3065–3127). Replace with:

```python
# _check_construct_host_coverage moved to scripts/s05/filter_d_altlocus.py
# (Issue #4 Session 2) and re-imported at the top of this module for
# backward compatibility with existing callers.

```

- [ ] **Step 4: Run DAG test**

Run:
```bash
pytest tests/test_s05_import_dag.py -v
```
Expected: **12 PASSED** (was 11 after Task 1; gains parametrize for `filter_d_altlocus.py`).

- [ ] **Step 5: Run full pytest**

Run:
```bash
pytest tests/ -q
```
Expected: 237 PASS + 1 skipped (235 baseline + 2 DAG parametrize entries from Tasks 1+2). Snapshot tests must remain green.

- [ ] **Step 6: Verify monolith CLI smoke**

Run:
```bash
python scripts/s05_insert_assembly.py --help > /dev/null && echo OK
```
Expected: `OK`.

- [ ] **Step 7: Commit**

```bash
git add scripts/s05/filter_d_altlocus.py scripts/s05/__init__.py scripts/s05_insert_assembly.py
git commit -m "$(cat <<'EOF'
Issue #4 Session 2 [2/2]: extract filter_d_altlocus.py from s05 monolith

_check_construct_host_coverage (Filter D — construct+host combined
coverage check, catches host genomic DNA mis-assembled via construct-
element homology) moves verbatim into scripts/s05/filter_d_altlocus.py.
Three thresholds (CONSTRUCT_HOST_MIN_COMBINED, CONSTRUCT_MIN_FRACTION,
CONSTRUCT_HOST_MIN_PIDENT) migrate with it.

DAG: filter_d_altlocus pre-registered at stage 5; parametrized DAG test
picks it up automatically. Closes Session 2 of the s05 split (filters
A/B/C/D all extracted).

Tests: 236 PASS + 1 skipped → 237 PASS + 1 skipped. No verdict logic
touched; snapshot fixtures unchanged.

Per docs/superpowers/specs/2026-04-19-s05-module-split-design.md §4 Session 2.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Session 2 Exit Gate

- [ ] **Step 1: Confirm Session 2 exit criteria**

Run all guards in one go:
```bash
pytest tests/ -q
pytest tests/test_s05_import_dag.py -v
python scripts/s05_insert_assembly.py --help > /dev/null && echo "CLI OK"
git log --oneline -3
wc -l scripts/s05_insert_assembly.py scripts/s05/*.py | tail -15
```

Expected:
- pytest: 237 PASS + 1 skipped
- DAG: 12 PASSED
- CLI: `CLI OK`
- git log shows the two new commits + recent history
- Monolith line count drops by ~150 lines (Filter A ~80 + Filter D ~63 + constants 3 + comment-lines added ~6 → net ~−138)

If any guard fails, do **not** push. Investigate, then either fix in a follow-up commit or revert with `git reset --hard HEAD~2` and re-run the plan.

- [ ] **Step 2: Update spec progress note (optional, low priority)**

If you want a paper trail in the spec itself, append to `docs/superpowers/specs/2026-04-19-s05-module-split-design.md` under §4 (or in a new "Progress log" section): `Session 2 completed YYYY-MM-DD; commits <sha1> <sha2>.` This is not required for correctness — git log is the authoritative record.

---

## Risks and rollback

| Risk | Mitigation |
|------|------------|
| Snapshot drift (cucumber_line225 / rice_G281) | Verbatim function move; no behaviour change. If it drifts, revert immediately — implies the move wasn't truly verbatim. |
| Constant import cycle (verdict.py uses `INSERT_HOST_FRACTION` etc.) | Verdict.py only references those names in **comments/docstrings**, not Python imports (verified by grep). No new edges added to the DAG. |
| `_check_chimeric_assembly` reusing `_chrom` BLAST output | Filter A still writes `_<stem>_vs_host_chrom.tsv` to the same workdir; Filter C still reads it. Cross-module file convention preserved. |
| Standalone `python scripts/s05_insert_assembly.py` invocation breaks | Step 7 of each task verifies `--help` smoke. The fallback `sys.path.insert(0, ...)` branch handles the no-package-context case. |

Rollback: each task is one commit. `git revert <sha>` or `git reset --hard HEAD~1` returns to baseline.
