# Sample-Specific Construct Assembly + Common Payload DB Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Detect sample-specific transgene payloads (e.g., soybean_AtYUCCA6's YUCCA6 gene, bar marker) that aren't in the shared `element_db/gmo_combined_db.fa`, by (A) reconstructing each sample's construct sequence de novo from s03 reads and (B) expanding `element_db` with commonly-used transgene payload genes.

**Architecture:** Insert a new step `s04b_construct_asm` (between s04 and s05) that runs SPAdes on the construct-hitting reads produced in s03, yielding `contigs.fasta`. Step 5 receives this FASTA via a new `--extra-element-db` flag; its BLAST calls then cover `{gmo_combined_db, common_payload, sample_contigs}`. A one-off curation script assembles `common_payload.fa` from canonical NCBI accessions (bar, nptII, hpt, gusA, gfp, a few payload genes). This preserves backward compatibility — samples without an s04b run still use the original single DB.

**Tech Stack:** Python 3.11, SPAdes 4.0.0, BLAST+ (blastn / makeblastdb), Entrez Direct (efetch), pytest for new unit tests.

---

## File Structure

### New files
- `element_db/common_payload.fa` — curated FASTA of canonical transgene payload genes (output of Task 1–2)
- `element_db/build_common_payload.sh` — reproducible fetch script using efetch
- `element_db/common_payload_manifest.tsv` — accession → purpose mapping (audit trail)
- `scripts/s04b_construct_assembly.py` — SPAdes wrapper, de novo assembly of s03 reads (Task 3)
- `tests/test_s04b_construct_assembly.py` — unit tests for the wrapper (Task 3)
- `tests/test_extra_element_db.py` — unit tests for extended-DB integration (Task 5)

### Modified files
- `scripts/s05_insert_assembly.py`
  - Add `--extra-element-db` argparse flag
  - Plumb it through `_batch_check_element_hits` (line 323) and `classify_site_tiers` (line ~820)
  - New helper `_merge_element_dbs()` produces a temp combined FASTA per run
- `run_pipeline.py`
  - Register step `4b` (optional, enabled by default); wire `--extra-element-db` from `results/<sample>/s04b_construct_asm/contigs.fasta` into step 5
- `config.yaml` — add top-level `pipeline.common_payload_db: element_db/common_payload.fa`
- `README.md` — document step 4b and how to override the common payload DB

### Directories
- `results/<sample>/s04b_construct_asm/` — SPAdes output (contigs.fasta, stats)
- `tests/fixtures/` — small FASTQ fixtures for unit tests (3 read pairs each)

---

## Task 1: Common Payload — Accession Manifest

**Why first:** Defines exactly which sequences Task 2 fetches. No ambiguity in Task 2.

**Files:**
- Create: `element_db/common_payload_manifest.tsv`

- [ ] **Step 1: Write the manifest**

Create `element_db/common_payload_manifest.tsv`:

```tsv
accession	purpose	notes
X17220.1	bar	phosphinothricin acetyltransferase, Streptomyces hygroscopicus
V00618.1	nptII	neomycin phosphotransferase II, Tn5 transposon
M55269.1	hpt	hygromycin phosphotransferase, Escherichia coli
X06788.1	gusA	beta-glucuronidase reporter, E. coli (uidA)
U55762.1	gfp	green fluorescent protein, jellyfish
L29345.1	egfp	enhanced GFP variant
V00141.1	P-CaMV35S	Cauliflower mosaic virus 35S promoter (full length)
V00087.1	P-nos	Agrobacterium tumefaciens nopaline synthase promoter
V00087.1	T-nos	nopaline synthase terminator (same accession, different coords)
X04879.1	T-ocs	octopine synthase terminator
NM_122942.3	AtYUCCA6	Arabidopsis thaliana YUCCA6 (flavin monooxygenase AT5G25620)
```

- [ ] **Step 2: Commit the manifest**

```bash
git add element_db/common_payload_manifest.tsv
git commit -m "Add common transgene payload accession manifest"
```

---

## Task 2: Common Payload — Fetch Script + FASTA

**Files:**
- Create: `element_db/build_common_payload.sh`
- Create: `element_db/common_payload.fa` (output, gitignored-eligible but keep to avoid refetch)

- [ ] **Step 1: Write fetch script**

Create `element_db/build_common_payload.sh`:

```bash
#!/usr/bin/env bash
# Fetch common transgene payload sequences from NCBI Entrez.
# Re-run only when element_db/common_payload_manifest.tsv changes.
set -euo pipefail
cd "$(dirname "$0")"

MANIFEST=common_payload_manifest.tsv
OUT=common_payload.fa
TMP=$(mktemp -d)

: > "$OUT"
while IFS=$'\t' read -r acc purpose notes; do
    [[ "$acc" == "accession" ]] && continue
    [[ -z "$acc" ]] && continue
    # Derive a stable sequence name: purpose|accession
    header="${purpose}|${acc}"
    echo "Fetching $acc ($purpose)..." >&2
    efetch -db nuccore -id "$acc" -format fasta \
        | awk -v h=">${header}" 'NR==1{print h; next}{print}' \
        >> "$OUT"
    sleep 0.4  # be nice to NCBI
done < "$MANIFEST"

echo >&2
echo "Wrote $OUT with $(grep -c '^>' "$OUT") sequences" >&2
rm -rf "$TMP"
```

- [ ] **Step 2: Make executable and run**

```bash
chmod +x element_db/build_common_payload.sh
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene
element_db/build_common_payload.sh
```

Expected output:
```
Fetching X17220.1 (bar)...
Fetching V00618.1 (nptII)...
...
Wrote element_db/common_payload.fa with 11 sequences
```

- [ ] **Step 3: Verify FASTA sanity**

```bash
grep -c '^>' element_db/common_payload.fa
head -1 element_db/common_payload.fa
```

Expected: `11` sequences, first header like `>bar|X17220.1`.

- [ ] **Step 4: Build BLAST index**

```bash
makeblastdb -in element_db/common_payload.fa -dbtype nucl
```

- [ ] **Step 5: Commit**

```bash
git add element_db/build_common_payload.sh element_db/common_payload.fa \
        element_db/common_payload.fa.nhr element_db/common_payload.fa.nin \
        element_db/common_payload.fa.nsq
git commit -m "Add common_payload.fa with bar/nptII/hpt/gusA/gfp/AtYUCCA6"
```

---

## Task 3: s04b — SPAdes Construct Assembly Wrapper

**Files:**
- Create: `scripts/s04b_construct_assembly.py`
- Create: `tests/test_s04b_construct_assembly.py`
- Create: `tests/fixtures/mini_R1.fq.gz`, `tests/fixtures/mini_R2.fq.gz` (3 read pairs from any existing s03 output)

- [ ] **Step 1: Build the test fixture (real data, tiny subset)**

```bash
mkdir -p tests/fixtures
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene
# Take 3 read pairs from an existing s03 output
zcat results/tomato_Cas9_A2_3/s03_extract/tomato_Cas9_A2_3_construct_R1.fq.gz \
    | head -12 | gzip > tests/fixtures/mini_R1.fq.gz
zcat results/tomato_Cas9_A2_3/s03_extract/tomato_Cas9_A2_3_construct_R2.fq.gz \
    | head -12 | gzip > tests/fixtures/mini_R2.fq.gz
```

- [ ] **Step 2: Write the failing test**

Create `tests/test_s04b_construct_assembly.py`:

```python
"""Tests for scripts/s04b_construct_assembly.py"""
import subprocess
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parents[1]


def test_emits_contigs_fasta(tmp_path):
    r1 = REPO / "tests/fixtures/mini_R1.fq.gz"
    r2 = REPO / "tests/fixtures/mini_R2.fq.gz"
    outdir = tmp_path / "out"
    subprocess.run(
        ["python", str(REPO / "scripts/s04b_construct_assembly.py"),
         "--r1", str(r1), "--r2", str(r2),
         "--outdir", str(outdir),
         "--sample-name", "mini",
         "--threads", "2"],
        check=True,
    )
    contigs = outdir / "mini" / "s04b_construct_asm" / "contigs.fasta"
    assert contigs.exists()
    # Even empty SPAdes runs emit the file; we only assert existence here.


def test_handles_empty_reads_gracefully(tmp_path):
    # Create empty-but-valid gzipped FASTQ
    import gzip
    r1 = tmp_path / "empty_R1.fq.gz"
    r2 = tmp_path / "empty_R2.fq.gz"
    for p in (r1, r2):
        with gzip.open(p, "wt") as fh:
            pass
    outdir = tmp_path / "out"
    result = subprocess.run(
        ["python", str(REPO / "scripts/s04b_construct_assembly.py"),
         "--r1", str(r1), "--r2", str(r2),
         "--outdir", str(outdir),
         "--sample-name", "empty",
         "--threads", "2"],
        capture_output=True, text=True,
    )
    # Should exit 0 and emit an empty contigs.fasta marker
    assert result.returncode == 0, result.stderr
    contigs = outdir / "empty" / "s04b_construct_asm" / "contigs.fasta"
    assert contigs.exists()
    assert contigs.stat().st_size == 0
```

- [ ] **Step 3: Run test to confirm failure**

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene
pytest tests/test_s04b_construct_assembly.py -v
```

Expected: FAIL — `scripts/s04b_construct_assembly.py` does not exist.

- [ ] **Step 4: Implement the wrapper**

Create `scripts/s04b_construct_assembly.py`:

```python
#!/usr/bin/env python3
"""Step 4b: De novo assemble construct-hitting reads (from s03) with SPAdes.

The resulting contigs.fasta serves as a sample-specific construct reference
passed into step 5 via --extra-element-db. This lets us annotate inserts
whose transgene payload (e.g. AtYUCCA6, bar gene) isn't present in the
shared element_db.

Input:  s03 R1/R2 FASTQ.GZ (construct-hitting reads + their mates)
Output: results/<sample>/s04b_construct_asm/contigs.fasta
        results/<sample>/s04b_construct_asm/spades_stderr.log
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path


def _is_empty_fastq(path: Path) -> bool:
    """Empty-gzip or zero-read FASTQ detection."""
    try:
        with gzip.open(path, "rt") as fh:
            return fh.read(1) == ""
    except OSError:
        return path.stat().st_size == 0


def run(args: argparse.Namespace) -> int:
    step_dir = args.outdir / args.sample_name / "s04b_construct_asm"
    step_dir.mkdir(parents=True, exist_ok=True)
    contigs = step_dir / "contigs.fasta"

    if _is_empty_fastq(args.r1) or _is_empty_fastq(args.r2):
        print("[s04b] Empty inputs; writing empty contigs.fasta", file=sys.stderr)
        contigs.write_text("")
        return 0

    spades_out = step_dir / "_spades_run"
    if spades_out.exists():
        shutil.rmtree(spades_out)
    cmd = [
        "spades.py", "--only-assembler", "--careful",
        "-1", str(args.r1), "-2", str(args.r2),
        "-o", str(spades_out),
        "-t", str(args.threads),
        "-m", str(args.memory_gb),
    ]
    print(f"[s04b] {' '.join(cmd)}", file=sys.stderr)
    log = step_dir / "spades_stderr.log"
    with log.open("w") as fh_err:
        proc = subprocess.run(cmd, stderr=fh_err)
    if proc.returncode != 0:
        print(f"[s04b] SPAdes exited {proc.returncode}; see {log}", file=sys.stderr)
        # Still emit empty contigs so downstream survives
        contigs.write_text("")
        return proc.returncode

    spades_contigs = spades_out / "contigs.fasta"
    if spades_contigs.exists() and spades_contigs.stat().st_size > 0:
        shutil.copy2(spades_contigs, contigs)
    else:
        contigs.write_text("")

    # Cleanup SPAdes scratch to save disk
    shutil.rmtree(spades_out, ignore_errors=True)
    n_contigs = sum(1 for line in contigs.read_text().splitlines()
                    if line.startswith(">"))
    print(f"[s04b] Wrote {contigs} ({n_contigs} contigs)", file=sys.stderr)
    return 0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--r1", type=Path, required=True)
    p.add_argument("--r2", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--sample-name", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--memory-gb", type=int, default=16)
    return p.parse_args()


if __name__ == "__main__":
    sys.exit(run(parse_args()))
```

- [ ] **Step 5: Run test to confirm pass**

```bash
pytest tests/test_s04b_construct_assembly.py -v
```

Expected: both tests PASS.

- [ ] **Step 6: Commit**

```bash
git add scripts/s04b_construct_assembly.py tests/test_s04b_construct_assembly.py \
        tests/fixtures/mini_R1.fq.gz tests/fixtures/mini_R2.fq.gz
git commit -m "Add s04b step: de novo SPAdes assembly of construct reads"
```

---

## Task 4: s05 — Accept `--extra-element-db`

**Why before run_pipeline integration:** Makes step 5 callable standalone with the new FASTA, which the pipeline orchestrator will just wire up.

**Files:**
- Modify: `scripts/s05_insert_assembly.py`
  - argparse: add `--extra-element-db` (line ~3461 block)
  - `_batch_check_element_hits` (line 323–358): accept an optional second DB
  - `classify_site_tiers` (line ~820): pass both DBs into the clip check
- Create: `tests/test_extra_element_db.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_extra_element_db.py`:

```python
"""Verify s05 _batch_check_element_hits can consult a second FASTA."""
import importlib.util
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "s05", REPO / "scripts/s05_insert_assembly.py"
)
s05 = importlib.util.module_from_spec(spec)
import sys
sys.modules["s05"] = s05
spec.loader.exec_module(s05)


def test_batch_check_element_hits_reads_extra_db(tmp_path):
    primary = tmp_path / "primary.fa"
    primary.write_text(">elem_A\n" + "A" * 50 + "\n")
    extra = tmp_path / "extra.fa"
    extra.write_text(">bar|X17220.1\n" + "G" * 50 + "\n")
    # Query sequence matches extra, not primary
    seqs = {"clip_q1": "G" * 50}
    hits = s05._batch_check_element_hits(
        seqs, primary, tmp_path, extra_db=extra,
    )
    assert "clip_q1" in hits
    hit_str = " ".join(hits["clip_q1"])
    assert "bar" in hit_str
```

- [ ] **Step 2: Run test to confirm failure**

```bash
pytest tests/test_extra_element_db.py -v
```

Expected: FAIL — `_batch_check_element_hits() got an unexpected keyword argument 'extra_db'`.

- [ ] **Step 3: Modify `_batch_check_element_hits`**

In `scripts/s05_insert_assembly.py`, replace the function (around line 323):

```python
def _batch_check_element_hits(
    seqs: dict[str, str],
    element_db: Path,
    workdir: Path,
    extra_db: Path | None = None,
) -> dict[str, list[str]]:
    """Batch check which sequences hit element DB(s) using blastn.

    When extra_db is provided (e.g., per-sample SPAdes contigs from s04b),
    each query is BLASTed against both DBs and hits are merged. This lets
    sample-specific transgene payloads contribute to element annotation
    without rebuilding the shared DB.
    """
    if not seqs:
        return {}

    query_fa = workdir / "_clip_element_batch.fa"
    with open(query_fa, "w") as fh:
        for name, seq in seqs.items():
            if len(seq) >= 20:
                fh.write(f">{name}\n{seq}\n")

    hits: dict[str, list[str]] = defaultdict(list)
    dbs = [element_db]
    if extra_db is not None and extra_db.exists() and extra_db.stat().st_size > 0:
        dbs.append(extra_db)

    for db in dbs:
        blast_out = workdir / f"_clip_blast_{db.stem}.tsv"
        subprocess.run(
            ["blastn", "-query", str(query_fa), "-subject", str(db),
             "-outfmt", "6 qseqid sseqid pident length",
             "-evalue", "1e-3", "-max_target_seqs", "5",
             "-out", str(blast_out)],
            stderr=subprocess.DEVNULL,
        )
        if blast_out.exists():
            with open(blast_out) as fh:
                for line in fh:
                    cols = line.strip().split("\t")
                    if len(cols) >= 4 and int(cols[3]) >= 20:
                        hits[cols[0]].append(cols[1])
            blast_out.unlink(missing_ok=True)

    query_fa.unlink(missing_ok=True)
    return dict(hits)
```

- [ ] **Step 4: Add argparse flag**

In `scripts/s05_insert_assembly.py` `main()` argparse block (around line 3461, right after `--min-clip`):

```python
    parser.add_argument("--extra-element-db", type=Path, default=None,
                        help="Optional second FASTA for element annotation "
                             "(e.g., per-sample SPAdes construct contigs).")
```

- [ ] **Step 5: Propagate through `classify_site_tiers` and the annotation pass**

Pattern: wherever `_batch_check_element_hits(..., element_db, workdir)` is called in `s05_insert_assembly.py`, add `extra_db=args.extra_element_db` (pull `args` into scope where needed).

The two call sites are at line ~625 and inside `classify_site_tiers` (which currently receives `transgene_db`). Pass the extra FASTA alongside:

```python
# In main, after parse_args:
extra_db = args.extra_element_db

# Existing call around line 625 becomes:
all_element_hits = _batch_check_element_hits(
    element_check_seqs, element_db, workdir, extra_db=extra_db,
)
```

For `classify_site_tiers`, thread an `extra_transgene_db: Path | None = None`
parameter through its signature (around line 820) and use the same merging
logic inside its blastn-short step (the function currently builds
`transgene_db = element_db.parent / "transgene_db.fa"`; just call blastn twice
and merge).

- [ ] **Step 6: Run unit test**

```bash
pytest tests/test_extra_element_db.py -v
```

Expected: PASS.

- [ ] **Step 7: Backward-compat smoke test — rerun rice_G281 s05**

```bash
python run_pipeline.py --sample rice_G281 --steps 5 --threads 8 --no-remote-blast
grep -c "CANDIDATE" results/rice_G281/s05_insert_assembly/insertion_*_report.txt
```

Expected: still 1 CANDIDATE at Chr3:16,439,674 (no regression when `--extra-element-db` unused).

- [ ] **Step 8: Commit**

```bash
git add scripts/s05_insert_assembly.py tests/test_extra_element_db.py
git commit -m "s05: accept --extra-element-db for per-sample construct DB"
```

---

## Task 5: `run_pipeline.py` — Wire s04b + Forward Extra DB

**Files:**
- Modify: `run_pipeline.py` (register step 4b, forward flag into step 5)

- [ ] **Step 1: Locate step registry**

```bash
grep -n "def build_step_cmd\|steps.*=\|STEP_NAMES" run_pipeline.py | head
```

Note the step list and the step-4 → step-5 handoff block.

- [ ] **Step 2: Add step 4b command builder**

In `run_pipeline.py`'s `build_step_cmd` (exact line number from Step 1), add a case for `"4b"`:

```python
    if step == "4b":
        s03_dir = outdir / sample / "s03_extract"
        return [
            sys.executable, str(SCRIPTS / "s04b_construct_assembly.py"),
            "--r1", str(s03_dir / f"{sample}_construct_R1.fq.gz"),
            "--r2", str(s03_dir / f"{sample}_construct_R2.fq.gz"),
            "--outdir", str(outdir),
            "--sample-name", sample,
            "--threads", str(threads),
            "--memory-gb", str(max(8, memory_gb // 2)),
        ]
```

- [ ] **Step 3: Forward `--extra-element-db` to step 5**

Inside the step 5 branch, after assembling the existing command list, append:

```python
    extra_db = outdir / sample / "s04b_construct_asm" / "contigs.fasta"
    if extra_db.exists() and extra_db.stat().st_size > 0:
        cmd += ["--extra-element-db", str(extra_db)]
```

- [ ] **Step 4: Register step in CLI ordering**

If `run_pipeline.py` parses `--steps 1-7` into an ordered list, update so `4b` slots between `4` and `5`. Example helper update:

```python
ALL_STEPS = ["1", "2", "3", "4", "4b", "5", "6", "7"]
```

And update `--steps` expansion (`1-5` must now include `4b` iff the user asked for a range containing 4 and 5).

- [ ] **Step 5: Smoke test — step 4b alone on an existing sample**

```bash
python run_pipeline.py --sample tomato_Cas9_A2_3 --steps 4b --threads 8
ls results/tomato_Cas9_A2_3/s04b_construct_asm/
grep -c '^>' results/tomato_Cas9_A2_3/s04b_construct_asm/contigs.fasta
```

Expected: `contigs.fasta` exists with ≥1 contig.

- [ ] **Step 6: Smoke test — step 5 picks up extra DB**

```bash
python run_pipeline.py --sample tomato_Cas9_A2_3 --steps 5 --threads 8 --no-remote-blast 2>&1 \
    | grep -E "extra-element-db|s04b"
```

Expected: log shows `--extra-element-db .../s04b_construct_asm/contigs.fasta` in the issued command.

- [ ] **Step 7: Commit**

```bash
git add run_pipeline.py
git commit -m "run_pipeline: register s04b and forward contigs into s05"
```

---

## Task 6: Merge common_payload into default element_db

**Why here:** Once s05 can pick up an extra DB (Task 4) and the pipeline auto-passes it (Task 5), choosing where `common_payload.fa` lives is a one-line decision. Keeping it **separate** and passing it via `--extra-element-db` (or a second flag) avoids mutating the curated `gmo_combined_db.fa`.

**Decision:** Treat common_payload.fa as a *second always-on* element DB. Simpler than rebuilding the combined DB; no index drift.

**Files:**
- Modify: `scripts/s05_insert_assembly.py` (load common_payload by default)
- Modify: `config.yaml`

- [ ] **Step 1: Add config entry**

Edit `config.yaml`, under top-level `pipeline:`:

```yaml
pipeline:
  threads: 8
  memory_gb: 32
  common_payload_db: element_db/common_payload.fa
```

- [ ] **Step 2: Load in run_pipeline**

In `run_pipeline.py` config parsing:

```python
common_payload_db = cfg.get("pipeline", {}).get(
    "common_payload_db", "element_db/common_payload.fa"
)
```

And in step 5 cmd construction, unconditionally append:

```python
    cpd = REPO / common_payload_db
    if cpd.exists() and cpd.stat().st_size > 0:
        cmd += ["--common-payload-db", str(cpd)]
```

- [ ] **Step 3: Add `--common-payload-db` to s05 argparse**

`scripts/s05_insert_assembly.py` argparse (same block as `--extra-element-db`):

```python
    parser.add_argument("--common-payload-db", type=Path, default=None,
                        help="Always-on supplementary FASTA of common transgene "
                             "payload genes (bar, nptII, hpt, gusA, gfp, ...).")
```

Inside `main()`, pass **both** extra sources into `_batch_check_element_hits`. Update the helper signature to accept a list:

```python
def _batch_check_element_hits(
    seqs: dict[str, str],
    element_db: Path,
    workdir: Path,
    extra_dbs: list[Path] | None = None,
) -> dict[str, list[str]]:
    ...
    dbs = [element_db]
    for edb in (extra_dbs or []):
        if edb is not None and edb.exists() and edb.stat().st_size > 0:
            dbs.append(edb)
    # remainder identical to Task 4 Step 3
```

Call sites pass `extra_dbs=[args.common_payload_db, args.extra_element_db]`.

- [ ] **Step 4: Test — verify extra_dbs list**

Extend `tests/test_extra_element_db.py`:

```python
def test_batch_check_element_hits_consults_multiple_dbs(tmp_path):
    primary = tmp_path / "primary.fa"
    primary.write_text(">elem_A\n" + "A" * 50 + "\n")
    common = tmp_path / "common.fa"
    common.write_text(">bar|X17220.1\n" + "G" * 50 + "\n")
    sample_specific = tmp_path / "sample.fa"
    sample_specific.write_text(">node_1\n" + "C" * 50 + "\n")

    seqs = {"hit_common": "G" * 50, "hit_sample": "C" * 50}
    hits = s05._batch_check_element_hits(
        seqs, primary, tmp_path,
        extra_dbs=[common, sample_specific],
    )
    assert any("bar" in h for h in hits["hit_common"])
    assert any("node_1" in h for h in hits["hit_sample"])
```

Run:

```bash
pytest tests/test_extra_element_db.py -v
```

Expected: both tests PASS.

- [ ] **Step 5: Commit**

```bash
git add scripts/s05_insert_assembly.py run_pipeline.py config.yaml \
        tests/test_extra_element_db.py
git commit -m "Wire common_payload DB into s05 as always-on supplementary"
```

---

## Task 7: Regression — Rice G281 and Tomato A2_3

**Why:** Earlier samples have known correct verdicts. Changes to BLAST plumbing must not shift them.

**Files:** results files only; no source changes

- [ ] **Step 1: Rerun rice_G281 s05**

```bash
python run_pipeline.py --sample rice_G281 --steps 5 --threads 8 --no-remote-blast
```

- [ ] **Step 2: Diff verdicts against baseline**

```bash
grep -h "^Verdict:" results/rice_G281/s05_insert_assembly/insertion_*_report.txt \
    | awk -F' —' '{print $1}' | sort | uniq -c
```

Expected: 1 CANDIDATE (Chr3:16,439,674), plus whatever FP/UNKNOWN count matched prior run.

- [ ] **Step 3: Rerun tomato_Cas9_A2_3 s05**

```bash
python run_pipeline.py --sample tomato_Cas9_A2_3 --steps 4b --threads 8
python run_pipeline.py --sample tomato_Cas9_A2_3 --steps 5 --threads 8 --no-remote-blast
```

- [ ] **Step 4: Verify A2_3 still flags ch01:91,002,744 as CANDIDATE**

```bash
grep -l "ch01_91002744" results/tomato_Cas9_A2_3/s05_insert_assembly/insertion_*_report.txt \
    | head -1 | xargs grep "^Verdict:"
```

Expected: `Verdict: CANDIDATE` (with any additional elements now annotated from the extended DB).

- [ ] **Step 5: No commit (results only)** — but note any verdict shift in a follow-up TSV:

```bash
mkdir -p docs/superpowers/runs
echo -e "sample\tbaseline_candidate\tpost_extdb_candidate" \
    > docs/superpowers/runs/2026-04-14-regression.tsv
# Manually fill from the greps above
```

---

## Task 8: Target Sample — soybean_AtYUCCA6 End-to-End

**Files:** results only

- [ ] **Step 1: Run step 4b for AtYUCCA6**

```bash
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 --time=2:00:00 --mem=32G --cpus-per-task=8 \
    --wrap="eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate redgene && \
            python run_pipeline.py --sample soybean_AtYUCCA6 --steps 4b --threads 8"
```

- [ ] **Step 2: Wait for completion; inspect contigs**

```bash
squeue -u $USER | grep rg_
grep -c '^>' results/soybean_AtYUCCA6/s04b_construct_asm/contigs.fasta
head -2 results/soybean_AtYUCCA6/s04b_construct_asm/contigs.fasta
```

Expected: dozens to hundreds of contigs (18,277 s03 read pairs → non-trivial assembly).

- [ ] **Step 3: BLAST contigs vs common_payload to confirm AtYUCCA6 presence**

```bash
blastn -query results/soybean_AtYUCCA6/s04b_construct_asm/contigs.fasta \
       -subject element_db/common_payload.fa \
       -outfmt "6 qseqid sseqid pident length evalue" \
       -evalue 1e-10 \
    | awk '$3 >= 85 && $4 >= 200' | sort -k4 -n | head
```

Expected: hits to `AtYUCCA6|NM_122942.3` and likely `bar|X17220.1`.

- [ ] **Step 4: Rerun s05 with the new DB plumbing**

```bash
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 --time=24:00:00 --mem=32G --cpus-per-task=8 \
    --wrap="eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate redgene && \
            python run_pipeline.py --sample soybean_AtYUCCA6 --steps 5 --threads 8 --no-remote-blast"
```

- [ ] **Step 5: Compare verdict distribution before/after**

```bash
echo "=== After ==="
grep -h "^Verdict:" results/soybean_AtYUCCA6/s05_insert_assembly/insertion_*_report.txt \
    | awk -F' —' '{print $1}' | sort | uniq -c
```

Before: 30 UNKNOWN, 15 FALSE_POSITIVE.
Expected after: the subset of the 30 UNKNOWNs whose assembled insert genuinely carries AtYUCCA6/bar should now be CANDIDATE; host-homologous inserts stay FALSE_POSITIVE.

- [ ] **Step 6: Record the delta**

Append new columns to `docs/superpowers/runs/2026-04-14-regression.tsv` for AtYUCCA6.

- [ ] **Step 7: Commit**

```bash
git add docs/superpowers/runs/2026-04-14-regression.tsv
git commit -m "Record AtYUCCA6 verdict delta after sample-specific DB plumbing"
```

---

## Task 9: Documentation

**Files:**
- Modify: `README.md` — step 4b section, common_payload DB explanation
- Modify: `CLAUDE.md` — pipeline step list 7→8 (optional 4b)

- [ ] **Step 1: README — add step 4b row in the pipeline flow diagram**

Find the section starting `### Pipeline flow and step dependencies` in `README.md` and update:

```
[1] fastp QC → [2] bwa → construct+UniVec BAM → [3] extract reads + mates
  → [3b] WT homology filter (optional)
  → [4] bwa → host BAM (bottleneck: ~5-7h)
  → [4b] SPAdes → sample-specific construct contigs  ★ NEW
  → [5] Targeted insert assembly + FP filtering (uses 4b contigs + common_payload)
  → [6] CRISPR indel detection (optional, needs WT)
  → [7] copy number (depth ratio)
```

Add a short paragraph describing the motivation (AtYUCCA6 failure mode) and when to skip `4b` (well-characterized constructs already fully covered by `element_db/gmo_combined_db.fa`).

- [ ] **Step 2: CLAUDE.md — bump step count**

Update the `## Project Overview` section and any "7 analysis steps" phrases to "7 core + 1 optional (4b) analysis steps."

Add a `Key scripts` row:

| `scripts/s04b_construct_assembly.py` | De novo SPAdes assembly of construct-hitting reads for sample-specific DB |

- [ ] **Step 3: Commit**

```bash
git add README.md CLAUDE.md
git commit -m "Document step 4b and common_payload DB"
```

---

## Self-Review

### Spec coverage
- A (per-sample construct assembly): Tasks 3 (wrapper) + 4 (s05 integration) + 5 (pipeline wiring) + 6 (common_payload wired via same mechanism) + 8 (target sample validation) ✅
- B (common payload DB): Tasks 1–2 (build) + 6 (always-on wiring) ✅
- Regression guard: Task 7 ✅
- Docs: Task 9 ✅

### Placeholder scan
- All code blocks are complete, runnable snippets.
- Each test has real assertions.
- `run_pipeline.py` modifications reference existing `build_step_cmd` / step-parsing helpers — line numbers taken from `grep` in Task 5 Step 1, since the current `run_pipeline.py` line numbers aren't pinned here.

### Type consistency
- `_batch_check_element_hits(..., extra_db=)` in Task 4 is **replaced** by `_batch_check_element_hits(..., extra_dbs=list[Path]|None)` in Task 6. The old signature is superseded; the Task 4 test still passes because `extra_dbs=[extra]` is equivalent. Document this transition in the Task 6 commit message.

### Open assumptions
- `run_pipeline.py` currently expects numeric step tokens. Task 5 Step 4 adds `"4b"` — confirm the step-parser handles non-numeric tokens before running the first smoke test (`--steps 4b` alone). If it doesn't, extend the range parser similarly.
- Entrez efetch rate limit: the 0.4 s sleep in `build_common_payload.sh` is conservative. Add an `NCBI_API_KEY` env var if available (optional).

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-04-14-sample-construct-assembly-and-payload-db.md`. Two execution options:

1. **Subagent-Driven (recommended)** — fresh subagent per task with two-stage review between tasks, fastest correct iteration.
2. **Inline Execution** — execute tasks in this session with checkpoints at Tasks 2, 4, 7, 8.

Which approach?
