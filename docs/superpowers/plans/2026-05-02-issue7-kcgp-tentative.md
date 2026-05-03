# Issue #7 KCGP Nomenclature (tentative_v0) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a placeholder KCGP site-ID format (`KCGP-{HOST3}-{EVENT}-{YY}-{SITE4}`) and a standalone post-hoc CLI that emits `kcgp_mapping.tsv` per sample, with the tentative ID surfaced on the PDF cover page under a "pending official spec" warning.

**Architecture:** A pure ID library (`scripts/util/kcgp_id.py`) holds `HOST_CODES`, `build_kcgp_id`, and `parse_kcgp_id`. A standalone CLI (`scripts/util/assign_kcgp_ids.py`) reads `insertion_*_report.txt` files, sorts CANDIDATE sites, and writes `kcgp_mapping.tsv`. `scripts/reports/insertion_pdf.py` gains a soft-fail `load_kcgp_mapping()` loader and the cover page prepends a tentative-ID block when the mapping exists.

**Tech Stack:** Python 3.11, pyyaml (already in env), matplotlib, pytest. No new dependencies.

**Spec:** [`docs/superpowers/specs/2026-05-02-issue7-kcgp-tentative-design.md`](../specs/2026-05-02-issue7-kcgp-tentative-design.md) (commit `307b861`).

**File Structure:**
- Create: `scripts/util/kcgp_id.py` (~80 lines) — `HOST_CODES`, `build_kcgp_id`, `parse_kcgp_id`, `_normalize_host_path`
- Create: `scripts/util/assign_kcgp_ids.py` (~120 lines) — CLI that produces `kcgp_mapping.tsv`
- Modify: `scripts/reports/insertion_pdf.py` — add `load_kcgp_mapping`, prepend cover page block
- Create: `tests/test_kcgp_id.py` — 8 tests
- Create: `tests/test_assign_kcgp_ids.py` — 6 tests
- Create: `tests/test_insertion_pdf_kcgp_cover.py` — 5 tests
- Modify: `config.yaml` — add `event:` field to `rice_G281`, `tomato_Cas9_A2_3`, `cucumber_line225`
- 255-test baseline must remain green; total grows to ~274 PASS.

---

## Task 1: `scripts/util/kcgp_id.py` — pure ID library

**Files:**
- Create: `scripts/util/kcgp_id.py`
- Create: `tests/test_kcgp_id.py`

- [ ] **Step 1.1: Write the test file with 8 tests**

Create `tests/test_kcgp_id.py`:

```python
"""Tests for scripts/util/kcgp_id.py — KCGP tentative_v0 ID format."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "util" / "kcgp_id.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("kcgp_id", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["kcgp_id"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def test_build_basic():
    mod = _load_module()
    out = mod.build_kcgp_id("db/Osativa_323_v7.0.fa", "G281", 2026, 1)
    assert out == "KCGP-OSA-G281-26-0001"


def test_build_all_5_hosts():
    mod = _load_module()
    cases = [
        ("db/Osativa_323_v7.0.fa", "OSA"),
        ("db/SLM_r2.0.pmol.fasta", "SLE"),
        ("db/CucSat_B10v3.fa",     "CSA"),
        ("db/Zm_B73_v5.fa",        "ZMA"),
        ("db/Gmax_v4.0.fa",        "GMA"),
    ]
    for host_path, expected_code in cases:
        out = mod.build_kcgp_id(host_path, "X", 2026, 1)
        assert f"KCGP-{expected_code}-X-26-0001" == out


def test_build_unknown_host_raises():
    mod = _load_module()
    with pytest.raises(KeyError) as exc:
        mod.build_kcgp_id("db/unknown_host.fa", "G281", 2026, 1)
    msg = str(exc.value)
    assert "unknown_host" in msg or "OSA" in msg  # message lists key or known codes


def test_round_trip():
    mod = _load_module()
    cases = [
        ("db/Osativa_323_v7.0.fa", "G281", 2026, 1),
        ("db/SLM_r2.0.pmol.fasta", "A2_3", 2026, 42),
        ("db/CucSat_B10v3.fa",     "line225", 2025, 9999),
    ]
    for host, event, year, site_ord in cases:
        kid = mod.build_kcgp_id(host, event, year, site_ord)
        parsed = mod.parse_kcgp_id(kid)
        assert parsed["host_code"] == mod.HOST_CODES[mod._normalize_host_path(host)]
        assert parsed["event"] == event
        assert parsed["year_2digit"] == year % 100
        assert parsed["site_ord"] == site_ord


def test_parse_invalid_format():
    mod = _load_module()
    with pytest.raises(ValueError):
        mod.parse_kcgp_id("KCGP-XX")
    with pytest.raises(ValueError):
        mod.parse_kcgp_id("not-a-kcgp-id")


def test_event_with_underscore():
    mod = _load_module()
    out = mod.build_kcgp_id("db/SLM_r2.0.pmol.fasta", "A2_3", 2026, 1)
    assert out == "KCGP-SLE-A2_3-26-0001"
    parsed = mod.parse_kcgp_id(out)
    assert parsed["event"] == "A2_3"


def test_normalize_host_path():
    mod = _load_module()
    assert mod._normalize_host_path("db/Osativa_323_v7.0.fa") == "Osativa_323_v7.0"
    assert mod._normalize_host_path("db/SLM_r2.0.pmol.fasta") == "SLM_r2.0.pmol"
    # Already-stripped name passes through
    assert mod._normalize_host_path("Osativa_323_v7.0") == "Osativa_323_v7.0"


def test_site_overflow():
    mod = _load_module()
    # 4-digit field: 9999 OK, 10000 overflow
    mod.build_kcgp_id("db/Osativa_323_v7.0.fa", "X", 2026, 9999)
    with pytest.raises(ValueError):
        mod.build_kcgp_id("db/Osativa_323_v7.0.fa", "X", 2026, 10000)
```

- [ ] **Step 1.2: Run tests, confirm 8 fail**

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
pytest tests/test_kcgp_id.py -v
```
Expected: 8 FAIL, all with `FileNotFoundError` or `AttributeError` (script does not exist yet).

- [ ] **Step 1.3: Implement `scripts/util/kcgp_id.py`**

Create the file with this exact content:

```python
#!/usr/bin/env python3
"""KCGP nomenclature (tentative_v0).

Placeholder ID format: KCGP-{HOST3}-{EVENT}-{YY}-{SITE4}.
Will be replaced with the official KCGP spec when delivered. The
mapping TSV schema is stable across spec versions; only the kcgp_id
values and spec_version tag change.
"""
from __future__ import annotations

import re
from pathlib import Path

HOST_CODES: dict[str, str] = {
    "Osativa_323_v7.0":   "OSA",   # rice
    "SLM_r2.0.pmol":      "SLE",   # tomato
    "CucSat_B10v3":       "CSA",   # cucumber
    "Zm_B73_v5":          "ZMA",   # corn
    "Gmax_v4.0":          "GMA",   # soybean
}

_KCGP_RE = re.compile(r"^KCGP-([A-Z]{3})-([A-Za-z0-9_]+)-(\d{2})-(\d{4})$")


def _normalize_host_path(host_ref_path: str) -> str:
    """Strip leading directory and a single trailing .fa/.fasta extension."""
    name = Path(host_ref_path).name
    for ext in (".fasta", ".fa"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def build_kcgp_id(
    host_ref_path: str,
    event: str,
    year: int,
    site_ord: int,
    spec: str = "tentative_v0",
) -> str:
    """Return the KCGP placeholder ID for one site.

    Raises KeyError if the host_ref normalized name is not registered in
    HOST_CODES. Raises ValueError if site_ord exceeds the 4-digit field.
    """
    norm = _normalize_host_path(host_ref_path)
    if norm not in HOST_CODES:
        raise KeyError(
            f"host '{norm}' (from {host_ref_path}) not in HOST_CODES; "
            f"known: {sorted(HOST_CODES)}"
        )
    if site_ord < 0 or site_ord > 9999:
        raise ValueError(f"site_ord {site_ord} out of 4-digit range (0-9999)")
    host3 = HOST_CODES[norm]
    yy = year % 100
    return f"KCGP-{host3}-{event}-{yy:02d}-{site_ord:04d}"


def parse_kcgp_id(kcgp_id: str) -> dict:
    """Parse a KCGP ID into its 4 fields. Raises ValueError on malformed input."""
    m = _KCGP_RE.match(kcgp_id)
    if not m:
        raise ValueError(f"not a KCGP-formatted ID: {kcgp_id!r}")
    host3, event, yy, site = m.groups()
    return {
        "host_code": host3,
        "event": event,
        "year_2digit": int(yy),
        "site_ord": int(site),
    }
```

- [ ] **Step 1.4: Run tests, confirm 8 PASS**

```bash
pytest tests/test_kcgp_id.py -v
pytest tests/ -q
```
Expected: 8 PASS in `test_kcgp_id.py`. Full suite: 263 PASS + 1 skipped (255 baseline + 8 new).

- [ ] **Step 1.5: Commit**

```bash
git add scripts/util/kcgp_id.py tests/test_kcgp_id.py
git commit -m "$(cat <<'EOF'
Issue #7 [1/4]: add scripts/util/kcgp_id.py (tentative_v0)

Pure ID library: HOST_CODES dict (5 hosts), build_kcgp_id, parse_kcgp_id,
_normalize_host_path. 8 unit tests cover the round-trip, all 5 hosts,
unknown-host KeyError, malformed-input ValueError, underscore-in-event,
and 4-digit overflow.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: `scripts/util/assign_kcgp_ids.py` — CLI that emits `kcgp_mapping.tsv`

**Files:**
- Create: `scripts/util/assign_kcgp_ids.py`
- Create: `tests/test_assign_kcgp_ids.py`

- [ ] **Step 2.1: Write the test file with 6 tests**

Create `tests/test_assign_kcgp_ids.py`:

```python
"""Tests for scripts/util/assign_kcgp_ids.py — KCGP mapping CLI."""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "util" / "assign_kcgp_ids.py"


def _write_config(path: Path, sample: str, host: str, event: str | None) -> None:
    """Minimal config.yaml with one sample entry."""
    body = f"""samples:
  {sample}:
    host_reference: {host}
"""
    if event is not None:
        body += f"    event: {event}\n"
    path.write_text(body)


def _write_report(s05_dir: Path, host_chr: str, pos: int, verdict: str) -> None:
    s05_dir.mkdir(parents=True, exist_ok=True)
    path = s05_dir / f"insertion_{host_chr}_{pos}_report.txt"
    path.write_text(
        "Header line\n"
        "More header\n"
        "\n"
        "Soft-clip junction stats\n"
        "  ...\n"
        "  ...\n"
        "  ...\n"
        f"Verdict: {verdict} — explanation text\n"
    )


def _run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True, text=True,
    )


def _read_mapping(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_no_s05_dir(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/Osativa_323_v7.0.fa", "EVT")
    sample_dir = tmp_path / "results" / "sampX"
    sample_dir.mkdir(parents=True)

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode == 0, res.stderr
    rows = _read_mapping(sample_dir / "kcgp_mapping.tsv")
    assert rows == []  # header-only TSV produces empty DictReader


def test_no_reports(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/Osativa_323_v7.0.fa", "EVT")
    sample_dir = tmp_path / "results" / "sampX"
    (sample_dir / "s05_insert_assembly").mkdir(parents=True)

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode == 0, res.stderr
    assert _read_mapping(sample_dir / "kcgp_mapping.tsv") == []


def test_assigns_ordinal_to_candidate_only(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/Osativa_323_v7.0.fa", "EVT")
    sample_dir = tmp_path / "results" / "sampX"
    s05 = sample_dir / "s05_insert_assembly"
    _write_report(s05, "Chr1", 100, "CANDIDATE")
    _write_report(s05, "Chr2", 200, "FALSE_POSITIVE")
    _write_report(s05, "Chr3", 300, "UNKNOWN")

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode == 0, res.stderr
    rows = _read_mapping(sample_dir / "kcgp_mapping.tsv")
    by_chr = {r["host_chr"]: r for r in rows}
    assert by_chr["Chr1"]["kcgp_id"] == "KCGP-OSA-EVT-26-0001"
    assert by_chr["Chr1"]["verdict"] == "CANDIDATE"
    assert by_chr["Chr2"]["kcgp_id"] == "-"
    assert by_chr["Chr2"]["verdict"] == "FALSE_POSITIVE"
    assert by_chr["Chr3"]["kcgp_id"] == "-"
    assert by_chr["Chr3"]["verdict"] == "UNKNOWN"
    for r in rows:
        assert r["spec_version"] == "tentative_v0"


def test_sort_order_natural(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/Osativa_323_v7.0.fa", "EVT")
    sample_dir = tmp_path / "results" / "sampX"
    s05 = sample_dir / "s05_insert_assembly"
    # Three CANDIDATE rows; insertion-order is Chr2, Chr10, Chr1
    _write_report(s05, "Chr2", 100, "CANDIDATE")
    _write_report(s05, "Chr10", 200, "CANDIDATE")
    _write_report(s05, "Chr1", 50, "CANDIDATE")

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode == 0, res.stderr
    rows = _read_mapping(sample_dir / "kcgp_mapping.tsv")
    candidate_rows = [r for r in rows if r["verdict"] == "CANDIDATE"]
    # Natural-sort: Chr1 < Chr2 < Chr10
    assert [r["host_chr"] for r in candidate_rows] == ["Chr1", "Chr2", "Chr10"]
    assert [r["kcgp_id"] for r in candidate_rows] == [
        "KCGP-OSA-EVT-26-0001",
        "KCGP-OSA-EVT-26-0002",
        "KCGP-OSA-EVT-26-0003",
    ]


def test_missing_event_in_config(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/Osativa_323_v7.0.fa", event=None)
    sample_dir = tmp_path / "results" / "sampX"
    sample_dir.mkdir(parents=True)

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode != 0
    assert "event" in res.stderr.lower()


def test_unknown_host_in_config(tmp_path):
    cfg = tmp_path / "config.yaml"
    _write_config(cfg, "sampX", "db/unknown_genome.fa", "EVT")
    sample_dir = tmp_path / "results" / "sampX"
    sample_dir.mkdir(parents=True)

    res = _run("--sample-dir", str(sample_dir), "--config", str(cfg), "--year", "2026")
    assert res.returncode != 0
    err = res.stderr.lower()
    assert "host" in err or "osa" in err
```

- [ ] **Step 2.2: Run tests, confirm 6 fail**

```bash
pytest tests/test_assign_kcgp_ids.py -v
```
Expected: 6 FAIL with `FileNotFoundError` (script does not exist yet).

- [ ] **Step 2.3: Implement `scripts/util/assign_kcgp_ids.py`**

Create the file with this exact content:

```python
#!/usr/bin/env python3
"""Assign KCGP tentative_v0 IDs to sites in a sample's s05 output.

Reads results/{sample}/s05_insert_assembly/insertion_*_report.txt, sorts
CANDIDATE sites first then by host_chr (natural) and pos_5p, and writes
results/{sample}/kcgp_mapping.tsv with one row per discovered report.
Only CANDIDATE rows receive a 4-digit ordinal; FALSE_POSITIVE / UNKNOWN
rows carry kcgp_id = "-".
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from datetime import datetime
from pathlib import Path

import yaml

# Package import: scripts/util/kcgp_id.py
sys.path.insert(0, str(Path(__file__).resolve().parent))
from kcgp_id import HOST_CODES, build_kcgp_id  # noqa: E402

REPORT_RE = re.compile(r"^insertion_(.+)_(\d+)_report\.txt$")
VERDICT_PRIORITY = {
    "CANDIDATE": 0,
    "FALSE_POSITIVE": 1,
    "UNKNOWN": 2,
}
_NATURAL_RE = re.compile(r"(\d+|\D+)")


def _natural_key(s: str):
    """Split into alphanumeric runs, convert digit runs to int for natural sort."""
    parts = []
    for tok in _NATURAL_RE.findall(s):
        parts.append((0, int(tok)) if tok.isdigit() else (1, tok))
    return parts


def log(msg: str) -> None:
    print(f"[assign_kcgp_ids] {msg}", file=sys.stderr, flush=True)


def parse_report(path: Path) -> tuple[str, int, str] | None:
    """Return (host_chr, pos_5p, verdict) or None if filename is unparseable."""
    m = REPORT_RE.match(path.name)
    if not m:
        log(f"skip {path.name}: filename does not match insertion_*_report.txt")
        return None
    host_chr, pos_str = m.group(1), m.group(2)
    pos = int(pos_str)
    verdict = "UNKNOWN"
    try:
        for line in path.read_text().splitlines():
            if line.startswith("Verdict:"):
                token = line[len("Verdict:") :].strip().split()[0]
                # Strip trailing punctuation (e.g., "CANDIDATE — ...") leaves CANDIDATE
                verdict = token.rstrip(".,—-")
                break
        else:
            log(f"warn {path.name}: no Verdict line, treating as UNKNOWN")
    except OSError as exc:
        log(f"warn {path.name}: read failed ({exc}), treating as UNKNOWN")
    return host_chr, pos, verdict


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--sample-dir", required=True, type=Path,
                   help="results/{sample}/")
    p.add_argument("--config", required=True, type=Path,
                   help="config.yaml")
    p.add_argument("--year", type=int, default=datetime.now().year,
                   help="default: current year")
    args = p.parse_args(argv)

    sample = args.sample_dir.name
    try:
        cfg = yaml.safe_load(args.config.read_text()) or {}
    except (OSError, yaml.YAMLError) as exc:
        log(f"FATAL: cannot read config: {exc}")
        return 1

    samples = cfg.get("samples", {}) or {}
    if sample not in samples:
        log(f"FATAL: sample {sample!r} not in {args.config} (samples: ...)")
        return 1

    entry = samples[sample] or {}
    event = entry.get("event")
    if not event:
        log(
            f"FATAL: sample {sample!r} has no 'event:' field in {args.config}; "
            "add 'event:' under samples.{sample} (e.g., event: G281)"
        )
        return 1

    host_ref = entry.get("host_reference")
    if not host_ref:
        log(f"FATAL: sample {sample!r} has no 'host_reference:' field")
        return 1

    # Validate host_ref is registered before doing any work
    from kcgp_id import _normalize_host_path
    norm = _normalize_host_path(host_ref)
    if norm not in HOST_CODES:
        log(
            f"FATAL: host '{norm}' (from {host_ref}) not in HOST_CODES; "
            f"known: {sorted(HOST_CODES)}"
        )
        return 1

    # Discover reports
    s05_dir = args.sample_dir / "s05_insert_assembly"
    reports: list[tuple[str, int, str, Path]] = []
    if s05_dir.exists():
        for report_path in sorted(s05_dir.glob("insertion_*_report.txt")):
            parsed = parse_report(report_path)
            if parsed is None:
                continue
            host_chr, pos, verdict = parsed
            reports.append((host_chr, pos, verdict, report_path))
    else:
        log(f"warn: {s05_dir} does not exist; emitting header-only mapping")

    # Sort
    reports.sort(key=lambda r: (
        VERDICT_PRIORITY.get(r[2], 99),
        _natural_key(r[0]),
        r[1],
    ))

    # Count CANDIDATEs to validate ordinal range before assigning
    n_cand = sum(1 for r in reports if r[2] == "CANDIDATE")
    if n_cand > 9999:
        log(f"FATAL: {n_cand} CANDIDATE sites exceeds 4-digit ordinal range")
        return 1

    # Assign
    out_path = args.sample_dir / "kcgp_mapping.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ordinal = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "kcgp_id", "internal_site", "verdict", "host_chr", "pos_5p", "spec_version",
        ])
        for host_chr, pos, verdict, report_path in reports:
            if verdict == "CANDIDATE":
                ordinal += 1
                kid = build_kcgp_id(host_ref, event, args.year, ordinal)
            else:
                kid = "-"
            internal = f"{sample}_{host_chr}_{pos}"
            writer.writerow([kid, internal, verdict, host_chr, pos, "tentative_v0"])

    log(f"wrote {out_path} ({len(reports)} rows, {n_cand} CANDIDATE)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2.4: Run tests, confirm 6 PASS**

```bash
pytest tests/test_assign_kcgp_ids.py -v
pytest tests/ -q
```
Expected: 6 PASS in `test_assign_kcgp_ids.py`. Full suite: 269 PASS + 1 skipped (263 from Task 1 + 6 new).

- [ ] **Step 2.5: Commit**

```bash
git add scripts/util/assign_kcgp_ids.py tests/test_assign_kcgp_ids.py
git commit -m "$(cat <<'EOF'
Issue #7 [2/4]: add scripts/util/assign_kcgp_ids.py CLI

Standalone post-hoc CLI: reads insertion_*_report.txt, sorts CANDIDATE
sites first then natural by chr/pos, emits kcgp_mapping.tsv with one
row per discovered report. Only CANDIDATE gets a 4-digit ordinal;
FALSE_POSITIVE / UNKNOWN carry kcgp_id = "-". Hard fails on missing
event field, missing host_reference, or unknown host code.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: PDF cover-page integration

**Files:**
- Modify: `scripts/reports/insertion_pdf.py` — add `load_kcgp_mapping`, edit `_page_cover`
- Create: `tests/test_insertion_pdf_kcgp_cover.py`

- [ ] **Step 3.1: Write the test file with 5 tests**

Create `tests/test_insertion_pdf_kcgp_cover.py`:

```python
"""Tests for KCGP cover-page integration in scripts/reports/insertion_pdf.py."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "reports" / "insertion_pdf.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("insertion_pdf", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["insertion_pdf"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


def _write_mapping(path: Path, rows: list[tuple[str, str, str]]) -> None:
    """rows: list of (kcgp_id, verdict, internal_site). Other columns stubbed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["kcgp_id\tinternal_site\tverdict\thost_chr\tpos_5p\tspec_version"]
    for kid, verdict, internal in rows:
        lines.append(f"{kid}\t{internal}\t{verdict}\tChr1\t100\ttentative_v0")
    path.write_text("\n".join(lines) + "\n")


def test_load_kcgp_mapping_absent(tmp_path):
    mod = _load_module()
    assert mod.load_kcgp_mapping(tmp_path) is None


def test_load_kcgp_mapping_with_candidate(tmp_path):
    mod = _load_module()
    _write_mapping(tmp_path / "kcgp_mapping.tsv", [
        ("KCGP-OSA-G281-26-0001", "CANDIDATE", "rice_G281_Chr3_16439674"),
        ("-", "FALSE_POSITIVE", "rice_G281_Chr1_8923"),
    ])
    out = mod.load_kcgp_mapping(tmp_path)
    assert out is not None
    assert out["kcgp_id"] == "KCGP-OSA-G281-26-0001"
    assert out["spec_version"] == "tentative_v0"


def test_load_kcgp_mapping_corrupt(tmp_path):
    mod = _load_module()
    p = tmp_path / "kcgp_mapping.tsv"
    p.write_bytes(b"\x00\x01garbage")
    # Should not raise; soft-fail returns None
    assert mod.load_kcgp_mapping(tmp_path) is None or isinstance(
        mod.load_kcgp_mapping(tmp_path), dict
    )


def test_cover_with_mapping(tmp_path, monkeypatch):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    _write_mapping(tmp_path / "kcgp_mapping.tsv", [
        ("KCGP-OSA-G281-26-0001", "CANDIDATE", "rice_G281_Chr3_16439674"),
    ])

    captured = {}
    def fake_page_text(pdf, title, lines, **kwargs):
        captured["title"] = title
        captured["lines"] = lines

    monkeypatch.setattr(mod, "_page_text", fake_page_text)

    out_pdf = tmp_path / "cover.pdf"
    with PdfPages(out_pdf) as pdf:
        mod._page_cover(pdf, "rice_G281", {"sample": "rice_G281"}, sample_dir=tmp_path)

    body = "\n".join(captured["lines"])
    assert "KCGP ID" in body
    assert "tentative_v0" in body
    assert "KCGP-OSA-G281-26-0001" in body


def test_cover_skips_dash_only_mapping(tmp_path, monkeypatch):
    mod = _load_module()
    from matplotlib.backends.backend_pdf import PdfPages

    # All rows have kcgp_id = "-" (no CANDIDATE site)
    _write_mapping(tmp_path / "kcgp_mapping.tsv", [
        ("-", "FALSE_POSITIVE", "rice_G281_Chr1_8923"),
        ("-", "UNKNOWN", "rice_G281_Chr1_99"),
    ])

    captured = {}
    def fake_page_text(pdf, title, lines, **kwargs):
        captured["lines"] = lines

    monkeypatch.setattr(mod, "_page_text", fake_page_text)

    out_pdf = tmp_path / "cover.pdf"
    with PdfPages(out_pdf) as pdf:
        mod._page_cover(pdf, "rice_G281", {"sample": "rice_G281"}, sample_dir=tmp_path)

    body = "\n".join(captured["lines"])
    assert "KCGP ID" not in body  # no CANDIDATE = no KCGP block
```

- [ ] **Step 3.2: Run tests, confirm 5 fail**

```bash
pytest tests/test_insertion_pdf_kcgp_cover.py -v
```
Expected: 5 FAIL — `load_kcgp_mapping` does not exist, and `_page_cover` does not accept `sample_dir`.

- [ ] **Step 3.3: Add `load_kcgp_mapping` to `scripts/reports/insertion_pdf.py`**

Insert the new loader immediately AFTER the existing `load_editing_data` function (added in Issue #6 Task 1) and BEFORE `build_sample_info`:

```python
def load_kcgp_mapping(sample_dir: Path) -> dict[str, Any] | None:
    """Return {'kcgp_id', 'spec_version'} for the lowest-ordinal CANDIDATE row, or None.

    Soft-fails on missing file, IO error, or malformed CSV — the PDF still renders.
    """
    path = Path(sample_dir) / "kcgp_mapping.tsv"
    if not path.exists():
        return None
    try:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("verdict") == "CANDIDATE" and row.get("kcgp_id", "-") != "-":
                    return {
                        "kcgp_id": row["kcgp_id"],
                        "spec_version": row.get("spec_version", "tentative_v0"),
                    }
    except (OSError, csv.Error):
        return None
    return None
```

- [ ] **Step 3.4: Modify `_page_cover` to accept `sample_dir` and prepend KCGP block**

Locate `_page_cover` in `scripts/reports/insertion_pdf.py` (around line 289). Current signature:

```python
def _page_cover(pdf: PdfPages, sample_name: str, audit: dict[str, Any]) -> None:
```

Replace the entire function body with:

```python
def _page_cover(
    pdf: PdfPages,
    sample_name: str,
    audit: dict[str, Any],
    sample_dir: Path | None = None,
) -> None:
    lines = [
        f"Sample:        {sample_name}",
        f"Generated:     {audit.get('generated_at', '—')}",
        f"Commit:        {audit.get('pipeline_commit', '—')}  (dirty={audit.get('pipeline_dirty', '—')})",
        "",
        "Input SHA-256:",
    ]
    for k, v in (audit.get("input_sha256") or {}).items():
        lines.append(f"  {k}: {v}")
    lines += ["", "DB manifest:"]
    for db in audit.get("db_manifest", []) or []:
        lines.append(
            f"  - {db.get('name')} md5={db.get('md5')} built={db.get('build_date')} "
            f"seqs={db.get('seq_count')}"
        )
    if not audit:
        lines.append("  (audit_header.json not present — placeholder page)")

    if sample_dir is not None:
        kcgp = load_kcgp_mapping(sample_dir)
        if kcgp:
            lines = [
                f"⚠ KCGP ID ({kcgp['spec_version']} — pending official spec):",
                f"    {kcgp['kcgp_id']}",
                "",
            ] + lines

    _page_text(pdf, "RedGene Insertion Report", lines)
```

- [ ] **Step 3.5: Update the call site inside `generate_pdf`**

Locate the line in `generate_pdf` that calls `_page_cover`:

```python
        _page_cover(pdf, sample_name, audit); pages += 1                # (1)
```

Replace it with:

```python
        _page_cover(pdf, sample_name, audit, sample_dir=sample_dir); pages += 1  # (1)
```

(`sample_dir` is already a local variable in `generate_pdf`.)

- [ ] **Step 3.6: Run new tests + full suite**

```bash
pytest tests/test_insertion_pdf_kcgp_cover.py -v
pytest tests/ -q
```
Expected: 5 PASS in `test_insertion_pdf_kcgp_cover.py`. Full suite: 274 PASS + 1 skipped (269 from Task 2 + 5 new). Zero regressions in `tests/test_insertion_pdf.py` (the existing scaffold tests still call `_page_cover(pdf, sample_name, audit)` without `sample_dir`, which works because the new arg has `None` default).

- [ ] **Step 3.7: Commit**

```bash
git add scripts/reports/insertion_pdf.py tests/test_insertion_pdf_kcgp_cover.py
git commit -m "$(cat <<'EOF'
Issue #7 [3/4]: surface KCGP tentative ID on PDF cover

New loader load_kcgp_mapping() finds the lowest-ordinal CANDIDATE row
in kcgp_mapping.tsv (soft-fails on missing or malformed file). The
cover page now accepts an optional sample_dir; when a KCGP mapping
exists with at least one CANDIDATE, the cover prepends a "⚠ KCGP ID
(tentative_v0 — pending official spec)" block. Mapping absence
preserves the existing cover layout exactly, protecting the 18-test
scaffold suite.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: 3-sample acceptance validation

**Files:**
- Modify: `config.yaml` — add `event:` field to 3 samples
- No code changes. Operational validation.

- [ ] **Step 4.1: Add `event:` field to `config.yaml`**

Edit `config.yaml`. Inside the `samples:` block, find each of the three samples and add an `event:` line indented to match other fields. Final state for each:

```yaml
  rice_G281:
    name: "Rice G281 transgenic event"
    host_reference: db/Osativa_323_v7.0.fa
    construct_reference: db/G281_construct.fa
    event: G281
    reads:
      ...
```

```yaml
  tomato_Cas9_A2_3:
    name: "Micro-Tom Cas9 1-A2-3 (Cas9 present, heterozygous edit)"
    host_reference: db/SLM_r2.0.pmol.fasta
    construct_reference: element_db/gmo_combined_db.fa
    wt_control: tomato_WT
    grna: db/tomato_grna.txt
    event: A2_3
    reads:
      ...
```

```yaml
  cucumber_line225:
    name: "Cucumber thaumatin line 225 (2 T-DNA copies + backbone, Chr2)"
    host_reference: db/CucSat_B10v3.fa
    ...
    event: line225
    reads:
      ...
```

Verify YAML is still valid:

```bash
python -c "import yaml; yaml.safe_load(open('config.yaml'))"
```

Expected: no output (no exception).

- [ ] **Step 4.2: Generate `kcgp_mapping.tsv` for all 3 samples**

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

for s in rice_G281 tomato_Cas9_A2_3 cucumber_line225; do
  python scripts/util/assign_kcgp_ids.py \
    --sample-dir results/$s \
    --config config.yaml
done
```

Verify each mapping TSV exists and has at least one `KCGP-` row:

```bash
for s in rice_G281 tomato_Cas9_A2_3 cucumber_line225; do
  echo "=== $s ==="
  head -3 results/$s/kcgp_mapping.tsv
  echo "CANDIDATE count: $(grep -c '^KCGP-' results/$s/kcgp_mapping.tsv)"
done
```

Expected:
- `rice_G281`: at least one `KCGP-OSA-G281-26-0001` row
- `tomato_Cas9_A2_3`: at least one `KCGP-SLE-A2_3-26-0001` row
- `cucumber_line225`: at least one `KCGP-CSA-line225-26-0001` row

If any sample produces zero `KCGP-` rows, STOP and inspect — possibly `s05_insert_assembly/` has no CANDIDATE-verdict reports for that sample. Report and ask before proceeding.

- [ ] **Step 4.3: Regenerate PDFs with cover-page KCGP block**

```bash
for s in rice_G281 tomato_Cas9_A2_3 cucumber_line225; do
  python run_pipeline.py --sample $s --steps 8
done
```

Verify each PDF cover now shows the tentative ID:

```bash
for s in rice_G281 tomato_Cas9_A2_3 cucumber_line225; do
  echo "=== $s ==="
  pdftotext results/$s/${s}_insertion_report.pdf - 2>/dev/null | head -20 | grep -A 2 "KCGP\|tentative"
done
```

Expected: each output contains the line `⚠ KCGP ID (tentative_v0 — pending official spec):` followed by the actual ID.

If `pdftotext` is unavailable, capture page count instead and confirm it is unchanged from the Issue #6 baseline (12 / 11 / 10 pages respectively) — the cover page got more lines, but page count stays the same.

- [ ] **Step 4.4: Comment on Issue #7 with results**

```bash
~/micromamba/bin/gh issue comment 7 --body "$(cat <<'EOF'
tentative_v0 implementation landed across 4 commits:

- [1/4] scripts/util/kcgp_id.py — pure ID library (HOST_CODES + build/parse, 8 tests)
- [2/4] scripts/util/assign_kcgp_ids.py — CLI emitting kcgp_mapping.tsv (6 tests)
- [3/4] scripts/reports/insertion_pdf.py — cover-page KCGP block (5 tests, soft-fail when mapping absent)
- [4/4] config.yaml — event: field added for rice_G281 / tomato_Cas9_A2_3 / cucumber_line225

Acceptance:
- rice_G281: KCGP-OSA-G281-26-0001
- tomato_Cas9_A2_3: KCGP-SLE-A2_3-26-0001
- cucumber_line225: KCGP-CSA-line225-26-0001

Test suite: 255 → 274 PASS (+19 across 3 new test files). Zero regressions.

⚠ Still BLOCKED on official KCGP spec from team lead. The mapping TSV
schema is stable; once the spec arrives, swap HOST_CODES + build_kcgp_id
body, bump spec_version to v1, and re-run assign_kcgp_ids.py against
each sample. No schema migration required.

Issue stays OPEN pending v1 spec.
EOF
)"
```

- [ ] **Step 4.5: Commit `config.yaml` changes**

```bash
git add config.yaml
git commit -m "$(cat <<'EOF'
Issue #7 [4/4]: add event: field to 3 acceptance samples

rice_G281 → G281, tomato_Cas9_A2_3 → A2_3, cucumber_line225 → line225.
Required by scripts/util/assign_kcgp_ids.py to derive the EVENT field
of the KCGP tentative_v0 ID. Other samples can be added incrementally
as they are processed.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Spec Coverage Matrix

| Spec section | Implemented in |
|--------------|---------------|
| `HOST_CODES` (5 hosts) | Task 1.3 |
| `build_kcgp_id`, `parse_kcgp_id`, `_normalize_host_path` | Task 1.3 |
| Host code tests (5 hosts, unknown KeyError) | Task 1.1 |
| Round-trip + parse-invalid + overflow tests | Task 1.1 |
| `assign_kcgp_ids.py` CLI + edge cases | Task 2.3 |
| Sort order (verdict priority + natural chr + pos) | Task 2.3 |
| Ordinal assignment to CANDIDATE only | Task 2.3 |
| `kcgp_mapping.tsv` schema | Task 2.3 |
| Edge-case tests (no s05, no reports, missing event, unknown host) | Task 2.1 |
| `load_kcgp_mapping` soft-fail loader | Task 3.3 |
| PDF cover prepend with `⚠ tentative_v0` warning | Task 3.4 |
| Cover-tests (with/without/corrupt mapping, dash-only) | Task 3.1 |
| `config.yaml` `event:` field | Task 4.1 |
| 3-sample acceptance | Tasks 4.2-4.4 |
| Migration path doc | Spec only — implementation is one-shot replacement of `HOST_CODES` + `build_kcgp_id` body, no plan task needed for v1 swap |

## Notes for the Implementer

- **Read the spec first.** The 3-state verdict semantics, sort key, and PDF integration model are documented there.
- **Verbatim policy.** This codebase prefers no docstrings unless WHY is non-obvious. The two new util files have meaningful docstrings on the public functions only.
- **No new dependencies.** `pyyaml` is already in the redgene env (used by the pipeline).
- **Test isolation.** All tests use `tmp_path` and write fixture TSVs / configs locally. No real BAM access.
- **Existing scaffold protection.** `tests/test_insertion_pdf.py` (18 tests) calls `_page_cover(pdf, name, audit)` without `sample_dir`. Default `sample_dir=None` keeps that signature backward-compatible — verify in Step 3.6.
- **Commit cadence.** Tasks 1-3 each ship as a single commit. Task 4 has a config commit (4.5) and the operational steps (4.2-4.4) emit artifacts but no commits.
