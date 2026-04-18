# Chain-of-custody log (Issue #8 scaffold)

Every RedGene pipeline step will eventually emit an append-only JSONL trail
to `results/<sample>/chain_of_custody.jsonl`. The scaffold lives at
[`scripts/util/coc_logger.py`](../../scripts/util/coc_logger.py); wiring
into `run_pipeline.py.build_step_cmd` is an explicit v1.1 follow-up PR and
is **not** part of this commit.

## Schema

One JSON object per line. Required keys are marked `*`.

| field              | type     | notes                                                     |
|--------------------|----------|-----------------------------------------------------------|
| `ts` *             | string   | ISO-8601 UTC timestamp (microsecond precision)            |
| `sample` *         | string   | Sample key matching `config.yaml`                         |
| `step` *           | string   | e.g. `s03`, `s04`, `s05`                                  |
| `event` *          | string   | `start` / `end` / `error` / arbitrary custom tag          |
| `user` *           | string   | `$USER` / `getpass.getuser()` fallback                    |
| `input_sha256`     | string   | Hex digest of a primary input (populated by caller)       |
| `output_sha256`    | string   | Hex digest of a primary output (populated by caller)      |
| `cmd`              | string   | Full command line (or exception mirror on error)          |
| `exit_code`        | int      | 0 on `end`, non-zero on `error`                           |
| `wall_time_sec`    | float    | Populated on `end` / `error`                              |
| `error`            | string   | Error message (only on `event == "error"`)                |
| `traceback`        | string   | Truncated to last 4 KB (only on `event == "error"`)       |
| `coc_chain_id`     | string   | **Proposed** — hash linking back to `audit_header.json`   |

### `coc_chain_id` cross-link proposal

Goal: let a reviewer run one `sha256` check to prove that a given JSONL trail
was produced by the same pipeline run described in `audit_header.json`.

Approach:

1. At step 1 entry, compute
   `coc_chain_id = sha256(audit_header.json contents)[:16]`.
2. Store it inside `audit_header.json` as a new field
   (`"coc_chain_id": "<hex>"`) and on every JSONL line.
3. Reviewer: `grep coc_chain_id chain_of_custody.jsonl | sort -u` must return
   exactly one value, and it must equal the value in `audit_header.json`.

This field is **not** emitted by the current scaffold (out of scope). The
v1.1 wire-in PR will add it to both the logger and the audit-header writer
(`scripts/_write_audit_header.py`).

## Usage

### As a context manager (primary interface)

```python
from pathlib import Path
from scripts.util.coc_logger import CocLogger, sha256_file

coc_path = Path(f"results/{sample}/chain_of_custody.jsonl")

with CocLogger(path=coc_path, sample=sample, step="s03") as cl:
    cl.log("input", {"input_sha256": sha256_file(r1)})
    run_subprocess(["bwa", "mem", ...])
    cl.log("output", {"output_sha256": sha256_file(out_bam)})
```

On clean exit this emits 4 lines: `start`, `input`, `output`, `end`.
On exception the `end` line is replaced with an `error` line that carries
the exception type, message, and a truncated traceback.

### From bash (secondary interface)

```bash
python scripts/util/coc_logger.py \
    --path results/S1/chain_of_custody.jsonl \
    --sample S1 --step s03 \
    --event checkpoint \
    --payload '{"note":"post-fastp","pairs":48213221}'
```

## Invariants (enforced by tests)

- **Append-only:** two successive `with CocLogger(...)` blocks produce 4
  lines, never 2. No truncation.
- **ISO-8601:** `ts` parses via `datetime.fromisoformat`.
- **User populated:** `user` is never empty.
- **Error path:** an exception inside the `with` block writes exactly one
  extra `error` line with `exit_code != 0`, then re-raises.
- **SHA256 helper:** `sha256_file()` hex-digest matches `hashlib.sha256()`.

## Wire-in plan (deferred to v1.1)

- Patch `run_pipeline.py.build_step_cmd` so every step is wrapped in
  `with CocLogger(...)`.
- Compute `coc_chain_id` once per run and persist into
  `results/<sample>/audit_header.json`.
- Add a `tools/verify_coc.py` that re-checks the chain_id and re-hashes
  the listed input/output files to verify nothing was tampered with
  after the fact.
- Integrate into `scripts/reports/insertion_pdf.py` section (7) so the
  PDF appendix links the CoC chain-id.
