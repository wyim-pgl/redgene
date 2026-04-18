#!/usr/bin/env python3
"""Chain-of-custody integrity verifier (Issue #8).

Reads ``results/<sample>/chain_of_custody.jsonl`` and checks:

  (a) Monotone timestamps — every record's ``ts`` must be >= its predecessor.
  (b) Pre/post pairing — every ``pre`` event for step S must have a matching
      ``post`` event for the same step S; no ``pre`` without ``post`` or vice
      versa.
  (c) Hash continuity — where step N's outputs are step N+1's inputs
      (``output_sha256`` entries appear as ``input_sha256`` entries in the
      next step's ``pre`` record), the digests must match.

Exits 0 on success, 1 on any integrity failure.

Usage::

    python tools/verify_coc.py results/rice_G281/chain_of_custody.jsonl
    python tools/verify_coc.py --sample rice_G281 --outdir results
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Core verifier
# ---------------------------------------------------------------------------

def load_records(jsonl_path: Path) -> list[dict[str, Any]]:
    """Parse all JSONL records from *jsonl_path*."""
    records: list[dict[str, Any]] = []
    with open(jsonl_path, encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            raw = raw.strip()
            if not raw:
                continue
            try:
                records.append(json.loads(raw))
            except json.JSONDecodeError as exc:
                raise ValueError(f"Line {lineno}: invalid JSON — {exc}") from exc
    return records


def verify_monotone_timestamps(records: list[dict[str, Any]]) -> list[str]:
    """Return a list of error messages for any non-monotone timestamp pair."""
    errors: list[str] = []
    prev_ts: datetime | None = None
    prev_idx: int = 0
    for i, rec in enumerate(records):
        raw_ts = rec.get("ts", "")
        try:
            ts = datetime.fromisoformat(raw_ts)
        except ValueError:
            errors.append(f"Record {i}: unparseable timestamp {raw_ts!r}")
            continue
        if prev_ts is not None and ts < prev_ts:
            errors.append(
                f"Record {i} (step={rec.get('step')}, event={rec.get('event')}): "
                f"timestamp {raw_ts} is earlier than record {prev_idx} ({prev_ts.isoformat()})"
            )
        prev_ts = ts
        prev_idx = i
    return errors


def verify_pre_post_pairs(records: list[dict[str, Any]]) -> list[str]:
    """Return errors for any step whose ``pre`` and ``post`` events are unpaired.

    A valid run for step S produces at least one ``pre`` record followed by
    at least one ``post`` record (the CocLogger start/end events bracket them).
    We also require that no ``post`` appears before its matching ``pre``.
    """
    errors: list[str] = []
    # Track open pre-entries per step
    open_steps: dict[str, int] = {}  # step -> count of unmatched pre records

    for i, rec in enumerate(records):
        event = rec.get("event", "")
        step = rec.get("step", "?")

        if event == "pre":
            open_steps[step] = open_steps.get(step, 0) + 1
        elif event == "post":
            if open_steps.get(step, 0) == 0:
                errors.append(
                    f"Record {i}: 'post' for step {step!r} has no preceding 'pre'"
                )
            else:
                open_steps[step] -= 1

    # Any step with unmatched pre entries
    for step, count in open_steps.items():
        if count > 0:
            errors.append(
                f"Step {step!r}: {count} 'pre' event(s) without a matching 'post'"
            )
    return errors


def verify_hash_continuity(records: list[dict[str, Any]]) -> list[str]:
    """Return errors where step N's output SHA-256 differs from step N+1's input.

    We collect (step, output_sha256) from all ``post`` records, and
    (step, input_sha256) from all ``pre`` records, then check that whenever
    the same filename appears in consecutive step outputs→inputs the digests
    match.

    Only filenames that appear in *both* the predecessor's output and the
    successor's input are checked; additional files in either set are ignored
    (steps may produce many outputs but only pass some downstream).
    """
    errors: list[str] = []

    # Ordered list of (step, event, sha_dict) for pre/post records
    timeline: list[tuple[str, str, dict[str, str]]] = []
    for rec in records:
        event = rec.get("event", "")
        step = rec.get("step", "?")
        if event == "pre":
            sha = rec.get("input_sha256") or {}
            if isinstance(sha, dict):
                timeline.append((step, "pre", sha))
        elif event == "post":
            sha = rec.get("output_sha256") or {}
            if isinstance(sha, dict):
                timeline.append((step, "post", sha))

    # Walk consecutive post→pre pairs
    for j in range(len(timeline) - 1):
        step_j, event_j, sha_j = timeline[j]
        step_k, event_k, sha_k = timeline[j + 1]
        if event_j == "post" and event_k == "pre":
            # sha_j = previous step's outputs, sha_k = next step's inputs
            for fname, digest_j in sha_j.items():
                if fname in sha_k:
                    digest_k = sha_k[fname]
                    if digest_j != digest_k:
                        errors.append(
                            f"Hash mismatch for '{fname}': "
                            f"step {step_j!r} output={digest_j[:12]}… "
                            f"but step {step_k!r} input={digest_k[:12]}…"
                        )
    return errors


def verify(jsonl_path: Path, *, verbose: bool = False) -> bool:
    """Run all three checks and print any errors.  Returns True if clean."""
    records = load_records(jsonl_path)
    if verbose:
        print(f"[coc-verify] {jsonl_path}: {len(records)} records")

    all_errors: list[str] = []
    all_errors.extend(verify_monotone_timestamps(records))
    all_errors.extend(verify_pre_post_pairs(records))
    all_errors.extend(verify_hash_continuity(records))

    if all_errors:
        print(f"[coc-verify] FAIL — {len(all_errors)} error(s) in {jsonl_path}:",
              file=sys.stderr)
        for err in all_errors:
            print(f"  • {err}", file=sys.stderr)
        return False

    if verbose:
        print(f"[coc-verify] OK — all checks passed.")
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "jsonl", nargs="?", type=Path, metavar="JSONL",
        help="Direct path to chain_of_custody.jsonl",
    )
    group.add_argument(
        "--sample", type=str,
        help="Sample key (used with --outdir to locate the JSONL file)",
    )
    p.add_argument(
        "--outdir", type=Path, default=Path("results"),
        help="Pipeline output root (default: results)",
    )
    p.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print summary even on success",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.sample:
        jsonl = args.outdir / args.sample / "chain_of_custody.jsonl"
    else:
        jsonl = args.jsonl

    if not jsonl.exists():
        print(f"[coc-verify] ERROR: file not found: {jsonl}", file=sys.stderr)
        return 2

    ok = verify(jsonl, verbose=args.verbose)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
