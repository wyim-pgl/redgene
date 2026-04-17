#!/usr/bin/env python3
"""Chain-of-custody logger scaffold (Issue #8).

Writes an append-only JSONL audit trail to
``results/<sample>/chain_of_custody.jsonl``. Each line is a single JSON
object following the schema documented in
``docs/architecture/chain_of_custody.md``:

    {"ts": "ISO8601", "sample": "...", "step": "s03",
     "event": "start|end|error|<custom>",
     "input_sha256": "...", "output_sha256": "...",
     "cmd": "...", "user": "...", "exit_code": 0,
     "wall_time_sec": 1234}

Usage (as a context manager wrapping any pipeline step):

    from scripts.util.coc_logger import CocLogger, sha256_file

    with CocLogger(path=Path("results/S1/chain_of_custody.jsonl"),
                   sample="S1", step="s03") as cl:
        cl.log("input", {"input_sha256": sha256_file(r1)})
        subprocess.run([...], check=True)
        cl.log("output", {"output_sha256": sha256_file(out_bam)})

The context manager auto-emits ``start`` on enter and ``end`` (or
``error`` on exception) on exit, with ``wall_time_sec`` populated. Wire-in
into ``run_pipeline.py.build_step_cmd`` is deferred to a separate v1.1 PR
(see the architecture doc for the proposal).
"""
from __future__ import annotations

import getpass
import hashlib
import json
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256_file(path: Path, *, chunk_size: int = 1 << 16) -> str:
    """Return the hex SHA-256 digest of a file, reading in chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _user() -> str:
    """Best-effort current-user string (env → getpass fallback)."""
    return os.environ.get("USER") or os.environ.get("LOGNAME") \
        or getpass.getuser() or "unknown"


def _now_iso() -> str:
    """UTC ISO-8601 timestamp compatible with ``datetime.fromisoformat``."""
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


class CocLogger:
    """Append-only chain-of-custody JSONL logger (context manager)."""

    def __init__(self, path: Path, sample: str, step: str) -> None:
        self.path = Path(path)
        self.sample = sample
        self.step = step
        self._start_ts: float | None = None
        self._fh = None  # lazily-opened file handle

    # -- context manager ----------------------------------------------------

    def __enter__(self) -> "CocLogger":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        # Open append-only so parallel steps cannot corrupt each other.
        self._fh = open(self.path, "a", encoding="utf-8")
        self._start_ts = time.time()
        self._write_event("start", {})
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        wall = None
        if self._start_ts is not None:
            wall = time.time() - self._start_ts
        if exc_type is None:
            self._write_event("end", {"exit_code": 0, "wall_time_sec": wall})
        else:
            err_msg = f"{exc_type.__name__}: {exc}"
            self._write_event(
                "error",
                {
                    "exit_code": 1,
                    "wall_time_sec": wall,
                    "error": err_msg,
                    "cmd": err_msg,  # mirrored so grep-by-cmd also catches errors
                    "traceback": "".join(traceback.format_exception(exc_type, exc, tb))[-4000:],
                },
            )
        if self._fh is not None:
            self._fh.close()
            self._fh = None
        # Propagate exceptions (do not swallow)
        return False

    # -- public API ---------------------------------------------------------

    def log(self, event: str, payload: dict[str, Any] | None = None) -> None:
        """Append a custom-event line with the given JSON-serializable payload."""
        self._write_event(event, payload or {})

    # -- internal -----------------------------------------------------------

    def _write_event(self, event: str, payload: dict[str, Any]) -> None:
        rec = {
            "ts": _now_iso(),
            "sample": self.sample,
            "step": self.step,
            "event": event,
            "user": _user(),
        }
        rec.update(payload)
        line = json.dumps(rec, sort_keys=False, separators=(",", ":"))
        if self._fh is None:  # direct call outside a with-block (should be rare)
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.path, "a", encoding="utf-8") as fh:
                fh.write(line + "\n")
        else:
            self._fh.write(line + "\n")
            self._fh.flush()


def main(argv: list[str] | None = None) -> int:
    """CLI shim so the logger can be invoked from bash (`coc_logger.py log ...`).

    Minimal; the primary interface is the Python context manager.
    """
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--path", type=Path, required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--step", required=True)
    p.add_argument("--event", required=True)
    p.add_argument("--payload", default="{}",
                   help="JSON payload (default: empty object)")
    args = p.parse_args(argv)

    try:
        payload = json.loads(args.payload)
    except json.JSONDecodeError as exc:
        print(f"[coc] bad --payload: {exc}", file=sys.stderr)
        return 2

    logger = CocLogger(path=args.path, sample=args.sample, step=args.step)
    logger.log(args.event, payload)  # type: ignore[arg-type]
    return 0


if __name__ == "__main__":
    sys.exit(main())
