"""Guard for issue #18: assemble_insert scratch dirs must be site-scoped.

The per-site SLURM fan-out (scripts/submit_s05_array.sh, --phase 2_3) runs one
process per site, all sharing a single `step_dir` passed to
`assemble_insert(workdir=step_dir)`. If any scratch path is derived directly
from `workdir` with a site-independent name (e.g. `workdir/_mm2_r1`), two
concurrent sites collide on it — one task's cleanup removes a file another is
writing, and worse, a task can map its reads against another site's contig.
This is not hypothetical: it crashed 26 of 106 tasks on soybean_UGT72E3
(SLURM array 5811054) with `_mm2_r1/_mm2ext.bam: No such file or directory`.

Scratch must instead be namespaced by site. The final `{site_id}_insert.fasta`
still belongs in `workdir` so Phase 4 can find it — that line is allowed to be
`workdir`-scoped and is excluded below.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

ASM = Path(__file__).resolve().parents[1] / "scripts" / "s05" / "assembly.py"

# Scratch basenames that concurrent sites would collide on if scoped to workdir.
_SCRATCH_TOKENS = ("_mm2", "_pilon", "_ssake", "_foreign_refine", "_host_term")


def _source_lines() -> list[tuple[int, str]]:
    return list(enumerate(ASM.read_text().splitlines(), start=1))


# Only these two functions receive the shared step_dir as `workdir` under the
# fan-out; the assembly helpers (_minimap2_extend, _ssake_extend,
# check_host_termination, pilon_fill) take an already-scoped scratch dir from
# their caller, so a fixed name off their local `workdir` is correct there.
_SHARED_WORKDIR_FUNCS = ("assemble_insert", "refine_with_foreign_reads")


def _shared_workdir_line_ranges() -> list[tuple[int, int]]:
    tree = ast.parse(ASM.read_text())
    ranges: list[tuple[int, int]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name in _SHARED_WORKDIR_FUNCS:
            ranges.append((node.lineno, node.end_lineno or node.lineno))
    return ranges


# Interpolations that namespace a path by site — a `workdir/...` scratch path is
# safe only if it carries one of these.
_SITE_NAMESPACE_TOKENS = ("site.site_id", "site_id", "insert_fasta.stem")


def test_no_workdir_scoped_scratch_paths():
    """A `workdir/<scratch>` path is allowed only if namespaced by site.

    A fixed name like `workdir/"_foreign_refine"` collides across concurrent
    sites; `workdir/f"_foreign_refine_{insert_fasta.stem}"` does not. The final
    insert FASTA (`workdir/f"{site.site_id}_insert.fasta"`) is site-namespaced by
    construction and therefore never flagged.
    """
    offenders: list[str] = []
    pat = re.compile(r"\bworkdir\s*/\s*f?[\"']([^\"']+)[\"']")
    ranges = _shared_workdir_line_ranges()
    for lineno, line in _source_lines():
        if not any(lo <= lineno <= hi for lo, hi in ranges):
            continue  # helper bodies get a scoped scratch dir from their caller
        for m in pat.finditer(line):
            literal = m.group(1)
            if not any(tok in literal for tok in _SCRATCH_TOKENS):
                continue
            if any(ns in literal for ns in _SITE_NAMESPACE_TOKENS):
                continue  # site-scoped → safe
            offenders.append(f"{ASM.name}:{lineno}: workdir/{literal!r}")
    assert not offenders, (
        "scratch paths must be site-scoped (issue #18), found workdir-scoped "
        "with no site namespace:\n" + "\n".join(offenders)
    )


def test_check_host_termination_gets_site_scoped_dir():
    """`check_host_termination(...)` must not be passed the bare `workdir`.

    It writes fixed-name `_host_term.fa`/`_host_term.paf`, so a shared workdir
    collides across concurrent sites.
    """
    tree = ast.parse(ASM.read_text())
    bad: list[int] = []
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "check_host_termination"
        ):
            last = node.args[-1] if node.args else None
            if isinstance(last, ast.Name) and last.id == "workdir":
                bad.append(node.lineno)
    assert not bad, (
        f"check_host_termination passed bare workdir at lines {bad} "
        "— must be a site-scoped scratch dir (issue #18)"
    )


def test_assemble_insert_creates_a_site_scoped_scratch_dir():
    """assemble_insert must materialise a scratch dir named by site_id."""
    src = ASM.read_text()
    # e.g. scratch = workdir / f"_scratch_{site.site_id}"
    assert re.search(r"workdir\s*/\s*f[\"']_scratch_\{site\.site_id\}[\"']", src), (
        "assemble_insert should derive a per-site scratch dir like "
        "`workdir / f\"_scratch_{site.site_id}\"` (issue #18)"
    )
