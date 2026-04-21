"""DAG no-cycle guard for the ``scripts/s05/`` package (Issue #4, Session 1).

The multi-session refactor documented in
``docs/superpowers/specs/2026-04-19-s05-module-split-design.md`` §2.2 fixes a
strict layering order:

    primitives  ->  {site_discovery, classify}  ->  read_extraction
                                                          |
                                                          v
                                                     assembly
                                                          |
                                                          v
                                                     annotation
                                                          |
            +-------+--------+--------+--------+----------+
            v       v        v        v        v
        filter_a filter_b filter_c filter_d
            +-------+--------+--------+----------+
                                                 |
                                                 v
                                             verdict  <-- config_loader
                                                 |
                                                 v
                                              report
                                                 |
                                                 v
                                       fanout_orchestrator

A later-stage module may import from an earlier (or same) stage; the reverse
is forbidden. This test walks the AST of every ``scripts/s05/*.py`` file and
flags any ``ImportFrom`` whose target module belongs to a strictly later
stage, which would imply a cycle once both sides are fully extracted.

Session 1 lands only a subset of the above modules (``primitives``,
``filter_b_flank``, ``filter_c_chimeric``, plus the already-extracted
``verdict``/``config_loader``). The four remaining shim files
(``annotate_report``, ``classify``, ``per_site``, ``site_discovery``) only
import back from the monolith ``scripts.s05_insert_assembly`` and are placed
at their final stage so later sessions can promote them in-place without
changing this test.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
PKG_DIR = REPO_ROOT / "scripts" / "s05"


# Stage ranking per spec §2.2. Lower rank = earlier in the pipeline.
# Modules at the same stage may freely import each other (e.g., the four
# filters consume the same FilterEvidence surface), but a module must never
# import from a strictly later stage.
STAGE: dict[str, int] = {
    "primitives": 0,
    "site_discovery": 1,
    "classify": 1,
    "read_extraction": 2,
    "per_site": 2,           # legacy shim; same stage as read_extraction (renames in Session 5)
    "assembly": 3,
    "annotation": 4,
    "annotate_report": 4,    # legacy shim; same stage as annotation (retired in Session 7)
    "filter_a_host": 5,
    "filter_b_flank": 5,
    "filter_c_chimeric": 5,
    "filter_d_altlocus": 5,
    "verdict": 6,
    "config_loader": 7,      # imports VerdictRules from verdict
    "report": 8,
    "fanout_orchestrator": 9,
}


def _iter_module_files() -> list[Path]:
    files = [p for p in sorted(PKG_DIR.glob("*.py")) if p.name != "__init__.py"]
    assert files, f"no modules found under {PKG_DIR}"
    return files


def _extract_s05_sibling_imports(src: str) -> list[tuple[str, int]]:
    """Return (sibling_module_name, lineno) for every ``from scripts.s05.X`` or
    ``from .X`` import found in src. Non-s05 imports are ignored.
    """
    tree = ast.parse(src)
    siblings: list[tuple[str, int]] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom):
            continue

        target: str | None = None
        # Relative import: ``from .sibling import X`` or ``from . import sibling``.
        if node.level and node.module:
            target = node.module.split(".")[0]
        elif node.level and not node.module:
            # ``from . import X``
            for alias in node.names:
                if alias.name in STAGE:
                    siblings.append((alias.name, node.lineno))
            continue
        elif node.module:
            # Absolute import: must start with scripts.s05.X.
            parts = node.module.split(".")
            if len(parts) >= 3 and parts[0] == "scripts" and parts[1] == "s05":
                target = parts[2]

        if target and target in STAGE:
            siblings.append((target, node.lineno))
    return siblings


@pytest.mark.parametrize("module_file", _iter_module_files(), ids=lambda p: p.name)
def test_module_imports_respect_dag(module_file: Path) -> None:
    """No module may import from a strictly later DAG stage (cycle guard)."""
    module_name = module_file.stem
    assert module_name in STAGE, (
        f"Unknown module {module_name}; add it to STAGE in "
        f"tests/test_s05_import_dag.py with the appropriate rank."
    )
    src = module_file.read_text()
    siblings = _extract_s05_sibling_imports(src)
    my_stage = STAGE[module_name]
    violations = [
        (target, lineno) for target, lineno in siblings
        if STAGE[target] > my_stage
    ]
    assert not violations, (
        f"{module_file.name} imports from later-stage sibling(s): "
        + ", ".join(
            f"{target} (stage {STAGE[target]}) at line {lineno}"
            for target, lineno in violations
        )
        + f"; this module is stage {my_stage}."
    )


def test_required_session1_modules_present() -> None:
    """Sanity check: every Session-1 ranked module file exists under scripts/s05/."""
    present = {p.stem for p in _iter_module_files()}
    # Only enforce existence for the Session-1 modules; later sessions will
    # add the rest. STAGE can hold forward-looking entries without failing.
    required = {
        "primitives", "filter_b_flank", "filter_c_chimeric",
        "verdict", "config_loader",
    }
    missing = required - present
    assert not missing, f"expected Session-1 modules missing: {missing}"
