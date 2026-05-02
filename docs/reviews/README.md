# Code Reviews — Archived

Historical code reviews captured as part of the project's review-as-you-go practice.

| Date | Review | Status |
|------|--------|--------|
| 2026-04-10 | [Pipeline reorganization (Step 9 → Step 5 centric)](2026-04-10-pipeline-reorg.md) | Triaged 2026-05-02 — all critical findings fixed by Issue #4 refactor; 1 minor item remained, cleaned up in commit `64e7d67`. |
| 2026-04-10 | [s09_targeted_assembly.py refactor (commit e11ef0c)](2026-04-10-s09-refactor.md) | Triaged 2026-05-02 — all crash bugs (NameError, tuple-arity) fixed during s09→s05 module split. |
| 2026-04-29 | [Issue #4 + #6 branch review (commits 46aaf40..39a123f)](2026-04-29-branch-review.md) | Triaged 2026-05-02 — `_apply_canonical_override` claim was already fulfilled in `scripts/s05/verdict.py:319`; 4 minor items addressed in commit `64e7d67`. |

Triage methodology: each finding is checked against current code (state ∈ {FIXED, STALE, ACTIVE, WONTFIX}). Active findings under MINOR severity are batched into polish commits. Anything CRITICAL/IMPORTANT becomes its own issue.
