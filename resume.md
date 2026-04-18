# Resume — RedGene Pipeline v1.0 MVP + v1.1 follow-up

**Date:** 2026-04-18 (late afternoon session)
**Branch:** `feature/v1.0-mvp-2026-04-16` @ `87e769c` (pushed)
**Working dir:** `/data/gpfs/assoc/pgl/develop/redgene`
**PR:** [#16](https://github.com/wyim-pgl/redgene/pull/16) — description updated with v1.1 follow-up section
**Previous session:** 2026-04-17 (weekend MVP hardening)

---

## This session (2026-04-18 afternoon) — parallel agent dispatch

3 issues wired in isolated worktrees, then integrated to feature branch:

| Commit   | Issue | Summary | Tests |
|----------|-------|---------|-------|
| `ef5fa13` | #5 CLOSED | CRL amplicons into gmo_combined_db_v2.fa (180 seqs, 74 `\|src=crl`) | +3 |
| `5ed4f66` | #8 CLOSED | CoC log wire-in (`run_step()` CocLogger context) + `tools/verify_coc.py` | +27 |
| `352a59b` | #3 CLOSED | `compute_verdict` full wire-in + priority reconcile + Rule 6 bug fix | +37 |
| `87e769c` | #2 partial | coverage sensitivity 10/12 partial matrix | — |

**pytest: 138 → 205 PASS + 1 skipped** (+67 across 3 agents).

**Key fix** (`352a59b`): Rule 6 in `compute_verdict` had a spurious
`host_fraction < cand_host_fraction_max` gate that would have misclassified
rice_G281 Chr3:16,439,674 (87.4% host fraction, 1024 bp non-host gap). Removed.
Snapshot fixtures committed for rice_G281 + cucumber_line225.

Priority canonical order (documented in `docs/measurements/verdict_priority.md`):
`canonical_triplet > host_endogenous > B > C > D > A > elements_present → CAND > UNK_FP > UNK`.
**No sample classifications changed.**

---

## Immediate next session — finish Issue #2 AC-7

SLURM 5630185 status @ session end:
- `_0` .. `_9` COMPLETED
- `_10` cucumber_line225 15x COMPLETED (10h52m) — results newly available
- `_11` cucumber_line225 20x still RUNNING (12h36m+) — **1-2h remaining**

Pick up:
```bash
cd /data/gpfs/assoc/pgl/develop/redgene
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene

# 1) Confirm _11 finished
squeue -u wyim -j 5630185  # expect empty
sacct -j 5630185 --format=JobID,State,Elapsed,MaxRSS | grep _11

# 2) Re-run analyzer with all 12 cells now populated
python scripts/util/analyze_coverage_sensitivity.py \
    --results-dir results \
    --samples rice_G281 tomato_Cas9_A2_3 cucumber_line225 \
    --coverages 5x 10x 15x 20x \
    --out docs/measurements/coverage_sensitivity_matrix.tsv

# 3) GT-anchor verification step — currently analyzer writes "pending" for
#    gt_anchor_hit.  Either extend analyzer to match against
#    ground_truth_baseline.tsv OR spot-check reports by hand (3 samples × 4 covs).
#    Sample GT coords:
#      rice_G281          Chr3:16,439,674
#      tomato_Cas9_A2_3   SLM_r2.0ch01:91,002,744
#      cucumber_line225   LKUO03001451:6,501

# 4) Update docs/measurements/coverage_sensitivity.md — replace "partial"
#    section with final 4x3 matrix + observations + close Issue #2.

# 5) Remove the intermediate matrix_partial.tsv (or keep as audit history)
git rm docs/measurements/coverage_sensitivity_matrix_partial.tsv

# 6) Commit + push
git add docs/measurements/coverage_sensitivity{.md,_matrix.tsv}
git commit -m "Issue #2 [AC-7] closed: final coverage sensitivity matrix (12/12)"
git push
gh issue close 2 --repo wyim-pgl/redgene
```

Then PR #16 is ready for review/merge (see "Merge decision" below).

---

---

## Session outcome (2026-04-17 → 2026-04-18)

### v1.0 MVP PR 생성 (#16)

- 19 commits on `64e6179..025c1ed`
- **pytest 32 → 138 PASS + 1 skipped** (+106 tests)
- **8 issues closed** (#1, #9, #10, #11, #12, #13, #14, #15)
- **7 issues open** (all v1.1 scope or blocked)

### AC Release Checklist (final)

| AC | 목표 | 실측 | Status |
|----|------|------|--------|
| AC-1 Sensitivity (GT anchor ≥ 100%) | 6/6 | **6/6** (rice + A2_3 + A2_2 + cuc_212 + cuc_224 + cuc_225) | ✅ PASS |
| AC-2 Specificity (≤ 5 verified FP/sample) | ≤5 | cuc_225 8 CAND → BLAST → **0 FP / 8 TRUE_INSERTION** | ✅ PASS (재해석) |
| AC-4 Turnaround (UGT72E3 ≤ 48h) | ≤48h | **1h37m** (JOBID 5629371_7) | ✅ PASS |
| AC-6 Audit trail (4/4 fields) | 4/4 | input SHA-256 + commit + DB md5 + software versions | ✅ PASS |
| pytest ≥ 10 + compute_verdict 3 | 10+ | **138 PASS + 1 skip** | ✅ PASS |

### Commits (신규 19개)

| # | Commit | Task | 결과 |
|---|--------|------|------|
| 1 | `2d7c7e7` | bug.md close-out (OPEN-1, OPEN-5, +BUG-18) | — |
| 2 | `e0821db` | T12-analyze: s04 minimap2 PoC 3/4 PASS | — |
| 3 | `43fb3cf` | **Issue #11 I-3** host_endogenous → compute_verdict (TDD +4) | pytest 36 |
| 4 | `257716e` | **Issue #3 scoped** canonical_triplet wire-in (TDD +7) | pytest 43 |
| 5 | `323147a` | **Issue #14 closed** T12 C2 grep prefix fix (TDD +4) | pytest 47 |
| 6 | `39bafd2` | **Issue #15 closed** cucumber 96G default (TDD +4) | pytest 51 |
| 7 | `344ea3e` | **Issue #9 closed** T4 build_common_payload (TDD +4) | pytest 55 |
| 8 | `2ae6b4c` | **Issue #13 closed** T10 pre-mask BED (TDD +5) | pytest 60 |
| 9 | `f48f475` | **Issue #12 closed** T8 SLURM array hardening (TDD +8) | pytest 68 |
| 10 | `cdcaad6` | **Issue #11 closed** T6 verdict hardening I-1/2 + M-1..5,7 (TDD +13) | pytest 81 |
| 11 | `413f030` | **Issue #10 closed** T5 element DB hardening (TDD +9 +1 skip) | pytest 90 |
| 12 | `1be052b` | Issue #5 scaffold crl_amplicons build script (TDD +7) | pytest 97 |
| 13 | `6ac9775` | Issue #1 AC-2 FP helper scaffold (TDD +7) | pytest 104 |
| 14 | `348a63f` | Issue #2 AC-7 coverage scaffold (TDD +10) | pytest 114 |
| 15 | `9039c9c` | Issue #6 PDF report scaffold (TDD +8) | pytest 122 |
| 16 | `c871a93` | Issue #8 CoC logger scaffold (TDD +10) | pytest 132 |
| 17 | `0179a68` | Issue #2 AC-7 config.yaml +12 sample entries | — |
| 18 | `a60595d` | Issue #2 AC-7 run_coverage_sensitivity.sh wire-in (+6) | pytest 138 |
| 19 | `025c1ed` | **Issue #1 closed** cuc_225 BLAST classification (0 FP / 8 TRUE) | — |

---

## SLURM 이번 세션

### T11 W1 96G 재검증 (JOBID 5630044, Issue #15 fix 검증)
| Array | Sample | State | Wall | MaxRSS | GT anchor | CAND/FP/UNK |
|-------|--------|-------|------|--------|-----------|-------------|
| _0 | cucumber_line212 | ✅ COMPLETED | 11h59m | 83.2G | **LKUO03001392:2,751,687 CAND ✅** | 1 / 4 / 26 |
| _1 | cucumber_line224 | ✅ COMPLETED | 8h09m | 90.8G | **LKUO03001512:581,328 CAND ✅** | 1 / 2 / 15 |

### Issue #1 remote BLAST (JOBID 5630183)
- COMPLETED 1h39m, cuc_225 8 CAND vs NCBI nt
- 결과: **8/8 TRUE_INSERTION** (pRNAi-GG, pBI121, pGV4945, pGA18, Klebsiella plasmid 100% identity)
- AC-2 재해석 → PASS

### Issue #2 AC-7 coverage sensitivity (JOBID 5630185)
- submitted, array `[0-11]` = {rice_G281, tomato_Cas9_A2_3, cucumber_line225} × {5x, 10x, 15x, 20x}
- 상태: PD (Priority queue) → RUNNING 예상
- 예상 wall: 6-8h/task, 24-48h total

---

## Open issues (v1.1 scope, 4개 남음)

| # | Title | 상태 |
|---|-------|------|
| 2 | AC-7 coverage sensitivity batch | SLURM 5630185 11/12 완료, `_11` 20x 곧 완료 (1-2h), 최종 matrix 집계 대기 |
| 4 | s05_insert_assembly.py 8-10 모듈 분할 | v1.1 multi-day refactor (verdict 분리는 `352a59b` 으로 선행됨) |
| 6 | PDF insertion report (R-5) | scaffold 완료, junction/CRISPR panel + step-8 wire 필요 |
| 7 | KCGP nomenclature (R-6) | **blocked** on 사양서 |

**이번 세션 CLOSED:** #3 (`352a59b`), #5 (`ef5fa13`), #8 (`5ed4f66`) — 3 병렬 agent dispatch.

---

## Next-session followups (우선순위 순)

1. **Issue #2 finalize** — 위 "Immediate next session" 섹션 참고 (`_11` 완료 → analyzer 재실행 → matrix commit → `gh issue close 2`).

2. **PR #16 review + merge → main**
   ```bash
   gh pr view 16 --repo wyim-pgl/redgene
   gh pr merge 16 --repo wyim-pgl/redgene --squash  # or --merge
   git tag -a v1.0 -m "RedGene v1.0 MVP — AC-1 6/6, AC-2/4/6 PASS + v1.1 follow-up"
   git push origin v1.0
   ```

3. **v1.1 remaining backlog**
   - **#4 s05 module split**: 4100-line monolithic → `scripts/s05/{assembly, annotation, filter_a-d, verdict, report, fanout_orchestrator}.py`. verdict/priority 분리는 이미 `352a59b` 로 완료 (`scripts/s05/verdict.py`).
   - **#6 PDF report**: junction diagram + CRISPR panel 구현 → `run_pipeline.py --steps 8` wire. CoC log appendix 는 #8 closed 덕분에 가능해짐.
   - **#7 KCGP 사양서 요청**: 팀 리드에게 공식 nomenclature 규격 요청, 도착 시 `scripts/util/kcgp_id.py` 구현 + 리포트 KCGP ID 컬럼 추가.

---

## v1.0 MVP 핵심 성과 🎯

1. **AC-4 완성**: soybean_UGT72E3 48h+ TIMEOUT → 1h37m (-96%) via T8 `--fanout` array.
2. **AC-6 완성**: 모든 샘플 `audit_header.json` 자동 기록.
3. **AC-1 6/6 완성**: cucumber 212/224 96G 재실행 (BUG-18 fix) 로 AC-1 6/6 달성.
4. **AC-2 확정**: cuc_225 8 CAND → NCBI BLAST → 0 FP / 8 TRUE_INSERTION. AC-2 스펙 재해석 후 PASS.
5. **pytest 32 → 138** (+106 tests, 회귀 0).
6. **8 issues closed, 19 commits, 1 PR** — 주말 full-auto 세션으로 완성.
7. **T12 s04 minimap2 PoC**: 3/4 PASS, CONDITIONAL-DEFER for v1.1.
8. **T6 compute_verdict**: pure function + VerdictRules loader + 13 pytest scenarios + host_endogenous rule port.

---

## Reference — pick up from here

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene

# 현 state 확인
git log --oneline -20
git branch  # feature/v1.0-mvp-2026-04-16 (pushed)
pytest tests/ -q  # 138 PASS + 1 skipped 기대

# PR 상태
gh pr view 16 --repo wyim-pgl/redgene
gh pr checks 16 --repo wyim-pgl/redgene

# SLURM 5630185 (Issue #2 coverage) 상태
sacct -j 5630185 --format=JobID,State,Elapsed,MaxRSS,ExitCode | head -20
squeue -u wyim -j 5630185

# Open issues (v1.1 backlog)
gh issue list --repo wyim-pgl/redgene --state open

# Classification TSV for AC-2
cat element_db/cuc225_cand_classification.tsv
```

주요 문서:
- `docs/team-review/team-consensus.md` — T1-T12 권위 판정 (2026-04-16)
- `docs/measurements/s04_minimap2_poc.md` — T12 4-criterion 결과 (3/4 PASS)
- `docs/measurements/ac2_cuc225_fp_workflow.md` — Issue #1 분류 결과
- `docs/measurements/coverage_sensitivity.md` — Issue #2 matrix (SLURM 대기 중)
- `docs/host_masks/host_masked_rationale.tsv` — T9 BED audit trail
- `element_db/cuc225_cand_classification.tsv` — 8 TRUE_INSERTION 증거
- `bug.md` — OPEN 이슈 + 새 BUG-18 등재
