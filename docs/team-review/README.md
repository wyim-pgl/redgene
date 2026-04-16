# RedGene Team Review — 2026-04-16

**Team:** `redgene-review-2026-04-16` (5 teammates + lead)
**Base commit:** `729835d` (main)
**Goal:** RedGene 파이프라인을 한국 격리소 GMO/LMO 검출 assay용 production-grade로 완성하기 위한 팀 리뷰 산출물 인덱스.

---

## 산출물 (6종)

| # | 파일 | 저자 | 줄수 | 역할 |
|---|------|------|------|------|
| 1 | [`team-consensus.md`](team-consensus.md) | won_yim (PI) | 267 | **최종 합의문.** MVP 12 task (T1-T12), Phase 2/v2.0 이관 목록, AC-1~8 현재 vs 목표 vs 도달, regulatory 제출 선언문 초안. Round 4 구현 시 권위 판정 기준. |
| 2 | [`gap-analysis.md`](gap-analysis.md) | compute_eng | 158 | 7 샘플 TP/FP/FN 매트릭스, AC-1~8 정량화, 병목 분해, gap closure 기여도 매트릭스, 회귀 리스크 매트릭스. MVP 후 production 도달률 **75% (6/8 AC)**. |
| 3 | [`refactor-roadmap.md`](refactor-roadmap.md) | bio_king | 138 | Phase 1 (MVP ≤2주, P1-1~P1-11) / Phase 2 (v1.1 1-2개월, P2-1~P2-9) / Phase 3 (v2.0 grant). 각 Phase task 표 + critical path + risk register + 10일 day-by-day. |
| 4 | [`test-strategy.md`](test-strategy.md) | bio_king | 475 | 현 10 pytest 분석, Phase 1 13건 추가 (compute_verdict 7 + run_blastn 3 + 기타 3), synthetic spike-in 설계, regression gate TSV, GitHub Actions CI 2-job, auto-revert 정책. |
| 5 | [`element-db-expansion.md`](element-db-expansion.md) | gmo_expert | 469 | MUST-HAVE 15-seq subset (P0 6 + P1 8 + P2 1), build_common_payload.sh subregion patch (BUG-9 완전 해결), canonical_triplets config schema, build_host_endogenous_bed.sh workflow, DB governance (cd-hit-est/md5/CHANGELOG). |
| 6 | [`v2-grant-scope.md`](v2-grant-scope.md) | haibao (옵셔널) | 540 | v1.0 7개 fundamental limitations, Level 2 Spec (dual-anchor + PE discordancy, v1.1 candidate), Level 3 Spec (GRIDSS + MCscan + long-read, v2.0 grant), k-mer 실험 설계, ROI vs Risk 분석, funding pitch (3 FTE-year / 18개월). |

---

## Task 명명 정합성 (T* / P1-* / W* 매핑)

세 문서에서 다른 명명을 씀 — Round 4 구현 시 이 mapping을 참조:

| team-consensus T* | refactor-roadmap P1-* | won_yim_round3 W* | 내용 |
|-------------------|------------------------|---------------------|------|
| T1 | P1-10 | W1 | AC-6 audit header 4-field |
| T2 | (measurement) | W5 (선결 측정) | `s05_stats.txt` round-3 수렴율 측정 |
| T3 | P1-7 | W5 | `assemble_insert(max_rounds=3 또는 5)` 수정 |
| T4 | (DB 작업) | W2 | SpCas9+sgRNA+AtYUCCA6+UGT72E3+P0/P1 seqs 추가 |
| T5 | P1-3 + P1-9 | W6 | cd-hit-est dedup + 4-way source tag |
| T6 | P1-1 + P1-5 | W3 | VerdictRules dataclass + compute_verdict 분리 |
| T7 | (부분 분할) | (T8 전제) | 최소 4-way s05 모듈 분할 boundary |
| T8 | (신규) | W4 | `run_pipeline.py` 레벨 per-site fan-out 임시안 |
| T9 | (신규 BED) | W7 | element 100%-ortholog pre-mask BED script + rationale |
| T10 | (site_discovery patch) | (W7 통합) | s05 Phase 1 BED intersect + `FALSE_NEGATIVE_MASKED` |
| T11 | (재검증) | W8 | UGT72E3 + AtYUCCA6 재제출 + 6 샘플 재검증 |
| T12 | (PoC) | (보너스) | rice_G281 s04 minimap2 PoC |

**공수 합산:** team-consensus 기준 ~8 person-day / refactor-roadmap 기준 11 pd (P1-7이 Phase 2 이관 시 10.5 pd). 실제 구현 시 team-consensus.md의 T1-T12가 권위 목록.

---

## 핵심 수치 (6 deliverable cross-reference)

| 지표 | 현재 | MVP 후 예상 | 근거 문서 |
|------|------|-------------|-----------|
| AC-1 Sensitivity (GT recall) | 5/6 (83%) | 6/8 ~ 7/8 | gap-analysis §2, team-consensus §6 |
| AC-2 Specificity (max FP/샘플) | 8 (cuc_225) | ≤ 5 전 샘플 | gap-analysis §2, team-consensus §6 |
| AC-4 Turnaround (UGT72E3) | 48h+ TIMEOUT | ~3-4h | gap-analysis §3, team-consensus §3 |
| AC-6 Audit fields | 1/4 (commit hash) | 4/4 | team-consensus §2.1, element-db §5 |
| s05 BLAST 호출/sample | ~735 (UGT72E3) | ~200-300 (DB 정규화 후) | compute_eng_round2 §Q1 |
| Per-site wall time (soy) | 28 min | ~22 min (round 축소) | compute_eng_round2 §ROI |
| BLAST refactor ROI | — | -5 ~ -7% | compute_eng_round2 §Q1 (실측 수정) |
| Array 병렬화 ROI | — | -92% (UGT72E3) | compute_eng_round1 §3.A |
| MVP 후 production 도달 | — | 75% (6/8 AC) | gap-analysis §6 결론 |

---

## Round 진행 요약

**Round 1 (독립 분석):** 5 teammate 병렬. 총 1,552줄 노트 (notes/ 참조).
**Round 2 (cross-talk):** 4건 합의 (C-1~C-4) + 1건 충돌 (minimap2 Phase 1 교체, compute_eng 반대) + haibao Level 2/3 v2.0 grant 이관.
**Round 3 (synthesis):** won_yim Round 3 판정 + 6 deliverable 작성.
**Round 4 (리뷰 + cleanup):** 리드가 본 index 작성 + 팀 shutdown.

---

## notes/ 인덱스 (Round 1-3 원본 노트, 13 파일)

| 저자 | Round 1 | Round 2 | Round 3 |
|------|---------|---------|---------|
| bio_king | `bio_king_round1.md` (263) | `bio_king_round2.md` (40) | (deliverable 직접 작성) |
| gmo_expert | `gmo_expert_round1.md` (221) | `gmo_expert_round2.md` (127) | (deliverable 직접 작성) |
| compute_eng | `compute_eng_round1.md` (350) | `compute_eng_round2.md` (254) | (deliverable 직접 작성) |
| haibao | `haibao_round1.md` (556) | `haibao_round2.md` (202) | (deliverable 직접 작성) |
| won_yim | `won_yim_round1.md` (162) | `won_yim_round2.md` (274) | `won_yim_round3.md` (230) |

**총 노트 볼륨:** 2,679줄.

---

## Unresolved Items (Round 4 이후 판정)

team-consensus §2.2 + won_yim §Unresolved 기반:

1. **rice_G281 s04 minimap2 PoC** — won_yim 조건부 승인 (T12). 승격 기준 4개: Chr3:16,439,674 CAND 유지 / Phase 1 site ±20% / s04 wall time -30% / MAPQ 분포 유지 (BUG-7 class 방어). 실패 시 자동 기각, v1.1 로드맵 이관.
2. **AC-7 Coverage robustness** — 9 subsampled variants 실행은 Phase 2 이관. v1.0에는 포함 안 됨.
3. **BUG-7 microhomology MAPQ 0-19 회귀 체크리스트** — T11 재검증 시 rice Chr3:16,439,674 / A2_3 ch01:91,002,744 보전 필수.
4. **soybean exploratory 2 샘플의 AC-2 판정** — GT 없이 FP budget ≤ 5 달성 가능한지 T11 결과 보고 판정.

---

## Sign-off

본 Team Review는 Round 1-4를 통해 5 teammate + lead의 합의로 종결. MVP v1.0 release 구현은 team-consensus.md의 T1-T12를 권위 목록으로 이번 주(2026-04-20~24) 수행. 결과는 `docs/team-review/implementation-log.md` (Round 4 이후 생성) 에 기록 예정.

**Lead:** team-lead (Claude via Claude Code)
**PI 승인:** won_yim
**Date:** 2026-04-16
