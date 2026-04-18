# Gap Analysis — RedGene 현 파이프라인의 Production 도달률 정량화

**Author:** compute_eng · **Date:** 2026-04-16 · **Team Round:** 3 · **Base commit:** `729835d`

**목적:** "외래 유전자 정체 + 위치" 규제 목표에 대한 현 파이프라인의 달성률을 실측으로 정량하고, won_yim MVP 결정에 따른 예상 도달률을 계산. 실측 근거는 `results/*/s05_insert_assembly/insertion_*_report.txt`, `sacct`, `resume.md`.

---

## 1. 샘플별 TP / FP / FN / UNKNOWN 매트릭스

### 1.1 Ground-truth anchor 샘플 (AC-1 직접 측정 대상)

| Sample | Host | GT locus | GT recall? | CAND total | FP (CAND - TP) | UNKNOWN | FN (missed GT) | 검증 상태 |
|---|---|---|---|---|---|---|---|---|
| **rice_G281** | Osativa 374 Mbp | Chr3:16,439,674 (2-copy H2H) | ✅ **Chr3:16439674 CANDIDATE** | 3 (16439674 + Chr3:29074002 + Chr11:2877409) | 2 (remote BLAST 검증 필요) | 33 | 0 | resume.md PASS |
| **tomato_Cas9_A2_3** | SLM_r2.0 833 Mbp | ch01:91,002,744 | ✅ **ch01 CANDIDATE** | 2 (91002744 + ch06:41781251) | 1 (ch06 review 필요) | 0 | 0 | resume.md PASS |
| **cucumber_line212** | B10v3 332 Mbp | Chr6 single T-DNA (thaumatin) | ✅ **LKUO03001392:2,751,687 CANDIDATE** | 1 | 0 | 26 | 0 | resume.md PASS |
| **cucumber_line224** | B10v3 332 Mbp | Chr2 single T-DNA (G6838 disrupt) | ✅ **LKUO03001512:581,328 CANDIDATE** | 1 | 0 | 15 | 0 | resume.md PASS |
| **cucumber_line225** | B10v3 332 Mbp | Chr2 2-copy + backbone (2 expected) | ✅ **LKUO03001451:6,501 CANDIDATE** | 8 (1 GT + 7 new) | 7 (remote BLAST 검증 필요) | 26 | 0 (추가 copy backbone은 8개 중 1-2개?) | **FP budget 8 > AC-2 MVP 5** ⚠️ |
| **tomato_Cas9_A2_2** | SLM_r2.0 833 Mbp | ch08:65,107,378 | ❓ **재검증 미완료** | ? | ? | ? | ? | resume.md에 항목 없음 (미재검증 GT) |

**GT Recall = 5/6 = 83.3% (tomato_A2_2 재검증 전 기준)**. AC-1 MVP threshold 80%는 **달성**, production 100%는 **미달**.

### 1.2 Exploratory 샘플 (GT 없음 — AC-2, AC-4만 측정)

| Sample | Host | GT status | CAND | FP count | UNKNOWN | Wall time | Verdict 해석 |
|---|---|---|---|---|---|---|---|
| **soybean_AtYUCCA6** | Gmax 1.1 Gbp | GT 미확립 | **0** | 15 | 30 | RUNNING 21h+/24h | 0-CAND 문제 (resume.md:96-105) |
| **soybean_UGT72E3** | Gmax 1.1 Gbp | GT 미확립 | **0** | 15 | 38 | **24h TIMEOUT** (partial 53/105) | 0-CAND + TIMEOUT (AC-4 FAIL) |

**AC-2 specificity:** 기존 WT-homology 필터 효과로 FP는 15/샘플 수준 (AC-2 MVP 5는 초과). 단 이 "FP"는 UNKNOWN 시뮬레이션으로 보아야 함 — GT가 없어서 true CAND와 구분 불가.

---

## 2. AC별 현재 / MVP 목표 / Production 목표 / MVP 채택 후 예상 도달

| AC | Metric | 현재 (2026-04-16) | MVP Min | Production | **MVP 채택 후 예상** |
|---|---|---|---|---|---|
| **AC-1** | Sensitivity (GT recall) | **5/6 = 83%** | ≥80% @ core 5 host | 100% core 6 @ ≥15x | **6/6 = 100%** (tomato_A2_2 재검증 완료 + canonical triplet → soybean 2 샘플도 CAND 1+개 기대) ✅ target |
| **AC-2** | Specificity (FP/sample) | **rice 2, A2_3 1, cuc_225 7, soybean 15(=UNK 상당부분)** | ≤ 5 CAND FP | ≤ 2 CAND FP | **rice 0-1 (remote BLAST 후), A2_3 0-1, cuc_225 0-2 (remote BLAST 확정)** — Q1a Cas9 DB + Q4a L1-5 cd-hit dedup으로 중복 제거 |
| **AC-3** | Annotation completeness | **CAND 중 element named 비율 ~60%** (rice Chr3:29M, cuc_225 new 7건 UNKNOWN 다수) | ≥70% | ≥90% | **~75-85%** (SpCas9+sgRNA DB, Cas9 샘플 annotation 복구. YUC6 유보라 soybean은 여전히 UNKNOWN 위험) |
| **AC-4** | Wall time | **rice 1h15m ✅, A2_3 29m ✅, cuc 9-14h ⚠️, soybean 24h+ ❌** | ≤48h | ≤24h (<1Gbp) / ≤48h (≥1Gbp) | **rice ~1h, A2_3 ~30m, cuc ~10h, soybean array+round축소 ~3-4h** — Q2 array + Q4a L1-3 assembly round 3 |
| **AC-5** | Operator UX | config.yaml + 1 sbatch | 1 sample + 1 sbatch | zero-config | **MVP 유지** (Q1b canonical triplet config 확장만, operator UX 크게 변하지 않음) |
| **AC-6** | Audit trail | **report에 commit hash ✗, input hash ✗, DB version ✗, sw version ✗** | 3/4 | 4/4 + SHA-256 lock | **4/4** (won_yim Q5 R-1~R-4 MVP 필수 채택) ✅ |
| **AC-7** | Coverage robustness | **측정 없음** (coverage matrix 미실행) | 15x만 OK | 15x+10x OK (5x/3x degrade) | **측정 예정** (현 MVP에 포함 안됨 — Phase 2) ⚠️ 미달성 |
| **AC-8** | Failure transparency | **현재 TIMEOUT 시 verdict 없음** (UGT72E3 52/105 보고 누락) | structured error | auto retry | **compute_verdict 분리로 partial verdict + LOD 로깅** (Q3b + R-7) ✅ |

### MVP vs Production 판정

- **현재 상태:** MVP 부분 통과 (AC-1/AC-5는 min 달성, AC-2/AC-4/AC-6 미달, AC-7 미측정)
- **MVP 채택 후 예상:** AC-1~AC-6 모두 MVP threshold 달성. **AC-7만 Phase 2로 연기** (coverage sensitivity batch 지연).
- **Production 도달률:** MVP 채택 후 약 **6/8 AC = 75%**. AC-7 coverage matrix 및 AC-2 production threshold(≤2 FP)는 v1.1 대상.

---

## 3. 병목 분해 (Round 1 §1 재사용, 최소 업데이트)

| Sample | s04 wall | s05 wall | s05 per-site | s05 BLAST 비중 | Phase 3 비중 | MaxRSS | 현 이슈 |
|---|---|---|---|---|---|---|---|
| rice_G281 | ~1h (cached) | 1:15 | 3.4min × 21 | 11-18% | 82-89% | 7.7 GB / 32 GB | OK |
| tomato_A2_3 | ~30m (cached) | 29m | 4min × 2 | 11-18% | 82-89% | 11.3 GB / 48 GB | OK |
| cucumber_line224 | ~6h | 9:21 | 31min × 18 | 11-18% | 82-89% | 87 GB / 96 GB ⚠️ | **OOM margin 9%** |
| cucumber_line212 | ~6h | 13:29 | 26min × 31 | 11-18% | 82-89% | 80 GB / 96 GB | OK |
| cucumber_line225 | ~6h | 13:42 | 21min × 38 | 11-18% | 82-89% | 60 GB / 96 GB | OK |
| soybean_AtYUCCA6 | ~7h | RUNNING 21h | 16min × 100 | 11-18% | 82-89% | 17 GB / 48 GB | **진행 중** |
| soybean_UGT72E3 | ~7h | **24h TIMEOUT** | **28min × 105** | 11-18% | 82-89% | 56 GB / 96 GB | **48h 필요 (105×28min)** ❌ |

**핵심 병목:** s05 **Phase 3 per-site k-mer 확장 + Pilon (82-89%)**. BLAST refactor로는 -5~7%만 개선. Array 병렬화만이 -92% 돌파구.

---

## 4. Gap closure 기여도 매트릭스 — MVP 채택 항목 × AC

| MVP 항목 | AC-1 sens | AC-2 spec | AC-3 annot | AC-4 wall | AC-6 audit | AC-8 fail | Rationale |
|---|---|---|---|---|---|---|---|
| **Q1a SpCas9 + sgRNA DB 추가** | +10% (Cas9 샘플 annotation) | +0% | **+15%** | +0-1% | +0% | +0% | Cas9 샘플 A2_2/A2_3 element 증거 확보 |
| **Q1b canonical-triplet rule** | +10-15% (soybean 2 샘플 CAND↑) | +0% | +5% | +0% | +5% (audit trail) | +0% | soybean 0-CAND 문제 잠재 해결 |
| **Q2 per-site SLURM array** | +0% | +0% | +0% | **−92%** (48h → 3-4h) | +0% | +0% | UGT72E3 TIMEOUT 해결 경로 |
| **Q3b compute_verdict() 분리** | +0% | +2% (verdict 투명성) | +0% | +0% | +10% | **+40%** | partial verdict + LOD 기록 |
| **Q3c run_blastn() wrapper** | +0% | +2% (BUG-2/15 방지) | +0% | −5~7% | +5% | +5% | atomic write, DB 일관성 |
| **Q3d extra-DB source 태깅** | +0% | +0% | +5% | +0% | +5% | +0% | Q1b 전제 조건 |
| **Q4a L1-3 assembly round 3** | +0-2% (측정 선행) | +0% | +0% | **−20~30%** | +0% | +0% | per-site 14→5min (haibao §3) |
| **Q4a L1-4 BLAST in-memory 공유** | +0% | +0% | +0% | −5~7% | +0% | +0% | Q3c와 중첩 |
| **Q4a L1-5 cd-hit DB dedup** | +2% | +3% (BUG-3 류 방지) | +0% | +0% | **+10%** | +0% | DB reproducibility |
| **Q5 R-1~R-4 regulatory hashes** | +0% | +0% | +0% | +0% | **+60%** | +5% | input/commit/DB/sw hash |
| **Q5 R-7 false-negative LOD** | +0% | +0% | +0% | +0% | +5% | **+20%** | "detected limit" 명시 |

**AC별 총 예상 개선:**
- AC-1: +22~27% → **83% (현재) → 100%+ (capped)** ✅
- AC-2: +7% → **FP budget 개선, remote BLAST 후 확정**
- AC-3: +25% → **60% → 75-85%** (MVP ≥70% 달성)
- AC-4: **−92% (UGT72E3 48h → 3-4h)** + s05 BLAST -10% ✅
- AC-6: +95% → 거의 0%(현재) → 100% (R-1~R-4 + Q3c atomic + Q4a cd-hit)
- AC-8: +75% → **partial verdict + LOD + exit code structure**

---

## 5. 회귀 리스크 매트릭스 — MVP 변경사항 × 기존 6 GT

| 변경 항목 | rice Chr3:16.4M | A2_3 ch01:91M | A2_2 ch08:65M (미검증) | cuc_212 Chr6 | cuc_224 Chr2 | cuc_225 Chr2 (2-copy) | 회귀 가능성 |
|---|---|---|---|---|---|---|---|
| Q1a Cas9/sgRNA DB 추가 | 🟢 영향 없음 | 🟡 ch01 annotation 변화 가능 | 🟡 미검증 | 🟢 | 🟢 | 🟢 | **LOW** (read-only annotation 확장) |
| Q1b canonical-triplet rule | 🟢 | 🟢 | 🟡 | 🟢 | 🟢 | 🟡 (8 CAND 평가 변화) | **LOW-MED** (verdict 승격만, 강등 없음) |
| Q2 per-site array | 🟡 race condition (unmapped cache) | 🟡 | 🟡 | 🟡 | 🟡 | 🟡 | **MED** (전체 6 GT 모두에 영향. 철저한 regression TDD 필요) |
| Q3b compute_verdict 분리 | 🟡 로직 경계 변경 | 🟡 | 🟡 | 🟡 | 🟡 | 🟡 | **MED** (7 TDD 시나리오 + 샘플 재실행으로 완화) |
| Q3c run_blastn wrapper | 🟢 | 🟢 | 🟢 | 🟢 | 🟢 | 🟢 | **LOW** (기계적 치환, atomic write 추가만) |
| Q4a L1-3 round 3 | 🟡 3-round 내 수렴 미확인 | 🟡 | 🟡 | 🟡 | 🟡 | 🟡 | **MED** (haibao §5 도전질문 1 — 측정 선행 필요) |
| Q4a L1-5 cd-hit dedup | 🟡 DB ID 변화로 element tag 변화 | 🟡 | 🟡 | 🟡 | 🟡 | 🟡 | **MED** (report element_name 변경 가능) |

**총합:** 개별 변경은 LOW-MED, 합쳐서 MVP 배포 전 **6 GT 전수 재검증 + regression TSV diff 필수**. won_yim MVP 판정 (Round 2 §Q3a에서 "refactor 도중 회귀가 GT 회복보다 클 수 있다")는 여기에서 입증됨.

### 회귀 방어 체계 (MVP 선결 조건)
1. 각 commit 직전 `pytest tests/` 10/10 PASS
2. 각 commit 직후 rice_G281 + A2_3 smoke test (~1.5h, 2 샘플만) → CAND 회복 확인
3. MVP 병합 직전 **6 GT 전수 재검증** (8-20h, array 병렬 시 단축 가능)
4. 기존 snapshot (`docs/superpowers/runs/2026-04-15-*-postfix.txt`)과 verdict diff 0 또는 설명 가능한 증가만 허용

---

## 6. Remaining Gaps to Production (MVP 이후)

### Phase 2 (v1.1, 1-3개월)
| 항목 | 출처 | 이유 | 예상 공수 |
|---|---|---|---|
| **AC-7 coverage matrix 실행** | won_yim R1 §2.3 | 15x/10x/5x/3x 9 variant batch | 2-3일 (array job) |
| **s05 8-10 모듈 분할** | bio_king R1 §2 | refactor 회귀 리스크 → production 후 | 2-3주 |
| **haibao Level 1 Phase 2: cluster_window adaptive + LAST switch** | haibao R1 Level 1 L1-1/L1-2 | GT regression 위험 있음 | 1주 |
| **PE discordancy signal** | haibao Level 2 부분 | BUG-7 재발 방지 | 1주 |
| **R-5 PDF-ready report** | won_yim R2 Q5 | regulatory UX 개선 | 1주 |
| **R-6 KCGP nomenclature 맵핑** | won_yim R2 Q5 | gmo_expert 협업 | 1주 |
| **R-8 chain of custody 표준화** | won_yim R2 Q5 | SLURM JobID + operator 집계 | 3일 |
| **Canola/감자 host reference 추가** | won_yim R1 §2.1 | 수입 우선순위 2차 | 1-2주 |
| **WT-free fallback algorithm** | haibao + won_yim R2 §Q4d 3 | soybean/corn exploratory | 2주 |
| **quality-degraded synthetic read 검증** | won_yim R2 §3 gmo_expert 3 | 격리소 input 품질 대응 | 1주 |

### v2.0 / Grant proposal (3-6개월)
| 항목 | 출처 | 이유 |
|---|---|---|
| **haibao Level 2 dual-anchor rewrite** | haibao Level 2 | 논문 투고 경쟁력 |
| **haibao Level 3 GRIDSS-style graph + macrosynteny** | haibao Level 3 | 알고리즘 재설계 |
| **ONT long-read readiness** | haibao Level 3 | USDA/APHIS 트렌드 대응 |
| **Nextflow/Snakemake 마이그레이션** | compute_eng R1 §3.C | resume capability, profile |
| **Real-time regulatory portal integration** | won_yim | 격리소 실시간 보고 |

---

## 7. 결론 — 정량 요약

- **현재 상태:** MVP 6/8 AC 부분 달성 (AC-1 83%, AC-5 / AC-3 min / AC-6 0%, AC-2/AC-4/AC-7 미달).
- **MVP 채택 후 예상:** **AC-1~AC-6 6개 모두 MVP threshold 달성** (≈ 75% production 도달).
- **남는 gap:** AC-7 coverage matrix (v1.1), AC-2 production threshold (≤2 FP) (v1.1), haibao Level 2/3 (v2.0).
- **최대 ROI 항목:** Q2 array 병렬화 (AC-4 −92%), Q5 regulatory hashes (AC-6 +95%), Q3b compute_verdict 분리 (AC-8 +75%).
- **최대 리스크 항목:** Q2 array race condition, Q3b verdict 로직 경계, Q4a round 3 (측정 선행 필수).
- **판정:** MVP 채택안이 **regulatory production deploy 가능 최소 기능**에 도달. v1.1/v2.0 로드맵은 논문·경쟁력 목적으로 분리.

---

**Round 3 compute_eng deliverable 완료. won_yim team-consensus.md synthesis에 참조 제공.**
