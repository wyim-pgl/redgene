# Round 3 — won_yim (PI) — 최종 Synthesis 및 권위 판정

**작성자:** won_yim (PI / end-user, 한국 격리소 GMO/LMO 검출 assay)
**Date:** 2026-04-16
**기반:** Round 1 (4명 + 본인) + Round 2 (4명 cross-talk + 본인 Q1-Q5 판정) + team-lead Round 2 brief
**이 문서의 권위:** Round 4에서 team-lead가 5개 deliverable 문서 통합 시 **최종 판정 근거**.

---

## 1. 최종 확정 — Round 2 4건 합의 (채택 확정)

| # | 합의 사항 | 확정 판정 | 근거 |
|---|----------|-----------|------|
| C-1 | **Canonical-triplet rule → `config.yaml` 외부화 (`VerdictRules` dataclass)** | **MVP 확정 채택** | bio_king + gmo_expert + won_yim 3자 합의. bio_king R2 Q2에서 "DEFAULT_TRIPLETS 하드코딩 fallback + config override"로 TDD 가능성 보존 확인. regulatory flexibility + audit 우위. |
| C-2 | **DB 정규화(cd-hit + source tag) + in-memory `BlastHits` sharing** (P0 #3 확장) | **MVP 확정 채택** | bio_king R2 Q1 수정 답변에서 "내 래퍼만으로 부족, DB 정규화까지 가야 BUG-2/3/4 class 근본 차단" 인정. haibao P2 dedup과 중첩 → 1회 작업으로 해결. |
| C-3 | **2-tier priority (element_db-family > univec)** | **확정 채택** | gmo_expert R2 Q3에서 3-tier는 BUG-3 재발 risk라고 스스로 2-tier 제안. common_payload / element_db / CRL을 같은 tier 2로 묶고 tie는 bitscore로 해결. 구현 단순성 + 도메인 정합. |
| C-4 | **haibao Level 2/3 → v2 grant 이관** | **확정 기각 (MVP out)** | haibao R2 Q1 "Level 2 rewrite는 MVP scope out 수용" 명시. PE discordancy / dual-anchor / MCscan 상위 워크플로는 v2.0 grant 제안서로 분리. |

**bio_king R2에서 추가된 핵심 제약:**
- 구현 순서는 **"모듈 분할 최소부 → verdict 순수화 → DB 정규화 → in-memory sharing"** 고정.
- CLI (`--phase 1,1.5` / `--site-id <id>` / `--phase 4`)는 **분할 완료 후** 추가 (CLI-first는 band-aid가 되어 나중에 뜯어야 함).
- 단, UGT72E3 즉시 재제출이 필요하면 `run_pipeline.py` 레벨 fan-out(별도 worktree N-way sbatch) **임시안**은 허용.

---

## 2. compute_eng 실측 ROI 기반 최종 채택 순서 (수용 확정)

compute_eng R2 §종합 권고의 ROI 기준을 **권위 순서**로 확정. team-lead 요약 그대로 채택:

| 우선순위 | 항목 | 실측 효과 | 최종 판정 |
|----------|------|-----------|-----------|
| **1 (MVP 필수)** | **s05 per-site SLURM array** | UGT72E3 48h → 4h (**-92%**) | **이번 주 구현** |
| 2 (MVP 조건부) | s04 BWA → minimap2 전환 | -40% 가능하나 BUG-7 회귀 리스크 | **PoC 먼저 (아래 §5)** |
| 3 (MVP 채택, 성능 보조) | BLAST 래퍼 refactor + in-memory sharing | -5 ~ -7% (code quality 주 가치) | **이번 주 구현** (C-2와 합체) |
| 4 (MVP 채택, 안전) | DB 확장 (P0 6 + P1 8 seqs) | **+1 ~ +1.3%만 증가** | **이번 주 구현** |

**AC-4 달성 경로:** array(-92%) + assembly round 8→3 (-15 ~ -20% 추가, haibao P1) = **UGT72E3 48h → ~3h** 보수적 기대. MVP AC-4 (≤ 48h) 강력 통과, production AC-4 (≤ 24h) 도달.

---

## 3. haibao P1/P2/P6 최종 판정

haibao R2는 Level 1에서 5개 patch 중 **3개만 사수** (P1, P2, P6). 나머지 2개(P4 cluster_window adaptive, P5 LAST 옵션)는 Phase 2 이관 자진 동의.

| patch | 판정 | Rationale |
|-------|------|-----------|
| **P1 — assembly `max_rounds=3`** | **MVP 채택** | 1줄 수정 (s05:2372). haibao R2 Q2에서 "1개만이라면 P1"로 스스로 지목. **선결 조건: `s05_stats.txt` aggregate로 round-3 수렴율 사전 측정 30분 작업** (resume.md "most converge by 6"이 관찰치일 뿐 optimal 아님). 측정에서 round-3 기여 ≤ 10%면 8→3, 10~30%면 8→5, ≥ 30%면 유지. |
| **P2 — `transgene_db` build-time cd-hit-est dedup** | **MVP 채택 (C-2와 통합)** | `element_db/build_common_payload.sh`에 `cd-hit-est -c 0.95` 1 step 추가. Makefile 1줄. BUG-3 근본 해결. C-2 DB 정규화 작업과 같은 PR로 묶음. |
| **P6 — element vs host 100%-ortholog pre-mask BED** | **MVP 채택** | **team-lead 요청 1차 질문.** gmo_expert R2 Q2와 haibao R2 Q4 양쪽이 **결합 권고**. 근거: <br> • AC-2 (cucumber line225 8 CAND 과잉) 해결에 직접 기여 <br> • 구현 비용 ≤ 1일 (BLAST 1회 + bedtools + 20줄 Python, haibao R2 Q4 Workflow 1 스케치 확정) <br> • 회귀 위험 **매우 낮음** (mask된 site는 drop이 아니라 low-confidence tag로만 downgrade, gmo_expert R2 §Regulatory 안전성 판정의 "FALSE_NEGATIVE_MASKED 태그 + rationale.tsv" 구조로 regulatory audit 요건 충족) <br> • **제약 1:** 100%-ortholog (corn × P-Ubi1, rice × P-Act1, soybean × AtYUCCA6 cultivar drift 포함)만 **pre-mask**. EPRV / pararetrovirus 잔재는 **사후 Filter D**에서 처리 (gmo_expert R2 §Regulatory matrix) <br> • **제약 2:** BED + `host_masked_rationale.tsv` 2-파일 쌍을 **모든 sample report에 첨부** (regulatory audit trail). |
| P3 — BLAST 호출 통합 (in-memory 3-filter) | **MVP 채택 (C-2에 흡수)** | compute_eng 실측 -5 ~ -7%. 독립 patch로는 ROI 약함이지만 C-2 in-memory `BlastHits` sharing 작업과 동일 코드 경로 → 별개 작업 아님. |
| P4 — cluster_window adaptive | **Phase 2 연기** | haibao R2 Q2에서 회귀 위험 "매우 높음" (BUG-7 revert class). MVP 후 재검토. |
| P5 — LAST 옵션 | **Phase 2 연기** | gmo_expert 131 element 대한 LAST regression 미검증. v2에서 검토. |
| k-mer k=15 → 21 변경 | **실험 측정 후 결정 (보류)** | haibao R2 Q3 본인 자인: "실험 검증 없이 default 변경 금지". Phase 2에서 실측 후 재판정. |

---

## 4. Q1-Q5 최종 판정 확정 (Round 2 판정 재확인 + 업데이트)

### Q1 (AC-1 sensitivity)
- **Q1a. SpCas9 CDS + sgRNA scaffold DB 추가 → MVP 확정 채택.**
  compute_eng R2 Q3 실측: +15 seqs로 BLAST 비용 +1%만 증가, 20% 한계 내 완전 안전. gmo_expert R2 Q1 MUST-HAVE subset 15개 (P0 6 + P1 8 + P2 1) 중 **P0 6개 + P1 8개 = 14개 이번 주 추가**. CRL subset(P2)은 gmo_expert가 raw TSV 공급 후 Phase 2.
  - **조건 추가 (gmo_expert R2):** SpCas9는 Phase 1.5 `--cas9-info-only` flag로 CANDIDATE 승격 트리거에서 제외 → tomato_Cas9_A2_* 샘플의 site count 폭증 방지.
  - **조건 추가 (compute_eng R2 부록 B):** 새 entries에 `_filter_host_endogenous` 확장 적용 (host cross-reaction pre-filter, per-sample +1분).
- **Q1b. Canonical-triplet rule (config-driven) → MVP 확정 채택** (C-1과 일치). bio_king의 `VerdictRules` dataclass + `DEFAULT_TRIPLETS` fallback 구현 방식 확정.
- **Q1c. AtYUCCA6 CDS 복원 → MVP 채택으로 승격.**
  - **Round 2 유보 → Round 3 채택 사유:**
    gmo_expert R2 Q1에서 AtYUCCA6/UGT72E3를 P1 must-have 14개에 포함시켰고, compute_eng가 +1% cost만 인정함. **Q1a(15개 추가)와 별도 작업이 아니라 같은 batch**. soybean 0-CAND 해결의 **direct path**로 채택.
  - **조건:** gmo_expert R2 §Matrix의 "soybean × AtYUCCA6 host ortholog (75% 유사도, systemic FP source)"는 **P6 pre-mask BED에 반드시 포함**. 이 두 작업이 세트로 묶여야 복원이 안전.

### Q2 (AC-4 turnaround)
**per-site SLURM array 병렬화 → MVP 필수 확정 지정** (§2 우선순위 1번). bio_king R2 Q3이 "분할 먼저, CLI는 분할 완료 후" 원칙을 제약으로 추가. **won_yim 확정 판정:**
- bio_king 원칙 수용: `--phase 1,1.5` / `--site-id <id>` / `--phase 4` CLI entrypoint는 **최소 모듈 분할 (site_discovery / classify / per_site / annotate 4-way boundary) 완료 후** 추가.
- **그러나** UGT72E3 소이빈 즉시 재제출이 Round 3 직후 필요 → bio_king이 제안한 **`run_pipeline.py` 레벨 fan-out 임시안** (site list split + N-way sbatch via 독립 worktree) 병행 허용. **임시안은 v1.0 production release 후 제거 조건.**

### Q3 (refactor scope)
**won_yim Round 2 판정 유지.** 단, C-2 (DB 정규화 + in-memory sharing) 추가로 MVP refactor 범위가 약간 확장:

| refactor 항목 | 최종 판정 | Rationale |
|---|-----------|-----------|
| s05 8-10 전체 모듈 분할 | **Phase 2 연기** | Round 2 판정 유지. |
| **최소 4-way 모듈 분할** (site_discovery / classify / per_site / annotate+report) | **MVP 채택** | Q2 per-site array CLI의 전제 조건 (bio_king R2 Q3). 8-10 모듈 분할의 subset, 회귀 위험 제한적. |
| `compute_verdict()` pure-function 분리 | **MVP 확정 채택** | AC-6 audit trail + bio_king R2 Q2 `VerdictRules` dataclass 주입 경로. |
| `run_blastn()` 공통 래퍼 | **MVP 확정 채택** | C-2에 흡수. |
| **DB 정규화 + in-memory `BlastHits` sharing** (bio_king R2 Q1 수정) | **MVP 추가 채택** | C-2. 단독 P3 BLAST 통합 흡수. |
| extra-DB source tag 3-way 분리 | **MVP 채택 → 4-way로 확장** | "element_db / common_payload / sample_contig / univec" + canonical_triplet virtual tag. gmo_expert R2 Q3의 2-tier priority와 정합 (tier 2 내부에서 source별 bitscore tie). |

### Q4 (haibao Level 2/3)
- **Level 2 rewrite (dual-anchor + PE discordancy + multi-signal) → v2.0 grant 이관 확정** (C-4).
  - haibao R2 Q1 "Level 2 기각 OK" 명시적 동의. PE discordancy 추가는 v1.1 로드맵 후보 (MVP 통과 이후).
- **Level 3 full redesign (GRIDSS + macrosynteny + long-read) → v2.0 grant proposal 자료로 분리 확정** (C-4). MVP 완전 out.
- **MCscan Workflow 2 (synteny-based EPRV mask) → v2 이관.** haibao R2 Q4 Workflow 2 "수동 큐레이션 필요, v2 infrastructure" 동의.

### Q5 (regulatory MUST 8개)

team-lead 요청 "3-4개 고정"에 따라 **Round 2 판정 5개 → 4개로 축소** (이번 주 구현 현실성 반영):

| # | 요건 | 최종 판정 | Rationale |
|---|------|-----------|-----------|
| R-1 | **Input SHA-256 hash 기록** | **MVP 필수 (이번 주 구현)** | `run_pipeline.py` 시작부 sha256sum + 로그. ≤ 0.5일. AC-6 최저선. |
| R-2 | **Pipeline commit hash lock + dirty/clean 체크** | **MVP 필수 (이번 주 구현)** | `git rev-parse HEAD` + `git status --porcelain` empty check. ≤ 0.2일. AC-6 최저선. |
| R-3 | **BLAST DB md5 + build date manifest** | **MVP 필수 (이번 주 구현)** | Q4a P2 cd-hit dedup 구현과 같은 PR에 `db/element_db_manifest.tsv` (name, md5, build_date, seq_count) 출력. 추가 비용 0일. |
| R-4 | **Software version manifest** | **MVP 필수 (이번 주 구현)** | `bwa --version`, `samtools --version`, `blastn -version`, `spades.py --version`, `minimap2 --version`, `cd-hit-est --version` 로그. ≤ 0.2일. |
| R-7 | **False-negative 명시 (coverage/MAPQ/LOD)** | **Phase 2 이관 (Round 2 판정 변경)** | Q3 `compute_verdict` 분리와 같은 PR 원했으나, **이번 주 구현 범위 초과**. v1.0에서 report에 "not detected"시 detection limit 문구 표준 삽입만 하고, 정량 LOD (coverage × MAPQ × BLAST evalue 통합) 은 Phase 2. |
| R-5 | PDF-ready report | Phase 2 | Round 2 판정 유지. |
| R-6 | KCGP element nomenclature | Phase 2 | Round 2 판정 유지. gmo_expert 추가 작업 필요. |
| R-8 | Chain of custody log | Phase 2 | Round 2 판정 유지. |

**MVP 필수 = 4개 (R-1, R-2, R-3, R-4)**. **R-7은 Phase 2로 조정** (이번 주 현실 범위). 이 4개는 통합 spec이 간단 — `run_pipeline.py` 시작부 `_write_audit_header()` 단일 함수에서 전부 출력.

---

## 5. 잔존 충돌 1건 최종 판정 — s04 minimap2 전환

**충돌 요약:**
- **haibao Round 1/2:** s04 BWA → minimap2 `-ax sr` 전환, -40% wall time, hybrid short+long preparation.
- **compute_eng R2 Q2 실측:** Phase 1 site scan을 minimap2로 교체 시 +2-3x 악화 (재매핑 비용), MAPQ 분포 변화로 BUG-7 회귀 가능. 단, **s04 host mapping 단계 전환은 긍정적**으로 명시 (R2 §종합 권고 2번).

**혼동 정리:**
- haibao Level 2 "Phase 1 site discovery minimap2 dual-anchor" = compute_eng가 반대한 것.
- haibao Level 1 언급한 "s04 host mapping BWA → minimap2" = compute_eng가 찬성한 것.

**team-lead 요청 1차 질문 (rice_G281 PoC만 승인 vs 완전 기각):**

### 최종 판정 — **rice_G281 PoC만 MVP 승인 (Conditional)**

**Rationale:**
1. compute_eng 부록 C 제안 실험 설계는 정밀(-1 rice 샘플, 기존 s01 재활용, 총 wall time ~4-6h, 승격 기준 명시: CAND 유지 + Phase 1 site ±20% + wall time -30% 이상).
2. rice_G281은 AC-1 GT anchor이며 verdict 회귀 판정 기준이 **객관적** (Chr3:16,439,674 CAND 유지 여부).
3. **PoC 실패 시 자동 기각** — haibao가 스스로 수용할 수 있는 기각 근거 (측정 기반, 논거 아님).
4. PoC 성공 시에도 **MVP 배포는 rice/tomato 두 host까지만**. cucumber / soybean / corn host에서의 추가 PoC는 **Phase 2로 이관** (host별 MAPQ 분포 이질성이 validation 범위 초과).

**PoC 승인 조건 (won_yim 명시):**
- **승격 기준 (all must pass):**
  1. rice_G281 Chr3:16,439,674 → CANDIDATE 유지 (verdict 동일)
  2. Phase 1 site 수가 BWA 대비 ±20% 이내
  3. s04 wall time -30% 이상 (BWA 대비)
  4. MAPQ 분포 비교: soft-clip reads의 MAPQ < 20 비율이 BWA 대비 증가하지 않음 (BUG-7 class 회귀 방지)
- **실패 시 처리:** s04 BWA 유지. haibao의 minimap2 전환 주장은 v1.1 로드맵으로 자동 이관.
- **실행 소유자:** compute_eng (측정) + haibao (해석) 공동. Round 4 이후 PoC 수행.
- **시한:** v1.0 production release 전.

---

## 6. 이번 주 구현 가능 범위 (team-lead 요청 — 최종 MVP 구현 목록)

**"이번 주 (2026-04-16 ~ 2026-04-22, 5 영업일)"** 구현 가능 범위를 공수 기반으로 확정.

### 6.1 이번 주 구현 (MVP core — 반드시 완료)

| # | 항목 | 담당 | 공수 | 의존성 | AC 기여 |
|---|------|------|------|--------|---------|
| W1 | **AC-6 audit header 4종** (R-1~R-4: input sha256, commit hash, DB md5, software versions) | compute_eng | 0.5일 | 없음 | AC-6 |
| W2 | **SpCas9 + sgRNA scaffold DB 추가** (gmo_expert P0 2개 먼저) | gmo_expert | 0.5일 | W1 (DB manifest) | AC-1, AC-3 |
| W3 | **canonical-triplet rule config 외부화** (`VerdictRules` dataclass + DEFAULT_TRIPLETS fallback) | bio_king | 1-2일 | 없음 | AC-1 |
| W4 | **per-site SLURM array 임시안** (`run_pipeline.py` fan-out, 모듈 분할 전) | compute_eng + bio_king | 1-2일 | 없음 | AC-4 |
| W5 | **assembly round 8→3** (P1, haibao, 선결 측정 30분 포함) | haibao | 0.5일 | `s05_stats.txt` aggregate | AC-4 |
| W6 | **transgene_db cd-hit-est dedup + source tag** (P2, C-2 일부) | bio_king | 0.5일 | W2 DB 확장 직후 | AC-1, AC-2 |
| W7 | **element 100%-ortholog pre-mask BED** (P6, haibao Workflow 1) | haibao + gmo_expert | 1일 | W2 DB 확장 | AC-2 |
| W8 | **UGT72E3 재제출 + AtYUCCA6 검증** (W4 + W5 + W7 결과 확인) | compute_eng | 0.5일 (submit) + ~8h wall | W4, W5, W7 | AC-1 (0-CAND 해결), AC-4 |

**주간 총 공수:** ~6-8 person-day. 4명 담당자 parallelize 가능 → 5 영업일 내 수용 가능.

**AC 목표 달성 예상:**
- AC-1: 5/6 → **6/6 또는 7/8** (Cas9 positive 검증, soybean 0-CAND 해결 기대)
- AC-2: rice/tomato/cucumber에서 FP ≤ 5 유지. P6 pre-mask로 cucumber 8 CAND → 3-4 예상.
- AC-4: UGT72E3 48h TIMEOUT → **≤ 4h**. MVP ≤ 48h 강력 통과.
- AC-6: 4/4 manifest 구비. Regulatory 제출 최저선 도달.

### 6.2 이번 주 OUT (MVP 이월 / Phase 2 / v2)

| 항목 | 분류 | Rationale |
|------|------|-----------|
| s05 8-10 전체 모듈 분할 | Phase 2 | 최소 4-way 분할만 이번 주 수용 가능 (W4와 함께), 전체 분할은 v1.1. |
| `compute_verdict()` pure 함수 분리 (full) | v1.1 (2주차) | W3 `VerdictRules` 구조 합의 선행. 합의 후 1주 작업. |
| DB 정규화 + in-memory `BlastHits` sharing (full C-2) | v1.1 (2주차) | 최소 cd-hit dedup (W6)만 이번 주, 전체 in-memory sharing은 모듈 분할 이후. |
| s04 minimap2 PoC | v1.0 이전 PoC | §5 조건부 승인, 별도 실험 branch. 이번 주 구현 OUT, PoC 결과 Round 4 이후 판정. |
| haibao P3/P4/P5 | Phase 2 | 본인 동의. |
| R-5/R-6/R-7/R-8 | Phase 2 | §Q5 판정대로. |
| haibao Level 2 (dual-anchor + PE) | v1.1 또는 v2 | PE discordancy만 v1.1, 전체는 v2 grant. |
| haibao Level 3 (GRIDSS + synteny + long-read) | v2.0 grant proposal | 완전 out. |
| k-mer k=15 → 21 | Phase 2 (측정 필요) | haibao 본인 자인. |
| `crl_amplicons.fa` 82개 통합 | Phase 2 | gmo_expert P2. |
| MCscan Workflow 2 (EPRV synteny mask) | v2 | haibao + gmo_expert 동의. |

### 6.3 이번 주 리스크 및 대응

| 리스크 | 완화 |
|--------|------|
| W4 임시안이 나중에 뜯어낼 band-aid | v1.0 release checklist에 "W4 fan-out 제거 및 CLI entrypoint 교체" 명시. |
| W5 `max_rounds=3` 측정 결과가 round-3 기여 ≥ 30%면 recall 회귀 | 측정 결과 ≥ 30%시 8→5로 타협. 측정 0.5h 포함된 공수. |
| W7 pre-mask BED에서 진짜 T-DNA junction 손실 | gmo_expert R2 "FALSE_NEGATIVE_MASKED 태그 + rationale.tsv" 구조 엄수. mask된 site는 drop 아니라 downgrade only. |
| W8 UGT72E3 재제출 시 여전히 0-CAND | fallback: AtYUCCA6 CDS 포함 (gmo_expert P1, Q1c 채택됨) + W3 canonical-triplet rule 동시 활성화 시 2중 보장. |

---

## 7. v1.0 Production Release 선언 기준 (이번 주 이후)

이번 주 W1-W8 완료 후 다음 주 (2026-04-23 ~) **v1.0 release candidate** 판정 기준:

1. **AC-1 sensitivity ≥ 100% (GT anchor 6 samples)**: rice_G281, tomato_A2_2, tomato_A2_3, cucumber_line212/224/225 전원 CANDIDATE 회수.
2. **AC-2 specificity ≤ 5 CAND FP / 샘플**: MVP 수용 기준 (≤2 production은 v1.1 목표).
3. **AC-4 turnaround ≤ 48h / 샘플**: soybean UGT72E3 포함 **모든 샘플**에서 통과.
4. **AC-6 audit trail 4/4**: report header에 input sha256 + commit hash + DB manifest + software versions.
5. **Regression test 10/10 PASS** (bio_king pytest baseline) + new tests for `compute_verdict` (최소 3 시나리오: CANDIDATE/FP/UNKNOWN).
6. **exploratory 샘플 (soybean AtYUCCA6/UGT72E3, corn ND207) 중 ≥ 1 CAND 회수** OR FP ≤ 5 유지 (AC-2 exploratory 판정).

**v1.0 통과 후:** 한국 격리소 **pilot assay 배포 가능** 수준. 실제 검역 사용 전 v1.1 (Phase 2 주요 항목) 및 외부 validation 필요.

---

## 8. teammate별 Round 4 실행 지시 (won_yim PI 권한)

| 담당 | Round 4 주 작업 | 책임 산출물 |
|------|-----------------|-------------|
| **bio_king** | W3 (`VerdictRules` config 외부화) + 최소 4-way 모듈 분할 boundary 정의 | PR 1개: `s05/verdict.py` + `s05/config_loader.py`, pytest 3+ test |
| **gmo_expert** | W2 (Cas9/sgRNA DB + P0 6 + P1 8 seqs) + W7 BED rationale 큐레이션 + `host_masked_rationale.tsv` 초안 | DB PR + rationale TSV |
| **compute_eng** | W1 audit header + W4 fan-out + W8 재제출 + rice_G281 s04 minimap2 PoC | run script 묶음 + audit header code + PoC 결과 |
| **haibao** | W5 assembly round 측정 + W7 pre-mask BED script (Workflow 1) + s04 PoC 해석 | measurement report + `scripts/build_element_mask_bed.sh` |

**Round 4 동기화 주기:** 일 1회 TaskList + 주 1회 broadcast. 주 후반(금요일) integration day — 전체 샘플 재검증 batch 실행.

---

## 9. 최종 한줄 요약 (team-lead synthesis용)

> **MVP v1.0은 "AC-1 6/6 GT recall + AC-2 ≤ 5 FP + AC-4 ≤ 48h + AC-6 4-field audit header" 4개 조건 달성으로 정의한다. 이번 주 W1-W8 (8 task, 6-8 person-day, 4명 parallelize)으로 도달 가능하며, s05 전체 refactor·Level 2 rewrite·minimap2 Phase1 교체·k-mer 변경·MCscan 상위 워크플로는 모두 v1.1 이후로 이관한다. s04 minimap2는 rice_G281 PoC 조건부 승인만 허용. regulatory 제출 최저선은 input sha256 + commit hash + DB md5 + software versions 4종으로 고정한다.**

---

**Round 3 완료. Round 4에서 team-lead가 5개 deliverable 통합 시 본 문서를 권위 판정 근거로 사용. won_yim 대기 상태 전환.**
