# Round 2 — won_yim (PI/end-user) — 1차 Scope 판정

**작성자:** won_yim · **Team:** redgene-review-2026-04-16 · **Date:** 2026-04-16
**기반:** 4명 Round 1 노트 (bio_king / gmo_expert / compute_eng / haibao) + 본인 Round 1 §5
**판정 기준:** Round 1 §1의 8개 Acceptance Criteria (AC-1~AC-8)에 대한 직접·간접 기여도.

---

## 요약 (결정 일람)

| Q | 핵심 제안 | 판정 | MVP / Phase 2 | AC 기여 |
|---|---------|------|----------|---------|
| Q1 | SpCas9 + sgRNA scaffold DB 추가 | **MVP 채택** | MVP | AC-1, AC-3, 규제 MUST #6 |
| Q1 | Soybean canonical-triplet rule (bar+P-35S+T-ocs) | **MVP 채택 (config-driven)** | MVP | AC-1, AC-3 |
| Q1 | AtYUCCA6 CDS 복원 | **유보 (Round 3)** | — | — |
| Q2 | per-site SLURM array 병렬화 | **MVP 필수 지정** | MVP | AC-4 |
| Q3 | s05 8-10 모듈 분할 | **Phase 2로 연기 (부분 예외)** | Phase 2 | — |
| Q3 | `compute_verdict()` pure-function 분리 | **MVP 예외 채택** | MVP | AC-6, AC-8 |
| Q3 | BLAST 공통 래퍼 `run_blastn()` | **MVP 채택** | MVP | AC-4, AC-6 |
| Q4 | haibao Level 1 patch 5개 | **부분 채택 (3/5)** | MVP 2개 / Phase 2 3개 | AC-1, AC-4 |
| Q4 | haibao Level 2 module rewrite (dual-anchor) | **Phase 2 연기** | Phase 2 | — |
| Q4 | haibao Level 3 full redesign (GRIDSS) | **기각 (v2.0 grant)** | out | — |
| Q5 | regulatory MUST 8개 | **3/5/0 분류** | 표 참조 | AC-6 |

---

## Q1 — AC-1 Sensitivity 관련 제안 (gmo_expert)

### Q1a. SpCas9 CDS + sgRNA scaffold DB 추가 — **MVP 채택 (P0)**

**Rationale:**
- gmo_expert 노트 §2, §6.3에서 지적한 "Cas9 검출 기능 부재"는 **regulatory MUST**.
  tomato_Cas9_A2_1/A2_2/A2_3 세 샘플이 `config.yaml: cas9_present: true`만으로 판정되고
  있고, 서열 기반 증거가 없다. 이는 Round 1 §4 regulatory #6 (KCGP nomenclature) 및
  #7 (false-negative 명시)에 직접 위반.
- AC-1 (sensitivity) 지표가 **"tomato_A2_2 ch08:65,107,378 Cas9-positive"** 를 포함한다면,
  sgRNA scaffold 없이는 recall 계산 자체가 불완전.
- gmo_expert 부록 B에 accession 후보 (NC_017053.1 SpCas9, pX330/pRGEB32 sgRNA scaffold
  80bp) 명시 → 구현 비용 ≤ 1일 (`build_common_payload.sh`에 efetch 2줄 추가).

**조건:**
- **gmo_expert가 우려한** "tomato_Cas9_A2_* 샘플에서 Phase 1.5 site count 폭증" 위험
  → 확장 직후 **세 샘플 재검증 필수** (배치 작업 ≤ 3h). 회귀 없을 경우에만 main 병합.
- compute_eng 노트 §7 질문 1 (DB 확장 시 BLAST 시간 선형 증가)을 고려, Cas9/sgRNA 2 seq만
  추가 → Phase 1.5 BLAST 시간 영향 ≤ 1%. 수용.

**예상 공수:** 1일 (DB update + efetch + 샘플 재검증).

### Q1b. Soybean canonical-triplet rule (bar + P-35S + T-ocs → auto-CANDIDATE) — **MVP 채택 (config-driven)**

**Rationale:**
- gmo_expert §5 결론에서 "bar + P-CaMV35S + T-ocs 같은 canonical triplet이 single contig에
  있으면 LMO-positive로 확정하기에 domain적으로 충분"이라 명시. 이는 **식약처 element-based
  screening 공식 해석과 일치** (Round 1 §1 AC-1/AC-3 목표 달성 경로).
- AC-1 현재 점수 **5/6 GT** (83%) → soybean이 exploratory이나 canonical-triplet 규칙이
  성립하면 bar+P-35S+T-nos가 포함된 AtYUCCA6/UGT72E3도 자동 CANDIDATE 승격 가능 → AC-1의
  **"0-CAND 문제" 우회 해결**.
- bio_king §5 smell #5 (`classify_site_tiers`가 extra_db를 "element_db"로 강제 태깅)과
  정합. "canonical_triplet" tag를 3rd source로 분리하면 tier 세분화 가능.

**Scope 제한 (gmo_expert vs bio_king 충돌 조정):**
- gmo_expert는 "per-sample white-list를 config로" 선호. bio_king은 "code에서 rule로".
  **won_yim 결정: config.yaml driven**. 이유:
  1. **규제 flexibility** — 새 marker 추가 시 코드 수정 없이 operator가 반영 가능
     (AC-5 operator UX 강화).
  2. **Auditability** — config.yaml이 보고서에 포함되면 "이 샘플은 어떤 triplet 규칙으로
     판정되었나"가 투명 (AC-6 audit trail).
- 구현 방식: `config.yaml` 샘플 entry에 `auto_cand_markers: [bar, P-CaMV35S, T-ocs]` 또는
  global default `canonical_triplets:` 섹션. s05 `compute_verdict`가 이를 읽어 "모든
  marker가 하나의 contig / insert 내 ≥90% identity로 발견되면 CANDIDATE로 승격".

**예상 공수:** 2-3일 (config schema + verdict 로직 + 3-host 재검증).

### Q1c. AtYUCCA6 CDS 복원 — **유보 (Round 3 판정)**

**Rationale:**
- gmo_expert §5 H1: "AtYUCCA6 제거 정책이 soybean 0-CAND 원인일 수 있음" — **가설 단계**.
- Q1b (canonical-triplet rule)가 통과하면 AtYUCCA6 복원 없이도 bar/P-35S/T-ocs 조합으로
  soybean CANDIDATE 확보 가능. **우선 Q1b로 해결 시도 → 실패 시에만 복원 검토.**
- 복원 시 rice/tomato/cucumber 호스트에 YUC6 ortholog (75% homology) 오탐 위험이 있어
  per-host blacklist가 필요 → 부담 증가.
- **결정:** Round 3에서 Q1b 재검증 결과 본 뒤 확정.

---

## Q2 — AC-4 Turnaround: per-site SLURM array 병렬화 (compute_eng) — **MVP 필수 지정**

**Rationale:**
- compute_eng §3.A의 실측 근거가 명확: UGT72E3 105 sites × 28min = **48h (24h TIMEOUT 기확인)**.
  array 25 concurrent × 4 waves → **~3-4h 벽시간**. AC-4 (≤ 48h, 목표 ≤ 24h) 달성에 **유일한
  cost-effective 경로**.
- CPU-hour 총량은 동일 → cluster cost 증가 0. queue 경쟁만 ↑.
- haibao §3 Level 1 제안 #3 (assembly round 8 → 3, per-site 14min → 5min)과 **독립적
  적용 가능** → 두 개 조합 시 UGT72E3 ≤ 2h.
- soybean 외에도 coverage sensitivity batch (AC-7 측정을 위한 9 subsampled 샘플) 실행에도
  array가 유리.

**조건:**
- **bio_king의 per-site CLI 계약** (compute_eng §7 질문 1, 2)과 동시 작업 필요:
  - `s05_insert_assembly.py --site-id-only <id>`, `--phase-4-only` 옵션 추가.
  - 이는 Q3의 refactor 중 **최소 분리만** 포함 (전체 8-10 모듈 분할 X, Phase 2+3 독립 CLI만 O).
- `_spades_run`, `step_dir` 공유 경로 격리 (compute_eng §3.A 언급). array task별 `workdir=step_dir/_array_<site_id>`.
- bin/rg-sbatch wrapper (compute_eng BUG-16 fix)와 함께 병합 시 SLURM env 오염 방지.

**예상 공수:** 3-5일 (per-site CLI + array wrapper + regression test).

---

## Q3 — s05 Refactor Scope (bio_king) — 원칙 적용 + 1 예외

### Q3 원칙 재확인 (Round 1 §5): **"MVP 통과 전 refactor 금지"**

8개 AC 중 refactor 자체가 기여하는 AC는 **AC-6 (audit/reproducibility)** 과 **AC-8 (failure
transparency)** 뿐. AC-1~AC-4에 직접 기여 없음. 따라서 **기본은 Phase 2로 연기**.

### Q3 세부 판정

**Q3a. s05 8-10 모듈 분할 (bio_king §2) — 기각 (Phase 2)**
- 비용: PR 8-10회, 대규모 diff, 회귀 위험. AC-1 현재 5/6 GT 상태에서 refactor 도중 발생할
  회귀 확률이 **GT 회복보다 크다**. bio_king 본인 §1 표에서도 "큰 diff, plan+TDD 필수"로
  리스크 명시.
- 재검토 시점: AC-1 6/6 + AC-2 ≤ 2 FP 달성 후 (v1.0 production 태그 이후).

**Q3b. `compute_verdict()` pure-function 분리 (bio_king §1 target #2, §4 smell #1) — MVP 예외 채택**
- **AC-6 audit trail에 직접 기여**: verdict 계산 로직이 318줄의 `generate_report` 안에 묻혀
  있어 보고서에 "왜 이 verdict로 판정되었나"를 체계적으로 기록 불가. pure function 분리 시
  `FilterEvidence` dataclass → verdict 테이블이 report 하단에 자동 첨부 가능.
- **AC-8 failure transparency**: TIMEOUT 시점에 partial verdict를 계산할 수 있게 됨 (현재는
  전체 BLAST 끝난 후에만 verdict 결정).
- bio_king §3 TDD 목록 #1 (7개 시나리오)이 즉시 부착 가능 → regression 방어.
- 비용: 1주 (1 모듈 + 7 test + regression validation).

**Q3c. BLAST 공통 래퍼 `run_blastn()` (bio_king target #3) — MVP 채택**
- 29곳 중복 제거, BUG-2/15 류 재발 방지. **AC-6 (atomic write, 0-byte DB 방어)에 직접 기여**.
- haibao §1.3 비판 ("BLAST가 4-5회 호출, 결과 간 일관성 없음")과도 정합. 공통 래퍼 + in-memory
  result caching이 공통 전제.
- 비용: 2-3일 (wrapper + 기계적 치환).

**Q3d. extra-DB source 태깅 3-way 분리 (bio_king smell #5) — MVP 채택**
- Q1b canonical-triplet rule 구현 시 **선행 조건**. "element_db" vs "common_payload" vs
  "sample_contig" vs "canonical_triplet" 4-way source 태깅이 verdict 정책의 기초.
- 비용: 1일.

---

## Q4 — haibao 제안 판정

### Q4a. Level 1 patch 5개 (haibao §3 Level 1) — **3/5 채택**

| # | 제안 | 판정 | Rationale |
|---|------|------|-----------|
| L1-1 | `cluster_window` adaptive (coverage 기반) | **Phase 2 연기** | BUG-7 리버트가 보여주듯 site-discovery 파라미터 변경은 GT 회귀 직결. MVP 안정성 우선. |
| L1-2 | BLAST → LAST 스위치 옵션 (`--aligner last|blast`) | **Phase 2 연기** | LAST는 blastn 대비 2-5x 빠름이나 element_db 검증 미완료. **gmo_expert가 131 element에 대한 LAST regression 검증 먼저.** |
| L1-3 | **Assembly round 8 → 3** | **MVP 채택** | UGT72E3 per-site 14min → 5min. Q2 array와 독립 효과. 단 **측정 선행**: s05_stats.txt 분석으로 round 3 이후 growth 기여 실측 (haibao §5 도전질문 1). 실측에서 기여 ≥ 10%이면 round 5로 조정. |
| L1-4 | **중복 BLAST 통합 (in-memory 3-filter)** | **MVP 채택** | haibao §1.3 / bio_king §1 target #3과 중첩. AC-4 직접 기여. Q3c의 `run_blastn()` wrapper와 같이 구현. |
| L1-5 | **transgene_db build-time dedup (cd-hit-est)** | **MVP 채택** | BUG-3 근본 해결. regulatory MUST #3 (DB versioning)과 정합. `element_db/` Makefile에 `cd-hit-est -c 0.95` 스텝. 비용 ≤ 1일. |

### Q4b. Level 2 module rewrite (dual-anchor + PE discordancy) — **Phase 2 연기**

**Rationale:**
- haibao §1.1 비판 ("SR-only는 2010년대 방식")은 **기술적으로 맞음**. GRIDSS/Manta는 multi-signal.
- 그러나 **MVP 통과 관점에서는 비용 과다**:
  - 비용: 2-3주 + 전면 regression 검증. 현 5/6 GT 회귀 위험 **GT 증가분보다 클 가능성**.
  - AC-1 개선분이 측정되지 않음. haibao §5 도전질문 3 ("PE 추가하면 BUG-7 regression 제거되나?")이
    **답변 없는 가설**.
- **대안:** MVP 통과 후 v1.1에서 **PE discordancy만 단독 추가** (dual-anchor rewrite 없이).
  cluster_window adaptive와 함께 Phase 2 묶음.

**재검토 조건:** v1.0 production 이후 논문 투고 전. haibao §7 우선순위 2 ("1개월 내")에 동의.

### Q4c. Level 3 full redesign (GRIDSS-style breakpoint graph + macrosynteny + long-read) — **기각 (v2.0 grant proposal)**

**Rationale:**
- haibao §3 Level 3 본인 평가: "2-3개월, 기존 95% 재작성, high risk". 격리소 production
  배포 일정(수개월) 및 MVP 통과 기준에 부합 안 함.
- **Macrosynteny signal (haibao §4)**: 기술적으로 매력적이나 현재 6 GT 중 **cucumber
  B10v3는 contig-level 조립 (haibao §4.1: "synteny가 약한 case")으로 적용 제한**. 나머지 3 host
  (rice/tomato/soybean)도 synteny block 구축 1회 추가 공수 필요. AC-1 직접 기여가 불명확.
- **Long-read readiness**: 한국 격리소 현재 **Illumina 150bp PE 표준**이 고정 (Round 1 §1
  AC-5 "표준 Illumina 150bp PE" 제약). ONT MinION 현장 전환은 USDA/APHIS 트렌드로 한국보다
  1-2년 뒤 도입 예상. **v2.0 grant proposal로 분리**.
- haibao §5 도전질문 `won_yim` #3 ("regulatory가 bp-resolution 정말 필요한가")은 **Yes**:
  KCGP/CRL event-specific assay가 LB/RB junction bp 좌표를 요구. ±10bp fallback 불가.
  단, 이것이 GRIDSS 필요로 이어지는 것은 아님. 현 방식으로 bp 수준 CAND 3/6 이미 달성.

### Q4d. haibao 반박 대응 논리 (Round 3 대비)

haibao가 "Level 2를 MVP에 넣어야 한다"고 재주장 시 예상 논점:
1. **"BUG-7 class regression이 reviewer에게 잡힐 것"** → 응답: "논문 리비전 주기는 ≥ 3개월,
   MVP 배포는 그 전. reviewer 대응은 v1.1 PE discordancy로 충분."
2. **"per-site 14min TIMEOUT은 알고리즘 문제"** → 응답: "Q2 array + L1-3 round 축소로
   48h TIMEOUT 해결. 알고리즘 교체 없이도 AC-4 달성. ROI 계산이 유리."
3. **"WT 없을 때 MAPQ=60 FP 무방비"** → 응답: "현 6 GT 중 WT 없는 샘플 = soybean/corn
   2 exploratory. AC-2 목표(≤ 2 FP / 샘플)에서 WT 없이 달성 가능한지가 **Round 3 Q1-bonus
   검증 대상**. 달성 불가 시 v1.1에서 WT-free fallback 구현."

**won_yim 최종 입장:** haibao Level 2/3은 기술적 완성도 추구. 그러나 **production
deployment = 규제 기관 제출 가능한 최소 기능** 이 MVP 정의. 기술적 완성도는 v1.1/v2.0에서.

---

## Q5 — Regulatory MUST 8개 분류

Round 1 §4 기반, MVP 필수 vs Phase 2 vs 폐기:

| # | 요건 | 분류 | Rationale |
|---|------|------|-----------|
| R-1 | **Input SHA-256 hash 기록** | **MVP 필수** | 법적 책임 소재 증거. 구현 비용 ≤ 0.5일 (`run_pipeline.py` 시작부 sha256sum + 로그). AC-6 직접. |
| R-2 | **Pipeline commit hash lock** | **MVP 필수** | `git rev-parse HEAD` + clean 체크. 비용 ≤ 0.2일. AC-6 직접. |
| R-3 | **BLAST DB versioning (md5 + build date)** | **MVP 필수** | Q4a L1-5 (cd-hit dedup) 구현 시 함께 처리. 재현성 핵심. AC-6. |
| R-4 | **Software version manifest** | **MVP 필수** | `bwa --version`, `samtools --version`, `blastn -version`, `spades.py --version` 등 로그. 비용 ≤ 0.2일. AC-6. |
| R-5 | **PDF-ready report (표준 4쪽 포맷)** | **Phase 2** | AC-5 (operator UX) target이나 MVP는 텍스트 report로 충분. regulatory 기관이 즉시 PDF를 요구하지 않음 (제출 시점에 manual 변환 가능). |
| R-6 | **KCGP element nomenclature 준수** | **Phase 2 (gmo_expert 승인 필요)** | element_db 현재 이름 체계가 EUginius. KCGP 표준 맵핑 테이블 필요. gmo_expert §2 부록 B의 accession과 KCGP 이름 교차검증 공수 ≤ 1주. Phase 2. |
| R-7 | **False-negative 명시 (coverage/MAPQ/LOD)** | **MVP 필수** | AC-8 (failure transparency) 직접. `compute_verdict`가 "not detected" 시 검출 한계 (min coverage, effective MAPQ) 계산 로그. Q3b와 함께 구현. |
| R-8 | **Chain of custody log (operator ID, SLURM JobID, timestamp)** | **Phase 2** | regulatory 요건이나 SLURM JobID + operator는 이미 `.err` 파일에 포함. 표준화된 수집은 Phase 2. |

**MVP 필수 = 5개 (R-1, R-2, R-3, R-4, R-7)**. **Phase 2 = 3개 (R-5, R-6, R-8)**. 폐기 = 0개.

**총 비용:** MVP 필수 5개 합산 ≤ 1.5일 (R-3은 Q4a L1-5와 중첩이므로 추가 비용 제로).

---

## 3명 teammate 미답변 질문에 대한 PI 입장

### gmo_expert가 won_yim에게 질문한 것 (gmo_expert §7 4번 항목)

1. **"Sensitivity vs 특이도 법적 허용치"**
   - 식약처 입장: **false-negative가 치명적** (수입 LMO 누락 = 검역 실패). 따라서 AC-1 sensitivity
     목표를 100% (core 6 GT)로 엄격하게. AC-2 specificity는 manual review로 보완 가능하므로
     ≤ 2 CAND FP / 샘플.
   - 수치 해석: **검출 실패 비용 >> FP 비용**. FN 1건 허용 시 국제 통상분쟁 가능.

2. **"Cas9-present 표기의 legal risk"**
   - Q1a 채택으로 해결. **MVP 배포 전까지 config flag 기반 표기 금지.** 현재 보고서에서
     "Cas9-present (by config)" 명기 + "서열 기반 미검증" 표시 필수.

3. **"Subsampled coverage가 격리소 input 품질과 일치하나"**
   - 부분 일치: 격리소는 MiSeq/NovaSeq 300bp/150bp 혼재, quality varies. 3x~15x subsampled는
     량 기반, quality는 반영 안 됨. **Phase 2에서 quality-degraded synthetic read 추가 검증
     필요.** MVP는 3-15x coverage matrix로 충분.

### haibao가 won_yim에게 질문한 것 (haibao §5 4번 항목)

1. **"WT 없을 때 production criterion에 포함"**
   - **동의**. AC-2 (specificity) 목표 달성을 WT-free 샘플 (soybean/corn exploratory)에서도
     검증해야 함. Round 3에서 해당 샘플의 FP budget 별도 측정 추가.

2. **"FP rate 숫자 목표 정의"**
   - Round 1 §1 AC-2 이미 정의: **≤ 5 CAND FP / 샘플 (MVP), ≤ 2 CAND FP / 샘플 (production)**.
     cucumber_line225 8 CAND (GT 2)는 **MVP 경계 + 1 over**. resume.md:146-149 (new CAND
     10건 remote BLAST 필요)이 제대로 수행되어야 AC-2 확정.

3. **"bp-resolution이 정말 필요한가"**
   - **Yes**. 위 Q4c에서 언급. CRL event-specific assay가 LB/RB junction bp 요구. 단 ±10bp
     tolerance는 허용 (short-read assembly의 typical breakpoint 불확실성).

---

## Round 3 Synthesis 준비 목록 (won_yim 담당)

다음 round에서 `docs/team-review/decisions.md` 작성 시 포함할 항목:

1. **MVP 범위 확정 리스트** (이 문서의 "MVP 채택" 전체 취합, 공수 합산)
2. **Phase 2 deferred list** (Phase 2 분류 전체)
3. **AC-1~AC-8 목표치 vs 예상 달성치** (이번 결정 후 도달 예정 점수)
4. **각 teammate R2 반론 수렴 결과** (bio_king: 왜 예외적으로 compute_verdict만 / haibao: 왜
   Level 2 연기 / compute_eng: 왜 array 우선 / gmo_expert: 왜 AtYUCCA6 유보)
5. **Implementation plan 12 task** (공수 총합 목표 ≤ 2주, single-owner 할당)
6. **Regulatory 제출 가능 선언문 초안** (MVP 완료 시점 기준 "이것은 격리소 production-ready"의
   공식 근거)

---

**Round 2 완료. team-lead에 경로 보고.**
