# Round 1 — won_yim (PI / end-user) — Production-Ready 기준과 검증 샘플 Rationale

**역할:** RedGene 프로젝트 PI 및 end-user (한국 GMO/LMO 격리소 검사 assay 담당자 대리)
**작성일:** 2026-04-16
**베이스 커밋:** `729835d` (main)
**관점:** 농림축산검역본부(APQA) / 국립종자원(KSVS) / 식품의약품안전처(MFDS) 제출용 분석 결과물 품질 기준

이 노트는 Round 2에서 bio_king·gmo_expert·compute_eng·haibao의 Round 1 결과물을 scope / cost / benefit으로 채택·기각하기 위한 **의사결정 기준**을 사전 고정합니다. Round 1 단계에서는 타 teammate 노트가 아직 도착 전이므로 일부 항목은 "Round 2 결정"으로 유보합니다.

---

## 1. Production-Ready 수용 기준 (Acceptance Criteria)

격리소 현장 배포 기준으로 8개 metric을 정의합니다. 각 항목은 **측정 가능한 숫자 threshold**와 **현재 수준 대비 gap**을 포함합니다.

| # | Metric | 정의 | MVP 최소 | 목표 (Production) | 측정 방법 |
|---|--------|------|---------|------|-----------|
| AC-1 | **Sensitivity (recall)** | Ground-truth insertion 중 CANDIDATE 또는 CANDIDATE_LOW_CONF로 보고되는 비율 | ≥ 80% (core 5 host), ≥ 1 CAND / site | 100% core 6 host @ ≥15x | 8개 GT 샘플 (rice·A2_3·cuc 212/224/225·tomato_A2_2·corn_ND207·+1 soybean GT 확보 시) 에 대해 per-site hit count |
| AC-2 | **Specificity (per-sample FP budget)** | 한 샘플에서 CANDIDATE로 보고되는 host-only 위양성 개수 | ≤ 5 CAND FP / 샘플 | ≤ 2 CAND FP / 샘플 | 수동 review + remote BLAST로 FP flag, rice·A2_3·cuc에서 count |
| AC-3 | **Annotation completeness** | CANDIDATE 중 element_name이 `UNKNOWN`/빈 값이 아닌 비율 | ≥ 70% | ≥ 90% | s05 report의 `element_annotation.tsv` 파싱 |
| AC-4 | **Turnaround (wall time)** | reads → final report (steps 1-7) 전체 시간, 16 CPU / 64G / SLURM 기준 | ≤ 48h / 샘플 | ≤ 24h / 샘플 (host ≤ 1 Gbp); ≤ 48h (≥ 1 Gbp: soybean/corn) | `sacct --format=Elapsed` |
| AC-5 | **Operator UX** | config.yaml 편집 분량 + 커맨드 수 | 1 sample entry 추가 + 1개 sbatch | zero-config 추가 (template 복사), 1 CLI | manual 실습 (외부 검사원 1인) |
| AC-6 | **Reproducibility / Audit trail** | 제출 보고서에 input hash, pipeline commit hash, DB hash, software version이 모두 포함 | 4종 중 3종 이상 | 4종 전부 + SHA-256 lockfile | report header 점검 |
| AC-7 | **Coverage robustness** | ≥ 15x (rice), ≥ 10x (cuc/soybean/corn), ≥ 10x (tomato)에서 GT insertion CAND 회수 | 15x 수준에서 only | 15x·10x 모두 (5x·3x degrade 허용) | coverage sensitivity 매트릭스 (config에 9 variant 이미 정의) |
| AC-8 | **Failure transparency** | 실패 모드(TIMEOUT, OOM, 0-CAND)에서 end-user가 원인 즉시 파악 가능 | structured error + 권장 재시도 | 자동 retry with adjusted mem/time + bug code 참조 | `run_pipeline.py` exit code + log tail |

### MVP vs Production 이진 판정 방식
- **MVP 통과 = AC-1~AC-5 min threshold 전부 만족 + AC-6 3/4 이상.**
- **Production 통과 = 8개 AC target 전부 만족.**
- 현재 상태는 "MVP 부분 통과"로 평가 (아래 gap 표 참조).

---

## 2. 검증 샘플 Rationale

### 2.1 왜 이 5 host인가 — 격리소 실제 유입 패턴과의 매핑

| Host | 격리소 관련성 (근거) | 검사 우선순위 | 현 샘플 커버리지 |
|------|-----------|----------------|----------------|
| **Rice** (*Oryza sativa*) | 한국 주식작물; GMO 쌀은 국내 상업재배 불허(식약처) → 수입 모니터링 최상위. 2019 IRRI Golden Rice 식용승인 이슈 재점화 | 1순위 | rice_G281 (lactoferrin RNAi, GT 있음) ✅ |
| **Tomato** (*Solanum lycopersicum*) | 유전자교정(CRISPR) 식품 사례 다수 (일본 고GABA 토마토 2021 승인). 교정·전이유전자 혼합형 검출 수요 | 1순위 | tomato_Cas9_A2_2/A2_3 (GT 있음), tomato_WT (음성 컨트롤) ✅ |
| **Cucumber** (*Cucumis sativus*) | 한국 종묘 수입 주요 품목; thaumatin(감미 단백질) 외래 원소 검사 연구 사례 | 2순위 | cucumber_212/224/225 (GT 3종, single/multi copy+backbone) ✅ |
| **Soybean** (*Glycine max*) | 수입량 최대 GM 작물(연 ~130만 톤). 다수 승인(Roundup Ready 등). 미승인 event 탐지가 실제 임무 | 1순위 | soybean_AtYUCCA6·UGT72E3 (GT 미확립, exploratory) ⚠️ |
| **Corn** (*Zea mays*) | 수입 사료용 GM 옥수수 연 ~1,000만 톤. MON810·NK603 등 이미 승인되어 있음. unapproved stacked event 탐지가 현업 수요 | 1순위 | corn_ND207 (Sci Rep 2025) ⚠️ Phase 1.5 |

**결론:** 5 host 선택은 한국 수입/재배 통계와 규제 우선순위를 기준으로 정당화됨. **누락 우선순위**는 Canola·감자(추후 Phase 2).

### 2.2 왜 이 전형적 vector 구성인가

| Vector 구성 | 대표 샘플 | 규제 이슈 | 검증 목적 |
|-------------|-----------|-----------|-----------|
| **T-DNA single copy + 식물 유래 promoter (Ubi1/Act1)** | rice_G281, cucumber_line212/224 | 가장 흔한 1세대 GMO 구조. 식물유래 promoter로 host cross-reaction FP 다수 | AC-2 (specificity), WT-homology filter 유효성 검증 |
| **T-DNA 2-copy head-to-head + backbone read-through** | rice_G281 (2-copy), cucumber_line225 | 의도되지 않은 백본 통합 규제 위반 사례 | Multi-copy assembly (s05) + alt-locus filter (Filter D) |
| **CRISPR/Cas9 editing (T-DNA 있음)** | tomato_A2_2, tomato_A2_3 | 전이유전자 제거 후 교정만 남긴 "null segregant" 식품 판정 | s06 indel detection + s04 host-mapping |
| **CRISPR edit only (T-DNA segregated out)** | tomato_A2_1 | null segregant 증명 부담 | s06 positive, s05 no-CAND 동시 충족 |
| **Bar marker / 박테리아-바이러스 공통 요소 (CaMV35S, nos)** | soybean, corn_ND207 | Element DB 일반 요소 검출 (핵심 regulatory signature) | common_payload.fa 매칭, AC-3 annotation |

**결론:** 현재 샘플은 한국 regulatory가 주목하는 vector 구성 대부분을 커버. **누락**: RNAi hairpin-only (non-T-DNA), gene-stack (3+ event combined). Phase 2 검토.

### 2.3 Ground-truth 샘플과 Exploratory 샘플의 역할 구분

| 구분 | 샘플 | 역할 | 수용 기준 |
|------|------|------|-----------|
| **GT anchor (AC-1 직접 측정)** | rice_G281 Chr3:16,439,674 / tomato_A2_3 ch01:91,002,744 / tomato_A2_2 ch08:65,107,378 / cuc_212 Chr6 / cuc_224 Chr2 / cuc_225 Chr2 (2-copy+backbone) | 파이프라인 recall 정량 측정 | AC-1 threshold 계산 기반 (6/6 CAND recovery 목표) |
| **Exploratory (no GT — FP 상한 + turnaround만 측정)** | soybean_AtYUCCA6 / soybean_UGT72E3 / corn_ND207 | 실제 운영 환경(GT 없음)에서 FP rate와 wall time 측정. "0 CAND"라도 FP budget ≤ AC-2 threshold면 OK | AC-2, AC-4 측정만 |
| **음성 컨트롤 (specificity floor)** | tomato_WT | T-DNA 없음. CAND = 0 (±1) 이어야 함 | AC-2 strict 검증 (WT에서 CAND ≥ 1 발생 시 production reject) |
| **Coverage stress test** | rice_G281_{15/10/5/3}x, tomato_A2_1_{15/10/5/3}x, tomato_A2_3_{5/3}x, cuc_224_{15/10/5/3}x, soybean_UGT72E3_{15/10/5/3}x, corn_ND207_{5/3/1}x | AC-7 측정 | ≥ 10x에서 GT 회수, 3x에서 graceful fail |

**결론:** soybean 0-CANDIDATE 문제(resume.md:96-105)는 GT anchor 샘플이 아니므로 **AC-1 실패로 간주하지 않음**. 단 AC-4 (wall time) 및 AC-2 (FP budget)만 적용. **따라서 Round 2에서 soybean-only 튜닝을 위해 s05 re-architecture 요구하는 제안은 우선순위 낮음 — scope out 유력.**

---

## 3. 현재 수준 vs 수용 기준 Gap

`resume.md`의 2026-04-16 검증 결과(`729835d` 기준) 기반:

| AC | 현재 점수 | 목표 | Gap | 우선순위 |
|----|-----------|------|-----|----------|
| **AC-1 Sensitivity** | 5/6 GT (rice·A2_3·cuc×3). tomato_A2_2 미재검증. Expected GT recall 83% | 100% | tomato_A2_2 재검증 (1 sample) | P1 |
| **AC-2 Specificity** | rice 3 CAND (1 true + 2 review 필요), cuc_225 8 CAND (2 expected), A2_3 2 CAND (1 true + 1 review). **FP 상한 미확정** | ≤ 2 CAND FP / 샘플 | New CANDIDATE 10건 remote BLAST 필요 (resume.md:146-149) | P1 |
| **AC-3 Annotation completeness** | 집계 미존재. 체감: rice 3/3 annotated, cuc_225 일부 UNKNOWN 유력 | ≥ 90% | 자동 집계 스크립트 추가 필요 | P2 |
| **AC-4 Turnaround** | rice 1h15m, A2_3 29m, cuc 9-13h, soybean 24h+ TIMEOUT | ≤ 48h 전부 | soybean 48h 재제출 or s05 per-site parallelize | P1 (compute_eng 담당) |
| **AC-5 Operator UX** | config.yaml에 host/construct/reads 수동 기재. `run_pipeline.py --sample X --steps 1-7` 1줄. 합격 직전 수준 | 현재 MVP 만족 | template 추가로 production 승격 가능 | P3 |
| **AC-6 Audit trail** | git commit hash는 로그에 있음. DB hash·input hash·version lock 부재 | 4/4 필수 | 보고서 header 추가 (s07 이후 단계) 필요 | P1 (regulatory 제출용 필수) |
| **AC-7 Coverage robustness** | 9 variant config 존재, 실행 미완 | ≥10x GT recall | subsampled 검증 batch 실행 필요 (resume.md deferred 8) | P2 |
| **AC-8 Failure transparency** | TIMEOUT = exit 0 + 불완전 report. OOM = SLURM kill, 진단 로그 부재 | structured error | run_pipeline.py wrapper 강화 | P2 |

**요약:** MVP 합격선에는 **AC-2 FP 검토 + AC-4 soybean 재제출 + AC-6 audit header** 3건 완료 시 도달. Production까지는 AC-3/7/8까지 추가 3건.

---

## 4. Regulatory 제출 관점에서 반드시 필요한 기능 (MUST)

한국 3개 기관(APQA / KSVS / MFDS) 제출 가능성을 전제로 **비타협 요건** 명세:

1. **Input SHA-256 hash 기록** — reads fastq.gz, host reference, construct reference, element DB 전부. 법적 책임 소재 증거.
2. **Pipeline commit hash lock** — 최종 보고서에 `git rev-parse HEAD` 및 `git diff --stat` (clean 여부). reproducibility 필수.
3. **BLAST DB versioning** — `nt` / `element_db` / `common_payload.fa`의 build date + md5. 동일 샘플을 6개월 후 재분석해도 같은 결과를 재현할 수 있어야 함.
4. **Software version manifest** — bwa / minimap2 / samtools / SPAdes / blastn / fastp 버전 `--version` 출력 저장.
5. **Electronic signature / PDF-ready report** — 보고서 한 장(또는 표준 4쪽) 포맷 고정. 검사관이 직접 수정할 수 없는 format (PDF/A 권장).
6. **KCGP (한국유전자검사평가원) 제안 element nomenclature 준수** — element_db 이름 체계를 국제 표준(bar·nptII·hpt)로 통일. internal code 절대 금지.
7. **False-negative 명시** — "검출 불가" 판정 시 coverage / MAPQ / 검출 한계를 명시. "not detected" ≠ "not present" 구분.
8. **Chain of custody log** — 샘플 수령→QC→분석→보고 각 단계 timestamp + operator ID. SLURM JobID도 함께.

**Round 2에서 충돌 예상:** haibao가 "아키텍처 교체(Nextflow·Snakemake 전환)"를 제안할 수 있음. 그러나 위 8개 요건은 **현 Python+subprocess 아키텍처로 1-2주 내 구현 가능**. 아키텍처 교체 ROI는 regulatory 요건을 이미 충족한 뒤에만 검토. → **현 시점 Scope out 유력.**

---

## 5. 다른 Teammate 제안에 대한 Scope 가이드 (1차 의견, Round 2에서 재검증)

| 제안 유형 | 채택/유보/기각 | Rationale |
|-----------|----------------|-----------|
| **s05 모듈 분할 (bio_king 예상 제안)** | **유보 → Round 2 결정** | 3806줄 단일 스크립트는 유지보수 부담이지만, MVP AC 충족에 직접 기여 아님. 단, test coverage 확보가 AC-6 "reproducibility"에 연결되면 채택. 분할 rationale이 "향후 기능 추가"뿐이면 Phase 2. |
| **Element DB governance (gmo_expert 예상 제안)** | **채택 유력** | AC-3 (annotation completeness) + regulatory 요건 #6 (KCGP nomenclature)에 직접 기여. 단, element DB를 샘플별 업데이트는 이미 기각됨 ([Element DB strategy memory](../../../.claude/projects/-data-gpfs-assoc-pgl-develop-redgene/memory/feedback_element_db_blast.md) 참조). 버전 관리와 md5 lockfile 방향 유력. |
| **SLURM 자원 최적화 (compute_eng 예상 제안)** | **채택 유력** | AC-4 (turnaround)에 직접 기여. 단, 4GB/thread 원칙 준수 ([memory_threads memory] 참조) 및 BWA -t 16 유지 ([bwa_threads memory] 참조). 소이빈 per-site parallelize는 수락, Nextflow 도입은 Round 2에서 ROI 별도 검증. |
| **알고리즘 교체 (haibao 예상 제안)** | **대부분 기각 예상** | MVP 통과 전 아키텍처 교체 금지. 단, s05 false positive filter 중 **alternative-locus filter (Filter D)** 개선안은 AC-2에 직접 기여하므로 채택 가능. |
| **Visualization 우선 구현** | **기각 (Phase 2)** | 한국 regulatory는 텍스트 보고서 기반. 시각화는 논문/학술용. MVP out. |
| **Coverage sensitivity 자동 batch** | **채택** | AC-7 측정에 필수. 9 variant 이미 config에 존재 → 실행 script만 필요 (ROI 높음). |

**결정권자 원칙:**
- **"MVP 통과 전 refactor 금지"** — AC-1/AC-2 미충족 상태에서 3806줄 분할은 시간 낭비.
- **"regulatory 요건 우선"** — AC-6 audit trail이 bio_king refactor보다 우선.
- **"GT anchor 6종 100% 우선"** — soybean/corn exploratory는 budget 남을 때.

---

## 6. 다른 Teammate에게 묻고 싶은 것

### bio_king (s05 refactor target 분석)
1. 3806줄 중 실제 CANDIDATE 판정 로직은 몇 줄인가? `_should_replace` / `classify_site_tiers` / `annotate_insert` / filter A-D 가 전체의 몇 %인가?
2. 리팩터링으로 얻는 **test coverage 증가분**이 측정 가능한가? 현재 pytest 10 test의 라인 coverage는?
3. s05 refactor가 AC-1 recall 개선에 직접 기여하는 경로가 있는가, 아니면 유지보수용인가?

### gmo_expert (element DB 정합성 + FP 케이스)
1. EUginius 131 elements 중 규제 기관이 "핵심 signature"로 요구하는 subset은 무엇인가? (예: bar, nptII, CaMV35S, nos = 필수 4종)
2. Element DB에 **KCGP / ISO 21569 / Codex Alimentarius 표준 명명**이 반영되어 있는가? (Regulatory 요건 #6)
3. Cucumber line225의 "8 CANDIDATE (GT=2)"는 backbone integration일 가능성이 있는가, 아니면 element DB cross-reaction인가?
4. Soybean 0-CANDIDATE는 element DB에 Gmax-specific promoter가 과다하게 포함되어 있을 가능성?

### compute_eng (HPC/SLURM 최적화)
1. soybean per-site 14분 × 105 sites를 24h 내로 맞추려면 s05 per-site를 SLURM array 분할하는 것이 실현 가능한가? (resume.md:108-112)
2. 현재 64G / 16 CPU / 24h가 AC-4 (≤ 48h) 달성에 병목이 어디인가? (step 4 host-mapping 5-7h vs step 5 per-site)
3. 4GB/thread rule과 BWA -t 16 유지하면서 추가 throughput 확보 방안이 있는가?

### haibao (알고리즘 재평가)
1. s05 Filter D (alt-locus minimap2 check)의 FP detection rate를 측정할 수 있는가? Cuc_225의 8 CAND 중 몇 개가 alt-locus에서 회수되는가?
2. Assembly-based junction detection 대안(de novo whole + graph-genome)의 **ROI**는? (MVP 통과 대비 필요 공수)
3. CRISPR s06 pileup-based indel detection을 대체할 알고리즘이 있는가? 현재 구현의 검출 한계(길이·MAPQ·heterozygosity)는?

---

## 7. Round 2 이후 최종 Synthesis 계획 (won_yim 담당)

Round 2에서 teammate 4명의 제안을 받은 뒤:
1. 각 제안을 **AC-1~AC-8 기여도**로 점수화 (직접·간접·무관).
2. **"MVP 통과에 필요한 최소 변경"** 목록 고정 (P1 tag). 추정 공수 ≤ 2주.
3. **Phase 2 deferred 목록** 별도 유지. 아키텍처 교체·visualization·exploratory host 추가 여기로.
4. 최종 결정 결과를 `docs/team-review/decisions.md` (Round 3)로 정리해 team-lead에게 제출.

---

**Round 1 종료. 대기 상태로 전환.** Round 2에서 taskUpdate 후 teammate 제안 수령 + scope 판정.
