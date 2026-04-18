# RedGene Team Review — Team Consensus (최종 결정문)

**문서 소유자:** won_yim (PI / end-user, 한국 격리소 GMO/LMO 검출 assay)
**Date:** 2026-04-16
**Team:** redgene-review-2026-04-16 (bio_king / gmo_expert / compute_eng / haibao / won_yim)
**기반:** Round 1-3 전 teammate 노트 + team-lead brief
**권위:** 본 문서는 MVP v1.0 release 기준 **최종 합의문**. Round 4 구현 시 판정 근거.

---

## §1. Executive Summary

**결정 한 줄:** RedGene MVP v1.0은 "AC-1 6/6 GT recall + AC-2 ≤ 5 FP / 샘플 + AC-4 ≤ 48h + AC-6 4-field audit header" 를 기준으로 정의하며, 12 task (≤ 10 person-day, 4명 parallelize로 1주)로 달성한다.

**MVP 범위 (이번 주 T1-T12):**
- **신규 기능 2**: canonical-triplet rule (config-driven), element vs host 100%-ortholog pre-mask BED
- **성능 개선 3**: per-site SLURM array 병렬화, assembly `max_rounds` 축소, `transgene_db` cd-hit dedup
- **DB 확장 1**: SpCas9 + sgRNA scaffold + AtYUCCA6 + UGT72E3 + P0/P1 element 14 seqs
- **Refactor 최소**: 최소 4-way 모듈 분할 + `compute_verdict()` 부분 분리 + `VerdictRules` dataclass
- **Regulatory audit 4**: input SHA-256 + commit hash + DB md5 + software versions

**공수 합산:** ~9.1 person-day, 4명 parallelize → **5 영업일** (2026-04-20 ~ 2026-04-24).

**도달 예정 AC 점수:**
| AC | 현재 | MVP 목표 | 도달 예상 |
|----|------|----------|-----------|
| AC-1 Sensitivity | 5/6 GT (83%) | 6/6 (100%) | **6/8 ~ 7/8** (soybean exploratory 복구 기대) |
| AC-2 Specificity | cucumber 8 CAND (GT=2) | ≤ 5 FP / 샘플 | **≤ 3-4 전 샘플** (pre-mask 효과) |
| AC-4 Turnaround | UGT72E3 48h TIMEOUT | ≤ 48h | **UGT72E3 ≤ 4h** (-92%) |
| AC-6 Audit trail | 0/4 | 4/4 필수 | **4/4 완전 충족** |

**Out-of-scope (연기):** s05 8-10 전체 refactor → v1.1. haibao Level 2 (dual-anchor + PE) → v1.1/v2. haibao Level 3 (GRIDSS + long-read + macrosynteny) → v2.0 grant.

**Regulatory 제출 가능성:** MVP v1.0 완료 시 **한국 격리소 pilot assay 수준** 도달. 실제 배포는 외부 validation + 법률 검토 후 (§7 선언문).

---

## §2. 5 Teammate Agree/Disagree 매트릭스

Round 1-2에서 제안된 23개 주요 항목에 대한 5명 입장 + 최종 결정.
표기: ✅ agree / ⚙ modify / ❌ disagree / ― neutral.

### §2.1 Refactor / Code Quality

| # | 제안 | bio_king | gmo_expert | compute_eng | haibao | won_yim | **최종 결정** |
|---|------|----------|------------|-------------|--------|---------|---------------|
| A | s05 8-10 전체 모듈 분할 | ✅ P0 | ― | ― | ⚙ (Level 2 함께) | ❌ Phase 2 | **Phase 2 연기** |
| B | 최소 4-way 분할 (site/classify/per_site/annotate) | ✅ | ― | ✅ (array CLI 전제) | ― | ✅ MVP | **MVP 채택 (T7)** |
| C | `compute_verdict()` pure 함수 분리 | ✅ P0 | ― | ✅ | ✅ | ✅ MVP 예외 | **MVP 부분 채택 (T6)** |
| D | `run_blastn()` 공통 래퍼 | ✅ P0 | ― | ⚙ (ROI -5~7%만) | ✅ | ✅ | **MVP 채택 (T5 합체)** |
| E | DB 정규화 + in-memory `BlastHits` sharing | ✅ (R2 수정) | ✅ | ✅ | ✅ P2 | ✅ | **MVP 채택 (T5)** |
| F | 2-tier priority (element_db-family > univec) | ✅ | ✅ (R2 자진 제안) | ― | ✅ | ✅ | **확정 채택** |

### §2.2 Feature / Scoring Rule

| # | 제안 | bio_king | gmo_expert | compute_eng | haibao | won_yim | **최종 결정** |
|---|------|----------|------------|-------------|--------|---------|---------------|
| G | canonical-triplet rule (config-driven) | ✅ (`VerdictRules`) | ✅ (config 선호) | ― | ― | ✅ MVP | **MVP 채택 (T6)** |
| H | SpCas9 + sgRNA scaffold DB 추가 | ― | ✅ P0 | ✅ (+1% 비용) | ― | ✅ MVP | **MVP 채택 (T4)** |
| I | AtYUCCA6 / UGT72E3 CDS 복원 | ― | ✅ P1 | ✅ | ― | ✅ MVP (R3 승격) | **MVP 채택 (T4)** |
| J | element 100%-ortholog pre-mask BED | ― | ✅ (100%-ortholog만, EPRV 제외) | ― | ✅ P6 | ✅ MVP | **MVP 채택 (T9)** |
| K | MCscan Workflow 2 (EPRV synteny mask) | ― | ⚙ (Workflow 1만) | ― | ⚙ (v2) | ❌ v2 | **v2 이관** |

### §2.3 Performance

| # | 제안 | bio_king | gmo_expert | compute_eng | haibao | won_yim | **최종 결정** |
|---|------|----------|------------|-------------|--------|---------|---------------|
| L | per-site SLURM array 병렬화 | ⚙ (분할 먼저) | ― | ✅ 강력 (-92%) | ✅ | ✅ MVP 필수 | **MVP 필수 (T8, 임시 fan-out)** |
| M | Assembly `max_rounds=8→3` | ― | ― | ✅ (복리) | ✅ P1 | ✅ (측정 조건부) | **MVP 채택 (T2 측정 + T3)** |
| N | `transgene_db` cd-hit-est dedup | ✅ | ✅ (효과 미미) | ✅ | ✅ P2 | ✅ | **MVP 채택 (T5)** |
| O | s04 BWA → minimap2 전환 | ― | ― | ⚙ (s04만, Phase 1 반대) | ✅ 강력 | ⚙ PoC 조건부 | **rice_G281 PoC만 승인 (T12)** |
| P | s05 Phase 1 minimap2 dual-anchor | ― | ― | ❌ (+2-3x, MAPQ 손실) | ✅ (Level 2) | ❌ | **기각 (v1.1 이상)** |
| Q | BLAST 호출 통합 (in-memory 3-filter) | ✅ | ― | ⚙ (ROI -5~7%) | ✅ P3 | ✅ (E에 흡수) | **MVP 채택 (T5 합체)** |

### §2.4 Algorithm / Architecture (haibao Level 1-3)

| # | 제안 | bio_king | gmo_expert | compute_eng | haibao | won_yim | **최종 결정** |
|---|------|----------|------------|-------------|--------|---------|---------------|
| R | cluster_window adaptive (P4) | ― | ― | ― | ⚙ (회귀 위험 자인) | ❌ Phase 2 | **Phase 2 연기** |
| S | LAST 옵션 (P5) | ― | ― | ― | ⚙ (검증 미완) | ❌ Phase 2 | **Phase 2 연기** |
| T | k-mer k=15 → 21 | ― | ― | ― | ⚙ (측정 필요 자인) | ⚙ (측정 선행) | **Phase 2 (측정 선행)** |
| U | haibao Level 2 (dual-anchor + PE) | ― | ― | ❌ | ⚙ (v2 수용) | ❌ v1.1 이상 | **v1.1 PE subset만 검토, 전체 v2** |
| V | haibao Level 3 (GRIDSS + long-read) | ― | ― | ❌ | ⚙ (v2 grant 수용) | ❌ | **v2.0 grant 분리** |

### §2.5 Regulatory / Audit

| # | 제안 | bio_king | gmo_expert | compute_eng | haibao | won_yim | **최종 결정** |
|---|------|----------|------------|-------------|--------|---------|---------------|
| W | AC-6 audit header R-1~R-4 (hash 4종) | ― | ✅ (식약처 필수) | ✅ (≤ 1일) | ― | ✅ MVP 필수 | **MVP 채택 (T1)** |

**23개 항목 결정 요약:** MVP 채택 **12개**, Phase 2 연기 **5개**, v1.1 이관 **2개**, v2.0 grant 분리 **2개**, 기각/조건부 **2개**.

---

## §3. MVP 12 Task List (owner + effort + dependency)

| # | Task | Owner | Effort (p-day) | AC 기여 | Dependency | 예상 PR 수 |
|---|------|-------|----------------|---------|------------|-----------|
| **T1** | AC-6 audit header 4-field (`_write_audit_header()`: input sha256 + git commit + DB md5 + software versions) | compute_eng | 0.5 | AC-6 | 없음 | 1 PR |
| **T2** | `s05_stats.txt` aggregate → round-3 수렴율 측정 (resume.md "most converge by 6" 검증) | haibao | 0.3 | (T3 전제) | 없음 | 측정 보고 only |
| **T3** | `assemble_insert(max_rounds=3 or 5)` 1-line + regression | haibao | 0.5 | AC-4 | T2 | 1 PR |
| **T4** | DB 확장: SpCas9 + sgRNA scaffold + AtYUCCA6 + UGT72E3 + gmo_expert P0 6 + P1 8 (총 14 seqs) + `element_db_manifest.tsv` | gmo_expert | 1.0 | AC-1, AC-3, AC-6 | T1 (manifest 포맷) | 1 PR |
| **T5** | `transgene_db` cd-hit-est -c 0.95 dedup + 4-way source tag loader + `BlastHits` dataclass 최소 + 중복 BLAST 통합 (megablast 3→1) | bio_king | 1.0 | AC-1, AC-2, AC-4 | T4 | 1-2 PR |
| **T6** | `VerdictRules` dataclass + config.yaml schema + `DEFAULT_TRIPLETS` fallback + `compute_verdict(evidence, rules)` 부분 분리 | bio_king | 1.5 | AC-1 | 없음 | 1 PR + 5+ pytest |
| **T7** | 최소 4-way s05 모듈 분할 boundary (site_discovery / classify / per_site / annotate+report) — 전체 분할 NOT now | bio_king | 1.0 | (T8/T10 전제) | 없음 | 1 PR (shim 경계) |
| **T8** | per-site SLURM array 임시안 (`run_pipeline.py --fanout`: site list split + N-way sbatch via 독립 worktree) | compute_eng | 1.5 | AC-4 | T7 | 1 PR + `submit_s05_array.sh` |
| **T9** | element vs host 100%-ortholog pre-mask BED (haibao Workflow 1, blastn + bedtools merge) + `host_masked_rationale.tsv` 큐레이션 | haibao (script) + gmo_expert (rationale) | 1.0 | AC-2 | T4 | 1 PR + per-host BED + TSV |
| **T10** | s05 Phase 1 직후 BED intersect + `FALSE_NEGATIVE_MASKED` 태그 로직 | bio_king | 0.5 | AC-2 | T7, T9 | 1 PR |
| **T11** | UGT72E3 + AtYUCCA6 재제출 + rice/tomato/cucumber 재검증 (6 GT + 2 exploratory) | compute_eng | 0.3 submit + ~8-16h wall | AC-1, AC-4 | T3, T4, T5, T8, T9, T10 | 재검증 표 |
| **T12** | rice_G281 s04 minimap2 PoC (조건부 승인 실험) | compute_eng + haibao | 0.5 setup + ~4-6h wall | (s04 PoC) | T11 | PoC 보고서 |

**총 공수:** 9.1 person-day (코드 작업 7.3 + submit/PoC 1.8). 4명 parallelize 시 **5 영업일 가능**.

### Critical Path

```
Day 1 (Mon): T1 ─┐                     (compute_eng, 0.5d)
                 ├─ T7 ──────────┐     (bio_king, 1.0d)
                 └─ T6 ──────────┤     (bio_king parallel, 1.5d)
Day 2 (Tue): T2 ─→ T3 ──────────┤     (haibao, 0.3+0.5d)
                   T4 ──────────┤     (gmo_expert, 1.0d)
Day 3 (Wed): T5 ─┬──────────────┤     (bio_king, 1.0d, T4 후)
                 └ T9 ──────────┤     (haibao+gmo_expert, 1.0d)
Day 4 (Thu): T10 + T8 ──────────┤     (bio_king 0.5d + compute_eng 1.5d)
Day 5 (Fri): T11 → T12            (compute_eng wall time 활용)
Day 6 (Sat): won_yim release 판정
```

---

## §4. Phase 2 Deferred List (v1.1 + v1.2)

### §4.1 v1.1 로드맵 (MVP 직후 4-6주)

| # | 항목 | 공수 (p-day) | AC 기여 | 담당 예상 |
|---|------|-------------|---------|-----------|
| D1 | `compute_verdict()` pure 함수 full 분리 + `FilterEvidence` dataclass | 5 | AC-6, AC-8 | bio_king |
| D2 | DB 정규화 full + in-memory `BlastHits` sharing (T5 확장) | 5 | AC-4 | bio_king |
| D3 | s05 전체 8-10 모듈 분할 (site_discovery / classify / reads / assemble / refine / annotate / filters / verdict / report / legacy) | 10 | code quality | bio_king |
| D4 | R-7 FN 명시 — 정량 LOD (coverage × MAPQ × e-value) 보고서 삽입 | 3 | AC-8 | bio_king + gmo_expert |
| D5 | haibao Level 2 subset — PE discordancy signal (dual-anchor 제외) | 7 | AC-1, AC-2 | haibao |
| D6 | s04 BWA → minimap2 전환 full deploy (T12 PoC 성공 조건) | 2 | AC-4 | compute_eng |
| D7 | `crl_amplicons.fa` 82개 통합 | 2 | AC-3 | gmo_expert |
| D8 | Coverage sensitivity batch 실행 (9 subsampled variants) | 2 submit + wall | AC-7 | compute_eng |
| D9 | Visualization 스크립트 QA (editing profile/effects/summary) | 3 | AC-5 | bio_king |
| D10 | BUG-14 worktree db/ symlink + BUG-15 atomic write | 2 | robustness | compute_eng |
| D11 | CI smoke test (synthetic 5-min pytest + GitHub Actions) | 3 | regression 방어 | compute_eng |

### §4.2 v1.2 / Phase 2 (v1.1 이후 ~3개월)

| # | 항목 | 공수 (p-day) | 담당 예상 |
|---|------|-------------|-----------|
| D12 | haibao P4 cluster_window adaptive | 3 | haibao |
| D13 | haibao P5 LAST 옵션 (`--aligner last`) | 5 | haibao |
| D14 | haibao k-mer k=15→21 실험 + 적용 (측정 선행) | 5 | haibao |
| D15 | R-5 PDF-ready report | 3 | bio_king |
| D16 | R-6 KCGP element nomenclature 맵핑 | 5 | gmo_expert |
| D17 | R-8 Chain-of-custody log 표준화 | 3 | compute_eng |
| D18 | MCscan Workflow 2 (EPRV synteny mask) | 7 | haibao + gmo_expert |
| D19 | cucumber / soybean / corn host s04 minimap2 PoC 확장 | 2 each | compute_eng |
| D20 | Nextflow / Snakemake 마이그레이션 평가 | 10 | compute_eng |

**Phase 2 총 공수 예상:** 약 60 person-day (~12주).

---

## §5. Phase 3 / v2.0 Grant Out-of-Scope

v1.x 코드베이스에서 분리되는 **장기 vision / grant proposal 분량**. 현 MVP와 무관.

| # | 항목 | 예상 scale | Rationale |
|---|------|------------|-----------|
| V1 | **haibao Level 2 full rewrite** — dual-anchor minimap2 + PE discordancy + breakpoint graph (GRIDSS/Manta-style) | 2-3 months | SR-only → multi-signal detection. BUG-7 class 해결. 학술 논문 리뷰 대응. |
| V2 | **haibao Level 3 full redesign** — GRIDSS-style BreakpointGraph + host/construct dual-BAM + macrosynteny native integration | 3-6 months | 2024-2026 production SV caller 수준 재구축. Rust/C++ migration 가능. |
| V3 | **ONT/HiFi long-read preset** — `data_type: short|long|hybrid` config 분기, per-site assembly 우회, `medaka`/`clair3` 교체 | 2 months | USDA/APHIS 2024-2025 guidance. 격리소 future-proofing. |
| V4 | **Macrosynteny full integration (`s05_synteny_boost.py`)** — rice/sorghum, tomato/potato, cucumber/melon, soybean/phaseolus MCscan block → synteny_score → verdict 가중치 | 2 months | haibao Axis D. Plant-specific GMO tool의 orthogonal signal. |
| V5 | **SvABA / GRIDSS integration or Rust prototype** — 외부 production SV caller wrapping 또는 scratch rewrite | 3-6 months | 학술 reference implementation 교체. Java 의존 회피 시 Rust/C++. |

**총 estimated grant scope:** 12-18 months, 2-3 FTE. haibao가 별도 `v2-grant-scope.md` (Task #16 옵셔널)에 정리.

---

## §6. AC-1 ~ AC-8 상태표

| AC | 항목 | 현재 (2026-04-16, `729835d`) | MVP 목표 | MVP 후 도달 예상 | Production 목표 |
|----|------|--------------------------------|----------|-------------------|------------------|
| **AC-1** | Sensitivity (GT CAND 회수율) | 5/6 GT (83%); tomato_A2_2 미재검증, soybean 0-CAND | ≥ 80% (MVP), 100% (production) | **6/8 ~ 7/8** (84-88%, T4+T6+T9 효과) | 100% core 6 host @ ≥15x |
| **AC-2** | Specificity (CAND FP / 샘플) | rice 3 (2 review), cuc_225 8 (GT=2), A2_3 2 (1 review) | ≤ 5 FP / 샘플 | **≤ 3-4 전 샘플** (T9 pre-mask) | ≤ 2 FP / 샘플 |
| AC-3 | Annotation completeness | 집계 부재 | ≥ 70% (MVP), ≥ 90% (production) | **측정 스크립트만 추가 (v1.1 정량)** | ≥ 90% |
| **AC-4** | Turnaround (wall time) | rice 1h15m, A2_3 29m, cuc 9-14h, soybean 24h+ **TIMEOUT** | ≤ 48h / 샘플 | **UGT72E3 ≤ 4h (-92%), cuc ≤ 14h, rice ≤ 1h** | ≤ 24h (host ≤ 1 Gbp); ≤ 48h (≥ 1 Gbp) |
| AC-5 | Operator UX | `--sample X --steps 1-7`; config 편집 수동 | 1 CLI + template config | **MVP 현재 수준 만족** | zero-config 추가 |
| **AC-6** | Audit trail | git commit은 로그에만 | 4/4 필수 (R-1~R-4) | **4/4 완전 충족** (T1) | 4/4 + SHA-256 lockfile |
| AC-7 | Coverage robustness | 9 subsampled variant config 존재, 미실행 | ≥ 10x GT recall | **Phase 2 D8 이관** | ≥ 10x (cuc/soy/corn), ≥ 15x rice |
| AC-8 | Failure transparency | TIMEOUT = exit 0, 진단 부재 | structured error + bug code 참조 | **부분 (T1만 포함)** | 자동 retry + LOD 명시 |

**MVP release 최저선:** 굵게 표시한 4 AC (AC-1/AC-2/AC-4/AC-6) 통과 필수.

---

## §7. Regulatory 제출 가능 선언문 초안

**중요:** 본 선언문은 MVP v1.0 completed + 모든 T1-T12 통과 후에만 사용 가능.

### §7.1 공식 제출 선언문 (격리소 assay 파일럿 배포 시 1-2 문단)

```
【RedGene v1.0 파이프라인 — 한국 격리소 GMO/LMO 검출 파일럿 assay 기술 선언】

RedGene v1.0 (commit hash: <HASH>, release date: 2026-04-DD)은 assembly-based
Illumina WGS 분석을 통해 수입 식물체의 transgene insertion 존재를 검출하는
파일럿 도구입니다. 본 도구는 한국 5 host (rice / tomato / cucumber / soybean /
corn)를 대상으로 6개 ground-truth 유전자변형 이벤트에 대해 검증되었으며, AC-1
sensitivity 100% (6/6 GT 회수), AC-2 specificity ≤ 5 CAND FP / 샘플, AC-4
turnaround ≤ 48h / 샘플 (16 CPU / 64G RAM SLURM 환경), AC-6 reproducibility
(입력 SHA-256 + commit hash + DB md5 + software version manifest 4-field 보고서
삽입) 기준을 충족합니다.

본 도구의 CANDIDATE 판정은 element-level screening 결과이며, 최종 GMO/LMO 법적
판정은 식품의약품안전처 고시 및 CRL-GMOMETHODS qPCR event-specific assay로
확정하여야 합니다. 파일럿 배포는 농림축산검역본부 / 국립종자원 /
식품의약품안전처 담당 부서 검토 하에서만 유효하며, Korean quarantine assay
공식 승인은 외부 validation (≥ 3개 독립 샘플 세트, CRL reference amplicon
호환성 검증) 완료 후 본 배포로 승격됩니다.
```

### §7.2 제출 문서 첨부 필수 항목 (Release checklist)

MVP v1.0 release 시 **보고서 첨부 파일**:

1. `run_report.txt` — 샘플별 verdict + input/DB/software manifest (AC-6 R-1~R-4)
2. `site_classification.tsv` — Phase 1.5 tier 분류 결과 (transgene-positive 모든 site)
3. `insertion_<site_id>_report.txt` — per-site verdict + evidence (CANDIDATE/FALSE_POSITIVE/UNKNOWN)
4. `element_annotation.tsv` — 각 insert의 element hits (annotation completeness 측정)
5. `host_masked_rationale.tsv` — pre-mask BED의 각 구간별 "왜 mask했는지" 법적 근거 (T9)
6. `element_db_manifest.tsv` — DB name / md5 / build date / seq count (T4)

### §7.3 법적 면책 (필수 삽입)

```
본 파일럿 도구의 "not detected" 판정은 "not present"와 동일하지 않습니다.
검출 한계는 coverage (≥ 10-15x), MAPQ 품질, host reference 완성도, element
DB 포함 여부에 따라 달라집니다. 전체 GMO/LMO 판정의 법적 효력은 본 도구가
아닌 qPCR event-specific assay 및 식약처 고시 기준에서 확정됩니다.
```

---

## §8. Unresolved Items (Round 4 결정 대기)

MVP 구현 중 또는 MVP 후 즉시 판정 필요한 미결 항목.

### §8.1 AtYUCCA6 / UGT72E3 CDS 복원 — Round 3 잠정 채택 → T11 검증 후 최종

- **현 상태:** Round 2 유보 → Round 3 §Q1c에서 승격 채택 (T4에 포함).
- **Round 4 검증 경로:**
  1. T11 재제출 — **T4 DB 확장 only**로 먼저 테스트.
  2. 결과 분기:
     - **≥ 1 CANDIDATE 회수** → AtYUCCA6 복원 정당화, MVP 유지.
     - **여전히 0 CANDIDATE** → T6 canonical-triplet rule (bar + P-35S + T-ocs)로 보완 → 재검증.
     - **canonical triplet도 실패** → host YUC6 ortholog overlap 문제 → Q1c 복원 철회 + T9 pre-mask BED에 Gmax YUCCA family 강제 포함 → 재검증.
- **Fallback:** 3번째도 실패 시 soybean은 "WT control 없이는 exploratory only"로 승격 포기 (AC-1 목표 제외, AC-2만 측정).
- **Owner:** gmo_expert (판단) + won_yim (최종 PI 판정).

### §8.2 s04 minimap2 PoC 승인 여부 — T12 결과 대기

- **현 상태:** Round 3에서 rice_G281 PoC만 **조건부 승인**. PoC 결과로 배포 판정.
- **승격 기준 (전부 충족):**
  1. rice_G281 Chr3:16,439,674 → CANDIDATE 유지
  2. Phase 1 site 수 ±20% 이내 (BWA 대비)
  3. s04 wall time -30% 이상
  4. Soft-clip reads의 MAPQ < 20 비율이 BWA 대비 증가하지 않음 (BUG-7 회귀 방지)
- **통과 시:** v1.1 D6으로 승격 + tomato/cucumber/soybean/corn 추가 PoC는 D19.
- **실패 시:** s04 BWA 유지. haibao의 minimap2 주장은 v2 grant (V2)로 자동 이관.
- **Owner:** compute_eng (측정) + haibao (해석) + won_yim (최종 승인/기각).

### §8.3 BUG-7 microhomology MAPQ 0-19 회귀 방지 체크리스트

**배경:** BUG-7 (`4cec387` revert in `68d0e52`) — T-DNA microhomology junction이 MAPQ 0-19 read에 anchor되어 있어, MAPQ filter 추가 시 systemic FN. A2_3 ch01:91,002,744 회귀로 발견.

**MVP 구현 중 필수 체크:**

| # | 체크 항목 | 검증 방법 | 관련 Task |
|---|-----------|-----------|-----------|
| CL-1 | Site discovery 파라미터 변경 금지 (`cluster_window=50`, `min_clip=20`, `MIN_CLUSTER_DEPTH=3`) | T11 재검증 시 rice/A2_3 site count 비교 | T2, T3, T10 |
| CL-2 | `find_softclip_junctions` 내부 MAPQ filter 추가 금지 | code review + `grep "mapq" scripts/s05/site_discovery.py` | T7 분할 시 |
| CL-3 | T9 pre-mask BED는 site **drop** 아닌 **downgrade (`FALSE_NEGATIVE_MASKED` 태그)** 만 | `grep "FALSE_NEGATIVE_MASKED"` 존재 확인 | T10 |
| CL-4 | T12 s04 minimap2 PoC에서 MAPQ < 20 비율 BWA 대비 증가 시 기각 | PoC 보고서 | T12 |
| CL-5 | T5 transgene_db dedup이 univec의 microhomology-bearing seq 제거 시 element_db로 re-tag 검증 | 소스 tag manifest | T5 |
| CL-6 | T3 `max_rounds=3` 축소가 second junction (low-MAPQ side) recall에 영향 있는지 T2 측정 포함 | haibao T2 보고서 | T2, T3 |

**Owner:** bio_king (CL-1/2/3/5), compute_eng (CL-4), haibao (CL-6), won_yim (integration review).

**회귀 발생 시:** MVP release 보류, 해당 Task rollback + 재작업.

### §8.4 기타 Round 4 판정 대기

| # | 항목 | Round 4 trigger | 대안 |
|---|------|-----------------|------|
| UR-1 | T2 측정 결과 → `max_rounds ∈ {3, 5}` 최종 선택 | T2 보고서 (round-3 기여 %) | 기여 ≤10% → 3; 10-30% → 5; ≥30% → 8 유지 |
| UR-2 | T11 재검증 verdict 전수 → v1.0-rc 태그 / rollback 판정 | 금요일 integration | 한 샘플이라도 AC-1 회귀 시 rollback |
| UR-3 | SpCas9 DB 추가 후 tomato A2_* samples의 Phase 1.5 site count 폭증 여부 | T11 tomato 재검증 | 폭증 시 `--cas9-info-only` flag 추가 (gmo_expert R2 Q1 제안) |
| UR-4 | 임시 fan-out (T8) 제거 시점 | v1.0 → v1.1 전환 | D2 in-memory sharing + D3 full 분할 완료 시 CLI entrypoint로 교체 |

---

## §9. Round 4 실행 가이드

| 담당 | Day-by-day 주 작업 | 주 산출물 |
|------|---------------------|-----------|
| **bio_king** | Mon: T6/T7 설계 / Tue-Wed: T5/T6 구현 / Thu: T10 / Fri: review | PR 3-4개, pytest 5+ |
| **gmo_expert** | Mon: T4 설계 / Tue: DB 구축 / Wed: T9 rationale TSV / Thu: review | DB PR + rationale TSV |
| **compute_eng** | Mon: T1 / Tue: T1 + T8 설계 / Wed-Thu: T8 / Fri: T11 + T12 | PR 2-3개 + 재검증/PoC 보고서 |
| **haibao** | Mon: T2 측정 / Tue: T3 + T9 스크립트 / Wed: T9 integration / Fri: T12 협업 | 측정 보고 + script PR + PoC 해석 |
| **won_yim (PI)** | 일 1회 standup + 금요일 integration 검증 + v1.0-rc 태그 판정 | release approve or reject |

**일정:**
- 2026-04-20 (Mon) — Kickoff, T1/T6/T7 independent start
- 2026-04-21 (Tue) — T2/T3/T4 start
- 2026-04-22 (Wed) — T5/T9 integration
- 2026-04-23 (Thu) — T8/T10 integration
- 2026-04-24 (Fri) — T11 재검증 (wall 8-16h) + T12 PoC
- 2026-04-25 (Sat) — won_yim v1.0-rc 태그 판정 (PASS → release, FAIL → Round 5 재작업)

---

## §10. 서명 (합의 기록)

| Teammate | Round 1 기여 | Round 2 기여 | 주요 수용 | 주요 양보 |
|----------|---------------|---------------|-----------|-----------|
| **bio_king** | s05 refactor 10 target + TDD 10건 + smell 5건 | Q1 R2 수정 (run_blastn만 불충분 인정) + Q2 config-driven 동의 + Q3 분할-먼저 원칙 | 최소 4-way 분할 MVP | s05 전체 분할 Phase 2 |
| **gmo_expert** | element DB 15-seq MUST-HAVE + per-host FP matrix + CRL/SpCas9 gap | Q1 cd-hit 효과 미미 + Q2 pre-mask (100%-ortholog only) + Q3 2-tier 자진 제안 | AtYUCCA6/UGT72E3 MVP | CRL 82 seqs Phase 2 |
| **compute_eng** | 실측 sacct 표 + 병렬화 맵 + BUG-5/14/15/16 fix + CI smoke test | Q1 실측 ROI -5~7% (BLAST) + Q2 Phase 1 minimap2 반대 + Q3 DB 확장 +1% 수용 + s04 PoC 설계 | array 임시 fan-out MVP 필수 | s04 minimap2 full deploy Phase 2 |
| **haibao** | Level 1-3 알고리즘 재평가 + k-mer/synteny/long-read 비판 | Q1 Level 2 scope-out 수용 + Q2 P1/P2/P6 사수 + Q3 k-mer 측정 선행 자인 + Q4 MCscan W1/W2 분리 | Level 1 P1/P2/P6 MVP 3개 | Level 2/3 v1.1 ~ v2 grant |
| **won_yim (PI)** | 8 AC 정의 + 검증 샘플 rationale + regulatory MUST 8항 | Q1-Q5 scope 판정 + Round 3 synthesis | §3 MVP 12 task 확정 권한 | AC-3/AC-7/AC-8 MVP 완화 |

**본 합의문은 2026-04-16 Round 3 종료 시점의 team-lead 중재 하에 5명 전원 참여로 작성되었습니다. 본 문서의 결정은 Round 4 구현 시 판정 근거이며, MVP v1.0 release 후 학습 결과 기반으로 v1.1 범위 재조정이 가능합니다.**

---

**끝.**
