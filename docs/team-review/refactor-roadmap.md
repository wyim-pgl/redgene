# RedGene Refactor Roadmap (2026-04-16)

**목적:** Korean quarantine GMO/LMO detection assay(AC)용 production-grade 파이프라인 완성.
**전제:** `main` @ `729835d`, 5/7 샘플 CANDIDATE 검증 PASS, pytest 10/10 PASS.
**설계 원칙:** TDD·증분 개선. 각 Phase 종료 시 (a) 새 테스트 green, (b) 기존 5개 CANDIDATE 회귀 없음, (c) bug.md에 결과 기록.

Phase 게이트:
- **Phase 1 (MVP, ≤2주 / ≤10 person-days):** regulatory-critical 결함 + per-site 속도 + TDD 인프라. 승인 기준 → `validation sample 5종 verdict 불변` + `audit trail fields 7개 모두 기록`.
- **Phase 2 (v1.1, 1-2개월):** 구조 개선 + 알고리즘 실험. 승인 기준 → `adaptive cluster_window A/B 검증` + `모듈 분할 완료 후 pytest 30+ green`.
- **Phase 3 (v2.0 grant scope, 6개월+):** 차세대 SV 알고리즘 + long-read ready. 별도 grant 이후.

---

## Phase 1 — MVP (≤ 10 person-days)

| # | Task | Owner | Effort (pd) | Risk | AC contribution | Dependency |
|---|------|-------|-------------|------|-----------------|------------|
| P1-1 | `compute_verdict()` pure 함수 분리 — `FilterEvidence` dataclass 입력, `(verdict, reason)` 출력. 현 `generate_report()` 318줄(s05:3197-3514)에서 BLAST I/O와 얽힌 tier 결정 로직 추출 | bio_king | 1.5 | Low — 기존 4-filter 로직을 그대로 옮김 | AC.4 (결정 근거 감사 가능), test가 verdict 자체를 커버 | 없음 |
| P1-2 | `run_blastn()` 공통 래퍼 — `s05` 12 blastn 호출 + `s04b:120` + `s06:74,93,143,473`를 단일 진입점으로. `task / query / subject_or_db / outfmt / opts / tag`. stem collision(BUG-2), 0-byte DB(BUG-15) 방어 포함 | bio_king | 1.5 | Low — 기계적 치환. tag로 filename 충돌 강제 차단 | AC.6 (도구 경로 통일 로그), BUG-2/15 재발 차단 | 없음 |
| P1-3 | Extra-DB 4-way source tag 정리 — `"element_db" / "payload" / "sample_contig" / "univec"`. `_should_replace`(s05:818)와 `classify_site_tiers` extra-DB 경로(s05:992) 동일 규칙. canonical_triplet matching은 tag별로 분기 | bio_king | 1.0 | Low — pytest 기존 `test_element_db_hit_wins_over_univec_at_tied_bitscore` 확장 | AC.3 (payload 추적), canonical triplet rule의 근거 | P1-1 |
| P1-4 | `_should_replace` element_db-family > univec 일관 적용 — family = {element_db, payload, sample_contig}. `annotate_insert`(s05:2768)·`classify_site_tiers`(s05:840) 양쪽에 동일 규칙 주입 | bio_king | 0.5 | Low | BUG-3 회귀 방어, rice Chr3:16,439,674 영구 보호 | P1-3 |
| P1-5 | canonical-triplet rule via `config.yaml` `VerdictRules` — Round 2 Q2 합의안. `DEFAULT_TRIPLETS = [{"bar","P-CaMV35S","T-ocs"}, ...]` fallback. `compute_verdict(ev, rules)` 시그니처 | bio_king | 1.0 | Low — config schema 확장만, 기존 샘플 config 변경 없음 | AC.2 (event-별 승격 근거), soybean 0-CAND 구제 | P1-1, gmo_expert H1 검토 |
| P1-6 | per-site CLI `--site-id-only <id>` / `--phase-4-only` — 최소 분리만. site list 파일을 먼저 쓰고 array task가 각자 site id로 재실행 | bio_king + compute_eng | 1.5 | Med — Phase 1.5 결과(`host_endo_ids`, `construct_flanking`)를 TSV로 dump/load해야 함. 디스크 재로드 포맷이 Phase 2 모듈 분할 시 바뀔 수 있음 — Round 2 Q3의 "band-aid" 주의 | UGT72E3/AtYUCCA6 24h 한계 돌파. array 진입점 제공 | P1-1, P1-2 |
| P1-7 | assembly max_rounds 8→3 (haibao L1-3) — `scripts/s05_insert_assembly.py:2372`의 `max_rounds=8` 기본값을 3으로. 수렴 히스토리(`growth_history`) 통계로 회귀 없음 증명 | haibao + bio_king 검증 | 0.5 | Med — cucumber line225(8 CANDIDATE)·A2_3 NODE_1+NODE_10 tiled 케이스에서 round 4-6에 도달한 site 있는지 s05_stats.txt 확인 필요 | per-site 14min→5min, soybean 24h 완화 | 검증용 rerun |
| P1-8 | 중복 BLAST 통합 (haibao L1-4 + 귀측 P0 #3) — `run_blastn` 래퍼 위에서 `_blast_insert_vs_host` 결과를 `_check_chimeric_assembly`(s05:3011)와 `_check_construct_host_coverage`가 공유. 현재 filename convention `_<stem>_vs_host_chrom.tsv`으로 일부 재사용 중이나 불완전 | bio_king | 1.0 | Low | per-site 3-filter 시간 1/3, 로그 일관 | P1-2 |
| P1-9 | `transgene_db` cd-hit-est -c 0.95 dedup (haibao L1-5) — build 시점(s05:868-881 `transgene_db` 조립부)에 `cd-hit-est -c 0.95 -n 10` 1회. 중복 엔트리(element_db ∩ univec) 제거 → BLAST space 축소 | gmo_expert + bio_king | 1.0 | Med — cd-hit-est 설치 필요. representative 선택에서 element_db 엔트리가 버려지지 않도록 `-g 1` + 이름 prefix로 순서 강제 | BLAST 5-10% 가속, AC 재현성(엔트리 수 명시) | makeblastdb pipeline, gmo_expert 승인 |
| P1-10 | SpCas9 + sgRNA scaffold DB 추가 (gmo_expert P0) — `element_db/cas9_scaffold.fa` 추가 (pSpCas9 backbone, tracrRNA, sgRNA scaffold ~100bp). CRISPR construct(A2 시리즈) 식별 보강 | gmo_expert | 0.5 | Low | CRISPR 계열 construct 자동 식별, AC.2 | 없음 |
| P1-11 | Audit trail R-1/R-2/R-3/R-4/R-7 통합 — reads/ref/DB SHA-256(R-1), git commit hash(R-2), DB md5 version(R-3), software manifest(R-4, micromamba env export), false-negative 명시(R-7, "no signal ≠ negative") 기록. `s05_stats.txt` 헤더 + `insertion_*_report.txt` footer | bio_king + compute_eng | 1.0 | Low — stat file 포맷 확장만 | AC.5 (규제 감사 대응), R-1~4 모두 | P1-2 (software 목록 수집) |
| **합계** | | | **11.0 pd** | | 5/7 CANDIDATE 유지 + soybean 3h 이하 + audit 7-field | |

**주의:** 11 pd로 목표 10 pd 초과. P1-7(assembly rounds 8→3)를 Phase 2로 이관하면 10.5 pd. 최종 조정은 team-lead 판단.

### Phase 1 Acceptance Criteria
- [ ] pytest 10 → 20+ green (Phase 1 테스트 strategy 참조)
- [ ] rice_G281 Chr3:16,439,674, A2_3 ch01:91,002,744, cucumber line212/224/225 primary site 모두 `CANDIDATE` 유지 (regression gate)
- [ ] soybean_UGT72E3 24h 내 완료 (per-site <10min)
- [ ] `insertion_*_report.txt`에 7-field audit trail 존재
- [ ] `compute_verdict()` 단독 호출 가능, BLAST 없이 verdict 산출

---

## Phase 2 — v1.1 (1-2개월, deferred)

| # | Task | Owner | Effort (pd) | Risk | AC contribution | Dependency |
|---|------|-------|-------------|------|-----------------|------------|
| P2-1 | s05 8-10 모듈 분할 — `scripts/s05/{models,io_fasta,blast,site_discovery,classify,reads,assemble,refine,annotate,filters,verdict,report,legacy}.py`. Round 1 §2 계획 | bio_king | 3.0 | Med — 대규모 import 재구성. Phase 1의 `run_blastn`·`compute_verdict` 분리가 pre-requisite라 Phase 2 착수 안전 | 유지보수성, 리뷰 가능성, TDD 확장 여지 | Phase 1 종료 |
| P2-2 | haibao L2 dual-anchor + PE discordancy site discovery — 현 soft-clip only(s05:373)에 PE read discordancy 추가. cluster_window adaptive 포함 | haibao + bio_king | 4.0 | High — site 수 증가로 FP 폭증 가능성. 기존 filter A/B/C/D 재조정 필요 | short/long T-DNA 모두 포괄, sensitivity↑ | P2-1 |
| P2-3 | haibao L1-1 cluster_window adaptive — coverage·read length 기반 동적 계산. BUG-7 회귀 방지 위해 fixture-based 결정적 테스트 필수 | haibao | 1.5 | Med — BUG-7 reversion 경험상 validation 샘플 5종 모두 통과해야 merge | coverage 5-15x 범위 강건성 | P2-1 |
| P2-4 | haibao L1-2 BLAST→LAST switch (평가) — clip-short alignment 정확도·속도 A/B. 채택 시 `run_blastn` 래퍼 대응 `run_last` 추가 | haibao | 2.0 | Med — LAST 설치 + fixture 결과 재측정. 실패 시 BLAST 유지 | 잠재적 per-site 20-30% 가속 | P2-1, P1-2 |
| P2-5 | Assembly step interface (Round 1 target #5) — `AssemblyStep` protocol로 kmer/mm2/Pilon/SSAKE 추상화. `assemble_insert()`(s05:2363-2647) 285줄 리팩토링 | bio_king | 2.0 | Med — 4 step growth 조기종료 조건 재설계 | 유지보수성, 새 assembler 추가 용이 | P2-1 |
| P2-6 | `extract_candidate_reads` tempfile 안전성 (Round 1 target #7) — 6 tempfile을 `TemporaryDirectory`로 격리 | bio_king | 0.5 | Low | 디스크 누수 방지 | P2-1 |
| P2-7 | R-5 PDF report — `report.py` 확장. jinja2 + weasyprint or matplotlib PDF. viz scripts(`plot_sample_summary.py`) 재사용 | bio_king + won_yim | 2.0 | Low | 규제 제출 포맷 | P2-1 |
| P2-8 | R-6 KCGP nomenclature — insertion 좌표·event 명명 규칙. `VerdictRules`의 `nomenclature_scheme` 키로 확장 | won_yim + gmo_expert | 1.0 | Low | AC.7 (공식 명칭) | P1-5 |
| P2-9 | R-8 chain of custody — sample 입수~보고서까지 단계별 hash chain. 각 step script의 `--lineage-file` 옵션 | bio_king + compute_eng | 1.5 | Med — step 간 통신 변경 | AC.5 확장 | P2-1 |
| **합계** | | | **17.5 pd** | | | |

### Phase 2 Acceptance Criteria
- [ ] pytest 20+ → 40+ green
- [ ] 모듈 분할 후 `from scripts.s05 import compute_verdict` 경로 동작
- [ ] Phase 1 validation 샘플 5종 + coverage subsample 9종 모두 regression 없음
- [ ] PDF report 1건 이상 sample-spot check 제출 가능

---

## Phase 3 — v2.0 grant scope (6개월+)

| # | Task | Owner | Effort (pd) | Risk | AC contribution | Dependency |
|---|------|-------|-------------|------|-----------------|------------|
| P3-1 | haibao L3 GRIDSS 기반 SV discovery — paired-end + split-read + assembly 통합. 현 find_softclip_junctions 대체 | haibao + 외부 bioinformatics 공동연구 | 15+ | High — assembly-based 접근과 완전 다른 알고리즘. 별도 validation 필요 | 대규모 SV, complex rearrangement 커버 | grant funding, Phase 2 종료 |
| P3-2 | Macrosynteny-aware site filter — host 근연종 synteny block 참조해 host-derived hit 제거. Filter B/D 보강 | haibao + gmo_expert | 8 | High — synteny DB 구축 + cultivar drift 대응 | host homology FP 재발 차단 | long-term |
| P3-3 | Long-read ready (ONT/PacBio HiFi) — s01~s05 pipeline을 long-read 입력 대응. `bwa` → `minimap2 -ax map-ont`, soft-clip 대신 long-read split | 전체 팀 | 20+ | High — 실험실 인프라 필요 | 1 kb 이상 complex T-DNA 해상 | grant |
| **합계** | | | **43+ pd** | | | |

### Phase 3 Acceptance Criteria (tentative)
- [ ] Illumina-only → Illumina + ONT hybrid 지원
- [ ] 5/7 샘플 + 신규 long-read 샘플 3건 end-to-end

---

## 전체 Dependency Graph

```
Phase 1 (≤ 10-11 pd, serial where noted):
  P1-2 (run_blastn) ─┬─> P1-8 (BLAST sharing) ─┐
                     ├─> P1-11 (audit trail) ──┤
                     └─> P1-6 (per-site CLI)   │
  P1-1 (compute_verdict) ─┬─> P1-3 (4-way tag) ─> P1-4 (merge rule) ─> P1-5 (triplet rule) ─┤
                          └─> P1-6 (per-site CLI)                                            │
  P1-7 (rounds 8→3) ─────── independent ────────────────────────────────────────────────────┤
  P1-9 (cd-hit-est) ─> makeblastdb build                                                     │
  P1-10 (Cas9 DB)   ─ independent                                                            │
                                                                                              ▼
                                                                                  Phase 1 Gate (AC)
                                                                                              │
Phase 2 (1-2 months):                                                                        ▼
  P2-1 (모듈 분할) ─┬─> P2-2 (dual-anchor) ─┐
                    ├─> P2-3 (adaptive window)
                    ├─> P2-4 (LAST)
                    ├─> P2-5 (assembly interface)
                    ├─> P2-6 (tempfile safety)
                    ├─> P2-7 (PDF report)
                    └─> P2-9 (chain of custody)
  P1-5 ─> P2-8 (KCGP nomenclature)
                                             ▼
                                    Phase 2 Gate (AC)
                                             │
Phase 3 (grant scope): P3-1/2/3
```

**Critical path (Phase 1):** P1-1 → P1-3 → P1-4 → P1-5 (≈ 4 pd). 병렬 가능 task는 P1-2, P1-7, P1-9, P1-10.

---

## Risk Register

| ID | Risk | Mitigation | Owner |
|----|------|------------|-------|
| R1 | P1-7 (rounds 8→3)이 cucumber line225 8 CANDIDATE 중 일부를 drop | s05_stats.txt에서 rounds 4-6에 수렴한 site 수 사전 집계. 그 수가 >0이면 기본값 5로 조정 | haibao + bio_king |
| R2 | P1-9 cd-hit-est가 element_db 엔트리를 univec representative로 대체 | `-g 1 -d 0 -n 10` + element_db prefix로 정렬해 항상 element_db가 rep 되도록 강제 | gmo_expert 리뷰 |
| R3 | P1-6 Phase 1.5 TSV dump 포맷이 P2-1 모듈 분할 시 바뀜 | dataclass → JSON schema로 버전 필드 포함. P2-1에서 loader만 교체 | bio_king |
| R4 | P1-5 VerdictRules config migration이 기존 7개 샘플 config 파일 깨뜨림 | 키 optional + fallback 제공, 모든 기존 config에 `grep -L canonical_triplets` → 미존재 시 default | bio_king |
| R5 | Phase 2 P2-2 dual-anchor가 FP 폭증 (cucumber line225 유사) | 별도 feature branch + worktree + coverage subsample 9종 전체 regression | haibao |
| R6 | Phase 2 모듈 분할 중 import-cycle (verdict ↔ filters ↔ report) | 분할 전 `pydeps`로 grep. cycle 발생 시 `models.py`에 공통 dataclass 격리 | bio_king |

---

## Phase 1 실행 순서 권고

1. **Day 1-2:** P1-2 (`run_blastn`) — 모든 후속 작업의 전제. 기계적 치환이라 TDD 안전.
2. **Day 3-4:** P1-1 (`compute_verdict`) + P1-11 (audit trail) 병렬. compute_verdict이 먼저 끝나면 P1-3/P1-4 착수 가능.
3. **Day 5:** P1-3 (source tag) + P1-4 (merge rule). 기존 `test_element_db_hit_wins_over_univec_at_tied_bitscore` 확장.
4. **Day 6:** P1-5 (VerdictRules config). soybean UGT72E3에서 triplet rule A/B 검증.
5. **Day 7:** P1-8 (BLAST sharing) + P1-10 (Cas9 DB).
6. **Day 8:** P1-9 (cd-hit-est) + P1-6 (per-site CLI) 착수.
7. **Day 9:** P1-6 완료, P1-7 (rounds 8→3) 검증.
8. **Day 10:** Regression suite 전체 실행, bug.md 업데이트, Phase 1 Gate 통과 판정.

Phase 1 종료 후 1-2주 관찰 구간 두고 Phase 2 착수.
