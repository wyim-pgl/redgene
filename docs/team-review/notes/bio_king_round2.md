# bio_king Round 2 — cross-talk 답변 (2026-04-16)

---

## Q1. haibao "DB 정규화 + BLAST 1회 + in-memory sharing" vs 내 `run_blastn` 공통 래퍼

겹치는 부분은 **호출 경로 통일**뿐이고 근본 처방은 다릅니다. 내 P0 #3은 29군데 `subprocess.run([...blastn...])`을 공통 래퍼로 치환해 BUG-2(stem collision)·옵션 불일치를 잡는 수준이라, BUG-3(source tie-break)·BUG-4(phase별 DB 누락) 같은 **결과 흐름 버그까지는 예방 못 합니다** — 그건 `_batch_check_element_hits`(s05:323), `classify_site_tiers`(s05:840), `annotate_insert`(s05:2824)가 **독립적으로 DB 리스트를 다시 구성**하는 구조 자체의 문제입니다(실제 BUG-4가 정확히 이 구조에서 터졌음).

따라서 haibao 방향(**DB 정규화 + 호출 1회 + 결과 공유**)에 동의하되 범위를 분리 제안:
- **수용:** cd-hit-est로 `element_db ∪ common_payload ∪ s04b/contigs`를 pre-merge + source tag(hierarchical `element_db > payload > sample_contig > univec`) 부여. query 1회, 결과를 `dict[site_id] → BlastHits`로 전달.
- **TDD 관점 강점:** `BlastHits`가 dataclass이면 `classify_site_tiers`·`annotate_insert`·`compute_verdict` 모두 가짜 `BlastHits` fixture로 독립 테스트 가능. 지금은 BLAST 파일 생성이 강제되어 unit test가 불가능(Round 1 smell #1과 같은 뿌리).
- **단, 순서:** 모듈 분할(P0 #1) + verdict 순수화(P0 #2) **먼저**, 그 다음 DB 정규화. 역순이면 in-memory sharing 경계가 어디서 끊기는지 알 수 없어 회귀 위험. BUG-4도 Phase 3(annotate)이 Phase 1.5(classify)의 extra_dbs 리스트를 잃어버려서 생겼으므로, 모듈 경계가 명확해진 뒤에 sharing 설계가 의미 있음.

**결론:** "내 래퍼만으로 충분" 아님. haibao 구조까지 가야 BUG-2/3/4 class 근본 차단. 단, 이행 순서는 **분할 → verdict 순수화 → DB 정규화 → sharing**.

---

## Q2. canonical-triplet rule: (a) `verdict.py` 하드코딩 vs (b) `config.yaml` 외부화

**(b) config 외부화에 동의.** 단, TDD 가능성을 잃지 않도록 loader를 분리.

이유:
- `config.yaml`은 이미 샘플별 `expected.insertion_chr/pos`, `construct_reference`, (계획된) `grna`를 담고 있어, 규제 대응용 event metadata의 자연스러운 귀속지. 새 event(예: thaumatin II 기반 cucumber)가 추가될 때마다 Python commit이 필요한 건 regulatory 작업 방식과 맞지 않음(won_yim 입장 타당).
- `compute_verdict()`는 여전히 pure 함수 유지 가능: 시그니처를 `compute_verdict(evidence: FilterEvidence, rules: VerdictRules)` 로 만들고, `VerdictRules`는 dataclass (`payload_whitelist: set[str]`, `canonical_triplets: list[frozenset[str]]`, threshold 상수). rule 로딩은 `load_verdict_rules(cfg: dict) -> VerdictRules`가 담당. 테스트는 dict fixture 주입으로 충분.
- 기본값(fallback): config에 `canonical_triplets` 키가 없으면 `s05/verdict.py`의 `DEFAULT_TRIPLETS = [frozenset({"bar", "P-CaMV35S", "T-ocs"}), ...]` 사용. 하드코딩 안전망 + config override는 양립 가능.

**결론:** (b) + loader 경계 + pure `compute_verdict(evidence, rules)`. won_yim 선호 수용, TDD 손실 없음.

---

## Q3. compute_eng의 `--site-id-only` / `--phase-4-only` CLI 순서

**즉시 시너지 쪽이지만, CLI-먼저는 반대 — 분할 먼저, CLI는 분할 완료 후.**

근거:
- 현재 `main()` (s05:3545-3805)은 Phase 1 → 1.5 → per-site loop → Phase 4 순으로 강하게 얽혀 있고(`host_endo_ids`, `construct_flanking`, `extra_dbs`가 phase간 공유), 중간 phase만 분리해 CLI에 노출하려면 **phase boundary가 함수 경계와 일치**해야 함. 지금 구조에서 `--phase-4-only` 추가하면 Phase 1.5 결과(`host_endo_ids`, `tier_results`)를 디스크에서 재로드하는 ad-hoc 직렬화 코드가 먼저 들어가고, 모듈 분할 시 그 코드를 또 뜯어야 함.
- 분할 계획에는 이미 phase별 모듈(`site_discovery.py`, `classify.py`, `assemble.py`, `annotate.py`, `report.py`)이 있고, 각자 `main()`에서 호출. `__main__.py`에 `--phase`, `--site-id` 플래그를 노출하면 자연스러움. 그 시점에 array job 진입점을 단순화할 수 있음.
- site 간 독립은 맞지만 Phase 1.5의 transgene-positive 필터(s05:1103)가 site 리스트를 **전역적으로** 결정하므로, array 병렬 전 Phase 1.5를 한 번 돌려 `assembly_sites.tsv`를 만들고 array task가 그걸 읽어 각자 site 1개 처리하는 fan-out 구조가 필요. 이 구조도 모듈 분할이 선결조건.

**결론:** 분할 → phase별 CLI entrypoint → array 병렬화. CLI를 먼저 넣어 band-aid 활성화는 soybean TIMEOUT 완화 효과는 있어도 나중에 뜯어내야 하므로 TDD상 권장 안 함. 다만 compute_eng가 **지금 당장 UGT72E3 resubmit을 살려야 한다면** `run_pipeline.py` 레벨에서 site list를 split해 독립 worktree로 N개 병렬 SLURM job을 쏘는 쪽이 s05 내부 변경 없이 가능 — 이 임시안은 별도 plan.
