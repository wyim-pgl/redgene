# RedGene Test Strategy (Phase 1 Focus, 2026-04-16)

**목표:** Phase 1 리팩토링 기간 동안 5/7 CANDIDATE 회귀 0건 보장 + `compute_verdict` / BLAST 공유 / VerdictRules loader를 TDD-first로 구현.
**현 상태:** `pytest tests/ -v` → 10 PASS. `tests/test_extra_element_db.py` 5건 + `tests/test_s04b_construct_assembly.py` 5건.
**목표 상태:** Phase 1 종료 시 pytest 30+ green (10 기존 + 10 pure-function + 3 spike-in + nightly regression).

---

## 1. 현 10 pytest 분석

### 커버되는 영역
| 테스트 | 대상 | 핵심 시나리오 |
|--------|------|----------------|
| `test_batch_check_element_hits_reads_extra_db` | `_batch_check_element_hits` (s05:323) | extra_dbs 단일 DB 읽기 |
| `test_batch_check_element_hits_without_extras_ignores_extra_seq` | 동 | extra_dbs=None guard |
| `test_batch_check_element_hits_consults_multiple_dbs` | 동 | common_payload + s04b contig 병합 |
| `test_annotate_insert_uses_extra_dbs` | `annotate_insert` (s05:2824) | BUG-4 회귀 방어 |
| `test_element_db_hit_wins_over_univec_at_tied_bitscore` | `_should_replace` (s05:818) | BUG-3 회귀 방어 |
| `test_spades_failure_on_tiny_input_yields_empty_contigs` | `s04b` | SPAdes 실패 시 clean 종료 |
| `test_handles_empty_reads_gracefully` | `s04b` | `_is_empty_fastq` |
| `test_filter_contigs_keeps_marker_positive_only` | `_filter_contigs_by_markers` (s04b:88) | BUG-5 fix 핵심 |
| `test_filter_contigs_handles_missing_marker_db` | 동 | 결측 DB graceful |
| `test_no_filter_keeps_raw_contigs` | 동 | `--no-filter` CLI |

### 커버되지 않는 영역 (gap)
| 영역 | 파일:라인 | 심각도 | 비고 |
|------|-----------|--------|------|
| verdict 결정 전체 (4 filter + UNKNOWN reclass) | s05:3197-3514 | **Critical** | BUG-6 계열(heuristic promotion) 재발 차단 불가능 |
| site clustering (cluster_window / min_depth) | s05:434-456 | High | BUG-7 cluster_window widening regression 재발 가능 |
| consensus 산정 | s05:228-278 | Med | 낮은 depth 샘플 정확도 |
| assembly merge / extend convergence | s05:2326, 1604 | Med | 각 assembler 조기종료 조건 |
| s03 error paths | s03:60-80 | Med | BAM 누락 시 sys.exit(1) 로직 |
| s06 indel parser | s06:243-375 | High | BUG-1 rename error 경험, 파서 pure라 TDD 최적 |
| BLAST 공통 래퍼 `run_blastn` (Phase 1 신규) | — | Critical | 19군데 치환의 안정성 |
| 3 filters BLAST 결과 sharing | s05:3115, 2999, 3050 | Med | Phase 1 P1-8 도입 시 필수 |
| VerdictRules loader (Phase 1 신규) | — | Med | config → dataclass 파싱 |

현재 coverage 수치(approx, 라인 기반): **s05 ~4%, s04b ~60%, s06 0%, s03 0%**. Phase 1 목표 → s05의 verdict/classify/filters만이라도 60% 이상.

---

## 2. Phase 1 추가 pytest 10건 (pure-function TDD)

Round 1 §3의 제안 목록을 Phase 1 delivery 우선순위로 재정렬.

### 2.1 `test_compute_verdict.py` (7 scenarios, Critical)

**목표:** `compute_verdict(evidence: FilterEvidence, rules: VerdictRules) -> tuple[str, str]` 결정 트리 완전 커버.

```python
# 7 scenarios
1. CANDIDATE:
   - foreign element ≥ 1, host_fraction < 0.80, largest_gap ≥ 500,
     no flanking, off_target_chrs < 2, construct_frac < 0.25
2. FALSE_POSITIVE A (host-fraction + small gap):
   - host_fraction = 0.85, largest_gap = 300  → FP, reason mentions 85%/300bp
3. FALSE_POSITIVE B (construct-flanking overlap):
   - site_chr=Chr11, site_pos=8800, flanking=[("Chr11", 8758, 8958)]
4. FALSE_POSITIVE C (multi-locus chimeric):
   - off_target_chrs = [("Chr5", 200), ("Chr7", 150)]
5. FALSE_POSITIVE D (construct+host explain insert):
   - construct_frac = 0.35, host_fraction = 0.60, combined = 0.90
6. UNKNOWN → FALSE_POSITIVE (host-only):
   - elements=[], host_fraction=0.90, construct_frac=0.02
7. UNKNOWN preserved:
   - elements=[], host_fraction=0.50  → verdict="UNKNOWN"
8. (bonus) CANDIDATE_LOW_CONF rule via VerdictRules:
   - canonical triplet {bar, P-CaMV35S, T-ocs} 모두 `sample_contig` source에서 매칭
     + host_fraction < 0.70  → CANDIDATE (P1-5 validation)
```

**특징:**
- BLAST I/O 없음. `FilterEvidence` dataclass만 조립해 호출.
- 각 scenario 4-8 줄. 총 60-80줄 테스트 모듈.
- Phase 1 P1-1/P1-3/P1-4/P1-5 작업 시 red-green-refactor 가능.

### 2.2 `test_run_blastn.py` (래퍼 3 cases, Critical)

```python
1. blastn-short with -subject: BUG-2 대응. 두 번 호출 시 output file 충돌 없음 확인.
2. blastn with -db: makeblastdb 자동 생성 + 0-byte DB 재빌드 (BUG-15)
3. blastn megablast with tag kwarg: filename에 tag 포함 (stem collision 방지)
```

fixture: 10bp query + 20bp subject FASTA. 실제 blastn 호출 필요(marks: `integration`).

### 2.3 `test_build_consensus.py` (2 cases)

```python
1. 3 seqs [AAAACG, AAATCG, AAAACG] → "AAAACG" (majority vote)
2. direction="left" vs "right" 처리 (s05:228-278)
```

### 2.4 `test_check_merge.py` (2 cases)

```python
1. contig_5p = "AAAATCGATCG", contig_3p = "CGATCGTTTT", min_overlap=5
   → "AAAATCGATCGTTTT" (5bp overlap merge)
2. 불일치: contig_5p="AAA", contig_3p="GGG" → None
```

### 2.5 `test_merge_annotations.py` (2 cases)

`_merge_annotations` (s05:2768).
```python
1. local hit (bitscore=300) + remote hit (bitscore=250) at same q_range
   → local kept (element_db priority)
2. Two hits, 85% overlap → dominated 제거; 79% overlap → 둘 다 유지
```

### 2.6 `test_site_overlaps_flanking.py` (3 cases)

`_site_overlaps_flanking` (s05:2986).
```python
1. flanking=(Chr3, 100, 200), site=(Chr3, 150) slop=500 → True
2. flanking=(Chr3, 100, 200), site=(Chr5, 150) → False (chrom 다름)
3. flanking=(Chr3, 100, 200), site=(Chr3, 800) slop=500 → False (slop 초과)
```

### 2.7 `test_extract_indels_from_pileup.py` (4 cases, High — s06 BUG-1 재발 방지)

`_extract_indels_from_pileup` (s06:338).
```python
1. "+3ACG" → [("insertion", "ACG")]
2. "-2TT" → [("deletion", "TT")]
3. "^F.+3ACG$" → [("insertion", "ACG")] (^F, $ skip)
4. "*..+1A-1T" → [("insertion","A"),("deletion","T")]
```

### 2.8 `test_strand_aware_extender.py` (2 cases)

`StrandAwareSeedExtender` (s05:1538).
```python
1. seed="AAATCG", reads=["AAATCGGGG", "AAATCGGGG", "AAATCGGGG"]
   → extend to "AAATCGGGG"
2. branch point ratio < 0.7 → stop
```

### 2.9 `test_s03_error_paths.py` (2 cases, Med)

```python
1. BAM 경로 존재하지 않음 → SystemExit with code 1
2. hit_count == 0 → SystemExit with code 1, "No reads mapped" stderr
```

`CalledProcessError`가 아니라 `sys.exit(1)` 직접 호출이므로 `pytest.raises(SystemExit)` 사용.

### 2.10 `test_verdict_rules_loader.py` (3 cases, Med — P1-5용)

```python
1. config.yaml 에 canonical_triplets 키 없음 → DEFAULT_TRIPLETS fallback
2. YAML에 단일 triplet → VerdictRules.canonical_triplets == [frozenset(...)]
3. payload_whitelist + threshold overrides → 모든 값 반영
```

---

## 3. Synthetic Insertion Spike-in Generator

**왜 필요한가:** 현재 모든 pytest는 function 단위. 전체 파이프라인 회귀를 확인하려면 실샘플(>1h) 재실행이 필요 — CI에 부적합. 합성 스파이크-인으로 5분 내 E2E 회귀 테스트를 만든다.

### 설계

```
tests/fixtures/spike_in/
├── synthetic_host.fa        # host slice 100kb (rice Chr3 부분 발췌)
├── synthetic_host.fa.fai
├── construct.fa             # 6-kb T-DNA construct (bar + P-CaMV35S + T-nos + ...)
├── spike_in_R1.fq.gz        # wgsim/pirs 합성: host + 1 CANDIDATE + 1 FP site
├── spike_in_R2.fq.gz
├── expected_sites.tsv       # chr, pos, verdict, insert_size, element_list
└── README.md
```

### 생성 스크립트 `tests/fixtures/spike_in/generate.sh`

```bash
#!/bin/bash
set -euo pipefail
# 1. host slice
samtools faidx db/Osativa_323_v7.0.fa Chr3:16000000-16100000 > synthetic_host.fa
# 2. T-DNA construct assemble from element_db
cat element_db/bar.fa element_db/P-CaMV35S.fa element_db/T-nos.fa > construct.fa
# 3. Build insertion: host[0..50000] + construct + host[50000..100000]
python make_spike.py --host synthetic_host.fa --construct construct.fa \
    --insert-pos 50000 --output spike_in_genome.fa
# 4. Simulate reads (15x coverage, 150bp PE)
wgsim -N 10000 -1 150 -2 150 -r 0.001 -e 0.001 \
    spike_in_genome.fa spike_in_R1.fq spike_in_R2.fq
gzip spike_in_R1.fq spike_in_R2.fq
# 5. expected_sites.tsv: Chr3:50000 CANDIDATE, bar+P-CaMV35S+T-nos
echo -e "Chr3\t50000\tCANDIDATE\t6000\tbar,P-CaMV35S,T-nos" > expected_sites.tsv
```

### E2E 회귀 테스트 `tests/test_spike_in_e2e.py`

```python
@pytest.mark.slow  # --runslow 플래그 필요. 기본 CI에서는 skip
def test_spike_in_produces_expected_verdict(tmp_path):
    # 1. run s01-s05 on spike_in_R1/R2 with synthetic_host.fa
    subprocess.run(
        ["python", "run_pipeline.py",
         "--sample", "spike_in", "--steps", "1-5",
         "--threads", "4", "--no-remote-blast"],
        check=True,
    )
    # 2. assert Chr3:50000 (±20) in insertion_*_report.txt with Verdict: CANDIDATE
    reports = list(Path("results/spike_in/s05_insert_assembly").glob("insertion_*_report.txt"))
    verdicts = parse_verdicts(reports)
    assert any(v["chr"] == "Chr3" and abs(v["pos"] - 50000) <= 20
               and v["verdict"] == "CANDIDATE" for v in verdicts)
    # 3. assert annotated elements include {bar, P-CaMV35S, T-nos}
```

**예상 실행 시간:** 100 kb host + 15x coverage = ~5분 (laptop), ~3분 (HPC).

**추가 시나리오 (Phase 1 필수):**
1. `spike_in_candidate` — 1 CANDIDATE site
2. `spike_in_false_positive` — construct-flanking 영역에 site 생성 → Filter B 작동 확인
3. `spike_in_chimeric` — 2 host chromosome에서 read mixed → Filter C 작동 확인

---

## 4. Regression Gate (block rules)

**원칙:** 5종 known positive site를 CANDIDATE에서 제거하는 변경은 merge 차단.

### Gate 정의 `tests/regression/known_positives.tsv`

| sample | chr | pos | must_be_verdict | tolerance_bp | min_insert_len |
|--------|-----|-----|-----------------|--------------|----------------|
| rice_G281 | Chr3 | 16439674 | CANDIDATE | 100 | 1000 |
| tomato_Cas9_A2_3 | SLM_r2.0ch01 | 91002744 | CANDIDATE | 100 | 5000 |
| cucumber_line212 | LKUO03001392.1 | 2751687 | CANDIDATE | 100 | 1000 |
| cucumber_line224 | LKUO03001512.1 | 581328 | CANDIDATE | 100 | 1000 |
| cucumber_line225 | LKUO03001451.1 | 6501 | CANDIDATE | 500 | 1000 |

### `tests/regression/check_regression.py`

```python
def check_regression(results_dir: Path, gate_tsv: Path) -> list[str]:
    """Compare current verdict against known_positives.tsv.
    Return list of regression messages; empty list = PASS.
    """
    failures = []
    for row in read_tsv(gate_tsv):
        reports = glob(f"{results_dir}/{row['sample']}/s05_insert_assembly/insertion_*_report.txt")
        hit = find_site_within_tolerance(reports, row['chr'], int(row['pos']),
                                          int(row['tolerance_bp']))
        if hit is None:
            failures.append(f"{row['sample']} {row['chr']}:{row['pos']} not found")
        elif hit['verdict'] != row['must_be_verdict']:
            failures.append(f"{row['sample']} verdict changed "
                            f"{row['must_be_verdict']} → {hit['verdict']}")
    return failures
```

### 실행 모드
- **Local:** `python tests/regression/check_regression.py` (수동)
- **Nightly:** cron job 0 2 * * * — Phase 1 완료 샘플 full rerun, gate 실패 시 Slack/메일 알림
- **Pre-merge:** PR CI는 spike_in만 돌리고, full regression은 merge 직전 수동

**예외:** 회귀가 발생했으나 *의도된* verdict 변화인 경우(예: soybean에서 false CANDIDATE를 FP로 정정), PR 설명에 `[regression-intentional]` + bug.md 업데이트 + gate_tsv 수정.

---

## 5. CI Integration (GitHub Actions)

### 5.1 CI 단계

```yaml
# .github/workflows/ci.yml
name: RedGene CI
on: [pull_request, push]
jobs:
  unit-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: mamba-org/setup-micromamba@v1
        with:
          environment-file: environment.yml  # Phase 1에서 추가
      - run: pytest tests/ -v -m "not slow"  # 단위 테스트만 (30+ tests, <2 min)
      - run: python tests/regression/check_spike_in_fixture_integrity.py

  spike-in-e2e:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: mamba-org/setup-micromamba@v1
      - run: pytest tests/test_spike_in_e2e.py -v --runslow  # ~5 min
```

### 5.2 pre-commit 훅 (compute_eng BUG-16 연계)

```yaml
# .pre-commit-config.yaml — Phase 1 추가
repos:
  - repo: local
    hooks:
      - id: no-shell-true-in-new-code
        name: "Reject new subprocess shell=True in scripts/"
        entry: python tools/check_no_shell_true.py
        language: system
        files: ^scripts/.*\.py$
      - id: no-sbatch-env-pollution
        name: "Ensure sbatch CLI flags override env (BUG-16)"
        entry: python tools/check_sbatch_cli_flags.py
        language: system
        files: ^run_.*\.sh$
```

### 5.3 BUG-14/15 와 통합

- **BUG-14 (worktree `db/` symlink):** CI 환경에서 `db/` 자동 생성 + `tests/fixtures/db/` 심볼릭. `tests/conftest.py`에 fixture 추가.
- **BUG-15 (0-byte transgene_db):** P1-2 `run_blastn` + `ensure_blast_db` 이 원자적 rename 구현하므로 CI에서도 재현 안 됨 확인.

---

## 6. Regression Auto-Revert 정책

**원칙:** 회귀는 빨리, 말없이 되돌린다. 수정은 별도 PR에서.

### 트리거

1. nightly regression gate FAIL (known_positives.tsv 중 1건 이상)
2. spike-in E2E 테스트 FAIL in main branch

### 자동 액션 (bot 또는 수동 on-call)

```bash
# Step 1: 회귀 commit 식별 (gate FAIL → git bisect)
git bisect start HEAD <last-known-good>
git bisect run tests/regression/check_regression.sh

# Step 2: revert + PR 자동 생성
git revert <bad-commit> --no-edit
gh pr create --title "revert: regression in <sample>/<site>" \
    --body "Auto-revert: known_positive <sample>:<chr>:<pos> verdict CANDIDATE → <new>. \
Bisected to <commit>. See tests/regression/last_fail.log"
```

### 판단 매트릭스

| 상황 | 액션 |
|------|------|
| 회귀 + `[regression-intentional]` 태그 부재 | 즉시 revert, on-call에 알림 |
| 회귀 + `[regression-intentional]` 태그 존재 | revert 안 함, bug.md 확인, 게이트 TSV 업데이트 PR로 요청 |
| spike-in E2E만 실패 (known_positives 통과) | main 유지, 48h 내 원인 분석 후 revert 판단 |
| 모든 regression + unit 실패 | 즉시 revert, 전체 팀 broadcast |

### 히스토리 보존

`docs/regression_log.md` — 시간순 revert 기록, 재착수 PR 링크. Phase 1 종료 시점에 주간 통계 리포트.

---

## 요약

- **Phase 1 테스트 추가 분량:** 10 unit + 3 spike-in E2E = 13건, 예상 작성 공수 **2.5 person-days** (roadmap P1-1~P1-5와 병행).
- **TDD 순서:** compute_verdict 7 scenarios를 **가장 먼저** (P1-1 red-green-refactor). run_blastn 래퍼 3 cases를 그 다음 (P1-2). 나머지 pure-function은 해당 로직 수정할 때 함께.
- **Regression gate는 지금 당장 추가해야** — Phase 1 구현 중 회귀를 즉시 탐지하기 위함. `known_positives.tsv`는 코드 변경 없이 오늘도 생성 가능.
- **CI 구성 Phase 1 필수:** unit + spike-in E2E만. regression nightly는 cron 기반.
- **auto-revert:** nightly + spike-in 실패 시만. local PR은 수동 판단.

---

## Appendix A. compute_verdict 3 핵심 시나리오 — Round 4 v1.0 checklist 복사용

team-lead Round 4에서 won_yim §7 "v1.0 release checklist"의 "Regression test 10/10 PASS + new compute_verdict 3 test" 항목에 바로 삽입 가능한 pytest snippet. 시나리오 선정 기준 → 하나는 happy path(CANDIDATE), 하나는 가장 흔한 FP 유형(FP-A host-fraction), 하나는 regulatory-critical edge case(UNKNOWN 유지).

**전제 시그니처 (P1-1 구현 후):**
```python
from scripts.s05.verdict import compute_verdict, DEFAULT_RULES
from scripts.s05.filters import FilterEvidence
```

### Scenario 1 — CANDIDATE (happy path)

```python
def test_compute_verdict_candidate_with_foreign_element():
    """1+ foreign element, low host-fraction, large foreign gap → CANDIDATE."""
    ev = FilterEvidence(
        elements=[(100, 600, "bar", "+", "element_db"),
                  (700, 1500, "P-CaMV35S", "+", "element_db")],
        host_endo_ids=set(),
        host_fraction=0.10,
        host_bp=600,
        insert_len=6000,
        n_count=0,
        largest_gap=5400,          # ≥ INSERT_MIN_FOREIGN_GAP (500)
        flanking_hit="",
        off_target_chrs=[],
        is_chimeric=False,
        construct_frac=0.08,
        combined_frac=0.18,
    )
    verdict, reason = compute_verdict(ev, DEFAULT_RULES)
    assert verdict == "CANDIDATE"
    assert reason == ""
```

### Scenario 2 — FALSE_POSITIVE A (host-fraction + small gap)

```python
def test_compute_verdict_fp_a_host_fraction_with_small_gap():
    """host_fraction ≥ 0.80 AND largest_gap < 500 → FP-A (most common FP pattern).
    Historical reference: insertion_22966 (96% host, 461bp gap) — see s05:72-74.
    """
    ev = FilterEvidence(
        elements=[(50, 450, "I-actin_rice", "+", "element_db")],   # foreign-tagged
        host_endo_ids=set(),
        host_fraction=0.96,
        host_bp=4800,
        insert_len=5000,
        n_count=0,
        largest_gap=200,          # < INSERT_MIN_FOREIGN_GAP (500)
        flanking_hit="",
        off_target_chrs=[],
        is_chimeric=False,
        construct_frac=0.02,
        combined_frac=0.97,
    )
    verdict, reason = compute_verdict(ev, DEFAULT_RULES)
    assert verdict == "FALSE_POSITIVE"
    assert "96%" in reason or "0.96" in reason
    assert "200" in reason
```

### Scenario 3 — UNKNOWN 유지 (regulatory critical)

```python
def test_compute_verdict_unknown_preserved_when_ambiguous():
    """No element annotations, moderate host-fraction → UNKNOWN (not auto-FP).

    Regulatory rationale: BUG-6 (f78912b reverted in 02e1579) showed that
    promoting UNKNOWN to CANDIDATE_LOW_CONF by heuristics produced host-DNA
    false positives. The inverse — auto-demoting every UNKNOWN to FP — would
    also mask genuine novel insertions (no annotation yet in element_db).
    UNKNOWN must be preserved when host_fraction is NOT high enough to claim
    host-only (≥ UNKNOWN_HOST_MIN_FRACTION = 0.85).
    """
    ev = FilterEvidence(
        elements=[],                          # no annotations
        host_endo_ids=set(),
        host_fraction=0.50,                   # < UNKNOWN_HOST_MIN_FRACTION (0.85)
        host_bp=2500,
        insert_len=5000,
        n_count=0,
        largest_gap=2500,
        flanking_hit="",
        off_target_chrs=[],
        is_chimeric=False,
        construct_frac=0.10,                  # > UNKNOWN_MAX_CONSTRUCT_FRAC (0.05) anyway
        combined_frac=0.60,
    )
    verdict, reason = compute_verdict(ev, DEFAULT_RULES)
    assert verdict == "UNKNOWN"
    assert reason == "no element annotations"
```

**실행:**
```bash
pytest tests/test_compute_verdict.py::test_compute_verdict_candidate_with_foreign_element \
       tests/test_compute_verdict.py::test_compute_verdict_fp_a_host_fraction_with_small_gap \
       tests/test_compute_verdict.py::test_compute_verdict_unknown_preserved_when_ambiguous \
       -v
```

**v1.0 release criterion (proposed):**
- [ ] 위 3 시나리오 PASS
- [ ] 기존 regression test 10/10 PASS
- [ ] `tests/regression/known_positives.tsv` 5/5 CANDIDATE 유지

나머지 4 시나리오(FP-B construct-flanking, FP-C chimeric, FP-D construct+host, UNKNOWN→FP host-only)는 nice-to-have, v1.0 release에 블로커 아님.
