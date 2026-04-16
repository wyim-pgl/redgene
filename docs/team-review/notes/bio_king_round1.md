# bio_king Round 1 — 코드 품질 검토 (2026-04-16)

**검토 대상:** `scripts/s01_qc.py` ~ `scripts/s07_copynumber.py`, 특히 `scripts/s05_insert_assembly.py` (3806줄).
**관점:** DRY · YAGNI · Pythonic · TDD · 리소스 안전성.
**전제:** `main` @ `729835d`, pytest 10/10 PASS, 5/7 샘플 검증 통과(resume.md).

---

## 1. 상위 10개 리팩토링 타겟

| # | 이름 | 파일 : 라인 | 우선순위 | 난이도 | 예상 효과 | 리스크 |
|---|------|-------------|---------|-------|-----------|-------|
| 1 | s05 모듈 분할 (site / assemble / annotate / filter / report) | `scripts/s05_insert_assembly.py:1-3806` | P0 | L | 리뷰/테스트 가능성 확보, import-cycle 없이 단위 테스트 부착 가능 | 큰 diff — 반드시 plan+TDD로. 회귀 위험은 `main()` 통합에 집중 |
| 2 | `generate_report()` 의 verdict 로직 분리 → 순수 함수 `compute_verdict()` | `s05:3197-3514` | P0 | M | verdict 단독 pytest 가능 (현재 0건). Filter A/B/C/D tier 결정이 side-effect-free 해짐 | 중 — 현재 `verdict` 계산이 중간에 BLAST를 호출해 파일을 만들기 때문에, BLAST 호출부 주입 필요 |
| 3 | BLAST 호출 공통 래퍼 `run_blastn(query, db_or_subject, outfmt, opts)` 도입 | `s05` 내 19회, `s04b:88-156`, `s06:73-100,142-147,473-476` 중복 | P0 | M | 29군데 중복 제거. BUG-2 같은 stem-collision, `check=True` vs 없음, `stderr=DEVNULL` 혼재를 정리 | 낮음. 시그니처 설계만 잘 잡으면 기계적 치환 |
| 4 | `_should_replace` 규칙을 `classify_site_tiers`의 extra-DB 경로에도 동일하게 | `s05:964-1002` vs `s05:928-951` | P1 | S | Step 3b에서 extra DB 결과를 강제로 `"element_db"` 소스로 넣는 로직(l.992)이 미묘함 — extra DB가 univec 엔트리를 포함할 때 버그 여지. BUG-3 회귀 방어 확장 | 낮음 |
| 5 | `assemble_insert()` 내 4종 assembler(kmer/mm2/Pilon/SSAKE) 순차 호출을 인터페이스화 | `s05:2363-2647` (285줄 단일 함수) | P1 | L | 각 assembler를 `AssemblyStep` protocol로 쪼개면 round 로그/성장 측정/cleanup이 단일화됨. `_pilon_r{rnd}` 등의 tempdir 생성·정리가 한 곳으로 | 중 — 현재 growth 조기종료 조건(`kmer==mm2==pilon==ssake==0`)이 ad-hoc이라 설계 변경 필요 |
| 6 | subprocess `shell=True` 6곳 제거 (pipeline 치환) | `s05:1519-1524,1756-1760,1930-1933,1958-1963,2034-2038,2146-2150`; `s03b:162` | P1 | M | shell 주입 표면적 제거, 경로에 공백/특수문자 안전. `minimap2 … | samtools sort`는 Popen+stdin chain로 | 낮음. 각 chain이 단순 mapping→sort 패턴 |
| 7 | `extract_candidate_reads`의 tempfile 스텝(6개 `_site_{id}_*.bam`)을 `TemporaryDirectory` 사용 + 예외 안전 정리 | `s05:1309-1499` | P1 | M | 현재 `check=True`로 중간 실패 시 tmp BAM이 남음. 대규모 sample 수천 site에서 누적. 16G 이상 디스크 여유에서 안전해지지만 soybean 사례(TIMEOUT)에서 눈에 띌 수 있음 | 낮음 |
| 8 | `main()` 내 per-site 순차 루프를 함수로 분리 + concurrent.futures 파라미터화 | `s05:3660-3728` | P2 | M | UGT72E3 24h TIMEOUT 문제의 구조적 해결 여지(site 간 독립). Round 2에서 compute_eng와 정책 합의 후 | 중 — BLAST/Pilon이 이미 내부 병렬이라 oversubscription 조심 |
| 9 | host BLAST DB 존재 체크 중복(3곳) 헬퍼화 `ensure_blast_db(fasta)` | `s05:727-745,1011-1015,883-889`; `s06:71-78` | P2 | S | 원자성(BUG-15: 0-byte DB) 결함 한 지점에서 수정 가능. `.tmp → mv` 패턴 부여 | 낮음 |
| 10 | `_extract_indels_from_pileup` / `_parse_pileup_indels` 유닛 테스트 부재 | `s06:243-375` | P2 | S | 현재 s06는 테스트가 0건. indel 파서는 pure string parsing이라 TDD에 이상적 (bug.md BUG-1 재발 방지) | 없음 |

---

## 2. s05 모듈화 제안

현재 `s05_insert_assembly.py` (3806줄)은 단일 파일에 전 파이프라인 단계가 섞여 있어 코드 리뷰·단위 테스트·변경 격리 모두 어렵다. 다음 구조로 분할 제안:

```
scripts/s05/
├── __init__.py
├── __main__.py              # 기존 main(), argparse, phase 오케스트레이션 (~120줄)
├── models.py                # JunctionCluster, InsertionSite, LegacyJunction, TierResult (~80줄)
├── io_fasta.py              # read_fasta, write_fasta, _read_fq_seqs, revcomp (~40줄)
├── blast.py                 # BLAST 공통 래퍼 (run_blastn, run_blastn_short, run_megablast, ensure_blast_db)
│                            # → _run_local_blast, _run_remote_blast, _batch_check_*, host_endo BLAST 모두 위임
├── site_discovery.py        # find_softclip_junctions, _build_consensus, _batch_check_maps_to_host (~350줄)
├── classify.py              # classify_site_tiers, _filter_host_endogenous, _should_replace,
│                            # write_tier_classification (~380줄) — Phase 1.5
├── reads.py                 # extract_candidate_reads, extract_unmapped_paired, _extract_seeds_at_positions (~350줄)
├── assemble.py              # StrandAwareSeedExtender, recruit_by_kmer, pilon_fill,
│                            # _minimap2_extend, _ssake_extend, _check_merge, _vote_extension (~700줄)
├── refine.py                # refine_with_foreign_reads, extract_foreign_reads, check_host_termination (~350줄)
├── annotate.py              # annotate_insert, _merge_annotations, _parse_blast6 (~200줄)
├── filters.py               # _blast_insert_vs_host, _check_chimeric_assembly,
│                            # _check_construct_host_coverage, _find_construct_flanking_regions,
│                            # _site_overlaps_flanking, 임계값 상수 (~300줄)
├── verdict.py               # compute_verdict() — pure 함수, BLAST 결과를 입력으로 받음 (~150줄)
├── report.py                # generate_report, write_stats (compute_verdict 호출) (~250줄)
└── legacy.py                # parse_legacy_junctions, legacy_junctions_to_sites (~100줄)
```

**설계 원칙:**
- `blast.py`는 "subprocess가 있는 유일한 통로"로 만든다. 나머지 모듈은 순수 파이썬 함수로 테스트 가능.
- `verdict.py`의 `compute_verdict()`는 입력(host_bp, largest_gap, construct_frac, off_target_chrs, elements, host_endo_ids)을 dataclass로 받아 `(verdict, reason)`을 리턴 — 현재 `generate_report()` 안에 묻혀 테스트 불가능.
- Phase 간 통신은 dataclass(`InsertionSite`, `AssemblyResult`, `AnnotationResult`, `VerdictInput`)로. dict는 금지.
- 전환 방식: 먼저 새 모듈을 만들고 기존 파일이 `from .s05 import *`로 래핑하는 compat shim 단계 → 테스트 통과 확인 → 기존 파일 제거. 1 commit 당 1 모듈 분리 권장.

**예상 PR 수:** 8-10회 (각 ≤400줄 diff). Round 2에서 plan 확정 필요.

---

## 3. 테스트 커버리지 gap

현재 테스트는 `_batch_check_element_hits`, `annotate_insert`의 extra_dbs, `_should_replace` tie-break, s04b 필터에만 몰려 있다. s05 핵심 3000줄 중 verdict·assembly·filter 경로는 단위 테스트 없음.

**최소 추가 TDD 목록:**

1. **`compute_verdict`** (모듈화 후 즉시 가능)
   - CANDIDATE: 1+ foreign element + host_fraction < 0.8
   - FALSE_POSITIVE A: host_fraction ≥ 0.8 AND largest_gap < 500
   - FALSE_POSITIVE B: site 좌표가 construct-flanking 슬롭 내
   - FALSE_POSITIVE C: off-target chr ≥ 2
   - FALSE_POSITIVE D: construct_frac ≥ 0.25 AND combined ≥ 0.85
   - UNKNOWN → FP: elements=[], host_fraction ≥ 0.85, construct_frac ≤ 0.05
   - UNKNOWN 유지: 위 조건 모두 불충족

2. **`_build_consensus`** (`s05:228`) — 다수결 투표, direction별 길이 처리, tie 처리

3. **`_check_merge`** (`s05:2326`) — min_overlap 경계, 불일치 시 None, 정확한 병합 결과

4. **`_extract_indels_from_pileup`** (`s06:338`) — `+3ACG`, `-2TT`, `^F.`, `$`, `*` 혼합 문자열

5. **`find_softclip_junctions` 클러스터링** — fixture BAM(tests/fixtures에 mini BAM 추가)으로 paired / single-direction / no-cluster 3케이스

6. **`_merge_annotations` greedy interval selection** — local/remote 80% overlap 규칙, bitscore tie에서 element_db 우선

7. **`_site_overlaps_flanking`** — slop=500 경계 내외 각 2건

8. **`StrandAwareSeedExtender`** — 간단한 contig+reads로 우측 extend 1개, branch point 정지 1개

9. **`s03` extract_reads error paths** — BAM missing, zero hits (현재 `sys.exit(1)` 경로 테스트 없음)

10. **`s06 call_variants_at_sites` WT subtraction** — mocked pileup으로 WT∩treatment 제거 확인

위 10건 중 1,2,3,4,7은 pure-function이라 fixture 없이 즉시 추가 가능. 가장 투자대비 효과가 높다.

---

## 4. 구체적 코드 smell 5개

### smell #1 — `generate_report()` 안에서 BLAST 호출 4번 + 순수 verdict 계산이 섞임
**위치:** `scripts/s05_insert_assembly.py:3197-3514` (318줄 단일 함수)

**문제:** verdict 결정이 Filter A (host BLAST) → Filter B (flanking) → Filter C (chimeric BLAST 재사용) → Filter D (construct BLAST) → UNKNOWN reclass BLAST 순으로 얽혀 있어, verdict 로직만 단독으로 테스트하거나 re-run하려면 전체 BLAST를 다시 돌려야 한다. 라인 3202의 TODO 주석 "If a verdict-only mode is needed later, extract into compute_verdict()"가 이미 그 필요를 인정함.

**개선 스케치:**
```python
# filters.py
@dataclass
class FilterEvidence:
    host_fraction: float
    host_bp: int
    largest_gap: int
    flanking_hit: str
    off_target_chrs: list[tuple[str, int]]
    is_chimeric: bool
    construct_frac: float
    combined_frac: float
    elements: list[tuple[int, int, str, str, str]]
    host_endo_ids: set[str]

def collect_filter_evidence(...) -> FilterEvidence:
    """All BLAST/IO happens here."""

# verdict.py
def compute_verdict(ev: FilterEvidence) -> tuple[str, str]:
    """Pure — no subprocess, no file IO."""
    ...

# report.py
def generate_report(..., evidence: FilterEvidence, verdict: str, reason: str):
    """Only formats."""
```

---

### smell #2 — subprocess `shell=True` 6건 (`bwa ... | samtools sort`)
**위치 (대표):** `scripts/s05_insert_assembly.py:1756-1760`
```python
subprocess.run(
    f"minimap2 -ax sr -t {threads} --secondary=no "
    f"{ref_fa} {r1_fq} {r2_fq} 2>/dev/null "
    f"| samtools sort -@ {threads} -o {bam}",
    shell=True, check=True,
)
```

**문제:** `ref_fa`/`r1_fq`/`r2_fq`가 Path 변수라 지금은 안전하지만 config.yaml-derived sample 이름이 경로에 섞이면 injection 표면. 또한 `2>/dev/null`로 stderr를 강제 버려 BUG-15(makeblastdb 0-byte DB)같은 실패 원인 추적이 어려움.

**개선:**
```python
p1 = subprocess.Popen(
    ["minimap2", "-ax", "sr", "-t", str(threads), "--secondary=no",
     str(ref_fa), str(r1_fq), str(r2_fq)],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
)
p2 = subprocess.run(
    ["samtools", "sort", "-@", str(threads), "-o", str(bam)],
    stdin=p1.stdout, check=True,
)
p1.stdout.close()
p1_err = p1.stderr.read().decode()
rc = p1.wait()
if rc:
    raise RuntimeError(f"minimap2 failed rc={rc}: {p1_err}")
```

---

### smell #3 — BLAST DB 존재 체크 3중 중복, 0-byte DB 방어 없음
**위치:** `s05:727-745`, `s05:883-889`, `s05:1011-1015`

각 지점마다 `host_ref.with_suffix(".fa.ndb").exists() or Path(str(host_ref) + ".ndb").exists()` 두 경로를 OR로 확인하는 동일한 블록이 복사되어 있다. BUG-15(db/transgene_db.fa 0-byte) 재발 방지 수단도 없음.

**개선 (`blast.py` 신규):**
```python
def ensure_blast_db(fasta: Path, *, force: bool = False) -> None:
    """Ensure a blastn-usable DB for `fasta` exists.

    Writes to <fasta>.{nhr,nin,nsq,ndb,...}. Atomic: builds into
    temp prefix first and only renames on success. Treats 0-byte
    existing .nhr as stale and rebuilds.
    """
    ndb = Path(str(fasta) + ".ndb")
    nhr = Path(str(fasta) + ".nhr")
    valid = (ndb.exists() and nhr.exists()
             and nhr.stat().st_size > 0)
    if valid and not force:
        return
    tmp_prefix = fasta.parent / f".{fasta.name}.tmp"
    subprocess.run(
        ["makeblastdb", "-in", str(fasta), "-dbtype", "nucl",
         "-out", str(tmp_prefix)],
        check=True, stderr=subprocess.PIPE,
    )
    for ext in (".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto", ".njs"):
        src = Path(f"{tmp_prefix}{ext}")
        if src.exists():
            src.rename(Path(f"{fasta}{ext}"))
```

---

### smell #4 — `extract_candidate_reads` 가 6개 tempfile을 만들고 마지막에만 cleanup
**위치:** `s05:1309-1499`

`region_bam`, `namelist`, `both_mates_bam`, `unmapped_bam`, `merged_bam`, `nsort_bam` 6개 파일을 step_dir에 나란히 만들고, 6개 subprocess 모두 `check=True`이라 중간 실패 시 앞서 쓴 BAM이 남는다. soybean처럼 1000+ site를 도는 sample에서 실패가 누적되면 step_dir가 수십 GB까지 증가.

**개선:**
```python
with tempfile.TemporaryDirectory(prefix=f"{site.site_id}_", dir=tmp_dir) as td:
    region_bam = Path(td) / "regions.bam"
    ...
    # 모든 중간 파일을 TemporaryDirectory 내부에 둠
    # 성공 시 out_r1/out_r2만 `tmp_dir` 밖으로 commit
```
BUG-11/12와 같은 계열. s04b에서 이미 `shutil.rmtree` pattern을 수정했으므로 동일 적용.

---

### smell #5 — `classify_site_tiers`의 extra-DB 경로가 source="element_db"로 강제
**위치:** `s05:992-999`
```python
if _should_replace(hits.get(qname), "element_db", bitscore):
    hits[qname] = {
        "element": cols[1],
        ...
        "source": "element_db",
    }
```

**문제:** `extra_transgene_dbs`에는 `common_payload.fa`(transgene 마커)도 들어가고 `s04b/contigs.fasta`(SPAdes 결과)도 들어간다. 둘 다 "element_db"로 태깅되어 downstream(`is_positive = hit.source == "element_db"`) 조건을 통과하지만, 실은 서로 다른 의미의 증거다. Round 2에서 gmo_expert가 "element_db / payload / sample_contig" 3-way source 분리를 요구할 여지가 크다. 최소한 `"source": "element_db_extra"`로 구분하고 `is_positive`는 그대로 통과시키면 추후 tier 세분화가 쉬워짐.

---

## 5. 다른 teammate에게 묻고 싶은 것

### haibao (알고리즘)
1. `find_softclip_junctions`의 `cluster_window=50`, `min_clip=20`, `MIN_CLUSTER_DEPTH=3`은 4x-15x coverage 전 범위에 고정값이다. BUG-7에서 300으로 넓히는 시도가 A2_3 회귀를 유발했는데, 샘플 coverage에 따라 자동 조정하는 알고리즘이 바람직한가 vs 고정값이 재현성 측면에서 나은가?
2. `assemble_insert`는 k-mer / mm2 / Pilon / SSAKE 4 assembler를 매 라운드 돌린다. minimap2 dual-anchor 기반 접근(= SPAdes 없이 read pair 양 끝 anchor를 host + insert 양쪽으로 직접 붙이기)이 대안으로 언급될 가능성이 높다. 현재 방식 vs dual-anchor 중 TDD·테스트 가능성 관점에서는 전자가 강하지만, 정확도 근거가 있는지 알고리즘 전문가 의견 필요.
3. `_merge_annotations`의 80% reciprocal overlap threshold는 어디서 나왔는가? element_db 내 동종 paralog(35S variant 다수)를 어떻게 처리할지 결정이 필요.

### gmo_expert (domain)
1. `CONSTRUCT_MIN_FRACTION = 0.25` (Filter D)는 "real T-DNA = 10% construct, FP = 50%" 경험치 기반이다. 샘플별(T-DNA 길이, 다중 copy) 분포를 실측해 threshold 보정 가능한지?
2. `element_db`(131종) + `common_payload.fa`(9종) 이중 구조에서 새 마커(예: thaumatin II 이미 추가됨) 추가 시 regression을 막는 버전 관리 규약이 없다 — git tag? checksum? 
3. BUG-3 fix(`element_db > univec`) 이후 tie-break 규칙이 도메인적으로 맞는가? 역으로 element_db 엔트리가 spurious hit일 때 univec이 이겼어야 하는 경우가 있을 수 있음.

### compute_eng (성능)
1. `assemble_insert` per-site ~14 min × 100 sites = 24h가 soybean TIMEOUT의 직접 원인. site 간 독립이므로 `ProcessPoolExecutor`로 N-way 병렬 가능. 단, 내부에서 `blastn -num_threads`, `minimap2 -t`, `samtools sort -@`을 이미 쓰므로 global concurrency 예산 설계 필요 (16 CPU를 site-parallel × tool-thread 로 어떻게 쪼갤지).
2. `_batch_check_element_hits`는 clip마다 blastn을 DB 수만큼 호출한다(`s05:352`). DB를 cat해서 한 번에 BLAST할 때와 per-DB 병렬(xargs -P / multiprocessing) 중 어느 쪽이 soybean-scale에서 빠른가?
3. `extract_unmapped_paired`가 step_dir에 caching(`s05:1502`). site 수천 개 중복 호출 시 디스크 vs 메모리 trade-off 점검 필요.

---

## 요약

- **P0 (blocker급 리팩토링):** s05 모듈 분할 → verdict 순수화 → BLAST 공통 래퍼 3건. 이 세 가지가 이후 모든 변경의 pre-requisite.
- **P1:** extra-DB source 태깅, assembler 인터페이스화, shell=True 제거, tempfile 안전성. 각각 S/M 난이도.
- **P2:** per-site 병렬화(Round 2 compute_eng 협의), BLAST DB helper, s06 파서 테스트.
- **TDD gap:** verdict 계산, consensus/merge pure 함수, pileup 파서는 fixture 없이 즉시 추가 가능 — 가장 먼저 부착해야 모듈 분할 중 회귀를 잡을 수 있다.
- **위험 포인트:** Round 2에서 haibao가 assembly 알고리즘 재설계를 주장할 가능성. 나는 증분 개선(P0→P1→P2 순서, 각 단계 pytest 녹색) 입장이며, 대안 알고리즘 채택은 기존 5/7 PASS 회귀 없음을 증명한 뒤에.
