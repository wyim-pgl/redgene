# Round 1 — compute_eng: HPC/SLURM 자원 최적화 분석

**Agent:** compute_eng · **Team:** redgene-review-2026-04-16 · **Date:** 2026-04-16
**Scope:** wall-time 병목 (s04 5-7h, s05 per-site 3-30min × 100+ sites), BUG-5/14/15/16 영구 fix, CI smoke test 설계

---

## 1. 샘플 크기별 실측 리소스 표

> 원천: `sacct -j 5628134,5628135,5628138,5628141,5628233`, `results/rg_*_<jobid>.err/.out` 페이즈 마커 grep, `ls -l results/<sample>/s04_host_map/<sample>_host.bam`

| Sample | Host | Host size | s04 BAM | s03 R1 construct | s04b raw/filt | Phase 1 sites | Phase 1.5 pos | Phase 4 verdicts | JobID | ReqMem | ReqCPU | Elapsed | MaxRSS | State | Per-site s05 | 추천 `--mem` / `--time` / `--cpus` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| rice_G281 | Osativa v7 | 374 Mbp | 6.3 GB | 448 KB | 122/6 | 33,411 | 21 | 3C/9F/9U | 5628134 | 32G | 8 | **1:15:20** | 7.7 GB | COMPLETED | **~3.4 min** | 32G / 6h / 8 |
| tomato_Cas9_A2_3 | SLM_r2.0 | 833 Mbp | 4.2 GB | 256 KB | 422/2 | 905 | 2 | 2C | 5628135 | 48G | 8 | **0:29:38** | 11.3 GB | COMPLETED | ~4 min | 32G / 2h / 8 |
| cucumber_line212 | B10v3 | 332 Mbp | 5.5 GB | 128 KB | 75/2 | — | 31 | 1C/4F/26U | 5628138_0 | 96G | 16 | **13:29:11** | 79.6 GB | COMPLETED | **~26 min** ⚠ | 96G / 18h / 16 |
| cucumber_line224 | B10v3 | 332 Mbp | 5.5 GB | 128 KB | 73/2 | — | 18 | 1C/2F/15U | 5628138_1 | 96G | 16 | **9:21:07** | 87.2 GB | COMPLETED | ~31 min ⚠ | 96G / 14h / 16 |
| cucumber_line225 | B10v3 | 332 Mbp | 5.5 GB | 192 KB | 123/2 | — | 38 | 8C/4F/26U | 5628138_2 | 96G | 16 | **13:42:22** | 59.6 GB | COMPLETED | ~21 min | 96G / 18h / 16 |
| soybean_AtYUCCA6 | Gmax v4 | 1.1 Gbp | 9.7 GB | 1.2 MB | 1,345/6 | — | 100 | 진행중 (78 processed) | 5628233 | 48G | 8 | **21:02 (running)** | 17.2 GB | RUNNING | **~16 min** | 64G / 36h / 16 |
| soybean_UGT72E3 | Gmax v4 | 1.1 Gbp | 9.6 GB | 1.3 MB | 1,641/6 | 28,256 | 105 | partial (48/53 processed) | 5628141 | 96G | 16 | **24:00:12** | 56.4 GB | **TIMEOUT** | **~28 min** ⚠⚠ | 128G / 48h / 16 또는 array 분할 |

**핵심 관측:**

- **MaxRSS vs ReqMem**: cucumber 계열이 96G 요청에 최대 **87 GB 사용**. 소이빈은 96G에 **56 GB / 48G에 17 GB**. 메모리는 그럭저럭 맞지만 cucumber line224는 여유 9%만 남아 **OOM 재발 위험**.
- **Per-site s05 time**는 host size보다 **Phase 1.5 positive count × host BAM size에 선형**:
  - rice 21 sites × (6.3 GB BAM) = 3.4 min/site
  - cucumber ~30 sites × (5.5 GB BAM) = 21-31 min/site
  - soybean 100+ sites × (9.7 GB BAM) = 16-28 min/site
- **s04 (host mapping)은 이번 검증 세트에서 재실행되지 않음** (`--steps 4b,5`만). s04 BAM은 이전에 생성되어 재사용. `bug.md` / `CLAUDE.md`에 기록된 5-7h 수치는 initial run 기준.
- **cpu-s2-core-0 파티션**만 사용 (사용자가 CLI로 명시 → BUG-16 회피). `SBATCH_PARTITION=cpu-s1-pgl-0` 환경변수가 여전히 세팅되어 있을 가능성.

---

## 2. 병목 분해

> 측정: `.err` 파일의 페이즈 마커 라인번호 + 파일 크기 + 프로세싱 완료 사이트수.

### s04 host mapping (5-7h, initial run)
- `scripts/s04_host_map.py:57-127`: `bwa mem -t N | samtools sort -o` 파이프 (직렬 한 스트림)
- `scripts/s04_host_map.py:149-183`: `samtools depth -a | python for-loop stream` → **대형 genome에서 실질적 병목 (soybean 1.1 Gbp 전 위치 스트리밍, 단일 스레드 Python 파싱)**
- 추정 분해(soybean 기준 7h):
  - bwa mem: ~55% (3.5-4h)
  - samtools sort: ~25% (1.5-2h, `-@ N` 멀티스레드)
  - samtools depth + Python stats: ~15% (~1h) — 이 부분은 Python side가 단일 코어 병목
  - bwa index / flagstat / index BAM: ~5%

### s05 insert assembly
**Per-site wall time 분해 (UGT72E3 28 min 기준, bio_king 리뷰와 교차검증 대상):**
- Phase 1 (scan host BAM soft-clips): **시간표상 전 site에 1회, 28,256 사이트 생성 → ~20min/total** (not per-site)
- Phase 1.5 (transgene-positive classification): **전 site에 1회 batch BLAST**, UGT72E3에서 `Transgene-positive (assemble): 105` 산출까지 ~15min
- Phase 2 (`extract_candidate_reads`): samtools region query, per-site ~30s–1min
- Phase 3 (`assemble_insert` k-mer 확장 + Pilon 갭필): **per-site ~15-25min — 전체 wall time의 90%**
  - `StrandAwareSeedExtender.load_paired_reads`: s03 1.3 MB (compressed) + junction-region 수백 KB 로드
  - `recruit_by_kmer` from unmapped pool: 매 라운드
  - up to 15 라운드 × (k-mer 확장 + minimap2 soft-clip 확장 + Pilon)
- Phase 3b (foreign read refinement): **per-site ~1-2min**
- Phase 4 (batch annotation): 전 site에 1회, BLAST + report 생성 ~3min

**주요 발견:**
- **s05의 per-site loop (line 3662-3728)는 완전 독립**: input은 `site`, `cand_r1/r2`, `host_bam` (read-only). output은 `<site>_insert.fasta`. 사이트간 공유 상태 없음 → **embarrassingly parallel**.
- Phase 4 annotation은 per-site loop 이후에 `combined_insert.fasta`로 일괄 수행. 이건 직렬 유지 가능 (BLAST 1회면 충분).
- **UGT72E3가 TIMEOUT 근본 원인**: 105 × 14min = 24.5h > 24h 한계. per-site 평균이 bug.md 추정(~10min)보다 더 긴 28min. 소이빈 host BAM이 큰 것과 관련.

---

## 3. 병렬화 가능 지점 맵

### A. s05 per-site array job (추천 — 가장 cost-effective)
- **구조:** `run_pipeline.py --steps 4b,5`로 Phase 1/1.5까지만 실행 후, Phase 2-3을 SLURM array로 분리. Phase 4는 array 완료 후 direct call.
- **코드 변경:**
  1. `scripts/s05_insert_assembly.py`에 `--site-id-only <id>` 옵션 추가 → 해당 사이트의 Phase 2+3만 실행.
  2. 새 wrapper `run_s05_array.sh`: Phase 1.5가 생성한 `site_tier_classification.tsv`에서 positive site 리스트 읽어 `#SBATCH --array=0-N` 동적 제출.
  3. `scripts/s05_insert_assembly.py`에 `--phase-4-only` 옵션 추가 → 모든 `<site>_insert.fasta`가 있으면 Phase 4만 실행.
- **wall time 예상:** 100 sites × 16min = 1600 min 직렬 → array 25 concurrent × 4 waves × 16min = **~64min 총 queue time**. UGT72E3 24h → **~3-4h** 가능.
- **Cluster cost:** +CPU hour 증가 0% (동일 총 CPU-hour, 벽시간만 단축), queue 경쟁 ↑.
- **리스크:** (a) array task 간 `_spades_run` 같은 공유 scratch가 있으면 충돌 — 확인 완료, per-site 파일명은 모두 `site_id` prefix라 안전. (b) `step_dir` 경로가 같아서 Python 프로세스가 같은 파일에 쓰면 충돌할 수 있음 — 각 array task의 workdir를 `step_dir/_array_<site_id>`로 격리 필요.

### B. GNU parallel (단일 노드, 48-core 활용)
- `parallel -j 16 python s05_single_site.py ::: ${SITES[@]}` 내부에서 실행.
- **장점:** SLURM array보다 스케줄 오버헤드 0.
- **단점:** 단일 노드 48 core 한계, host BAM random I/O 경합 심함. 소이빈 9.7 GB BAM에 16 프로세스 동시 읽기 → IOPS 포화 가능. 실제 wall time 개선 측정 필요.

### C. Nextflow/Snakemake 마이그레이션 (장기, 높은 cost)
- **cost:** 엔진 선택 + 파이프라인 재작성 1-2주, 테스트 1주, docs 업데이트.
- **benefit:**
  - 네이티브 SLURM array 지원, per-task resource 분리 (s04 32G vs s05 8G).
  - resume capability (`--resume` / `-resume`), 실패 사이트만 재실행.
  - profile 기반 파티션 스위칭 (`-profile pronghorn`).
- **추천:** Round 2 이후 장기 로드맵으로 별도 plan. 현재 구조로 array job 병렬화(A)부터 시도.

### D. s04 `samtools depth` 병렬화
- `scripts/s04_host_map.py:149-183` Python stream을 `mosdepth` (멀티스레드, ~10x 빠름)로 교체.
- **cost:** `mosdepth` 의존성 추가 (conda-forge에 있음), 50줄 대체.
- **benefit:** soybean s04가 7h → ~5.5-6h (약 15% 단축).

---

## 4. BUG-5 / 14 / 15 / 16 영구 해결책

### BUG-5: Phase 1.5 positive site 하드 상한 + alert
**현재 상태:** `s04b_construct_assembly.py`가 raw contig를 90%/200bp 필터로 1,345→6로 줄여 1,116→100 수준으로 축소했지만, **여전히 UGT72E3 105, AtYUCCA6 100 = 24h TIMEOUT 경계**.

**제안 fix (`scripts/s05_insert_assembly.py:3644-3656`):**
```python
# After classify_site_tiers
n_pos = len(assembly_sites)
MAX_POS_SITES = int(os.environ.get("REDGENE_MAX_POS_SITES", 60))
if n_pos > MAX_POS_SITES:
    log(f"WARNING: {n_pos} transgene-positive sites exceeds MAX_POS_SITES={MAX_POS_SITES}")
    log(f"  Likely causes: (a) element DB over-reacts with host, "
        f"(b) s04b filter too permissive, (c) actual multi-insertion sample.")
    log(f"  To proceed despite risk: REDGENE_MAX_POS_SITES={n_pos}")
    log(f"  To investigate: inspect site_tier_classification.tsv, s04b/contigs.fasta")
    sys.exit(2)  # distinct from runtime errors
```
- Wall time envelope: 60 sites × 30min = 30h < 48h time budget (여유 있게).
- Override path (`REDGENE_MAX_POS_SITES=200`) 로 수동 승인 케이스 지원.

### BUG-14: worktree `db/` symlink 자동화
**현재 상태:** `db/`는 real directory, main worktree 외에서 일부 reference 파일이 빠짐.

**제안 fix (새 파일 `scripts/setup_worktree_db.sh`):**
```bash
#!/bin/bash
# Usage: bash scripts/setup_worktree_db.sh [<main_repo_path>]
set -euo pipefail
MAIN="${1:-/data/gpfs/assoc/pgl/develop/redgene}"
WT_DB="$(git rev-parse --show-toplevel)/db"

if [ -d "$WT_DB" ] && [ ! -L "$WT_DB" ]; then
    echo "ERROR: $WT_DB is a real directory. Move to $WT_DB.bak and re-run."
    exit 1
fi
ln -sfn "$MAIN/db" "$WT_DB"
echo "Symlinked $WT_DB -> $MAIN/db"

# Validate each host ref used by config.yaml samples
python -c "
import yaml, pathlib
cfg = yaml.safe_load(open('config.yaml'))
for s, sc in cfg.get('samples', {}).items():
    hr = pathlib.Path(sc.get('host_reference', ''))
    if not hr.exists():
        print(f'MISSING host_ref for {s}: {hr}')
"
```
- `.claude/` 또는 README에 **worktree 생성 직후 반드시 실행**이라 명시.
- 더 깔끔한 대안은 `git ls-files db/ | grep -v '^db/.*\.bwt$\|\.fa\.amb$\|...'`로 tracked db/ 파일만 커밋, binary는 모두 gitignored → worktree clone 시 `db/`가 empty로 생기고 symlink 쉬워짐.

### BUG-15: atomic write for `db/transgene_db.fa`
**Root cause:** `scripts/s05_insert_assembly.py:867-887` — `transgene_db.fa`를 빌드하는 코드가 `.exists()`만 체크, size=0 체크 없음. 중단된 run이 zero-byte 남기면 `blastn`이 "Database memory map file error"로 실패.

**제안 fix (`scripts/s05_insert_assembly.py:867` 주변):**
```python
transgene_db = element_db.parent / "transgene_db.fa"

def _valid_file(p: Path, min_size: int = 100) -> bool:
    return p.exists() and p.is_file() and p.stat().st_size >= min_size

if not _valid_file(transgene_db):
    log(f"  Building {transgene_db} (atomic write)...")
    tmp = transgene_db.with_suffix(".fa.tmp")
    with open(tmp, "w") as fout:
        # ... existing concat logic ...
        pass
    os.replace(tmp, transgene_db)  # atomic
    subprocess.run(
        ["makeblastdb", "-in", str(transgene_db), "-dbtype", "nucl"],
        check=True,
    )
```
- 같은 패턴을 `transgene_db_clean.fa` (scripts/s05_insert_assembly.py:781), `element_db` rebuild path (line 2691)에 모두 적용.
- `makeblastdb` 출력 파일 (`.ndb/.nhr/.nin/...`)도 함께 size-check 필요: 이들 중 일부가 zero-byte면 `blastn`이 segfault할 수 있음.

### BUG-16: `sbatch` wrapper로 env override
**Root cause:** `SBATCH_PARTITION=cpu-s1-pgl-0` env var가 `#SBATCH --partition` 지시자보다 우선.

**제안 fix (새 파일 `bin/rg-sbatch`, chmod +x):**
```bash
#!/bin/bash
# Wrapper around sbatch that strips SBATCH_* env vars that would override
# script directives. Use for all redgene SLURM submissions.
set -eu
for var in SBATCH_PARTITION SBATCH_ACCOUNT SBATCH_QOS SBATCH_TIMELIMIT; do
    unset "$var"
done
exec sbatch "$@"
```
- README / CLAUDE.md 업데이트: "모든 SLURM 제출은 `bin/rg-sbatch`로" 표준화.
- `run_rerun_*.sh` 내부의 `#SBATCH` 지시자가 이제 실제로 우선 적용.
- 대안: 각 script 상단에 `unset SBATCH_PARTITION ...` 삽입. 덜 우아하지만 포괄적.

---

## 5. CI smoke test 설계

**목표:** 5분 내 완료, pipeline 전체 step 커버, PR마다 실행.

### Synthetic data 생성 (`tests/smoke/gen_data.py`)
- host: rice Chr1 앞 **100 kbp slice** (NCBI `Osativa_323_v7.0.fa` 50000-150000)
- construct: `element_db/common_payload.fa`의 bar + P-CaMV35S (이미 레포에 있음, ~2 kb 합본)
- inserted reads:
  - host 전체에 **50x coverage** synthetic Illumina reads (wgsim or Mason)
  - position 75000에 construct 5 kb 삽입 → ground truth
  - **총 FASTQ ~20 MB, 100bp PE, 50k pairs**
- Ground truth 사이트: Chr1_slice:75000 → 기대 verdict: CANDIDATE

### pytest integration (`tests/smoke/test_pipeline_smoke.py`)
```python
@pytest.mark.smoke
def test_pipeline_end_to_end(tmp_path):
    synthetic_dir = Path(__file__).parent / "data"
    r1, r2 = synthetic_dir / "smoke_R1.fq.gz", synthetic_dir / "smoke_R2.fq.gz"

    # Step 1-5 in-process
    subprocess.run([
        "python", "run_pipeline.py",
        "--sample", "smoke_test",
        "--steps", "1-5",
        "--threads", "4",
        "--config", str(synthetic_dir / "smoke_config.yaml"),
        "--outdir", str(tmp_path),
        "--no-remote-blast",
    ], check=True, timeout=300)  # 5 min hard cap

    reports = list(tmp_path.glob("smoke_test/s05_insert_assembly/insertion_*_report.txt"))
    assert len(reports) >= 1
    verdicts = [r.read_text().splitlines()[0] for r in reports]
    assert any("CANDIDATE" in v for v in verdicts), f"No CANDIDATE verdict: {verdicts}"
```

### GitHub Actions (`.github/workflows/smoke.yml`)
```yaml
name: Pipeline smoke test
on: [pull_request]
jobs:
  smoke:
    runs-on: ubuntu-latest
    timeout-minutes: 10
    steps:
      - uses: actions/checkout@v4
      - uses: mamba-org/setup-micromamba@v1
        with: { environment-file: environment.yml, environment-name: redgene }
      - run: pytest tests/smoke/ -v --timeout=300
```
- **제약:** GitHub Actions runner는 ubuntu-latest, 4 CPU, 7 GB RAM. Synthetic data를 100kbp로 제한해야 메모리 내 들어감.
- **SPAdes / Pilon 설치 시간**: micromamba 캐시 사용하면 2-3분, 초기 setup은 5-6분. 전체 워크플로 wall clock ~10분 예상.

### Baseline regression suite (`tests/regression/`, weekly cron)
- rice_G281 real data (이미 레포에 있음, ~ 6.3 GB BAM → 아마 공개 불가, internal SLURM으로만 실행)
- Expected: Chr3:16,439,674 → CANDIDATE, FP rate ≤ 95%.
- 매주 SLURM cron으로 수동 실행하고 결과를 `docs/regression/YYYY-WW.tsv`로 저장.

---

## 6. SLURM 제출 best practice 표준화

### Template (`run_rerun_*.sh` 개선안)
```bash
#!/bin/bash
#SBATCH --job-name=rg_<sample>_<step>
#SBATCH --partition=cpu-s2-core-0
#SBATCH --account=cpu-s2-pgl-0
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G            # 4 GB/thread per feedback_memory_threads.md
#SBATCH --time=24:00:00
#SBATCH --output=results/%x_%j.out
#SBATCH --error=results/%x_%j.err
#SBATCH --export=NONE        # clean env, avoid BUG-16

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene
cd /data/gpfs/assoc/pgl/develop/redgene

# Self-check
python -c "import sys; assert sys.version_info >= (3, 11)"
python -c "import pysam, Bio, yaml"
[ -f db/$(yq ".samples.<sample>.host_reference" config.yaml) ] || { echo "host_ref missing"; exit 1; }

python run_pipeline.py --sample <sample> --steps <range> --threads $SLURM_CPUS_PER_TASK --no-remote-blast
```
- `--export=NONE`: cleanest env. `SBATCH_PARTITION` override 자동 회피 (BUG-16).
- `#SBATCH --output=results/%x_%j.out` 패턴: job name이 파일명에 포함되어 grep 편리.
- `--threads $SLURM_CPUS_PER_TASK`: SLURM이 준 thread 수를 Python에 그대로 전달. `feedback_memory_threads.md`와 정합.

### 제출 CLI 통일
- 모든 `sbatch` 호출은 `bin/rg-sbatch` wrapper 통과. CLI flag로 sample 차이 흡수:
  ```bash
  bin/rg-sbatch --mem=128G --time=48:00:00 run_rerun_ugt72e3.sh
  ```
- `bin/rg-sbatch`가 sbatch 호출 전에 `--mem`, `--cpus-per-task`가 `4 GB × cpus` 규칙을 따르는지 검사, 안 맞으면 경고.

### 자원 테이블 (config로 관리)
`config/slurm_profiles.yaml` 새로 추가:
```yaml
rice_G281:        { mem: 32G, time: 6h, cpus: 8 }
tomato_Cas9_*:    { mem: 32G, time: 4h, cpus: 8 }
cucumber_line*:   { mem: 96G, time: 18h, cpus: 16 }
soybean_AtYUCCA6: { mem: 64G, time: 36h, cpus: 16 }
soybean_UGT72E3:  { mem: 128G, time: 48h, cpus: 16, array: per-site }
```
`run_pipeline.py --submit-slurm` 플래그 추가해서 이 config에서 자동으로 sbatch script 생성+제출.

---

## 7. 다른 teammate에게 묻고 싶은 것

### → bio_king (`s05` 3806줄 refactor 담당)
1. **Phase 2+3 로직을 독립 CLI (`s05_single_site.py`)로 분리할 때**, 어떤 helper (`StrandAwareSeedExtender`, `extract_unmapped_paired`, `recruit_by_kmer`)가 site 경계를 넘어 mutable state를 공유하나요? 제가 `extract_unmapped_paired`는 cached 파일 기반이라 안전해 보이는데, `extender` state는 per-site 새로 생성하는지 확인 부탁.
2. Refactor 후 per-site CLI의 input/output 계약(`InsertionSite` serialization 포맷)이 정해지면, 제가 만들 SLURM array wrapper(`run_s05_array.sh`)와 인터페이스 맞춰야 함.

### → gmo_expert (element DB 정합성)
1. **DB를 `gmo_combined_db.fa` (192 seqs)에서 가령 1000 seqs로 확장**하면 Phase 1.5 `classify_site_tiers`의 BLAST wall time이 몇 배 길어지나요? 현재 UGT72E3에서 batch BLAST가 Phase 1.5 전체 ~15분의 주요 비중.
2. 역으로 현재 UGT72E3/AtYUCCA6 100+ positive sites의 절대값을 줄이려면 DB의 어떤 element가 soybean 호스트와 과도 교차하나요? (bar? 35S?). 이들만 masking하는 접근이 더 정확도-친화적일 수 있음.

### → haibao (알고리즘 대안)
1. **BWA → minimap2 전환** 시 wall time 예상? s04 host mapping에서 1.1 Gbp soybean + 2×50 MB reads 기준으로, minimap2 `-ax sr`이 BWA-MEM 대비 1.5-2x 빠르다는 paper 있지만 실제 측정 필요.
2. **k-mer 기반 (mash / dashing) pre-screening**으로 Phase 1.5 positive site를 먼저 shortlist하면 BLAST 부하를 줄일 수 있나요? 정확도 trade-off가 어떻게 되는지?
3. s05 per-site assembly 자체가 필요한가? `minimap2 --splice` 로 insert junction + construct 직접 정렬 + SV caller로 대체하는 안이 있다면 속도-단순성 양쪽 이득.

---

## 8. Round 2 예상 갈등 포인트 (기록용)

- **vs gmo_expert — 정확도 vs 속도:** DB 확장이 recall↑, BLAST 시간↑, Phase 1.5 positive site 증가 → wall time ↑. 제 입장: "48h cap 안에서 돌아가는 정확도가 production-ready". gmo_expert는 "FN 1개가 치명적, 시간은 hardware로 해결". → per-sample profile로 DB 크기 config화 제안.
- **vs haibao — 알고리즘 교체 비용:** minimap2 전환은 wall time 15-30% 개선 예상, 코드 200줄 수정. 큰 refactor는 test coverage 부족 상태에서 리스크. 제 입장: smoke test 먼저 → 알고리즘 교체. haibao는 단번에 교체 선호할 수 있음.

---

## 9. 액션 요약 (Round 1 산출)

1. **즉시 적용 가능 (1 commit, <1일)**:
   - BUG-15 fix: `s05_insert_assembly.py:867-887`의 atomic write + size check
   - BUG-16 fix: `bin/rg-sbatch` wrapper 추가 + README 갱신
   - BUG-5 guard: `MAX_POS_SITES` env var 가드
   - `scripts/setup_worktree_db.sh` 추가 (BUG-14)

2. **Round 2 합의 필요**:
   - s05 per-site array job (bio_king refactor 결과 대기)
   - SLURM submit wrapper 표준화
   - CI smoke test pytest/GHA

3. **장기 로드맵**:
   - mosdepth 전환 (s04 15% 단축)
   - Nextflow/Snakemake 마이그레이션 (2-3주 공사)

---

**다음:** TaskUpdate Task #3 → completed, team-lead에 한 줄 보고.
