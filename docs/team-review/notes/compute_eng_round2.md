# Round 2 — compute_eng 답변: Q1/Q2/Q3

**Agent:** compute_eng · **Date:** 2026-04-16 · **Scope:** HPC 영향 정량화

---

## Q1 — BLAST 호출 감축 (bio_king + haibao 통합)

### makeblastdb per-sample 검증 (s05:867 근처)

**확인 결과 — per-sample 여러 번 호출됨 (비효율):**
| 위치 | 용도 | 빈도 |
|---|---|---|
| `s05:735` | host ref BLAST DB | per-sample 1회 (idempotent, `.ndb` 체크) |
| `s05:802` | `transgene_db_clean.fa` (host-endo 제거판) | per-sample 1회 |
| `s05:879+887` | `transgene_db.fa` (element + UniVec) | per-sample 1-2회 (체크 2곳) |
| `s05:2691` | `_run_local_blast` 내부 element_db temp index | **per-site 1회 × N sites** |

**즉, s05:2691의 makeblastdb가 진짜 반복**. UGT72E3 105 sites × ~1초 makeblastdb = ~2분. 전체 wall time 비중 낮지만, CLI가 stdout DEVNULL이라 log로 안 보임 (검증 불가는 아님, code path 분석으로 확정).

### 예상 단축 (UGT72E3 기준 실측 기반)

**Phase 1.5 batch BLAST:**
- 현재: `classify_site_tiers` 내부에서 `_batch_check_element_hits`가 element DB BLAST 1회, `_filter_host_endogenous`가 host BLAST 1회 = **2회 호출, ~15분 측정치**.
- bio_king 통합 래퍼 + 결과 캐시 공유: 2회 → 1회 호출 병합은 어려움 (DB가 다름). 대신 **DB index rebuild 중복 제거**만 가능 → **약 30초 단축 (무시 가능)**.
- 단축 예상: **Phase 1.5 15분 → 14분 (~7%)**.

**Phase 3/4 per-site BLAST (실측 call-site 카운트 — team-lead 요청):**

`grep -nE "(_run_local_blast|_run_remote_blast|_batch_check|_blast_insert_vs_host|\"blastn\")" scripts/s05_insert_assembly.py`로 확인한 **per-site 루프 내부 호출** (line 3662-3728 context):

| Line | 호출 | 용도 | per-site 빈도 |
|---|---|---|---|
| 2849 | `_run_local_blast(insert, element_db)` via `annotate_insert` | element annotation (primary) | **1회** |
| 2855 | `_run_local_blast(insert, extra_db)` via `annotate_insert` | s04b contig annotation | **1회** (extra_db 있으면) |
| 2865 | `_run_remote_blast(insert)` via `annotate_insert` | NCBI nt (실험에서는 `--no-remote-blast`로 skip) | 0회 (현재) |
| 2898 | `blastn border_fa vs insert_fasta` | border detection | **1회** |
| 2937 | `blastn megablast` (construct-flanking) | Filter B | **1회** |
| 3015 | `blastn megablast` (reuse cached) | Filter C | **1회** (캐시 hit 시 skip) |
| 3075 | `blastn megablast` | Filter D (alt-locus) | **1회** |
| 3139 | `_blast_insert_vs_host` megablast | host-fraction (Filter A) | **1-2회** (line 3315+3375) |

**Per-site BLAST 총합: 6-8회** (remote 제외). UGT72E3 105 sites × 7회 = **~735 BLAST 호출 per sample**.

**단일 BLAST 호출 시간 (실측 추정, insert 수 kb × DB 155 kb × host 1.1 Gbp):**
- Element DB BLAST (작은 DB): per-site ~10-20초 × 3회 = 30-60초
- Host megablast (1.1 Gbp): per-site ~30-60초 × 3-4회 = **2-4분**
- Border/flanking BLAST: per-site ~10-20초 × 2회 = 20-40초
- **per-site BLAST 합계: 3-5분** (28min의 **11-18%**)

**Phase 3 k-mer 확장 + Pilon (나머지 23-25min, 82-89%)가 진짜 병목**. BLAST는 **사실상 부차적**.

**통합 (bio_king 래퍼 + haibao in-memory 공유) 예상 단축:**
- In-memory 공유로 3-4회 중복 제거 → per-site **1-2분 단축**
- UGT72E3: per-site 28min → **26-27min (-5 ~ -7%)**
- 총 wall time: 48h → **45-46h (-4 ~ -6%)**

**통합 결론 (UGT72E3):**
- Phase 1.5: **-7%** (15 → 14min, DB rebuild 중복 제거)
- Phase 3/4 per-site: **-5 ~ -7%** (28min → 26-27min, BLAST 중복 제거)
- **전체 s05: -5 ~ -7%** (48h → 45-46h).

**중요한 재평가:** Round 1 prep에서 저는 per-site BLAST 영향을 -10%로 추정했지만, 실측 call-site 카운트 기반으로는 **-5 ~ -7%가 정확**. Phase 3 k-mer/Pilon이 82-89%를 차지하므로 BLAST refactor의 walltime ROI는 이 정도가 상한.

**결론**: BLAST 래퍼 refactor의 값은 **walltime보다 code quality / maintainability**에 있음. 29곳 중복 제거는 적극 지지하지만 성능 논거는 약함. **진짜 단축은 array job 병렬화(48h → 4h, -92%)**.

---

## Q2 — minimap2 dual-anchor 전환 (haibao)

### Phase 1 site scan 전환 영향

**현재 `find_softclip_junctions` (Phase 1):**
- Input: host BAM의 soft-clip reads iterate (pysam fetch)
- UGT72E3: 28,256 candidate sites 생성. err line 24 → 56546 = **~20분** (log rate 기준)
- BAM I/O 바인드, single-threaded pysam

**minimap2 `--split-prefix` dual-anchor 전환:**
- minimap2가 chimeric alignment를 SA:Z tag + primary 양쪽 모두 출력 → junction 직접 검출, Python 후처리 최소
- Input은 **원본 reads**를 재매핑해야 함 (s04 BAM 재사용 불가) — **비용 증가**
- 대안: s04 BAM의 unmapped + supplementary만 minimap2에 re-feed. 이 경우 ~5-10 GB 데이터 재매핑

**예상 수치 (soybean 1.1 Gbp, 47M reads 기준):**
| 항목 | 현재 (pysam iterate) | minimap2 dual-anchor | 변화 |
|---|---|---|---|
| Phase 1 wall time | ~20min | **~40-60min** (재매핑 cost) | **+2-3x** ❌ |
| Phase 1 MaxRSS | ~2 GB (pysam) | **~5 GB** (minimap2 index) | +2.5x |
| Site 수 | 28,256 | 예상 유사 (±10%) | ~동일 |
| Sensitivity | 현재 baseline | **MAPQ 낮은 microhomology anchor 손실 위험** | ⚠ BUG-7 회귀 가능 |

**통합 UGT72E3 비교:**
- 현재: Phase 1 (20min) + Phase 1.5 (15min) + Phase 2-3 (48h) = ~49h → TIMEOUT
- Array job만 적용 (baseline): 20+15+ (28min × 4 waves of 25 parallel) = **~2.2h**
- Array + minimap2 Phase 1: 50+15+ (28min × 4 waves) = **~2.6h** (+20% 악화)
- Array + minimap2 Phase 1 + haibao BLAST refactor: 50+14+ (25min × 4 waves) = **~2.4h**

**결론:** minimap2 dual-anchor 전환은 **Phase 1에 negative ROI**. 현재 pysam iterate는 20분이고 병목 아님. 48h의 실질 구성요소는 per-site loop(93%)다. minimap2로 교체할 지점은 s04 host mapping (→ -40%, 긍정적)이지, Phase 1 scan 아님. haibao 주장 재검토 요청 — BUG-7 회귀 리스크 + wall time 악화 조합.

---

## Q3 — DB 확장 비용 (gmo_expert P0/P1 15-30 seqs)

### 실측 baseline
- 현재: `element_db/gmo_combined_db.fa` 131 seqs + `common_payload.fa` 9 seqs = **140 seqs, 155,071 bytes**
- Phase 1.5 batch BLAST wall time: UGT72E3 ~15분 (err line 56546 → 58409)
- Phase 3/4 per-site BLAST 합: per-site ~2-4분 × 105 sites = **~5h** (전체 48h의 10%)

### 확장 시나리오 (P0/P1 15-30 seqs 추가, CD-HIT 후)
- **+15 seqs × 평균 1 kb = +15 kb → 155 → 170 KB (+10%)**
- **+30 seqs × 평균 1 kb = +30 kb → 155 → 185 KB (+20%)**
- 보수적: SpCas9 4.1 kb 포함 시 +30 kb, sgRNA/element 혼합 시 평균 700 bp → +21 kb

### 예상 wall time 영향
| 항목 | 현재 (UGT72E3) | +15 seqs | +30 seqs |
|---|---|---|---|
| makeblastdb (per-sample 2회) | ~2초 | ~2.2초 | ~2.4초 |
| Phase 1.5 batch BLAST | 15분 | **16분 (+7%)** | **17분 (+13%)** |
| Phase 3/4 per-site BLAST (5h 합계) | 5h | 5.4h (+8%) | 5.8h (+16%) |
| **총 s05 wall time** | **48h** | **48.3h (+1%)** | **48.6h (+1.3%)** |

**BLAST 시간 자체는 20% 한계 내.** 진짜 리스크는 **positive_count 증가** (Round 1 prep에서 지적):
- 새 element가 host와 교차하면 sites 100 → 120-130 → per-site loop +5-7h (+10-15%)
- SpCas9 (Streptococcus pyogenes origin)은 plant host와 homology 낮음 → 안전
- sgRNA는 synthetic, host와 무관 → 안전
- CRL 82 amplicons는 "amplicons" 자체는 GMO assay용이라 homology 낮을 것으로 예상 (gmo_expert 확인 필요)

**결론 (BLAST 비용만):**
- **+15 seqs: s05 wall time +1% (16분 → 16분)** ✅ 안전
- **+30 seqs: s05 wall time +1.3% (48h → 48.6h)** ✅ 안전 (array 병렬화 적용 후 4h → 4.1h)
- **20% 한계 내 충분히 유지 가능**.

**단서:** 새 entries의 **host cross-reaction 사전 검사** (`scripts/s05_insert_assembly.py:_filter_host_endogenous` 로직을 확장 DB에도 적용) 필수. 이 검사 자체가 BLAST 1회 추가 — per-sample ~1분, 허용 범위.

---

## 종합 권고 (Round 2 compute_eng 입장)

1. **BLAST 래퍼 refactor (bio_king+haibao)** — code quality 지지, 성능 ROI는 **-5 ~ -7%** (48h → 45-46h). 값은 유지보수성.
2. **minimap2 dual-anchor Phase 1 교체 (haibao)** — **반대**. Wall time +2x, MAPQ 손실로 BUG-7 회귀 위험. s04는 minimap2 전환 긍정적, Phase 1은 NO.
3. **DB 확장 (gmo_expert)** — **찬성**. BLAST 비용 +1-1.3%, 20% 한계 내. 선결 조건: host cross-reaction pre-filter를 새 entries에 적용.

**진짜 wall time 해결책 순서 (ROI 기준):**
- (a) **s05 per-site array job**: UGT72E3 48h → 4h (-92%), 코드 변경 ~200줄
- (b) s04 minimap2 전환: 5-7h → 3-4h (-40%), 15줄 변경 (rice_G281 PoC 선행 검증 필요 — 부록 D 참조)
- (c) BLAST 래퍼 refactor: -5 ~ -7%, 29곳 정리 (유지보수 중심)
- (d) DB 확장: +1% (비용, 정확도 개선)

---

## 부록 A — bio_king ↔ compute_eng: s05 refactor 인터페이스 책임 경계

### 제안 CLI 계약 (compute_eng 소비자 관점)

**Phase 분리 CLI 옵션 (bio_king 책임):**
```
python s05_insert_assembly.py \
  --phase 1,1.5 \
  --host-bam ... --host-ref ... --element-db ... \
  --outdir ... --sample-name ... \
  --site-list-out <step_dir>/positive_sites.json
```
- 출력: `positive_sites.json` — `[{"site_id", "host_chr", "pos_5p", "pos_3p", "seed_5p", "seed_3p", "confidence"}, ...]`
- Phase 1.5 batch BLAST 결과는 기존 `site_tier_classification.tsv` 유지.

```
python s05_insert_assembly.py \
  --phase 2,3 \
  --site-id <site_id> \
  --positive-sites <step_dir>/positive_sites.json \
  ...
```
- 출력: `<step_dir>/insertion_<site_id>_insert.fasta` (단일 사이트). Per-site 파일은 `site_id` prefix로 격리.

```
python s05_insert_assembly.py \
  --phase 4 \
  --positive-sites <step_dir>/positive_sites.json \
  --insert-fasta-dir <step_dir>
```

### compute_eng 책임 (SLURM array wrapper `bin/run_s05_array.sh`)
```bash
J1=$(rg-sbatch --parsable --cpus=8 --mem=32G --time=2h \
      --wrap="python ... --phase 1,1.5 ...")
J2=$(rg-sbatch --parsable --array=0-$((N-1))%25 --dependency=afterok:$J1 \
      --cpus=8 --mem=16G --time=1h \
      --wrap='SITE=$(jq -r ".[$SLURM_ARRAY_TASK_ID].site_id" ...); python ... --phase 2,3 --site-id $SITE ...')
J3=$(rg-sbatch --parsable --dependency=afterok:$J2 \
      --cpus=4 --mem=8G --time=30m \
      --wrap="python ... --phase 4 ...")
```
- `%25` = 동시 최대 25 task (GPFS I/O + BLAST DB 공유 경합 방지).
- Per-site 16G: cucumber 실측 87/96G ÷ 30 sites ≈ 3G/site + 여유 → 16G 안전.

### 미해결 질의 (bio_king에게)
1. `extract_unmapped_paired`가 `step_dir/unmapped_R{1,2}.fastq.gz` 캐시를 per-site loop에서 생성. **array task 동시 실행 시 race** → Phase 1 단계로 캐시 생성 끌어올림 필요.
2. `StrandAwareSeedExtender`가 per-site 새 인스턴스인지 확인. 인스턴스화 cost (s03 1.3 MB FASTQ 로드)가 site마다 반복되면 array task 당 ~30s 오버헤드 — 100 sites × 30s = 50min 누적.

---

## 부록 B — gmo_expert ↔ compute_eng: DB 확장 정량 모델 (Q3 확장)

### BLAST 비용 모델

`_batch_check_element_hits`: query = clip seqs (수천~수만), DB = transgene_db.fa (element_db + UniVec + extra). `blastn` wall time:
```
T(blastn) ∝ (query_len) × (db_len) / (threads × cache_efficiency)
```

**DB 크기 변화:**
- 현재: ~140 seqs, 155 KB
- +15 seqs (P0/P1 minimum): +15 kb → 170 KB (+10%)
- +30 seqs (P0/P1 full + Round 1 prep의 95 seqs 중간): +30 kb → 185 KB (+20%)
- Round 1 prep 최대치 (+95 seqs): +100-200 kb → 255-355 KB (+65-130%) — 비합의 초기 추정

Q3 답변의 +15/+30 seqs 시나리오는 **실제 gmo_expert CD-HIT 후 범위**와 정합. Round 1 prep의 "1.7x DB" 수치는 지나친 상한 — Q3 답변을 정설로 채택.

### 미해결 질의 (gmo_expert에게)
1. CD-HIT 후 P0/P1 entries의 **fragment-size distribution**: 평균 길이, min/max?
2. SpCas9 (4.1 kb CDS)가 soybean/rice host 유전자와 **homology 사전 검사 결과**? `_filter_host_endogenous` 확장 적용 전 단독 BLAST로 사전 선별 가능.
3. CRL 82 amplicon set의 구체 서열 (`element_db/crl_amplicons_raw.tsv` 참조 가능한지).

---

## 부록 C — haibao ↔ compute_eng: minimap2 s04 전환 PoC 설계 (Q2 확장)

### rice_G281 PoC 실험 설계

**목적:** s04 BWA → minimap2 전환의 wall time + MaxRSS + downstream verdict 영향 실측.

**프로토콜:**
1. `scripts/s04_host_map.py:86-99`만 교체 (bwa mem → minimap2 -ax sr):
   ```python
   mmi = host_ref.with_suffix(".mmi")
   if not mmi.exists():
       subprocess.run(["minimap2", "-d", str(mmi), str(host_ref), "-k", "21"], check=True)
   cmd = ["minimap2", "-ax", "sr", "-t", str(threads), "-R", rg_tag, str(mmi), str(r1), str(r2)]
   ```
2. rice_G281에 대해 s01 결과 재사용 + s04만 재실행.
3. s05 재실행 (`--steps 5 --no-remote-blast`).
4. 비교 항목:
   - Wall time (s04), MaxRSS (sacct)
   - Chr3:16,439,674 verdict (CANDIDATE 유지?)
   - Phase 1 site 수 비교 (28,256 range 유지?)
   - MAPQ 분포 (samtools view | awk for soft-clip reads)

**승격 기준:** Chr3:16,439,674 CANDIDATE 유지 + Phase 1 site 수 ±20% + wall time -30% 이상 → s04 교체 승격.

**소요:** 총 wall time ~4-6h (s04 ~1h + s05 rerun ~1-2h + 분석). won_yim Round 3 synthesis의 Level 2 rewrite scope 판정 evidence로 활용 가능.

### 미해결 질의 (haibao에게)
1. minimap2 supplementary MAPQ 분포가 BWA와 다르면 Phase 1 `find_softclip_junctions`의 `MAPQ 0-19 microhomology anchor read` 손실 위험. minimap2는 low-MAPQ supplementary 출력 방식이 달라 주의 필요.
2. `minimap2 --splice` + short-read SV caller (Sniffles/SVIM)로 s05 대체 가능? paper 공통 "short read SV caller는 T-DNA 100-500 bp 검출 sensitivity 낮음" → 포기 권장인지 재검토.
3. mash/dashing pre-screening: 28,256 candidate → 105 positive에서 BLAST 대체. 예상 wall time 개선 vs sensitivity loss 수치?
