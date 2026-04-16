# Haibao — Round 1 알고리즘 재평가 (2026-04-16)

**Reviewer 포지션:** JCVI/MCscan/LAST/ALLMAPS 관점의 bioinformatics veteran. assembly-based
transgene insertion detection은 short-read 시대의 레거시 타협이라는 전제에서 시작함.
**원칙:** 비판은 반드시 대안과 함께. short-read lock-in은 깨야 하지만 rewrite는 단계적.

읽은 자료: `CLAUDE.md`, `resume.md`, `bug.md`, `docs/superpowers/plans/2026-04-15-validation-and-cleanup.md`,
`scripts/s05_insert_assembly.py` 전체 흐름 (site discovery → classify → assemble →
annotate → 4-filter → verdict), `scripts/s03b_homology_filter.py`, `scripts/s06_indel.py` 요지.

---

## 1. 현 방식 비판 — 무엇이 근본적으로 취약한가

### 1.1 "assembly-based junction detection"은 split-read detection을 다시 발명한 것이다

`find_softclip_junctions`(s05_insert_assembly.py:373)는 SR(split-read)을 "soft-clip
consensus 클러스터"로 재구현하고 있다. 이 접근은:

- **정확히 2004~2010년대 SV caller들이 한 일**이다. Delly, Manta, GRIDSS, SvABA가
  이미 SR + discordant PE + breakpoint assembly를 production 품질로 풀어놓았다.
  RedGene는 그 중 **SR만** 쓰고 있고, **PE discordancy 신호를 완전히 버리고 있다**.
- `min_clip=20`, `MIN_CLUSTER_DEPTH=3`, `cluster_window=50`은 hard-coded.
  BUG-7 (2026-04-15 리버트)가 정확히 이 문제 — 파라미터 한 벌로 모든 host/construct
  조합을 덮으려고 하다가 microhomology 기반 junction을 통째로 잃었다. 이것은
  **튜닝 문제가 아니라 설계 문제**다. 원본 GRIDSS는 breakpoint-graph를 만들어서
  여러 signal을 동시에 요구하는 방식으로 이 class의 false negative를 체계적으로 회피한다.
- "bidirectional soft-clip at same position + consensus different + not in host"라는
  3-조건 validation은 깔끔하지만, 실제로는 **condition 3 ("not in host")만이 discriminative** 이고
  나머지 둘은 coverage가 약한 샘플(5x, BUG-7 대상 A2_3 같은 case)에서는
  systematic false negative source로 작용한다. Paired-direction 요구는 T-DNA의
  절반이 host repeat에 걸려 있는 경우(rice Chr3 head-to-head 2-copy가 바로 그런 예)
  `junction_3p`가 모호해져서 single-direction fallback으로 떨어지고, 거기서 `≥5 reads`
  hard cut에 걸린다.

### 1.2 "assembly"가 assembly가 아니다 — k-mer + minimap2 + Pilon + SSAKE 콜라주

`assemble_insert`(s05_insert_assembly.py:2363)는 per-site에 **4개의 서로 다른
extension 엔진**(k-mer, minimap2 soft-clip, Pilon, SSAKE)을 라운드마다 순차 호출하고
`max_rounds=8` (기본값), `s05_stats`에 convergence 기록만 남기는 구조다. 문제:

- **per-site 14분 × 105 sites = 24.5h** (resume.md OPEN-1, soybean_UGT72E3 TIMEOUT).
  이것은 "파이프라인을 빠르게 만드는 엔지니어링 문제"가 아니라 **알고리즘 선택이 잘못된
  것**이다. per-site 국소 assembly를 105번 도는 대신 **전체 construct-hitting read
  pool 한 번에 global de novo assembly** (SPAdes/MEGAHIT) 후 host에 map back 하는
  것이 정답이다. 실제로 s04b가 이미 그 구조로 1,345 → 6 contigs 필터링으로 잘
  작동하고 있다. **s04b의 아이디어를 s05로 올려야** 한다.
- 4개 엔진 stack은 서로 다른 데이터 표현에서 같은 읽기를 반복 처리한다. Pilon은
  **polish**용이지 extension 용이 아니다 — `gap_size=1000`의 N-run을 정답으로 fill하는
  것은 Pilon의 공식 usage가 아니다. 이것은 hack이다.
- k-mer extension은 `k=15, min_depth=2, min_ratio=0.7` — **15-mer는 plant genome
  크기 (1Gbp급 soybean)에서 unique하지 않다**. 4^15 = 10억, soybean 1.1Gbp면
  random하게 ≈1번 나오지만 repetitive sequence (plant genome의 50%+)에서 k-mer path가
  cycle로 들어가는 것은 예견된 실패다. "cycle detection"이 코드에 들어가 있다는
  사실 자체가 이 설계의 취약성을 증명한다 (s05_insert_assembly.py:2594-2609).

### 1.3 BLAST 의존과 minimap2 의존의 혼재 — 결정을 안 내린 설계

한 파일 안에서:
- `find_softclip_junctions`: minimap2 (`_batch_check_maps_to_host`)
- `classify_site_tiers`: **blastn-short** (s05:914)
- `annotate_insert`: **blastn** + `_run_remote_blast` (NCBI nt)
- `_check_chimeric_assembly`, `_blast_insert_vs_host`, `_find_construct_flanking_regions`,
  `_check_construct_host_coverage`: **megablast** (4개 다 별도 호출)

**같은 insert sequence에 대해 BLAST를 4~5번 돌린다**. 각각 `-db host_ref`,
`-subject element_db` 등으로 꾸려서. 이것은:

1. I/O + makeblastdb overhead 누적 (soybean 1.1Gbp host DB makeblastdb만 수 분)
2. **결과 간 일관성 보장이 없다** — 같은 insert의 host hit이 Filter C (chimeric)에선
   ≥98% identity만, Filter A (host-fraction)에선 ≥90%만, 다른 곳에선 threshold가
   또 다르다. 같은 raw alignment table에서 filter 별로 다른 threshold를
   적용해야지 **BLAST를 다시 부르면 안 된다**.
3. BUG-2 (blast stem collision), BUG-3 (bitscore tie), BUG-4 (annotate_insert never
   used extra_dbs) 모두 이 **중복 BLAST 호출의 arm-length coordination** 실패의
   증상이다. 근본 수정은 **BLAST 호출을 한 번으로 합치고, 결과를 in-memory로
   여러 filter가 공유**하는 구조.

### 1.4 WT-based filter (s03b)는 옳지만 **WT가 없을 때는 파이프라인이 반쯤 실명**이다

`s03b_homology_filter.py`는 construct-vs-host minimap2 `asm5`(identity 80%, min 50bp)
로 "host에서 construct랑 닮은 영역"을 호출하고, 거기에 떨어지는 read pair를 버린다.
논리는 타당하다. 그러나:

- **WT sample이 없으면 이 step은 건너뛴다**. 현실의 GMO quarantine assay에서 WT는
  cultivar별로 구하기 어려울 수 있다 (Korean Micro-Tom WT는 있지만 corn B73 WT short
  reads는? soybean Wm82 WT short reads는?). WT가 없을 때 MAPQ=60 false positive가
  어떻게 걸러지는지에 대한 대안이 **없다**.
- minimap2 `asm5`는 assembly-to-assembly preset이라 divergence 5% (= 95% identity)를
  가정한다. `min_identity=0.80`으로 내렸지만, PAF에서 반환되는 block 안의 identity
  계산은 `matching / block_len`이다. asm5 mode는 chaining penalty 구조 때문에
  **80% identity 영역을 애초에 chain으로 안 묶는다**. 이건 `-x asm10` (최대 10% divergence)
  또는 `-k 15 -w 10 -A 1 -B 2 -O 6 -E 2 -s 30 -N 50`처럼 LAST-style sensitive 파라미터로
  바꿔야 한다.
- 더 근본적으로, **homology search에 minimap2는 잘못된 도구다**. LAST (Kiełbasa et al.
  2011)는 custom score matrix + frequency-based seeding으로 **divergent homology
  (60-80% identity)**까지 잡아낸다. Plant T-DNA promoter element들 (Ubi1, Act1, Gt1)은
  cultivar 차이로 75-85% identity 구간에 자주 떨어진다. LAST로 옮기면 `s03b`의
  recall이 올라갈 것이다.

### 1.5 Macrosynteny signal을 완전히 무시하고 있다

RedGene는 **host genome이 잘 annotated 된 것을 전제**로 돌아간다 (rice MSU v7, tomato
SLM r2.0, cucumber B10v3, soybean Wm82). 이 genome들은 전부 closely-related species에
대해 **high-quality synteny block**이 이미 공개되어 있다 (Phytozome, Ensembl Plants).
T-DNA insertion이 host genome에 들어가면:

1. **Macrosynteny block이 깨진다** — insertion이 ≥5kb이면 adjacent gene order가
   교란된다. MCscan (JCVI/jcvi.compara.synteny)으로 rice_G281 Chr3:16.44M 주변 200kb를
   sorghum/maize reference에 대해 synteny block 그려보면, **그 지점에서만 synteny
   block이 끊기거나 inversion처럼 보인다**. 이것은 매우 강력한 orthogonal signal이다.
2. **Read depth의 local anomaly** — insertion site에서는 split reads가 많이 생기면서
   bwa가 map하지 않은 base pair가 국소적으로 쌓인다. `samtools depth`로 mean±3σ 벗어난
   구간을 찾으면 candidate 예측이 가능하다.

현 파이프라인은 (1), (2) 둘 다 signal로 사용하지 않고 오직 soft-clip cluster만 본다.
이건 **single-signal detection**이고 short-read SV caller의 교훈 (multi-signal 요구)을
역행한다.

### 1.6 Long-read readiness가 전혀 없다

GMO regulatory 분야는 **ONT MinION field-deployable**로 이동 중이다 (2024-2025 USDA/
APHIS guidance). RedGene의 코드 구조는:

- `pysam.AlignmentFile(host_bam)`에서 `fetch()`로 alignment를 읽는다 — 여기까지는 long-read도 OK.
- 그러나 `soft-clip consensus` 로직이 `min_clip=20, MIN_CLUSTER_DEPTH=3, cluster_window=50`
  에 고정. ONT raw read의 soft-clip은 1-2kb 단위로 지저분하게 clip되고 position drift이
  5-10bp 정도 흔하다. **이 파라미터로는 ONT read를 아예 사용하지 못한다**.
- `_build_consensus`는 majority vote (≥51%) — ONT/HiFi raw error rate (5-10%)에서
  majority vote는 작동하지만 **polishing이 필수**인데 현재 코드는 short-read 가정만
  있다.

full redesign 시 long-read (bp-resolution junction)와 short-read (coverage/depth)를
**complementary signal**로 받도록 인터페이스를 열어야 한다. 현 구조는 short-read lock-in.

### 1.7 Element DB의 "bitscore tie-break" 버그는 **database design 문제의 증상**이다

BUG-3(bitscore tie에서 univec이 element_db를 이김)이 왜 일어났는가? `element_db`와
`univec_vectors.fa`가 **서로 중복 sequence를 담고 있다**는 사실 때문. `_should_replace`
패치(s05_insert_assembly.py:818)는 source label을 우선하지만, 이는 **database를
정규화했어야 할 것을 runtime에서 땜빵**한 것이다.

정답은:

1. `transgene_db.fa`를 빌드할 때 **element_db 항목과 중복되는 univec 항목은 미리
   제거**한다 (`cd-hit-est -c 0.95` 같은 clustering으로).
2. 각 sequence에 **hierarchical source tag** (e.g., `>BAR__PMID12345__element_db` vs
   `>pCambia1300_bar__univec`) 부여.
3. BLAST 결과에서 같은 "semantic identity"의 hit을 dedupe할 때 source priority는
   **DB level에서** (symlink/build script) 이미 해결.

Runtime `_should_replace` patch는 symptom만 잡는다. `element_db/build_common_payload.sh`
(2026-04-10 신설)가 build time dedup을 할 수 있는 자리다.

---

## 2. 5개 재평가 축 — 현재 방식 vs 대안

### Axis A: Site discovery (soft-clip vs split-read vs discordant PE vs dual-anchor)

| 방식 | 현 RedGene | GRIDSS/Manta | dual-anchor (mm2 --split-prefix) | depth-gap |
|------|-----------|--------------|----------------------------------|-----------|
| Sensitivity | SR only, single-side 가능 | SR+PE+local assembly | SR + split alignment | depth anomaly |
| FP rate | host repeat에 약함 (BUG-7) | multi-signal로 낮음 | split alignment의 MAPQ로 제어 | coverage 외에 약함 |
| 단일 파라미터 민감도 | **매우 높음** (cluster_window, min_clip, MIN_CLUSTER_DEPTH) | 각 signal별 계층적 threshold | minimap2 기본 세팅으로 robust | 통계 프레임워크 | 
| Microhomology 처리 | **실패** (BUG-7) | breakpoint graph가 처리 | split alignment가 자연스럽게 처리 | 못 봄 |
| Short + long read | short만 | Manta는 short only, GRIDSS는 hybrid | minimap2는 both | both |
| Implementation cost | existing (복잡) | 외부 tool import | **가벼움** | 쉬움 |

**추천:** **dual-anchor (minimap2 `--split-prefix` + `--sr` 또는 `-ax sr`)**로 site
discovery를 교체. rationale:
- minimap2 `--split-prefix`는 construct reference를 host에 **chimeric read as one
  alignment**로 그대로 보고한다. T-DNA가 host에 붙어 있으면 read가 construct 20bp +
  host 130bp로 split되고, 두 쪽의 좌표를 동시에 반환한다.
- 이것은 soft-clip cluster 로직의 super-set이면서 **cluster_window tuning이 필요 없다**
  (read length 자체가 cluster 경계).
- GRIDSS급 full breakpoint graph는 overkill이고, minimap2 dual-anchor는 같은 tool
  체인으로 해결된다. PE discordancy는 phase 2에서 추가 가능 (bwa-mem raw BAM의 `TLEN`
  outlier).

### Axis B: Assembly (per-site targeted vs global de novo vs assembly-free breakpoint graph)

| 방식 | 현 RedGene (k-mer+mm2+Pilon+SSAKE) | SPAdes/MEGAHIT global | GRIDSS breakpoint assembly | Flye (long-read) |
|------|-----------------------------------|------------------------|----------------------------|-------------------|
| Wall time | 105 sites × 14 min (OPEN-1) | single SPAdes run | breakpoint-localized, 빠름 | long-read 필요 |
| Coverage 요구 | ≥10x | ≥10x | junction 근처 depth | ≥10x long-read |
| Ploidy/multi-copy 처리 | **palindromic insert에 약함** (head-to-head T-DNA) | meta-mode가 잘 처리 | Δ-path로 해결 | native |
| 외부 tool | Pilon (polishing 오남용) | SPAdes | gridss.jar | Flye |
| Fast iteration | 각 site가 독립이라 병렬 | single global | per-region | per-region |
| Quality | per-site 4-engine stack의 예측불가 품질 | 일관된 quality | breakpoint-focused quality | highest |

**추천:**

1. **단기 (patch)**: 현 `assemble_insert`를 유지하되 per-site round를 **max 3으로 축소**
   (현재 8). 실험적으로 resume.md "most converge by 6"이라 했지만, Pilon iteration은
   **diminishing return**가 급격하다. round 1-3에서 잡히는 junction이 이후 라운드에도
   남는다. 8→3으로 바꾸면 per-site 14min → ~5min, 105 sites가 24h 안에 들어온다.
2. **중기 (rewrite-module)**: s04b의 SPAdes global assembly를 주 엔진으로 격상.
   s05는 global contig를 host에 map back 해서 break/chimera junction을 찾는 쪽으로 바꾼다.
   site → assemble → annotate의 loop를 뒤집어서 **assemble → contigs → map → sites**
   로 바꾼다. s04b contigs.fasta (6개)가 이미 그 역할을 할 수 있다.
3. **장기 (full redesign)**: GRIDSS-style breakpoint graph. breakpoint edge는 SR +
   discordant PE + local assembly 각각의 evidence를 scoring으로 합친다. 이것이 production
   SV caller가 converge한 답이다.

### Axis C: Alignment search (BLAST vs minimap2 vs LAST)

| 용도 | 현재 도구 | 대안 | 추천 |
|-----|---------|------|------|
| clip → host (short, highly similar) | blastn-short | minimap2 `-ax sr` | minimap2 (속도 10x+, plant host 1Gbp DB makeblastdb 불필요) |
| clip → element_db (short, divergent OK) | blastn-short | **LAST** (custom matrix) | LAST — construct element cultivar drift 75-85% identity 대응 |
| insert → host (long-ish, repeat-heavy) | megablast | minimap2 `asm10` | minimap2 (동일 이유) |
| insert → element_db (full annotation) | blastn + NCBI remote | LAST + local curated DB + remote | LAST local + NCBI fallback |
| construct-vs-host flanking | megablast | minimap2 `asm5` | minimap2 (이미 s03b에 있음, 재사용) |

**핵심 메시지:** RedGene는 **BLAST와 minimap2를 정책 없이 섞어 쓴다**. 원칙:
- **short query (≤100bp)**: minimap2 `-ax sr` 또는 blastn-short (둘 다 OK, speed vs
  sensitivity trade-off 명확)
- **medium query (100bp-10kb), divergent OK**: LAST
- **long query (≥10kb), assembly-level**: minimap2 `asm5/asm10/asm20` tier
- **production annotation (element DB)**: LAST + curated DB

BLAST는 `-remote` 할 때 외에는 권장하지 않는다. LAST는:
- JCVI 레퍼런스 구축에서 이미 입증된 tool (MCscanX, LAST pipeline).
- Custom scoring matrix로 GC-rich (maize)와 AT-rich (rice) host 모두 같은 로직으로 다룸.
- DB build 1회, query는 multi-thread, blastn 대비 2-5x 빠름.

### Axis D: Structural signal boost — macrosynteny

**제안:** `scripts/s05_synteny_boost.py` 신설, `generate_report` 직전에 호출.

입력:
1. `site.host_chr`, `site.pos_5p` (insertion 후보 좌표)
2. Pre-computed MCscan synteny block (rice vs sorghum, tomato vs potato, cucumber vs
   melon, soybean vs common bean 등 closely-related reference 쌍으로 1회 구축)

출력:
- `synteny_score`: insertion 후보 좌표 ±50kb 구간 내 synteny block의 **break evidence**
  (inversion, gap, out-of-order gene).
- `disrupted_gene`: insertion이 gene body 안에 떨어지면 gene ID (예: cucumber_line224의
  G6838는 이 로직으로 자동 annotation).

**Scoring 방식 (Haibao 제안):**

```
synteny_score = w1 * (1 if breaks_adjacent_synteny_block else 0)
              + w2 * (1 if inside_gene_body else 0)
              + w3 * (1 if near_TE_boundary else 0)   # TE junction은 natural SV
              - w4 * (1 if inside_known_SV_hotspot else 0)
```

이 score를 verdict에 어떻게 녹일지:
- `CANDIDATE_HIGH_CONF`: 기존 조건 + synteny_score ≥ 2
- `CANDIDATE`: 기존 조건만
- `UNKNOWN → CANDIDATE_LOW_CONF` 승격: element hit이 부족해도 synteny_score ≥ 3이면
  fallback 승격 가능 (BUG-6에서 reverted된 heuristic의 올바른 교체재).

**구현 비용:** MCscan은 `pip install jcvi`로 설치. synteny block 구축은 host 당 1회
(~10분). query는 bedtools intersect급 속도 (O(log N)).

### Axis E: Long-read readiness

현재 구조가 ONT/HiFi를 **흡수 가능한 지점**:

| Step | 현재 | Long-read 대응 |
|------|-----|---------------|
| s04 host mapping | bwa mem | minimap2 `-ax map-ont`/`map-hifi`로 교체 (이미 dependency) |
| s05 find_softclip_junctions | MIN_CLUSTER_DEPTH=3, cluster_window=50 | 같은 로직 재사용 가능하나 window=500, min_clip=100으로 preset 분기 필요 |
| s05 assemble_insert | 4-engine stack | **long-read가 있으면 이 단계가 필요 없음** — read 자체가 full insert를 span |
| s06 indel | `mpileup -Q 0` | ONT는 noisy이므로 pileup 대신 `medaka` 또는 `clair3`로 교체 |

**추천:** config.yaml에 `data_type: short|long|hybrid` 필드 추가. `run_pipeline.py`에서
data_type에 따라 preset 분기. s05는 long-read preset이 주어지면 per-site assembly를
건너뛰고 read 자체를 annotate로 보냄. 이는 **`long-read preset = production-grade`**
로 가는 자연스러운 경로.

---

## 3. 대안 아키텍처 제안 — 3 레벨

### Level 1: Patch (2-5 days, minimal risk)

**목표:** 현 구조 유지, 가장 심각한 FP/FN/TIMEOUT 포인트만 타겟.

1. **cluster_window adaptive** — `find_softclip_junctions`가 `cluster_window`를
   per-chromosome read depth 기반으로 자동 조정. Low coverage (≤5x): 100. Normal
   (≥10x): 50. BUG-7의 리버트를 이 방식으로 대체.
2. **BLAST → LAST 스위치 옵션** — `--aligner last|blast` 플래그 추가 (기본은 blast로
   호환성 유지). classify_site_tiers와 annotate_insert에만 적용.
3. **Assembly round 8 → 3** — `assemble_insert(max_rounds=3)`로 축소. UGT72E3/AtYUCCA6
   의 24h wall을 ~6h로 단축. resume.md "most converge by 6"은 경험칙인데, per-site
   evidence file (s05_stats.txt)을 재분석해서 실제 round-3 convergence rate를
   측정해야 한다 (Haibao 검증 요청).
4. **중복 BLAST 통합** — `_blast_insert_vs_host`, `_check_chimeric_assembly`,
   `_check_construct_host_coverage`가 megablast를 3번 호출한다. **한 번 호출 후 pandas
   DataFrame으로 in-memory 3-filter**. 수 분 절약 × 105 sites.
5. **transgene_db build-time dedup** — `element_db/` Makefile에 `cd-hit-est -c 0.95`
   스텝 추가. BUG-3 근본 해결.

### Level 2: Rewrite-module (2-3 weeks, medium risk)

**대상 모듈:** `find_softclip_junctions` + `classify_site_tiers` (s05의 앞 1200줄).

**설계:**

```
# 신설: scripts/s05a_site_discovery.py
# input: host BAM (from s04), construct ref
# output: sites.tsv (chr, pos, strand, evidence_type, evidence_score, clip_5p, clip_3p)

signals = [
    extract_split_reads(host_bam),            # SR (현 soft-clip cluster)
    extract_discordant_PE(host_bam),          # 새로운 signal
    dual_anchor_mm2(construct_ref, host_bam), # minimap2 --split-prefix
]
sites = merge_signals_with_breakpoint_graph(signals, window=150)
# 각 site는 multi-evidence score
```

다중 signal을 요구하므로 **BUG-7 class의 regression이 알고리즘 레벨에서 차단됨**.
cluster_window 같은 single hard-coded parameter가 전체 sensitivity를 좌우하지 않게 된다.

`classify_site_tiers`는 multi-evidence site를 받으므로 "is_positive = has_element_hit"의
boolean이 아니라 **Bayesian-style scoring** (element_hit weight + PE support weight +
dual_anchor weight + host_exclusion weight)로 바뀐다. UNKNOWN이 줄어든다 (resume.md에서
AtYUCCA6 30 UNKNOWN / 15 FP / 0 CANDIDATE는 bool 분류의 실패; score 기반이면 "low-conf
candidate" tier가 생긴다).

### Level 3: Full redesign (2-3 months, high risk, high reward)

**목표:** assembly-free T-DNA junction graph + long-read ready.

**핵심 idea (GRIDSS-style):**

```
# redgene v2 architecture
reads → minimap2(-ax sr|-ax map-ont) → host BAM
       → minimap2(construct as host) → construct BAM
       → BreakpointGraph(host BAM, construct BAM)
            nodes: genome positions (host or construct)
            edges: read-pair support (SR, PE, split)
BreakpointGraph
   → find_edges(score >= threshold)
   → for each edge (host_pos ↔ construct_pos):
        local_assemble_bubble(edge) if assembly needed
        annotate_with_LAST(element_DB + common_payload + per_sample_s04b)
        synteny_boost(host_pos, macrosynteny_DB)
        verdict = compute_verdict(evidence_vector)
```

장점:
- **short + long read를 같은 pipeline으로 처리** (read-type 분기가 minimap2 preset에
  국한).
- **per-site assembly 24h 문제가 사라짐** — assembly는 breakpoint 주변만, 전체 pool은
  global SPAdes 1회.
- **macrosynteny가 설계 레벨에 통합** — verdict가 local evidence + structural context.
- **production SV caller와 동일한 pattern** — 10년치 SV detection 교훈을 상속.

비용:
- 기존 95% 재작성. test data 검증 (rice G281, tomato A2_*, cucumber line*, soybean)
  전면 재실시.
- GRIDSS처럼 Java 의존을 피하고 싶다면 SvABA (C++) 또는 **structural-variant-caller in
  Rust** (최근 academic prototype) 참고 가능.

**Haibao의 입장:** Level 1은 **다음 주 안에 해야 하는 것** (OPEN-1, OPEN-5 해결에 필수).
Level 2는 **논문 투고 전에 해야 하는 것** (현 BUG-7 class의 regression이 reviewer에게
잡힐 것). Level 3는 **v2로 분리해서 grant proposal**에 담을 것. production GMO detection은
이대로는 regulator가 안 받는다 — FP rate가 사람이 매번 manual review해야 하는 수준.

---

## 4. Macrosynteny boost 설계 — 구체안

### 4.1 어느 케이스에 signal이 강한가

resume.md의 3개 ground truth case를 synteny 관점에서 재검토:

**rice G281 Chr3:16,439,674 (2-copy head-to-head T-DNA, 36bp del, lactoferrin RNAi)**
- Rice MSU v7 Chr3:16.4M 주변 100kb는 **sorghum Sb03 ~15-16Mb**와 high synteny block.
- Head-to-head 2-copy T-DNA (~15kb total + 36bp del)는 **synteny gap**을 만든다.
- MCscan output에서 anchors (ortholog gene pairs) 사이의 "unaligned block" 길이가
  local 평균의 ≥3σ인 지점으로 자동 검출 가능.
- 예측: synteny_score ≥ 3 (block break + inside/near gene + ≥5kb anomaly).

**tomato A2_3 SLM_r2.0ch01:91,002,744 (Cas9 vector, heterozygous)**
- Tomato SLM r2.0 ch01:91M은 potato (Solanum tuberosum) v6.1 chr01 및 S. pennellii와
  synteny block 형성.
- T-DNA size ~15kb가 ortholog gene interval 안에 떨어지면 synteny_score 기여.
- Heterozygous이면 read depth는 보통 수준 유지 — synteny signal이 **depth로 잡히지
  않는 heterozygous 삽입에서 특히 유용**.

**cucumber line224 LKUO03001512.1:581,328 (G6838 disruption, single T-DNA)**
- Cucumber B10v3은 contig-level 조립 (8,035 contigs)이라 chromosomal synteny block이
  melon (Cucumis melo) 대비 제한적.
- 그러나 **G6838 gene disruption**은 synteny의 gene-level event — G6838이 melon ortholog에
  대한 1:1 anchor인지 여부가 signal.
- 예측: 8,035 contig 환경에서 synteny는 per-chromosome 블록 대신 per-contig gene
  collinearity로 재정의해야 함. 이 샘플은 **synteny가 약한 case**이고 element DB signal이
  여전히 주도권.

### 4.2 구현 스케치

`scripts/s05_synteny_boost.py`:

```python
# pre-computed inputs (build 1회):
#   db/synteny/rice_vs_sorghum.anchors.simple
#   db/synteny/tomato_vs_potato.anchors.simple
#   db/synteny/cucumber_vs_melon.anchors.simple
#   db/synteny/soybean_vs_phaseolus.anchors.simple
# 형식: JCVI MCscanX anchors.simple (블록 시작/끝 + 평균 Ks)

def score_insertion_site(
    host_chr: str, site_pos: int,
    synteny_anchors: Path,
    host_gff: Path,            # db/Osativa_323_v7.0.gene_exons.gff3
    window: int = 50000,
) -> dict:
    # 1. 해당 좌표가 synteny block 내부인지, 블록 경계 ±window 내인지 판단
    # 2. 블록 내부면 "내부 + 블록 평균 gene interval 대비 local anomaly"를 측정
    # 3. GFF로 nearest gene / inside gene body 판정
    # 4. synteny_score 산출
    return {
        "inside_synteny_block": bool,
        "breaks_block": bool,
        "local_gene_interval_anomaly": float,  # log-ratio
        "inside_gene": str | None,
        "disrupted_gene_ortholog": str | None,  # sorghum/melon ortholog
        "synteny_score": int,  # 0-5
    }
```

`generate_report`에서 호출 → verdict compute 전에 `synteny_score`를 evidence에 추가.

### 4.3 기여도 평가 방법

이 signal이 실제로 얼마나 discriminative한지는 **test data를 ground truth vs FP로
나눈 뒤 ROC**로 본다:

| 샘플 | ground truth 좌표 | FP 좌표 (현 파이프라인 FP) | synteny_score 차이 |
|------|-----------------|-------------------------|-------------------|
| rice G281 | Chr3:16,439,674 | Chr3:29M, Chr11:2.87M (new CAND) | 측정 필요 |
| A2_3 | ch01:91,002,744 | ch06:41,781,251 (new CAND) | 측정 필요 |
| cucumber line225 | LKUO03001451.1:6,501 | 7개 new CAND | 측정 필요 |

Round 2에서 bio_king이 "element DB가 discriminative"라고 주장하면, 이 ROC 비교로
**synteny가 orthogonal signal임**을 입증해야 한다.

---

## 5. 팀 cross-talk에 던질 도전 질문

### To `bio_king` (s05 3806-line refactor target)

1. "**per-site 4-engine stack (k-mer+mm2+Pilon+SSAKE)** 는 실제로 round별 marginal
   contribution이 있나? s05_stats.txt에서 각 engine별 growth bp를 누적 집계해라.
   Pilon이 실제로 +bp 기여하는 라운드 비율이 어떻게 되나? SSAKE는?"
   — 만약 SSAKE contribution이 10% 이하라면 **SSAKE 제거**가 refactor의 첫 수순.
2. "`max_rounds=8`은 아무 근거가 없다. resume.md에 'most converge by 6'이라 했지만,
   이건 outcome 관찰이지 optimal 선택이 아니다. **max_rounds=3으로 돌렸을 때 CANDIDATE
   복원율이 얼마나 떨어지나?** 이걸 측정하지 않고는 refactor가 무의미하다."
3. "`find_softclip_junctions`는 SR만 쓴다. **PE discordancy signal을 추가하면 BUG-7
   클래스의 regression이 제거되는가?** 내 주장은 yes. 당신의 intuition은?"

### To `gmo_expert` (element DB + FP 케이스)

1. "**element_db (131 entries)와 common_payload (9 entries) 사이 중복**이 얼마나 있나?
   BUG-3이 중복 때문이었다면, cd-hit-est -c 0.95로 dedup 했을 때 N entries가 줄어드나?
   이게 실제 문제의 크기다."
2. "**FP 케이스들 (A2_3 ch06:41.78M, rice Chr11:2.87M, cucumber line225 +7)**에 대해
   synteny block context를 봤나? 이들이 synteny break 시그널이 없다면 FP 확률이 높다는
   내 가설을 검증해 달라."
3. "**host-endogenous exclusion Tier 1/2 threshold (90%/50% vs 75%/30%)**는 어디서
   왔나? cultivar drift 실험적 데이터 기반인가, 임의 선택인가? **rice Nipponbare vs
   Xiushui identity 분포**를 BLAST로 재측정하고 percentile 기반으로 재설정해야 한다."

### To `compute_eng` (HPC/SLURM)

1. "**per-site 14분 × 105 sites = 24.5h TIMEOUT은 알고리즘 문제, 리소스 문제가 아니다**.
   시간을 48h로 늘리는 것은 band-aid. 당신이 리소스 최적화를 얘기할 때, 이 point를
   빠뜨리면 엔지니어링이 엉뚱한 곳을 파게 된다."
2. "**s05 안에서 BLAST가 per-site 3~5회 호출된다** (host-frac, chimeric, construct-host,
   flanking, annotate). 이들은 host DB를 같은 insert에 대해 반복 query하는 중복
   I/O. single-pass batch로 뭉치면 soybean 1.1Gbp host에서 **wall time이 얼마나
   단축되나?** 측정 가능한가?"
3. "`run_pipeline.py`의 step orchestration에서 **s05 per-site loop를 SLURM array로
   병렬화**할 수 있나? s04b contigs.fasta가 사전 가능하면 per-site 독립이라 array
   병렬화는 trivially 가능해야 한다. 왜 안 되어 있나?"

### To `won_yim` (PI, production-ready 기준)

1. "**현재 파이프라인은 WT sample 없이 작동하지 않는다**. Korean quarantine assay에서
   WT 확보가 항상 가능한가? 없으면 s03b를 bypass하는 fallback은 **MAPQ=60 FP에
   대해 완전히 무방비**다. 이것을 production criterion에 포함시켜야 한다."
2. "**FP rate 목표를 숫자로 정의**하자. resume.md에서 cucumber line225가 1 CAND ground
   truth에 8 CAND 보고를 냈다. FP rate를 '사람이 매번 manual review해야 하는 수준
   (>50% FP)'로 남기면 production이 아니다. 목표: per-sample FP ≤1. 이게 spec에
   들어가야 validation이 의미가 있다."
3. "**regulatory 관점에서 bp-resolution junction이 정말 필요한가?** ONT이면 trivial,
   short-read이면 assembly로 쥐어짜고 있다. Korean regulation (lab_QIA)이 요구하는
   resolution이 실제로 몇 bp인가? ±10bp로 충분하면 assembly 부하가 확 줄어든다."

---

## 6. 참고 논문 (Round 2 준비)

> Round 1 단계에서는 아직 외부 웹 검증을 돌리지 않았음. 아래는 기억 기반 citation, Round
> 2에서 WebFetch로 확인 예정.

1. **GRIDSS:** Cameron DL, et al. *GRIDSS: sensitive and specific genomic rearrangement
   detection using positional de Bruijn graph assembly.* Genome Research 2017.
   → 현재 RedGene 대비 "multi-signal breakpoint graph"의 production 구현.

2. **minimap2 split-read dual-anchor:** Li H. *Minimap2: pairwise alignment for nucleotide
   sequences.* Bioinformatics 2018. `--split-prefix`, `-ax sr` preset 문서.
   → Axis A / Axis C의 근거.

3. **MCscan / macrosynteny signal:** Tang H, et al. *Synteny and collinearity in
   plant genomes.* Science 2008; jcvi.compara.synteny module (JCVI github).
   → Haibao 본인 authorship (과장 ok). Axis D의 backing.

4. **LAST for divergent homology:** Kiełbasa SM, et al. *Adaptive seeds tame genomic
   sequence comparison.* Genome Research 2011.
   → Axis C element DB 75-85% identity 구간 처리.

5. **Long-read T-DNA characterization:** Jupe F, et al. *The complex architecture and
   epigenomic impact of plant T-DNA insertions.* PLOS Genetics 2019 (ONT data).
   → Axis E long-read readiness 근거.

Round 2에서 각 논문의 **실제 결과 수치**를 확인 후 반박 무기로 사용.

---

## 7. Round 1 정리 — Haibao의 결론

1. **현 파이프라인은 short-read SR-only SV caller의 2010년대 방식을 2026년에 재발명**
   하고 있다. multi-signal (SR+PE+dual-anchor) + breakpoint graph로 옮겨야 한다.
2. **BUG-3, BUG-4, BUG-7은 전부 설계 취약성의 증상**. patch로는 땜빵이고, Level 2
   rewrite가 실질적 해결.
3. **macrosynteny signal은 plant host의 high-quality annotation을 활용하는 orthogonal
   infrastructure**. RedGene가 plant-specific GMO tool인 이상, 안 쓰는 것은 낭비.
4. **per-site 14분 × 105 sites TIMEOUT은 리소스 문제가 아니라 per-site assembly
   설계 문제**. global s04b SPAdes + breakpoint graph가 정답.
5. **long-read readiness 없이는 production-grade이라 부를 수 없다**. regulatory 분야
   트렌드를 역행.

**추천 우선순위:**
1. (이번 주) Level 1 patch 5개 — BUG regression 방어.
2. (1개월) Level 2 module rewrite: s05a (site discovery multi-signal) + MCscan boost.
3. (v2) Level 3 full redesign + long-read preset.

Level 1만 해도 OPEN-1, OPEN-5 (soybean TIMEOUT, 0 CANDIDATE) 해결 가능성 높다. Round 2
에서 bio_king이 "증분 개선"을 주장하면 나는 Level 2를 밀고, `won_yim` PI 결정에
따름.

— Haibao, 2026-04-16
