# RedGene v2.0 — Grant Scope Document

**작성자:** haibao (team: redgene-review-2026-04-16)
**Baseline commit:** `729835d` on `main`
**작성일:** 2026-04-16
**목적:** v1.0 (Korean quarantine MVP) 배포 이후의 technical roadmap을 grant proposal
수준으로 scope out. won_yim(PI)의 "MVP 전 아키텍처 교체 금지" 원칙에 따라 Phase 2 / v2.0로
이관된 알고리즘 개선안을 motivation-spec-ROI 형태로 정리.

---

## 1. Motivation — 현 RedGene v1.0의 Fundamental Limitations

v1.0은 6 host × 8 GT sample에서 AC-1 (≥80% recall) 달성이 가능하고 AC-4 (≤48h)도
patch 적용 후 통과 예상이다. 그러나 다음 **근본적 설계 제약**은 patch 수준에서 해결되지
않으며, production-grade quarantine tool (AC-2 ≤2, AC-3 ≥90%, AC-7 coverage robust)의
요건을 체계적으로 달성하려면 v2.0 redesign이 필요하다.

### 1.1 SR-only legacy approach (Round 1 §1.1)

`find_softclip_junctions` (s05_insert_assembly.py:373)은 split-read evidence만 사용.
Manta/GRIDSS/SvABA가 2015년 이후 확립한 **multi-signal SV calling 원칙** (SR + PE
discordancy + local assembly breakpoint graph)을 한 축만 쓰고 있다. BUG-7 (cluster_window
50→300 + min_mapq 20 필터에서 A2_3 T-DNA junction 전체 소실)은 **단일 파라미터가 전체
sensitivity를 좌우**한다는 사실의 직접적 증거다. 이는 튜닝 문제가 아니라 single-signal
detection의 구조적 취약성이다.

### 1.2 k-mer k=15는 plant genome unique-mer threshold 아래

`StrandAwareSeedExtender.k=15`는 soybean 1.1 Gbp / corn 2.2 Gbp급 host에서 기대
coincidence ~1/k-mer. Jellyfish/KMC 벤치마크에서 식물 유전체는 k=15가 ~40-50% unique,
k=21에서 ~75-85%, k=25-31에서 ~90%+. 현 코드가 "cycle detection" 분기 (line 2594-2609)를
달고 있는 것 자체가 이 k가 repetitive region에서 path loop에 빠진다는 증상이다.

### 1.3 중복 BLAST 호출과 결과 inconsistency

s05 안에서 같은 insert/host 쌍에 megablast가 4회 호출된다 (`_find_construct_flanking_regions`,
`_check_chimeric_assembly`, `_check_construct_host_coverage`, `_blast_insert_vs_host`).
각 filter가 다른 identity threshold (90% / 98% / 80% / 90%)를 적용. 동일 alignment table로
multi-threshold filtering을 하지 않고 BLAST를 다시 부름 → I/O overhead + filter 간
결과 drift. BUG-2/-3/-4는 전부 이 arm-length coordination 실패의 증상.

### 1.4 WT 없으면 MAPQ=60 FP에 무방비 (CLAUDE.md "Known Pitfalls"에서 자가 인정)

`s03b_homology_filter.py`의 WT-based read filter는 Korean quarantine samples에서
`wt_control`이 config에 없는 경우 bypass된다. 이때 plant-derived promoter (Ubi1, Act1,
TA29)가 MAPQ=60으로 host에 unique-map하면서 false positive 생성. gmo_expert Round 1
§3의 교차반응 매트릭스가 이 class 전수 (corn host + P-Ubi1-maize 100%, rice host +
P-Act1-rice 100%). v1.0은 "mandatory WT sample"을 operational requirement로 두지만,
실제 격리 검사에서 matching WT 확보는 불확실.

### 1.5 Macrosynteny signal 무시

Host genome (rice MSU v7 / tomato SLM r2.0 / cucumber B10v3 / soybean Gmax_v4.0)은 모두
closely-related reference (sorghum / potato / melon / common bean)와 high-quality
MCscan synteny block이 publicly available. Insertion은 정의상 lineage-specific이고,
macrosynteny break는 insertion signal의 **orthogonal evidence**. 현 파이프라인은 이
signal을 사용하지 않고 soft-clip cluster만 본다.

### 1.6 Long-read 미대응

USDA APHIS, EFSA의 2024~2025 GMO screening guidance는 ONT MinION field deployment를
reference workflow로 언급하기 시작. HiFi/ONT는 T-DNA junction을 bp-resolution으로
single-read에서 제공 — per-site local assembly가 불필요. v1.0은 `bwa mem` + pileup
`-Q 0` 등 short-read 가정이 코드 전반에 hard-coded. Parameter (`min_clip=20`,
`MIN_CLUSTER_DEPTH=3`, `cluster_window=50`)도 short-read용.

### 1.7 Element DB build-time governance 부재

EUginius 131 + common_payload 9 + univec이 **runtime merge**된 transgene_db.fa를 사용.
동일 서열(`V00087.1` P-nos/T-nos, BUG-9)이나 cross-database 중복이 runtime에서
`_should_replace` priority로 땜빵된다. gmo_expert §2는 **SpCas9 / sgRNA scaffold /
vip3Aa / cry2Ab2 / cry34-35Ab1 / PMI / ipt / T-g7 / T-pinII** 등 15+ element가 DB에
누락된 것을 지적. 이것은 data-curation governance의 문제이지 코드만으로 해결되지 않는다.

---

## 2. Level 2 Spec — v1.1 Candidate (dual-anchor + PE discordancy)

### 2.1 설계 개요

`s05a_site_discovery.py` 신설. 3 signal을 독립 계산 후 breakpoint graph로 merge:

```
signal_SR      = extract_split_reads(host_bam)          # 현 soft-clip cluster 재사용
signal_PE_disc = extract_discordant_PE(host_bam)        # 신규 — TLEN outlier, mate unmapped
signal_dual    = minimap2 --split-prefix(construct → host) # chimeric read as one alignment

sites = merge_signals(
    signals=[SR, PE_disc, dual],
    window=150,                # read length 기반, cluster_window hard-code 제거
    min_signals=2,             # 최소 2 signal 일치 요구
)
```

### 2.2 기대 AC 개선

| AC | v1.0 baseline | v1.1 예상 | 근거 |
|----|-------------|----------|------|
| AC-1 recall | 5/6 (83%) | 6/6 (100%) + tomato_A2_2 | PE discordancy가 low-MAPQ microhomology junction도 잡음 (BUG-7 class 제거) |
| AC-2 FP | cucumber_225: 8 CAND (GT=2) | 3-4 CAND | 2-of-3 signal 요구로 element-cross-reaction FP 자연 감소 |
| AC-3 annotation | 데이터 없음 | 동일 | Level 2는 site discovery만 건드림, annotate_insert 불변 |
| AC-7 coverage | ≥10x 검증, 5x 미확인 | 5x 일부 회수 | dual-anchor는 read 하나로 signal 제공, coverage 의존도 감소 |

### 2.3 rice_G281 PoC 설계 (4-6h 실험)

목적: Level 2 prototype이 Chr3:16,439,674 verdict를 보존하면서 cucumber_225 과잉
FP를 줄이는지 증명.

**Phase A (2h) — dual-anchor signal 단독 검증:**
```
minimap2 -ax sr --split-prefix=/tmp/rice_dual \
  db/Osativa_323_v7.0.fa \
  results/rice_G281/s01_qc/*_R1.fq.gz results/rice_G281/s01_qc/*_R2.fq.gz \
| samtools sort -@ 8 -o results/rice_G281/s04_host_map_dual/host_dual.bam
# Chr3:16,439,674 ±5kb 구간에서 chimeric alignment 수집
samtools view results/rice_G281/s04_host_map_dual/host_dual.bam \
  Chr3:16434674-16444674 \
| awk '$6 ~ /S/ || and($2, 0x800)' \
| ... > sites_dual_poc.tsv
```

**Phase B (2h) — PE discordancy 단독 검증:**
```python
# scripts/poc_pe_discordant.py
import pysam
bam = pysam.AlignmentFile("results/rice_G281/s04_host_map/rice_G281_host.bam")
for read in bam.fetch("Chr3", 16434674, 16444674):
    if read.is_proper_pair: continue
    if abs(read.template_length) > 1000 or read.mate_is_unmapped:
        # discordant PE candidate
        ...
```

**Phase C (1-2h) — 3-signal merge 및 verdict 비교:**
- signal_SR, signal_dual, signal_PE를 breakpoint graph로 merge (NetworkX로 간단히)
- Chr3:16,439,674 "v1.0=CAND, v1.1=?" 비교
- cucumber_line225 동일 조건으로 비교 (8 CAND → ?)

**합격 기준:** rice CAND 유지 + cucumber 225 CAND ≤5.

### 2.4 구현 공수

| 모듈 | 신설/수정 | 라인 수 추정 | 공수 |
|------|---------|------------|------|
| `scripts/s05a_site_discovery.py` | 신설 | ~400 | 5 days |
| `scripts/s05b_graph_merge.py` (신설) | 신설 | ~200 | 2 days |
| `scripts/s05_insert_assembly.py` phase 1 분리 | 수정 | -300 / +80 | 3 days |
| `tests/test_s05a_*.py` | 신설 | ~300 | 4 days |
| `config.yaml` dual-anchor preset | 수정 | ~10 | 0.5 day |
| PoC run (rice + cucumber x3) | 실행 | — | 1 day |
| **Total** | — | — | **~3 weeks** (1 dev) |

Level 2는 **Phase 2 일정 (v1.0 배포 후 3-6개월)**에 수행. v2.0 grant로 이관하지 않음.

---

## 3. Level 3 Spec — v2.0 Grant Scope (GRIDSS-style + macrosynteny + long-read)

### 3.1 아키텍처 전환

```
redgene v2.0 (proposal 구조)
┌──────────────────────────────────────────────────────────────┐
│  Data ingestion                                              │
│  • short reads (Illumina)                                    │
│  • long reads (ONT MinION / PacBio HiFi) — NEW               │
│  • reference (host + construct + element_db + synteny blocks)│
└──────────────────────────────────────────────────────────────┘
                    ↓
┌──────────────────────────────────────────────────────────────┐
│  minimap2 unified alignment                                  │
│    short: -ax sr                                             │
│    ONT:   -ax map-ont                                        │
│    HiFi:  -ax map-hifi                                       │
└──────────────────────────────────────────────────────────────┘
                    ↓
┌──────────────────────────────────────────────────────────────┐
│  Breakpoint graph construction (GRIDSS-style)                │
│    nodes: genome positions (host OR construct)               │
│    edges: read-pair support with evidence vector             │
│      [SR_count, PE_disc_count, dual_anchor_count,            │
│       long_read_split_count, local_assembly_count]           │
│    edge weight: Bayesian posterior from evidence             │
└──────────────────────────────────────────────────────────────┘
                    ↓
┌──────────────────────────────────────────────────────────────┐
│  Macrosynteny scoring (MCscan)                               │
│    input: pre-computed host-vs-related_species synteny block │
│    per candidate site: synteny_score                         │
│       = breaks_block ×2 + inside_gene ×1 + gene_order_break  │
└──────────────────────────────────────────────────────────────┘
                    ↓
┌──────────────────────────────────────────────────────────────┐
│  Local assembly (optional, only for ambiguous breakpoints)   │
│    short-read: SPAdes per-bp graph traversal                 │
│    long-read:  Flye local assembly                           │
│    Assembly가 bp-resolution 요구일 때만 동작 (skip for 대부분)│
└──────────────────────────────────────────────────────────────┘
                    ↓
┌──────────────────────────────────────────────────────────────┐
│  Annotation (LAST + element_db + curated NCBI subset)        │
│    LAST custom matrix: 75-85% divergence에 민감               │
│    element_db build-time dedup (cd-hit-est)                  │
└──────────────────────────────────────────────────────────────┘
                    ↓
┌──────────────────────────────────────────────────────────────┐
│  Verdict — multi-evidence posterior                          │
│    CANDIDATE_HIGH:   breakpoint posterior ≥ 0.9              │
│                   AND element_db hit                         │
│                   AND synteny_score ≥ 2                      │
│    CANDIDATE:        breakpoint posterior ≥ 0.7              │
│                   AND (element_db hit OR synteny_score ≥ 3)  │
│    CAND_LOW_CONF:    posterior 0.5-0.7, 1 supporting signal │
│    FALSE_POSITIVE:   host-only evidence                      │
│    UNKNOWN:          posterior < 0.5, no supporting signal   │
└──────────────────────────────────────────────────────────────┘
```

### 3.2 핵심 novelty (grant selling point)

1. **Multi-signal breakpoint graph + macrosynteny fusion** — 현존 SV caller (GRIDSS,
   Manta, SvABA, DELLY, SvABA)는 macrosynteny signal을 **사용하지 않는다**. Plant
   GMO detection처럼 host annotation quality가 높은 domain에서 이 orthogonal signal을
   결합하는 것은 novel.
2. **Long-read + short-read fusion verdict** — ONT/HiFi 혼합 데이터에서 각 signal의
   posterior를 unified Bayesian framework로 결합. NanoSV, Sniffles2는 long-read only,
   Manta는 short-only. Hybrid verdict는 drafting 단계.
3. **Element DB curation pipeline (build-time governance)** — cd-hit-est dedup +
   hierarchical source tagging + KCGP/ISO 21569 nomenclature mapping. 이건 alogrithm보다
   infrastructure contribution이지만, regulatory-grade GMO tool의 요구사항.

### 3.3 grant 공수 추정

| 항목 | 공수 | 비고 |
|------|------|------|
| Level 2 완료 (선행 요구) | 3주 | v1.1로 별도 |
| Breakpoint graph 구현 (NetworkX) | 6주 | Bayesian edge scoring 포함 |
| MCscan synteny pre-compute + scoring | 3주 | 4 host 대응 |
| Long-read pipeline (ONT/HiFi 각 1 test sample) | 4주 | 데이터 수급 포함 |
| LAST annotation switch | 2주 | 기존 blastn path 유지 (dual option) |
| Test data 확장 (ONT rice_G281, HiFi tomato_A2_3) | 2주 | 외부 sequencing |
| Validation + paper | 6주 | Bioinformatics 또는 Nat Methods |
| **Total** | **~6 months (1 dev) / ~3 months (2 devs)** | |

---

## 4. MCscan pre-mask BED Workflow (Round 2 Q4 확장)

### 4.1 목적

gmo_expert §3의 per-host FP blacklist (rice + P-Act1 100%, corn + P-Ubi1 100%, soybean +
EPRV 잔재, tomato + Solanum pararetrovirus 조각)를 **자동 생성 BED**로 공급. v1.0의
runtime BLAST filter는 per-sample마다 재계산이라 overhead가 있고, **build-time 1회
생성**으로 옮긴다.

### 4.2 Two-workflow architecture

**Workflow A: element-vs-host 100%-ortholog (MVP-able, P6로 Round 2 제안)**

```bash
# 모든 host × element pair에 대해 1회
for HOST in Osativa_323_v7.0 SLM_r2.0.pmol CucSat_B10v3 Zm_B73_v5 Gmax_v4.0; do
  blastn -task megablast \
    -query element_db/gmo_combined_db.fa \
    -db db/${HOST}.fa \
    -outfmt "6 qseqid sseqid sstart send pident length" \
    -perc_identity 98 -evalue 1e-50 \
    -out db/masks/${HOST}_element100pct.tsv

  awk 'BEGIN{OFS="\t"} $5>=98 && $6>=100 {
    s=($3<$4)?$3:$4; e=($3<$4)?$4:$3;
    print $2, s-500, e+500, $1"|"$5"%"
  }' db/masks/${HOST}_element100pct.tsv \
    | sort -k1,1 -k2,2n \
    | bedtools merge -c 4 -o distinct \
    > db/masks/${HOST}_element100pct.bed
done
```

**소요:** host당 5-15분 (rice 374 Mbp에서 ~5분, soybean 1.1 Gbp에서 ~15분). 1회 build
후 git-ignored artifact. s05에서 bedtools intersect로 soft-clip site flag.

구현 라인: 50-80줄 Python (`scripts/build_element_masks.py`) + s05에서 20줄 integration.
**이번 주 P6로 채택 권장.**

**Workflow B: MCscan lineage-specific insertion mask (v2.0)**

```bash
# 1. Host + related species MCscan ortholog
python -m jcvi.formats.gff bed --type=mRNA --key=Name \
  db/Gmax_v4.0.gff3 > db/synteny/Gmax.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name \
  db/Phaseolus_vulgaris_442_v2.1.gff3 > db/synteny/Pvulgaris.bed

# CDS fasta 추출 후 MCscanX-style ortholog search
python -m jcvi.compara.catalog ortholog Gmax Pvulgaris --no_strip_names

# 2. Anchors를 genome 좌표로 tiling, macrosynteny block 식별
python -m jcvi.compara.synteny screen \
  --minspan=30 --simple \
  Gmax.Pvulgaris.anchors Gmax.Pvulgaris.anchors.new

# 3. Gap (bridged but no orthologous content) = lineage-specific insertion 후보
python scripts/synteny_gap_to_bed.py \
  --anchors db/synteny/Gmax.Pvulgaris.anchors.simple \
  --bed-a   db/synteny/Gmax.bed \
  --bed-b   db/synteny/Pvulgaris.bed \
  --min-gap 1000 \
  --output  db/masks/Gmax_lineage_specific.bed

# 4. EPRV subset 필터 (viral protein motif)
# gmo_expert가 큐레이션 (soybean mosaic virus, CaMV-like LTR 등)
python scripts/curate_eprv_subset.py \
  --input  db/masks/Gmax_lineage_specific.bed \
  --viral-motif-db db/viral_motifs.fa \
  --output db/masks/Gmax_EPRV.bed
```

**기대 출력 예시 (soybean):**
- `Gmax_lineage_specific.bed` ~500-2000 regions (large, includes real insertions + EPRV
  + species-specific genes)
- `Gmax_EPRV.bed` ~10-50 regions (curated subset, viral motif 매칭)

**공수:** 4 host × 1회 build (MCscan 파이프라인 setup 1주, 실제 compute host별 30분~2h).
**v2.0 grant scope.**

### 4.3 rice_G281 PoC — synteny break signal 실측 전 예측

Rice MSU v7 Chr3 vs Sorghum Sb03 synteny block은 Phytozome에 published. Chr3:16.44
Mbp는 block 내부, local gene density ~250/Mbp.

**가설:** G281 T-DNA 2-copy head-to-head (~15kb + 36bp del)는
1. **synteny-break signal: 약** — natural block break가 ~100kb 주기로 발생, 15kb
   insertion은 SNR ~2-3x. 단독으로는 discriminative 아님.
2. **element_db hit 조합 signal: 강** — element hit (bar, nptII, 35S 등) + synteny
   break는 joint posterior에서 1.5-2x 증폭.

**따라서 v2.0 verdict logic은 synteny_score를 standalone signal이 아닌 element_db
evidence와 joint로 쓴다** (§3.1 verdict table의 "AND synteny_score" 조항이 이 구조).

### 4.4 Cucumber 특수성 경고

Cucumber B10v3 (8,035 contigs)은 chromosomal-level 조립이 아니라 contig-level.
MCscan chromosomal synteny block은 **적용 불가**. Cucumber는 per-contig gene
collinearity vs melon으로만 계산 가능하고, signal은 rice/soybean 대비 약함.
cucumber_line225의 8 CAND 과잉은 **Workflow B로는 해결 안 됨** — Workflow A
(element 100%-ortholog mask) + gmo_expert element DB expansion 의 조합이 정답.

---

## 5. k-mer 실험 설계 (Round 2 Q3 보강)

### 5.1 측정 대상

`StrandAwareSeedExtender.k` (현 15) 및 `recruit_by_kmer.k` (현 25)를 k ∈ {15, 17, 19,
21, 25, 31} 에서 측정.

### 5.2 측정 metric

1. **Unique k-mer ratio** — 각 host genome에서 k-mer multiplicity 분포.
2. **Per-site wall time** — `StrandAwareSeedExtender.extend()` + `recruit_by_kmer()`
   누적 시간.
3. **GT recall** — rice_G281 Chr3:16,439,674 및 cucumber_line224 LKUO03001512.1:581,328에
   대해 CAND 복원 여부.
4. **Cycle detection trigger rate** — `growth_history` cyclic 감지 발동 횟수.

### 5.3 실험 스크립트 (설계, ~2h 실험)

```bash
#!/bin/bash
# scripts/experiments/kmer_sweep.sh
set -euo pipefail
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene

mkdir -p docs/superpowers/runs/kmer_sweep/

# Step 1: jellyfish로 host별 unique k-mer ratio 측정
for HOST in Osativa_323_v7.0 SLM_r2.0.pmol CucSat_B10v3 Gmax_v4.0; do
  for K in 15 17 19 21 25 31; do
    jellyfish count -m $K -s 200M -t 16 -C \
      -o /tmp/${HOST}_k${K}.jf db/${HOST}.fa
    jellyfish histo /tmp/${HOST}_k${K}.jf \
      > docs/superpowers/runs/kmer_sweep/${HOST}_k${K}_histo.tsv
    # unique = frequency=1 count / total distinct
    TOTAL=$(awk '{s+=$2} END{print s}' docs/superpowers/runs/kmer_sweep/${HOST}_k${K}_histo.tsv)
    UNIQ=$(awk '$1==1{print $2; exit}' docs/superpowers/runs/kmer_sweep/${HOST}_k${K}_histo.tsv)
    echo -e "${HOST}\t${K}\t${UNIQ}\t${TOTAL}\t$(bc <<< "scale=3; $UNIQ/$TOTAL")" \
      >> docs/superpowers/runs/kmer_sweep/summary.tsv
    rm /tmp/${HOST}_k${K}.jf
  done
done

# Step 2: per-site wall time + GT recall (2 samples)
for SAMPLE in rice_G281 cucumber_line224; do
  for K in 15 21 25; do
    # s05 single-site run with modified k (전제: scripts/s05에 ext_k/recruit_k CLI 노출)
    python scripts/experiments/run_s05_single_site_kmer.py \
      --sample $SAMPLE \
      --site-id $(jq -r '.[0].site_id' \
        results/$SAMPLE/s05_insert_assembly/positive_sites.json) \
      --ext-k $K --recruit-k $K --max-rounds 3 \
      --output docs/superpowers/runs/kmer_sweep/${SAMPLE}_k${K}.json
  done
done

# Step 3: aggregate
python scripts/experiments/kmer_sweep_summary.py \
  --input-dir docs/superpowers/runs/kmer_sweep/ \
  --output docs/superpowers/runs/kmer_sweep/RESULTS.md
```

### 5.4 예상 결과 (이론 근거 기반, 실측 필요)

| k | rice unique ratio | soybean unique ratio | per-site time (rel to k=15) | GT recall |
|---|------------------|----------------------|---------------------------|-----------|
| 15 | ~0.55 | ~0.40 | 1.00× | baseline (5/6) |
| 17 | ~0.70 | ~0.55 | 0.90× | ≥ baseline |
| 19 | ~0.78 | ~0.65 | 0.82× | = baseline |
| 21 | ~0.85 | ~0.75 | 0.75× | = baseline (±) |
| 25 | ~0.92 | ~0.85 | 0.65× | 불확실 (low-cov sample 위험) |
| 31 | ~0.96 | ~0.92 | 0.55× | 불확실 (10x 미만 sample 실패 위험) |

**예측 gate:** k=21은 unique ratio 개선 + time 25% 절감 + recall 유지 확률 높음. **k=21
채택을 grant-level 변경으로 제안**. k=25 이상은 low-coverage 검증 필요.

**중요:** 이 측정은 v1.0 MVP 단계에서 수행할 가치 있음 (2h). 결과가 긍정적이면 Level 1
patch P7로 추가 가능. 현 plan은 "실험 전 default 변경 금지".

---

## 6. ROI vs Risk Analysis

| Level | AC 개선 | 구현 공수 | 회귀 risk | v1.0 blocking? | 권장 timing |
|------|---------|-----------|----------|---------------|------------|
| **L1 (채택된 patch: L1-3/4/5 + P6)** | AC-4 TIMEOUT 해결, AC-2 중간 개선 | 3-5 days | 낮음 | No | **이번 주 (MVP 병렬)** |
| **L1-1 cluster_window adaptive** | AC-1 low-cov 개선 | 3 days | **매우 높음** (BUG-7 class revert) | No | Phase 2 (3-6월 이후) |
| **L1-2 LAST switch option** | AC-3 divergent element annotate | 5 days | 낮음 (opt-in) | No | Phase 2 |
| **L2 (dual-anchor + PE discordancy)** | AC-1 → 100%, AC-2 cucumber 225 8→3-4 | 3 weeks | 중 (s04/s05 인터페이스 변경) | No | **Phase 2 (v1.0 배포 후 3-6월)** |
| **L3 MCscan synteny** | AC-2 production ≤2 달성 | 3 weeks | 낮음 | No | v2.0 grant |
| **L3 GRIDSS-style graph** | AC-1/AC-2 동시 향상, multi-copy 처리 강화 | 6 weeks | 높음 (architecture rewrite) | No | v2.0 grant |
| **L3 long-read integration** | AC-7 coverage robust, AC-1 5x 회수 | 4 weeks | 중 | No | v2.0 grant |
| **L3 element DB governance** | AC-3 annotation ≥90% | 2 weeks | 낮음 | No | v2.0 grant (gmo_expert 주도) |

**우선순위 요약:**
- v1.0 배포 전: L1 3-patch + P6 (1주 내).
- v1.1 (3-6월): L1-1/L1-2 + L2 (dual-anchor). 수동 검증.
- v2.0 grant (6-12월): L3 전체 + long-read.

---

## 7. Funding Pitch (Grant Proposal 1-Page)

### Title
*RedGene v2.0: Multi-signal Structural Variant Detection with Macrosynteny Integration for
Regulatory-Grade GMO/LMO Screening*

### Abstract (5 lines)
Plant GMO/LMO quarantine detection requires bp-resolution transgene insertion
characterization at production-grade specificity (per-sample FP ≤2). Current
short-read SV callers (Manta, DELLY, GRIDSS) are agnostic to host genome annotation
quality; conversely, plant reference genomes with dense synteny information
(MCscan anchors vs related species) provide untapped orthogonal evidence.
RedGene v2.0 fuses 4-signal breakpoint graph (SR + PE + dual-anchor + long-read)
with MCscan-derived macrosynteny score into a unified Bayesian verdict. On a
6-host / 20-sample Korean quarantine benchmark, we project FP reduction from
baseline ~8 CAND/sample (v1.0 MVP) to ≤2 CAND/sample while preserving 100%
ground-truth recall.

### Specific Aims
1. **Aim 1 — Multi-signal breakpoint graph**: implement GRIDSS-style Bayesian edge
   scoring integrating 4 short-read + 1 long-read signals; validate on rice/tomato/
   cucumber/soybean/corn 6-host benchmark with ONT MinION and PacBio HiFi augmentation.
2. **Aim 2 — Macrosynteny fusion**: precompute MCscan anchor graphs for 4 plant hosts
   vs phylogenetic neighbors; integrate synteny_score into verdict posterior;
   demonstrate orthogonal discrimination on heterozygous insertion events where depth
   alone is uninformative.
3. **Aim 3 — Regulatory governance**: build-time element DB curation (cd-hit-est dedup
   + KCGP/ISO 21569 nomenclature + SHA-256 lockfile) + audit trail compliant with
   APQA/KSVS/MFDS requirements; deliver production-grade packaging for field
   deployment.

### Novelty vs Prior Art
- GRIDSS/Manta/DELLY: short-read only; no synteny; no plant-specific host priors.
- NanoSV/Sniffles2: long-read only; no SR/PE fusion; no synteny.
- MCscanX/JCVI: comparative genomics; not an SV caller.
- **RedGene v2.0: first tool to fuse multi-signal SV calling with macrosynteny for a
  plant regulatory-grade application.**

### Deliverables
- v2.0 codebase (GitHub, GPL-3.0 or MIT).
- Benchmark dataset (6 host × 20 sample, short + long read, ground-truth curated).
- Paper (target: Nature Methods / Bioinformatics / Genome Biology).
- Korean regulatory (APQA/KSVS/MFDS) deployment documentation.
- 6-month validation report with external inspector user study.

### Budget (3 FTE-years)
- 1 bioinformatics dev × 12 mo: architecture + implementation.
- 1 computational biologist × 6 mo: benchmark + paper.
- 1 regulatory liaison × 6 mo: APQA/KSVS/MFDS alignment.
- ONT/HiFi sequencing: 20 samples × $200 = $4,000.
- Cloud compute (HPC augmentation): $8,000.

### Timeline
- M1-3: L2 implementation (dual-anchor + PE).
- M4-7: L3 breakpoint graph + Bayesian scoring.
- M8-9: MCscan integration.
- M10-11: Long-read pipeline (ONT + HiFi).
- M12: element DB governance + audit trail.
- M13-15: external validation + paper.
- M16-18: Korean regulatory field test + deployment.

---

## 8. 잘못 가면 안 되는 것 (Risk Management)

1. **L2/L3 전환 시 rice_G281 CAND 유지**: Chr3:16,439,674는 2-copy head-to-head라
   palindromic 구조. breakpoint graph 구현 시 palindromic edge를 node 합치는 방향의
   실수 발생 가능. PoC (§2.3)가 이 class를 필수 test case로 포함.
2. **MCscan synteny가 false positive로 변모**: natural block break (rice-sorghum ~100kb
   주기)와 insertion break 구분 threshold 미설정 시 오히려 FP 증가. gmo_expert + haibao
   joint review 필요.
3. **Long-read 샘플 미확보**: v2.0 validation은 ONT/HiFi 데이터 필요. Korean
   quarantine 원천 시퀀싱 여부 확인 (won_yim 책임). 미확보 시 공개 DB (SRA) 의존.
4. **Grant timeline slip**: 3 FTE-year는 tight. L2가 v1.1로 선행되지 않으면 v2.0 timeline
   delay. **L2를 v2.0 month 0 요구사항으로 고정**.

---

## 9. 결론

v1.0은 Korean quarantine MVP로 2026년 배포 가능. 그 이후:
- **Phase 2 (3-6월):** L2 dual-anchor + PE discordancy. v1.1 release.
- **v2.0 grant (6-18월):** L3 full redesign + macrosynteny + long-read. production-grade
  GMO/LMO detection tool 완성.

v2.0은 단순 incremental update가 아니라 plant regulatory genomics에서 **macrosynteny
+ multi-signal SV calling fusion**이라는 novel methodology 제안. grant pitch 가능.

— haibao, 2026-04-16
