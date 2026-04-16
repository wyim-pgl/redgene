# Haibao — Round 2 Cross-talk 답변 (2026-04-16)

**상대:** team-lead (won_yim PI 원칙 "MVP 통과 전 아키텍처 교체 금지" 적용)
**입장:** Level 2 rewrite는 scope out 수용. Level 1 patch 중 2개(round 축소 + DB dedup)만 사수. Level 3는 v2 grant 자료로 분리 동의.

---

## Q1. MVP vs Level 2 ROI — Level 2 기각에 동의하나

**답변: 부분 동의. Level 2 rewrite는 MVP scope out 수용.** 단 "AC-2 달성에 Level 2 필수"는 아니지만, **AC-2 근본 해결**은 Level 2 없이 어렵다는 점은 기록 남긴다.

AC-1 (≥80% recall, 현 5/6=83%):
- BUG-3/4 regression 방지는 Level 1 patch(round 축소·DB dedup)로 충분. 이미 `4d3f2bb`/`d56b18f`가 main에 들어있고 resume.md의 5/6 복구가 증거. **Level 2 없이 AC-1 통과 가능.**

AC-2 (≤5 CAND FP/sample):
- cucumber 225의 8 CAND(GT=2)는 **element_db cross-reaction + alt-locus filter(Filter D)의 threshold 문제**가 1차 원인이다. gmo_expert가 지적한 EPRV/ortholog pre-mask BED가 들어가면 6→2~3 수준으로 떨어진다고 본다. **Level 1 + gmo_expert mask만으로 AC-2 ≤5 달성 가능.** 단, "≤2 production target"까지는 Level 2의 PE discordancy signal이 있어야 확실. MVP 합격선(≤5)은 Level 1만으로 OK.
- AC-2 근본 risk: BUG-7 class(microhomology junction in low-MAPQ read)는 Level 1으로 **방어만** 되고 **탐지 개선은 안 된다**. MVP 샘플에서 이 class가 다시 터지면 patch로 막을 수 있지만, production 단계(AC-2 ≤2)에서는 Level 2 필요. → v2 grant로 미룸.

AC-4 (≤48h):
- compute_eng의 per-site SLURM array(16G/1h × 25 concurrent)가 병목 해결. **Level 2 불필요.** 단, array 병렬화 후에도 per-site 14분이 **5분으로 떨어지면** 전체가 극적으로 개선되므로 Level 1 중 `assembly round 8→3`은 **array 병렬화와 복리 효과**가 있다.

**결론:** Level 2 기각 OK. **사수: (a) assembly round 축소 (b) DB dedup**. Level 2는 v2 문서로 이관.

---

## Q2. Level 1 5-patch 우선순위 ranking

AC-1·AC-2·AC-4 기여도 × 구현 난이도 × 회귀 위험 matrix:

| # | Patch | AC-1 | AC-2 | AC-4 | 구현 난이도 | 회귀 위험 | **우선순위** |
|---|------|------|------|------|-----------|----------|-------------|
| P1 | **assembly round 8→3** (`max_rounds=3`) | 낮음 (recall 약간 하락 가능) | 중 (per-site 시간↓ → array 효과) | **매우 높음** (14min→5min, UGT72E3 24h→8h 기대) | **낮음** (1줄 수정 + stats 재측정) | 중 (s05_stats.txt로 round-3 수렴율 사전 측정 필수) | **★ 1 (이번 주)** |
| P2 | **transgene_db build-time dedup** (`cd-hit-est -c 0.95`) | 중 (BUG-3 class 근본 해결) | 중 (univec FP 자동 제거) | 낮음 | **낮음** (Makefile 1스텝) | 매우 낮음 (build time 변경, runtime 그대로) | **★ 2 (이번 주)** |
| P3 | **BLAST 호출 통합** (megablast 3→1, in-memory filter) | 낮음 | 낮음 | 중 (per-site BLAST 시간 30% 절감 추정) | **높음** (3 함수 인터페이스 변경, DataFrame 도입, test 확장) | 높음 (filter threshold 간 결과 drift 검증 필요) | 3 (Phase 2 이후) |
| P4 | **cluster_window adaptive** (coverage 기반 per-chr) | 중 (BUG-7 방어 + low-cov recall↑) | 낮음 (FP 증가 risk) | 낮음 | **중** (coverage estimation + unit test 추가) | **매우 높음** (BUG-7 revert가 바로 이 class) | 4 (MVP 후로 보류) |
| P5 | **LAST 옵션 추가** (`--aligner last`) | 낮음 | 낮음 (divergent host element 소폭↑) | 낮음 (LAST는 blastn보다 약간 빠름) | 중 (LAST install + config + wrapper) | 중 (기본값 BLAST 유지하므로 우회 옵션) | 5 (v2로 이관) |

**이번 주 채택 권장: P1 + P2.**
- P1은 AC-4 직접. `scripts/s05_insert_assembly.py:2372`의 `max_rounds=8` → `max_rounds=3` 1줄. 단 prerequisite으로 resume.md "most converge by 6"이 round-3 수렴율로 얼마나 떨어지는지 `s05_stats.txt` aggregate 측정 (30분 작업).
- P2는 AC-1/AC-2 근본. `element_db/build_common_payload.sh`에 `cd-hit-est` 1줄 추가 (이미 bioconda 내 있음).

**1개만이라면 P1.** AC-4 TIMEOUT이 MVP 합격 최대 장애물. P2는 이미 `_should_replace` runtime patch가 band-aid로 동작하므로 후속.

---

## Q3. k-mer critique 정량화 (k=15 → 21/25/31)

**솔직한 한계 고백:** 실험적 숫자는 이 session에서 측정하지 않았다. 아래는 **이론 근거 + 비슷한 plant genome 벤치마크에서의 경험치**이며 "실험 검증 필요"를 병기한다.

### 이론 (4^k 대비 genome size)

| k | 4^k | soybean 1.1 Gbp 기준 coincidence 기대치 | plant repeat(50%) 포함 unique ratio 이론치 |
|---|------|---------------------------------------|-----------------------------------------|
| 15 | 1.07e9 | 1.03× random (≈1 hit/kmer) | **~40-50%** (repetitive region dominant) |
| 21 | 4.40e12 | 0.00025× (희소) | ~75-85% |
| 25 | 1.13e15 | <0.001× | ~90% |
| 31 | 4.61e18 | 근사 unique | ~95% (실제 SPAdes 기본값 유사 범위) |

Marçais & Kingsford (Jellyfish) 및 BCALM 벤치마크에 따르면, 식물 1-2 Gbp 유전체에서
k=21 이상은 repeat 영역에서도 "mostly unique" 수준(~80%+). 현재 k=15는 **plant
context에서 거의 무조건 multi-hit이 된다**.

### per-site time 영향 추정 (불확실성 높음)

- k=15 → k=21로 올리면 `recruit_by_kmer`의 candidate read 수가 **50-70% 감소** (false
  recruitment 제거). 단, `StrandAwareSeedExtender.extend`는 k-mer path가 unique해져서
  cycle detection 분기가 줄어듦 → **per-site 시간 10-20% 감소 가능**.
- 반면 **recall은 감소할 위험**: k=21 seed가 low-coverage (≤5x) read에서는 안
  잡힘. min_depth=2 유지 시 k=21이 실제로 extension을 시작할 수 있을지는 sample별.
- **결론:** k=21은 "repeat FP 감소"에는 유리, "low-coverage recall"에는 불리. MVP
  target이 ≥10x coverage이므로 **k=21 시도 가치 있음**. k=25/31은 SPAdes-level 요구.

### 실험 계획 (권장, ~2h 작업)

```
# 각 k에 대해 rice_G281, cucumber_line225 (multi-copy), soybean_AtYUCCA6 (genome 큰
# host) 세 샘플로 single-site rerun
for k in 15 17 19 21 25; do
  python -c "from scripts.s05_insert_assembly import assemble_insert; \
             assemble_insert(..., ext_k=$k, recruit_k=$k, max_rounds=3)"
done
# metric: round-3 converged 비율, per-site wall time, GT CAND 복원율
```

**실험 결과 없이 default 변경 금지.** 현 k=15는 regression risk. 측정 전에는 "권장 변경"
포지션만 유지.

---

## Q4. MCscan synteny BED 생성 workflow (gmo_expert 연계)

gmo_expert의 제안: **EPRV / 100%-ortholog pre-mask BED**를 s05 post-assembly filter
입력으로 추가. MCscan으로 이 BED를 자동 생성하는 건 **두 가지 level**로 나뉜다.

### Workflow 1: 100%-ortholog mask (gmo_expert core request)

**목적:** host genome에 construct element와 100% 동일한 endogenous copy가 있는 좌표를
BED로 미리 빼둔다. 예: rice host에서 P-Act1 construct를 쓰면 `LOC_Os11g06390` 좌표를
사전 mask → s05 soft-clip cluster가 거기서 발견되면 즉시 FP 처리.

**MCscan을 써야 하는가?** **아니다. 이 부분은 단순 BLAST로 충분.**

```bash
# 1회 pre-compute (host + construct 조합별)
blastn -task megablast -query element_db/gmo_combined_db.fa -db db/Osativa_323_v7.0.fa \
  -outfmt "6 qseqid sseqid sstart send pident length" -perc_identity 98 -evalue 1e-50 \
  -out db/masks/rice_element_100pct.tsv
# → BED 변환 (±500bp slop 추가해서 site detection 좌표와 매칭)
awk 'BEGIN{OFS="\t"} $5>=98 && $6>=100 {
  s=($3<$4)?$3:$4; e=($3<$4)?$4:$3;
  print $2, s-500, e+500, $1"|"$5"%"
}' db/masks/rice_element_100pct.tsv | sort -k1,1 -k2,2n | bedtools merge \
  > db/masks/rice_element_100pct.bed
```

s05 수정 지점: `find_softclip_junctions` 직후 `bedtools intersect`로 site를 이 mask와
대조, hit된 site는 `confidence=low_endogenous_ortholog`로 flag. **구현 < 200줄.**

### Workflow 2: MCscan macrosynteny-based EPRV / integrated viral mask (상위 레벨)

**목적:** soybean의 soybean mosaic virus-like EPRV, tomato의 Solanum pararetrovirus 같은
**integrated pararetrovirus**는 BLAST로 잡기 어렵다 (partial, diverged). MCscan으로
closely-related species(soybean vs common bean, tomato vs potato) synteny block을
만들면, **"soybean에만 있고 common bean에는 없는 block"** = 후보 EPRV / lineage-specific
integration region이다.

```bash
# MCscan pipeline (JCVI), host별 1회 구축
python -m jcvi.formats.gff bed --type=mRNA --key=Name \
  db/Gmax_v4.0.gff3 > db/synteny/Gmax.bed
python -m jcvi.compara.catalog ortholog Gmax Pvulgaris --no_strip_names
# → Gmax.Pvulgaris.anchors 생성 (syntenic gene pair)
python -m jcvi.compara.synteny screen --minspan=30 --simple \
  Gmax.Pvulgaris.anchors Gmax.Pvulgaris.anchors.new
# Anchors를 host 좌표로 tiling, synteny block "bridged" 영역 찾기
# gap 영역 = lineage-specific insertion candidate (EPRV 후보 포함)
python scripts/synteny_gap_to_bed.py \
  db/synteny/Gmax.Pvulgaris.anchors.simple \
  db/Gmax.bed \
  --min-gap 1000 \
  > db/masks/Gmax_lineage_specific.bed
```

이 BED는 "EPRV도 포함될 수 있으나 진짜 insertion도 포함"이라 **수동 큐레이션이 필요**한
intermediate 산출물. gmo_expert가 이 BED를 검토해서 "EPRV subset만" 분리하는 게
이상적. MVP scope out 유력 — v2 infrastructure로 이관.

### rice Chr3:16,439,674 synteny break 실제 signal 있나 (추측)

**Rice MSU v7 Chr3 ↔ Sorghum bicolor Sb03 synteny block (Phytozome data)**

Rice Chr3: ~0-40 Mbp 구간은 Sorghum Sb03 ~50-58 Mbp와 synteny block으로 대응됨 (JCVI
grape/sorghum/rice comparison paper 2008의 알려진 블록). Chr3:16.44 Mbp는 **블록의
내부**, 주변 gene density 중간 (rice annotation에서 protein-coding gene ~250 per
Mbp 수준).

**예측:**
- G281의 T-DNA 2-copy head-to-head (~15kb + 36bp del)는 이 block 안에 떨어진다.
- Sorghum ortholog region에 대해 **±20kb synteny continuity**를 측정하면, insertion
  바로 양쪽 gene pair는 여전히 Sorghum에 대응하지만 **insertion 지점에서 gene
  order의 local "gap" (Sorghum 쪽 ortholog 없는 15kb foreign sequence)**가 남는다.
- MCscan anchors.simple의 block score가 해당 좌표에서 **국소 break**를 보이는 것이
  기대되는 signal.

**정량적으로 얼마나 강한가:** Rice와 Sorghum은 **50 MYA 분화 + polyploid WGD 영향**으로
synteny block의 natural break가 ~100kb마다 존재한다. 15kb insertion break는 natural
background 대비 **signal-to-noise 2-3x** 정도. 즉 synteny_score가 discriminative이긴
하되 standalone signal로는 약함. **element_db hit + synteny break 조합**으로만 유용.

### MVP scope 권고

MCscan full pipeline은 MVP에서 제외. 단 **Workflow 1 (element 100%-ortholog BLAST
pre-mask)** 은 gmo_expert의 core 요구이고, **BLAST 1회 + BED conversion 50줄로 이번
주 안에 구현 가능**. 이것만 Level 1 patch 목록에 "P6"로 추가 권장:

**P6 (신설): element vs host 100%-ortholog pre-mask BED**
- AC-2 기여도: **높음** (cucumber 8 CAND 과잉의 상당 부분 제거 기대)
- 구현 난이도: 낮음 (BLAST + bedtools + 20줄 Python)
- 회귀 위험: 매우 낮음 (mask된 site는 low-confidence로만 downgrade, 완전 drop 아님)
- 우선순위: **P2(DB dedup)와 공동 2위**

gmo_expert가 EPRV mask까지 원하면 그 부분은 Workflow 2로 v2 이관.

---

## 정리 — Haibao가 Round 2에서 방어하는 scope

| 항목 | 입장 |
|------|------|
| Level 1 P1 (assembly round 8→3) | **이번 주 채택 주장** (P1 → AC-4 직접 해결) |
| Level 1 P2 (transgene_db dedup) | **이번 주 채택 주장** (P2 → AC-1/AC-2 근본) |
| Level 1 P6 신설 (element 100%-ortholog BED) | **이번 주 채택 주장** (gmo_expert 연계, AC-2 기여) |
| Level 1 P3/P4/P5 | MVP 후 Phase 2로 이관 수락 |
| Level 2 rewrite (dual-anchor + PE + multi-signal) | **scope out 수용** — v2 grant 문서로 분리 |
| Level 3 full redesign + MCscan synteny | v2 grant 자료로 이관. MVP scope out 수용 |
| k-mer k=15 → 21 | **실험 측정 후 결정** 보류. default 변경 금지 |

won_yim의 "MVP 통과 전 아키텍처 교체 금지" 원칙과 충돌하지 않는 범위에서 **P1 + P2 +
P6**만 Level 1 채택 주장. Round 3 결정 회의에서 이 3개 patch에 대해 bio_king/compute_eng
입장 확인 요청.

— Haibao, 2026-04-16
