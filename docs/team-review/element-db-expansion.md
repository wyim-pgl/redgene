# Element DB Expansion Plan

**작성자:** gmo_expert
**날짜:** 2026-04-16
**Cycle:** redgene-review-2026-04-16 Round 3 (MVP 확정, Phase 2 로드맵)
**관련:** `bug.md` BUG-4/8/9, `resume.md` soybean 0-CANDIDATE, Round 1/2 notes, won_yim `team-consensus.md`
**Status:** MVP deliverable, Round 4 검토 대기

**won_yim Round 3 synthesis 반영 (2026-04-16 업데이트):**
- **AtYUCCA6 MVP 승격 확정** (Round 2 유보 → Round 3 채택): P1 batch에 정식 포함, W2 workstream
- **P6 pre-mask BED MVP 채택** + **soybean × AtYUCCA6 75% ortholog 필수 포함** (AtYUCCA6 복원의 안전장치)
- **Cas9 `--cas9-info-only` flag 확정** (site count 폭증 방지, CANDIDATE 승격 트리거 제외)
- **gmo_expert 담당 workstream:** W2 (P0 6 + P1 8 seqs 추가) + W7 (pre-mask BED rationale 수동 큐레이션)

---

## 목차

1. [MVP 필수 확장 리스트 (15 MUST-HAVE)](#1-mvp-필수-확장-리스트-15-must-have)
2. [build_common_payload.sh 수정 patch (BUG-9 완전 해결)](#2-build_common_payloadsh-수정-patch-bug-9-완전-해결)
3. [Canonical-triplet rule config 스키마](#3-canonical-triplet-rule-config-스키마)
4. [Per-host FP blacklist BED workflow](#4-per-host-fp-blacklist-bed-workflow)
5. [DB governance (regulatory R-3 MUST)](#5-db-governance-regulatory-r-3-must)
6. [Phase 2 로드맵](#6-phase-2-로드맵)

---

## 1. MVP 필수 확장 리스트 (15 MUST-HAVE)

Round 2 Q1 답변에서 확정한 subset. CD-HIT 절감 효과가 낮으므로 (서열이 기능적으로 독립), regulatory MUST-HAVE 필터로 자체 축소. 총 추가 bp ≈ 17 kb (현 DB 121 kb 대비 +17%).

### P0 — CRISPR 및 승인 이벤트 marker (6개)

| Name | Accession | 좌표 | Size | Efetch 한 줄 |
|---|---|---|---|---|
| **SpCas9 (plant codon-opt)** | Addgene pX330 합성 CDS (MH782573.1 참조) | 1-4104 | ~4.1 kb | `efetch -db nuccore -id MH782573.1 -format fasta -seq_start 1 -seq_stop 4104` |
| **sgRNA scaffold (tracrRNA-derived, universal 80bp)** | synthetic — `GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT` | — | 80 bp | `printf ">regulatory\|sgRNA_scaffold\|synthetic\|SpCas9_universal\n<seq>\n" > sgRNA_scaffold.fa` |
| **PMI (manA, *E. coli*)** | Z11946.1 | 3-1177 | ~1.2 kb | `efetch -db nuccore -id Z11946.1 -format fasta -seq_start 3 -seq_stop 1177` |
| **vip3Aa** | KX344421.1 | 1-2370 | ~2.4 kb | `efetch -db nuccore -id KX344421.1 -format fasta` |
| **cry2Ab2** | AY007234.1 | 1-1899 | ~1.9 kb | `efetch -db nuccore -id AY007234.1 -format fasta` |
| **cry34Ab1 + cry35Ab1** | AY126452.1 + AY126453.1 | full | ~0.45+1.1 kb | `for A in AY126452.1 AY126453.1; do efetch -db nuccore -id $A -format fasta; done` |

**주의:** SpCas9는 Phase 1.5 positive site 수를 증가시킬 가능성 (tomato_Cas9_A2_* 샘플에서 construct 내부 중복 hit). 대응: `--cas9-info-only` flag로 CANDIDATE 승격 트리거에서 제외, 대신 `element_annotation.tsv`에 informational tag로만 기록 — won_yim Round 3에서 MVP 확정. bio_king s05 refactor에서 구현 (W8 workstream).

**Workstream mapping (won_yim Round 3):**
- **W2** = 이 섹션(§1) 15 seq 확장 실행 + §2 build_common_payload patch + §5 governance Makefile. **owner: gmo_expert.**
- **W7** = §4 pre-mask BED rationale 수동 큐레이션 (자동 생성된 BED에 rationale_code 및 citation 부여). **owner: gmo_expert.**
- **W8** = `--cas9-info-only` flag 구현 + §3 canonical triplet rule 코드화. **owner: bio_king.**

### P1 — 격리소 수입 이벤트 및 soybean 0-CAND 해결 (8개)

| Name | Accession | 좌표 | Size | 비고 |
|---|---|---|---|---|
| **T-pinII** | D00063.1 | 1-250 | ~250 bp | `element_summary.tsv`의 "T-PinII curated" entry에 실체 서열 채움 |
| **T-g7** | AF242881.1 | 1-250 | ~250 bp | Debode 2013 reference |
| **CTP2 (AtEPSPS transit peptide)** | AT1G48860 (TAIR11) or X64538.1:1-216 | 1-216 | 216 bp | GTS 40-3-2 RR soybean 필수 |
| **2mEPSPS (maize double-mutant)** | X63374.1 또는 이와 동등한 event-derived subregion | 1-1340 | ~1.3 kb | GA21 이벤트 |
| **AtYUCCA6** | NM_122473.3 | 1-1100 | ~1.1 kb | **MVP 채택 확정** (won_yim Round 3 synthesis). soybean × AtYUCCA6 host ortholog 75% 매칭 위험은 §4 P6 pre-mask BED로 격리 — 복원 선결 조건 |
| **UGT72E3** | 적절한 Arabidopsis/합성 CDS (KM009065.1 등) | full | ~1.4 kb | AtYUCCA6와 동일 정책. Gm UGT family 교차반응도 §4 BED 큐레이션 범위 |
| **P-nos subregion** | V00087.1 | 1847-2113 | 267 bp | **BUG-9 완전 해결** — `element_summary.tsv` 좌표 |
| **T-nos subregion** | V00087.1 | 1277-1536 | 260 bp | 동일 |

### P2 — CRL 호환성 최소 (1 항목, 최대 12 seqs)

| Name | Source | 비고 |
|---|---|---|
| **CRL unique subset** | `element_db/crl_amplicons.fa` 82개 중 EUginius 131 서열과 non-redundant한 ~10-12개 | Section 5.1 dedup Makefile로 자동 추출 |

**왜 subset만?** crl_amplicons.fa 82개 중 QL/QT prefix가 이미 EUginius와 겹치는 것이 대다수 (AX033493, CRL-GMOMETHODS 표기의 동일 amplicon). 실질 unique는 10~15개 추정. 완전 통합은 Phase 2.

---

## 2. build_common_payload.sh 수정 patch (BUG-9 완전 해결)

### 문제
현행 `build_common_payload.sh`는 `efetch -db nuccore -id V00087.1 -format fasta` 로 **full accession 8kb** 서열을 받아옴. `common_payload_manifest.tsv`에 `V00087.1	P-nos	...` 한 줄만 있으므로, 받아온 서열은 **nos operon 전체(P-nos + T-nos + CDS 포함)**. P-nos로 태깅되었으나 T-nos region의 BLAST hit도 같은 entry로 반환 → silent 오분류. BUG-9의 `b61b7f1` commit은 T-nos row 제거만 했고 P-nos의 region 정제는 미완.

### Patch 설계
`common_payload_manifest.tsv`를 4-column으로 확장하여 optional `seq_start` / `seq_stop` 지원.

**수정된 manifest (`common_payload_manifest.tsv`):**
```tsv
accession	purpose	notes	subregion
X17220.1	bar	phosphinothricin acetyltransferase, Streptomyces hygroscopicus	260-811
V00618.1	nptII	neomycin phosphotransferase II, Tn5 transposon	1-795
M55269.1	hpt	hygromycin phosphotransferase, Escherichia coli	113-1138
X06788.1	gusA	beta-glucuronidase reporter, E. coli (uidA)	1-1812
U55762.1	gfp	green fluorescent protein, jellyfish	
L29345.1	egfp	enhanced GFP variant	
V00141.1	P-CaMV35S	Cauliflower mosaic virus 35S promoter	6909-7439
V00087.1	P-nos	Agrobacterium tumefaciens nopaline synthase promoter	1847-2113
V00087.1	T-nos	Agrobacterium tumefaciens nopaline synthase terminator	1277-1536
X04879.1	T-ocs	octopine synthase terminator	5300-5700
```

**수정된 `build_common_payload.sh` (핵심 부분):**
```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

MANIFEST=common_payload_manifest.tsv
OUT=common_payload.fa

TMPOUT=$(mktemp)
trap 'rm -f "$TMPOUT"' EXIT

# Skip header row via NR>1 in awk feed
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r acc purpose notes subregion; do
    [[ -z "$acc" ]] && continue
    header="${purpose}|${acc}"

    # Build efetch args; include subregion only when specified
    efetch_args=(-db nuccore -id "$acc" -format fasta)
    if [[ -n "${subregion:-}" ]]; then
        seq_start="${subregion%%-*}"
        seq_stop="${subregion##*-}"
        efetch_args+=(-seq_start "$seq_start" -seq_stop "$seq_stop")
        header="${header}:${subregion}"
    fi

    echo "Fetching $acc ($purpose) ${subregion:+region $subregion}..." >&2
    efetch "${efetch_args[@]}" \
        | awk -v h=">${header}" 'NR==1{print h; next}{print}' \
        >> "$TMPOUT"
    sleep 0.4
done

mv "$TMPOUT" "$OUT"
echo "Wrote $OUT with $(grep -c '^>' "$OUT") sequences"
```

### 검증 체크리스트
- [ ] 생성된 `common_payload.fa`에 11개 entry 존재 (기존 9 + T-nos 복원 + T-ocs subregion)
- [ ] P-nos entry 길이 = 267 bp (= 2113 - 1847 + 1)
- [ ] T-nos entry 길이 = 260 bp
- [ ] Headers: `>P-nos|V00087.1:1847-2113`, `>T-nos|V00087.1:1277-1536`
- [ ] 기존 rice/tomato validation 재실행해서 Phase 3 `element_annotation.tsv`에서 P-nos vs T-nos 분리 정확도 확인

---

## 3. Canonical-triplet rule config 스키마

### 배경
Round 1 §5 (soybean 0-CAND), Round 2 Q3 (bio_king priority) 연계. won_yim/bio_king 합의: "bar + P-CaMV35S + T-ocs canonical triplet이 insert 내부에 함께 존재하면 CANDIDATE로 강제 승격" rule을 config-driven으로.

### config.yaml 추가 스키마

```yaml
# 파일: config.yaml (기존 pipeline: 섹션 아래)
pipeline:
  threads: 8
  memory_gb: 32
  common_payload_db: element_db/common_payload.fa
  # NEW: canonical triplet 기반 CANDIDATE 승격 rule
  canonical_triplets:
    # 기본 정책: bar + P-35S + T-nos/T-ocs 3종이 insert에 같이 있으면 무조건 CANDIDATE
    default:
      required_any:
        - [bar, P-CaMV35S, T-nos]        # any-of 조합 1
        - [bar, P-CaMV35S, T-ocs]        # any-of 조합 2
        - [nptII, P-CaMV35S, T-nos]      # any-of 조합 3
      min_identity: 0.90
      min_coverage: 0.80
      min_hits_per_element: 1

  # 샘플별 override (construct-specific canonical marker set)
  sample_canonical_triplets:
    rice_G281:
      required_any:
        - [hLF1, P-Gt1, T-nos]           # G281 고유 payload
      min_identity: 0.90
      min_coverage: 0.80
    soybean_AtYUCCA6:
      required_any:
        - [AtYUCCA6, bar, T-nos]         # K-GMO 승인 이벤트 spec
        - [bar, P-CaMV35S, T-ocs]        # 위 실패 시 fallback canonical
      min_identity: 0.85                 # plant-derived는 85%로 완화
      min_coverage: 0.70
    tomato_Cas9_A2_3:
      required_any:
        - [SpCas9, P-CaMV35S, T-nos]     # Cas9-present LMO spec
      min_identity: 0.90
      min_coverage: 0.80
      cas9_info_only: false              # A2_3에서는 Cas9를 CANDIDATE 트리거로 인정

  # Verdict 정책
  verdict_rules:
    auto_cand_via_triplet: true          # triplet 만족 시 FALSE_POSITIVE → CANDIDATE 승격
    auto_cand_via_element_db_only: true  # 기존 behavior 유지 (element_db hit 있으면 CAND)
    chimeric_filter: true                 # Filter C 유지
    host_fraction_max: 0.80              # Filter A
```

### 구현 스케치 (bio_king에 위임)
```python
# scripts/s05_insert_assembly.py 일부
def apply_canonical_triplet(site, element_hits, config, sample_name):
    """Return True if site should be promoted to CANDIDATE by triplet rule."""
    sample_rule = config.get("sample_canonical_triplets", {}).get(sample_name)
    default_rule = config["pipeline"]["canonical_triplets"]["default"]
    rule = sample_rule or default_rule
    
    present_markers = {h["name"] for h in element_hits
                       if h["identity"] >= rule["min_identity"]
                       and h["coverage"] >= rule["min_coverage"]}
    
    for required_set in rule["required_any"]:
        if set(required_set).issubset(present_markers):
            site["triplet_matched"] = required_set
            return True
    return False
```

### 안전장치
- Triplet rule은 **기존 element_db hit 기반 CAND 판정과 OR 결합**, AND 아님. 즉 rule이 failed여도 element_db 방식으로 CAND 될 수 있음.
- `cas9_info_only: true`이면 SpCas9 hit은 존재는 기록하되 triplet 구성원으로 인정 안 함 (false discovery 억제).
- 샘플별 override가 우선, 없으면 default.

---

## 4. Per-host FP blacklist BED workflow

### 설계 원칙 (Round 2 Q2 + won_yim Round 3 synthesis 확정)

| 대상 | Pre-mask (s02 전) | 사후 filter (Filter D) | 비고 |
|---|---|---|---|
| 100% ortholog (corn × P-Ubi1, rice × P-Act1, tomato × P-TA29) | ✅ Pre-mask | — | 기본 100% 매칭군 |
| **soybean × AtYUCCA6 (~75% ortholog, GmYUC family 10+ copies)** | ✅ **Pre-mask (MUST)** | — | **AtYUCCA6 MVP 복원의 필수 안전장치** (won_yim Round 3) — mask 없이 복원 시 systemic FP |
| **soybean × UGT72E3 (Gm UGT family)** | ✅ Pre-mask | — | UGT72E3 MVP 복원의 동일 안전장치 |
| EPRV / soybean pararetrovirus 잔재 | — | ✅ Filter D 유지 | **pre-mask 기각** (insertion hotspot 가능성, team-lead 지시) |
| 60~80% 유사 영역 | — | ✅ 기존 s03b + Filter D | 기존 로직 유지 |

**Pre-mask MUST 조건 요약:**
1. **100% self-ortholog** (corn×P-Ubi1, rice×P-Act1): 무조건 pre-mask
2. **Plant-derived MVP payload × host family** (soybean×AtYUCCA6, soybean×UGT72E3): AtYUCCA6/UGT72E3 **복원과 pre-mask는 패키지 딜**. 둘 중 하나만 실시 금지. 70-80% ortholog를 mask하지 않으면 Round 1 §5 H1이 재발.

**기각 항목:** EPRV의 pre-mask는 Filter D로 전담. 나머지는 기존 s03b+Filter D로 유지.

### 자동 BED 생성 스크립트 설계

**`scripts/build_host_endogenous_bed.sh`** (신규):
```bash
#!/usr/bin/env bash
# Generate host-endogenous FP blacklist BED from element_db × host BLAST
# Output: results/host_masks/{host}_endogenous.bed + rationale.tsv
set -euo pipefail

HOST_FA="$1"            # e.g. db/Osativa_323_v7.0.fa
HOST_NAME="$2"          # e.g. rice
ELEMENT_DB="${3:-element_db/gmo_combined_db.fa}"
OUTDIR="${4:-results/host_masks}"
MIN_PIDENT="${5:-98}"   # 98% identity threshold (near-100% ortholog only)
MIN_LEN="${6:-100}"     # ≥100 bp alignment

mkdir -p "$OUTDIR"
TSV="$OUTDIR/${HOST_NAME}_endogenous.tsv"
BED="$OUTDIR/${HOST_NAME}_endogenous.bed"
RATIONALE="$OUTDIR/${HOST_NAME}_endogenous_rationale.tsv"

# Step 1: BLAST element_db → host, retain near-perfect matches
blastn -query "$ELEMENT_DB" -subject "$HOST_FA" \
       -outfmt "6 qseqid sseqid pident length qlen slen sstart send evalue bitscore" \
       -perc_identity "$MIN_PIDENT" -evalue 1e-20 \
  | awk -v m="$MIN_LEN" '$4 >= m' > "$TSV"

# Step 2: Convert to BED (sorted, 0-based)
awk 'BEGIN{OFS="\t"} {
       start = ($7 < $8) ? $7-1 : $8-1;
       end   = ($7 < $8) ? $8   : $7;
       print $2, start, end, $1"|pident="$3"|len="$4, $10, "+"
     }' "$TSV" | sort -k1,1 -k2,2n > "$BED.raw"

# Step 3: Merge adjacent intervals within 100bp
bedtools merge -i "$BED.raw" -d 100 -c 4,5 -o distinct,mean > "$BED"
rm -f "$BED.raw"

# Step 4: Rationale TSV (for regulatory audit trail)
{
  echo -e "host_region\telement_matched\tpident\talignment_len\trationale_code\tcitation"
  awk 'BEGIN{OFS="\t"} {
         split($4, a, "|");
         print $1":"$2"-"$3, a[1], a[2], a[3], "E-01", "auto-detected via element_db BLAST ≥98% identity"
       }' "$BED"
} > "$RATIONALE"

echo "Wrote $BED ($(wc -l < $BED) regions) and $RATIONALE"
```

### 수동 큐레이션 단계 (regulatory audit trail)

`host_masked_rationale.tsv`는 자동 생성되나, **각 구간에 대해 rationale_code를 수동으로 업그레이드**:

| rationale_code | 의미 | 근거 문헌 |
|---|---|---|
| E-01 | Auto-detected 98%+ (unverified) | (auto) |
| E-02 | Known endogenous homolog (curated) | e.g. CLAUDE.md Known Pitfalls, Debode 2013 |
| E-03 | 100% match — same organism gene (self) | e.g. corn × P-Ubi1-maize |
| E-04 | Published reference-gene coordinate | e.g. Christensen & Quail 1996 |
| X-01 | Reviewer rejected mask | (인접 insertion hotspot 가능성으로 제외) |

### 4대 host별 priority pair

| Host | Element | 좌표 출처 | rationale_code | 예상 BED 레코드 수 | MVP 상태 |
|---|---|---|---|---|---|
| **corn (Zm_B73_v5)** | P-Ubi1-maize (S94464.1) | `build_host_endogenous_bed.sh corn` 자동 + curation | E-03 | 1-3 (ZmUbi1 locus만) | MVP |
| **rice (Osativa_323_v7.0)** | P-Act1-rice (S44221.1) | 동일 | E-03 | 1-2 (LOC_Os11g06390) | MVP |
| **rice** | P-Ubi1-maize | 자동 | E-02 | ~5-10 (Os Ubi family) | MVP |
| **tomato (SLM_r2.0)** | P-TA29 (X52283.1) | 자동 | E-02 | 2-5 (Solanum tapetum) | MVP |
| **soybean (Gmax_v4.0)** | **AtYUCCA6 (NM_122473.3)** | 자동 + 수동 큐레이션 | E-02 | 10+ (Gm YUC family) | **MVP-BLOCKING** — AtYUCCA6 복원의 선행 조건. 미완료 시 AtYUCCA6를 element_db에 재투입 금지 |
| **soybean** | **UGT72E3 (Gm UGT family)** | 자동 + 수동 큐레이션 | E-02 | 15+ (Gm UGT paralog) | **MVP-BLOCKING** — UGT72E3 복원의 선행 조건 |
| **soybean** | GmUbi family × P-Ubi1 | 자동 | E-02 | ~5 | MVP |
| **cucumber (CucSat_B10v3)** | (낮은 우선순위, dicot 공통만) | 자동 | E-02 | ~3 | Phase 2 |

**MVP-BLOCKING 의미:** 해당 BED 구간이 큐레이션·validate 되기 전까지 §1 P1 `AtYUCCA6` / `UGT72E3` entry를 `element_db/gmo_combined_db.fa`에 merge하지 않음. 순서: (1) `build_host_endogenous_bed.sh soybean` 실행 → (2) rationale.tsv 수동 큐레이션 → (3) s02 integration 검증 → (4) 그 후에야 AtYUCCA6/UGT72E3 efetch + DB merge. gmo_expert 단독 workstream (W2+W7 동시 진행, 내부 순서 고정).

### s02 BWA 단계 integration

```python
# run_pipeline.py — step 2 (host BWA)
host_mask_bed = Path(f"results/host_masks/{sample.host_name}_endogenous.bed")
if host_mask_bed.exists():
    # Post-BWA: filter out reads mapped entirely within mask
    cmd_pipe = (
        f"bwa mem -t {threads} {host_ref} {r1} {r2} | "
        f"samtools view -b -L <(complement {host_mask_bed}) - | "
        f"samtools sort -@ {threads} -o {host_bam}"
    )
    # And record masked-region soft-clip reads separately as FALSE_NEGATIVE_MASKED
    mark_masked_junctions(host_bam, host_mask_bed, 
                          out=f"{outdir}/s02b_false_negative_masked.tsv")
```

`FALSE_NEGATIVE_MASKED` 태그 별도 보고 → regulatory audit trail 유지.

---

## 5. DB governance (regulatory R-3 MUST)

### 5.1 dedup Makefile 스케치

**`element_db/Makefile`** (신규):
```makefile
# Combined DB build + dedup
CD_HIT_EST := /data/gpfs/home/wyim/scratch/bin/cd-hit-v4.6.8-2017-1208/cd-hit-est
DEDUP_ID   := 0.95
MIN_LEN    := 50

GMO_COMBINED := gmo_combined_db.fa
COMMON       := common_payload.fa
CRL          := crl_amplicons.fa
NEW_P0P1     := new_p0p1_elements.fa        # Round 3 §1 결과물

UNIFIED_RAW  := gmo_all_unified_raw.fa
UNIFIED      := gmo_all_unified.fa
UNIFIED_MD5  := gmo_all_unified.fa.md5

.PHONY: all dedup verify clean

all: $(UNIFIED) $(UNIFIED_MD5)

$(UNIFIED_RAW): $(GMO_COMBINED) $(COMMON) $(CRL) $(NEW_P0P1)
	cat $^ > $@

$(UNIFIED): $(UNIFIED_RAW)
	$(CD_HIT_EST) -i $< -o $@ -c $(DEDUP_ID) -n 10 -M 4000 -T 4 -l $(MIN_LEN)
	@echo "Dedup: $$(grep -c '^>' $<) → $$(grep -c '^>' $@) sequences"

$(UNIFIED_MD5): $(UNIFIED)
	md5sum $< | awk '{print $$1}' > $@
	date -u +"%Y-%m-%dT%H:%M:%SZ" > $@.build_date

verify: $(UNIFIED)
	makeblastdb -in $< -dbtype nucl
	bwa index $<
	samtools faidx $<

clean:
	rm -f $(UNIFIED_RAW) $(UNIFIED) $(UNIFIED).clstr $(UNIFIED_MD5) $(UNIFIED_MD5).build_date
```

### 5.2 Version pinning 과 audit

- `element_db/CHANGELOG.md` (신규): 빌드 날짜, md5, 엔트리 수, 이번 Round의 diff (추가된 P0 P1 15개).
- run_pipeline.py는 DB 경로 + md5를 `results/{sample}/run_manifest.json`에 기록하여 "어느 DB version으로 돌았는지"를 추적.

### 5.3 `element_db/README.md` 템플릿 (엔트리별 메타)

**`element_db/README.md`** 의 각 엔트리:

```markdown
### SpCas9 (plant codon-optimized)

| Field | Value |
|---|---|
| FASTA header | `>cds\|SpCas9\|Streptococcus_pyogenes\|MH782573.1\|plant_codon_optimized_CDS\|mvp_2026-04-16` |
| Source | GenBank MH782573.1 (Xiang et al. 2017 PLoS ONE) |
| Subregion | 1-4104 (CDS only, NLS 제외) |
| Length | 4,104 bp |
| Owner | gmo_expert (Round 3 MVP) |
| Rationale | CRISPR LMO 법적 판정, CLAUDE.md Known Pitfalls §Cas9-present |
| Cross-react | SpCas9 has no plant homolog — safe |
| Regulatory tag | K-LMO-2024-Cas9-required |
| Updated | 2026-04-16 |
```

엔트리 15개 × 각 섹션. regulatory 제출 시 이 메타가 evidence package.

---

## 6. Phase 2 로드맵 (MVP 이후)

Team-lead/won_yim 분류에 따라 Phase 2로 deferred된 항목을 명세:

### 6.1 KCGP (Korean Codex GMO Nomenclature) mapping

- **목적:** 식약처 보고서의 한글 명칭 ↔ EUginius 영문 tag ↔ genbank accession 3자 매핑.
- **산출물:** `element_db/kcgp_nomenclature.tsv` (국가 표준 명칭 → 코드 → RedGene FASTA header 매핑).
- **의존성:** 농촌진흥청 2024/2025 고시 PDF 크롤, OCR 필요. 수동 작업 2-3인일.

### 6.2 CRL 82 seq full integration

- MVP (§1)에서는 unique ~10개만 dedup 후 포함.
- Phase 2에서 나머지 70여개 amplicon도 `crl_amplicons`로 별도 tag 후 unified DB에 포함.
- 정책: amplicon source tag (`|crl_v1`) 보존 — BUG-3 재발 방지용 priority 비교에서 `element_db` tier로 처리 (§ Round 2 Q3 2-tier 정책).

### 6.3 Cas9-present flag override

- 현재 `config.yaml: expected.cas9_present`는 ground-truth 메타데이터. 
- Phase 2: SpCas9 BLAST hit이 confirmed되면 런타임 flag 자동 설정, ground-truth flag와 cross-validation → mismatch 시 warning.
- 구현: `scripts/s05_insert_assembly.py`에서 Cas9 hit 발견 시 `site["cas9_detected"] = True` 기록 후 `run_manifest.json`에 요약.

### 6.4 기타 deferred

- **FLP/Cre, T-HSP17.5, P-PcUbi:** 국내 승인 이벤트 거의 부재. 2026 하반기 재검토.
- **EPRV pre-mask:** 팀 합의로 기각. Filter D 유지.
- **Chain of custody (regulatory R-3):** won_yim Phase 2 분류. 별도 deliverable.

---

## 7. Round 4 검토 체크리스트

이 MVP 완결성 확인을 위해 Round 4에서 검증 요망:

- [ ] §1 P0 6개 efetch 스크립트 실제 실행 → FASTA 생성 가능한지
- [ ] §2 patch 적용 후 `common_payload.fa`가 11개 entry, 각 subregion 길이 검증
- [ ] §3 config 스키마가 기존 config.yaml parse를 깨지 않는지 (run_pipeline.py config load 테스트)
- [ ] §3 `apply_canonical_triplet` 로직이 rice_G281 ground truth에 대해 `[hLF1, P-Gt1, T-nos]` 조합을 실제 인식하는지 (G281_construct.fa + g281_elements.fa로 spot-check)
- [ ] §4 `build_host_endogenous_bed.sh`을 rice에 시범 실행 → BED 레코드 수가 합리적 (1-10개) 범위인지
- [ ] §4 **`build_host_endogenous_bed.sh soybean` 실행 + rationale 수동 큐레이션 완료** (MVP-BLOCKING — AtYUCCA6/UGT72E3 DB merge 선행 조건)
- [ ] §4 soybean BED에 GmYUC family 10+ 좌표, GmUGT family 15+ 좌표 포함 검증
- [ ] §1 AtYUCCA6/UGT72E3 DB merge는 §4 soybean BED curation 이후 순서 준수
- [ ] §5 Makefile로 dedup 실행 → ≤5% seq 감소 확인 (CD-HIT 예상 효과)
- [ ] §5 `run_manifest.json` 기록 구조 설계 확인
- [ ] bio_king W8과 연결: `--cas9-info-only` flag가 s05에서 CANDIDATE 승격 경로에 실제 반영되는지 (tomato_Cas9_A2_3 spot-check)

---

## 8. 참조

- Round 1: `docs/team-review/notes/gmo_expert_round1.md` — DB 충분성, per-host FP, soybean 0-CAND
- Round 2: `docs/team-review/notes/gmo_expert_round2.md` — CD-HIT 분석, pre-mask + MCscan 결합, 2-tier priority
- 관련 bug: BUG-3 (element_db > univec), BUG-4 (extra_dbs in annotate_insert), BUG-8 (AtYUCCA6 제거), BUG-9 (V00087.1 subregion 미분리)
- Regulatory: 식약처 고시 제2024-, 농촌진흥청 GMO 검사방법, EU CRL-GMOMETHODS
- Citation: Debode et al. 2013 (Food Analytical Methods), Christensen & Quail 1996 (Transgenic Res), McElroy et al. 1990 (Plant Cell), Klee et al. 1987 (MGG), Kim et al. 2023 (Mol Breeding) — soybean AtYUCCA6 event

---

**Status:** Round 3 MVP deliverable 완료. Round 4 team-lead 검토 대기.
