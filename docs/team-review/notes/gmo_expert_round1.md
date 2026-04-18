# GMO Domain Expert — Round 1 Notes

**작성자:** gmo_expert (team: redgene-review-2026-04-16)
**날짜:** 2026-04-16
**범위:** Element/transgene DB 정합성, host-endogenous FP, Korean quarantine 요건, soybean 0-CANDIDATE의 domain 해석

---

## 1. 현 element DB 충분성 평가

### 요약 (한 줄)
| DB 파일 | 서열수 | 평가 | 근거 |
|---|---|---|---|
| `element_db/gmo_combined_db.fa` | 131 | **부분적으로 충분** | EUginius canonical 14 CDS + 10 promoter + 4 terminator는 20세기 말~2010년대 이벤트는 커버. 그러나 **SpCas9 CDS, sgRNA scaffold, PMI, ipt, vip3Aa, cry2Ab2, cry34/35Ab1, T-g7, T-pinII, T-HSP17.5, T-rbcSE9(독립서열 X), P-PcUbi, P-MTL, CTP2, 2mEPSPS** 등이 canonical FASTA 엔트리로 빠져 있음 (`element_summary.tsv`에는 "curated / full / 0" 로만 명시되어 실제 서열이 없음) |
| `element_db/common_payload.fa` | 9 | **충분** (한계 명확) | bar/nptII/hpt/gusA/gfp/egfp/P-CaMV35S/P-nos/T-ocs. BUG-9 (V00087.1 중복)는 `b61b7f1`로 해결. P-nos와 T-nos를 같은 accession에서 subregion으로 분리하는 efetch -seq_start/-seq_stop 추가는 여전히 TODO. |
| `element_db/gmo_corn_combined_db.fa` → `db/gmo_corn_combined_db.fa` | 192 | **옥수수 한정 충분** | EUginius 130 + corn LB/RB border 62 (31 events × LB+RB). CRL amplicon 82 seqs는 빠짐. 옥수수 전용이므로 GM 옥수수 이벤트 ID는 잘 커버. |
| `element_db/transgene_db.fa` | 수천? | **UniVec 그대로** | 1.76 MB, Phase 1.5의 primary DB (파이프라인이 `univec`로 참조). Vector backbone 스크리닝에는 적합하지만 T-DNA element 식별 해상도 낮음 — 그래서 element_db가 univec보다 우선되도록 BUG-3에서 `_should_replace` 추가됨 (커밋 `4d3f2bb`). |
| `element_db/crl_amplicons.fa` | 82 | **미통합** | CRL-GMOMETHODS amplicon 82개. `element.md`에는 언급되어 있으나 파이프라인이 현재 공급하는 `extra-element-db`에는 빠져 있음. 식약처/농진청 인증 assay 호환성을 위해 포함 권장. |
| `element_db/g281_elements.fa` | 4 | **샘플 전용** | hLF1(NM_002343.4), G6-EPSPS(AY945950.1), Gt1(X54314.1), Pepc(X15239.1). 수동 작성으로 **rice G281 이벤트에만 적용**되는 custom payload. 다른 rice line 검증에는 역으로 부담이 될 수 있음. |

**핵심 결론:** 기존 131 EUginius는 **한국 격리소 2010년대 이전 이벤트 스크리닝에는 OK**이지만, 2016년 이후 **CRISPR-based gene editing 기반 LMO (Cas9 vector), 신형 marker (PMI/ipt), 신형 terminator (T-HSP17.5, T-rbcSE9)** 검출에는 공백이 존재. `construct_reference`를 `element_db/gmo_combined_db.fa` 그 자체로 쓰는 tomato/cucumber/soybean 샘플들은 이 공백이 직접 Phase 2(host mapping) 입력에 영향을 준다.

---

## 2. 누락 또는 확장 권장 element (우선순위 표)

| 이름 | 분류 | 왜 필요 | 우선순위 | 제안 accession |
|---|---|---|---|---|
| **SpCas9 CDS** | cds | CRISPR LMO 탐지에 필수. `tomato_Cas9_A2_*` 샘플의 construct_reference는 `gmo_combined_db.fa`지만, **DB에 SpCas9가 없음**. Cas9 존재 판정은 현재 `config.yaml: expected.cas9_present` flag에만 의존. | **P0** | NC_017053.1 (*S. pyogenes* SF370 Cas9) 또는 addgene pX330 주형 (plant codon-optimized) |
| **sgRNA scaffold** (80bp tracrRNA-derived) | regulatory | CRISPR 이벤트의 universal marker (sequence-defined, ~80bp). `tomato_grna.txt` 존재하나 DB에 scaffold 서열 없음. | **P0** | pMDC123, pHEE401E, pRGEB32 등 (식물 CRISPR 대표 backbone) |
| **Cas9 NLS (nucleoplasmin NLS)** | regulatory | Plant Cas9 vector의 ~80% 는 C-terminal npNLS 사용. Cas9 존재 간접 증거. | P1 | pX330: AAAAAGAAGAAGAGGAAGGTTGATCCA |
| **PMI (phosphomannose isomerase)** | cds | 신세대 positive-selection marker (대체로 Syngenta MIR604). `element_summary.tsv`에 없음. | P1 | Z11946.1 (*E. coli* manA) |
| **ipt (isopentenyl transferase)** | cds | Cytokinin self-excision marker. MAT vector 시스템. | P2 | X00639.1 |
| **FLP recombinase / FRT sites** | regulatory | Marker excision 시스템, GLA21 등. | P2 | X03772.1 (FLP), 48bp FRT consensus |
| **Cre / loxP** | regulatory | 유사 marker excision. | P2 | X03453.1 |
| **vip3Aa** | cds | MIR162 maize (신형 Bt). `element_summary.tsv`는 "curated/full/0" = 실제 서열 없음. | P1 | KX344421.1 또는 DOE database |
| **cry2Ab2** | cds | MON89034. 동일하게 서열 누락. | P1 | AY007234.1 |
| **cry34Ab1 / cry35Ab1** | cds | DAS-59122 binary toxin. 서열 누락. | P1 | AY126452.1 / AY126453.1 |
| **T-g7** (Agrobacterium gene 7 terminator) | terminator | Debode 2013 등에서 자주 등장. 서열 누락. | P1 | AF242881.1:1-250 |
| **T-pinII** (potato proteinase inhibitor II) | terminator | Bt 이벤트에 자주 사용. KP784699.1 within QL-CON-gat_T-pinII에만 간접 포함. 단독 canonical 엔트리 없음. | P1 | D00063.1 또는 X04118.1 |
| **T-HSP17.5** (heat shock terminator) | terminator | Newer Corteva/Syngenta vectors. | P2 | Z11575.1 |
| **T-rbcSE9-pea** (full seq) | terminator | T-E9는 존재(X00806.1)하나 pCAMBIA 변종·Dow 변종과의 분화 필요. | P2 | PvRbcS 변형들 |
| **GOX (glyphosate oxidoreductase)** | cds | GT73 (canola). `element_summary.tsv`에 curated 표기만. | P2 | L48604.1 |
| **2mEPSPS (maize)** | cds | GA21. "host-derived modified" — host 오탐 매우 높음. 서열 없음. | P2 | X63374.1 |
| **CTP2 (Arabidopsis EPSPS)** | regulatory | RR soybean GTS 40-3-2. `element.md`에 "accession unclear" 표기됨. | P1 | AT1G48860 (TAIR11) 또는 X64538.1 subregion |
| **P-PcUbi (parsley ubi)** | promoter | BASF 이벤트 (LibertyLink). 서열 없음. | P2 | ? (BASF 특허도면) |
| **P-MTL / P-PLD3** | promoter | 옥수수 haploid induction (DH-lines) — Korean quarantine과는 무관하나 확인 필요. | P3 | |

### 서열 부재 위험도가 특히 높은 Top-3
1. **SpCas9 + sgRNA scaffold**: tomato A2_* 세 샘플이 Cas9-present 이벤트인데, 현재 s05 `annotate_insert`로 "Cas9"가 직접 검출되지 않음. Cas9의 존재·유무 판정이 config flag에 의존 = 임상·규제 proof 용도로는 불충분.
2. **PMI/ipt/FLP**: 최근 국내 승인 LMO (Syngenta Golden Rice 2 포함) 에서 사용되는 marker이므로 격리소 업무 보장에 영향.
3. **vip3Aa/cry2Ab2/cry34-35Ab1**: Bt 계통 신형 toxin. 수입 옥수수·면화 검출에 필수.

---

## 3. Per-host endogenous FP blacklist

아래는 construct element 서열이 host genome의 endogenous ortholog와 sequence similarity를 공유하여 MAPQ=60 false positive를 만드는 실제 교차반응 매트릭스. `CLAUDE.md` Known Pitfalls와 `element_summary.tsv` 기반.

| Host 유전체 | Element DB 서열 | Endogenous homolog / 위치 | 교차반응 유형 | 현재 대응 |
|---|---|---|---|---|
| **Rice** (Osativa_323_v7.0) | P-Ubi1-maize (S94464.1) | OsUbi1/2/3 (LOC_Os03g13170, LOC_Os02g06640) | Ubiquitin 5'UTR+intron1 ~70% 유사 | s03b (WT homology filter) |
| Rice | P-Act1-rice (S44221.1) | Act1 자기자신 (LOC_Os11g06390) | **100% 동일** — P-Act1이 rice construct에 쓰이면 host reference에서 그대로 매핑 | s03b만으로 불충분, WT 필수 |
| Rice | QT-TAX-OS-002 (XM_052303253.1) | taxon-specific 서열 자체가 rice ID용 | 설계상 rice에서 양성 | Phase 1.5에서 taxon-specific는 construct hit로 처리 금지 필요 |
| Rice | QT-TAX-SPS (KT225496.1) | sucrose phosphate synthase | rice 자체 유전자 | Phase 1.5 제외 필요 |
| **Tomato** (SLM_r2.0) | P-TA29 (X52283.1) | Nicotiana homolog의 Solanum ortholog (SlTA29) | Solanaceae 공통 tapetum 유전자 | s03b |
| Tomato | P-SSuAra (X13611.1) | SlRbcS (multi-copy) | Calvin cycle 유전자 ~65% 유사 | s03b 필요 |
| Tomato | Endogenous pararetrovirus (EPRV) | Solanum pararetrovirus 조각이 genome에 integrated | CaMV/FMV 계열과 LTR 유사성 | **현재 미대응** — EPRV mask 필요 |
| **Cucumber** (CucSat_B10v3) | P-Ubi1-maize, P-Act1-rice | CsUBQ10, CsACT2/7 | monocot-based 서열이지만 actin은 dicot과 보존도 높음 | s03b 필요 |
| **Corn** (Zm_B73_v5) | **P-Ubi1-maize (S94464.1)** | ZmUBI1 자기자신 (100%) | **100% 동일** — 옥수수 host일 때 자기자신 매칭 | `--wt_control` 필수, s03b로도 bypass 안 됨 |
| Corn | QT-TAX-ZM-003/004, zSSIIb, wx012 | 해당 유전자 자체 | 옥수수 ID 용 탐지 타깃 | Phase 1.5 제외 필요 |
| Corn | P-Act1-rice | ZmActin family | 단자엽 actin 보존 | s03b |
| Corn | I-hsp70-maize (custom 추가분) | 옥수수 자기 자신 (X03714.1) | **100% 동일** | corn host일 때 construct hit로 처리하면 오탐 |
| **Soybean** (Gmax_v4.0) | P-CaMV35S, P-FMV | EPRV (soybean mosaic virus family integrated copies) | pararetrovirus 잔재 | EPRV mask 필요 |
| Soybean | P-Ubi1-maize | GmUBI family | ubiquitin 보존도 | s03b |
| Soybean | cp4-epsps (L29358.1) | GmEPSPS (Glyma01G179200 등) | CP4 ≠ Gm EPSPS 이지만 보존 motif ~45% | 일반적으로 문제 아님 |
| Soybean | AtYUCCA6 (NM_122473.3) — **Arabidopsis 식물 유전자** | GmYUCCA family 10+ copies | **YUC6 ↔ Glyma06g12890/19g11570 등 ~75% 동일** | **BUG-8 기 식별** — 식물 유전자를 marker로 쓰면 host cross-react 최악. `0c5a9c1`에서 제거됨 |

### 특히 위험한 페어 3개
1. **Corn host + P-Ubi1-maize (100%)**: corn ND207 샘플은 construct 내부에 P-Ubi1이 있으면 host에서도 무조건 match. 현재 ND207은 `db/gmo_corn_combined_db.fa`를 사용하나 여기에도 P-Ubi1-maize 포함.
2. **Rice host + P-Act1-rice (100%)**: rice line이 P-Act1을 사용하는 경우 host-fraction 필터만으로는 불충분. G281 construct는 P-Pepc/P-Gt1 기반이라 현재 샘플은 안전하지만, 향후 rice 이벤트 확장 시 위험.
3. **Soybean + AtYUCCA6 (식물 유전자)**: `0c5a9c1`로 제거되었으나 **이 제거가 오히려 현재 soybean 샘플의 0-CANDIDATE 원인 가능성 높음** (아래 §5 참조).

---

## 4. common_payload.fa 9개 marker 평가

| marker | accession | 실질 영향 | 판정 |
|---|---|---|---|
| bar | X17220.1 | 356bp CDS. rice/tomato/cucumber/soybean 이벤트 다수에 포함. 잘 선택됨. | OK |
| nptII | V00618.1 | 795bp CDS. kanamycin selection marker canonical. | OK |
| hpt | M55269.1 | 1026bp CDS. Hygromycin marker. | OK |
| gusA | X06788.1 | uidA reporter. | OK |
| gfp | U55762.1 | A. victoria 원본. | OK (한계: 식물 codon-optimized 변종은 sequence drift 있음) |
| egfp | L29345.1 | Enhanced GFP. | OK |
| P-CaMV35S | V00141.1 | CaMV 완전 게놈 → 35S region 전체 (full 8kb accession). Efetch할 때 subregion 자르지 않았다면 P-35S + T-35S 둘 다 포함 → **이중 검출 혼동 위험**. | **확인 필요** |
| P-nos | V00087.1 | Agrobacterium Ti octopine plasmid nos operon. P-nos와 T-nos가 같은 accession 안에 존재. BUG-9에서 T-nos 중복 발견. **현재 P-nos가 실제로 promoter region만 담고 있는지, 아니면 nos operon 전체인지 verify 필요**. | **확인 필요** |
| T-ocs | X04879.1 | A. tumefaciens octopine synthase. ocs operon. 마찬가지로 subregion 검증 필요. | 확인 필요 |

**BUG-9 실제 영향 재평가:**
- Fix `b61b7f1`로 T-nos 중복 row는 제거됨.
- 그러나 P-nos entry 자체가 V00087.1 full sequence이면 T-nos region까지 포함되어, BLAST hit이 P-nos 태그로만 보고되는 silent 오분류가 여전히 발생 가능.
- **권장 조치:** `build_common_payload.sh`에 `efetch -seq_start 1847 -seq_stop 2113` (P-nos) / `-seq_start 1277 -seq_stop 1536` (T-nos) 추가하여 subregion 전용 FASTA 생성. `element_summary.tsv` 좌표와 일치.

---

## 5. Soybean 0-CANDIDATE 근본 원인 hypothesis (domain 관점)

resume.md §"Soybean 0-CANDIDATE problem" 에 4가지 가설이 있었음. Domain 관점으로 재분석:

### H1: AtYUCCA6는 식물 유래 → 본질적 host cross-react
- **근거 (BUG-8):** `NM_122942.3` 최초 typo → `NM_122473.3` (YUC6 / AT5G25620) 정정 → 결국 `0c5a9c1`에서 전면 제거.
- **결과:** AtYUCCA6의 CDS가 element_db에서 사라지면, s04b가 construct reads에서 AtYUCCA6 서열을 포함한 contig를 만들어도, s05 `annotate_insert`가 그것을 element로 인식 못 함 → **s05 Phase 3에서 NODE_X가 element hit을 얻지 못하고 UNKNOWN으로 남음**. 이것이 30 UNKNOWN/0 CANDIDATE의 직접 원인일 가능성 매우 높음.
- **해결안:** AtYUCCA6 CDS (NM_122473.3, 약 1.1kb)를 `gmo_combined_db.fa`에 다시 추가하되, **soybean host에만 한정된 filter**를 파이프라인 단에 둠. 혹은 AtYUCCA6 대신 **flanking non-plant 서열**(bar marker, P-CaMV35S, T-nos 등 canonical 요소)만으로 CANDIDATE 판정하는 mode를 추가.
  - 실제 Kim et al. 2023 (Mol Breeding) paper의 construct는 `35S::AtYUCCA6:T-ocs` + `35S::bar:T-nos` — bar와 P-35S는 이미 common_payload에 있음.
  - 즉 **"bar marker가 insert 내부에 존재하면 CANDIDATE로 promote"** 하는 rule이 있다면 0 → N+ 가능.

### H2: s04b contig filter가 너무 aggressive
- resume.md: soybean 1,345 → 6 filtered (99.6%). 이 6개 contig가 insert 영역을 전부 커버한다고 장담 못함.
- 특히 T-DNA가 단독 copy + 가까운 bacterial backbone 혼합 integration인 경우, SPAdes contig가 구간별로 쪼개지면 filter가 "marker 포함" 기준을 통과하는 건 한두 개뿐일 수 있음.
- **검증 방법:** `results/soybean_AtYUCCA6/s04b_construct_assembly/contigs.fasta` 6개와 contigs_all.fasta 1,345개를 각각 s05에 집어넣고 Phase 1.5 positive site수 / 최종 CANDIDATE 차이 비교.

### H3: Gmax v4.0 genome 내부 retroviral elements + CaMV35S-like region
- Soybean은 SMV/BPMV와의 오랜 공진화로 endogenous pararetrovirus-like sequences가 존재. `P-CaMV35S` BLAST hit이 host region에서 발생하면 soft-clip 기반 junction detection 시 양쪽 모두 MAPQ=60 가능.
- Filter D (alternative-locus, minimap2)가 이를 잡는다면 FP 쪽으로 떨어짐. 이것이 15 FP의 상당 부분일 가능성.

### 결론 및 권장
1. **BUG-8 제거 정책 재고.** AtYUCCA6/UGT72E3처럼 식물 유래 CDS를 element_db에 넣는 것은 host cross-react 위험이 크지만, construct-specific은 정확히 그 식물 유전자에 기반한 전사체. 해결책은 제거가 아니라 **per-sample white-list**: `config.yaml`에 `payload_cds: [AtYUCCA6, UGT72E3]`을 두고, s05가 host 매칭 시 이 flag를 infer해서 CANDIDATE 판정.
2. **CANDIDATE 판정 rule 확장.** 현재 "element_db hit이 있어야 CANDIDATE"라는 정책으로 보임 (`_should_replace`의 element_db > univec). 이 정책이 construct-specific payload CDS가 DB에 없는 soybean 같은 샘플에 불리. bar/P-35S/T-ocs 같은 canonical marker만 있어도 CANDIDATE로 승격하는 rule 추가 권장.
3. **식약처/농진청 관점**: 규제 목적으로는 "bar + P-35S-CaMV + T-ocs" 같은 canonical triplet이 single contig에 있으면 LMO-positive로 확정하기에 domain적으로 충분함. AtYUCCA6 존재 여부는 construct ID 정도.

---

## 6. 한국 규제 요건 체크리스트 (식약처/농진청/환경부)

### 6.1 법적 근거
- **LMO법** (유전자변형생물체의 국가 간 이동 등에 관한 법률, 2008~) + 시행령
- **식약처 고시 제2024-** 유전자재조합식품 안전성 심사: qPCR event-specific assay 요구
- **농촌진흥청 고시 GMO 검사방법**: CRL-GMOMETHODS 기반 element screening
- **환경부**: 방출 LMO 환경영향평가 (비의도 유출 검정)

### 6.2 detection assay가 만족해야 할 최소 요건
| 요건 | 현 RedGene 상태 | 개선 필요 |
|---|---|---|
| P-CaMV35S / T-nos 2-plex screening | 있음 (gmo_combined_db + common_payload) | OK |
| Event-specific (LB/RB junction) 확인 | 부분 (corn_border_db only) | tomato/cucumber/soybean event junction 서열 추가 |
| Copy number (단일 vs 다중 T-DNA) | s07_copynumber 있음 | bp-level quantification 검증 필요 |
| Zygosity (호모 vs 헤테로) | s07 내부 구현 | 공식 validation 필요 |
| Cas9/sgRNA 존재 확인 | **없음** (config flag로만) | SpCas9 + scaffold 서열 추가 필수 |
| CRISPR editing (indel) 검출 | s06_indel | OK |
| LOD (limit of detection) — 5x coverage | CLAUDE.md: "≥10x 권장, 5x 부분" | 5x 실험 검증 필요 (이미 subsampled 샘플 있음) |
| chain-of-custody / audit trail | (없음) | run ID, DB version, md5 기록 파이프라인화 |
| CRL reference amplicon 일치도 | crl_amplicons.fa 82개 **미통합** | DB 파이프라인에 합류 |

### 6.3 특히 미흡한 세 가지
1. **Cas9 검출 기능 부재**: tomato_Cas9_A2_* 샘플이 "Cas9 presence"를 식별한다고 주장하지만 서열 기반 판정 아님.
2. **Event-junction 서열 부족**: 옥수수는 62 corn-border + 29 corn-event이지만 벼·토마토·오이·콩은 event-specific 서열 없음.
3. **CRL amplicon 미통합**: 식약처가 실무적으로 참조하는 EU CRL 표준 amplicon DB가 파이프라인에 안 들어가고 있음.

---

## 7. 다른 teammate에게 묻고 싶은 것

### bio_king (s05 refactor 담당)
1. `annotate_insert` 함수가 `extra_dbs` 루프에서 BLAST 결과를 merge할 때, **primary hit vs secondary hit의 tagging 규칙**을 어디서 정하고 있나? BUG-3의 `_should_replace`와 유사하게 "element_db hit이 common_payload hit보다 우선" 또는 그 반대의 정책이 필요. soybean 샘플 재시도 시 AtYUCCA6를 element_db에 다시 넣으면 common_payload와의 충돌 가능성.
2. "bar + P-35S + T-ocs triplet 존재 시 CANDIDATE 승격" 같은 **canonical-marker override** rule을 code에 넣을지, config의 `payload_whitelist`로 뺄지 의견?
3. Cas9 CDS를 `element_db`에 추가하면 `tomato_Cas9_A2_*` 샘플에 대한 Phase 1.5 site count가 폭증할 가능성 (host에 유사 endonuclease 없으므로 FP는 적지만 construct 내부 중복 hit). 적용 시 배치 재검증 필요.

### compute_eng (HPC/SLURM 최적화 담당)
1. DB를 10개 → 20개로 확장 시 BLAST 1회 비용이 거의 선형 증가. 현재 UGT72E3가 24h TIMEOUT인데 DB를 2배로 하면 48h로 늘어남. 해결 방안은 (a) 서열 길이가 짧은 amplicon은 별도 low-evalue cut-off로 빠른 스크리닝, (b) DB 클러스터링 (CD-HIT 95%)으로 redundancy 제거.
2. `common_payload.fa`를 한 번만 makeblastdb하면 되는데, s05가 매 샘플마다 재build하는지 확인 필요 (`.nhr/.nin/.nsq` cache 공유).
3. crl_amplicons.fa (82 seq, 8kb)를 추가했을 때 Phase 1.5 속도 측정 제안.

### haibao (알고리즘 재평가 담당)
1. Sample-specific DB (s04b contig filter) 접근 vs **full-catalog 단일 DB + synteny-based reclassification** 접근 중 어느 쪽을 선호? 예: 전체 element + corn-border + CRL 모두 포함한 "universal DB" + 각 host의 endogenous ortholog에 대한 contra-DB를 BLAST reciprocal로 돌려서 FP를 사후에 빼는 방식.
2. **Synteny-based host endogenous masking**: rice P-Act1, corn P-Ubi1 같은 100% 동일성 문제는 minimap2 alt-locus 필터(Filter D)로 대부분 잡을 수 있다고 보는데, 실제 validation 결과가 있는지?
3. Soybean EPRV 같은 pararetrovirus-in-host는 host reference에서 pre-mask 하는 것이 답인가, 아니면 s05의 FP filter에서 처리하는 것이 답인가?

### won_yim (production 기준/검증 담당)
1. 한국 격리소 production 환경에서 "CANDIDATE 판정 민감도 vs 특이도"의 법적 허용치는? 식약처는 false-negative를 특히 싫어함 (수입 LMO 누락).
2. Cas9 검출 기능이 지금 없는데 tomato A2_* 샘플 검증 리포트에 "Cas9-present" 라고 표기해도 되는 legal risk 범위는?
3. 각 host별 subsampled coverage (3x/5x/10x/15x) validation이 실제 격리소 input data 품질과 일치하는가?

---

## 8. Round 2 갈등 예상

- **vs compute_eng:** DB에 20개 element 추가 시 BLAST 비용 +50% — compute_eng는 반대할 가능성. 대안: CD-HIT 클러스터링으로 redundancy 제거 후 추가.
- **vs bio_king:** "canonical-marker override" 를 코드에 넣을지 config로 뺄지. 나는 **config로 빼야** 한다고 주장할 예정 (host별/이벤트별 white-list가 regulatory flexibility 확보).
- **vs haibao:** "host pre-mask vs 사후 filter". 나는 **host pre-mask(EPRV, 100%-ortholog 구간)를 미리 BED로 빼고 s02 BWA 단계에서 제외**하자고 주장할 예정 (SLURM 시간 절감 + 알고리즘 단순화).

---

## 부록 A — 현재 DB sequence 카운트 요약

```
element_db/gmo_combined_db.fa     : 131 (EUginius + 1 thaumatin II)
element_db/common_payload.fa      :   9 (bar, nptII, hpt, gusA, gfp, egfp, P-CaMV35S, P-nos, T-ocs)
element_db/gm_new_elements.fa     :   3 (CTP4, I-hsp70-maize, P-FMV34S-full)
element_db/g281_elements.fa       :   4 (hLF1, G6-EPSPS, Gt1, Pepc) — rice G281 전용
element_db/crl_amplicons.fa       :  82 (CRL-GMOMETHODS amplicons) — 미통합
element_db/euginius_fullseq.fa    : ~98 (EUginius full sequences)
element_db/transgene_db.fa        : ≈7k (UniVec full)
db/gmo_corn_combined_db.fa        : 192 (EUginius 130 + corn LB/RB 62)
db/Cas9_construct.fa              : 131 (실체는 EUginius 복제, Cas9 서열 없음)
db/G281_construct.fa              : 132 (EUginius 131 + custom rice G281 payload)
```

## 부록 B — 현 DB에 반드시 추가해야 할 5개 서열 (P0~P1)

| Priority | Name | Accession | Rationale |
|---|---|---|---|
| P0 | SpCas9 CDS (plant codon-opt) | (pX330 또는 pRGEB32) | tomato A2_* sample Cas9 검출 |
| P0 | sgRNA scaffold (80bp) | (보편 서열) | CRISPR universal marker |
| P1 | AtYUCCA6 CDS | NM_122473.3 | soybean sample 0-CANDIDATE 해결 |
| P1 | UGT72E3 CDS | (Soybean-opt) | soybean UGT72E3 sample 0-CANDIDATE 해결 |
| P1 | crl_amplicons.fa 병합 | (82 seq) | CRL 표준 호환성 |

---

**상태:** Round 1 분석 완료. Round 2 Discussion에서 bio_king/compute_eng/haibao와 DB 확장 및 per-host FP 전략 정합.
