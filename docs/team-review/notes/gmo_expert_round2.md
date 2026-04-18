# GMO Domain Expert — Round 2 Cross-talk 답변

**작성자:** gmo_expert
**날짜:** 2026-04-16
**수신자:** team-lead (Q1: compute_eng 반박, Q2: haibao 반박, Q3: bio_king 제안)

---

## Q1 — DB 확장 비용 vs Regulatory MUST-HAVE subset

### 답변 요지
CD-HIT-est `-c 0.95` 클러스터링은 **실제로 큰 도움 안 됨** — Round 1에서 제안한 서열들의 상당수가 이미 다른 organism·다른 region 유래로 95% 이하 유사성. 유일하게 겹치는 대 후보는:
- `crl_amplicons.fa` 82개 ↔ `gmo_combined_db.fa` 131개: QL/QT prefix가 같고 accession도 공유 (AX033493, CRL-GMOMETHODS 등 중복) → 예상 redundancy 제거 후 **unique seq ~40개**. 그러나 이들은 amplicon 길이(60~500bp)가 다르고 primer 위치가 다른 partial sequence라 CD-HIT-est -c 0.95로도 독립 cluster 유지될 확률 높음 → 실질 감소분은 10~15개에 그칠 것.
- SpCas9 vs P-35S, PMI vs nptII 등 기능적으로 완전 독립 → redundancy 0.

**결론:** CD-HIT로는 목록 크기를 크게 줄이지 못하므로, **"regulatory MUST-HAVE 필터"로 자체 축소**가 답. 아래 15개 subset 제안.

### Minimal MUST-HAVE subset (15개 이하, 식약처·농진청 기준)

| # | Priority | Name | 분류 | 근거 (regulatory) | 예상 size |
|---|---|---|---|---|---|
| 1 | **P0** | SpCas9 CDS (plant codon-opt) | cds | CRISPR LMO 법적 판정의 유일한 molecular evidence | ~4.2 kb |
| 2 | **P0** | sgRNA scaffold (80bp tracrRNA-derived) | regulatory | CRISPR universal marker, host-free | 80 bp |
| 3 | **P0** | PMI (manA) | cds | Syngenta/Corteva 신형 positive selection; 격리소 샘플 점유율 증가 | ~1.2 kb |
| 4 | **P0** | vip3Aa | cds | MIR162 (수입 옥수수 주요 이벤트) | ~2.4 kb |
| 5 | **P0** | cry2Ab2 | cds | MON89034 (국내 수입 승인 옥수수) | ~1.9 kb |
| 6 | **P0** | cry34Ab1 + cry35Ab1 | cds ×2 | DAS-59122 (승인 이벤트, binary toxin) | ~0.45 + 1.1 kb |
| 7 | **P1** | T-pinII | terminator | Bt crops 보편, 현재 KP784699.1 내부 amplicon으로만 간접 존재 | ~250 bp |
| 8 | **P1** | T-g7 | terminator | Debode 2013 reference, 다수 이벤트 | ~250 bp |
| 9 | **P1** | CTP2 (AtEPSPS) | regulatory | GTS 40-3-2 (국내 최다 수입 soybean) | ~216 bp |
| 10 | **P1** | 2mEPSPS (maize) | cds | GA21 (승인 옥수수) | ~1.3 kb |
| 11 | **P1** | AtYUCCA6 | cds | soybean 0-CAND 해결 직결 | ~1.1 kb |
| 12 | **P1** | UGT72E3 | cds | soybean 0-CAND 해결 직결 | ~1.4 kb |
| 13 | **P1** | P-nos subregion (V00087.1:1847-2113) | promoter | BUG-9 잔존 silent 오분류 해결 | 267 bp |
| 14 | **P1** | T-nos subregion (V00087.1:1277-1536) | terminator | BUG-9 동일 | 260 bp |
| 15 | **P2** | CRL subset (unique 10~12개) | amplicons | 식약처 참조 CRL 호환성 최소 | ~1 kb 총합 |

**추정 비용 영향:** 위 15개 총 bp ≈ 17~20 kb. 현 element_db 121 kb에 추가해도 **DB 크기 +17%**. BLAST 비용은 hit candidate 수 증가가 주 요인이지 DB size linear에 그치지 않음. 단, **SpCas9 CDS (4.2 kb)는 tomato 샘플의 Phase 1.5 positive site 수를 증가**시킬 가능성 — 이는 compute_eng의 우려가 맞음. 대응책: Phase 1.5 Cas9 hit은 informational tag로만 처리하고, CANDIDATE 승격 트리거에서 제외하는 flag (`--cas9-info-only`).

**P3로 제외한 이유:** FLP/Cre (marker excision, 국내 승인 이벤트에 거의 부재), T-HSP17.5 (신형이나 국내 샘플 아직 없음), P-PcUbi (BASF 전용), GOX (시대적 표본 부족). 이들은 Phase 2에서 추가 검토 권장하지 Round 1 P0/P1에서 제외.

---

## Q2 — MCscan synteny + pre-mask 결합 가능성 / regulatory 안전성

### 답변 요지
**결합 권장. 단, "pre-mask-with-audit-trail" 방식으로.** 양자택일이 아님.

### 구조
```
1. MCscan (host vs outgroup related genome)으로 synteny block 후보 자동 도출
   → EPRV, Ubi/Act 등 conserved region의 BED 초안
2. 도메인 전문가 수동 큐레이션 (식약처 제출용 evidence 문서화)
   → host_masked.bed + host_masked_rationale.tsv (각 구간별 "왜 mask했는지" 법적 근거)
3. s02 BWA: 원본 host.bam 생성 후, s03b 또는 신설 s02b에서 BED-intersect로 "mask된 영역에 걸린 soft-clip junction"을 **FALSE_NEGATIVE_MASKED** 태그로 분리
4. 최종 report에 mask된 영역 인접 junction은 FALSE_NEGATIVE_MASKED로 명시 — 은폐 X, 명시적 비판정
```

### Regulatory 안전성 판정
식약처·농진청 acceptance criteria는 "false negative는 명시되어야 한다"가 핵심. **pre-mask로 reads를 버리는 것 자체가 위반은 아니다** — 문제는 "왜 버렸는지"의 audit trail. 위 3-단계가 "BED rationale.tsv + FALSE_NEGATIVE_MASKED 태그"를 제공하면 규제 요건 충족.

**다만 한 가지 중요 단서:** Pre-mask BED에는 **~100% ortholog 영역 (corn × P-Ubi1, rice × P-Act1)만 포함**시킬 것. EPRV처럼 "host에 integrated pararetrovirus 조각"은 종종 실제 T-DNA insertion과 인접하거나 host 내 hotspot이 될 수 있음. 이를 pre-mask하면 일부 진짜 insertion까지 놓칠 위험 — EPRV는 pre-mask가 아니라 사후 Filter D (alt-locus) 쪽으로 처리.

### 요약 판정 matrix
| 대상 | Pre-mask (s02 전) | 사후 filter (Filter D) | Rationale |
|---|---|---|---|
| 100% ortholog (corn × P-Ubi1, rice × P-Act1) | ✅ | — | 100% 동일은 read mapping 자체가 의미 없음 |
| Soybean AtYUCCA6 host ortholog | ✅ | — | 75% 유사도지만 systemic FP source |
| EPRV / host pararetrovirus 잔재 | — | ✅ | insertion hotspot 가능성, pre-mask 위험 |
| 다른 60~80% 유사 영역 | — | ✅ | 기존 s03b + Filter D로 충분 |

haibao의 MCscan 제안은 **pre-mask BED 후보 자동 생성**에 완벽히 적합. 내 Round 1 주장과 충돌 없이 결합. 두 방법 결합이 최적.

---

## Q3 — 3-tier priority vs "element_db-family > univec"

### 답변 요지
**3-tier (common_payload > element_db > CRL > univec)가 domain 관점에서 정당하지만, 구현상 2-tier (canonical > univec)로 충분**. 이유: 계층이 깊어질수록 BUG-3 class 재발 위험이 다시 올라감 (tie-break 정책이 복잡해짐).

### Domain justification
| Tier | 구성 | Sequence 특성 | 신뢰도 |
|---|---|---|---|
| 1 | common_payload | canonical CDS, GenBank fetched, single-organism, single-accession | **최고** (well-characterized marker) |
| 2 | element_db (EUginius) | amplicon ~80-500bp, qPCR primer pair 기반, regulatory-validated | 높음 |
| 3 | CRL amplicons | EU 표준 amplicon, 일부 CRL-GMOMETHODS 자체 curation (no accession) | 중간 |
| 4 | univec | NCBI UniVec, vector backbone 전체, redundant, old | 낮음 |

Domain상 common_payload가 element_db보다 우선하는 것은 타당함 — common_payload는 "이 marker가 있다"는 **definitive 증거**, element_db amplicon은 "이 element family를 target으로 설계된 primer 증폭 영역" (diagnostic 증거). 즉 **CDS full hit > amplicon hit > vector backbone hit** 순이 domain적으로 옳음.

### 그러나 구현 권장은 2-tier
BUG-3은 "동일 tier 내부에서 tie를 어떻게 깰 것인가"가 아니라 "element_db 계열을 univec보다 명시적으로 우선시킬 것인가"가 본질. 따라서:

```python
SOURCE_PRIORITY = {
    "common_payload": 2,
    "element_db":     2,   # 동일 tier로 합침
    "crl_amplicons":  2,
    "univec":         1,
}
# Tie (같은 tier 내부) → bitscore 비교
# Tier 차이 있으면 → 상위 tier 승
```

이렇게 하면:
- `_should_replace(existing, new)` 로직이 단순 (priority 비교 1회 + bitscore 비교 1회)
- 3개 element DB-family 서로 간에는 bitscore가 tie-break → 데이터 특성 따라 자연스러운 선택
- BUG-3 재발 방지 핵심인 "element_db-family 서열이 있는 한 univec hit으로 덮이지 않음"이 보장

### 추가 방어선
CDS-level full-length hit은 amplicon hit보다 거의 항상 bitscore가 높음 (서열 길이가 10배 이상). 따라서 domain상 우선순위가 동일 tier 내부 bitscore 비교로 자연스럽게 드러남 — 별도 계층 분리 없어도 OK.

**bio_king 제안 수용 + 경계 조건:** "element_db-family 모두 univec보다 우선"으로 단순화. Tie-break은 bitscore. 단, `common_payload`의 서열은 canonical accession full CDS이므로 bitscore 우위가 자연스럽게 보장 → 명시적 tier 2-tier로 충분.

---

## 요약 matrix (3 답변 한 줄 요약)

| Q | 결론 | 핵심 근거 |
|---|---|---|
| Q1 | CD-HIT로 절감 불가, 대신 **regulatory MUST-HAVE 15개 subset** 제안 (P0 6개 + P1 8개 + P2 1개) | 기능적 독립 서열이라 redundancy 낮음; 격리소 법적 요구에 맞춰 자체 필터링이 답 |
| Q2 | **MCscan + pre-mask 결합 권장** (100% ortholog만 pre-mask, EPRV는 사후 Filter D) + FALSE_NEGATIVE_MASKED 태그 + rationale.tsv | 양자택일 아님. audit trail로 regulatory 요건 충족 |
| Q3 | **2-tier 충분** (element_db-family > univec), tie는 bitscore. 3-tier 구현은 BUG-3 class 재발 위험 | Domain 우선순위는 자연스러운 bitscore 차이로 표현됨 |

---

**상태:** Round 2 답변 완료. Round 3 또는 종합 투표 대기.
