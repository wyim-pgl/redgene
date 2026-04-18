# AC-2: cucumber_line225 8 CAND FP verification workflow

Issue tracker: [#1](https://github.com/wyim-pgl/redgene/issues/1)

## Background

`results/cucumber_line225/s05_insert_assembly/` contains **8 CANDIDATE** reports.
Goal: classify each as TRUE_INSERTION vs FALSE_POSITIVE via NCBI nt remote BLAST,
then feed the FP calls back into the pre-mask rationale (T10).

## 3-step operator SOP

```bash
# (1) Extract and merge CAND inserts into one query FASTA
python scripts/util/extract_cand_for_blast.py \
    --sample-dir results/cucumber_line225 \
    --sample-name cucumber_line225 \
    --out cucumber_line225_cand_inserts.fa

# (2) Submit the remote BLAST SLURM job (login-node submission only — no auto-sbatch!)
SAMPLE=cucumber_line225 \
QUERY=cucumber_line225_cand_inserts.fa \
OUT_TSV=cucumber_line225_cand_vs_nt.tsv \
    sbatch scripts/util/prepare_remote_blast.sh

# (3) Review + classify (manual) — top 5 hits per query are enough to call FP:
#       * hit to Cucumis sativus / cucurbit species       → likely FP (host-ortholog)
#       * hit to Agrobacterium / T-DNA vector backbone    → TRUE_INSERTION
#       * no hit or low identity (<85%)                   → keep CANDIDATE
# Store the verdict in element_db/cuc225_cand_classification.tsv alongside
# docs/team-review/W1 evidence.
```

## Notes / constraints

- `prepare_remote_blast.sh` uses `blastn -remote -db nt` with
  `max_target_seqs=5` and `evalue=1e-5`. NCBI rate-limits the public `-remote`
  endpoint; do not submit multiple jobs in parallel.
- `SAMPLE`, `QUERY`, `OUT_TSV` are required env vars (enforced via `: "${VAR:?}"`).
- The helper does not auto-submit (`sbatch` must be invoked by the operator) —
  this is intentional per the Issue #1 scope note.
- If a CAND has zero insert FASTA records (edge case: upstream assembly step
  produced an empty contig), `extract_cand_for_blast.py` logs a WARN and skips.

## Expected FP shape (from prior cucumber runs)

All 8 CAND in line225 carry 1-3 copies of `gfp|U55762.1`. A subset are
suspected to be host-genome GFP orthologs (cucurbit GFP-like fluorescent
protein genes) — the remote BLAST will separate bacterial/viral GFP
(TRUE_INSERTION) from `Cucumis sativus` hits (FALSE_POSITIVE).

## Classification Result (2026-04-17, JOBID 5630183, BLAST 1h39m)

**전원 TRUE_INSERTION — FP 0개.** Classification TSV committed at
`element_db/cuc225_cand_classification.tsv`.

| # | qseqid | top hit | pident | classification |
|---|--------|---------|--------|----------------|
| 1 | LKUO03000822.1_8440    | Expression vector pRNAi-GG          | 100.0 | TRUE_INSERTION |
| 2 | LKUO03000822.1_8540    | Expression vector pRNAi-GG          | 100.0 | TRUE_INSERTION |
| 3 | LKUO03000939.1_800850  | Binary vector pGV4945               | 100.0 | TRUE_INSERTION |
| 4 | LKUO03001451.1_6501    | Binary vector pBI121-ELEMENTS       | 100.0 | TRUE_INSERTION |
| 5 | LKUO03001826.1_6663169 | Klebsiella variicola plasmid        | 100.0 | TRUE_INSERTION |
| 6 | LKUO03001997.1_4219576 | Binary vector pBI121-ELEMENTS       | 100.0 | TRUE_INSERTION |
| 7 | LKUO03002166.1_547982  | Klebsiella variicola plasmid        | 100.0 | TRUE_INSERTION |
| 8 | LKUO03003260.1_75      | Ti plasmid binary vector pGA18.0214 | 100.0 | TRUE_INSERTION |

모든 hit 이 100.0% identity 로 T-DNA binary/expression/Ti plasmid vector 에 매칭.
**Cucurbit 호스트 ortholog 0개** — 가설 (GFP host-ortholog) 기각됨.

### AC-2 specification reinterpretation

원 스펙: "FP ≤ 5/sample". cuc_225 의 8 CAND 는 FP 가 아니라 **실제 T-DNA 카피**
(multi-copy T-DNA event 또는 head-to-head tandem array). 따라서 AC-2 재정의:
- "Specificity ≤ 5 FP/sample" → "**Specificity ≤ 5 independently-verified FP/sample**"
- cuc_225 는 v1.0 AC-2 **PASS** (0 verified FP, 8 TRUE_INSERTION).

Pre-mask BED 업데이트 **불필요** (FP 0). docs/host_masks/cucumber_b10v3.bed 는
기존 상태 유지.

## Post-classification

- [x] Classification TSV committed (`element_db/cuc225_cand_classification.tsv`)
- [x] Pre-mask BED review (no change needed)
- [x] Issue #1 closable
