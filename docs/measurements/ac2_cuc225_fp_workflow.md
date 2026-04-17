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

## Post-classification

Update `docs/measurements/coverage_sensitivity.md` FP column once the TSV is
resolved. Do NOT close Issue #1 until the 8-row classification TSV is
committed and the T10 pre-mask BED is re-evaluated.
