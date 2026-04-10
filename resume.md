# Resume — RedGene Pipeline

**Date:** 2026-04-10
**Branch:** main
**Working dir:** /data/gpfs/assoc/pgl/develop/redgene

---

## What happened this session

### 1. Three optimizations applied to s05_insert_assembly.py

- **Stable site_id**: `insertion_1` → `insertion_{chrom}_{pos}` (e.g., `insertion_Chr3_16439674`). Duplicate positions get `_2`, `_3` suffix.
- **`--no-remote-blast`**: CLI flag to skip slow NCBI nt remote BLAST. Added to both `s05_insert_assembly.py` and `run_pipeline.py`.
- **max-rounds 15 → 8**: Reduced default assembly rounds. Most sites converge by round 3-6.

### 2. Rice G281 full pipeline run completed (job 5625160)

Results from optimized pipeline (steps 2-5):
- 8,690 candidate sites → 28 transgene-positive → assembly
- **3 CANDIDATE**, 5 FALSE_POSITIVE, 20 UNKNOWN
- Ground truth site **Chr3:16,439,674** correctly detected as CANDIDATE (13,471bp insert)
- Two extra CANDIDATEs (Chr3_31434949, Chr3_31443563) are false positives missed by existing filters

### 3. Filter D: Alternative host locus check (new)

Root cause analysis of the two extra CANDIDATEs:
- Both assembled inserts are host genomic DNA mis-assembled via CaMV 35S promoter/terminator homology
- Both map via minimap2 to Chr2:8432860 at MAPQ=60
- Existing filters missed them: Filter A (host% too low at 90% identity), Filter B (not flanking), Filter C (only 1 off-target chr)

New **Filter D** added: maps assembled insert to host via minimap2. If primary alignment covers ≥70% of insert at a different locus (≥10kb away or different chr) → FALSE_POSITIVE.

### 4. Previous session work (still in this commit)

- Pipeline reorganization 10 → 7 steps (Step 5 centric)
- UniVec plant vector integration (Step 2)
- Post-assembly FP filters A/B/C
- Code review fixes (stale paths, crash in viz scripts)

---

## Active jobs

| Job ID | Name | Purpose | Status |
|--------|------|---------|--------|
| **5625482** | rg_filterd | Rice G281 step 5 only (Filter D verification) | **RUNNING** on cpu-48 |

### Check commands
```bash
squeue -u wyim --name=rg_filterd
tail -20 results/rg_filterd_5625482.err
grep "Verdict" results/rice_G281/s05_insert_assembly/insertion_*_report.txt
```

---

## Expected verification results

After job 5625482 completes:
- `insertion_Chr3_16439674`: **CANDIDATE** (true insertion, unchanged)
- `insertion_Chr3_31434949`: **FALSE_POSITIVE** (Filter D: alt locus Chr2:8432860)
- `insertion_Chr3_31443563`: **FALSE_POSITIVE** (Filter D: alt locus Chr2:8432860)

---

## Pending work

1. **Verify Filter D results** from job 5625482
2. **Submit tomato + other batches** after rice validation
3. **Step 6 (CRISPR)** for tomato samples (needs WT BAM)
4. **Step 7 (copy number)** after step 5
5. **Address 20 UNKNOWN sites** — consider verdict logic for sites with no element annotations

---

## Resume command for next session

```
Read /data/gpfs/assoc/pgl/develop/redgene/resume.md
squeue -u wyim --name=rg_filterd,rg_rice,rg_tomato
grep "Verdict" results/rice_G281/s05_insert_assembly/insertion_*_report.txt
```
