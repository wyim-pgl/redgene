# Resume — RedGene Pipeline

**Date:** 2026-04-15
**Branch:** `feature/sample-construct-db` (HEAD `10e3b9b`), 20 commits ahead of main
**Working dir:** `/data/gpfs/assoc/pgl/develop/redgene/.worktrees/sample-construct-db` (isolated worktree; `db/`, `test_data/`, `data/`, `results/` symlink back to the main path)

---

## Feature summary

Extends the pipeline so that **sample-specific transgene payloads (AtYUCCA6 gene, bar marker, etc.) that aren't in the shared `element_db/gmo_combined_db.fa` can still be annotated**. Two mechanisms:

1. **`element_db/common_payload.fa`** — curated FASTA of 9 bacterial/viral-origin transgene markers (bar, nptII, hpt, gusA, gfp, egfp, P-CaMV35S, P-nos, T-ocs), fetched from NCBI via `element_db/build_common_payload.sh`. Wired as an always-on `--common-payload-db` flag on s05.
2. **New step 4b** (`scripts/s04b_construct_assembly.py`) — de novo SPAdes assembly of construct-hitting reads from s03. The resulting contigs are passed to s05 via `--extra-element-db`, giving s05 a per-sample reference it can BLAST against. After Task 11 the wrapper also **filters the contigs** by marker-DB hits (≥90% / ≥200 bp) so only payload-bearing contigs are forwarded — prevents 25× positive-site inflation seen on raw contigs.

`run_pipeline.py` now registers step 4b between 4 and 5; `--steps 1-5` expands to `1,2,3,4,4b,5`. `pipeline.common_payload_db` in `config.yaml` names the always-on DB. All three extras are threaded through `_batch_check_element_hits` (Phase 1.5) and `annotate_insert` (Phase 3) via `extra_dbs: list[Path] | None`.

---

## Commits on the branch (oldest → newest)

```
a152d7f  gitignore: add .worktrees/
0c90fe6  Add common transgene payload accession manifest
76c15ee  Add common_payload.fa with bar/nptII/hpt/gusA/gfp/AtYUCCA6
b61b7f1  common_payload: drop duplicate T-nos, route summary to stdout
b02be31  common_payload: atomic publish + trap cleanup
c8d9ec4  Add s04b step: de novo SPAdes assembly of construct reads
34fe8e7  s04b: robustify empty detection + scratch cleanup + test name
ec7cf10  s05: accept --extra-element-db for per-sample construct DB
8fc9057  s05: disambiguate per-DB blast output file to prevent stem collision
6de9f75  run_pipeline: register s04b and forward contigs into s05
5ca9e44  run_pipeline: document step 4b and correct memory comment
c3c96b0  Wire common_payload DB into s05 as always-on supplementary
ae62ea2  test(s05): exercise unrelated-clip negative case in multi-DB test
0c5a9c1  common_payload: drop AtYUCCA6 (host-origin plant gene)
66fb4e1  Record rice_G281 + A2_3 regression deltas after extra-DB plumbing
0e3cbc0  regression: diagnose rice Chr3:16,439,674 bitscore-tie exclusion
3998bdf  Document step 4b and common_payload DB
d56b18f  s05: extend annotate_insert (Phase 3) with extra_dbs list
01f4654  regression: document AtYUCCA6 25x transgene-positive inflation
4d3f2bb  s05: element_db beats univec regardless of bitscore in classify_site_tiers
10e3b9b  s04b: filter contigs by marker-DB hits to curb extra_db inflation
```

20 commits. Two reverts already landed upstream (`4cec387` site-discovery fix → reverted in `68d0e52`; `f78912b` borders-heuristic reclass → reverted in `02e1579`) before this feature branch forked.

---

## Tests

`tests/test_extra_element_db.py` + `tests/test_s04b_construct_assembly.py` — **10 passing**:

- `_batch_check_element_hits` with single extra DB, multiple extra DBs, and no extras
- `annotate_insert` with extra_dbs
- `_should_replace` priority (element_db > univec regardless of bitscore; same-source bitscore-best with incumbent ties)
- s04b SPAdes-failure-on-tiny-input, empty-input short-circuit
- `_filter_contigs_by_markers` keeps marker-positive, drops random, handles missing marker DB, `--no-filter` keeps raw

Run: `cd <worktree> && eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene && pytest tests/ -v`.

---

## End-to-end verification status

| Sample | Result |
|--------|--------|
| rice_G281 (pre-Task 10) | Regression: Chr3:16,439,674 lost — both clips univec-only at equal bitscore, s04b contig covered the region but tied. **Task 10 (commit `4d3f2bb`) fixes the tie-break; needs rerun to confirm.** |
| tomato_Cas9_A2_3 | Positive signal from partial run: Phase 1 VALIDATED `SLM_r2.0ch01:91,002,744` with 3p=NODE_1 from s04b contigs at 100%/124bp (strictly longer than univec 76bp, so the existing strict `>` merge already worked here). Full rerun still pending. |
| soybean_AtYUCCA6 | Phase 1.5 inflated to 1,116 positive sites (25× baseline 45) on raw s04b contigs; run killed. **Task 11 (commit `10e3b9b`) filters contigs 1,345 → 6 (99.6% reduction); needs rerun to confirm final verdict distribution.** |
| cucumber_line212/224/225 | Pre-existing OOM at 32GB — separate issue, not touched by this branch. |
| soybean_UGT72E3 | Pre-existing 24h timeout on 64GB — separate issue, not touched. |

Regression snapshots are in `docs/superpowers/runs/2026-04-14-regression.tsv`.

---

## Deferred work (not blockers for finishing)

- Rerun s05 on rice_G281, tomato_Cas9_A2_3, soybean_AtYUCCA6 end-to-end with the current HEAD to confirm the two fixes land as expected. Each takes 30 min–several hours on SLURM.
- Cucumber + UGT72E3 reruns at 64GB (pre-existing issues, tracked in the main worktree's task list).
- No PR yet — awaiting `finishing-a-development-branch` decision (merge / push / keep / discard).

---

## Reference — how to test this branch quickly

```bash
cd /data/gpfs/assoc/pgl/develop/redgene/.worktrees/sample-construct-db
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

# Unit tests (under 10s)
pytest tests/ -v

# End-to-end on a known sample (requires SLURM or long local run)
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 --time=24:00:00 \
    --mem=64G --cpus-per-task=8 --chdir="$PWD" \
    --wrap="eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate redgene \
            && python run_pipeline.py --sample soybean_AtYUCCA6 --steps 4b,5 --threads 8 --no-remote-blast"

# Inspect verdicts after completion
grep -h "^Verdict:" results/soybean_AtYUCCA6/s05_insert_assembly/insertion_*_report.txt \
    | awk -F' —' '{print $1}' | sort | uniq -c
```

---

## Plan document

`docs/superpowers/plans/2026-04-14-sample-construct-assembly-and-payload-db.md` — full 9-task plan + accepted deviations. Tasks 10 and 11 were added mid-flight as follow-ups after Task 7 and Task 8 surfaced two design gaps; both are now committed on this branch.
