# Resume — RedGene Pipeline

**Date:** 2026-04-15
**Branch:** `main` @ `8bb8b0a` — `feature/sample-construct-db` merged (fast-forward)
**Working dir:** `/data/gpfs/assoc/pgl/develop/redgene`
**Worktree:** removed after merge (was at `.worktrees/sample-construct-db`)

---

## Session outcome — what changed on `main`

### New pipeline step `4b`
- `scripts/s04b_construct_assembly.py` — de novo SPAdes assembly of the construct-hitting reads produced by s03. Emits `results/<sample>/s04b_construct_asm/contigs.fasta` (filtered by marker-DB hits) + `contigs_all.fasta` (unfiltered debug copy).
- Registered in `run_pipeline.py`: `--steps 1-5` now expands to `1,2,3,4,4b,5`. Step 5 auto-appends `--extra-element-db <contigs.fasta>` when it exists and is non-empty.

### Shared `element_db/common_payload.fa`
- 9 canonical bacterial/viral transgene markers: `bar`, `nptII`, `hpt`, `gusA`, `gfp`, `egfp`, `P-CaMV35S`, `P-nos`, `T-ocs`.
- Fetched from NCBI by `element_db/build_common_payload.sh` (manifest in `common_payload_manifest.tsv`).
- Wired as the always-on `--common-payload-db` flag on s05 via `pipeline.common_payload_db` in `config.yaml`.

### s05 extended to consult multiple element DBs
- `_batch_check_element_hits` (Phase 1.5 positive classification) and `annotate_insert` (Phase 3 element annotation) now both take `extra_dbs: list[Path] | None`.
- `classify_site_tiers` merge rule uses a new `_should_replace` helper: **element_db hit beats univec regardless of bitscore**; same-source ties keep the incumbent. This restores sites whose clips match vector backbone at the same identity as a sample-specific contig (rice Chr3:16,439,674 was the motivating case).

### Pre-filter for s04b contigs
- Raw s04b contigs can cross-react with generic host features (AtYUCCA6 spike test: 1,345 contigs inflated Phase 1.5 positives from 45 → 1,116). The wrapper now keeps only contigs with at least one `blastn` hit ≥90 % / ≥200 bp against `common_payload.fa` + `gmo_combined_db.fa`. AtYUCCA6 spot-check: 1,345 → 6 contigs (99.6 % reduction).
- Escape hatch: `--no-filter` preserves legacy raw contigs.

### Tests
`tests/test_extra_element_db.py` + `tests/test_s04b_construct_assembly.py` — **10/10 passing** (pytest 9.0.3 in `redgene` micromamba env).

### Bug + regression log
- `bug.md` — 17 bugs + 5 open issues + session observations (committed `8bb8b0a`).
- `docs/superpowers/runs/2026-04-14-regression.tsv` — rice / A2_3 / AtYUCCA6 verdict deltas and diagnostic notes.

---

## Branch/commit map

`main` now contains everything from the former `feature/sample-construct-db` branch. Highlights (newest → oldest):

```
8bb8b0a  bug: log 17 bugs + 5 open issues from sample-construct-db cycle
de83e60  resume: capture sample-construct-db branch state for finish
10e3b9b  s04b: filter contigs by marker-DB hits to curb extra_db inflation          (Task 11)
4d3f2bb  s05: element_db beats univec regardless of bitscore in classify_site_tiers (Task 10)
01f4654  regression: document AtYUCCA6 25x transgene-positive inflation
d56b18f  s05: extend annotate_insert (Phase 3) with extra_dbs list                  (bug-4 fix)
3998bdf  Document step 4b and common_payload DB
0e3cbc0  regression: diagnose rice Chr3:16,439,674 bitscore-tie exclusion
66fb4e1  Record rice_G281 + A2_3 regression deltas after extra-DB plumbing
0c5a9c1  common_payload: drop AtYUCCA6 (host-origin plant gene)
ae62ea2  test(s05): exercise unrelated-clip negative case in multi-DB test
c3c96b0  Wire common_payload DB into s05 as always-on supplementary                 (Task 6)
5ca9e44  run_pipeline: document step 4b and correct memory comment
6de9f75  run_pipeline: register s04b and forward contigs into s05                   (Task 5)
8fc9057  s05: disambiguate per-DB blast output file to prevent stem collision
ec7cf10  s05: accept --extra-element-db for per-sample construct DB                 (Task 4)
34fe8e7  s04b: robustify empty detection + scratch cleanup + test name
c8d9ec4  Add s04b step: de novo SPAdes assembly of construct reads                  (Task 3)
b02be31  common_payload: atomic publish + trap cleanup
b61b7f1  common_payload: drop duplicate T-nos, route summary to stdout
76c15ee  Add common_payload.fa with bar/nptII/hpt/gusA/gfp/AtYUCCA6                 (Task 2)
0c90fe6  Add common transgene payload accession manifest                            (Task 1)
a152d7f  gitignore: add .worktrees/
1fe8bfe  Re-apply s06 indel_seq NameError fix after f78912b revert
02e1579  Revert "Reclassify UNKNOWN inserts as CANDIDATE_LOW_CONF ..."
```

Two in-session reverts (`68d0e52`, `02e1579`) kept two flawed heuristics out of `main`.

---

## End-to-end validation status

Unit tests all pass. Full-pipeline reruns to confirm the two late fixes are **not yet done** and are the first things to do next session.

| Sample | Expected outcome | Status |
|--------|-----------------|--------|
| rice_G281 | Chr3:16,439,674 restored as CANDIDATE after Task 10 (`_should_replace`) | ⏳ rerun needed |
| tomato_Cas9_A2_3 | ch01:91,002,744 remains CANDIDATE; s04b rescue already validated via partial run | ⏳ full rerun needed |
| soybean_AtYUCCA6 | Positive-site count drops from 1,116 → ~baseline-scale after Task 11 filter; UNKNOWN sites with real T-DNA promoted to CANDIDATE via `annotate_insert` extra_dbs | ⏳ rerun needed |

Suggested rerun command (each sample as a separate SLURM job):

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene

sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 --time=24:00:00 \
  --mem=64G --cpus-per-task=8 --chdir="$PWD" \
  --wrap="eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate redgene \
          && python run_pipeline.py --sample <sample> --steps 4b,5 --threads 8 --no-remote-blast"
```

---

## Open items carried over

Pre-existing issues outside this feature's scope, still pending:

- **cucumber_line212 / 224 / 225**: s05 OOM at 32 GB. Needs 64 GB rerun (job 5626608 was submitted earlier in session; final state not re-verified).
- **soybean_UGT72E3**: s05 TIMEOUT at 64 GB / 24 h (job 5626560). Task 11 contig filter may help reduce assembly queue.
- **CRISPR indel validation** for A2_1/A2_2/A2_3 already confirmed earlier in the broader session (SlPHD_MS1 9 bp del in A2_2 matches Seol et al. 2025 ground truth).
- **Step 7 copy number** done for all four tomato samples.
- **Visualization scripts** (`plot_editing_profile.py`, `plot_editing_effects.py`, `plot_sample_summary.py`) and coverage sensitivity batch still untouched.

---

## Reference — how to pick up from here

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene

git log --oneline -10                   # confirm state
pytest tests/ -v                        # 10 tests, <10 s
cat bug.md                              # open items + deferred fixes
cat docs/superpowers/runs/2026-04-14-regression.tsv   # verdict deltas

# Submit a validation run (rice first — smallest genome, fastest round-trip)
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 --time=6:00:00 \
  --mem=32G --cpus-per-task=8 --chdir="$PWD" \
  --wrap="eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate redgene \
          && python run_pipeline.py --sample rice_G281 --steps 4b,5 --threads 8 --no-remote-blast"

# Inspect verdicts after completion
grep -h "^Verdict:" results/rice_G281/s05_insert_assembly/insertion_*_report.txt \
    | awk -F' —' '{print $1}' | sort | uniq -c
grep "^Verdict:" results/rice_G281/s05_insert_assembly/insertion_Chr3_16439674_report.txt
```

Plan document and in-progress notes:

- `docs/superpowers/plans/2026-04-14-sample-construct-assembly-and-payload-db.md`
- `docs/superpowers/runs/2026-04-14-regression.tsv`

Stash left from the merge (main's pre-session edits to `README.md`, `resume.md`, `run_batch_other.sh`, `run_batch_tomato.sh`) is still on `stash@{0}` as `main-local-session-edits-pre-merge`. `git stash show` to inspect; `git stash drop` to discard if unneeded.
