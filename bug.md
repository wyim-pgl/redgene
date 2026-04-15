# RedGene Session Bug Log

All bugs found and their resolution status across the `feature/sample-construct-db` development cycle (2026-04-12 → 2026-04-15). Format per entry: **severity**, **what broke**, **root cause**, **fix / status**, and a commit or artifact pointer.

Only bugs that actually affected runs or reviews are listed — cosmetic / nitpick feedback that went straight to a commit isn't tracked here.

---

## Code bugs (fixed on `main`)

### BUG-1 — `s06_indel.py` NameError on every editing site
- **Severity:** high (crashes step 6 on first site)
- **Symptom:** `NameError: name 'indel_seq' is not defined. Did you mean: 'indels'?` at `scripts/s06_indel.py:298`
- **Root cause:** a rename of `indel_seq` → `indel_seq_upper` missed one call site.
- **Fix:** commit `1fe8bfe` (re-applied after an unrelated revert reverted it a second time).
- **How found:** job 5626101 failed in 8 s during tomato step 6.

### BUG-2 — `_batch_check_element_hits` stem-collision overwrites BLAST output
- **Severity:** medium (silent data loss when primary + extra DBs share a filename stem)
- **Symptom:** second iteration's blast output file overwrote the first's before it was parsed.
- **Root cause:** `blast_out = workdir / f"_clip_blast_{db.stem}.tsv"` — not unique across iterations.
- **Fix:** commit `8fc9057` — prefix filename with loop index (`_clip_blast_{i}_{db.stem}.tsv`).
- **How found:** Task 4 quality review.

### BUG-3 — `classify_site_tiers` bitscore-tie silently drops real T-DNAs
- **Severity:** high (cost us rice_G281 Chr3:16,439,674 in Task 7 regression)
- **Symptom:** real T-DNA site excluded from transgene-positive classification even though s04b contig covered the same region at 100% identity.
- **Root cause:** bitscore-best merge condition was a strict `>`. When the primary `transgene_db` (univec) and an extra DB (element_db) tied, the first iteration won; univec-tagged hit passed, and the element_db-required positive policy dropped the site.
- **Fix:** commit `4d3f2bb` — `_should_replace(existing, new_src, new_bit)` helper enforces element_db > univec regardless of bitscore; same-source ties keep the incumbent.
- **How found:** Task 7 rice rerun after step 4b and diagnostic BLAST showing the contig covered AF234315.1:6222-6323 at 100%/33bp.

### BUG-4 — `annotate_insert` (Phase 3) never used `extra_dbs`
- **Severity:** high (AtYUCCA6 UNKNOWN sites stayed UNKNOWN even after s04b contigs were wired in)
- **Symptom:** `element_annotation.tsv` came out empty on AtYUCCA6 despite s04b NODE_1 containing bar / gfp / P-CaMV35S / P-nos / T-ocs at ≥99% identity.
- **Root cause:** Tasks 4 and 6 extended `_batch_check_element_hits` and `classify_site_tiers` (Phase 1.5) to accept `extra_dbs`, but left `annotate_insert` calling `_run_local_blast(insert, element_db)` against the primary DB only.
- **Fix:** commit `d56b18f` — `annotate_insert(..., extra_dbs: list[Path] | None = None)` runs one BLAST per DB and merges hits; `_run_local_blast` gained a `tag` kwarg to disambiguate intermediate filenames.
- **How found:** Task 8 AtYUCCA6 run showed 30 UNKNOWN / 15 FP / 0 CANDIDATE (= baseline), empty `element_annotation.tsv`.

### BUG-5 — s04b raw contigs inflate Phase 1.5 positive sites 25×
- **Severity:** high (AtYUCCA6 job 5627124 aborted at 2h28m; 1,116 sites vs baseline 45)
- **Symptom:** `Transgene-positive (assemble): 1116` on soybean_AtYUCCA6 when passing the raw 1,345-contig SPAdes output as `--extra-element-db`. At ~10 min/site this exceeds 24 h SLURM limit.
- **Root cause:** SPAdes contigs from construct-hitting reads contain generic plant/viral sequences (P-CaMV35S, generic NOS-terminator regions, etc.) that cross-react with broad host features in soybean/rice/tomato genomes.
- **Fix:** commit `10e3b9b` — `s04b_construct_assembly.py` now filters contigs by marker-DB hits (≥90% / ≥200 bp against `common_payload.fa` + `gmo_combined_db.fa`) before writing `contigs.fasta`; unfiltered copy preserved as `contigs_all.fasta`. AtYUCCA6 spot-check: 1,345 → 6 contigs (99.6 % reduction).
- **How found:** Task 8 end-to-end rerun.

### BUG-6 — `UNKNOWN → CANDIDATE_LOW_CONF` heuristic produces false positives
- **Severity:** high (would silently promote FP sites to candidates if merged)
- **Symptom:** heuristic `borders ≥ 4 AND largest_gap ≥ 1000 bp` promoted two AtYUCCA6 UNKNOWN sites that NCBI remote BLAST later showed were 100% host DNA (NC_038245.2:43,783,085 matched *Glycine soja* chr1 at 99.67 % / 1,522 bp; NC_016089.4:14,246,877 was a multi-chrom host repeat).
- **Root cause:** `_blast_insert_vs_host`'s `largest_gap` compares against a single host reference. Reference-assembly gaps / cultivar-specific regions appear as "non-host" even when the same sequence is present in another soybean assembly in NCBI.
- **Fix:** the f78912b heuristic was reverted in `02e1579`. The cross-reference check is not automatic — users now rely on remote BLAST or curated DBs. The accompanying `scripts/postprocess_unknown_reclass.py` is kept as a template for future better-gated heuristics.
- **How found:** NCBI remote BLAST of the two promoted gaps during in-session validation.

### BUG-7 — `find_softclip_junctions` cluster-window widening drops real T-DNAs
- **Severity:** medium (would have caused rice and A2_3 regressions if merged)
- **Symptom:** proposed `cluster_window = 50 → 300` + `min_mapq = 20` filter removed A2_3's `ch01:91,002,744` candidate entirely.
- **Root cause:** T-DNA junctions sometimes rely on microhomology anchors that register as MAPQ 0-19. The new MAPQ filter dropped those reads; the wider cluster window by itself didn't compensate.
- **Fix:** commit `4cec387` reverted in `68d0e52`. Current site-discovery parameters stay at `cluster_window=50`, no MAPQ filter.
- **How found:** A2_3 rerun showed 0 candidates; BAM inspection revealed the junction's anchor reads lived below MAPQ 20.

---

## Documentation / data bugs (fixed in same session)

### BUG-8 — `NM_122942.3` is not AtYUCCA6
- **Severity:** low (plan typo; not in any production script)
- **Symptom:** the original plan/manifest listed AtYUCCA6 as `NM_122942.3`, which actually resolves to a TIR-NBS-LRR disease-resistance gene.
- **Fix:** Task 2 implementer caught it and corrected to `NM_122473.3` (Arabidopsis YUC6 / AT5G25620). Later removed entirely in `0c5a9c1` because its plant origin causes host cross-reaction in every tested species.
- **How found:** Task 2 implementer review.

### BUG-9 — duplicate `V00087.1` row for P-nos and T-nos
- **Severity:** low (produced two identical FASTA entries with different headers)
- **Symptom:** `common_payload.fa` had 11 sequences but P-nos and T-nos were identical bodies under two names — same accession fetched twice.
- **Fix:** commit `b61b7f1` removed the T-nos row. A proper solution would sub-region the V00087.1 nos operon (promoter vs. terminator) via `efetch -seq_start/-seq_stop`, left as a potential future improvement.
- **How found:** Task 2 spec review.

### BUG-10 — `build_common_payload.sh` leaves partial output + leaks TMP on abort
- **Severity:** low (no user-facing failure yet; robustness bug)
- **Symptom:** `set -e` abort in the efetch loop left a partial `common_payload.fa` on disk with no indication of truncation, and `TMP=$(mktemp -d)` was never cleaned on non-zero exit.
- **Fix:** commit `b02be31` — single `TMPOUT=$(mktemp)` + `trap 'rm -f "$TMPOUT"' EXIT`, `mv $TMPOUT $OUT` only on successful loop completion.
- **How found:** Task 2 code-quality review.

### BUG-11 — `s04b_construct_assembly.py` crashes on missing input file
- **Severity:** low (intended behaviour was a clean empty output, not a traceback)
- **Symptom:** `_is_empty_fastq(Path('/no/such/path'))` raised `FileNotFoundError` instead of returning `True`.
- **Fix:** commit `34fe8e7` — nested `try/except OSError` so missing / unreadable is treated as empty.
- **How found:** Task 3 code-quality review.

### BUG-12 — s04b scratch dir `_spades_run/` leaks on SPAdes failure
- **Severity:** low (disk usage; multi-GB per sample when SPAdes aborts)
- **Symptom:** on SPAdes non-zero exit the wrapper wrote empty `contigs.fasta` and returned 0, but `shutil.rmtree(spades_out)` only ran on the success path.
- **Fix:** commit `34fe8e7` — cleanup called in both branches.
- **How found:** Task 3 code-quality review.

### BUG-13 — `_batch_check_element_hits` does not tolerate `extra_dbs = []`
- **Severity:** low (edge case; `main()` always passes a list, but the nullable-list contract was ambiguous)
- **Symptom:** no bug in practice; review flagged that `None` and `[]` take different code paths that happen to converge.
- **Fix:** left as-is but covered by the negative-case pytest `test_batch_check_element_hits_without_extras_ignores_extra_seq` (commit `ae62ea2`).

---

## Infrastructure / environment bugs (workarounds applied)

### BUG-14 — worktree `db/` is a real directory, not a symlink
- **Severity:** medium (s05 failures at runtime until ref files are linked)
- **Symptom:** initial worktree setup symlinked `test_data/`, `data/`, `results/` back to the main path, but `db/` was a real directory that shared only a handful of tracked files. Worktree runs failed at host-BWA step when the tomato / soybean / cucumber references were missing.
- **Workaround:** symlinked `SLM_r2.0.pmol.*`, `Gmax_v4.0.fa*`, etc. individually in-session. A future hardening would either make `db/` a symlink or checkpoint a list of required files for each sample.
- **How found:** Task 4 rice smoke test; Task 7 A2_3 rerun; Task 8 AtYUCCA6 step 5 crash (commit context in regression TSV).

### BUG-15 — `db/transgene_db.fa` can be left zero-byte from an earlier failed run
- **Severity:** medium (blastn exits with "Database memory map file error")
- **Symptom:** `transgene_db.fa` is built by the pipeline if missing, but a previously-interrupted run left a 0-byte file that passed the `exists()` check.
- **Workaround:** deleted the 0-byte file so the regeneration path fired on next run.
- **Suggested real fix:** size-check the file before trusting it, or write to `<path>.tmp` then `mv` atomically.

### BUG-16 — SLURM `SBATCH_PARTITION` env vars override `#SBATCH` directives
- **Severity:** medium (jobs routed to a saturated partition despite the script saying otherwise)
- **Symptom:** two `rg_s05fix_*` submissions landed on `cpu-s1-pgl-0` even though their scripts specified `cpu-s2-core-0`. SLURM priority rules mean env > script-level directive.
- **Workaround:** pass partition/account/time explicitly on the `sbatch` CLI. User env `SBATCH_PARTITION=cpu-s1-pgl-0` is the culprit and should probably be cleared before submitting to other partitions.
- **How found:** session observation on jobs 5626028/5626029 before resubmitting as 5626030/5626031.

### BUG-17 — `pytest` not installed in `redgene` micromamba env
- **Severity:** low (blocks running tests from that env)
- **Symptom:** fresh `redgene` env didn't have `pytest`; tests were failing to import `pysam` because of a wrong env activation, masking the real issue.
- **Workaround:** `micromamba install -n redgene -c conda-forge pytest`. Should be added to the environment manifest for future creates.

---

## Known pre-existing issues (not from this session, still open)

### OPEN-1 — `soybean_UGT72E3` s05 TIMEOUT at 64 GB / 24 h
- Job 5626100 hit 12 h limit, resubmitted as 5626560 at 64 G/24 h — still running at session end.
- Root cause probably same as AtYUCCA6: Phase 1.5 over-selects transgene-positive sites. Task 11 filter should help once applied.

### OPEN-2 — cucumber line212/224/225 s05 OOM at 32 GB
- Batch job 5626023 lost three samples to exit-code-9 kills.
- Resubmitted at 64 G (5626608) earlier in session; status unverified.

### OPEN-3 — rice_G281 end-to-end rerun with Task 10 + 11 fixes
- Chr3:16,439,674 should now survive Phase 1.5 after `4d3f2bb`. Not yet verified by a full rerun on the merged `main`.

### OPEN-4 — tomato_Cas9_A2_3 end-to-end rerun
- Partial run confirmed the s04b rescue mechanism works (NODE_1 at 100%/124 bp). Full rerun to Phase 4 verdict still pending.

### OPEN-5 — soybean_AtYUCCA6 end-to-end rerun with Task 11 filter
- Spot-check showed contigs 1,345 → 6. Full s05 rerun should now finish in time and give a real verdict distribution.

---

## Session-level observations

- Every major bug in this session was surfaced by running the pipeline on a real sample, not by static review alone. Reviews caught cosmetic / robustness issues (BUG-2, -10, -11, -12) but missed the merge-priority bug (BUG-3) and the Phase 3 annotation gap (BUG-4) until end-to-end runs exposed them.
- The two reverts in-session (`68d0e52`, `02e1579`) demonstrate the value of always validating threshold / heuristic changes against a known-good sample before committing defaults.
- Worktree isolation was worth the setup cost: two different s05 reruns interfered because they shared the `results/` symlink, but no code regression escaped back to `main` until the explicit fast-forward merge.
