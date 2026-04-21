# Validation & Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `feature/sample-construct-db` 머지 (`8bb8b0a`) 이후의 세 가지 회귀(rice Chr3:16,439,674 / tomato A2_3 ch01:91,002,744 / soybean AtYUCCA6 positive-site count)를 end-to-end 재실행으로 확인하고, OOM/TIMEOUT로 보류된 네 샘플(cucumber_line212/224/225, soybean_UGT72E3)을 완료하고, bug.md OPEN-1~5 항목을 닫는다.

**Architecture:** 새 코드 변경 없음 — 이미 `main`에 머지된 Task 10/11 fix (commits `4d3f2bb`, `10e3b9b`)가 실제로 동작하는지 검증. `run_pipeline.py --steps 4b,5` 재실행 → verdict 파일 비교 → 결과에 따라 bug.md 업데이트. cucumber/soybean은 메모리/타임만 조정해서 재제출.

**Tech Stack:** Bash + SLURM (partition `cpu-s2-core-0`, account `cpu-s2-pgl-0`), Python 3.11 (`run_pipeline.py`, `scripts/s04b_construct_assembly.py`, `scripts/s05_insert_assembly.py`), micromamba env `redgene`.

---

## Context Snapshot (2026-04-15)

**Branch:** `main` @ `729835d`. Worktree `.worktrees/sample-construct-db`는 머지 후 제거됨.

**유닛 테스트:** `pytest tests/ -v` → 10/10 PASS (`tests/test_extra_element_db.py`, `tests/test_s04b_construct_assembly.py`).

**보류 중인 stash:** `stash@{0}: On main: main-local-session-edits-pre-merge` — pre-merge 시점의 `README.md`, `resume.md`, `run_batch_other.sh`, `run_batch_tomato.sh` 편집. 머지 이후 최신 상태에서 재검토 필요.

**재검증 대상 회귀 (bug.md OPEN-3~5, resume.md status 표):**

| Sample | Expected outcome | 관련 commit |
|--------|-----------------|-------------|
| `rice_G281` | Chr3:16,439,674 restored as `CANDIDATE` after `_should_replace` merge rule | `4d3f2bb` |
| `tomato_Cas9_A2_3` | ch01:91,002,744 remains `CANDIDATE`; s04b rescue already spot-checked | `6de9f75`, `d56b18f` |
| `soybean_AtYUCCA6` | Transgene-positive site count `1,116 → ~baseline-scale` (spot-check: contigs 1,345 → 6) | `10e3b9b` |

**주의 — SLURM 환경 변수 (BUG-16):** 사용자 쉘에 `SBATCH_PARTITION=cpu-s1-pgl-0`이 설정되어 있을 수 있음. 이 환경변수는 `#SBATCH --partition` 지시자보다 우선순위가 높아서, 의도한 파티션으로 보내지 않음. 이 plan의 모든 `sbatch` 호출은 **CLI 플래그로 명시적 전달** — `sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 ...`.

---

## Phase 1: Housekeeping Before Reruns

### Task 1: Inspect & resolve pre-merge stash

**의존성:** 없음. 제일 먼저 실행해서 작업 트리를 깨끗하게 만들고 시작한다.

**Files:**
- Inspect: `stash@{0}` (paths `README.md`, `resume.md`, `run_batch_other.sh`, `run_batch_tomato.sh`)

- [ ] **Step 1: stash 내용을 요약 diff로 본다**

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
git stash show stash@{0} --stat
git stash show -p stash@{0}
```

Expected: 4개 파일 편집 내역. 현재 파일과 비교해서 이미 반영된 것인지, 놓친 편집인지 판단.

- [ ] **Step 2: 현재 파일과 stash 간 3-way diff**

```bash
# 각 파일이 이미 동일한 내용을 갖고 있는지 빠르게 확인
for f in README.md resume.md run_batch_other.sh run_batch_tomato.sh; do
  echo "=== $f (stash vs. HEAD) ==="
  git diff stash@{0} -- "$f" | head -40
done
```

- [ ] **Step 3: 판단 & 액션**

- 모든 내용이 이미 반영되어 있거나 더 이상 유효하지 않으면 → `git stash drop stash@{0}`.
- 반영되지 않은 유효한 편집이 있으면 → 해당 파일만 `git checkout stash@{0} -- <file>`로 꺼낸 뒤 `git diff` 확인하고 commit.
- 판단이 애매하면 drop하지 말고 그대로 두고 다음 task로 넘어간다 (plan 마지막 Task 10에서 다시 본다).

- [ ] **Step 4: 결정 근거 기록**

`bug.md`의 "Known pre-existing issues" 섹션 끝에 한 줄 추가:

```markdown
### OPEN-6 — pre-merge stash `main-local-session-edits-pre-merge` 처리 결과
- 2026-04-15: [drop / partial cherry-pick / still pending] — [한 줄 근거].
```

---

### Task 2: 유닛 테스트 재확인 (clean baseline)

**의존성:** Task 1 완료.

- [ ] **Step 1: pytest 실행**

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene
pytest tests/ -v 2>&1 | tee /tmp/pytest_validation_plan.log
```

Expected:
```
tests/test_extra_element_db.py ........ [N passed]
tests/test_s04b_construct_assembly.py ... [N passed]
============ 10 passed in <N>s ============
```

- [ ] **Step 2: 실패 시 중단**

테스트가 실패하면 Phase 2 이후 재실행에 의미가 없다. 원인 파악 후 수정 commit을 내고 다시 실행.

---

## Phase 2: Core Validation Reruns

모든 재실행은 **`--steps 4b,5 --no-remote-blast`**로 제한. s01~s04는 이미 완료되었고 bottleneck인 s04 host mapping을 재실행하지 않는다. `--no-remote-blast`는 NCBI 원격 BLAST를 건너뛴다 (회귀 검증 목적이라 필요 없음).

각 task의 sbatch 스크립트는 **heredoc** 대신 `run_rerun_*.sh` 파일로 저장한다 (재현성 + 실패 시 resubmit 용이).

### Task 3: rice_G281 rerun — restore Chr3:16,439,674

**의존성:** Task 2 완료.

**Files:**
- Create: `run_rerun_rice_validation.sh`
- Read: `results/rice_G281/s05_insert_assembly/insertion_Chr3_16439674_report.txt` (재실행 후)

**Context:** Chr3:16,439,674는 2-copy head-to-head T-DNA. `4d3f2bb` 이전에는 univec-element_db bitscore tie에서 univec이 이겨서 transgene-positive classification에서 drop됨.

- [ ] **Step 1: rerun 스크립트 작성**

```bash
cat > run_rerun_rice_validation.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=rg_rice_revalidate
#SBATCH --output=results/rg_rice_revalidate_%j.out
#SBATCH --error=results/rg_rice_revalidate_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene
python run_pipeline.py --sample rice_G281 --steps 4b,5 --threads 8 --no-remote-blast
EOF
chmod +x run_rerun_rice_validation.sh
```

- [ ] **Step 2: baseline (현재) verdict 스냅샷을 남긴다**

```bash
mkdir -p docs/superpowers/runs
grep -h "^Verdict:" results/rice_G281/s05_insert_assembly/insertion_*_report.txt 2>/dev/null \
    | sort | uniq -c | sort -rn \
    > docs/superpowers/runs/2026-04-15-rice-baseline.txt
cat docs/superpowers/runs/2026-04-15-rice-baseline.txt
```

이 파일은 나중에 Task 9에서 diff 기준으로 사용.

- [ ] **Step 3: SLURM 제출 (CLI로 partition 명시 — BUG-16 회피)**

```bash
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 \
       --time=6:00:00 --mem=32G --cpus-per-task=8 \
       --chdir="$PWD" ./run_rerun_rice_validation.sh
```

Expected output: `Submitted batch job <JOBID>` — JOBID를 기록.

- [ ] **Step 4: 진행 모니터링**

```bash
# 대략 30~90 min 소요. 완료 신호는 squeue에서 사라지는 것.
squeue -u $USER --name=rg_rice_revalidate
# 스트리밍 확인:
tail -f results/rg_rice_revalidate_<JOBID>.out
```

- [ ] **Step 5: Chr3:16,439,674 복원 확인**

```bash
grep "^Verdict:" results/rice_G281/s05_insert_assembly/insertion_Chr3_16439674_report.txt
```

Expected: `Verdict: CANDIDATE` (또는 CANDIDATE_HIGH_CONF류; FALSE_POSITIVE / UNKNOWN이면 실패).

- [ ] **Step 6: 전체 verdict 분포 스냅샷**

```bash
grep -h "^Verdict:" results/rice_G281/s05_insert_assembly/insertion_*_report.txt \
    | awk -F' —' '{print $1}' | sort | uniq -c | sort -rn \
    > docs/superpowers/runs/2026-04-15-rice-postfix.txt
diff docs/superpowers/runs/2026-04-15-rice-baseline.txt \
     docs/superpowers/runs/2026-04-15-rice-postfix.txt
```

- [ ] **Step 7: 실패 시 디버깅 경로**

만약 Chr3:16,439,674가 여전히 CANDIDATE가 아니면:
  1. `results/rice_G281/s04b_construct_asm/contigs.fasta`가 비어 있지 않은지 확인.
  2. `results/rice_G281/s05_insert_assembly/element_annotation.tsv`에 Chr3:16,439,674 엔트리와 어떤 DB hit가 기록되었는지 확인.
  3. `scripts/s05_insert_assembly.py`의 `_should_replace` 호출부에 로그 추가 후 단일 사이트만 재실행.

---

### Task 4: tomato_Cas9_A2_3 rerun — preserve ch01:91,002,744

**의존성:** Task 3 제출 후 (병렬 가능 — 다른 sample, 별도 SLURM job).

**Files:**
- Create: `run_rerun_a2_3_validation.sh`
- Read: `results/tomato_Cas9_A2_3/s05_insert_assembly/insertion_*_report.txt`

**Context:** A2_3의 T-DNA는 SLM_r2.0ch01:91,002,744 (config.yaml `expected.tdna_insertion_pos: 65107378`와 위치가 다르다는 점은 이미 세션 중 확인됨 — config 쪽이 outdated). s04b rescue가 이 사이트에서 NODE_1 100%/124bp 매치를 만들어내는 게 spot-check로 확인되었으나 end-to-end verdict는 미검증.

- [ ] **Step 1: rerun 스크립트 작성**

```bash
cat > run_rerun_a2_3_validation.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=rg_a2_3_revalidate
#SBATCH --output=results/rg_a2_3_revalidate_%j.out
#SBATCH --error=results/rg_a2_3_revalidate_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene
python run_pipeline.py --sample tomato_Cas9_A2_3 --steps 4b,5 --threads 8 --no-remote-blast
EOF
chmod +x run_rerun_a2_3_validation.sh
```

- [ ] **Step 2: baseline verdict 스냅샷**

```bash
grep -h "^Verdict:" results/tomato_Cas9_A2_3/s05_insert_assembly/insertion_*_report.txt 2>/dev/null \
    | awk -F' —' '{print $1}' | sort | uniq -c | sort -rn \
    > docs/superpowers/runs/2026-04-15-a2_3-baseline.txt
cat docs/superpowers/runs/2026-04-15-a2_3-baseline.txt
```

- [ ] **Step 3: SLURM 제출**

```bash
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 \
       --time=8:00:00 --mem=48G --cpus-per-task=8 \
       --chdir="$PWD" ./run_rerun_a2_3_validation.sh
```

토마토 genome은 833Mbp — rice (374Mbp)보다 약 2.2배. s05 자체는 host genome size의 영향이 적지만, BLAST step에서 일부 delay 있을 수 있으므로 wall time을 6h → 8h로 늘림.

- [ ] **Step 4: ch01:91,002,744 검증**

```bash
# 정확한 파일명이 tolerance window 때문에 약간 다를 수 있음
ls results/tomato_Cas9_A2_3/s05_insert_assembly/insertion_SLM_r2.0ch01_9100*_report.txt
grep "^Verdict:" results/tomato_Cas9_A2_3/s05_insert_assembly/insertion_SLM_r2.0ch01_9100*_report.txt
```

Expected: 최소 1개 `CANDIDATE` (BUG-7 회귀: 0 candidates면 실패).

- [ ] **Step 5: 전체 verdict 분포**

```bash
grep -h "^Verdict:" results/tomato_Cas9_A2_3/s05_insert_assembly/insertion_*_report.txt \
    | awk -F' —' '{print $1}' | sort | uniq -c | sort -rn \
    > docs/superpowers/runs/2026-04-15-a2_3-postfix.txt
```

---

### Task 5: soybean_AtYUCCA6 rerun — Task 11 filter validation

**의존성:** Task 3 제출 후 (병렬 가능).

**Files:**
- Create: `run_rerun_atyucca6_validation.sh`
- Read: `results/soybean_AtYUCCA6/s04b_construct_asm/contigs.fasta` (필터된 contig 수)
- Read: `results/soybean_AtYUCCA6/s05_insert_assembly/insertion_*_report.txt`

**Context:** BUG-5 fix 검증. 이전 job 5627124는 2h28m에 중단 (1,116 positive sites × 10 min/site > 24h limit). s04b contig 필터(≥90%/≥200bp vs `common_payload.fa` + `gmo_combined_db.fa`)가 1,345 → 6 contigs로 줄였으므로, Phase 1.5 positive site 수도 baseline(45)에 근접해야 함.

- [ ] **Step 1: rerun 스크립트 작성**

```bash
cat > run_rerun_atyucca6_validation.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=rg_atyucca6_revalidate
#SBATCH --output=results/rg_atyucca6_revalidate_%j.out
#SBATCH --error=results/rg_atyucca6_revalidate_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene
python run_pipeline.py --sample soybean_AtYUCCA6 --steps 4b,5 --threads 8 --no-remote-blast
EOF
chmod +x run_rerun_atyucca6_validation.sh
```

- [ ] **Step 2: s04b 재실행 전 현재 contig 수 기록**

```bash
if [ -s results/soybean_AtYUCCA6/s04b_construct_asm/contigs.fasta ]; then
  n_filt=$(grep -c "^>" results/soybean_AtYUCCA6/s04b_construct_asm/contigs.fasta)
  n_all=$(grep -c "^>" results/soybean_AtYUCCA6/s04b_construct_asm/contigs_all.fasta 2>/dev/null || echo NA)
  echo "AtYUCCA6 pre-rerun: filtered=$n_filt all=$n_all" \
       > docs/superpowers/runs/2026-04-15-atyucca6-baseline.txt
else
  echo "AtYUCCA6 pre-rerun: no s04b output yet" \
       > docs/superpowers/runs/2026-04-15-atyucca6-baseline.txt
fi
cat docs/superpowers/runs/2026-04-15-atyucca6-baseline.txt
```

- [ ] **Step 3: SLURM 제출 (메모리 48G, time 12h — 이전 24h TIMEOUT 경험상 여유 있게)**

```bash
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 \
       --time=12:00:00 --mem=48G --cpus-per-task=8 \
       --chdir="$PWD" ./run_rerun_atyucca6_validation.sh
```

- [ ] **Step 4: s04b 필터 결과 확인 (job 완료 전 중간에도 확인 가능)**

```bash
grep -c "^>" results/soybean_AtYUCCA6/s04b_construct_asm/contigs.fasta
grep -c "^>" results/soybean_AtYUCCA6/s04b_construct_asm/contigs_all.fasta
```

Expected: filtered ≈ 6 (±몇 개는 SPAdes stochasticity로 허용), `contigs_all`은 1,300~1,500 범위.

- [ ] **Step 5: Phase 1.5 transgene-positive site count 확인**

```bash
# s05가 stderr로 찍는 "Transgene-positive (assemble): N" 라인
grep "Transgene-positive" results/rg_atyucca6_revalidate_*.err results/rg_atyucca6_revalidate_*.out 2>/dev/null
```

Expected: 수십~수백 대 (baseline 45 대비 크게 벗어나지 않음). ≥1,000이면 필터 실패로 판단.

- [ ] **Step 6: 최종 verdict 분포**

```bash
grep -h "^Verdict:" results/soybean_AtYUCCA6/s05_insert_assembly/insertion_*_report.txt \
    | awk -F' —' '{print $1}' | sort | uniq -c | sort -rn \
    > docs/superpowers/runs/2026-04-15-atyucca6-postfix.txt
cat docs/superpowers/runs/2026-04-15-atyucca6-postfix.txt
```

Expected: 최소 1개 이상 `CANDIDATE` (BUG-4 fix로 `annotate_insert`의 extra_dbs가 element annotation을 채워주므로).

---

## Phase 3: OOM / TIMEOUT Retries

bug.md OPEN-1 (soybean_UGT72E3), OPEN-2 (cucumber line212/224/225) 해결.

### Task 6: cucumber_line212 / 224 / 225 rerun @ 64G

**의존성:** Phase 2의 검증들이 **진행 중이어도 병렬 실행 가능**. 세 샘플은 서로 독립이고 64G/16 CPU로 분리 제출한다.

**Files:**
- Create: `run_rerun_cucumber_oom.sh`

**Context:** 직전 batch job 5626023에서 exit code 9 (OOM kill) 세 번. 5626608을 64G로 재제출했으나 상태 미확인. 이번에는 확실한 flag 세트로 다시 돌린다.

- [ ] **Step 1: 이전 job 상태 먼저 확인**

```bash
sacct -j 5626608 --format=JobID,State,ExitCode,Elapsed,MaxRSS -n 2>/dev/null || echo "job not in accounting"
# 이미 성공했으면 Task 6 skip
ls -la results/cucumber_line212/s05_insert_assembly/ 2>/dev/null | head
ls -la results/cucumber_line224/s05_insert_assembly/ 2>/dev/null | head
ls -la results/cucumber_line225/s05_insert_assembly/ 2>/dev/null | head
```

세 샘플 모두 `insertion_*_report.txt` 파일이 존재하고 job 5626608이 COMPLETED면 → Step 6 (verdict 요약)으로 skip.

- [ ] **Step 2: rerun 스크립트 작성 (샘플별 개별 제출)**

```bash
cat > run_rerun_cucumber_oom.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=rg_cuc_oom
#SBATCH --output=results/rg_cuc_oom_%A_%a.out
#SBATCH --error=results/rg_cuc_oom_%A_%a.err
#SBATCH --array=0-2

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene
SAMPLES=(cucumber_line212 cucumber_line224 cucumber_line225)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "=== Running $SAMPLE ==="
python run_pipeline.py --sample "$SAMPLE" --steps 4b,5 --threads 16 --no-remote-blast
EOF
chmod +x run_rerun_cucumber_oom.sh
```

array job을 쓰면 세 개가 동시에 queue에 들어가고 각각 독립적으로 실패/성공을 보고함.

- [ ] **Step 3: SLURM 제출**

```bash
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 \
       --time=16:00:00 --mem=64G --cpus-per-task=16 \
       --chdir="$PWD" ./run_rerun_cucumber_oom.sh
```

- [ ] **Step 4: OOM 여부 모니터링 (job 완료 후)**

```bash
sacct -j <JOBID> --format=JobID,State,ExitCode,Elapsed,MaxRSS -n
```

Expected: 세 array task 모두 `COMPLETED`, MaxRSS는 64G 이하. exit code 9 재발 시 → 96G로 한 번 더 재시도하되, BUG-5 스타일의 positive site inflation이 원인일 가능성을 먼저 점검 (Task 5와 같은 로직).

- [ ] **Step 5: 라인별 T-DNA 검출 확인**

config.yaml의 expected:
- line212 → Chr6, single T-DNA
- line224 → Chr2, single T-DNA (disrupts G6838)
- line225 → Chr2, 2 T-DNA copies + backbone

```bash
for sample in cucumber_line212 cucumber_line224 cucumber_line225; do
  echo "=== $sample ==="
  grep -l "^Verdict: CANDIDATE" results/${sample}/s05_insert_assembly/insertion_*_report.txt 2>/dev/null \
    | sed 's|.*/insertion_||; s|_report.txt$||'
done
```

- [ ] **Step 6: verdict 요약 기록**

```bash
for sample in cucumber_line212 cucumber_line224 cucumber_line225; do
  echo "=== $sample ==="
  grep -h "^Verdict:" results/${sample}/s05_insert_assembly/insertion_*_report.txt 2>/dev/null \
      | awk -F' —' '{print $1}' | sort | uniq -c | sort -rn
done > docs/superpowers/runs/2026-04-15-cucumber-postfix.txt
cat docs/superpowers/runs/2026-04-15-cucumber-postfix.txt
```

---

### Task 7: soybean_UGT72E3 rerun — leverage Task 11 filter

**의존성:** Task 5 (AtYUCCA6)가 **완료되어 결과가 나와 있으면** 필터 효과 예상치를 조정할 수 있음. 대기할 필요는 없고, Phase 2와 병렬로 제출 가능.

**Files:**
- Create: `run_rerun_ugt72e3.sh`

**Context:** OPEN-1. job 5626560이 64G/24h TIMEOUT. BUG-5 fix (Task 11 contig filter)가 적용된 상태에서 재실행하면 Phase 1.5 positive site 수가 줄어 assemble queue를 통과할 가능성이 높다. 단, host가 soybean (1.1 Gbp, AtYUCCA6와 동일 genome)이라 여전히 시간이 오래 걸릴 수 있으므로 time은 24h로 시작.

- [ ] **Step 1: rerun 스크립트 작성**

```bash
cat > run_rerun_ugt72e3.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=rg_ugt72e3_retry
#SBATCH --output=results/rg_ugt72e3_retry_%j.out
#SBATCH --error=results/rg_ugt72e3_retry_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene
python run_pipeline.py --sample soybean_UGT72E3 --steps 4b,5 --threads 16 --no-remote-blast
EOF
chmod +x run_rerun_ugt72e3.sh
```

- [ ] **Step 2: SLURM 제출**

```bash
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 \
       --time=24:00:00 --mem=64G --cpus-per-task=16 \
       --chdir="$PWD" ./run_rerun_ugt72e3.sh
```

- [ ] **Step 3: s04b 필터 효과 중간 점검**

```bash
# job 시작 후 ~1h 이내에 s04b 완료됨
grep -c "^>" results/soybean_UGT72E3/s04b_construct_asm/contigs.fasta
grep -c "^>" results/soybean_UGT72E3/s04b_construct_asm/contigs_all.fasta
```

`contigs.fasta`에 100개 이상이 있으면 필터가 soybean/UGT72E3 construct에 대해 충분히 엄격하지 않다는 뜻 — Task 9 디버깅 대상에 추가.

- [ ] **Step 4: Phase 1.5 site count**

```bash
grep "Transgene-positive" results/rg_ugt72e3_retry_*.err results/rg_ugt72e3_retry_*.out 2>/dev/null
```

Expected: 수백 대 이하. 1,000+이면 24h 내 완료 불가 가능성 높음 → 현재 job은 그대로 두고 Task 9 디버깅으로 넘어가서 필터 임계값(현재 90%/200bp) 조정 검토.

- [ ] **Step 5: 완료 후 CANDIDATE 확인**

```bash
grep -h "^Verdict:" results/soybean_UGT72E3/s05_insert_assembly/insertion_*_report.txt 2>/dev/null \
    | awk -F' —' '{print $1}' | sort | uniq -c | sort -rn \
    > docs/superpowers/runs/2026-04-15-ugt72e3-postfix.txt
cat docs/superpowers/runs/2026-04-15-ugt72e3-postfix.txt
```

Expected: 최소 1개 `CANDIDATE` (bar marker + 35S promoter + UGT72E3 CDS 조합).

---

## Phase 4: Tracking Document Updates

### Task 8: bug.md 의 OPEN-1 ~ OPEN-5 확정

**의존성:** Task 3, 4, 5, 6, 7 모두 완료. 즉 Phase 2 + Phase 3 전체.

**Files:**
- Modify: `bug.md` ("Known pre-existing issues" 섹션)
- Modify: `docs/superpowers/runs/2026-04-14-regression.tsv` (또는 `2026-04-15-regression.tsv`를 새로 만듦)
- Read: `docs/superpowers/runs/2026-04-15-*.txt` (Phase 2-3에서 만든 snapshot들)

- [ ] **Step 1: 각 OPEN 항목에 결과 append**

`bug.md`의 OPEN-1 ~ OPEN-5 각각에 결과 줄을 추가. 템플릿:

```markdown
### OPEN-1 — `soybean_UGT72E3` s05 TIMEOUT at 64 GB / 24 h
- Job 5626100 hit 12 h limit, resubmitted as 5626560 at 64 G/24 h — still running at session end.
- Root cause probably same as AtYUCCA6: Phase 1.5 over-selects transgene-positive sites. Task 11 filter should help once applied.
- **2026-04-15 rerun (job <JOBID>):** [COMPLETED / TIMEOUT / OOM]. Phase 1.5 positives: <N>. Final verdicts: <CANDIDATE x / FP x / UNKNOWN x>. [closed / still open — reason].
```

OPEN-2 (cucumber), OPEN-3 (rice), OPEN-4 (A2_3), OPEN-5 (AtYUCCA6) 각각 동일 스타일.

- [ ] **Step 2: regression TSV 업데이트**

`docs/superpowers/runs/2026-04-14-regression.tsv`의 스키마를 확인 후, 2026-04-15 재실행 결과 row를 append. 스키마가 rice/A2_3/AtYUCCA6 verdict 델타만 담고 있으면, cucumber/UGT72E3는 별도 파일(`2026-04-15-oom-retry.tsv`)로 분리해도 됨 — 기존 파일의 컬럼과 맞지 않으면 덮어쓰지 말 것.

```bash
head -1 docs/superpowers/runs/2026-04-14-regression.tsv
# 같은 컬럼 체계에 맞춰서 append
```

- [ ] **Step 3: commit**

```bash
git add bug.md docs/superpowers/runs/
git commit -m "bug.md: close OPEN-1..5 with 2026-04-15 revalidation outcomes"
```

---

### Task 9: Diagnostic follow-up — 실패한 검증 사례만

**의존성:** Task 8 (OPEN 항목 결정 결과).

**적용 조건:** Task 3-7 중 **기대와 다른 결과**가 나온 항목이 있을 때만 수행. 모두 expected대로 나왔으면 Task 9는 skip하고 Phase 5로 진행.

- [ ] **Step 1: 실패한 케이스를 목록화**

```bash
# Task 8 commit 시점의 bug.md를 재읽어서 "still open"이나 예상 밖 결과를 가진 항목만 수집
grep -A2 "2026-04-15 rerun" bug.md | grep -B1 "still open\|unexpected\|FAILED"
```

- [ ] **Step 2: 각 실패 케이스에 대해 가설 → 1단계 디버깅**

흔한 케이스와 1차 조치 (새 Task로 분리할지는 이 단계에서 판단):

| 실패 패턴 | 1차 점검 지점 |
|----------|--------------|
| rice Chr3:16,439,674가 여전히 CANDIDATE 아님 | `element_annotation.tsv`에서 해당 사이트가 univec만 매칭하는지, `contigs.fasta`가 비었는지 |
| A2_3 ch01:91,002,744 0 candidates | `find_softclip_junctions` 재호출 로그에서 cluster_window/MAPQ 필터 확인 (BUG-7 회귀) |
| AtYUCCA6 positive count 여전히 >500 | `contigs.fasta`에 남아 있는 6개 contig의 염기서열을 BLAST로 수동 검사, 필터 임계값(90%/200bp) 재평가 |
| cucumber OOM 재발 | 어느 단계인지 (`Phase 1.5` vs `targeted assembly`) slurmstepd 로그로 확인 |
| UGT72E3 TIMEOUT 재발 | `contigs.fasta` 건수 확인. >100이면 soybean host에 대한 필터 세밀화 필요 |

각 조사 건은 **별도 commit**으로 결과를 `bug.md` 새 엔트리 (BUG-18+, OPEN-7+)에 기록.

- [ ] **Step 3: 코드 수정이 필요하면 별도 plan으로 분리**

이 plan은 validation-first이므로 fix 구현은 별도 plan으로 쓴다. 새 plan 파일 경로 예시:
`docs/superpowers/plans/2026-04-15-<specific-fix>.md`

---

## Phase 5: Deferred Visualization & Coverage Sensitivity

resume.md "Open items carried over"의 viz / coverage 항목은 검증이 끝난 뒤에 합쳐서 처리한다. Phase 1-4가 깨끗하게 닫히지 않으면 Phase 5는 이 plan 밖으로 밀어내는 것이 안전.

### Task 10: Tomato CRISPR editing visualization

**의존성:** Phase 2 완료 (A2_3 재실행 후 s04/s06 결과가 모두 최신).

**Files:**
- Script: `scripts/viz/plot_editing_profile.py`
- Script: `scripts/viz/plot_editing_effects.py`
- Input: `results/tomato_Cas9_A2_*/s06_indel/editing_sites.tsv`, `results/tomato_Cas9_A2_*/s06_indel/grna_targets.tsv`, `results/tomato_WT/s04_host_map/tomato_WT_host.bam`
- Output: `results/tomato_Cas9_A2_{1,2,3}_editing_profile.png`, `results/tomato_Cas9_A2_{1,2,3}_editing_effects.png`

- [ ] **Step 1: s06 결과 존재 여부 확인**

```bash
for sample in tomato_Cas9_A2_1 tomato_Cas9_A2_2 tomato_Cas9_A2_3; do
  echo "=== $sample ==="
  ls results/${sample}/s06_indel/ 2>/dev/null
done
```

`editing_sites.tsv`가 존재하지 않는 샘플이 있으면 해당 샘플에 대해 `run_pipeline.py --sample <sample> --steps 6` 먼저 실행 (s06_indel은 짧으므로 login node에서 실행 가능).

- [ ] **Step 2: editing profile plot 생성 (세 샘플)**

```bash
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene
cd /data/gpfs/assoc/pgl/develop/redgene

for sample in tomato_Cas9_A2_1 tomato_Cas9_A2_2 tomato_Cas9_A2_3; do
  python scripts/viz/plot_editing_profile.py \
    --treatment-bam results/${sample}/s04_host_map/${sample}_host.bam \
    --wt-bam       results/tomato_WT/s04_host_map/tomato_WT_host.bam \
    --host-ref     db/SLM_r2.0.pmol.fasta \
    --grna-targets   results/${sample}/s06_indel/grna_targets.tsv \
    --editing-sites  results/${sample}/s06_indel/editing_sites.tsv \
    --sample-name ${sample} --outdir results
done
```

Expected output: `results/${sample}_editing_profile.png` 세 개.

- [ ] **Step 3: editing effect annotation plot 생성**

```bash
for sample in tomato_Cas9_A2_1 tomato_Cas9_A2_2 tomato_Cas9_A2_3; do
  python scripts/viz/plot_editing_effects.py \
    --editing-sites  results/${sample}/s06_indel/editing_sites.tsv \
    --gff            db/SLM_r2.0.gff3.gz \
    --host-ref       db/SLM_r2.0.pmol.fasta \
    --sample-name ${sample} --outdir results
done
```

- [ ] **Step 4: Ground truth 대조 (config.yaml expected)**

Expected (Seol et al. 2025):
- **SlAMS** (Solyc08g062780, chr08): 6bp deletion
- **SlPHD_MS1** (Solyc04g008420, chr04): 9bp deletion
- A2_1: heterozygous (Cas9 removed)
- A2_2: homozygous
- A2_3: heterozygous

`editing_sites.tsv`의 해당 locus 행과 frameshift 호출이 플롯에 표시되었는지 PNG 육안 확인.

---

### Task 11: Sample summary figures

**의존성:** Phase 2-3 완료 (각 샘플의 s05 / s06 / s07 결과 존재).

**Files:**
- Script: `scripts/viz/plot_sample_summary.py`
- Output: `results/<sample>_summary.png`

- [ ] **Step 1: 대상 샘플 목록**

```bash
SAMPLES=(
  rice_G281
  tomato_Cas9_A2_1 tomato_Cas9_A2_2 tomato_Cas9_A2_3
  cucumber_line212 cucumber_line224 cucumber_line225
  soybean_AtYUCCA6 soybean_UGT72E3
  corn_ND207
)
```

`corn_ND207`은 resume.md 시점에 상태 미확인이므로 s05가 끝난 것만 포함. 확인:

```bash
for s in "${SAMPLES[@]}"; do
  if [ -d "results/$s/s05_insert_assembly" ] && \
     ls results/$s/s05_insert_assembly/insertion_*_report.txt &>/dev/null; then
    echo "$s: READY"
  else
    echo "$s: skip (no s05)"
  fi
done
```

- [ ] **Step 2: summary plot 생성 (READY 샘플만)**

```bash
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene
cd /data/gpfs/assoc/pgl/develop/redgene

for s in "${SAMPLES[@]}"; do
  if ls results/$s/s05_insert_assembly/insertion_*_report.txt &>/dev/null; then
    echo "=== $s ==="
    python scripts/viz/plot_sample_summary.py --sample-name "$s" --outdir results \
      || echo "  plot failed — continue"
  fi
done
```

`|| continue`로 하나가 실패해도 나머지 샘플은 계속 처리.

- [ ] **Step 3: 결과 파일 확인**

```bash
ls -la results/*_summary.png
```

---

### Task 12: Coverage sensitivity 실행 (optional — 이번 plan에서 끝까지 안 해도 됨)

**의존성:** Phase 2-4 모두 완료, Phase 5 Task 10-11도 완료. 여기까지 여력이 남았을 때만.

**Context:** config.yaml에는 이미 subsample variant가 정의되어 있음 (`rice_G281_{15,10,5,3}x`, `cucumber_line224_{15,10,5,3}x`, `soybean_UGT72E3_{15,10,5,3}x`, `tomato_A2_1_{15,10,5,3}x`, `tomato_A2_3_{5,3}x`). 각각 s01-s05 전체를 돌려야 해서 시간/리소스 소모가 크다.

- [ ] **Step 1: 대상 범위를 좁힌다**

모든 subsample을 돌리지 말고, **coverage cliff 검증에 핵심적인 것만**:
- `rice_G281_{10,5,3}x`: 10x/5x 경계 확인 (CLAUDE.md에 "≥15x for rice" 기준 기재되어 있음 — 10x는 회색지대)
- `cucumber_line224_{10,5,3}x`: 10x 기준 검증
- `soybean_UGT72E3_{10,5,3}x`: 10x 기준 검증

위 9개 샘플 외에는 이 plan에서 돌리지 않는다.

- [ ] **Step 2: batch 스크립트 작성**

```bash
cat > run_coverage_sensitivity_2026_04_15.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=rg_cov_sens
#SBATCH --output=results/rg_cov_sens_%A_%a.out
#SBATCH --error=results/rg_cov_sens_%A_%a.err
#SBATCH --array=0-8

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene
SAMPLES=(
  rice_G281_10x rice_G281_5x rice_G281_3x
  cucumber_line224_10x cucumber_line224_5x cucumber_line224_3x
  soybean_UGT72E3_10x soybean_UGT72E3_5x soybean_UGT72E3_3x
)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "=== $SAMPLE ==="
python run_pipeline.py --sample "$SAMPLE" --steps 1-5 --threads 16 --no-remote-blast
EOF
chmod +x run_coverage_sensitivity_2026_04_15.sh
```

Coverage subsample은 `--steps 1-5` (s01/s02 새로 돌림 — 이 샘플들은 아직 실행 안 됨).

- [ ] **Step 3: 제출 — memory 64G, time 24h**

```bash
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 \
       --time=24:00:00 --mem=64G --cpus-per-task=16 \
       --chdir="$PWD" ./run_coverage_sensitivity_2026_04_15.sh
```

- [ ] **Step 4: 결과 정리 TSV**

```bash
OUT=docs/superpowers/runs/2026-04-15-coverage-sensitivity.tsv
echo -e "sample\tcoverage\tcandidate_sites\tground_truth_detected" > $OUT

for host in rice_G281 cucumber_line224 soybean_UGT72E3; do
  for cov in 10x 5x 3x; do
    sub="${host}_${cov}"
    if [ -d "results/${sub}/s05_insert_assembly" ]; then
      n=$(grep -l "^Verdict: CANDIDATE" results/${sub}/s05_insert_assembly/insertion_*_report.txt 2>/dev/null | wc -l)
      echo -e "${host}\t${cov}\t${n}\tTBD" >> $OUT
    fi
  done
done
cat $OUT
```

ground_truth_detected 컬럼은 각 host의 config.yaml `expected.insertion_chr/pos` 정보를 보고 수동으로 채운다 (자동화는 이 plan 밖).

---

## Dependency Graph

```
Phase 1: Housekeeping
  Task 1: stash 처리 ─┐
  Task 2: 유닛 테스트 ─┤── (Phase 1 완료)
                       │
  ┌────────────────────┼─────────────────────────────┐
  │                    │                             │
  ▼                    ▼                             ▼
Phase 2 (validation) Phase 3 (OOM/TIMEOUT retry)  (둘은 병렬)
  Task 3: rice         Task 6: cucumber x3
  Task 4: A2_3         Task 7: UGT72E3
  Task 5: AtYUCCA6         │
         │                 │
         └──── Task 8: bug.md OPEN close-out
                    │
                    ▼
             Task 9: diagnostic follow-up (실패 케이스만)
                    │
                    ▼
Phase 5 (deferred — optional):
  Task 10: CRISPR viz    (A2 series s06 완료 전제)
  Task 11: summary figs  (Phase 2-3 완료 전제)
  Task 12: coverage sensitivity (여력 있을 때만)
```

**Critical path:** Task 1 → Task 3 (rice 검증이 가장 빠른 signal) → Task 8.
Task 4/5/6/7은 모두 Task 3와 병렬 제출 가능 (별도 SLURM job).

---

## Time Estimates

| Task | 실행 시간 | Compute | 비고 |
|------|----------|---------|------|
| Task 1-2 | ~10 min | login node | stash inspection + pytest |
| Task 3 | ~1h | 1 sbatch | rice `--steps 4b,5` (cached s04 BAM 전제) |
| Task 4 | ~2h | 1 sbatch | tomato 833Mbp |
| Task 5 | ~4h | 1 sbatch | soybean + 필터된 contig로 진행 |
| Task 6 | ~4h (array) | 3 array tasks | cucumber x3 병렬 |
| Task 7 | ~12h | 1 sbatch | soybean UGT72E3 |
| Task 8-9 | ~30 min | login node | doc 업데이트 + 필요 시 디버깅 |
| Task 10 | ~10 min | login node | matplotlib viz |
| Task 11 | ~15 min | login node | 9-10개 summary plots |
| Task 12 | ~24h | 1 array job | 9 coverage variants (optional) |

**총 wall time:** Task 3-7이 병렬이면 ~12-16h (UGT72E3가 longest pole). Task 12를 포함하면 +24h.

---

## Self-Review Checklist

- [x] **Spec coverage:** resume.md "Validation status" 표의 3개 샘플 모두 Phase 2 Task 3-5에 대응. "Open items carried over"의 cucumber/UGT72E3 → Task 6-7. viz → Task 10-11. coverage sensitivity → Task 12.
- [x] **Placeholder scan:** 모든 sbatch 명령어에 명시적 partition/account/time/mem. JOBID는 "제출 후 기록" 형태로 명시.
- [x] **SLURM BUG-16 회피:** 모든 `sbatch` 호출이 CLI로 `--partition=cpu-s2-core-0 --account=cpu-s2-pgl-0`를 전달.
- [x] **파일 경로:** validation output은 `docs/superpowers/runs/2026-04-15-*.txt`로 통일. rerun shell은 `run_rerun_*.sh`로 통일.
- [x] **스텝 단위:** 각 task가 2-5 분 단위 스텝 (스크립트 작성 → 스냅샷 → 제출 → 모니터 → 검증 → 기록).

---

## Execution Handoff

Plan 저장 위치: `docs/superpowers/plans/2026-04-15-validation-and-cleanup.md`

실행 옵션:

1. **Subagent-Driven (권장)** — Phase 1부터 task별로 fresh subagent 디스패치. Task 간에 결과 확인 후 다음 task 진행. 장점: rice/A2_3/AtYUCCA6가 동시에 queue에 들어가는 동안 main session은 다른 일 할 수 있음.
2. **Inline Execution** — 현재 세션에서 executing-plans skill로 순차 진행. 장점: sbatch 제출/확인이 한 세션에 끝남.

어느 방식으로 갈까요?
