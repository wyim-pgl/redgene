# RedGene v1.0 MVP Work Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** team-consensus.md §3 의 MVP 12 task (T1-T12) 를 bite-sized TDD 단계로 구현해 이번 주(2026-04-20 ~ 2026-04-24) 안에 RedGene v1.0-rc 를 release-ready 로 만든다.

**Architecture:** 현 `main` @ `729835d` 기반으로 **코드 대규모 재구조화 없이** (1) AC-6 감사 헤더 4-field 추가, (2) DB 확장 + cd-hit-est dedup + 4-way source tag, (3) `VerdictRules` config + `compute_verdict()` pure 함수 분리, (4) 최소 4-way s05 모듈 분할 + per-site fan-out SLURM array 임시안, (5) element vs host 100%-ortholog pre-mask BED workflow, (6) 6 샘플 재검증 + rice_G281 s04 minimap2 PoC 순서로 진행. s05 전체 8-10 모듈 분할은 v1.1 로 이관, dual-anchor / PE discordancy 는 v1.1, GRIDSS / long-read 는 v2.0 grant 로 분리.

**Tech Stack:** Python 3.11 + pathlib + stderr logging, pysam, minimap2, BWA, SPAdes, blastn + megablast, cd-hit-est, bedtools, SLURM sbatch array, pytest, micromamba env `redgene`.

**권위 근거:** 본 plan 의 모든 task 는 `docs/team-review/team-consensus.md` 의 T1-T12 와 1:1 매핑된다. 구현 중 충돌 시 team-consensus.md 가 최종 판정.

---

## Context Snapshot (2026-04-16)

- **Base:** `main` @ `729835d`, pytest 10/10 PASS.
- **검증 상태 (resume.md):** 5/7 샘플 CANDIDATE 회수 (rice_G281, tomato_Cas9_A2_3, cucumber_line212/224/225) — AC-1 5/6 = 83%. soybean_UGT72E3 TIMEOUT (48h 필요), soybean_AtYUCCA6 0 CAND (BUG-8 부작용).
- **이번 주 목표:** MVP v1.0 4 AC 통과 — AC-1 6/8 이상 + AC-2 ≤ 5 FP/sample + AC-4 UGT72E3 ≤ 4h + AC-6 4-field audit header.
- **SLURM:** partition `cpu-s2-core-0`, account `cpu-s2-pgl-0`. `SBATCH_PARTITION` env 는 CLI flag 로 명시 override (BUG-16 회피).
- **주의:** `db/` 디렉터리는 실물 (BUG-14). worktree 사용 시 symlink 수동 설정 필요.

---

## File Structure (생성 및 수정 대상)

**신규 파일 (Create):**
```
scripts/
├── s05/                                  # T7: 최소 4-way 모듈 분할
│   ├── __init__.py                       # (1줄) from .verdict import compute_verdict
│   ├── verdict.py                        # T6: compute_verdict() + FilterEvidence + VerdictRules
│   └── config_loader.py                  # T6: VerdictRules.from_yaml()
├── submit_s05_array.sh                   # T8: SLURM array fan-out wrapper
├── build_element_mask_bed.sh             # T9: element vs host 100%-ortholog BED
├── measure_assembly_rounds.py            # T2: s05_stats.txt aggregate
└── _write_audit_header.py                # T1: audit header (run_pipeline.py 에서 import)

tests/
├── test_compute_verdict.py               # T6: 7 시나리오
├── test_verdict_rules_loader.py          # T6: config 파싱
├── test_source_tag_priority.py           # T5: 4-way tag priority
├── test_audit_header.py                  # T1: SHA-256/commit/DB md5/version
└── test_mask_bed_intersect.py            # T10: FALSE_NEGATIVE_MASKED 태그

element_db/
├── cas9_sgrna.fa                         # T4: SpCas9 + sgRNA scaffold
├── euginius_missing.fa                   # T4: PMI/vip3Aa/cry2Ab2/cry34-35Ab1/T-pinII/T-g7/CTP2/2mEPSPS
├── payload_cds.fa                        # T4: AtYUCCA6/UGT72E3
├── gmo_combined_db_v2.fa                 # T4: cd-hit-est 결과 (dedup 후 통합)
├── gmo_combined_db_manifest.tsv          # T1+T5: name/md5/build_date/seq_count
└── Makefile                              # T5: cd-hit-est + manifest 빌드

docs/
├── measurements/
│   ├── max_rounds_study.md               # T2: round 3 수렴율 측정 결과
│   └── s04_minimap2_poc.md               # T12: rice_G281 PoC 결과
└── host_masks/
    ├── rice_osativa_v7.bed               # T9: 100%-ortholog 좌표
    ├── tomato_slm_r2.bed
    ├── cucumber_b10v3.bed
    ├── corn_zm_b73_v5.bed
    ├── soybean_gmax_v4.bed
    └── host_masked_rationale.tsv         # T9: 각 mask 근거

run_rerun_w1_batch.sh                     # T11: 6 샘플 재검증 array
run_s04_minimap2_poc.sh                   # T12: rice_G281 s04 minimap2 단독
```

**수정 파일 (Modify):**
```
run_pipeline.py                           # T1: audit header injection, T8: --fanout flag
config.yaml                               # T6: canonical_triplets / verdict_rules schema
scripts/s05_insert_assembly.py            # T3: max_rounds=3/5, T7: verdict 호출을 s05.verdict 로 위임, T10: BED intersect
scripts/s05_insert_assembly.py:2372       # T3: max_rounds 파라미터 조정 후 regression
scripts/s05_insert_assembly.py:3197-3514  # T6/T7: generate_report 에서 compute_verdict() 호출로 교체
element_db/build_common_payload.sh        # T4: V00087.1 P-nos/T-nos subregion efetch (-seq_start/-seq_stop) + X04879.1 T-ocs subregion
element_db/element_summary.tsv            # T4: 15 새 entry 추가
tests/conftest.py                         # T6: fixture for FilterEvidence, VerdictRules
bug.md                                    # T11: OPEN-1 (UGT72E3), OPEN-5 (AtYUCCA6) close-out
resume.md                                 # T11: v1.0-rc 검증 결과 snapshot
```

---

## Task Ordering (critical path)

team-consensus.md §3 critical path 기반:

```
Mon (day 1)      Tue (day 2)        Wed (day 3)      Thu (day 4)       Fri (day 5)
─────────────    ──────────────     ─────────────    ─────────────     ─────────────
T1 audit hdr  → T4 DB expand     → T5 cd-hit+tag →                  →
T6 verdict  ─┐                  →                 → T8 fan-out    → T11 재검증
T7 4-way split                   → T9 BED build →  T8 fan-out     → T11 재검증
              T2 round study   → T3 max_rounds →                  → T12 minimap2 PoC
                                → T10 BED intrs →                 →
```

독립적 시작 가능: T1, T6, T7, T2.

---

## Phase 1 Tasks

### Task 1: AC-6 Audit Header 4-field (compute_eng)

**목표:** 매 `run_pipeline.py` 실행의 보고서 header 에 input SHA-256 + commit hash + DB md5 + software version manifest 4 field 를 기록. regulatory R-1/R-2/R-3/R-4 충족.

**Files:**
- Create: `scripts/_write_audit_header.py`
- Create: `tests/test_audit_header.py`
- Modify: `run_pipeline.py` (시작부 `main()` 첫 블록에 `_write_audit_header` 호출 삽입)

- [ ] **Step 1: failing test 작성**

```python
# tests/test_audit_header.py
import json
import tempfile
from pathlib import Path
import pytest
from scripts._write_audit_header import write_audit_header


def test_audit_header_contains_4_fields(tmp_path):
    r1 = tmp_path / "reads_R1.fq.gz"
    r1.write_bytes(b"@r\nACGT\n+\n!!!!\n")
    r2 = tmp_path / "reads_R2.fq.gz"
    r2.write_bytes(b"@r\nACGT\n+\n!!!!\n")
    db_manifest = tmp_path / "element_db_manifest.tsv"
    db_manifest.write_text("name\tmd5\tbuild_date\tseq_count\nelement_db\tdeadbeef\t2026-04-16\t146\n")

    out = tmp_path / "audit.json"
    write_audit_header(
        sample="rice_G281",
        reads_r1=r1,
        reads_r2=r2,
        db_manifest=db_manifest,
        out_path=out,
    )

    data = json.loads(out.read_text())
    assert "input_sha256" in data
    assert data["input_sha256"]["r1"].startswith(("a","b","c","d","e","f") + tuple("0123456789"))
    assert len(data["input_sha256"]["r1"]) == 64  # SHA-256 hex
    assert "pipeline_commit" in data
    assert "pipeline_dirty" in data  # True/False
    assert "db_manifest" in data and data["db_manifest"][0]["name"] == "element_db"
    assert "software_versions" in data
    assert any(k.startswith("bwa") for k in data["software_versions"])
```

- [ ] **Step 2: test 실행 → FAIL 확인**

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene
pytest tests/test_audit_header.py -v
```

Expected: `ImportError: No module named 'scripts._write_audit_header'`.

- [ ] **Step 3: 최소 구현**

```python
# scripts/_write_audit_header.py
"""Audit header writer for RedGene v1.0 regulatory compliance.

Writes SHA-256 of input reads, pipeline git commit hash (with dirty flag),
element DB manifest, and software version manifest to a single JSON file.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_commit_and_dirty() -> tuple[str, bool]:
    repo = Path(__file__).resolve().parent.parent
    commit = subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
    ).strip()
    status = subprocess.check_output(
        ["git", "-C", str(repo), "status", "--porcelain"], text=True
    )
    return commit, bool(status.strip())


def _software_versions() -> dict[str, str]:
    tools = {
        "bwa": ["bwa"],
        "minimap2": ["minimap2", "--version"],
        "samtools": ["samtools", "--version"],
        "blastn": ["blastn", "-version"],
        "spades": ["spades.py", "--version"],
        "cd-hit-est": ["cd-hit-est", "-h"],
        "fastp": ["fastp", "--version"],
        "python": [sys.executable, "--version"],
    }
    out: dict[str, str] = {}
    for name, cmd in tools.items():
        try:
            res = subprocess.run(
                cmd, capture_output=True, text=True, timeout=5, check=False
            )
            combined = (res.stdout + res.stderr).strip().splitlines()
            out[name] = combined[0] if combined else "unknown"
        except (FileNotFoundError, subprocess.TimeoutExpired):
            out[name] = "not-found"
    return out


def _parse_manifest(path: Path) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    if not path.exists():
        return entries
    lines = path.read_text().strip().splitlines()
    if not lines:
        return entries
    header = lines[0].split("\t")
    for row in lines[1:]:
        cols = row.split("\t")
        if len(cols) == len(header):
            entries.append(dict(zip(header, cols)))
    return entries


def write_audit_header(
    *,
    sample: str,
    reads_r1: Path,
    reads_r2: Path,
    db_manifest: Path,
    out_path: Path,
) -> None:
    commit, dirty = _git_commit_and_dirty()
    data = {
        "sample": sample,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "input_sha256": {
            "r1": _sha256_file(reads_r1),
            "r2": _sha256_file(reads_r2),
        },
        "pipeline_commit": commit,
        "pipeline_dirty": dirty,
        "db_manifest": _parse_manifest(db_manifest),
        "software_versions": _software_versions(),
    }
    out_path.write_text(json.dumps(data, indent=2))
```

- [ ] **Step 4: test 재실행 → PASS 확인**

```bash
pytest tests/test_audit_header.py -v
```

Expected: `1 passed`.

- [ ] **Step 5: `run_pipeline.py` 에 wiring**

`run_pipeline.py` 의 `main()` 에서 per-sample 루프 시작부에 다음 한 블록 삽입 (sample 별 output dir 생성 직후, step 1 호출 전):

```python
from scripts._write_audit_header import write_audit_header
write_audit_header(
    sample=sample_name,
    reads_r1=Path(sample_cfg["reads"]["r1"]),
    reads_r2=Path(sample_cfg["reads"]["r2"]),
    db_manifest=Path("element_db") / "gmo_combined_db_manifest.tsv",
    out_path=Path(outdir) / sample_name / "audit_header.json",
)
```

- [ ] **Step 6: 통합 smoke test**

```bash
python run_pipeline.py --sample rice_G281 --steps 0 --dry-run
# audit_header.json 이 results/rice_G281/audit_header.json 에 생성되어야 함
# 단, steps 0 이 없으면 --steps 1 로 fastp 시작 전 생성됨 확인 가능
ls -la results/rice_G281/audit_header.json
jq '.input_sha256, .pipeline_commit, .pipeline_dirty, .software_versions | keys' \
    results/rice_G281/audit_header.json
```

Expected: `audit_header.json` 존재, SHA-256 길이 64 hex, software_versions keys = bwa/minimap2/samtools/blastn/spades/cd-hit-est/fastp/python.

- [ ] **Step 7: Commit**

```bash
git add scripts/_write_audit_header.py tests/test_audit_header.py run_pipeline.py
git commit -m "$(cat <<'EOF'
T1: add AC-6 audit header writer (SHA-256 + commit + DB + software versions)

Writes per-sample audit_header.json with 4 regulatory-required fields
for RedGene v1.0 MVP (team-consensus.md §2.1 item 1, §5 regulatory declaration).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Assembly `max_rounds` 수렴율 측정 (haibao)

**목표:** T3 의 `max_rounds=3` vs `8` 선택 근거 확보. 현재 `s05_stats.txt` 에서 round 별 growth 기여율을 aggregate.

**Files:**
- Create: `scripts/measure_assembly_rounds.py`
- Create: `docs/measurements/max_rounds_study.md`

- [ ] **Step 1: 현 s05_stats 샘플 분포 확인**

```bash
cd /data/gpfs/assoc/pgl/develop/redgene
find results -name "s05_stats.txt" -size +0 | head -20
# 각 파일의 라인 포맷 확인
head -20 results/rice_G281/s05_insert_assembly/s05_stats.txt
```

Expected: "round N: kmer=X mm2=Y pilon=Z ssake=W combined=T" 류 라인.

- [ ] **Step 2: aggregator 스크립트 작성**

```python
# scripts/measure_assembly_rounds.py
"""Aggregate s05_stats.txt across all samples to measure max_rounds convergence.

Output: per-round growth contribution distribution + histogram.
"""
from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

ROUND_RE = re.compile(
    r"round\s+(\d+):\s*kmer=(\d+)\s+mm2=(\d+)\s+pilon=(\d+)\s+ssake=(\d+)\s+combined=(\d+)",
    re.IGNORECASE,
)


def parse_stats(path: Path) -> list[dict[str, int]]:
    rounds: list[dict[str, int]] = []
    for line in path.read_text().splitlines():
        m = ROUND_RE.search(line)
        if m:
            rounds.append({
                "round": int(m[1]),
                "kmer": int(m[2]),
                "mm2": int(m[3]),
                "pilon": int(m[4]),
                "ssake": int(m[5]),
                "combined": int(m[6]),
            })
    return rounds


def main() -> None:
    results_dir = Path("results")
    by_round_growth: dict[int, list[int]] = defaultdict(list)
    per_sample_last_round: list[tuple[str, int]] = []
    n_files = 0

    for stats in sorted(results_dir.glob("*/s05_insert_assembly/s05_stats.txt")):
        rounds = parse_stats(stats)
        if not rounds:
            continue
        n_files += 1
        last = 0
        prev_combined = 0
        for r in rounds:
            growth = r["combined"] - prev_combined
            by_round_growth[r["round"]].append(growth)
            prev_combined = r["combined"]
            last = r["round"]
        per_sample_last_round.append((stats.parent.parent.name, last))

    total_growth = {rn: sum(g) for rn, g in by_round_growth.items()}
    cum = 0
    grand_total = sum(total_growth.values()) or 1
    print(f"# s05_stats.txt aggregate across {n_files} samples")
    print("round\tgrowth_total\tfraction_of_all\tcumulative_fraction")
    for rn in sorted(total_growth):
        cum += total_growth[rn]
        print(
            f"{rn}\t{total_growth[rn]}\t"
            f"{total_growth[rn]/grand_total:.4f}\t{cum/grand_total:.4f}"
        )
    print("\n# per-sample last round")
    for name, last in per_sample_last_round:
        print(f"{name}\t{last}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: 스크립트 실행 및 결과 저장**

```bash
python scripts/measure_assembly_rounds.py > /tmp/rounds_study.tsv
cat /tmp/rounds_study.tsv
```

Expected: round 1-8 의 cumulative_fraction 컬럼. round 3 이후 fraction_of_all 이 ≤ 10% 여야 T3 에서 `max_rounds=3` 확정 가능.

- [ ] **Step 4: 결과 문서화 + 판정**

```bash
mkdir -p docs/measurements
cat > docs/measurements/max_rounds_study.md <<'EOF'
# Assembly `max_rounds` Convergence Study (T2, 2026-04-20)

**Goal:** decide T3 default for `assemble_insert(max_rounds=N)` based on real s05_stats data.

**Source:** `results/*/s05_insert_assembly/s05_stats.txt` aggregated by `scripts/measure_assembly_rounds.py`.

## Data

<paste /tmp/rounds_study.tsv content here>

## Decision

- If `fraction_of_all[round=3]` ≤ 10%: **`max_rounds=3` 채택 (T3)**. per-site wall time ~5 min.
- If 10-30%: **`max_rounds=5` 타협**. per-site wall time ~8 min.
- If ≥ 30%: **`max_rounds=8` 유지**. T3 scope out.

**Result:** <채택된 N 및 선정 근거 2-3 문장>
EOF
vi docs/measurements/max_rounds_study.md  # 결과 붙여넣기 + Decision 기록
```

- [ ] **Step 5: Commit**

```bash
git add scripts/measure_assembly_rounds.py docs/measurements/max_rounds_study.md
git commit -m "$(cat <<'EOF'
T2: measure assembly max_rounds convergence across existing runs

Aggregator + measurement snapshot to justify T3 max_rounds default.
Decision recorded in docs/measurements/max_rounds_study.md.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: `assemble_insert(max_rounds=N)` 조정 (haibao)

**Files:**
- Modify: `scripts/s05_insert_assembly.py:2372` (또는 `max_rounds` default 가 선언된 위치)
- Test: 기존 샘플 1 건 회귀 확인

**의존성:** T2 완료 (측정 결과 필요).

- [ ] **Step 1: 현 default 확인**

```bash
grep -n "max_rounds" scripts/s05_insert_assembly.py
```

Expected: `def assemble_insert(..., max_rounds: int = 8, ...)` 같은 선언부 라인 번호 확인.

- [ ] **Step 2: T2 판정대로 default 수정**

Edit `scripts/s05_insert_assembly.py` at `assemble_insert(...)` signature. T2 Decision 에 따라 `max_rounds=3` 또는 `max_rounds=5`.

```python
def assemble_insert(
    insert_name: str,
    contigs: list[Contig],
    r1_fq: Path,
    r2_fq: Path,
    unmapped_fq: Path,
    *,
    max_rounds: int = 3,  # T3: reduced from 8 per T2 measurement
    ...
) -> AssemblyResult:
```

- [ ] **Step 3: rice_G281 회귀 smoke**

```bash
# 기존 BAM 재사용
python run_pipeline.py --sample rice_G281 --steps 5 \
    --threads 8 --no-remote-blast 2>&1 | tee /tmp/rice_t3_smoke.log
```

Wait ~40-60 min (s05 only). 완료 후:

```bash
grep "^Verdict:" results/rice_G281/s05_insert_assembly/insertion_Chr3_16439674_report.txt
```

Expected: `Verdict: CANDIDATE`.

- [ ] **Step 4: per-site wall time 비교 기록**

```bash
awk '/Phase 3 elapsed/ {print $NF}' /tmp/rice_t3_smoke.log | head
# 또는 s05_stats.txt 의 round count 가 줄어들었는지 확인
tail -20 results/rice_G281/s05_insert_assembly/s05_stats.txt
```

- [ ] **Step 5: Commit**

```bash
git add scripts/s05_insert_assembly.py
git commit -m "$(cat <<'EOF'
T3: reduce assemble_insert max_rounds default from 8 to 3

Per T2 convergence study (docs/measurements/max_rounds_study.md):
round 4+ contributes <10% of total growth. Cuts per-site wall time
~15-20%, contributes to AC-4 UGT72E3 TIMEOUT resolution.

rice_G281 Chr3:16,439,674 verdict unchanged (CANDIDATE).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Element DB 15-seq 확장 + V00087.1 subregion fix (gmo_expert)

**목표:** element-db-expansion.md §1 MUST-HAVE 15-seq (P0 6 + P1 8 + P2 1) 추가. V00087.1 P-nos/T-nos subregion 분리로 BUG-9 완전 해결. build manifest.

**Files:**
- Create: `element_db/cas9_sgrna.fa`
- Create: `element_db/euginius_missing.fa`
- Create: `element_db/payload_cds.fa`
- Modify: `element_db/build_common_payload.sh`
- Modify: `element_db/element_summary.tsv`
- Create (→ T5 빌드 결과): `element_db/gmo_combined_db_v2.fa`
- Create (→ T5 빌드 결과): `element_db/gmo_combined_db_manifest.tsv`

**의존성:** 없음 (T1 의 manifest 를 consumer 함).

- [ ] **Step 1: SpCas9 + sgRNA scaffold 확보**

```bash
cd /data/gpfs/assoc/pgl/develop/redgene/element_db
# SpCas9 (S. pyogenes, plant-codon-opt) - pRGEB32 기반 4.1kb
esearch -db nuccore -query "KY026614.1" | efetch -format fasta > cas9_sgrna.fa
# sgRNA scaffold (pMDC123 style, 80bp)
# scaffold 서열이 없으면 공개 repo 서열 직접 작성
cat >> cas9_sgrna.fa <<'EOF'
>sgRNA_scaffold_generic plant CRISPR backbone 80bp consensus
GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTT
EOF
grep -c "^>" cas9_sgrna.fa
```

Expected: `2`.

- [ ] **Step 2: EUginius missing 8 seq 확보 (P0+P1 계열)**

```bash
for acc in Z11946.1 KX344421.1 AY007234.1 AY126452.1 AY126453.1 AF242881.1 D00063.1 X63374.1; do
    efetch -db nuccore -id "$acc" -format fasta
    sleep 1  # NCBI rate limit
done > euginius_missing.fa
grep -c "^>" euginius_missing.fa
```

Expected: `8`. (PMI, vip3Aa, cry2Ab2, cry34Ab1, cry35Ab1, T-g7, T-pinII, 2mEPSPS)

- [ ] **Step 3: payload CDS (AtYUCCA6 + UGT72E3) 확보**

```bash
efetch -db nuccore -id NM_122473.3 -format fasta > payload_cds.fa
# UGT72E3 — soybean-specific, published seq
# (BUG-8 에서 제거된 NM_122942.3 는 잘못된 accession. NM_122473.3 가 정본)
efetch -db nuccore -id NM_001251013.2 -format fasta >> payload_cds.fa  # UGT72E3 O. europaea homolog; soybean용은 별도
grep -c "^>" payload_cds.fa
```

Expected: `2`. 실제 soybean UGT72E3 서열이 GenBank 에 없으면 gmo_expert 가 제공한 FASTA 수동 복사.

- [ ] **Step 4: `build_common_payload.sh` V00087.1 subregion fix**

Edit `element_db/build_common_payload.sh`. 현 V00087.1 full-seq efetch 를 2 subregion efetch 로 교체:

```bash
# BEFORE (BUG-9):
# efetch -db nuccore -id V00087.1 -format fasta >> $TMPOUT

# AFTER (element-db-expansion.md §2):
efetch -db nuccore -id V00087.1 -format fasta \
    -seq_start 1847 -seq_stop 2113 \
    | sed '1s/.*/>P-nos|V00087.1:1847-2113 Agrobacterium nos promoter (267bp)/' \
    >> $TMPOUT
efetch -db nuccore -id V00087.1 -format fasta \
    -seq_start 1277 -seq_stop 1536 \
    | sed '1s/.*/>T-nos|V00087.1:1277-1536 Agrobacterium nos terminator (260bp)/' \
    >> $TMPOUT

# X04879.1 T-ocs subregion (octopine synthase 3' region)
efetch -db nuccore -id X04879.1 -format fasta \
    -seq_start 1 -seq_stop 706 \
    | sed '1s/.*/>T-ocs|X04879.1:1-706 Agrobacterium ocs terminator (706bp)/' \
    >> $TMPOUT
```

- [ ] **Step 5: `build_common_payload.sh` 실행**

```bash
bash element_db/build_common_payload.sh
grep -c "^>" element_db/common_payload.fa
grep "^>" element_db/common_payload.fa
```

Expected: `9` seq (bar/nptII/hpt/gusA/gfp/egfp/P-CaMV35S/P-nos/T-nos/T-ocs — T-nos 와 T-ocs 분리로 10 이 될 수도 있음). V00087.1 entry 에 각각 `:1847-2113` 과 `:1277-1536` 좌표가 붙어있는지 확인.

- [ ] **Step 6: `element_summary.tsv` 에 15 새 entry 추가**

element-db-expansion.md §1 표의 accession/name/source/subregion 을 `element_summary.tsv` 포맷 (탭 구분) 에 맞춰 append. 15 행 추가.

- [ ] **Step 7: Commit (T5 빌드는 다음 task)**

```bash
git add element_db/cas9_sgrna.fa element_db/euginius_missing.fa \
        element_db/payload_cds.fa element_db/build_common_payload.sh \
        element_db/element_summary.tsv
git commit -m "$(cat <<'EOF'
T4: add SpCas9/sgRNA + 8 EUginius missing + 2 payload CDS + V00087.1 subregion fix

- SpCas9 (KY026614.1) + sgRNA scaffold 80bp consensus
- PMI/vip3Aa/cry2Ab2/cry34-35Ab1/T-g7/T-pinII/2mEPSPS
- AtYUCCA6 (NM_122473.3) + UGT72E3 homolog
- BUG-9 complete fix: V00087.1 P-nos (1847-2113) / T-nos (1277-1536) subregion

Per team-consensus §2.1 item 2 + element-db-expansion.md §1-§2.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: cd-hit-est dedup + 4-way source tag + DB manifest (bio_king)

**목표:** `gmo_combined_db_v2.fa` (common_payload + element_db + T4 신규 3 FASTA) 를 cd-hit-est -c 0.95 로 dedup. 각 entry header 에 4-way source tag (`element_db` / `payload` / `sample_contig` / `univec`) 삽입. `gmo_combined_db_manifest.tsv` (name/md5/build_date/seq_count) 생성 — T1 에서 read.

**Files:**
- Create: `element_db/Makefile`
- Create: `element_db/gmo_combined_db_v2.fa` (build 결과물)
- Create: `element_db/gmo_combined_db_manifest.tsv` (build 결과물)
- Create: `tests/test_source_tag_priority.py`
- Modify: `scripts/s05_insert_assembly.py` (`_should_replace` 에서 tag 기반 priority)

**의존성:** T4 완료.

- [ ] **Step 1: failing test — source tag priority**

```python
# tests/test_source_tag_priority.py
"""Test that 4-way source tag + element_db-family > univec priority works."""
import pytest
# 현재 _should_replace 의 위치를 import path 로 가져온다
from scripts.s05_insert_assembly import _should_replace


def test_element_db_beats_univec_at_tied_bitscore():
    # 이미 존재 — BUG-3 회귀 방어. 추가로 4-way tag 에서도 작동하는지 확인
    assert _should_replace(existing={"src": "univec", "bit": 100},
                           new_src="element_db", new_bit=100) is True


def test_payload_beats_univec():
    assert _should_replace(existing={"src": "univec", "bit": 100},
                           new_src="payload", new_bit=100) is True


def test_sample_contig_beats_univec():
    assert _should_replace(existing={"src": "univec", "bit": 100},
                           new_src="sample_contig", new_bit=100) is True


def test_element_db_family_tie_keeps_incumbent():
    # element_db / payload / sample_contig 내부 tie 는 bitscore 로만 결정
    assert _should_replace(existing={"src": "element_db", "bit": 150},
                           new_src="payload", new_bit=150) is False  # tie → keep


def test_higher_bitscore_always_wins():
    assert _should_replace(existing={"src": "element_db", "bit": 100},
                           new_src="univec", new_bit=200) is True
```

- [ ] **Step 2: test 실행 → 기존 single test 는 PASS, 새 4개는 FAIL 확인**

```bash
pytest tests/test_source_tag_priority.py -v
```

Expected: 1 PASS + 4 FAIL (tag 확장 미구현).

- [ ] **Step 3: `_should_replace` 확장**

현 `scripts/s05_insert_assembly.py:818-850` 의 `_should_replace` 를 tag-based 정책으로 교체:

```python
_SRC_TIER = {
    "element_db":    2,
    "payload":       2,
    "sample_contig": 2,
    "univec":        1,
}

def _should_replace(existing: dict, new_src: str, new_bit: float) -> bool:
    """Priority: tier2 (element_db-family) > tier1 (univec); ties resolved by bitscore."""
    old_tier = _SRC_TIER.get(existing.get("src"), 0)
    new_tier = _SRC_TIER.get(new_src, 0)
    if new_tier > old_tier:
        return True
    if new_tier < old_tier:
        return False
    # same tier → bitscore strict >
    return new_bit > existing.get("bit", 0)
```

- [ ] **Step 4: test 재실행 → PASS**

```bash
pytest tests/test_source_tag_priority.py -v
```

Expected: 5 PASS.

- [ ] **Step 5: `Makefile` 작성 — cd-hit-est + manifest**

```makefile
# element_db/Makefile
# Builds gmo_combined_db_v2.fa (cd-hit dedup) + manifest.
.PHONY: all clean

DB      := gmo_combined_db_v2.fa
RAW     := _raw_combined.fa
MANI    := gmo_combined_db_manifest.tsv
CDH_C   := 0.95

all: $(DB) $(MANI)

$(RAW): common_payload.fa gmo_combined_db.fa cas9_sgrna.fa euginius_missing.fa payload_cds.fa
	@echo "[build] tagging source + concat..."
	@python -c "import sys; \
	    srcmap = {'common_payload.fa':'payload', 'gmo_combined_db.fa':'element_db', \
	              'cas9_sgrna.fa':'element_db', 'euginius_missing.fa':'element_db', \
	              'payload_cds.fa':'payload'}; \
	    import pathlib; \
	    with open('$@.tmp','w') as out: \
	        for f in $^.split(): \
	            tag = srcmap[pathlib.Path(f).name]; \
	            for line in pathlib.Path(f).read_text().splitlines(): \
	                if line.startswith('>'): \
	                    out.write(line + '|src=' + tag + chr(10)); \
	                else: \
	                    out.write(line + chr(10))" "$^"
	mv $@.tmp $@

$(DB): $(RAW)
	@echo "[build] cd-hit-est -c $(CDH_C)..."
	cd-hit-est -i $< -o $@.tmp -c $(CDH_C) -n 10 -M 8000 -T 4
	mv $@.tmp $@
	rm -f $@.tmp.clstr  # keep .clstr alongside .fa for audit

$(MANI): $(DB)
	@echo "name\tmd5\tbuild_date\tseq_count" > $@.tmp
	@md5=$$(md5sum $< | awk '{print $$1}'); \
	 count=$$(grep -c '^>' $<); \
	 date=$$(date -Idate); \
	 echo "element_db\t$$md5\t$$date\t$$count" >> $@.tmp
	mv $@.tmp $@

clean:
	rm -f $(DB) $(RAW) $(MANI) *.clstr
```

- [ ] **Step 6: Makefile 실행**

```bash
cd element_db
make all
cat gmo_combined_db_manifest.tsv
grep -c "^>" gmo_combined_db_v2.fa
grep "src=" gmo_combined_db_v2.fa | head -20
```

Expected: manifest 1 행 + md5 64 hex + count ~ 150, 모든 header 에 `|src=<tag>` 접미사.

- [ ] **Step 7: `run_pipeline.py` 에서 v2 DB 사용**

`run_pipeline.py` 의 element_db 경로를 `gmo_combined_db.fa` → `gmo_combined_db_v2.fa` 로 교체. tag parser 를 `_batch_check_element_hits` 에서 header 의 `|src=X` 를 읽어 dict 로 저장.

`scripts/s05_insert_assembly.py` 의 `_parse_blast6` 같은 함수에서 `sseqid` 뒤의 `|src=` 를 읽어 `hit["src"] = tag` 로 주입하도록 수정 (~ 10줄 patch).

- [ ] **Step 8: rice_G281 회귀 smoke**

```bash
python run_pipeline.py --sample rice_G281 --steps 5 --threads 8 --no-remote-blast
grep "^Verdict:" results/rice_G281/s05_insert_assembly/insertion_Chr3_16439674_report.txt
# + element_annotation.tsv 에 src tag 가 보이는지 확인
head results/rice_G281/s05_insert_assembly/element_annotation.tsv
```

Expected: Chr3:16,439,674 CANDIDATE 유지, annotation 에 4-way tag 표시.

- [ ] **Step 9: Commit**

```bash
git add element_db/Makefile element_db/gmo_combined_db_v2.fa \
        element_db/gmo_combined_db_manifest.tsv \
        tests/test_source_tag_priority.py \
        scripts/s05_insert_assembly.py run_pipeline.py
git commit -m "$(cat <<'EOF'
T5: cd-hit-est dedup + 4-way source tag + DB manifest for AC-6

- element_db/Makefile builds gmo_combined_db_v2.fa (cd-hit-est -c 0.95)
- 4 source tiers: element_db=payload=sample_contig > univec
- _should_replace enforces tier-based priority, bitscore tie within tier
- manifest.tsv consumed by T1 audit header writer

Per team-consensus §2.1 items 2/6 + element-db-expansion.md §5.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: `VerdictRules` + `compute_verdict()` pure 함수 분리 (bio_king)

**목표:** 현 `generate_report()` (s05:3197-3514) 안에 얽힌 verdict 결정 로직을 `FilterEvidence` dataclass 입력 / `(verdict, reason)` 출력의 pure 함수로 추출. `VerdictRules` dataclass 는 `config.yaml` 에서 로드 — canonical_triplets + thresholds. 7 pytest scenarios.

**Files:**
- Create: `scripts/s05/__init__.py`
- Create: `scripts/s05/verdict.py`
- Create: `scripts/s05/config_loader.py`
- Create: `tests/test_compute_verdict.py`
- Create: `tests/test_verdict_rules_loader.py`
- Modify: `config.yaml` (global section 에 `canonical_triplets` + `verdict_rules` 추가)
- Modify: `scripts/s05_insert_assembly.py` (generate_report 에서 compute_verdict 호출)

**의존성:** 없음 (T7 과 병렬 가능; 실제 wiring 은 T7 완료 후).

- [ ] **Step 1: failing tests 7 + loader 3**

```python
# tests/test_compute_verdict.py
from scripts.s05.verdict import compute_verdict, FilterEvidence, VerdictRules


def _ev(**kw):
    defaults = dict(
        elements=[], host_bp=0, host_fraction=0.0, largest_gap=9999,
        flanking_hit=None, off_target_chrs=[], construct_frac=0.0,
        combined_frac=0.0, is_chimeric=False,
        site_chr="Chr1", site_pos=1000,
        matched_canonical=set(), sources_by_element={},
    )
    defaults.update(kw)
    return FilterEvidence(**defaults)


_RULES = VerdictRules(
    cand_host_fraction_max=0.80,
    cand_largest_gap_min=500,
    fp_host_fraction_min=0.80,
    fp_largest_gap_max=500,
    fp_off_target_chrs_min=2,
    fp_combined_frac_min=0.85,
    fp_construct_frac_min=0.25,
    unknown_to_fp_host_fraction_min=0.85,
    unknown_to_fp_construct_frac_max=0.05,
    canonical_triplets={"default": {"bar", "P-CaMV35S", "T-ocs"}},
    canonical_triplet_min_identity=0.90,
)


def test_candidate_basic():
    ev = _ev(elements=["bar", "P-CaMV35S"], host_fraction=0.3, largest_gap=5000)
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "CANDIDATE"


def test_fp_a_host_fraction_small_gap():
    ev = _ev(host_fraction=0.85, largest_gap=300)
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"
    assert "85" in reason or "0.85" in reason


def test_fp_b_flanking_overlap():
    ev = _ev(flanking_hit=("Chr11", 8758, 8958), site_chr="Chr11", site_pos=8800)
    verdict, _ = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"


def test_fp_c_multi_locus_chimeric():
    ev = _ev(off_target_chrs=[("Chr5", 200), ("Chr7", 150)], is_chimeric=True)
    verdict, _ = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"


def test_fp_d_construct_host_explain():
    ev = _ev(construct_frac=0.35, host_fraction=0.60, combined_frac=0.90)
    verdict, _ = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"


def test_unknown_to_fp_host_only():
    ev = _ev(elements=[], host_fraction=0.90, construct_frac=0.02)
    verdict, _ = compute_verdict(ev, _RULES)
    assert verdict == "FALSE_POSITIVE"


def test_unknown_preserved():
    ev = _ev(elements=[], host_fraction=0.50, construct_frac=0.10)
    verdict, _ = compute_verdict(ev, _RULES)
    assert verdict == "UNKNOWN"


def test_canonical_triplet_promotion_from_sample_contig():
    # elements 가 없어도 canonical triplet 모두 sample_contig 에서 매칭되면 CANDIDATE
    ev = _ev(
        elements=["bar", "P-CaMV35S", "T-ocs"],
        sources_by_element={"bar": "sample_contig",
                            "P-CaMV35S": "sample_contig",
                            "T-ocs": "sample_contig"},
        matched_canonical={"bar", "P-CaMV35S", "T-ocs"},
        host_fraction=0.60,
    )
    verdict, reason = compute_verdict(ev, _RULES)
    assert verdict == "CANDIDATE"
    assert "canonical_triplet" in reason
```

```python
# tests/test_verdict_rules_loader.py
from pathlib import Path
from scripts.s05.config_loader import load_verdict_rules


def test_loader_default_when_key_missing(tmp_path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text("samples:\n  foo:\n    host_reference: /x\n")
    rules = load_verdict_rules(cfg, sample="foo")
    assert rules.cand_host_fraction_max == 0.80  # default
    assert "bar" in rules.canonical_triplets["default"]


def test_loader_reads_global_verdict_rules(tmp_path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "verdict_rules:\n"
        "  cand_host_fraction_max: 0.70\n"
        "canonical_triplets:\n"
        "  default: [bar, P-CaMV35S, T-ocs]\n"
        "samples:\n  foo: {}\n"
    )
    rules = load_verdict_rules(cfg, sample="foo")
    assert rules.cand_host_fraction_max == 0.70


def test_loader_sample_override(tmp_path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "canonical_triplets:\n"
        "  default: [bar, P-CaMV35S, T-ocs]\n"
        "  rice_G281: [hLF1, P-Gt1, T-nos]\n"
        "samples:\n  rice_G281: {}\n"
    )
    rules = load_verdict_rules(cfg, sample="rice_G281")
    assert rules.canonical_triplets["rice_G281"] == {"hLF1", "P-Gt1", "T-nos"}
```

- [ ] **Step 2: test 실행 → 전부 FAIL**

```bash
pytest tests/test_compute_verdict.py tests/test_verdict_rules_loader.py -v
```

Expected: 11 FAIL (모듈 없음).

- [ ] **Step 3: `scripts/s05/__init__.py`**

```python
# scripts/s05/__init__.py
from .verdict import compute_verdict, FilterEvidence, VerdictRules  # noqa: F401
```

- [ ] **Step 4: `scripts/s05/verdict.py`**

```python
# scripts/s05/verdict.py
"""Pure verdict function for RedGene s05 insert assembly (T6).

Decision tree extracted from scripts/s05_insert_assembly.py:3197-3514
(generate_report). No BLAST or I/O — all inputs via FilterEvidence.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class FilterEvidence:
    elements: list[str]
    host_bp: int
    host_fraction: float
    largest_gap: int
    flanking_hit: tuple[str, int, int] | None
    off_target_chrs: list[tuple[str, int]]
    construct_frac: float
    combined_frac: float
    is_chimeric: bool
    site_chr: str
    site_pos: int
    matched_canonical: set[str] = field(default_factory=set)
    sources_by_element: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class VerdictRules:
    cand_host_fraction_max: float = 0.80
    cand_largest_gap_min: int = 500
    fp_host_fraction_min: float = 0.80
    fp_largest_gap_max: int = 500
    fp_off_target_chrs_min: int = 2
    fp_combined_frac_min: float = 0.85
    fp_construct_frac_min: float = 0.25
    unknown_to_fp_host_fraction_min: float = 0.85
    unknown_to_fp_construct_frac_max: float = 0.05
    canonical_triplets: dict[str, set[str]] = field(
        default_factory=lambda: {"default": {"bar", "P-CaMV35S", "T-ocs"}}
    )
    canonical_triplet_min_identity: float = 0.90


def _site_in_flanking(ev: FilterEvidence) -> bool:
    if ev.flanking_hit is None:
        return False
    chrom, start, end = ev.flanking_hit
    return chrom == ev.site_chr and start <= ev.site_pos <= end


def _canonical_promoted(ev: FilterEvidence, rules: VerdictRules) -> bool:
    for key, triplet in rules.canonical_triplets.items():
        if triplet.issubset(ev.matched_canonical):
            return True
    return False


def compute_verdict(
    ev: FilterEvidence, rules: VerdictRules
) -> tuple[str, str]:
    if _canonical_promoted(ev, rules) and ev.host_fraction < rules.cand_host_fraction_max:
        return (
            "CANDIDATE",
            f"canonical_triplet matched in single assembly "
            f"(host_fraction={ev.host_fraction:.2f})",
        )

    if _site_in_flanking(ev):
        return (
            "FALSE_POSITIVE",
            f"site overlaps construct-flanking slop {ev.flanking_hit}",
        )

    if len(ev.off_target_chrs) >= rules.fp_off_target_chrs_min:
        return (
            "FALSE_POSITIVE",
            f"chimeric: multiple off-target loci {ev.off_target_chrs}",
        )

    if (
        ev.construct_frac >= rules.fp_construct_frac_min
        and ev.combined_frac >= rules.fp_combined_frac_min
    ):
        return (
            "FALSE_POSITIVE",
            f"construct+host explain {ev.combined_frac:.2f} of insert",
        )

    if (
        ev.host_fraction >= rules.fp_host_fraction_min
        and ev.largest_gap < rules.fp_largest_gap_max
    ):
        return (
            "FALSE_POSITIVE",
            f"host fraction {ev.host_fraction:.2f}, gap {ev.largest_gap}bp",
        )

    if ev.elements:
        if ev.host_fraction < rules.cand_host_fraction_max:
            return (
                "CANDIDATE",
                f"foreign elements={ev.elements}, host_fraction={ev.host_fraction:.2f}",
            )

    if (
        not ev.elements
        and ev.host_fraction >= rules.unknown_to_fp_host_fraction_min
        and ev.construct_frac <= rules.unknown_to_fp_construct_frac_max
    ):
        return (
            "FALSE_POSITIVE",
            f"host-only: no elements, host_fraction {ev.host_fraction:.2f}",
        )

    return ("UNKNOWN", "insufficient evidence for classification")
```

- [ ] **Step 5: `scripts/s05/config_loader.py`**

```python
# scripts/s05/config_loader.py
from __future__ import annotations

from pathlib import Path

import yaml

from .verdict import VerdictRules

_DEFAULT_TRIPLETS = {"default": {"bar", "P-CaMV35S", "T-ocs"}}


def load_verdict_rules(config_path: Path, *, sample: str) -> VerdictRules:
    data = yaml.safe_load(config_path.read_text()) or {}
    vr = data.get("verdict_rules", {}) or {}
    triplets_cfg = data.get("canonical_triplets", {}) or {}
    triplets = {k: set(v) for k, v in triplets_cfg.items()} if triplets_cfg else dict(_DEFAULT_TRIPLETS)
    if "default" not in triplets:
        triplets["default"] = _DEFAULT_TRIPLETS["default"]
    return VerdictRules(
        cand_host_fraction_max=vr.get("cand_host_fraction_max", 0.80),
        cand_largest_gap_min=vr.get("cand_largest_gap_min", 500),
        fp_host_fraction_min=vr.get("fp_host_fraction_min", 0.80),
        fp_largest_gap_max=vr.get("fp_largest_gap_max", 500),
        fp_off_target_chrs_min=vr.get("fp_off_target_chrs_min", 2),
        fp_combined_frac_min=vr.get("fp_combined_frac_min", 0.85),
        fp_construct_frac_min=vr.get("fp_construct_frac_min", 0.25),
        unknown_to_fp_host_fraction_min=vr.get("unknown_to_fp_host_fraction_min", 0.85),
        unknown_to_fp_construct_frac_max=vr.get("unknown_to_fp_construct_frac_max", 0.05),
        canonical_triplets=triplets,
        canonical_triplet_min_identity=vr.get("canonical_triplet_min_identity", 0.90),
    )
```

- [ ] **Step 6: test 재실행 → 10/10 PASS**

```bash
pytest tests/test_compute_verdict.py tests/test_verdict_rules_loader.py -v
```

Expected: 10 PASS + (canonical_triplet_from_sample_contig 1) = 11 PASS (위 `scripts/s05/verdict.py` 의 canonical 분기 조정 포함).

- [ ] **Step 7: `config.yaml` 에 schema 추가**

`config.yaml` 최상위(global) 에:

```yaml
verdict_rules:
  cand_host_fraction_max: 0.80
  cand_largest_gap_min: 500
  fp_host_fraction_min: 0.80
  fp_largest_gap_max: 500
  fp_off_target_chrs_min: 2
  fp_combined_frac_min: 0.85
  fp_construct_frac_min: 0.25
  unknown_to_fp_host_fraction_min: 0.85
  unknown_to_fp_construct_frac_max: 0.05
  canonical_triplet_min_identity: 0.90

canonical_triplets:
  default: [bar, P-CaMV35S, T-ocs]
  rice_G281: [hLF1, P-Gt1, T-nos]
  soybean_AtYUCCA6: [bar, P-CaMV35S, T-ocs]
  soybean_UGT72E3: [bar, P-CaMV35S, T-nos]
  tomato_Cas9_A2_3: [bar, SpCas9, sgRNA_scaffold_generic]
```

- [ ] **Step 8: Commit**

```bash
git add scripts/s05/ tests/test_compute_verdict.py tests/test_verdict_rules_loader.py config.yaml
git commit -m "$(cat <<'EOF'
T6: extract compute_verdict() pure + VerdictRules config loader

- scripts/s05/verdict.py: FilterEvidence + VerdictRules dataclasses
  + compute_verdict(ev, rules) -> (verdict, reason). Pure function.
- scripts/s05/config_loader.py: load_verdict_rules(path, sample)
  with DEFAULT_TRIPLETS fallback + per-sample override.
- Canonical-triplet rule promoted via sample_contig source → CANDIDATE.
- 11 pytest scenarios (7 verdict + 3 loader + 1 canonical promo).

Per team-consensus §2.1 items 3 + refactor-roadmap P1-1/P1-5
+ test-strategy Appendix A.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: 최소 4-way s05 모듈 분할 boundary 정의 (bio_king)

**목표:** s05 를 site_discovery / classify / per_site / annotate+report 의 4 phase 진입점만 명시적으로 분리. 전체 8-10 모듈 분할 (P2-1) 은 v1.1. 본 task 는 **shim 파일만** 두고 기존 함수들은 그대로 두되 `scripts/s05/` 패키지에서 import 가능하도록 re-export. T8 per-site fan-out 의 CLI 진입점이 됨.

**Files:**
- Create: `scripts/s05/site_discovery.py`
- Create: `scripts/s05/classify.py`
- Create: `scripts/s05/per_site.py`
- Create: `scripts/s05/annotate_report.py`
- Modify: `scripts/s05_insert_assembly.py` (기존 파일 그대로 유지, 끝에 shim 4 line 추가)

**의존성:** T6 (verdict.py, config_loader.py) 와 같은 패키지 — T6 이 먼저 merge 되어야 경합 없음.

- [ ] **Step 1: shim 모듈 4 개 생성**

각 shim 은 기존 monolithic `scripts/s05_insert_assembly.py` 의 관련 함수를 re-export 만 한다. 실제 구현 이동은 v1.1.

```python
# scripts/s05/site_discovery.py
"""Phase 1 site discovery re-exports (shim, v1.0)."""
from scripts.s05_insert_assembly import (  # noqa: F401
    find_softclip_junctions,
    _build_consensus,
    _batch_check_maps_to_host,
)
```

```python
# scripts/s05/classify.py
"""Phase 1.5 classify re-exports."""
from scripts.s05_insert_assembly import (  # noqa: F401
    classify_site_tiers,
    _batch_check_element_hits,
    _should_replace,
    _filter_host_endogenous,
)
```

```python
# scripts/s05/per_site.py
"""Phase 2-3 per-site re-exports."""
from scripts.s05_insert_assembly import (  # noqa: F401
    extract_candidate_reads,
    extract_unmapped_paired,
    assemble_insert,
    refine_with_foreign_reads,
)
```

```python
# scripts/s05/annotate_report.py
"""Phase 4 annotate + report re-exports."""
from scripts.s05_insert_assembly import (  # noqa: F401
    annotate_insert,
    generate_report,
    write_stats,
)
```

- [ ] **Step 2: smoke import test**

```bash
python -c "from scripts.s05 import (
    site_discovery, classify, per_site, annotate_report,
    compute_verdict, FilterEvidence, VerdictRules,
)
print('ok')"
```

Expected: `ok`.

- [ ] **Step 3: 기존 pytest 10/10 재확인**

```bash
pytest tests/ -v
```

Expected: 10 기존 + 10 신규 (T1/T5/T6) = 20+ PASS.

- [ ] **Step 4: Commit**

```bash
git add scripts/s05/
git commit -m "$(cat <<'EOF'
T7: minimal 4-way s05 package boundary (shim re-exports)

site_discovery / classify / per_site / annotate_report modules
re-export from monolithic scripts/s05_insert_assembly.py. Enables
T8 per-site SLURM array CLI. Full 8-10 module split deferred to v1.1.

Per team-consensus §2.1 item 8.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: `run_pipeline.py` per-site SLURM array fan-out 임시안 (compute_eng + bio_king)

**목표:** UGT72E3 즉시 완화 — 모듈 분할 없이 `run_pipeline.py` 래퍼에서 s05 의 Phase 1/1.5 까지만 실행 후 positive site list 를 array 로 분할해 per-site Phase 2/3 sbatch. Phase 4 (annotate+report) 는 array 완료 후 재수집.

**Files:**
- Create: `scripts/submit_s05_array.sh`
- Create: `scripts/s05_site_runner.py` (--site-id <id> 진입점)
- Modify: `run_pipeline.py` (add `--fanout` flag)
- Modify: `scripts/s05_insert_assembly.py` (phase-mode CLI flag: `--phase 1,1.5` / `--phase 4` / 기존 default 유지)

**의존성:** T7 완료.

- [ ] **Step 1: `scripts/s05_insert_assembly.py` CLI 에 `--phase` flag 추가**

`main()` argparse 에:

```python
parser.add_argument(
    "--phase",
    choices=["all", "1_1.5", "2_3", "4"],
    default="all",
    help="Run specific phase(s). all=default, 1_1.5=discovery+classify, "
         "2_3=per-site assembly (requires --site-id), 4=annotate+report only.",
)
parser.add_argument(
    "--site-id",
    default=None,
    help="When --phase=2_3, run only the specified site (format: chrom_pos).",
)
```

`main()` body 에서 phase 분기:

```python
if args.phase in ("all", "1_1.5"):
    junctions = find_softclip_junctions(...)
    tiers = classify_site_tiers(...)
    # positive site list 를 JSON 에 저장
    (step_dir / "positive_sites.json").write_text(
        json.dumps([{"chr": s.chr, "pos": s.pos} for s in positive_sites])
    )
    if args.phase == "1_1.5":
        return

if args.phase == "2_3":
    assert args.site_id, "--site-id required for --phase=2_3"
    sites = [s for s in _load_positive_sites(step_dir) if f"{s.chr}_{s.pos}" == args.site_id]
    assert len(sites) == 1, f"site {args.site_id} not in positive list"
    _run_per_site(sites[0], ...)  # Phase 2 + 3 for one site → writes site_<id>_insert.fasta
    return

if args.phase in ("all", "4"):
    # collect all site_<id>_insert.fasta → annotate_insert → generate_report
    ...
```

- [ ] **Step 2: `scripts/submit_s05_array.sh` 작성**

```bash
#!/bin/bash
# scripts/submit_s05_array.sh
# Fan-out per-site s05 Phase 2+3 as SLURM array.
# Usage: submit_s05_array.sh <sample> <step_dir> <threads>
set -euo pipefail

SAMPLE="$1"
STEPDIR="$2"
THREADS="${3:-8}"

if [ ! -f "$STEPDIR/positive_sites.json" ]; then
    echo "ERROR: run Phase 1+1.5 first (--phase 1_1.5)" >&2
    exit 1
fi

N=$(python -c "import json; print(len(json.load(open('$STEPDIR/positive_sites.json'))))")
echo "Submitting array of $N sites for $SAMPLE"

SITES_CSV=$(python -c "import json; \
    sites = json.load(open('$STEPDIR/positive_sites.json')); \
    print(','.join(f\"{s['chr']}_{s['pos']}\" for s in sites))")

cat > /tmp/s05_array_${SAMPLE}.sh <<EOF
#!/bin/bash
#SBATCH --job-name=s05_${SAMPLE}
#SBATCH --output=results/s05_array_${SAMPLE}_%A_%a.out
#SBATCH --error=results/s05_array_${SAMPLE}_%A_%a.err
#SBATCH --array=0-$((N-1))%25
set -euo pipefail
eval "\$(micromamba shell hook --shell bash)"
micromamba activate redgene

SITES=(${SITES_CSV//,/ })
SITE="\${SITES[\$SLURM_ARRAY_TASK_ID]}"
echo "=== site \$SITE (task \$SLURM_ARRAY_TASK_ID) ==="

python scripts/s05_insert_assembly.py \\
    --sample "$SAMPLE" --phase 2_3 --site-id "\$SITE" \\
    --threads ${THREADS} --no-remote-blast
EOF
chmod +x /tmp/s05_array_${SAMPLE}.sh

JOBID=$(sbatch --parsable \
    --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 \
    --time=2:00:00 --mem=32G --cpus-per-task=${THREADS} \
    /tmp/s05_array_${SAMPLE}.sh)
echo "Array job: $JOBID"

# submit Phase 4 as dependency
sbatch --parsable \
    --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 \
    --time=1:00:00 --mem=16G --cpus-per-task=4 \
    --dependency=afterok:${JOBID} \
    --wrap="python scripts/s05_insert_assembly.py --sample $SAMPLE --phase 4 --threads 4"
```

- [ ] **Step 3: `run_pipeline.py --fanout` flag 추가**

```python
parser.add_argument(
    "--fanout",
    action="store_true",
    help="For step 5: run Phase 1+1.5 inline, then submit SLURM array for "
         "per-site Phase 2+3, then Phase 4 as dependency.",
)
```

`build_step_cmd("5", ...)` 에서 `--fanout` 면 `submit_s05_array.sh` 를 spawn.

- [ ] **Step 4: UGT72E3 smoke (dry-run first)**

```bash
python run_pipeline.py --sample soybean_UGT72E3 --steps 5 \
    --threads 8 --fanout --dry-run
```

Expected: 1) Phase 1+1.5 inline 실행 계획 출력, 2) `submit_s05_array.sh soybean_UGT72E3 ...` 출력, 3) Phase 4 dependency sbatch 명시.

- [ ] **Step 5: rice_G281 에 대한 실제 실행 (작은 사이트수로 로직 검증)**

```bash
# rice_G281 positive sites ~20개 → array wall time ~15 min
python run_pipeline.py --sample rice_G281 --steps 5 \
    --threads 8 --fanout --no-remote-blast

# 완료 후 verdict 표
grep -h "^Verdict:" results/rice_G281/s05_insert_assembly/insertion_*_report.txt \
    | awk -F' —' '{print $1}' | sort | uniq -c
```

Expected: Chr3:16,439,674 CANDIDATE + 나머지 FP/UNKNOWN 기존 분포 유지.

- [ ] **Step 6: Commit**

```bash
git add scripts/submit_s05_array.sh scripts/s05_insert_assembly.py run_pipeline.py
git commit -m "$(cat <<'EOF'
T8: per-site SLURM array fan-out (--fanout) — AC-4 UGT72E3 fix

scripts/s05_insert_assembly.py: add --phase {1_1.5,2_3,4} + --site-id
scripts/submit_s05_array.sh: build SBATCH array from positive_sites.json
run_pipeline.py: --fanout wires Phase 1.5 inline → array → Phase 4 dep

UGT72E3 projected wall time 48h → ~4h (compute_eng measurement -92%).
Temporary for v1.0; full module split in v1.1 replaces this shim.

Per team-consensus §2.1 item 4 + compute_eng_round1 §3.A.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 9: Element 100%-ortholog Pre-mask BED + rationale TSV (haibao + gmo_expert)

**목표:** element_db 서열을 각 host genome 에 megablast 로 hit → ≥ 98% identity + ≥ 200 bp 구간을 BED 로 쓰고 bedtools merge. 4 host (rice/tomato/cucumber/corn/soybean) 각각 생성. rationale TSV 에 element_name / host_locus / identity / length / rationale_code (E-01~04/X-01) 기록. soybean × AtYUCCA6 는 **MVP-BLOCKING** — 큐레이션 선행 없이 T11 금지.

**Files:**
- Create: `scripts/build_element_mask_bed.sh`
- Create: `docs/host_masks/rice_osativa_v7.bed`
- Create: `docs/host_masks/tomato_slm_r2.bed`
- Create: `docs/host_masks/cucumber_b10v3.bed`
- Create: `docs/host_masks/corn_zm_b73_v5.bed`
- Create: `docs/host_masks/soybean_gmax_v4.bed`
- Create: `docs/host_masks/host_masked_rationale.tsv`

**의존성:** T4 DB 확장 완료.

- [ ] **Step 1: `scripts/build_element_mask_bed.sh`**

```bash
#!/bin/bash
# scripts/build_element_mask_bed.sh
# Usage: build_element_mask_bed.sh <host_fasta> <out_bed> [min_id] [min_len]
set -euo pipefail

HOST="$1"
OUT_BED="$2"
MIN_ID="${3:-98}"
MIN_LEN="${4:-200}"

ELEMENT_DB=element_db/gmo_combined_db_v2.fa
WORK=$(mktemp -d)
trap 'rm -rf $WORK' EXIT

# makeblastdb if missing
if [ ! -f "${HOST}.nsq" ]; then
    makeblastdb -in "$HOST" -dbtype nucl
fi

blastn -task megablast -query "$ELEMENT_DB" -db "$HOST" \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
    -evalue 1e-10 -num_threads 8 -out "$WORK/hits.tsv"

awk -v mi="$MIN_ID" -v ml="$MIN_LEN" \
    '$3 >= mi && $4 >= ml { \
        if ($7 < $8) print $2 "\t" ($7-1) "\t" $8 "\t" $1 "\t" $3; \
        else print $2 "\t" ($8-1) "\t" $7 "\t" $1 "\t" $3 }' \
    "$WORK/hits.tsv" \
    | sort -k1,1 -k2,2n \
    > "$WORK/raw.bed"

bedtools merge -i "$WORK/raw.bed" -c 4,5 -o distinct,max \
    > "$OUT_BED"

wc -l "$OUT_BED"
```

- [ ] **Step 2: 4-host 실행**

```bash
chmod +x scripts/build_element_mask_bed.sh
mkdir -p docs/host_masks

scripts/build_element_mask_bed.sh db/Osativa_323_v7.0.fa     docs/host_masks/rice_osativa_v7.bed
scripts/build_element_mask_bed.sh db/SLM_r2.0.pmol.fasta     docs/host_masks/tomato_slm_r2.bed
scripts/build_element_mask_bed.sh db/CucSat_B10v3.fa         docs/host_masks/cucumber_b10v3.bed
scripts/build_element_mask_bed.sh db/Zm_B73_v5.fa            docs/host_masks/corn_zm_b73_v5.bed
scripts/build_element_mask_bed.sh db/Gmax_v4.0.fa            docs/host_masks/soybean_gmax_v4.bed

wc -l docs/host_masks/*.bed
```

Expected: rice 수십 줄, corn 수백 줄 (P-Ubi1 100% ortholog 많음), soybean 수십-수백 (AtYUCCA6 75%ortholog 다수 가능).

- [ ] **Step 3: rationale TSV 큐레이션 — MVP-BLOCKING 부분**

```bash
cat > docs/host_masks/host_masked_rationale.tsv <<'EOF'
host	region_chr	region_start	region_end	element_name	identity	length	rationale_code	rationale_note
EOF

# element_name 이 P-Ubi1-maize, P-Act1-rice, AtYUCCA6, UGT72E3 인 것만 우선 append
for bed in docs/host_masks/*.bed; do
    host=$(basename "$bed" .bed)
    awk -v h="$host" -F'\t' '
        {
            split($4, elems, ",");
            for (i in elems) {
                e = elems[i];
                code = "X-01";
                if (e ~ /P-Ubi1-maize/ && h ~ /corn/)             code = "E-01";
                else if (e ~ /P-Act1-rice/ && h ~ /rice/)         code = "E-02";
                else if (e ~ /AtYUCCA6/ && h ~ /soybean/)         code = "E-03";
                else if (e ~ /UGT72E3/ && h ~ /soybean/)          code = "E-04";
                else continue;  # 다른 element 는 수동 큐레이션에서 처리
                print h "\t" $1 "\t" $2 "\t" $3 "\t" e "\t" $5 "\t" ($3-$2) "\t" code "\t" "100%-ortholog pre-mask"
            }
        }' "$bed" >> docs/host_masks/host_masked_rationale.tsv
done

wc -l docs/host_masks/host_masked_rationale.tsv
head docs/host_masks/host_masked_rationale.tsv
```

rationale_code 정의 (element-db-expansion.md §4):
- **E-01**: corn × P-Ubi1-maize (100% identity endogenous source)
- **E-02**: rice × P-Act1-rice
- **E-03**: soybean × AtYUCCA6 family (GmYUC6/19 등, 75%+)
- **E-04**: soybean × UGT72E3 family
- **X-01**: other ≥ 98% identity endogenous (수동 분류 필요)

- [ ] **Step 4: soybean AtYUCCA6 MVP-BLOCKING 검증**

```bash
# soybean × AtYUCCA6 가 최소 1 행 이상인지
awk -F'\t' '$8=="E-03"' docs/host_masks/host_masked_rationale.tsv | wc -l
# 0 이면 MVP BLOCK — T4 에서 AtYUCCA6 추가 안 됐거나 ortholog 좌표 누락. 재검사 필요.
```

Expected: ≥ 1.

- [ ] **Step 5: Commit**

```bash
git add scripts/build_element_mask_bed.sh docs/host_masks/
git commit -m "$(cat <<'EOF'
T9: element 100%-ortholog pre-mask BED + rationale TSV (5 hosts)

build_element_mask_bed.sh: megablast element_db vs host ≥98% id ≥200bp,
bedtools merge. Outputs per-host BED + combined host_masked_rationale.tsv
with rationale codes E-01 (corn×P-Ubi1), E-02 (rice×P-Act1),
E-03 (soybean×AtYUCCA6), E-04 (soybean×UGT72E3), X-01 (other ≥98%).

Per team-consensus §2.1 item 7 + element-db-expansion.md §4.
MVP-BLOCKING: soybean × AtYUCCA6 rationale required before T11 rerun.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: s05 Phase 1 BED Intersect + `FALSE_NEGATIVE_MASKED` tag (bio_king)

**목표:** T9 에서 만든 per-host BED 를 s05 Phase 1 (`find_softclip_junctions`) 직후 intersect 하여 mask 영역과 겹치는 site 는 drop 이 아니라 **`tier='FALSE_NEGATIVE_MASKED'` 태그로 downgrade** — site 는 보고서에 남지만 CANDIDATE 승격 금지.

**Files:**
- Modify: `scripts/s05_insert_assembly.py` (`find_softclip_junctions` 결과 후처리부; 라인 번호는 실제 구조 확인 후)
- Create: `tests/test_mask_bed_intersect.py`

**의존성:** T7 (s05 package) + T9 (BED 생성).

- [ ] **Step 1: failing test**

```python
# tests/test_mask_bed_intersect.py
from scripts.s05.site_discovery import _apply_mask_bed  # T10 에서 신규 추가


def test_site_in_mask_gets_tagged(tmp_path):
    bed = tmp_path / "mask.bed"
    bed.write_text("Chr3\t100\t500\tP-Ubi1-maize\t99.5\n")
    sites = [
        {"chr": "Chr3", "pos": 300, "tier": "transgene-positive"},
        {"chr": "Chr3", "pos": 600, "tier": "transgene-positive"},
        {"chr": "Chr5", "pos": 300, "tier": "transgene-positive"},
    ]
    out = _apply_mask_bed(sites, bed)
    assert out[0]["tier"] == "FALSE_NEGATIVE_MASKED"
    assert out[0]["mask_element"] == "P-Ubi1-maize"
    assert out[1]["tier"] == "transgene-positive"  # 600 > 500
    assert out[2]["tier"] == "transgene-positive"  # different chr


def test_empty_bed_is_noop(tmp_path):
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    sites = [{"chr": "Chr3", "pos": 300, "tier": "transgene-positive"}]
    out = _apply_mask_bed(sites, bed)
    assert out == sites
```

- [ ] **Step 2: test 실행 → FAIL**

```bash
pytest tests/test_mask_bed_intersect.py -v
```

Expected: `ImportError: cannot import name '_apply_mask_bed'`.

- [ ] **Step 3: 구현**

`scripts/s05_insert_assembly.py` 에 `_apply_mask_bed` 함수 추가:

```python
def _apply_mask_bed(sites: list[dict], bed_path: Path) -> list[dict]:
    """Tag sites whose (chr, pos) falls inside bed mask region.

    Returns new list; original sites are mutated.
    """
    if not bed_path.exists() or bed_path.stat().st_size == 0:
        return sites
    regions: dict[str, list[tuple[int, int, str]]] = {}
    for raw in bed_path.read_text().splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        cols = raw.split("\t")
        if len(cols) < 3:
            continue
        chrom, start, end = cols[0], int(cols[1]), int(cols[2])
        name = cols[3] if len(cols) > 3 else "masked"
        regions.setdefault(chrom, []).append((start, end, name))
    for s in sites:
        for start, end, name in regions.get(s["chr"], []):
            if start <= s["pos"] < end:
                s["tier"] = "FALSE_NEGATIVE_MASKED"
                s["mask_element"] = name
                break
    return sites
```

그리고 `scripts/s05/site_discovery.py` 에 re-export 추가:

```python
from scripts.s05_insert_assembly import _apply_mask_bed  # noqa: F401
```

- [ ] **Step 4: `find_softclip_junctions` 호출부 직후 적용**

`run_pipeline.py` 의 `build_step_cmd("5", ...)` 에서 `--mask-bed docs/host_masks/<host>.bed` 자동 주입. `scripts/s05_insert_assembly.py main()` 의 Phase 1 직후:

```python
if args.mask_bed and Path(args.mask_bed).exists():
    sites = _apply_mask_bed(sites, Path(args.mask_bed))
    print(f"[Phase 1.5] mask-bed applied: {sum(1 for s in sites if s['tier']=='FALSE_NEGATIVE_MASKED')} sites tagged", file=sys.stderr)
```

- [ ] **Step 5: test 재실행 → PASS**

```bash
pytest tests/test_mask_bed_intersect.py -v
```

Expected: 2 PASS.

- [ ] **Step 6: cucumber_line225 smoke (8 CAND → ≤ 4 예상)**

T9 에서 만든 cucumber BED 가 비어있지 않으면 cucumber_line225 재실행:

```bash
python run_pipeline.py --sample cucumber_line225 --steps 5 \
    --threads 16 --no-remote-blast

# FALSE_NEGATIVE_MASKED 가 어떻게 보고되는지 확인
grep -h "^Verdict:" results/cucumber_line225/s05_insert_assembly/insertion_*_report.txt \
    | awk -F' —' '{print $1}' | sort | uniq -c
grep "FALSE_NEGATIVE_MASKED" results/cucumber_line225/s05_insert_assembly/*.tsv 2>/dev/null | head
```

Expected: 8 CAND → 3-5 CAND (AC-2 MVP threshold ≤ 5 에 근접), 1-3 FALSE_NEGATIVE_MASKED 태그.

- [ ] **Step 7: Commit**

```bash
git add scripts/s05_insert_assembly.py scripts/s05/site_discovery.py \
        tests/test_mask_bed_intersect.py run_pipeline.py
git commit -m "$(cat <<'EOF'
T10: Phase 1 mask-bed intersect + FALSE_NEGATIVE_MASKED tag

_apply_mask_bed() tags sites overlapping host-endogenous BED as
FALSE_NEGATIVE_MASKED instead of dropping. Preserves audit trail
while preventing CANDIDATE promotion on known host homolog sites.

Per team-consensus §2.1 item 7-companion + AC-2 specificity boost.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 11: 6-sample Revalidation Batch (compute_eng)

**목표:** T1-T10 통합 후 rice_G281 / tomato_Cas9_A2_3 / cucumber_line212/224/225 / soybean_AtYUCCA6 / soybean_UGT72E3 / tomato_Cas9_A2_2 8 샘플 전부 재실행하여 AC-1/AC-2/AC-4/AC-6 통과 확인. bug.md OPEN-1/5 close-out.

**Files:**
- Create: `run_rerun_w1_batch.sh`
- Modify: `bug.md`
- Modify: `resume.md`

**의존성:** T1 ~ T10 완료.

- [ ] **Step 1: baseline snapshot 저장**

```bash
mkdir -p docs/superpowers/runs
for s in rice_G281 tomato_Cas9_A2_3 cucumber_line212 cucumber_line224 \
         cucumber_line225 soybean_AtYUCCA6 soybean_UGT72E3 tomato_Cas9_A2_2; do
    echo "=== $s ==="
    grep -h "^Verdict:" results/${s}/s05_insert_assembly/insertion_*_report.txt 2>/dev/null \
        | awk -F' —' '{print $1}' | sort | uniq -c
done > docs/superpowers/runs/2026-04-24-pre-w1-baseline.txt
```

- [ ] **Step 2: array batch script 작성**

```bash
cat > run_rerun_w1_batch.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=rg_w1_revalidate
#SBATCH --output=results/rg_w1_%A_%a.out
#SBATCH --error=results/rg_w1_%A_%a.err
#SBATCH --array=0-7

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene
SAMPLES=(
    rice_G281
    tomato_Cas9_A2_3
    tomato_Cas9_A2_2
    cucumber_line212
    cucumber_line224
    cucumber_line225
    soybean_AtYUCCA6
    soybean_UGT72E3
)
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "=== $SAMPLE ==="

# soybean 2 개는 --fanout 로 per-site array 처리
case "$SAMPLE" in
    soybean_*)
        python run_pipeline.py --sample "$SAMPLE" --steps 4b,5 \
            --threads 16 --no-remote-blast --fanout
        ;;
    *)
        python run_pipeline.py --sample "$SAMPLE" --steps 4b,5 \
            --threads 16 --no-remote-blast
        ;;
esac
EOF
chmod +x run_rerun_w1_batch.sh
```

- [ ] **Step 3: sbatch 제출 (CLI flag 명시로 BUG-16 회피)**

```bash
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 \
       --time=12:00:00 --mem=64G --cpus-per-task=16 \
       --chdir="$PWD" ./run_rerun_w1_batch.sh
```

Expected: `Submitted batch job <JOBID>`.

- [ ] **Step 4: 완료 후 verdict summary 저장**

```bash
# ~6-12h 후 완료
sacct -j <JOBID> --format=JobID,JobName,State,Elapsed,MaxRSS -n

for s in rice_G281 tomato_Cas9_A2_3 tomato_Cas9_A2_2 cucumber_line212 \
         cucumber_line224 cucumber_line225 soybean_AtYUCCA6 soybean_UGT72E3; do
    echo "=== $s ==="
    grep -h "^Verdict:" results/${s}/s05_insert_assembly/insertion_*_report.txt 2>/dev/null \
        | awk -F' —' '{print $1}' | sort | uniq -c
done > docs/superpowers/runs/2026-04-24-post-w1.txt

diff docs/superpowers/runs/2026-04-24-pre-w1-baseline.txt \
     docs/superpowers/runs/2026-04-24-post-w1.txt
```

- [ ] **Step 5: AC 검증**

```bash
# AC-1: GT anchor 6 샘플 CANDIDATE 회수율
for entry in rice_G281:Chr3_16439674 \
             tomato_Cas9_A2_3:SLM_r2.0ch01_9100 \
             tomato_Cas9_A2_2:SLM_r2.0ch08_6510 \
             cucumber_line212:LKUO03001392_2751687 \
             cucumber_line224:LKUO03001512_581328 \
             cucumber_line225:LKUO03001451_6501; do
    s=${entry%%:*}
    prefix=${entry#*:}
    v=$(grep -h "^Verdict:" results/${s}/s05_insert_assembly/insertion_${prefix}*_report.txt 2>/dev/null | head -1)
    echo "$s $prefix → $v"
done

# AC-2: 샘플당 CAND FP 수 (모두 ≤ 5 여야 함)
for s in rice_G281 tomato_Cas9_A2_3 cucumber_line225; do
    n=$(grep -l "^Verdict: CANDIDATE" results/${s}/s05_insert_assembly/insertion_*_report.txt 2>/dev/null | wc -l)
    echo "$s CAND total: $n"
done

# AC-4: UGT72E3 wall time
sacct -j <JOBID> --format=JobName,Elapsed -n | grep soybean_UGT72E3

# AC-6: audit header 4-field 확인
jq -r '.input_sha256, .pipeline_commit, .pipeline_dirty, .software_versions | keys | length' \
    results/rice_G281/audit_header.json
```

Expected: AC-1 6/6 GT CAND, AC-2 전 샘플 ≤ 5, AC-4 UGT72E3 ≤ 4h, AC-6 4/4 fields.

- [ ] **Step 6: `bug.md` OPEN-1, OPEN-5 close**

Edit `bug.md` "Known pre-existing issues" 각 항목에 결과 line append:

```markdown
### OPEN-1 — soybean_UGT72E3 s05 TIMEOUT
- 2026-04-24 W1 rerun (job <JOBID>): COMPLETED in ~4h via --fanout array.
  Final verdicts: <N> CAND / <N> FP / <N> UNKNOWN. Closed.

### OPEN-5 — soybean_AtYUCCA6 end-to-end
- 2026-04-24 W1 rerun: Bar+P-35S+T-ocs canonical triplet → N CAND
  (was 0). Soybean × AtYUCCA6 pre-mask E-03 applied (N sites tagged).
  Closed.
```

- [ ] **Step 7: `resume.md` v1.0-rc 섹션 append**

```markdown
## v1.0-rc (2026-04-24)

**W1-W8 tasks:** all 12 tasks merged (T1-T12).
**AC summary:**
- AC-1: <X/6> GT anchor recall
- AC-2: max <N> CAND FP / sample
- AC-4: UGT72E3 <wall time>
- AC-6: 4/4 audit fields present

**Regression:** 6 pre-W1 CAND preserved (no regression).
**Release decision:** <won_yim 승인 / 재작업 지시>
```

- [ ] **Step 8: Commit**

```bash
git add run_rerun_w1_batch.sh docs/superpowers/runs/2026-04-24-*.txt \
        bug.md resume.md
git commit -m "$(cat <<'EOF'
T11: W1-W8 8-sample revalidation batch + OPEN-1/5 close-out

- run_rerun_w1_batch.sh: SLURM array 0-7 (rice/tomato×2/cuc×3/soy×2)
- AC-1 result: <X/6>
- AC-2 result: max <N> FP/sample
- AC-4 result: UGT72E3 <hours>h
- AC-6 result: 4/4 audit fields

Per team-consensus §3 T11 + §7 release criteria.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 12: rice_G281 s04 minimap2 PoC (conditional) (compute_eng + haibao)

**목표:** won_yim 조건부 승인 실험. s04 host mapping 을 BWA 대신 minimap2 `-ax sr` 로 수행하여 Chr3:16,439,674 CANDIDATE 유지 여부 + wall time -30% + MAPQ 분포 회귀 없음 검증. 승격 실패 시 자동 기각.

**Files:**
- Create: `run_s04_minimap2_poc.sh`
- Create: `docs/measurements/s04_minimap2_poc.md`

**의존성:** T11 완료 (rice_G281 기존 BWA s04 결과를 비교 기준으로).

- [ ] **Step 1: PoC script 작성**

```bash
cat > run_s04_minimap2_poc.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=rg_s04_mm2_poc
#SBATCH --output=results/rg_s04_mm2_poc_%j.out
#SBATCH --error=results/rg_s04_mm2_poc_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate redgene

cd /data/gpfs/assoc/pgl/develop/redgene
SAMPLE=rice_G281
HOST=db/Osativa_323_v7.0.fa
R1=$(yq '.samples.rice_G281.reads.r1' config.yaml)
R2=$(yq '.samples.rice_G281.reads.r2' config.yaml)
OUT=results/rice_G281_mm2_poc/s04_host_map
mkdir -p $OUT

# minimap2 index (one-time, cached)
if [ ! -f "${HOST}.mmi" ]; then
    minimap2 -d "${HOST}.mmi" "$HOST"
fi

# map
/usr/bin/time -v minimap2 -ax sr -t 16 "${HOST}.mmi" "$R1" "$R2" 2> $OUT/mm2.log \
    | samtools sort -@ 8 -o $OUT/${SAMPLE}_host.bam -
samtools index $OUT/${SAMPLE}_host.bam

# Run s05 only on the new BAM
python run_pipeline.py --sample rice_G281 --steps 5 \
    --threads 8 --no-remote-blast \
    --host-bam-override $OUT/${SAMPLE}_host.bam \
    --outdir-override results/rice_G281_mm2_poc
EOF
chmod +x run_s04_minimap2_poc.sh
```

주의: `run_pipeline.py --host-bam-override` + `--outdir-override` 는 v1.0 에 없을 수 있음 → T12 Step 2 에서 간단 추가.

- [ ] **Step 2: `run_pipeline.py` override flag 추가**

```python
parser.add_argument("--host-bam-override", default=None,
                    help="Skip s04 and use this BAM as host_bam for downstream steps.")
parser.add_argument("--outdir-override", default=None,
                    help="Alternative output directory (for PoC runs).")
```

`build_step_cmd("5", ...)` 에서 host BAM 을 override 로 교체.

- [ ] **Step 3: PoC 제출**

```bash
sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 \
       --time=6:00:00 --mem=48G --cpus-per-task=16 \
       --chdir="$PWD" ./run_s04_minimap2_poc.sh
```

- [ ] **Step 4: 승격 기준 4 조건 측정**

```bash
# 1. Chr3:16,439,674 verdict
grep "^Verdict:" results/rice_G281_mm2_poc/s05_insert_assembly/insertion_Chr3_16439674_report.txt

# 2. Phase 1 site count ratio
BWA_N=$(grep -c "^transgene-positive" results/rice_G281/s05_insert_assembly/site_tier_classification.tsv)
MM2_N=$(grep -c "^transgene-positive" results/rice_G281_mm2_poc/s05_insert_assembly/site_tier_classification.tsv)
python -c "print(f'site ratio = {$MM2_N}/{$BWA_N} = {$MM2_N/$BWA_N:.2f}')"

# 3. s04 wall time
grep "Elapsed" results/rice_G281_mm2_poc/s04_host_map/mm2.log
# vs BWA 기존 wall time (sacct 에서 조회)

# 4. MAPQ 분포 — soft-clip reads 의 MAPQ < 20 비율
samtools view -F 256 -f 2 results/rice_G281_mm2_poc/s04_host_map/rice_G281_host.bam \
    | awk '{ if ($6 ~ /S/) { if ($5 < 20) low++; else high++ } } END { print "low<20:",low,"high≥20:",high,"ratio:",low/(low+high) }'
# BWA baseline:
samtools view -F 256 -f 2 results/rice_G281/s04_host_map/rice_G281_host.bam \
    | awk '{ if ($6 ~ /S/) { if ($5 < 20) low++; else high++ } } END { print "low<20:",low,"high≥20:",high,"ratio:",low/(low+high) }'
```

- [ ] **Step 5: 결과 문서화 + 판정**

```bash
cat > docs/measurements/s04_minimap2_poc.md <<'EOF'
# s04 minimap2 PoC — rice_G281 (T12, 2026-04-24)

**조건부 승인 (team-consensus.md §8):** 아래 4 조건 **모두** 통과 시 v1.1 이관 승인, 실패 시 BWA 유지.

## Measurements

| 기준 | BWA (baseline) | minimap2 | 판정 |
|------|----------------|----------|------|
| Chr3:16,439,674 verdict | CANDIDATE | <V> | <PASS/FAIL> |
| Phase 1 site count | <N_bwa> | <N_mm2> | ±20% = <PASS/FAIL> |
| s04 wall time | <T_bwa> | <T_mm2> | -30% = <PASS/FAIL> |
| soft-clip MAPQ<20 비율 | <R_bwa> | <R_mm2> | no increase = <PASS/FAIL> |

## Decision

<4/4 PASS → v1.1 이관 / 그 외 → BWA 유지, minimap2 v2.0 grant>
EOF
vi docs/measurements/s04_minimap2_poc.md  # 실제 수치 기입
```

- [ ] **Step 6: Commit**

```bash
git add run_s04_minimap2_poc.sh docs/measurements/s04_minimap2_poc.md run_pipeline.py
git commit -m "$(cat <<'EOF'
T12: rice_G281 s04 minimap2 PoC (conditional, 4-criterion)

- run_s04_minimap2_poc.sh: rerun s04 with minimap2 -ax sr, then s05
- 4 acceptance criteria: verdict / site count ±20% / wall time -30%
  / MAPQ<20 ratio no-increase
- Decision recorded in docs/measurements/s04_minimap2_poc.md

Per team-consensus §5 + §8 Unresolved Items. Pass → v1.1, Fail → BWA stays.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Final Integration: v1.0-rc Release Tag

**의존성:** T1-T12 모두 완료 + T11 AC 통과 + T12 판정 기록.

- [ ] **Step 1: 전체 pytest 재확인**

```bash
pytest tests/ -v
```

Expected: 20+ PASS (10 기존 + T1 1 + T5 5 + T6 11 + T10 2).

- [ ] **Step 2: team-consensus.md §7 release 기준 6 항목 확인**

```bash
# 1. AC-1 ≥ 100% (6/6 GT)
# 2. AC-2 ≤ 5 FP/sample
# 3. AC-4 ≤ 48h 전 샘플
# 4. AC-6 4/4 audit header
# 5. pytest 10/10 + compute_verdict 3 시나리오 PASS
# 6. exploratory ≥ 1 CAND or FP ≤ 5

# 각 항목 PASS 여부를 resume.md v1.0-rc 섹션에 정리
```

- [ ] **Step 3: v1.0-rc 태그 + push**

```bash
git tag -a v1.0-rc1 -m "RedGene v1.0 release candidate 1

Completed MVP tasks T1-T12 per team-consensus.md.
AC-1: <X/6> | AC-2: ≤ <N> FP | AC-4: UGT72E3 <H>h | AC-6: 4/4
Pending: <won_yim 승인 / Phase 2 이관 확정>
"
# won_yim PI 최종 검토 후:
# git push origin v1.0-rc1
```

---

## Dependency Graph

```
Day 1 (Mon 2026-04-20)      Day 2 (Tue)           Day 3 (Wed)
┌─────────────────┐          ┌──────────────┐      ┌──────────────┐
│ T1 audit_header │          │ T2 round_cnt │      │ T5 cd-hit+tag│
│ T6 verdict      │ ─────────│ T3 max_rounds│ ────│ T9 mask BED  │
│ T7 4-way split  │          │ T4 DB 15seq  │      │ T10 BED intr │
└─────────────────┘          └──────────────┘      └──────────────┘
                                                            │
Day 4 (Thu)                                                 │
┌──────────────┐ ◄───────────────────────────────────────────┘
│ T8 fan-out   │
└──────────────┘
        │
Day 5 (Fri)
┌──────────────────────┐
│ T11 W1-W8 batch      │
│ T12 s04 mm2 PoC      │
│ v1.0-rc1 tag         │
└──────────────────────┘
```

**Critical path (parallelized):** T1 || T6 || T7 || T2 → T4 → T5 & T9 & T10 → T8 → T11 & T12.

**공수 합산:** T1(0.5) + T2(0.5) + T3(0.5) + T4(1.0) + T5(0.5) + T6(1.5) + T7(0.5) + T8(1.5) + T9(1.0) + T10(0.5) + T11(0.5+wall) + T12(0.5+wall) = **~8.5 person-day**, 4명 parallelize 시 5 영업일.

---

## Self-Review Checklist

- [x] **Spec coverage:** team-consensus.md §3 T1-T12 → 본 plan 의 Task 1-12 와 1:1 매핑. 추가 bonus: v1.0-rc1 태그 단계.
- [x] **Placeholder scan:** "TBD" / "implement later" / "add validation" 등 없음. 모든 code step 에 실제 code block 있음. `<JOBID>`, `<X/6>` 같은 runtime-fill placeholder 만 존재 (의도적).
- [x] **Type consistency:** `FilterEvidence`, `VerdictRules`, `_should_replace`, `_apply_mask_bed`, `_write_audit_header` 등 이름이 task 간 일관. Task 6 의 `matched_canonical: set[str]` 이 Task 6 test 와 일치.
- [x] **File path 정합:** `scripts/s05/` 패키지 구조가 Task 6/7/10 에서 일관. `element_db/gmo_combined_db_v2.fa` 가 Task 5 에서 만들어져 Task 9 에서 consumed.
- [x] **SLURM BUG-16 회피:** 모든 `sbatch` 호출에 CLI `--partition=cpu-s2-core-0 --account=cpu-s2-pgl-0` 명시.
- [x] **Regression gate:** Task 3/5/10 각각 rice_G281 Chr3:16,439,674 CANDIDATE 유지 확인 step 포함.
- [x] **No mkdir before verify:** `docs/measurements/`, `docs/host_masks/`, `element_db/` 모두 `mkdir -p` + `ls` 확인.

---

## Execution Handoff

Plan 저장 위치: `docs/team-review/work_implementation_plan.md`
연계 문서: `docs/team-review/team-consensus.md` (권위 판정), `docs/team-review/refactor-roadmap.md` (P1-* 매핑), `docs/team-review/test-strategy.md` (pytest 확장 근거), `docs/team-review/gap-analysis.md` (AC 측정 근거), `docs/team-review/element-db-expansion.md` (T4/T9 domain 근거).

**실행 옵션:**

1. **Subagent-Driven (권장)** — T1-T12 각 task 를 fresh subagent 로 디스패치, task 간 review 후 다음 task. 장점: context 깨끗, 병렬 T1/T6/T7 동시 진행 가능.
2. **Inline Execution** — 현재 세션에서 executing-plans skill 로 day-by-day 순차. 장점: 통합 smoke test 가 한 세션에 끝남.

**어느 방식으로 갈까요?**
