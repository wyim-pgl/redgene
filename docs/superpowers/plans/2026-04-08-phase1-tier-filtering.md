# Phase 1 Tiered Filtering — Eliminate SV Assembly Waste

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce Phase 3 assembly targets from ~9,600 to <50 by adding a tiered element-screening filter between Phase 1 junction detection and Phase 3 assembly, without losing any ground truth sites.

**Architecture:** Phase 1 detects soft-clip junctions (unchanged logic), then a new Tier classifier (BLAST clip sequences against element_db at two sensitivity levels) categorizes sites as Tier A (element-positive, assemble), Tier B (ambiguous, assemble with flag), or Tier C (element-negative, skip assembly). A ground truth baseline is recorded first to prevent regressions.

**Tech Stack:** Python 3.11, blastn, minimap2, pysam, pathlib

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `scripts/s09_targeted_assembly.py` | Modify | Add tier classification between Phase 1 and Phase 3 |
| `scripts/s09_ground_truth_baseline.py` | Create | Standalone script to extract ground truth clip stats |
| `ground_truth_baseline.tsv` | Create (output) | Safety net — ground truth clip metrics |

All changes are in `scripts/s09_targeted_assembly.py` except the one-shot baseline script.

---

### Task 0: Record Ground Truth Baseline

**Files:**
- Create: `scripts/s09_ground_truth_baseline.py`
- Output: `ground_truth_baseline.tsv`

This is the safety net. We extract clip-level metrics for every known ground truth junction site so we can verify no parameter change drops them.

- [ ] **Step 0.1: Create the baseline extraction script**

```python
#!/usr/bin/env python3
"""Extract soft-clip metrics for known ground truth junction sites.

Scans the host BAM at each ground truth position for soft-clipped reads,
measures clip lengths, read counts, element DB identity, and host mapping.
Outputs ground_truth_baseline.tsv.
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import pysam

# Known ground truth sites (chr, position, sample, description)
GROUND_TRUTH = [
    # Rice G281
    ("Chr3", 16439674, "rice_G281", "T-DNA insertion, 2-copy head-to-head"),
    ("Chr3", 16439710, "rice_G281", "T-DNA insertion, second junction"),
    # Tomato Cas9 A2_3
    ("SLM_r2.0ch08", 65107378, "tomato_Cas9_A2_3", "T-DNA RB junction"),
    # Cucumber line 224
    ("LKUO03001512.1", 580628, "cucumber_line224", "T-DNA LB junction"),
    ("LKUO03001512.1", 581332, "cucumber_line224", "T-DNA RB junction"),
    # Cucumber line 212
    ("LKUO03001392.1", 2751693, "cucumber_line212", "T-DNA junction"),
    # Cucumber line 225
    ("LKUO03002166.1", 547987, "cucumber_line225", "T-DNA junction"),
    # Corn ND207
    ("NC_050098.1", 181367276, "corn_ND207", "T-DNA LB junction"),
    # Soybean UGT72E3
    ("NC_038254.2", 51882860, "soybean_UGT72E3", "T-DNA junction in Glyma.18g226800"),
    # Soybean AtYUCCA6 Site II
    ("NC_038255.2", 49789752, "soybean_AtYUCCA6", "T-DNA Site II junction"),
]


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def extract_clips_at_position(
    bam_path: Path, chrom: str, pos: int, window: int = 100, min_clip: int = 15,
) -> dict:
    """Extract soft-clip stats at a ground truth position."""
    right_clips: list[str] = []  # clips extending past mapped region (5' seed)
    left_clips: list[str] = []   # clips before mapped region (3' seed)
    total_reads = 0

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(chrom, max(0, pos - window), pos + window):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if not read.cigartuples or not read.query_sequence:
                    continue
                total_reads += 1
                cigar = read.cigartuples
                seq = read.query_sequence

                # Right soft-clip (5' junction seed)
                if cigar[-1][0] == 4 and cigar[-1][1] >= min_clip:
                    if read.reference_end and abs(read.reference_end - pos) <= window:
                        clip_len = cigar[-1][1]
                        right_clips.append(seq[-clip_len:])

                # Left soft-clip (3' junction seed)
                if cigar[0][0] == 4 and cigar[0][1] >= min_clip:
                    if read.reference_start and abs(read.reference_start - pos) <= window:
                        clip_len = cigar[0][1]
                        left_clips.append(seq[:clip_len])
    except ValueError:
        # Chromosome not in BAM
        return {"found": False}

    if not right_clips and not left_clips:
        return {"found": False, "total_reads": total_reads}

    # Build consensus for longest clips
    right_consensus = max(right_clips, key=len) if right_clips else ""
    left_consensus = max(left_clips, key=len) if left_clips else ""

    return {
        "found": True,
        "total_reads": total_reads,
        "right_clip_count": len(right_clips),
        "left_clip_count": len(left_clips),
        "right_clip_max_len": max((len(c) for c in right_clips), default=0),
        "left_clip_max_len": max((len(c) for c in left_clips), default=0),
        "right_clip_median_len": sorted(len(c) for c in right_clips)[len(right_clips)//2] if right_clips else 0,
        "left_clip_median_len": sorted(len(c) for c in left_clips)[len(left_clips)//2] if left_clips else 0,
        "right_consensus": right_consensus,
        "left_consensus": left_consensus,
    }


def check_element_hit(seq: str, element_db: Path, evalue: str = "1e-5") -> dict:
    """BLAST a clip sequence against element_db, return best hit."""
    if not seq or len(seq) < 15:
        return {"hit": False}
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
        fh.write(f">query\n{seq}\n")
        query_fa = Path(fh.name)
    out_path = query_fa.with_suffix(".blast")
    subprocess.run(
        ["blastn", "-query", str(query_fa), "-subject", str(element_db),
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
         "-evalue", evalue, "-max_target_seqs", "1",
         "-out", str(out_path)],
        stderr=subprocess.DEVNULL,
    )
    result = {"hit": False}
    if out_path.exists() and out_path.stat().st_size > 0:
        cols = open(out_path).readline().strip().split("\t")
        if len(cols) >= 6:
            result = {
                "hit": True,
                "element": cols[1],
                "identity": float(cols[2]),
                "aln_length": int(cols[3]),
                "evalue": cols[4],
                "bitscore": float(cols[5]),
            }
    query_fa.unlink(missing_ok=True)
    out_path.unlink(missing_ok=True)
    return result


def check_host_mapping(seq: str, host_ref: Path) -> bool:
    """Check if clip sequence maps to host genome."""
    if not seq or len(seq) < 15:
        return False
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
        fh.write(f">query\n{seq}\n")
        query_fa = Path(fh.name)
    paf_path = query_fa.with_suffix(".paf")
    subprocess.run(
        ["minimap2", "-c", "--secondary=no", str(host_ref), str(query_fa)],
        stdout=open(paf_path, "w"), stderr=subprocess.DEVNULL,
    )
    maps = False
    if paf_path.exists():
        for line in open(paf_path):
            cols = line.strip().split("\t")
            if len(cols) >= 12:
                match_bp = int(cols[9])
                if match_bp >= 30:
                    maps = True
                    break
    query_fa.unlink(missing_ok=True)
    paf_path.unlink(missing_ok=True)
    return maps


def main():
    parser = argparse.ArgumentParser(description="Extract ground truth clip baselines")
    parser.add_argument("--host-bam", required=True)
    parser.add_argument("--host-ref", required=True)
    parser.add_argument("--element-db", required=True)
    parser.add_argument("--sample", required=True,
                        help="Sample name (must match GROUND_TRUTH entries)")
    parser.add_argument("--output", default="ground_truth_baseline.tsv")
    parser.add_argument("--min-clip", type=int, default=10,
                        help="Min clip for scanning (use low value to catch everything)")
    args = parser.parse_args()

    host_bam = Path(args.host_bam)
    host_ref = Path(args.host_ref)
    element_db = Path(args.element_db)

    # Filter ground truth to this sample
    sample_gt = [gt for gt in GROUND_TRUTH if gt[2] == args.sample]
    if not sample_gt:
        print(f"No ground truth entries for sample '{args.sample}'", file=sys.stderr)
        print(f"Available: {sorted(set(g[2] for g in GROUND_TRUTH))}", file=sys.stderr)
        sys.exit(1)

    # Check if output exists — append mode
    write_header = not Path(args.output).exists()
    with open(args.output, "a") as fout:
        if write_header:
            fout.write("sample\tchrom\tposition\tdescription\t"
                       "right_clip_count\tleft_clip_count\t"
                       "right_clip_max_len\tleft_clip_max_len\t"
                       "right_clip_median_len\tleft_clip_median_len\t"
                       "element_hit_5p\telement_id_5p\telement_identity_5p\telement_aln_len_5p\t"
                       "element_hit_3p\telement_id_3p\telement_identity_3p\telement_aln_len_3p\t"
                       "maps_to_host_5p\tmaps_to_host_3p\n")

        for chrom, pos, sample, desc in sample_gt:
            print(f"Checking {sample} {chrom}:{pos} ({desc})...", file=sys.stderr)
            clips = extract_clips_at_position(
                host_bam, chrom, pos, min_clip=args.min_clip)

            if not clips.get("found"):
                print(f"  WARNING: No clips found at {chrom}:{pos}", file=sys.stderr)
                fout.write(f"{sample}\t{chrom}\t{pos}\t{desc}\t"
                           f"0\t0\t0\t0\t0\t0\t"
                           f"False\t-\t0\t0\tFalse\t-\t0\t0\t"
                           f"False\tFalse\n")
                continue

            # Element hits
            elem_5p = check_element_hit(clips.get("right_consensus", ""), element_db)
            elem_3p = check_element_hit(clips.get("left_consensus", ""), element_db)

            # Host mapping
            host_5p = check_host_mapping(clips.get("right_consensus", ""), host_ref)
            host_3p = check_host_mapping(clips.get("left_consensus", ""), host_ref)

            print(f"  5p: {clips['right_clip_count']} clips, "
                  f"max {clips['right_clip_max_len']}bp, "
                  f"element={'YES' if elem_5p['hit'] else 'no'}, "
                  f"host={'YES' if host_5p else 'no'}", file=sys.stderr)
            print(f"  3p: {clips['left_clip_count']} clips, "
                  f"max {clips['left_clip_max_len']}bp, "
                  f"element={'YES' if elem_3p['hit'] else 'no'}, "
                  f"host={'YES' if host_3p else 'no'}", file=sys.stderr)

            fout.write(
                f"{sample}\t{chrom}\t{pos}\t{desc}\t"
                f"{clips['right_clip_count']}\t{clips['left_clip_count']}\t"
                f"{clips['right_clip_max_len']}\t{clips['left_clip_max_len']}\t"
                f"{clips['right_clip_median_len']}\t{clips['left_clip_median_len']}\t"
                f"{elem_5p.get('hit', False)}\t{elem_5p.get('element', '-')}\t"
                f"{elem_5p.get('identity', 0)}\t{elem_5p.get('aln_length', 0)}\t"
                f"{elem_3p.get('hit', False)}\t{elem_3p.get('element', '-')}\t"
                f"{elem_3p.get('identity', 0)}\t{elem_3p.get('aln_length', 0)}\t"
                f"{host_5p}\t{host_3p}\n"
            )

    print(f"Baseline written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
```

- [ ] **Step 0.2: Run baseline for rice G281**

```bash
eval "$(micromamba shell hook --shell bash)" && micromamba activate redgene
python scripts/s09_ground_truth_baseline.py \
  --host-bam results/rice_G281/s07_host_map/rice_G281_host.bam \
  --host-ref db/Osativa_323_v7.0.fa \
  --element-db db/G281_construct.fa \
  --sample rice_G281 \
  --output ground_truth_baseline.tsv
```

Expected: TSV with 2 rows (Chr3:16439674 and Chr3:16439710). Both should show:
- `right_clip_count >= 3` and `left_clip_count >= 3`
- `right_clip_max_len >= 30`
- `element_hit_5p = True` or `element_hit_3p = True`
- `maps_to_host_5p = False` and `maps_to_host_3p = False`

**STOP CHECK:** If either ground truth position shows 0 clips or no element hit, investigate before proceeding. This is the baseline safety net.

- [ ] **Step 0.3: Commit baseline**

```bash
git add scripts/s09_ground_truth_baseline.py ground_truth_baseline.tsv
git commit -m "Add ground truth baseline extraction script and rice G281 baseline"
```

---

### Task 1: Tighten Phase 1 Parameters

**Files:**
- Modify: `scripts/s09_targeted_assembly.py` (function `find_softclip_junctions`, line 301)

We test min_clip and min_support_reads thresholds against the baseline.

- [ ] **Step 1.1: Parameter sweep — min_clip**

Write and run a quick parameter sweep. Do NOT modify s09 yet — just analyze.

```bash
# Extract current Phase 1 stats from the running job's log
# Count candidates at different clip thresholds from the log
grep "VALIDATED" results/rice_G281/s09_targeted_assembly/slurm_5624223.err | \
  awk -F'[=,)]' '{print $2, $4}' | \
  awk '{gsub(/bp/,""); print}' | \
  sort -n | head -20
```

Alternatively, write an inline analysis of the ground_truth_baseline.tsv:

```python
#!/usr/bin/env python3
"""Analyze ground truth baselines to find safe parameter thresholds."""
import csv
import sys

with open("ground_truth_baseline.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)

print("Ground truth clip metrics:")
print(f"{'Sample':<25} {'Chr:Pos':<25} {'5p_clips':>8} {'5p_max':>6} "
      f"{'3p_clips':>8} {'3p_max':>6} {'elem_5p':>8} {'elem_3p':>8}")
print("-" * 110)

min_clip_len = 999
min_read_count = 999
for r in rows:
    r5 = int(r["right_clip_count"])
    l5 = int(r["left_clip_count"])
    r5_max = int(r["right_clip_max_len"])
    l5_max = int(r["left_clip_max_len"])
    min_clip_len = min(min_clip_len, r5_max if r5 > 0 else 999, l5_max if l5 > 0 else 999)
    min_read_count = min(min_read_count, r5 if r5 > 0 else 999, l5 if l5 > 0 else 999)
    print(f"{r['sample']:<25} {r['chrom']}:{r['position']:<15} "
          f"{r5:>8} {r5_max:>6} {l5:>8} {l5_max:>6} "
          f"{r['element_hit_5p']:>8} {r['element_hit_3p']:>8}")

print(f"\nMinimum clip length across all ground truth: {min_clip_len}bp")
print(f"Minimum read count across all ground truth: {min_read_count}")
print(f"\nSafe min_clip threshold: {min(min_clip_len - 5, 25)}bp")
print(f"Safe min_support threshold: {max(min_read_count - 1, 2)}")
```

Run: `python /dev/stdin < above_script`

- [ ] **Step 1.2: Build parameter sweep table**

Using the ground truth analysis, construct this table:

| min_clip | min_support | total_candidates_rice | ground_truth_retained | ground_truth_lost |
|----------|-------------|----------------------|----------------------|-------------------|

Test values: min_clip in {20, 25, 30}, min_support in {3, 5, 8}.

To get candidate counts without running the full pipeline, add a `--dry-run-phase1` flag temporarily:

```python
# In find_softclip_junctions, after line 509 (clip-difference filter):
# Add temporary logging:
log(f"  [SWEEP] min_clip={min_clip}, total_paired={len(paired_sites)}, "
    f"cond2_passed={len(cond2_passed)}")
```

Pick the most aggressive thresholds that retain ALL ground truth. If ground truth clips are all >40bp with >5 reads, we can safely set min_clip=25 and min_support=3.

- [ ] **Step 1.3: Apply chosen parameters**

Update `find_softclip_junctions` default parameters based on Step 1.2 results.

In `scripts/s09_targeted_assembly.py`, line 306-307:

```python
# BEFORE:
def find_softclip_junctions(
    host_bam: Path,
    host_ref: Path,
    element_db: Path | None,
    workdir: Path,
    min_clip: int = 20,
    cluster_window: int = 50,

# AFTER (example — use actual values from Step 1.2):
def find_softclip_junctions(
    host_bam: Path,
    host_ref: Path,
    element_db: Path | None,
    workdir: Path,
    min_clip: int = 25,      # was 20; ground truth min is Xbp
    cluster_window: int = 50,
```

Also update `MIN_CLUSTER_DEPTH` at line 359 if the sweep supports it:

```python
# BEFORE:
MIN_CLUSTER_DEPTH = 3

# AFTER (if ground truth min support >= 5):
MIN_CLUSTER_DEPTH = 3  # Keep conservative — tighter filter is in Tier classification
```

- [ ] **Step 1.4: Commit**

```bash
git add scripts/s09_targeted_assembly.py
git commit -m "Tighten Phase 1 min_clip based on ground truth baseline (was 20bp)"
```

---

### Task 2: Add Tiered Element Screening (Tier A/B/C)

**Files:**
- Modify: `scripts/s09_targeted_assembly.py` — add `classify_site_tiers()` function and `site_tier_classification.tsv` output

This is the core change. Add a new function between `find_softclip_junctions()` and the Phase 3 assembly loop.

- [ ] **Step 2.1: Add the tier classification function**

Insert after the `find_softclip_junctions` function (after line ~599), before the legacy junction parser:

```python
# ---------------------------------------------------------------------------
# Tier classification: element-based filtering before assembly
# ---------------------------------------------------------------------------

@dataclass
class TierResult:
    """Tier classification for one insertion site."""
    site_id: str
    tier: str           # "A", "B", "C"
    reason: str         # human-readable reason
    element_hit_std: str    # standard BLAST best hit (or "")
    element_identity_std: float
    element_aln_len_std: int
    element_hit_sens: str   # sensitive BLAST best hit (or "")
    element_identity_sens: float
    element_aln_len_sens: int
    host_map_5p: bool
    host_map_3p: bool
    annotation: str     # "T-DNA", "endogenous_SV", "TE_insertion", ""


def classify_site_tiers(
    sites: list[InsertionSite],
    element_db: Path,
    host_ref: Path,
    workdir: Path,
) -> tuple[list[InsertionSite], list[InsertionSite], list[TierResult]]:
    """Classify validated sites into Tier A/B/C based on element DB screening.

    Tier A: Standard BLAST hit >=80% identity over >=20bp → assemble
    Tier B: Only sensitive BLAST hit, or low-complexity clip, or short clip → assemble (flagged)
    Tier C: No element hit at all → skip assembly

    Returns (tier_a_b_sites, tier_c_sites, all_tier_results).
    """
    if not sites:
        return [], [], []

    tier_dir = workdir / "_tier_classification"
    tier_dir.mkdir(parents=True, exist_ok=True)

    # Collect all clip sequences for batch BLAST
    clip_seqs: dict[str, str] = {}
    for site in sites:
        if site.seed_5p and len(site.seed_5p) >= 15:
            clip_seqs[f"{site.site_id}_5p"] = site.seed_5p
        if site.seed_3p and len(site.seed_3p) >= 15:
            clip_seqs[f"{site.site_id}_3p"] = site.seed_3p

    if not clip_seqs:
        log("  Tier classification: no clip sequences to classify")
        results = [TierResult(s.site_id, "B", "no_clips", "", 0, 0, "", 0, 0,
                              False, False, "") for s in sites]
        return sites, [], results

    # Write all clips to one FASTA
    clip_fa = tier_dir / "all_clips.fa"
    with open(clip_fa, "w") as fh:
        for name, seq in clip_seqs.items():
            fh.write(f">{name}\n{seq}\n")

    # ---- Standard BLAST (evalue 1e-5) ----
    std_out = tier_dir / "std_blast.tsv"
    subprocess.run(
        ["blastn", "-query", str(clip_fa), "-subject", str(element_db),
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
         "-evalue", "1e-5", "-max_target_seqs", "5",
         "-out", str(std_out)],
        stderr=subprocess.DEVNULL,
    )

    # ---- Sensitive BLAST (evalue 1e-3, word_size 7) ----
    sens_out = tier_dir / "sens_blast.tsv"
    subprocess.run(
        ["blastn", "-query", str(clip_fa), "-subject", str(element_db),
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
         "-evalue", "1e-3", "-word_size", "7", "-max_target_seqs", "5",
         "-out", str(sens_out)],
        stderr=subprocess.DEVNULL,
    )

    # Parse standard hits: best hit per query
    std_hits: dict[str, dict] = {}
    if std_out.exists():
        with open(std_out) as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 6:
                    continue
                qname = cols[0]
                identity = float(cols[2])
                aln_len = int(cols[3])
                if qname not in std_hits or float(cols[5]) > std_hits[qname].get("bitscore", 0):
                    std_hits[qname] = {
                        "element": cols[1], "identity": identity,
                        "aln_length": aln_len, "bitscore": float(cols[5]),
                    }

    # Parse sensitive hits: best hit per query
    sens_hits: dict[str, dict] = {}
    if sens_out.exists():
        with open(sens_out) as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 6:
                    continue
                qname = cols[0]
                if qname not in sens_hits or float(cols[5]) > sens_hits[qname].get("bitscore", 0):
                    sens_hits[qname] = {
                        "element": cols[1], "identity": float(cols[2]),
                        "aln_length": int(cols[3]), "bitscore": float(cols[5]),
                    }

    # ---- Host mapping for Tier C clips ----
    # (done later, only for Tier C sites)

    # ---- Classify each site ----
    tier_results: list[TierResult] = []
    tier_ab: list[InsertionSite] = []
    tier_c: list[InsertionSite] = []

    for site in sites:
        key_5p = f"{site.site_id}_5p"
        key_3p = f"{site.site_id}_3p"

        # Standard BLAST: check for >=80% identity over >=20bp
        std_5p = std_hits.get(key_5p, {})
        std_3p = std_hits.get(key_3p, {})
        has_std_hit = (
            (std_5p.get("identity", 0) >= 80 and std_5p.get("aln_length", 0) >= 20) or
            (std_3p.get("identity", 0) >= 80 and std_3p.get("aln_length", 0) >= 20)
        )

        # Sensitive BLAST: any hit at all
        sens_5p = sens_hits.get(key_5p, {})
        sens_3p = sens_hits.get(key_3p, {})
        has_sens_hit = bool(sens_5p) or bool(sens_3p)

        # Low-complexity check
        def is_low_complexity(seq: str) -> bool:
            if not seq or len(seq) < 10:
                return False
            freq = max(seq.upper().count(b) for b in "ACGT") / len(seq)
            return freq > 0.60

        clip_low_complexity = (
            is_low_complexity(site.seed_5p) or is_low_complexity(site.seed_3p)
        )

        # Short clip check
        clip_short = (
            (site.seed_5p and len(site.seed_5p) < 30) or
            (site.seed_3p and len(site.seed_3p) < 30)
        )

        # Assign tier
        if has_std_hit:
            tier = "A"
            reason = "standard_blast_hit"
            best = std_5p if std_5p.get("bitscore", 0) >= std_3p.get("bitscore", 0) else std_3p
        elif has_sens_hit:
            tier = "B"
            reason = "sensitive_blast_only"
            best = sens_5p if sens_5p.get("bitscore", 0) >= sens_3p.get("bitscore", 0) else sens_3p
        elif clip_low_complexity:
            tier = "B"
            reason = "low_complexity_clip"
            best = {}
        elif clip_short:
            tier = "B"
            reason = "short_clip"
            best = {}
        else:
            tier = "C"
            reason = "no_element_hit"
            best = {}

        # Best standard hit info
        best_std = std_5p if std_5p.get("bitscore", 0) >= std_3p.get("bitscore", 0) else std_3p
        best_sens = sens_5p if sens_5p.get("bitscore", 0) >= sens_3p.get("bitscore", 0) else sens_3p

        tr = TierResult(
            site_id=site.site_id,
            tier=tier,
            reason=reason,
            element_hit_std=best_std.get("element", ""),
            element_identity_std=best_std.get("identity", 0),
            element_aln_len_std=best_std.get("aln_length", 0),
            element_hit_sens=best_sens.get("element", ""),
            element_identity_sens=best_sens.get("identity", 0),
            element_aln_len_sens=best_sens.get("aln_length", 0),
            host_map_5p=False,  # filled later for Tier C
            host_map_3p=False,
            annotation="",
        )
        tier_results.append(tr)

        if tier in ("A", "B"):
            tier_ab.append(site)
        else:
            tier_c.append(site)

    # ---- Host mapping for Tier C sites (batch) ----
    if tier_c:
        tierc_seqs: dict[str, str] = {}
        for site in tier_c:
            if site.seed_5p:
                tierc_seqs[f"{site.site_id}_5p"] = site.seed_5p
            if site.seed_3p:
                tierc_seqs[f"{site.site_id}_3p"] = site.seed_3p

        tierc_host_set = _batch_check_maps_to_host(tierc_seqs, host_ref, tier_dir)

        for tr in tier_results:
            if tr.tier == "C":
                maps_5p = f"{tr.site_id}_5p" in tierc_host_set
                maps_3p = f"{tr.site_id}_3p" in tierc_host_set
                tr.host_map_5p = maps_5p
                tr.host_map_3p = maps_3p
                if maps_5p or maps_3p:
                    tr.annotation = "endogenous_SV"
                else:
                    tr.annotation = "unknown_non_element"

    # Cleanup
    shutil.rmtree(tier_dir, ignore_errors=True)

    # Log summary
    n_a = sum(1 for t in tier_results if t.tier == "A")
    n_b = sum(1 for t in tier_results if t.tier == "B")
    n_c = sum(1 for t in tier_results if t.tier == "C")
    log(f"  Tier classification: {n_a} Tier A, {n_b} Tier B, {n_c} Tier C")
    for tr in tier_results:
        if tr.tier in ("A", "B"):
            log(f"    {tr.site_id}: TIER {tr.tier} ({tr.reason}) "
                f"— {tr.element_hit_std or tr.element_hit_sens}")

    return tier_ab, tier_c, tier_results
```

- [ ] **Step 2.2: Add tier TSV output function**

Insert right after `classify_site_tiers`:

```python
def write_tier_classification(
    tier_results: list[TierResult],
    output_path: Path,
) -> None:
    """Write site_tier_classification.tsv."""
    with open(output_path, "w") as fh:
        fh.write("site_id\ttier\treason\t"
                 "element_hit_std\telement_identity_std\telement_aln_len_std\t"
                 "element_hit_sens\telement_identity_sens\telement_aln_len_sens\t"
                 "host_map_5p\thost_map_3p\tannotation\n")
        for tr in tier_results:
            fh.write(f"{tr.site_id}\t{tr.tier}\t{tr.reason}\t"
                     f"{tr.element_hit_std}\t{tr.element_identity_std}\t{tr.element_aln_len_std}\t"
                     f"{tr.element_hit_sens}\t{tr.element_identity_sens}\t{tr.element_aln_len_sens}\t"
                     f"{tr.host_map_5p}\t{tr.host_map_3p}\t{tr.annotation}\n")
    log(f"  Tier classification written: {output_path}")
```

- [ ] **Step 2.3: Commit tier classification code**

```bash
git add scripts/s09_targeted_assembly.py
git commit -m "Add tiered element screening (Tier A/B/C) before Phase 3 assembly"
```

---

### Task 3: Wire Tier Classification into Main Pipeline

**Files:**
- Modify: `scripts/s09_targeted_assembly.py` — `main()` function (line ~2489)

- [ ] **Step 3.1: Insert tier classification between Phase 1 and Phase 3**

Replace the current block at lines 2489-2516 (Phase 1 → direct Phase 3) with tiered flow:

```python
    # ---- Phase 1: Soft-clip junction detection ----
    sites = find_softclip_junctions(
        host_bam, host_ref, element_db, step_dir,
        min_clip=args.min_clip,
    )

    # Fallback to step 6 junctions if no sites found
    if not sites and args.junctions:
        junctions_path = Path(args.junctions)
        if junctions_path.exists() and junctions_path.stat().st_size > 0:
            log("No soft-clip sites found, falling back to step 6 junctions...")
            legacy_juncs = parse_legacy_junctions(junctions_path)
            if legacy_juncs:
                sites = legacy_junctions_to_sites(legacy_juncs, host_bam)
                log(f"  Converted {len(legacy_juncs)} legacy junctions → "
                    f"{len(sites)} sites")

    if not sites:
        log("No insertion sites found. Nothing to assemble.")
        write_stats(step_dir / "s09_stats.txt", args.sample_name, 0, 0, [])
        return

    # ---- Phase 1.5: Tier classification ----
    log(f"\n{'=' * 60}")
    log(f"=== Phase 1.5: Tier Classification ({len(sites)} sites) ===")
    tier_ab_sites, tier_c_sites, tier_results = classify_site_tiers(
        sites, element_db, host_ref, step_dir,
    )

    # Write tier classification TSV (always — includes all sites)
    write_tier_classification(tier_results, step_dir / "site_tier_classification.tsv")

    # Summary
    log(f"  Assembly targets: {len(tier_ab_sites)} (Tier A+B)")
    log(f"  Skipped (Tier C): {len(tier_c_sites)}")
    n_sv = sum(1 for t in tier_results if t.annotation == "endogenous_SV")
    log(f"  Of Tier C: {n_sv} endogenous SVs, "
        f"{len(tier_c_sites) - n_sv} unknown non-element")

    # Use tier_ab_sites for Phase 3 (replacing `sites`)
    sites = tier_ab_sites

    if not sites:
        log("No Tier A/B sites to assemble.")
        write_stats(step_dir / "s09_stats.txt", args.sample_name,
                    len(tier_ab_sites) + len(tier_c_sites), 0, [])
        return

    log(f"\nProcessing {len(sites)} insertion site(s):")
    for site in sites:
        log(f"  {site.site_id}: {site.host_chr}:{site.pos_5p}"
            f"{f'-{site.pos_3p}' if site.pos_3p else ''} "
            f"({site.confidence}, 5p={len(site.seed_5p)}bp, "
            f"3p={len(site.seed_3p)}bp)")
```

- [ ] **Step 3.2: Remove the old MAX_NON_ELEMENT cap**

In `find_softclip_junctions`, remove lines 589-597 (the `MAX_NON_ELEMENT = 5` cap):

```python
# DELETE this block:
    # Limit: all element-hit sites + top N non-element sites
    MAX_NON_ELEMENT = 5
    element_sites = [s for s in validated_sites if s.has_element_hits]
    non_element_sites = [s for s in validated_sites if not s.has_element_hits]
    validated_sites = element_sites + non_element_sites[:MAX_NON_ELEMENT]
    if len(non_element_sites) > MAX_NON_ELEMENT:
        log(f"  Keeping {len(element_sites)} element-hit + "
            f"{MAX_NON_ELEMENT} top non-element sites "
            f"(dropped {len(non_element_sites) - MAX_NON_ELEMENT})")
```

The tier classifier now handles this filtering properly — all validated sites pass to the tier classifier, which decides what to assemble.

- [ ] **Step 3.3: Add Tier C sites to final output**

After the Phase 3/4 assembly loop, before `write_stats`, add Tier C site info to the stats:

```python
    # ---- Write Tier C skipped sites to stats ----
    stats_results = all_results.copy()
    for site in tier_c_sites:
        tier_info = next((t for t in tier_results if t.site_id == site.site_id), None)
        stats_results.append({
            "site_id": site.site_id,
            "insert_length": 0,
            "remaining_ns": 0,
            "rounds": 0,
            "status": f"skipped_{tier_info.annotation if tier_info else 'tier_c'}",
        })
```

- [ ] **Step 3.4: Syntax check**

```bash
python -c "import ast; ast.parse(open('scripts/s09_targeted_assembly.py').read()); print('OK')"
```

Expected: `OK`

- [ ] **Step 3.5: Commit**

```bash
git add scripts/s09_targeted_assembly.py
git commit -m "Wire tier classification into main pipeline, skip Tier C assembly"
```

---

### Task 4: Validation on Rice Data

**Files:**
- No code changes — run pipeline and verify results

- [ ] **Step 4.1: Run step 9 on rice G281 with new code**

```bash
sbatch run_rice_s09.sh
```

Wait for completion. Check log output.

- [ ] **Step 4.2: Verify ground truth is Tier A**

```bash
grep "Chr3:16439674\|insertion_18908\|TIER A" results/rice_G281/s09_targeted_assembly/slurm_*.err
cat results/rice_G281/s09_targeted_assembly/site_tier_classification.tsv | \
  awk -F'\t' '$2 == "A" || NR == 1'
```

**STOP CHECK:** `insertion_18908` (Chr3:16,439,674) MUST be Tier A. If it's Tier B or C, the element DB search thresholds need loosening.

- [ ] **Step 4.3: Compare before/after**

```bash
# Count assembly targets
echo "=== Before ==="
echo "Total validated: ~9600 + 12 element-hit"
echo "Assembly targets: 12 + 5 = 17 (MAX_NON_ELEMENT cap)"

echo "=== After ==="
grep -c "TIER A" results/rice_G281/s09_targeted_assembly/site_tier_classification.tsv
grep -c "TIER B" results/rice_G281/s09_targeted_assembly/site_tier_classification.tsv
grep -c "TIER C" results/rice_G281/s09_targeted_assembly/site_tier_classification.tsv
```

Expected: Tier A ~12, Tier B ~5-20, Tier C ~9,500+. Assembly targets drop from ~9,600 to <50.

- [ ] **Step 4.4: Verify runtime improvement**

Compare wall-clock time from SLURM output:
```bash
sacct -j <NEW_JOB_ID> --format=JobID,Elapsed,MaxRSS
sacct -j 5624223 --format=JobID,Elapsed,MaxRSS
```

- [ ] **Step 4.5: Final commit with validation results**

```bash
git add results/rice_G281/s09_targeted_assembly/site_tier_classification.tsv
git commit -m "Validate tier filtering on rice G281: Tier A/B/C classification verified"
```

---

## Constraints Checklist

- [x] TIER C sites are NEVER deleted — written to `site_tier_classification.tsv` and `s09_stats.txt`
- [x] Phase 3 assembly logic is untouched — only input site list changes
- [x] Ground truth baseline recorded BEFORE any parameter changes
- [x] Every parameter choice justified by `ground_truth_baseline.tsv`
- [x] If ground truth site is classified as Tier C → STOP and report error
