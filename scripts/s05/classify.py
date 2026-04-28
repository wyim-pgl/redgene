"""Phase 1.5 — Classify (Issue #4, Session 4).

Extracted verbatim from ``scripts/s05_insert_assembly.py``. Owns the
tier-based classification logic that turns Phase 1 soft-clip junctions
+ Phase 4 annotation hits into TierResult / InsertionSite verdicts.

Stage 1 in the s05 DAG (depends only on ``primitives``). Imports
``_parse_src_tag`` from ``.primitives`` (Session 3 migrated it there
to break a runtime cycle with ``annotation``).

Replaces the v1.0 shim that re-exported from the monolith.

Module-level constants mirrored from monolith (HOST_ENDO_*, CLIP_HOST_*).
After Session 4, the monolith imports these back via
``from s05.classify import ...`` for backward compatibility.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from .primitives import (
    log,
    _parse_src_tag,
    InsertionSite,
    TierResult,
)


# ---------------------------------------------------------------------------
# Host-endogenous exclusion thresholds (mirrored from monolith for classify)
# ---------------------------------------------------------------------------
# DB-level: BLAST transgene_db entries vs host genome to find host-derived
# elements (e.g., P-Act1-rice in a rice genome). Two tiers handle cultivar
# drift — P-Act1 hits Nipponbare at ~78%/39% due to Xiushui vs Nipponbare.
HOST_ENDO_T1_PIDENT = 90.0   # Tier 1 (clean host): min % identity
HOST_ENDO_T1_QCOVS = 50.0    # Tier 1: min query coverage %
HOST_ENDO_T2_PIDENT = 75.0   # Tier 2 (divergent host): min % identity
HOST_ENDO_T2_QCOVS = 30.0    # Tier 2: min query coverage %

# Per-clip host verification: stricter than DB-level because individual clips
# are short (~20-50 bp), so even moderate-identity host hits are significant.
# DB-level Tier 2 uses 75% because full-length elements (500-3000 bp) may
# diverge across cultivars; per-clip uses 95% because a 30 bp clip hitting
# host at <95% is too noisy to trust as evidence of host origin.
CLIP_HOST_PIDENT = 95.0       # Per-clip host check: min % identity
CLIP_HOST_MIN_LEN = 30        # Per-clip host check: min alignment length


def _batch_check_element_hits(
    seqs: dict[str, str],
    element_db: Path,
    workdir: Path,
    extra_dbs: list[Path] | None = None,
) -> dict[str, list[str]]:
    """Batch check which sequences hit element DB(s) using blastn.

    When ``extra_dbs`` is provided (e.g., the always-on ``common_payload.fa``
    and/or per-sample SPAdes contigs from s04b), each query is BLASTed
    against every DB in turn and hits are merged. This lets shared
    transgene payloads and sample-specific assemblies contribute to
    element annotation without rebuilding the shared DB.
    """
    if not seqs:
        return {}

    query_fa = workdir / "_clip_element_batch.fa"
    with open(query_fa, "w") as fh:
        for name, seq in seqs.items():
            if len(seq) >= 20:
                fh.write(f">{name}\n{seq}\n")

    hits: dict[str, list[str]] = defaultdict(list)
    dbs = [element_db]
    for edb in (extra_dbs or []):
        if edb is not None and edb.exists() and edb.stat().st_size > 0:
            dbs.append(edb)

    for i, db in enumerate(dbs):
        blast_out = workdir / f"_clip_blast_{i}_{db.stem}.tsv"
        subprocess.run(
            ["blastn", "-query", str(query_fa), "-subject", str(db),
             "-outfmt", "6 qseqid sseqid pident length",
             "-evalue", "1e-3", "-max_target_seqs", "5",
             "-out", str(blast_out)],
            stderr=subprocess.DEVNULL,
        )
        if blast_out.exists():
            with open(blast_out) as fh:
                for line in fh:
                    cols = line.strip().split("\t")
                    if len(cols) >= 4 and int(cols[3]) >= 20:
                        hits[cols[0]].append(cols[1])
            blast_out.unlink(missing_ok=True)

    query_fa.unlink(missing_ok=True)
    return dict(hits)


def _filter_host_endogenous(
    transgene_db: Path,
    host_ref: Path,
    tier_dir: Path,
    workdir: Path,
    threads: int,
) -> tuple[Path, set[str]]:
    """BLAST transgene_db vs host genome to remove host-derived entries.

    Returns (blast_db, exclude_ids): blast_db is the filtered DB path
    (or original if nothing excluded), exclude_ids are the removed entries.
    """
    blast_db = transgene_db
    exclude_ids: set[str] = set()
    exclude_details: dict[str, str] = {}

    # Ensure host BLAST DB exists
    host_blast_db_exists = (
        host_ref.with_suffix(".fa.ndb").exists()
        or Path(str(host_ref) + ".ndb").exists()
    )
    if not host_blast_db_exists:
        log("  Host BLAST DB not found — creating with makeblastdb...")
        result = subprocess.run(
            ["makeblastdb", "-in", str(host_ref), "-dbtype", "nucl"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        if result.returncode == 0:
            host_blast_db_exists = True
        else:
            log("  WARNING: makeblastdb failed, skipping host-endogenous exclusion")

    if not host_blast_db_exists:
        log("  WARNING: Host BLAST DB unavailable, skipping host-endogenous exclusion")
        return blast_db, exclude_ids

    host_hits_out = tier_dir / "host_endogenous_hits.tsv"
    log("  BLAST transgene_db vs host genome to find host-endogenous entries...")
    result = subprocess.run(
        ["blastn", "-task", "blastn",
         "-query", str(transgene_db), "-db", str(host_ref),
         "-outfmt", "6 qseqid qlen sseqid pident length qcovs",
         "-evalue", "1e-10", "-max_target_seqs", "1",
         "-num_threads", str(threads),
         "-out", str(host_hits_out)],
        stderr=subprocess.DEVNULL,
    )
    if result.returncode != 0:
        log(f"  WARNING: host-endogenous BLAST failed (rc={result.returncode})")
        return blast_db, exclude_ids

    if host_hits_out.exists():
        with open(host_hits_out) as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 6:
                    continue
                qseqid = cols[0]
                pident = float(cols[3])
                qcovs = float(cols[5])
                tier1 = pident >= HOST_ENDO_T1_PIDENT and qcovs >= HOST_ENDO_T1_QCOVS
                tier2 = pident >= HOST_ENDO_T2_PIDENT and qcovs >= HOST_ENDO_T2_QCOVS
                if tier1 or tier2:
                    tier_label = "T1" if tier1 else "T2"
                    exclude_ids.add(qseqid)
                    exclude_details[qseqid] = (
                        f"[{tier_label}] pident={pident:.1f}%, qcovs={qcovs:.0f}%"
                    )

    if exclude_ids:
        filtered_db = tier_dir / "transgene_db_clean.fa"
        n_total = 0
        n_excluded = 0
        with open(transgene_db) as fin, open(filtered_db, "w") as fout:
            write = True
            for line in fin:
                if line.startswith(">"):
                    n_total += 1
                    seq_id = line[1:].split()[0]
                    if seq_id in exclude_ids:
                        n_excluded += 1
                        write = False
                    else:
                        write = True
                        fout.write(line)
                elif write:
                    fout.write(line)
        log(f"  Host-endogenous exclusion: removed {n_excluded}/{n_total} entries")
        for eid in sorted(exclude_ids):
            log(f"    excluded: {eid} ({exclude_details[eid]})")
        subprocess.run(
            ["makeblastdb", "-in", str(filtered_db), "-dbtype", "nucl"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        blast_db = filtered_db
    else:
        log("  Host-endogenous exclusion: no entries matched host genome")

    # Write for debugging/auditability (Phase 4 receives exclude_ids directly)
    persist_path = workdir / "host_endogenous_ids.txt"
    with open(persist_path, "w") as fh:
        for eid in sorted(exclude_ids):
            fh.write(f"{eid}\t{exclude_details.get(eid, '')}\n")

    return blast_db, exclude_ids


_SRC_TIER = {
    "element_db":    2,
    "payload":       2,
    "sample_contig": 2,
    "univec":        1,
}

# Tier-2 sources collectively count as "element hits" for transgene_positive
# classification. Derived from _SRC_TIER so adding a new tier-2 source (e.g.,
# future 'crispr_ref') stays in sync automatically. Required because cd-hit
# -c 0.95 clustering of gmo_combined_db_v2.fa absorbed multiple element_db
# amplicons into payload-tagged representatives (cluster 0 P-CaMV35S, 35
# nptII, 54 hpt, 71 T-ocs, 77 P-nos in element_db/gmo_combined_db_v2.fa.clstr);
# a literal 'element_db' equality check at classify_site_tiers would drop
# rice_G281 Chr3:16,439,674 (CaMV35S-containing) -> AC-1 regression.
_TIER2_SRCS = frozenset(k for k, v in _SRC_TIER.items() if v >= 2)

# Module-level cache so we only warn once per unknown tag per run.
_UNKNOWN_SRC_WARNED: set[str] = set()


# _parse_src_tag moved to scripts/s05/primitives.py
# (Issue #4 Session 3) and re-imported at the top of this module for
# backward compatibility with existing callers.


def _should_replace(existing: dict | None, new_src: str, new_bit: float) -> bool:
    """Tag-tiered merge priority for classify_site_tiers `hits` map (Task T5).

    Tier-2 sources (element_db / payload / sample_contig) collectively beat
    tier-1 (univec). Within a tier, a strictly greater bitscore wins;
    ties go to the incumbent (matches the original `>` semantics and the
    BUG-3 regression guard in test_extra_element_db).

    The `existing` dict may use either the new `src`/`bit` schema or the
    legacy `source`/`bitscore` schema - both are supported so the historic
    callers in classify_site_tiers keep working while migration proceeds.
    """
    # I-1: warn once per run when an unknown tag shows up. Treated as tier 0
    # so it will never beat a known tag, which is the safe default.
    if new_src and new_src not in _SRC_TIER and new_src not in _UNKNOWN_SRC_WARNED:
        _UNKNOWN_SRC_WARNED.add(new_src)
        print(f"[s05] warn: unknown src tag {new_src!r}, treating as tier 0",
              file=sys.stderr)
    if existing is None:
        return True
    old_src = existing.get("src", existing.get("source", ""))
    old_bit = existing.get("bit", existing.get("bitscore", 0.0))
    new_tier = _SRC_TIER.get(new_src, 0)
    old_tier = _SRC_TIER.get(old_src, 0)
    if new_tier > old_tier:
        return True
    if new_tier < old_tier:
        return False
    # same tier -> strict '>' bitscore
    return new_bit > old_bit


def classify_site_tiers(
    sites: list[InsertionSite],
    element_db: Path,
    host_ref: Path,
    workdir: Path,
    threads: int = 4,
    min_identity: float = 80.0,
    min_aln_len: int = 20,
    extra_transgene_dbs: list[Path] | None = None,
) -> tuple[list[InsertionSite], list[InsertionSite], list[TierResult], set[str]]:
    """Classify sites by transgene-positive identification.

    BLASTs all clip sequences against transgene_db (element DB + UniVec)
    using blastn-short (optimized for 20-80bp queries).

    TRANSGENE-POSITIVE: at least one clip hits transgene_db → assemble
    TRANSGENE-NEGATIVE: no clip hits transgene_db → skip

    Returns (assembly_sites, skip_sites, all_tier_results, host_endo_ids).
    """
    if not sites:
        return [], [], [], set()

    tier_dir = workdir / "_tier_classification"
    tier_dir.mkdir(parents=True, exist_ok=True)

    # Locate transgene_db (element_db + UniVec combined)
    transgene_db = element_db.parent / "transgene_db.fa"
    if not transgene_db.exists():
        # Build transgene_db by combining element_db + UniVec
        univec_db = Path(__file__).resolve().parent.parent.parent / "db" / "univec_vectors.fa"
        log("  Building transgene_db (element DB + UniVec)...")
        with open(transgene_db, "w") as out:
            if element_db.exists():
                out.write(element_db.read_text())
            if univec_db.exists():
                out.write(univec_db.read_text())
        # Build BLAST index
        subprocess.run(
            ["makeblastdb", "-in", str(transgene_db), "-dbtype", "nucl"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )

    # Ensure BLAST DB index exists
    if not transgene_db.with_suffix(".fa.ndb").exists() and \
       not Path(str(transgene_db) + ".ndb").exists():
        subprocess.run(
            ["makeblastdb", "-in", str(transgene_db), "-dbtype", "nucl"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )

    blast_db, exclude_ids = _filter_host_endogenous(
        transgene_db, host_ref, tier_dir, workdir, threads
    )

    # Step 1: Collect all clip sequences
    clip_seqs: dict[str, str] = {}
    for site in sites:
        if site.seed_5p and len(site.seed_5p) >= 15:
            clip_seqs[f"{site.site_id}_5p"] = site.seed_5p
        if site.seed_3p and len(site.seed_3p) >= 15:
            clip_seqs[f"{site.site_id}_3p"] = site.seed_3p

    log(f"  BLASTing {len(clip_seqs)} clip sequences against transgene_db "
        f"(blastn-short, {blast_db.name})...")

    # Step 2: Write clips and BLAST against transgene_db
    clip_fa = tier_dir / "all_clips.fa"
    with open(clip_fa, "w") as fh:
        for name, seq in clip_seqs.items():
            fh.write(f">{name}\n{seq}\n")

    blast_out = tier_dir / "transgene_blast.tsv"
    result = subprocess.run(
        ["blastn", "-task", "blastn-short",
         "-query", str(clip_fa), "-db", str(blast_db),
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
         "-evalue", "1e-5", "-max_target_seqs", "3",
         "-num_threads", str(threads),
         "-out", str(blast_out)],
        stderr=subprocess.DEVNULL,
    )
    if result.returncode != 0:
        log(f"  WARNING: transgene BLAST failed (rc={result.returncode})")

    # Step 3: Parse hits — best hit per query with identity/length filter
    hits: dict[str, dict] = {}
    if blast_out.exists():
        with open(blast_out) as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 6:
                    continue
                qname = cols[0]
                pident = float(cols[2])
                aln_len = int(cols[3])
                bitscore = float(cols[5])

                if pident < min_identity or aln_len < min_aln_len:
                    continue

                sseqid = cols[1]
                # T5: 4-way source tag. v2 DB headers carry a |src=<tag>
                # suffix (payload / element_db / sample_contig); legacy
                # univec|... headers get mapped to 'univec', everything
                # else to 'element_db'.
                default_main = (
                    "univec" if sseqid.startswith("univec|") else "element_db"
                )
                element_id, src_tag = _parse_src_tag(sseqid, default_main)
                if _should_replace(hits.get(qname), src_tag, bitscore):
                    hits[qname] = {
                        "element": element_id,
                        "identity": pident,
                        "aln_length": aln_len,
                        "bitscore": bitscore,
                        "source": src_tag,
                    }

    log(f"  {len(hits)} clips hit transgene_db (identity>={min_identity}%, "
        f"aln>={min_aln_len}bp)")

    # Step 3b: Optional BLAST against extra transgene DBs. These include
    # the always-on shared payload DB (common_payload.fa) and per-sample
    # SPAdes contigs from s04b. Merged into `hits` so both shared and
    # sample-specific transgene references contribute to transgene-positive
    # classification.
    extra_dbs_list = [
        edb for edb in (extra_transgene_dbs or [])
        if edb is not None and edb.exists() and edb.stat().st_size > 0
    ]
    for i, edb in enumerate(extra_dbs_list):
        extra_blast_out = tier_dir / f"transgene_blast_extra_{i}_{edb.stem}.tsv"
        # T5 fallback tag when a header lacks |src=<tag>: per-sample SPAdes
        # contig FASTAs from s04b are labelled 'sample_contig'; everything
        # else (legacy common_payload, ad-hoc FASTAs) defaults to
        # 'element_db' so tier-2 priority still holds.
        default_tag = "sample_contig" if "contigs" in edb.name else "element_db"
        log(f"  BLASTing clips against extra transgene DB "
            f"({edb.name}, default_src={default_tag}, blastn-short)...")
        result_extra = subprocess.run(
            ["blastn", "-task", "blastn-short",
             "-query", str(clip_fa), "-subject", str(edb),
             "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
             "-evalue", "1e-5", "-max_target_seqs", "3",
             "-out", str(extra_blast_out)],
            stderr=subprocess.DEVNULL,
        )
        if result_extra.returncode != 0:
            log(f"  WARNING: extra transgene BLAST failed "
                f"(rc={result_extra.returncode})")
        n_extra = 0
        if extra_blast_out.exists():
            with open(extra_blast_out) as fh:
                for line in fh:
                    cols = line.strip().split("\t")
                    if len(cols) < 6:
                        continue
                    qname = cols[0]
                    pident = float(cols[2])
                    aln_len = int(cols[3])
                    bitscore = float(cols[5])
                    if pident < min_identity or aln_len < min_aln_len:
                        continue
                    sseqid = cols[1]
                    element_id, src_tag = _parse_src_tag(sseqid, default_tag)
                    if _should_replace(hits.get(qname), src_tag, bitscore):
                        hits[qname] = {
                            "element": element_id,
                            "identity": pident,
                            "aln_length": aln_len,
                            "bitscore": bitscore,
                            "source": src_tag,
                        }
                        n_extra += 1
        log(f"  Extra transgene DB {edb.name} contributed/updated "
            f"{n_extra} clip hits (total now {len(hits)})")

    # Step 3.5: Per-clip host verification
    # Some divergent host elements (e.g., I-actin rice intron) escape DB-level
    # exclusion because their full-entry qcovs is low, but the SHORT clip we
    # match against them happens to be a 100% conserved window. To catch this,
    # BLAST each transgene-hit clip back against the host genome with
    # blastn-short. If the clip ALSO hits the host at >=95% identity over
    # >=30bp, it's an endogenous conserved sequence — drop it.
    host_blast_db_exists = (
        host_ref.with_suffix(".fa.ndb").exists()
        or Path(str(host_ref) + ".ndb").exists()
    )
    if hits and host_blast_db_exists:
        hit_clip_fa = tier_dir / "hit_clips.fa"
        with open(hit_clip_fa, "w") as fh:
            for qname in hits:
                fh.write(f">{qname}\n{clip_seqs[qname]}\n")

        host_clip_blast = tier_dir / "hit_clips_vs_host.tsv"
        log(f"  Per-clip host verification: BLAST {len(hits)} hit clips "
            f"vs host genome (blastn-short)...")
        result = subprocess.run(
            ["blastn", "-task", "blastn-short",
             "-query", str(hit_clip_fa), "-db", str(host_ref),
             "-outfmt", "6 qseqid sseqid pident length evalue",
             "-evalue", "1e-5", "-max_target_seqs", "1",
             "-num_threads", str(threads),
             "-out", str(host_clip_blast)],
            stderr=subprocess.DEVNULL,
        )
        if result.returncode != 0:
            log(f"  WARNING: per-clip host BLAST failed (rc={result.returncode})")

        host_endogenous_clips: dict[str, tuple[float, int]] = {}
        if host_clip_blast.exists():
            with open(host_clip_blast) as fh:
                for line in fh:
                    cols = line.strip().split("\t")
                    if len(cols) < 5:
                        continue
                    qname = cols[0]
                    pident = float(cols[2])
                    length = int(cols[3])
                    if pident >= CLIP_HOST_PIDENT and length >= CLIP_HOST_MIN_LEN:
                        # Keep best (longest) hit per clip
                        prev = host_endogenous_clips.get(qname)
                        if prev is None or length > prev[1]:
                            host_endogenous_clips[qname] = (pident, length)

        if host_endogenous_clips:
            log(f"  Per-clip host filter: dropping {len(host_endogenous_clips)} "
                f"clips that match host (pident>={CLIP_HOST_PIDENT}%, length>={CLIP_HOST_MIN_LEN}bp)")
            for qname, (pid, ln) in sorted(host_endogenous_clips.items()):
                hit = hits[qname]
                log(f"    dropped: {qname} → {hit['element']} "
                    f"(host hit: {pid:.0f}%/{ln}bp)")
                del hits[qname]
            log(f"  After host verification: {len(hits)} transgene-only clips remain")
        else:
            log("  Per-clip host filter: no clips matched host")

    # Step 4: Classify each site
    tier_results: list[TierResult] = []
    assembly_sites: list[InsertionSite] = []
    skip_sites: list[InsertionSite] = []

    for site in sites:
        key_5p = f"{site.site_id}_5p"
        key_3p = f"{site.site_id}_3p"

        hit_5p = hits.get(key_5p, {})
        hit_3p = hits.get(key_3p, {})

        # Require at least one tier-2 element hit (not univec-only).
        # UniVec-only matches at this step are typically short (20-31bp)
        # alignments that match native plant DNA by chance. Real T-DNA
        # insertions always have at least one characteristic element
        # (promoter, selection marker, terminator) matching a tier-2 source.
        #
        # C-1 fix: tier-2 sources (element_db / payload / sample_contig) all
        # count. After cd-hit -c 0.95 clustering, former element_db amplicons
        # may be absorbed into payload-tagged representatives (e.g., CaMV35S
        # cluster 0 in gmo_combined_db_v2.fa.clstr). A literal "element_db"
        # check would miss rice_G281 Chr3:16,439,674 CaMV35S hits -> AC-1
        # regression. Use the _TIER2_SRCS module constant derived from
        # _SRC_TIER for a single source of truth.
        has_element_hit = (hit_5p.get("source") in _TIER2_SRCS) or \
                           (hit_3p.get("source") in _TIER2_SRCS)
        is_positive = has_element_hit

        tr = TierResult(
            site_id=site.site_id,
            chrom=site.host_chr,
            pos=site.pos_5p or site.pos_3p or 0,
            transgene_positive=is_positive,
            clip_5p_len=len(site.seed_5p) if site.seed_5p else 0,
            clip_3p_len=len(site.seed_3p) if site.seed_3p else 0,
            hit_5p=hit_5p.get("element", ""),
            hit_5p_identity=hit_5p.get("identity", 0),
            hit_5p_aln_len=hit_5p.get("aln_length", 0),
            hit_5p_source=hit_5p.get("source", ""),
            hit_3p=hit_3p.get("element", ""),
            hit_3p_identity=hit_3p.get("identity", 0),
            hit_3p_aln_len=hit_3p.get("aln_length", 0),
            hit_3p_source=hit_3p.get("source", ""),
        )
        tier_results.append(tr)

        if is_positive:
            assembly_sites.append(site)
        else:
            skip_sites.append(site)

    shutil.rmtree(tier_dir, ignore_errors=True)

    n_pos = len(assembly_sites)
    n_neg = len(skip_sites)
    log(f"  Transgene-positive (assemble): {n_pos}")
    log(f"  Transgene-negative (skip):     {n_neg}")

    # Log positive site details
    for tr in tier_results:
        if tr.transgene_positive:
            parts = []
            if tr.hit_5p:
                parts.append(f"5p={tr.hit_5p} ({tr.hit_5p_identity:.0f}%/{tr.hit_5p_aln_len}bp)")
            if tr.hit_3p:
                parts.append(f"3p={tr.hit_3p} ({tr.hit_3p_identity:.0f}%/{tr.hit_3p_aln_len}bp)")
            log(f"    {tr.site_id} {tr.chrom}:{tr.pos}: {', '.join(parts)}")

    return assembly_sites, skip_sites, tier_results, exclude_ids


def write_tier_classification(
    tier_results: list[TierResult],
    output_path: Path,
) -> None:
    """Write site_tier_classification.tsv."""
    with open(output_path, "w") as fh:
        fh.write("site_id\tchrom\tpos\ttransgene_positive\t"
                 "clip_5p_len\tclip_3p_len\t"
                 "hit_5p\thit_5p_identity\thit_5p_aln_len\thit_5p_source\t"
                 "hit_3p\thit_3p_identity\thit_3p_aln_len\thit_3p_source\n")
        for tr in tier_results:
            fh.write(f"{tr.site_id}\t{tr.chrom}\t{tr.pos}\t{tr.transgene_positive}\t"
                     f"{tr.clip_5p_len}\t{tr.clip_3p_len}\t"
                     f"{tr.hit_5p}\t{tr.hit_5p_identity}\t{tr.hit_5p_aln_len}\t{tr.hit_5p_source}\t"
                     f"{tr.hit_3p}\t{tr.hit_3p_identity}\t{tr.hit_3p_aln_len}\t{tr.hit_3p_source}\n")
    log(f"  Classification written: {output_path}")
