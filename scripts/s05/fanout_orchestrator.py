"""Phase 6 — Fanout Orchestrator (Issue #4, Session 6).

Extracted from ``scripts/s05_insert_assembly.py``. Owns the CLI dispatch
for s05's four phases (1_1.5 site discovery, 2_3 per-site assembly,
4 annotation, 'all' end-to-end). Stage 9 in the s05 DAG (top of pipeline).

main() splits the original monolith ``main()`` into:
- main: argparse + phase dispatch (~80 lines)
- _run_phase_1_1_5: Phase 1 + 1.5 (site discovery + classify)
- _run_phase_2_3: Phase 2 + 3 (read extraction + assembly per site)
- _run_phase_4: Phase 4 (annotation, FP filter eval, report writing)

The split is mechanical — variables and control flow preserved verbatim
from the monolith. Snapshot fixtures + DAG test verify no behavioral drift.
"""
from __future__ import annotations

import argparse
import gzip
import json
import pickle
import sys
from collections import Counter
from pathlib import Path

from .primitives import (
    STEP,
    log,
    read_fasta,
    InsertionSite,
    TierResult,
)
from .site_discovery import (
    find_softclip_junctions,
    _apply_mask_bed,
    MASKED_SOURCE_TAG,
    parse_legacy_junctions,
    legacy_junctions_to_sites,
)
from .classify import (
    classify_site_tiers,
    write_tier_classification,
)
from .read_extraction import (
    extract_candidate_reads,
)
from .assembly import (
    assemble_insert,
)
from .annotation import (
    annotate_insert,
)
from .filter_b_flank import (
    _find_construct_flanking_regions,
)
from .report import (
    generate_report,
    write_stats,
)
from .config_loader import load_verdict_rules


def _run_phase_1_1_5(
    args: argparse.Namespace,
    step_dir: Path,
    host_bam: Path,
    host_ref: Path,
    element_db: Path,
    extra_dbs: list[Path],
    positive_sites_json: Path,
    positive_sites_pkl: Path,
) -> tuple[list[InsertionSite], set[str]] | None:
    """Phase 1 + 1.5: site discovery + transgene-positive classify.

    Returns (sites, host_endo_ids) on success, or None if called for a
    phase that requires loading from pickle (2_3 / 4).  Early returns
    (no sites found, no transgene-positive sites) write stats and return
    ([], set()) sentinel rather than None.
    """
    sites: list[InsertionSite] = []
    host_endo_ids: set[str] = set()

    if args.phase in ("all", "1_1.5"):
        # ---- Phase 1: Soft-clip junction detection ----
        sites = find_softclip_junctions(
            host_bam, host_ref, element_db, step_dir,
            min_clip=args.min_clip,
            extra_dbs=extra_dbs,
        )

        # Fallback to step 6 junctions if no sites found
        if not sites and args.junctions:
            junctions_path = Path(args.junctions)
            if junctions_path.exists() and junctions_path.stat().st_size > 0:
                log("No soft-clip sites found, falling back to legacy junctions file...")
                legacy_juncs = parse_legacy_junctions(junctions_path)
                if legacy_juncs:
                    sites = legacy_junctions_to_sites(legacy_juncs, host_bam)
                    log(f"  Converted {len(legacy_juncs)} legacy junctions → "
                        f"{len(sites)} sites")

        if not sites:
            log("No insertion sites found. Nothing to assemble.")
            write_stats(step_dir / "s05_stats.txt", args.sample_name, 0, 0, [])
            # Still emit an empty positive_sites.json so the array wrapper
            # can detect 'nothing to fan out' cleanly.
            positive_sites_json.write_text("[]\n")
            positive_sites_pkl.write_bytes(pickle.dumps(
                {"sites": [], "host_endo_ids": set()}
            ))
            return [], set()

        log(f"\nPhase 1 found {len(sites)} insertion site(s):")
        for site in sites:
            log(f"  {site.site_id}: {site.host_chr}:{site.pos_5p}"
                f"{f'-{site.pos_3p}' if site.pos_3p else ''} "
                f"({site.confidence}, 5p={len(site.seed_5p)}bp, "
                f"3p={len(site.seed_3p)}bp)")

        # ---- T10: pre-mask BED intersect (host-endogenous ortholog regions) ----
        # Sites inside a T9 mask region are tagged FALSE_NEGATIVE_MASKED and
        # excluded from Phase 1.5 CANDIDATE classification. They still appear
        # in site_tier_classification.tsv for regulatory audit trail.
        masked_sites: list[InsertionSite] = []
        if args.mask_bed:
            bed_path = Path(args.mask_bed)
            _apply_mask_bed(sites, bed_path)
            masked_sites = [
                s for s in sites
                if getattr(s, "tier", None) == MASKED_SOURCE_TAG
            ]
            sites = [
                s for s in sites
                if getattr(s, "tier", None) != MASKED_SOURCE_TAG
            ]
            print(
                f"[Phase 1] mask-bed {bed_path.name}: "
                f"{len(masked_sites)}/{len(masked_sites) + len(sites)} "
                f"sites tagged {MASKED_SOURCE_TAG}",
                file=sys.stderr,
            )
            for s in masked_sites:
                log(f"  MASKED: {s.site_id} {s.host_chr}:{s.pos_5p} "
                    f"→ {getattr(s, 'mask_element', 'masked')}")

        # ---- Phase 1.5: Transgene-positive identification ----
        log("\n--- Phase 1.5: Transgene-positive identification ---")
        assembly_sites, skip_sites, tier_results, host_endo_ids = classify_site_tiers(
            sites, element_db, host_ref, step_dir,
            threads=args.threads,
            extra_transgene_dbs=extra_dbs,
        )
        # T10: append FALSE_NEGATIVE_MASKED entries to tier_results for audit
        for s in masked_sites:
            tier_results.append(TierResult(
                site_id=s.site_id,
                chrom=s.host_chr,
                pos=s.pos_5p or s.pos_3p or 0,
                transgene_positive=False,
                clip_5p_len=len(s.seed_5p) if s.seed_5p else 0,
                clip_3p_len=len(s.seed_3p) if s.seed_3p else 0,
                hit_5p=getattr(s, "mask_element", ""),
                hit_5p_source=MASKED_SOURCE_TAG,
                hit_3p=getattr(s, "mask_element", ""),
                hit_3p_source=MASKED_SOURCE_TAG,
            ))
        write_tier_classification(
            tier_results, step_dir / "site_tier_classification.tsv",
        )

        if not assembly_sites:
            log("No transgene-positive sites found. Nothing to assemble.")
            write_stats(step_dir / "s05_stats.txt", args.sample_name, len(sites), 0, [])
            positive_sites_json.write_text("[]\n")
            positive_sites_pkl.write_bytes(pickle.dumps(
                {"sites": [], "host_endo_ids": host_endo_ids}
            ))
            return [], host_endo_ids

        sites = assembly_sites  # Only assemble transgene-positive sites

        # ---- T8: persist Phase 1/1.5 state for array fan-out -----------------
        # positive_sites.json is the *public* contract for submit_s05_array.sh
        # to enumerate site_ids for --array indexing.  positive_sites.pkl is
        # the *internal* rehydration file for Phase 2_3 / 4 workers so they
        # can run against the exact same InsertionSite objects (including
        # seed_5p/3p consensus clips and JunctionCluster detail).
        positive_sites_json.write_text(json.dumps(
            [
                {
                    "site_id": s.site_id,
                    "chr": s.host_chr,
                    "pos": s.pos_5p,
                    "pos_3p": s.pos_3p,
                    "confidence": s.confidence,
                }
                for s in sites
            ],
            indent=2,
        ) + "\n")
        # SECURITY: pickle is consumed only from the same sample's step_dir
        # (same user writes and reads). Never rehydrate from untrusted sources.
        # v1.1 module split will replace with typed JSON + dataclasses loader.
        positive_sites_pkl.write_bytes(pickle.dumps(
            {"sites": sites, "host_endo_ids": host_endo_ids}
        ))
        log(f"  [T8] wrote {len(sites)} positive sites → "
            f"{positive_sites_json.name} (+ .pkl)")

        if args.phase == "1_1.5":
            log(f"\n[Phase 1.5] early return — "
                f"{len(sites)} positive sites saved for per-site fan-out.")
            return sites, host_endo_ids

    else:
        # phase in {"2_3", "4"} — rehydrate from pickle written by a prior
        # 1_1.5 / all run.  Fail loud if the pickle is missing so the user
        # notices they skipped Phase 1.5.
        if not positive_sites_pkl.exists():
            sys.exit(
                f"ERROR: --phase {args.phase} requires a prior "
                f"--phase 1_1.5 (or --phase all) run to have produced "
                f"{positive_sites_pkl} — not found."
            )
        loaded = pickle.loads(positive_sites_pkl.read_bytes())
        sites = loaded["sites"]
        host_endo_ids = loaded["host_endo_ids"]
        log(f"  [T8] loaded {len(sites)} positive sites + "
            f"{len(host_endo_ids)} host-endogenous IDs from {positive_sites_pkl.name}")

    return sites, host_endo_ids


def _run_phase_2_3(
    args: argparse.Namespace,
    step_dir: Path,
    host_bam: Path,
    host_ref: Path,
    element_db: Path,
    sites: list[InsertionSite],
) -> list[dict] | None:
    """Phase 2 + 3: per-site candidate read extraction + iterative assembly.

    Returns all_results list, or None if phase was skipped (phase=4 only).
    Returns early (None) when --phase=2_3 finishes a single site.
    """
    # ---- Phase 2 + 3: Per-site candidate extraction + iterative assembly ----
    all_results: list[dict] = []
    if args.phase in ("all", "2_3"):
        if args.phase == "2_3":
            # Filter to single site for array task.
            matching = [s for s in sites if s.site_id == args.site_id]
            if not matching:
                sys.exit(
                    f"ERROR: --site-id {args.site_id!r} not in "
                    f"positive_sites.json (have: "
                    f"{[s.site_id for s in sites]})"
                )
            sites_to_run = matching
        else:
            sites_to_run = sites

        for site in sites_to_run:
            log(f"\n{'=' * 60}")
            log(f"=== Processing {site.site_id}: "
                f"{site.host_chr}:{site.pos_5p} ===")

            # Phase 2: Extract candidate reads
            cand_r1 = step_dir / f"{site.site_id}_candidate_R1.fastq.gz"
            cand_r2 = step_dir / f"{site.site_id}_candidate_R2.fastq.gz"
            if not cand_r1.exists():
                n_cand = extract_candidate_reads(
                    host_bam, site, cand_r1, cand_r2,
                    flank=args.flank, threads=args.threads,
                    min_clip=args.min_clip,
                    junction_window=args.junction_window,
                )
                log(f"  Candidate reads: {n_cand:,} pairs")
            else:
                n_cand = 0
                opener = gzip.open if str(cand_r1).endswith(".gz") else open
                with opener(str(cand_r1), "rt") as fh:
                    for _ in fh:
                        n_cand += 1
                n_cand = n_cand // 4
                log(f"  Candidate reads (cached): {n_cand:,} pairs")

            # Phase 3: Iterative assembly
            s03_r1 = Path(args.s03_r1) if args.s03_r1 else None
            s03_r2 = Path(args.s03_r2) if args.s03_r2 else None
            assembly_fa, rounds, status = assemble_insert(
                site=site,
                candidate_r1=cand_r1,
                candidate_r2=cand_r2,
                host_bam=host_bam,
                host_ref=host_ref,
                element_db=element_db,
                workdir=step_dir,
                threads=args.threads,
                max_rounds=args.max_rounds,
                ext_k=args.seed_k,
                recruit_k=args.recruit_k,
                gap_size=args.gap_size,
                s03_r1=s03_r1,
                s03_r2=s03_r2,
            )

            # Record result
            if assembly_fa and assembly_fa.exists():
                contigs = read_fasta(assembly_fa)
                if contigs:
                    longest_name = max(contigs, key=lambda k: len(contigs[k]))
                    longest_seq = contigs[longest_name]
                    insert_len = len(longest_seq)
                    insert_ns = longest_seq.upper().count("N")
                else:
                    insert_len, insert_ns = 0, 0
                    status = "no_assembly"
            else:
                insert_len, insert_ns = 0, 0
                status = "no_assembly"

            result = {
                "site_id": site.site_id,
                "insert_length": insert_len,
                "remaining_ns": insert_ns,
                "rounds": rounds,
                "status": status,
            }
            all_results.append(result)

            # T8: per-site runs drop a small JSON sidecar so Phase 4
            # can rebuild all_results without re-running assembly.
            if args.phase == "2_3":
                (step_dir / f"{site.site_id}_result.json").write_text(
                    json.dumps(result, indent=2) + "\n"
                )

        if args.phase == "2_3":
            log(f"\n[Phase 2_3] done for site {args.site_id}. "
                f"Phase 4 collects results after the full array completes.")
            return None

    return all_results


def _run_phase_4(
    args: argparse.Namespace,
    step_dir: Path,
    host_bam: Path,
    host_ref: Path,
    element_db: Path,
    extra_dbs: list[Path],
    sites: list[InsertionSite],
    host_endo_ids: set[str],
    all_results: list[dict],
) -> None:
    """Phase 4: annotation + FP filter eval + report writing.

    Modifies all_results in place to attach verdict per site.
    """
    # ---- Phase 4: Annotation & Report ----
    # For --phase=4 we rehydrate all_results from per-site sidecars dropped
    # by each --phase=2_3 worker.  For --phase=all we already built
    # all_results above.
    if args.phase == "4":
        for site in sites:
            sidecar = step_dir / f"{site.site_id}_result.json"
            if sidecar.exists():
                all_results.append(json.loads(sidecar.read_text()))
            else:
                log(f"  [T8 Phase 4] missing sidecar for {site.site_id}, "
                    f"marking as no_assembly")
                all_results.append({
                    "site_id": site.site_id,
                    "insert_length": 0,
                    "remaining_ns": 0,
                    "rounds": 0,
                    "status": "no_assembly",
                })

    # ---- Combine inserts ----
    combined_insert = step_dir / "insert_only.fasta"
    with open(combined_insert, "w") as fout:
        for site in sites:
            site_fa = step_dir / f"{site.site_id}_insert.fasta"
            if site_fa.exists():
                fout.write(open(site_fa).read())

    # ---- Phase 4: Annotation & Report ----
    if combined_insert.exists() and combined_insert.stat().st_size > 0:
        log(f"\n{'=' * 60}")
        log("=== Phase 4: Annotation ===")
        ann_tsv, border_tsv = annotate_insert(
            combined_insert, element_db, step_dir, args.sample_name,
            no_remote_blast=args.no_remote_blast,
            extra_dbs=extra_dbs,
        )

        # host_endo_ids returned directly from classify_site_tiers (no file I/O)
        log(f"  Using {len(host_endo_ids)} host-endogenous element IDs "
            f"for post-filter verdict")

        # Detect construct-flanking regions (host DNA in construct reference)
        construct_flanking: list[tuple[str, int, int]] = []
        if args.construct_ref:
            construct_path = Path(args.construct_ref)
            if construct_path.exists():
                log("  Detecting construct-flanking host regions...")
                construct_flanking = _find_construct_flanking_regions(
                    construct_path, host_ref, step_dir, threads=args.threads,
                )

        # Generate report for each site
        verdict_counts: Counter = Counter()
        verdict_by_site: dict[str, str] = {}
        # Issue #3 wire-in: load canonical-triplet + threshold rules from config
        # (falls back to DEFAULT_TRIPLETS when --config is omitted / missing).
        rules = load_verdict_rules(
            Path(args.config) if args.config else Path("config.yaml"),
            args.sample_name,
        )
        for site in sites:
            site_fa = step_dir / f"{site.site_id}_insert.fasta"
            if site_fa.exists():
                result_info = next(
                    (r for r in all_results if r["site_id"] == site.site_id), None)
                if result_info:
                    _, verdict = generate_report(
                        site_fa, ann_tsv, border_tsv,
                        site, result_info["rounds"], result_info["status"],
                        step_dir,
                        host_endogenous_ids=host_endo_ids,
                        host_ref=host_ref,
                        element_db=element_db,
                        threads=args.threads,
                        construct_flanking=construct_flanking,
                        rules=rules,
                    )
                    verdict_counts[verdict] += 1
                    verdict_by_site[site.site_id] = verdict

        # Attach verdict onto the per-site results so the final summary reflects it.
        for r in all_results:
            r["verdict"] = verdict_by_site.get(r["site_id"], "NO_ASSEMBLY")
        log(f"  Phase 4 verdicts: "
            f"{dict(verdict_counts) if verdict_counts else '(none)'}")


def main() -> None:
    """CLI entrypoint — argparse + phase dispatch."""
    parser = argparse.ArgumentParser(
        description="Step 5: Targeted insert assembly "
                    "(soft-clip detection + k-mer extension + Pilon gap fill)")
    parser.add_argument("--junctions", default=None,
                        help="DIRECT-INVOCATION ESCAPE HATCH (run_pipeline.py never sets this). "
                             "junctions.tsv from step 6 (fallback if soft-clip "
                             "detection finds nothing)")
    parser.add_argument("--host-bam", required=True,
                        help="Host-mapped BAM from step 7")
    parser.add_argument("--host-ref", required=True,
                        help="Host reference FASTA")
    parser.add_argument("--element-db", required=True,
                        help="Element database FASTA for annotation")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--max-rounds", type=int, default=5,
                        help="Max alternating extension+Pilon rounds (T3: reduced from 8 per max_rounds_study.md)")
    parser.add_argument("--seed-k", type=int, default=15,
                        help="k-mer size for seed extension")
    parser.add_argument("--recruit-k", type=int, default=25,
                        help="k-mer size for unmapped read recruitment")
    parser.add_argument("--gap-size", type=int, default=1000,
                        help="Initial N gap size for Pilon scaffold")
    parser.add_argument("--flank", type=int, default=5000,
                        help="Flank size for candidate read extraction (bp)")
    parser.add_argument("--min-clip", type=int, default=20,
                        help="Minimum soft-clip length for junction detection")
    parser.add_argument("--extra-element-db", type=Path, default=None,
                        help="Optional second FASTA for element annotation "
                             "(e.g., per-sample SPAdes construct contigs).")
    parser.add_argument("--common-payload-db", type=Path, default=None,
                        help="Always-on supplementary FASTA of common transgene "
                             "payload genes (bar, nptII, hpt, gusA, gfp, ...).")
    parser.add_argument("--junction-window", type=int, default=20,
                        help="Window (bp) around junction for allele phasing "
                             "(reads within this distance classified as junction-proximal)")
    parser.add_argument("--construct-ref", default=None,
                        help="Construct reference FASTA (for flanking detection)")
    parser.add_argument("--no-remote-blast", action="store_true",
                        help="Skip remote NCBI nt BLAST (use local element_db only)")
    parser.add_argument(
        "--mask-bed",
        default=None,
        help="BED file of host-endogenous ortholog regions (T9/T10). "
             "Sites overlapping a region are tagged FALSE_NEGATIVE_MASKED "
             "(not dropped) after Phase 1 discovery, preventing CANDIDATE "
             "promotion while preserving the audit trail.",
    )
    parser.add_argument("--s03-r1", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--s03-r2", default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        "--config",
        default=None,
        help="Path to config.yaml for verdict rules / canonical triplets "
             "(Issue #3 wire-in). Missing file falls back to DEFAULT_TRIPLETS.",
    )
    # --- T8: per-site SLURM array fan-out (v1.0 임시안) -----------------------
    # v1.1 에서 full module split 완료 후 이 flag 제거 예정.
    parser.add_argument(
        "--phase",
        choices=["all", "1_1.5", "2_3", "4"],
        default="all",
        help="Run specific phase(s). all=default sequential (backwards compat). "
             "1_1.5=site discovery + transgene-positive classify, dumps "
             "positive_sites.json and exits. 2_3=per-site candidate read "
             "extraction + iterative assembly (requires --site-id). "
             "4=annotate + per-site report (consumes per-site "
             "<site_id>_insert.fasta dropped by 2_3 runs).",
    )
    parser.add_argument(
        "--site-id",
        default=None,
        help="For --phase=2_3: run assembly for only the given site_id "
             "(as stored in positive_sites.json).",
    )
    args = parser.parse_args()

    if args.phase == "2_3" and not args.site_id:
        parser.error("--site-id is required when --phase=2_3")

    step_dir = Path(args.outdir) / args.sample_name / STEP
    step_dir.mkdir(parents=True, exist_ok=True)

    host_bam = Path(args.host_bam)
    host_ref = Path(args.host_ref)
    element_db = Path(args.element_db)
    extra_db = args.extra_element_db
    common_payload_db = args.common_payload_db

    # Order: common_payload (shared) first, then per-sample extra.
    # A None entry is dropped naturally.
    extra_dbs = [p for p in [common_payload_db, extra_db] if p is not None]

    log(f"=== Step 5: Targeted Insert Assembly for {args.sample_name} ===")
    log(f"  Phase: {args.phase}"
        + (f" (site_id={args.site_id})" if args.site_id else ""))
    log(f"  Host BAM: {host_bam}")
    log(f"  Host ref: {host_ref}")
    log(f"  Element DB: {element_db}")
    if common_payload_db is not None:
        log(f"  Common payload DB: {common_payload_db}")
    if extra_db is not None:
        log(f"  Extra element DB: {extra_db}")

    # Paths used by T8 fan-out contract between phases.
    positive_sites_json = step_dir / "positive_sites.json"
    positive_sites_pkl = step_dir / "positive_sites.pkl"

    # ---- Phase 1 + 1.5 (site discovery + transgene-positive classify) ----
    phase1_result = _run_phase_1_1_5(
        args, step_dir, host_bam, host_ref, element_db, extra_dbs,
        positive_sites_json, positive_sites_pkl,
    )
    if phase1_result is None:
        return
    sites, host_endo_ids = phase1_result

    # Early exits: no sites or phase=1_1.5 with positive sites found
    if not sites or args.phase == "1_1.5":
        return

    # ---- Phase 2 + 3 ----
    phase23_result = _run_phase_2_3(
        args, step_dir, host_bam, host_ref, element_db, sites,
    )
    if phase23_result is None:
        # phase=2_3 single-site run completed; Phase 4 runs separately
        return
    all_results = phase23_result

    # ---- Phase 4: Annotation & Report ----
    _run_phase_4(
        args, step_dir, host_bam, host_ref, element_db, extra_dbs,
        sites, host_endo_ids, all_results,
    )

    # ---- Write stats ----
    write_stats(
        step_dir / "s05_stats.txt", args.sample_name,
        len(sites), 0, all_results,
    )

    log(f"\n{'=' * 60}")
    log(f"=== Step 5 complete ===")
    log(f"Output: {step_dir}")
    for r in all_results:
        v = r.get("verdict", "NO_ASSEMBLY")
        log(f"  {r['site_id']}: {r['insert_length']:,}bp, "
            f"{r['remaining_ns']} Ns, {r['rounds']} rounds → {r['status']} [{v}]")


if __name__ == "__main__":
    main()
