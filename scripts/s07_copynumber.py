#!/usr/bin/env python3
"""Step 7: Estimate transgene copy number from read depth ratio.

Compares median read depth on the construct BAM (Step 2) to the genome-wide
median depth from the host BAM (Step 4). The ratio estimates copy number:

    copy_number = construct_median_depth / genome_median_depth
    ~0.5 = hemizygous single-copy
    ~1.0 = homozygous single-copy

Cross-validates with candidate site count from Step 5 (1 site = single
locus, >1 may indicate multi-locus).

Output:
    {outdir}/{sample}/s07_copynumber/copynumber.tsv
"""

import argparse
import subprocess
import sys
from pathlib import Path


def log(msg: str) -> None:
    """Print message to stderr."""
    print(f"[s07_copynumber] {msg}", file=sys.stderr, flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Estimate transgene copy number from read depth ratio"
    )
    parser.add_argument("--construct-bam", type=Path, required=True,
                        help="Construct-mapped BAM (from s02_construct_map)")
    parser.add_argument("--host-bam", type=Path, required=True,
                        help="Host-mapped BAM (from s04_host_map)")
    parser.add_argument("--construct-ref", type=Path, required=True,
                        help="Construct reference FASTA")
    parser.add_argument("--host-ref", type=Path, required=True,
                        help="Host reference genome FASTA")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="Base output directory")
    parser.add_argument("--sample-name", type=str, required=True,
                        help="Sample name for output file naming")
    parser.add_argument("--site-stats", type=Path, default=None,
                        help="s05 stats file for site count cross-validation (optional)")
    return parser.parse_args()


def get_construct_median_depth(bam: Path) -> tuple[float, int, dict[str, list[int]]]:
    """Calculate median depth across covered construct elements.

    When using element_db (many reference sequences), most positions have
    0 coverage. We calculate per-element depth and use the element with
    the highest median coverage for copy number estimation.

    Returns:
        (median_depth, total_covered_positions, per_element_depths)
    """
    result = subprocess.run(
        ["samtools", "depth", "-a", str(bam)],
        capture_output=True, text=True, check=True,
    )

    # Group depths by reference element
    element_depths: dict[str, list[int]] = {}
    for line in result.stdout.strip().splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            ref_name = parts[0]
            depth = int(parts[2])
            if ref_name not in element_depths:
                element_depths[ref_name] = []
            element_depths[ref_name].append(depth)

    if not element_depths:
        return 0.0, 0, {}

    # Find element with highest median depth (best marker for copy number)
    best_element = None
    best_median = 0.0
    for elem, depths in element_depths.items():
        # Only consider elements with >50% positions covered
        covered = sum(1 for d in depths if d > 0)
        if covered < len(depths) * 0.5:
            continue
        depths_sorted = sorted(depths)
        median = depths_sorted[len(depths_sorted) // 2]
        if median > best_median:
            best_median = float(median)
            best_element = elem

    # If no element has >50% coverage, use all covered positions
    if best_element is None:
        all_covered = [d for depths in element_depths.values()
                       for d in depths if d > 0]
        if not all_covered:
            return 0.0, 0, element_depths
        all_covered.sort()
        median = all_covered[len(all_covered) // 2]
        return float(median), len(all_covered), element_depths

    total_covered = sum(1 for d in element_depths[best_element] if d > 0)
    log(f"  Best marker element: {best_element} "
        f"(median={best_median:.1f}x, {total_covered} covered positions)")
    return best_median, len(element_depths[best_element]), element_depths


def get_host_median_depth(bam: Path, max_positions: int = 1_000_000) -> tuple[float, int]:
    """Calculate median depth for the host genome BAM.

    Streams samtools depth output and samples up to max_positions for
    the median calculation to handle large genomes efficiently.

    Returns:
        (median_depth, total_positions_sampled)
    """
    depth_proc = subprocess.Popen(
        ["samtools", "depth", "-a", str(bam)],
        stdout=subprocess.PIPE, text=True,
    )

    depths: list[int] = []
    total_positions = 0
    sample_every = 1

    for line in depth_proc.stdout:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        total_positions += 1

        if total_positions % sample_every == 0:
            depths.append(int(parts[2]))

        # Once we have enough samples, increase sampling interval
        # to keep memory bounded while still getting a representative median
        if len(depths) == max_positions and sample_every == 1:
            # Estimate total genome size and set appropriate interval
            # At this point we've seen max_positions positions
            sample_every = max(1, total_positions // max_positions * 10)

    depth_proc.wait()

    if not depths:
        return 0.0, 0

    depths.sort()
    median = depths[len(depths) // 2]
    return float(median), len(depths)


def count_candidate_sites(site_stats: Path | None) -> int:
    """Count CANDIDATE sites from s05 stats file.

    Reads s05_stats.txt and counts lines with verdict=CANDIDATE.
    Returns 0 if no stats file provided or not found.
    """
    if site_stats is None or not site_stats.exists():
        return 0

    count = 0
    with open(site_stats) as fh:
        for line in fh:
            parts = line.strip().split("\t", 1)
            if len(parts) == 2 and parts[0].endswith("_verdict") and parts[1] == "CANDIDATE":
                count += 1
    return count


def classify_copy_number(depth_ratio: float) -> tuple[float, str]:
    """Classify the estimated copy number from depth ratio.

    Returns:
        (estimated_copies, description)
    """
    if depth_ratio < 0.25:
        return round(depth_ratio, 2), "low-level / mosaic"
    elif depth_ratio < 0.75:
        return 1.0, "hemizygous single-copy"
    elif depth_ratio < 1.25:
        return 1.0, "homozygous single-copy (or 2 hemizygous copies)"
    elif depth_ratio < 1.75:
        return 2.0, "likely 2-3 copies"
    else:
        return round(depth_ratio, 1), f"multi-copy (~{depth_ratio:.1f}x)"


def assess_site_validation(site_count: int, estimated_copies: float) -> str:
    """Cross-validate copy number estimate with candidate site count from s05.

    1 site = single insertion locus (consistent with single/tandem copy)
    >1 sites = possible multi-locus insertion
    0 sites = no candidate sites found

    Returns:
        Validation status string.
    """
    if site_count == 0:
        return "no_candidate_sites"
    elif site_count == 1:
        if estimated_copies <= 1.5:
            return "consistent_single_locus"
        else:
            return "depth_suggests_multicopy_single_locus"
    else:
        # >1 sites
        if estimated_copies > 1.5:
            return "consistent_multi_locus"
        else:
            return "multiple_sites_but_low_depth"


def determine_confidence(
    construct_median: float,
    host_median: float,
    site_count: int,
    site_validation: str,
) -> str:
    """Determine overall confidence in the copy number estimate."""
    if host_median < 1.0:
        return "Low"  # Insufficient host coverage for reliable ratio
    if construct_median < 1.0:
        return "Low"  # No construct coverage

    if site_count >= 1 and "consistent" in site_validation:
        return "High"
    elif site_count >= 1:
        return "Medium"
    else:
        return "Low"


def run_copynumber(
    construct_bam: Path,
    host_bam: Path,
    construct_ref: Path,
    host_ref: Path,
    outdir: Path,
    sample_name: str,
    site_stats: Path | None = None,
) -> None:
    """Estimate copy number and write results."""
    step_dir = outdir / sample_name / "s07_copynumber"
    step_dir.mkdir(parents=True, exist_ok=True)

    copynumber_tsv = step_dir / "copynumber.tsv"

    # Validate inputs
    for label, path in [
        ("Construct BAM", construct_bam),
        ("Host BAM", host_bam),
        ("Construct ref", construct_ref),
        ("Host ref", host_ref),
    ]:
        if not path.exists():
            log(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    log(f"Sample: {sample_name}")
    log(f"Construct BAM: {construct_bam}")
    log(f"Host BAM: {host_bam}")
    if site_stats:
        log(f"Site stats: {site_stats}")

    # Calculate construct median depth
    log("Calculating construct median depth...")
    construct_median, construct_positions, _ = get_construct_median_depth(construct_bam)
    log(f"  Construct median depth: {construct_median:.1f}x "
        f"({construct_positions:,} positions)")

    # Calculate host genome median depth (sample up to 1M positions)
    log("Calculating host genome median depth (sampling up to 1M positions)...")
    host_median, host_positions = get_host_median_depth(host_bam)
    log(f"  Host median depth: {host_median:.1f}x "
        f"({host_positions:,} positions sampled)")

    # Calculate depth ratio
    if host_median > 0:
        depth_ratio = construct_median / host_median
    else:
        log("WARNING: Host median depth is 0; cannot compute depth ratio")
        depth_ratio = 0.0

    log(f"  Depth ratio (construct/host): {depth_ratio:.3f}")

    # Classify copy number
    estimated_copies, copy_description = classify_copy_number(depth_ratio)
    log(f"  Estimated copies: {estimated_copies} ({copy_description})")

    # Count candidate sites and cross-validate
    site_count = count_candidate_sites(site_stats)
    log(f"  Candidate site count: {site_count}")

    site_validation = assess_site_validation(site_count, estimated_copies)
    log(f"  Site validation: {site_validation}")

    # Determine confidence
    confidence = determine_confidence(
        construct_median, host_median, site_count, site_validation,
    )
    log(f"  Confidence: {confidence}")

    # Write output TSV
    header = "\t".join([
        "sample",
        "construct_median_depth",
        "host_median_depth",
        "depth_ratio",
        "estimated_copies",
        "copy_description",
        "candidate_sites",
        "site_validation",
        "confidence",
    ])
    row = "\t".join([
        sample_name,
        f"{construct_median:.2f}",
        f"{host_median:.2f}",
        f"{depth_ratio:.4f}",
        f"{estimated_copies:.1f}",
        copy_description,
        str(site_count),
        site_validation,
        confidence,
    ])

    with open(copynumber_tsv, "w") as fh:
        fh.write(header + "\n")
        fh.write(row + "\n")

    log(f"Wrote copy number results to: {copynumber_tsv}")

    # Summary to stderr
    log("=== Copy Number Summary ===")
    log(f"  Construct median depth: {construct_median:.1f}x")
    log(f"  Host median depth:      {host_median:.1f}x")
    log(f"  Depth ratio:            {depth_ratio:.3f}")
    log(f"  Estimated copies:       {estimated_copies:.1f} ({copy_description})")
    log(f"  Candidate sites:        {site_count}")
    log(f"  Site validation:        {site_validation}")
    log(f"  Confidence:             {confidence}")
    log("Done.")


def main() -> None:
    args = parse_args()
    run_copynumber(
        args.construct_bam,
        args.host_bam,
        args.construct_ref,
        args.host_ref,
        args.outdir,
        args.sample_name,
        site_stats=args.site_stats,
    )


if __name__ == "__main__":
    main()
