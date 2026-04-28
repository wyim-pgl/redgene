"""RedGene s05 package (Issue #4 multi-session split).

Pure-function modules extracted from monolithic ``s05_insert_assembly.py``.
Session 1 (2026-04-19) added ``primitives``, ``filter_b_flank``, and
``filter_c_chimeric`` alongside the already-extracted ``verdict`` and
``config_loader``. Remaining decomposition (site_discovery, classify,
read_extraction, assembly, annotation, report, fanout_orchestrator) is
tracked in ``docs/superpowers/specs/2026-04-19-s05-module-split-design.md``.
Session 2 (2026-04-28) added ``filter_a_host`` and ``filter_d_altlocus``
filter modules.
Session 3 (2026-04-28) added ``annotation`` (Phase 4 BLAST-driven
element annotation), retiring the ``annotate_report.py`` shim's
re-export of ``annotate_insert``.
Session 4 (2026-04-28) added ``assembly`` (Phase 3 k-mer extension +
Pilon gap fill + minimap2/SSAKE refinement), the highest-risk extraction
in the s05 split. ``per_site.py`` shim updated to source
``assemble_insert`` / ``refine_with_foreign_reads`` from ``.assembly``.
"""
from .primitives import (  # noqa: F401
    STEP,
    log,
    revcomp,
    read_fasta,
    write_fasta,
    _read_fq_seqs,
    _parse_src_tag,
    JunctionCluster,
    InsertionSite,
    LegacyJunction,
    TierResult,
)
from .filter_a_host import (  # noqa: F401
    INSERT_HOST_MIN_PIDENT,
    _blast_insert_vs_host,
)
from .filter_b_flank import (  # noqa: F401
    CONSTRUCT_FLANK_PIDENT,
    CONSTRUCT_FLANK_MIN_LEN,
    CONSTRUCT_FLANK_SLOP,
    _find_construct_flanking_regions,
    _site_overlaps_flanking,
    filter_b_flanking_hit,
)
from .filter_c_chimeric import (  # noqa: F401
    CHIMERIC_MIN_PIDENT,
    CHIMERIC_MIN_OFFTARGET_BP,
    _check_chimeric_assembly,
)
from .filter_d_altlocus import (  # noqa: F401
    CONSTRUCT_HOST_MIN_COMBINED,
    CONSTRUCT_MIN_FRACTION,
    CONSTRUCT_HOST_MIN_PIDENT,
    _check_construct_host_coverage,
)
from .annotation import (  # noqa: F401
    _parse_blast6,
    _run_local_blast,
    _run_remote_blast,
    _merge_annotations,
    annotate_insert,
)
from .verdict import (  # noqa: F401
    compute_verdict,
    FilterEvidence,
    VerdictRules,
    _apply_canonical_override,
    REPORT_INTERESTING_VERDICTS,
)
from .config_loader import load_verdict_rules  # noqa: F401
from .assembly import (  # noqa: F401
    StrandAwareSeedExtender,
    assemble_insert,
    check_host_termination,
    extract_foreign_reads,
    pilon_fill,
    recruit_by_kmer,
    refine_with_foreign_reads,
)
