#!/usr/bin/env python3
"""KCGP nomenclature (tentative_v0).

Placeholder ID format: KCGP-{HOST3}-{EVENT}-{YY}-{SITE4}.
Will be replaced with the official KCGP spec when delivered. The
mapping TSV schema is stable across spec versions; only the kcgp_id
values and spec_version tag change.
"""
from __future__ import annotations

import re
from pathlib import Path

HOST_CODES: dict[str, str] = {
    "Osativa_323_v7.0":   "OSA",   # rice
    "SLM_r2.0.pmol":      "SLE",   # tomato
    "CucSat_B10v3":       "CSA",   # cucumber
    "Zm_B73_v5":          "ZMA",   # corn
    "Gmax_v4.0":          "GMA",   # soybean
}

_KCGP_RE = re.compile(r"^KCGP-([A-Z]{3})-([A-Za-z0-9_]+)-(\d{2})-(\d{4})$")


def _normalize_host_path(host_ref_path: str) -> str:
    """Strip leading directory and a single trailing .fa/.fasta extension."""
    name = Path(host_ref_path).name
    for ext in (".fasta", ".fa"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def build_kcgp_id(
    host_ref_path: str,
    event: str,
    year: int,
    site_ord: int,
    spec: str = "tentative_v0",
) -> str:
    """Return the KCGP placeholder ID for one site.

    Raises KeyError if the host_ref normalized name is not registered in
    HOST_CODES. Raises ValueError if site_ord exceeds the 4-digit field.
    """
    norm = _normalize_host_path(host_ref_path)
    if norm not in HOST_CODES:
        raise KeyError(
            f"host '{norm}' (from {host_ref_path}) not in HOST_CODES; "
            f"known: {sorted(HOST_CODES)}"
        )
    if site_ord < 0 or site_ord > 9999:
        raise ValueError(f"site_ord {site_ord} out of 4-digit range (0-9999)")
    host3 = HOST_CODES[norm]
    yy = year % 100
    return f"KCGP-{host3}-{event}-{yy:02d}-{site_ord:04d}"


def parse_kcgp_id(kcgp_id: str) -> dict:
    """Parse a KCGP ID into its 4 fields. Raises ValueError on malformed input."""
    m = _KCGP_RE.match(kcgp_id)
    if not m:
        raise ValueError(f"not a KCGP-formatted ID: {kcgp_id!r}")
    host3, event, yy, site = m.groups()
    return {
        "host_code": host3,
        "event": event,
        "year_2digit": int(yy),
        "site_ord": int(site),
    }
