"""RedGene s05 package (v1.0 MVP minimal split).

Pure-function modules extracted from monolithic s05_insert_assembly.py.
v1.1 will continue the decomposition; v1.0 stops at compute_verdict + config.
"""
from .verdict import compute_verdict, FilterEvidence, VerdictRules  # noqa: F401
from .config_loader import load_verdict_rules  # noqa: F401
