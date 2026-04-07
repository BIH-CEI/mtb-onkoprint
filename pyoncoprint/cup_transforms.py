"""CUP-specific transformation logic for the oncoprint pipeline.

Maps raw CUP data values (variant strings, IHC scores, etc.) to the shared
vocabulary defined in ``vocabulary.py``.
"""

import re
import numpy as np
import pandas as pd

from .vocabulary import INTERESTING_FUSION_PARTNERS


# ---------------------------------------------------------------------------
# SNV classification
# ---------------------------------------------------------------------------

# BRAF variant → mutation class mapping
_BRAF_CLASS_MAP = {
    "V600E": "mut class I",
    "V600K": "mut class I",
    "V600D": "mut class I",
    "V600R": "mut class I",
    # Class III: kinase-impaired, activating via RAS-dependent mechanisms
    "D594G": "mut class III",
    "D594N": "mut class III",
    "G464V": "mut class III",
    "G466V": "mut class III",
    "T241M": "mut class III",  # per collaborator classification
}

# Regex to extract amino acid change from variant strings
# Matches patterns like: p.V600E, p.G12C, pG12D (no dot), p.H2417Qfs*3
_AA_CHANGE_RE = re.compile(r"p\.?([A-Z]\d+[A-Z*])")
# Also match short forms without "p." prefix: V600E, G12C
_AA_SHORT_RE = re.compile(r"\b([A-Z]\d+[A-Z])\b")
# Match AF annotations to strip: "AF 56%", "AF 11%"
_AF_RE = re.compile(r"\s*AF\s+[\d.]+%")


def _extract_aa_change(variant_str):
    """Extract amino acid change from a variant annotation string.

    Examples:
        'BRAF V600E' → 'V600E'
        'TP53 p.M246V' → 'M246V'
        'KRAS p.G12C' → 'G12C'
        'ATM p.L263fs AF 56%' → 'L263fs'
    """
    # Strip AF annotations first
    cleaned = _AF_RE.sub("", variant_str)
    # Try p.XXX format first
    match = _AA_CHANGE_RE.search(cleaned)
    if match:
        return match.group(1)
    # Try short format (V600E without p.)
    match = _AA_SHORT_RE.search(cleaned)
    if match:
        return match.group(1)
    return None


def classify_snv(value, gene):
    """Classify a CUP SNV value into the shared vocabulary.

    Parameters
    ----------
    value : str or float
        Raw SNV value from the transposed CSV (e.g. 'TP53 p.M246V',
        'BRAF V600E', 'KRAS p.G12C|KRAS p.G12D').
    gene : str
        Gene name (row label without '_SNV' suffix).

    Returns
    -------
    str
        Vocabulary token: 'mut pathog/ likely pathog', 'mut class I',
        'mut class III', 'mut G12C', 'mut nonG12C', or empty string
        if no value.
    """
    if pd.isna(value) or str(value).strip() == "":
        return ""

    value = str(value).strip()

    # Multi-hit: take first variant for classification (details in supplement)
    if "|" in value:
        value = value.split("|")[0].strip()

    aa_change = _extract_aa_change(value)

    # Gene-specific classification
    if gene == "BRAF" and aa_change:
        braf_class = _BRAF_CLASS_MAP.get(aa_change)
        if braf_class:
            return braf_class
        return "mut pathog/ likely pathog"

    if gene == "KRAS" and aa_change:
        if aa_change == "G12C":
            return "mut G12C"
        return "mut nonG12C"

    # Default: any non-empty variant = pathogenic/likely pathogenic
    return "mut pathog/ likely pathog"


# ---------------------------------------------------------------------------
# Deletion / Amplification / Signature classification
# ---------------------------------------------------------------------------

def classify_cnv(value, assessed_status):
    """Classify a deletion or amplification value.

    Parameters
    ----------
    value : str or float
        Raw value (gene name if detected, empty otherwise).
    assessed_status : str
        'assessed', 'not assessed', or '' for the patient in this modality.

    Returns
    -------
    str
        'CNV loss', 'CNV amp', 'CNV assessed', 'not assessed', or ''.
    """
    if pd.isna(value) or str(value).strip() == "":
        if assessed_status == "assessed":
            return "CNV assessed"
        elif assessed_status == "not assessed":
            return "not assessed"
        return ""
    return ""  # placeholder — caller decides loss vs amp


def classify_del(value, assessed_status):
    """Classify a deletion value."""
    if pd.isna(value) or str(value).strip() == "":
        if assessed_status == "assessed":
            return "CNV assessed"
        elif assessed_status == "not assessed":
            return "not assessed"
        return ""
    return "CNV loss"


def classify_amp(value, assessed_status):
    """Classify an amplification value."""
    if pd.isna(value) or str(value).strip() == "":
        if assessed_status == "assessed":
            return "CNV assessed"
        elif assessed_status == "not assessed":
            return "not assessed"
        return ""
    return "CNV amp"


def classify_sig(value, assessed_status):
    """Classify a signature/overexpression value."""
    if pd.isna(value) or str(value).strip() == "":
        if assessed_status == "not assessed":
            return "not assessed"
        return ""
    return "mRNAexp overexp"


def classify_her2(value):
    """Classify HER2 IHC score into vocabulary tokens.

    Parameters
    ----------
    value : str or float
        Raw IHC value (0, 1, 2, 2+, 3, or empty/NaN).

    Returns
    -------
    str
        'HER2 0', 'HER2 1+', 'HER2 2+', 'HER2 3+', or 'not assessed'.
    """
    if pd.isna(value) or str(value).strip() == "":
        return "not assessed"
    s = str(value).strip().rstrip("+")
    try:
        score = int(float(s))
    except (ValueError, TypeError):
        return "not assessed"
    if score == 0:
        return "HER2 0"
    elif score == 1:
        return "HER2 1+"
    elif score == 2:
        return "HER2 2+"
    elif score >= 3:
        return "HER2 3+"
    return "not assessed"


# ---------------------------------------------------------------------------
# TMB classification
# ---------------------------------------------------------------------------

def classify_tmb(value, high=10, intermediate=6):
    """Classify TMB into high/intermediate/low/not assessed.

    Parameters
    ----------
    value : str or float
        Raw TMB value (numeric or 'na'/empty).
    high : float
        Threshold for 'high' (>=).
    intermediate : float
        Threshold for 'intermediate' (>=).
    """
    if pd.isna(value) or str(value).strip() in ("", "na"):
        return "not assessed"
    try:
        tmb = float(value)
    except (ValueError, TypeError):
        return "not assessed"
    if tmb >= high:
        return "high"
    if tmb >= intermediate:
        return "intermediate"
    return "low"


# ---------------------------------------------------------------------------
# Fusion parsing
# ---------------------------------------------------------------------------

def parse_fusions(fusions_str, cohort="CUP"):
    """Parse fusion annotations and return interesting partner genes.

    Parameters
    ----------
    fusions_str : str
        Raw value from 'Fusions detected' row (e.g.
        'FGFR3::KIAA0203', 'FKBP15::RET Fusion (e23::e12)',
        'APC::ARB2A, ARID4B::DNM3', 'assessed', 'not assessed', 'na').
    cohort : str
        Cohort name for looking up interesting partners.

    Returns
    -------
    list of str
        List of interesting partner gene names that should get 'fus det.'
        entries. Empty list if no fusions detected.
    """
    if pd.isna(fusions_str):
        return []
    fusions_str = str(fusions_str).strip()
    if fusions_str in ("assessed", "not assessed", "na", ""):
        return []

    interesting = INTERESTING_FUSION_PARTNERS.get(cohort, set())
    partners = []

    # Handle multiple fusions separated by comma
    for part in fusions_str.split(","):
        part = part.strip()
        # Remove annotations like "Fusion (e23::e12)"
        # First extract gene pair from GENE1::GENE2 pattern
        match = re.match(r"(\w[\w/.-]*)::(\w[\w/.-]*)", part)
        if match:
            gene1, gene2 = match.group(1), match.group(2)
            if gene1 in interesting:
                partners.append(gene1)
            if gene2 in interesting:
                partners.append(gene2)

    return partners
