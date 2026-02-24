"""Single source of truth for all valid marker values in the oncoprint pipeline.

This module defines the shared vocabulary between data preparation (NB01)
and visualization (NB02). Any value that appears in the data must be listed
here, otherwise validation will flag it.
"""

import warnings

# ---------------------------------------------------------------------------
# Base marker values — recognised by the visualisation layer
# ---------------------------------------------------------------------------

BASE_MARKER_VALUES = {
    # Meta
    "not assessed",
    "assessed",

    # Mutation
    "mut pathog/ likely pathog",
    "mut VUS",
    "mut assessed",

    # Fusion
    "fus det.",
    "fus assessed",

    # CNV
    "CNV amp",
    "CNV dup",
    "CNV loss",
    "CNV assessed",

    # mRNA expression
    "mRNAexp overexp",
    "mRNAexp underexp",
    "mRNAexp assessed",

    # IHC
    "IHC exp",
    "IHC neg.",
    "IHC assessed",

    # Complex biomarkers (full-cell)
    "high",
    "intermediate",
    "low",
    "positive",
}

# ---------------------------------------------------------------------------
# Gene-specific marker values — cohort-dependent extensions
# ---------------------------------------------------------------------------

GENE_SPECIFIC_MARKERS = {
    # BRAF variant classes
    "mut class I",
    "mut class II",
    "mut class III",

    # KRAS variant annotations
    "mut G12C",
    "mut nonG12C",
}

# All valid marker values (union)
ALL_MARKER_VALUES = BASE_MARKER_VALUES | GENE_SPECIFIC_MARKERS

# ---------------------------------------------------------------------------
# Post-transform replacements (NB01 Step 10)
# ---------------------------------------------------------------------------
# After transform_value() prepends the modality prefix, these fix
# nonsensical combinations like "mut normal" -> "mut assessed".

POST_TRANSFORM_REPLACEMENTS = {
    "mut normal": "mut assessed",
    "fus fus": "fus det.",
    "fus normal": "fus assessed",
    "CNV normal": "CNV assessed",
    "CNV del": "CNV loss",
    "IHC normal": "IHC assessed",
    "mRNAexp normal": "mRNAexp assessed",
}

# ---------------------------------------------------------------------------
# Complex biomarker value map (NB01 Step 11)
# ---------------------------------------------------------------------------
# Maps raw Excel values to display values for complex biomarkers
# (TMB, dMMR/MSI-h, HRD, signatures).

COMPLEX_BIOMARKER_VALUE_MAP = {
    "det": "positive",
    "det.": "positive",
    "yes": "positive",
    "n.e.": "not assessed",
    "e.": "assessed",
    "mutational signature det.": "positive",
}

# ---------------------------------------------------------------------------
# Interesting fusion partners per cohort
# ---------------------------------------------------------------------------
# For CUP: certain fusion partners are shown individually, the rest are
# collapsed into "fus det." on the partner gene row.

INTERESTING_FUSION_PARTNERS = {
    "CUP": {"FGFR3", "FGFR2", "RET", "EWSR1", "APC", "ARID4B"},
}

# ---------------------------------------------------------------------------
# Validation functions
# ---------------------------------------------------------------------------


def validate_marker_values(df, known_values=None, strict=False):
    """Check all cell values in *df* against the known vocabulary.

    Parameters
    ----------
    df : pandas.DataFrame
        Gene-data matrix (genes × patients) with string marker values.
    known_values : set, optional
        Accepted values. Defaults to ``ALL_MARKER_VALUES``.
    strict : bool
        If *True*, raise ``ValueError`` when unknown values are found.

    Returns
    -------
    set
        Unknown values (empty if everything is valid).
    """
    if known_values is None:
        known_values = ALL_MARKER_VALUES

    unique_values = set()
    for col in df.columns:
        unique_values.update(df[col].unique())

    unknown = unique_values - known_values
    if unknown:
        msg = (
            f"Unknown marker values found ({len(unknown)}): "
            + ", ".join(repr(v) for v in sorted(unknown))
        )
        if strict:
            raise ValueError(msg)
        warnings.warn(msg, stacklevel=2)

    return unknown


def validate_markers_dict(markers_dict, expected_values=None):
    """Check that *markers_dict* has an entry for every expected value.

    Parameters
    ----------
    markers_dict : dict
        The MARKERS dictionary mapping value -> style definition.
    expected_values : set, optional
        Values that must have a style entry. Defaults to ``BASE_MARKER_VALUES``.

    Returns
    -------
    set
        Missing keys (empty if the dict is complete).
    """
    if expected_values is None:
        expected_values = BASE_MARKER_VALUES

    missing = expected_values - set(markers_dict.keys())
    if missing:
        msg = (
            f"MARKERS dict is missing entries for: "
            + ", ".join(repr(v) for v in sorted(missing))
        )
        warnings.warn(msg, stacklevel=2)

    return missing
