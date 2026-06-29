"""Shared transformation logic for the oncoprint pipeline.

Functions in this module are called by NB01 (data preparation) to convert
raw Excel values into the vocabulary defined in ``vocabulary.py``.
"""

from .vocabulary import (
    COMPLEX_BIOMARKER_VALUE_MAP,
    INTERESTING_FUSION_PARTNERS,
    POST_TRANSFORM_REPLACEMENTS,
)


def transform_value(value, modality):
    """Prepend modality prefix to a raw cell value.

    - ``"n.e."`` → ``"not assessed"``
    - ``"e."``   → ``"<modality> normal"``  (later fixed by post-transform)
    - anything else → ``"<modality> <value>"``
    """
    if value == "n.e.":
        return "not assessed"
    elif value == "e.":
        return f"{modality} normal"
    else:
        return f"{modality} {value}"


def apply_post_transform_replacements(df):
    """Apply the standard post-transform replacements in-place.

    Fixes combinations like ``"mut normal"`` → ``"mut assessed"`` using
    the rules from :data:`vocabulary.POST_TRANSFORM_REPLACEMENTS`.

    Returns the modified DataFrame (same object).
    """
    for old_val, new_val in POST_TRANSFORM_REPLACEMENTS.items():
        df = df.replace(old_val, new_val)
    return df


def apply_complex_biomarker_map(df):
    """Map raw complex-biomarker values to display values.

    Uses :data:`vocabulary.COMPLEX_BIOMARKER_VALUE_MAP`.

    Returns the modified DataFrame (same object).
    """
    for old_val, new_val in COMPLEX_BIOMARKER_VALUE_MAP.items():
        df = df.replace(old_val, new_val)
    return df


def filter_fusion_partners(df, cohort):
    """Collapse uninteresting fusion partner annotations.

    For cohorts listed in :data:`vocabulary.INTERESTING_FUSION_PARTNERS`,
    fusion cells that contain a ``"fus "`` prefix followed by a partner name
    **not** in the interesting set are replaced with ``"fus det."``.

    This keeps specific partner annotations (e.g. ``"fus FGFR3"``) for
    interesting partners while collapsing the rest.

    Parameters
    ----------
    df : pandas.DataFrame
        Gene-data matrix after ``transform_value`` + post-transform.
    cohort : str
        Cohort name (e.g. ``"CUP"``).

    Returns the modified DataFrame.
    """
    interesting = INTERESTING_FUSION_PARTNERS.get(cohort)
    if interesting is None:
        return df

    def _collapse(val):
        if not isinstance(val, str) or not val.startswith("fus "):
            return val
        # Values like "fus det.", "fus assessed" are already canonical
        suffix = val[4:]  # strip "fus " prefix
        if suffix in ("det.", "assessed"):
            return val
        # Check if the partner gene is interesting
        if suffix not in interesting:
            return "fus det."
        return val

    return df.map(_collapse)
