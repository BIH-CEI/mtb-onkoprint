# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # MTB Oncoprint — Visualization (CUP Cohort)
#
# Reads intermediate CSVs produced by `01_data_preparation_CUP.py` and
# generates an oncoprint plot with IHC heatmap.
#
# **Edit only the cells marked "EDIT FOR YOUR COHORT".**

# %% [markdown]
# ## Configuration — EDIT FOR YOUR COHORT

# %%
# === Cell 1: Study metadata (match data prep) ===  EDIT FOR YOUR COHORT

STUDY_DATE = "2026-02-09"
COHORT_NAME = "CUP"

# %%
# === Cell 2: Histology configuration ===  EDIT FOR YOUR COHORT

HISTOLOGY_ROW = "Histology"

HISTOLOGY_COLORS = {
    "ADC": "cornflowerblue",
    "SCC": "salmon",
    "Poorly differentiated carcinoma": "orange",
    "Undifferentiated carcinoma": "gold",
    "Sarcomatoid carcinoma": "orchid",
    "Carcinoma, not otherwise specified": "silver",
    "NEC": "seagreen",
}

# %%
# === Cell 3: Annotation configuration ===  EDIT FOR YOUR COHORT

SEX_ROW = "Sex (M/F)"
SEX_COLORS = {"M": "lightblue", "F": "lightpink"}

FAVOURABLE_RISK_ROW = "Favourable risk (y/n)"
FAVOURABLE_RISK_DISPLAY = {"yes": "yes", "no": "no"}
FAVOURABLE_RISK_COLORS = {"yes": "limegreen", "no": "white"}

DIAGNOSIS_ADAPTED_ROW = "Diagnosis adapted (y/n)"
DIAGNOSIS_ADAPTED_DISPLAY = {"y": "yes", "n": "no"}
DIAGNOSIS_ADAPTED_COLORS = {"yes": "dodgerblue", "no": "white"}

MMT_COUNT_ROW = "MMT recommended (number; na/not applicable)"

RECOMMENDATION_BAR_COLOR = "green"

# %%
# === Cell 4: Manual corrections ===  EDIT FOR YOUR COHORT

# Override values for specific patients if needed
HISTOLOGY_OVERRIDES = {}       # e.g. {"Patient 3": "ADC"}
RISK_OVERRIDES = {}            # e.g. {"Patient 5": "yes"}
MMT_COUNT_OVERRIDES = {}       # e.g. {"Patient 10": 2}

GENE_DISPLAY_OVERRIDES = {}    # e.g. {"ERBB2/Her2": "ERBB2/HER2"}

# %%
# === Cell 5: Plot settings ===  EDIT FOR YOUR COHORT

FIGSIZE = (18, 8)
GAP = 0.15
TOP_PLOT = False
RIGHT_PLOT = False

# %%
# === Cell 6: Marker definitions ===  EDIT FOR YOUR COHORT

from matplotlib.patches import Polygon

MARKERS = {
    # --- Base ---
    "not assessed": dict(marker="fill", color="lightgrey", zindex=0),
    "assessed":     dict(marker="fill", color="darkgrey", zindex=1),

    # --- Mutation band [0.8 - 1.0] ---
    "mut pathog/ likely pathog": dict(
        marker=Polygon([[0,0.8],[0,1.0],[1.0,1.0],[1,0.8]]),
        color="red", zindex=5),
    "mut VUS": dict(
        marker=Polygon([[0,0.8],[0,1.0],[1.0,1.0],[1,0.8]]),
        color="orange", zindex=4),
    "mut assessed": dict(
        marker=Polygon([[0,0.8],[0,1.0],[1.0,1.0],[1,0.8]]),
        color="darkgrey", zindex=1),

    # Gene-specific variant annotations
    "mut class I": dict(
        marker=Polygon([[0,0.8],[0,1.0],[1,1.0],[1,0.8],
                         [0.4,0.8],[0.4,0.55],[0.2,0.55],[0.2,0.8]]),
        color="red", zindex=5),
    "mut class II": dict(
        marker=Polygon([[0,0.8],[0,1.0],[1,1.0],[1,0.8],
                         [0.8,0.8],[0.8,0.55],[0.6,0.55],[0.6,0.8]]),
        color="red", zindex=5),
    "mut class III": dict(
        marker=Polygon([[0,0.8],[0,1.0],[1,1.0],[1,0.8],
                         [0.4,0.8],[0.4,0.55],[0.2,0.55],[0.2,0.8],
                         [0.8,0.8],[0.8,0.55],[0.6,0.55],[0.6,0.8]]),
        color="darkorange", zindex=5),
    "mut G12C": dict(
        marker=Polygon([[0,1.0],[1.0,1.0],[1,0.8]]),
        color="red", zindex=5),
    "mut nonG12C": dict(
        marker=Polygon([[0,0.8],[0,1.0],[1.0,1.0]]),
        color="red", zindex=5),

    # --- Fusion band [0.6 - 0.8] ---
    "fus det.":     dict(marker=Polygon([[0,0.6],[0,0.8],[1,0.8],[1,0.6]]),
                         color="yellow", zindex=5),
    "fus assessed": dict(marker=Polygon([[0,0.6],[0,0.8],[1,0.8],[1,0.6]]),
                         color="darkgrey", zindex=1),

    # --- CNV band [0.4 - 0.6] ---
    "CNV amp":      dict(marker=Polygon([[0,0.4],[0,0.6],[1,0.6],[1,0.4]]),
                         color="chartreuse", zindex=3),
    "CNV dup":      dict(marker=Polygon([[0,0.4],[0,0.6],[1,0.6],[1,0.4]]),
                         color="cyan", zindex=5),
    "CNV loss":     dict(marker=Polygon([[0,0.4],[0,0.6],[1,0.6],[1,0.4]]),
                         color="blue", zindex=2),
    "CNV assessed": dict(marker=Polygon([[0,0.4],[0,0.6],[1,0.6],[1,0.4]]),
                         color="darkgrey", zindex=1),

    # --- mRNA expression band [0.2 - 0.4] ---
    "mRNAexp overexp":  dict(marker=Polygon([[0,0.2],[0,0.4],[1,0.4],[1,0.2]]),
                              color="magenta", zindex=5),
    "mRNAexp underexp": dict(marker=Polygon([[0,0.2],[0,0.4],[1,0.4],[1,0.2]]),
                              color="olive", zindex=3),
    "mRNAexp assessed": dict(marker=Polygon([[0,0.2],[0,0.4],[1,0.4],[1,0.2]]),
                              color="darkgrey", zindex=1),

    # --- IHC band [0.0 - 0.2] ---
    "IHC exp":      dict(marker=Polygon([[0,0],[0,0.2],[1,0.2],[1,0]]),
                         color="black", zindex=5),
    "IHC neg.":     dict(marker=Polygon([[0,0],[0,0.2],[1,0.2],[1,0]]),
                         color="gold", zindex=3),
    "IHC assessed": dict(marker=Polygon([[0,0],[0,0.2],[1,0.2],[1,0]]),
                         color="darkgrey", zindex=1),

    # --- Complex biomarkers (full-cell fill) ---
    "high":         dict(marker="fill", color="red", zindex=5),
    "intermediate": dict(marker="fill", color="orange", zindex=5),
    "low":          dict(marker="fill", color="yellow", zindex=5),
    "positive":     dict(marker="fill", color="red", zindex=5),

    # --- HER2 IHC score (full-cell, 4 levels) ---
    "HER2 0":  dict(marker="fill", color="whitesmoke", zindex=1),
    "HER2 1+": dict(marker="fill", color="lightskyblue", zindex=3),
    "HER2 2+": dict(marker="fill", color="dodgerblue", zindex=4),
    "HER2 3+": dict(marker="fill", color="darkblue", zindex=5),
}

# %%
# === Cell 7: Legend configuration ===  EDIT FOR YOUR COHORT

LEGEND_ENABLED = True

LEGEND_GROUPS = [
    ("Mutation", [
        ("not assessed", "not assessed"),
        ("mut assessed", "assessed"),
        ("mut pathog/ likely pathog", "pathogenic / likely pathogenic"),
        ("mut class I", "class I (BRAF only)"),
        ("mut class III", "class III (BRAF only)"),
        ("mut G12C", "G12C (KRAS only)"),
        ("mut nonG12C", "non-G12C (KRAS only)"),
    ]),
    ("Fusion", [
        ("fus assessed", "assessed"),
        ("fus det.", "fusion"),
    ]),
    ("Copy Number Variation", [
        ("CNV assessed", "assessed"),
        ("CNV amp", "amplification"),
        ("CNV loss", "deletion"),
    ]),
    ("mRNA expression", [
        ("mRNAexp assessed", "assessed"),
        ("mRNAexp overexp", "overexpression"),
    ]),
    ("TMB", [
        ("assessed", "assessed"),
        ("high", "high"),
        ("intermediate", "intermediate"),
        ("low", "low"),
    ]),
    ("HER2 (IHC)", [
        ("not assessed", "not assessed"),
        ("HER2 0", "score 0"),
        ("HER2 1+", "score 1+"),
        ("HER2 2+", "score 2+"),
        ("HER2 3+", "score 3+"),
    ]),
    ("Other Biomarker", [
        ("assessed", "assessed"),
        ("positive", "positive (dMMR)"),
    ]),
]

LEGEND_COLUMNS = 4  # only used as fallback if LEGEND_LAYOUT is None
LEGEND_FONT_SIZE = 9
LEGEND_SWATCH_SIZE = 0.6
LEGEND_SCALER_STYLE = "stepped"

# Declarative legend layout: each inner list = one row of groups placed side by side.
# Group names reference LEGEND_GROUPS titles, annotation names, or heatmap names.
LEGEND_LAYOUT = [
    ["Mutation", "Fusion", "Copy Number Variation", "mRNA expression",
     "TMB", "HER2 (IHC)", "Other Biomarker"],
    ["Sex", "Histology", "Favourable Risk", "Rec.", "Dx adapted"],
    ["PD-L1", "ER", "PR", "AR"],
]

# %%
# === Cell 8: IHC heatmap configuration ===  EDIT FOR YOUR COHORT

# Which IHC markers to display as heatmap rows, and their colormaps
IHC_HEATMAP_CONFIG = {
    "PD-L1":  {"cmap": "YlOrRd", "vmin": 0, "vmax": 100},
    "ER":     {"cmap": "PuRd",   "vmin": 0, "vmax": 100},
    "PR":     {"cmap": "PuRd",   "vmin": 0, "vmax": 100},
    "AR":     {"cmap": "PuRd",   "vmin": 0, "vmax": 100},
    # HER2 now shown as 4-level biomarker (0/1+/2+/3+), not heatmap
    # MMR shown as dMMR complex biomarker
}

# %% [markdown]
# ---
# ## Processing Logic (do not edit below unless customizing)

# %%
import numpy as np
import pandas as pd
import pyoncoprint as pop
import matplotlib.pyplot as plt
import os
# %matplotlib inline

# %%
# --- Step 1: Load CSVs ---

clinical = pd.read_csv(os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_clinical.csv"),
                        index_col=0)
clinical.columns = clinical.columns.astype(str)

gene_data_single = pd.read_csv(os.path.join("output", f"{STUDY_DATE}_gene-data-single-rows.csv"),
                                index_col=0)
gene_data_single.columns = gene_data_single.columns.astype(str)

gene_data_aggregated = pd.read_csv(os.path.join("output", f"{STUDY_DATE}_gene-data-aggregated-rows.csv"),
                                    index_col=0)
gene_data_aggregated.columns = gene_data_aggregated.columns.astype(str)

ihc_heatmap = pd.read_csv(os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_ihc_heatmap.csv"),
                           index_col=0)
ihc_heatmap.columns = ihc_heatmap.columns.astype(str)

print(f"Loaded: clinical {clinical.shape}, gene_single {gene_data_single.shape}")
print(f"  gene_aggregated {gene_data_aggregated.shape}, ihc_heatmap {ihc_heatmap.shape}")

# %%
# --- Step 2: Apply manual corrections ---

for patient_id, val in HISTOLOGY_OVERRIDES.items():
    if patient_id in clinical.columns and HISTOLOGY_ROW in clinical.index:
        clinical.loc[HISTOLOGY_ROW, patient_id] = val

for old_name, new_name in GENE_DISPLAY_OVERRIDES.items():
    gene_data_single.index = pd.Index([s.replace(old_name, new_name) for s in gene_data_single.index])

# %%
# --- Step 2b: Validate data and MARKERS dict ---

from pyoncoprint.vocabulary import validate_marker_values, validate_markers_dict

combined_data = pd.concat([gene_data_single, gene_data_aggregated])
unknown_values = validate_marker_values(combined_data, strict=True)
if not unknown_values:
    print("All data values are valid.")

missing_markers = validate_markers_dict(MARKERS)
if missing_markers:
    raise ValueError(f"MARKERS dict is incomplete — missing: {missing_markers}")
else:
    print("MARKERS dict covers all required base values.")

# %%
# --- Step 3: Sort patients by histology frequency + favourable risk ---

combined_oncoprint = pd.concat([gene_data_single, gene_data_aggregated])

# Remove patients with no data at all
empty_columns = combined_oncoprint.columns[
    combined_oncoprint.apply(lambda x: all(x.isin(['not assessed', 'normal'])))]
print(f"Removing {len(empty_columns)} empty patients: {list(empty_columns)}")

needed_patients = [c for c in combined_oncoprint.columns if c not in empty_columns]

# Build sorting frame from clinical data
sort_df = pd.DataFrame(index=needed_patients)
if HISTOLOGY_ROW in clinical.index:
    sort_df["histology"] = clinical.loc[HISTOLOGY_ROW, needed_patients].values
if FAVOURABLE_RISK_ROW in clinical.index:
    sort_df["risk"] = clinical.loc[FAVOURABLE_RISK_ROW, needed_patients].values

# Rank by histology frequency
if "histology" in sort_df.columns:
    hist_freq = sort_df["histology"].value_counts().rank(method="first", ascending=True)
    sort_df["_freq_rank"] = sort_df["histology"].map(hist_freq)
else:
    sort_df["_freq_rank"] = 0

sort_df = sort_df.sort_values(by=["_freq_rank", "risk"], ascending=False)
sorted_patient_index = sort_df.index

print(f"\n{len(sorted_patient_index)} patients with data, sorted by histology frequency")

# %%
# --- Step 4: Build annotations ---

annotations = {}

# Sex annotation
if SEX_ROW in clinical.index:
    sex_annot = pd.DataFrame(
        clinical.loc[SEX_ROW, sorted_patient_index]).T
    annotations["Sex"] = {
        "annotations": sex_annot,
        "colors": SEX_COLORS,
        "order": 0,
    }

# Histology annotation
if HISTOLOGY_ROW in clinical.index:
    histology_annot = pd.DataFrame(
        clinical.loc[HISTOLOGY_ROW, sorted_patient_index]).T
    # Strip trailing whitespace from histology values
    histology_annot.loc[HISTOLOGY_ROW] = histology_annot.loc[HISTOLOGY_ROW].map(
        lambda x: x.strip() if isinstance(x, str) else x)
    # Sort legend entries by frequency (most common first)
    histology_freq_order = (
        histology_annot.loc[HISTOLOGY_ROW]
        .value_counts()
        .index.tolist()
    )
    annotations["Histology"] = {
        "annotations": histology_annot,
        "colors": HISTOLOGY_COLORS,
        "legend_order": histology_freq_order,
        "order": 1,
    }

# Favourable risk annotation
if FAVOURABLE_RISK_ROW in clinical.index:
    risk_annot = pd.DataFrame(
        clinical.loc[FAVOURABLE_RISK_ROW, sorted_patient_index]).T
    risk_annot.loc[FAVOURABLE_RISK_ROW] = risk_annot.loc[FAVOURABLE_RISK_ROW].map(
        FAVOURABLE_RISK_DISPLAY)
    annotations["Favourable Risk"] = {
        "annotations": risk_annot,
        "colors": FAVOURABLE_RISK_COLORS,
        "order": 2,
    }

# Recommendation count annotation (MMT count)
if MMT_COUNT_ROW in clinical.index:
    mmt_data = clinical.loc[MMT_COUNT_ROW, sorted_patient_index].copy()
    mmt_data = mmt_data.replace("na", "0").replace("not applicable", "0")
    mmt_data = pd.to_numeric(mmt_data, errors="coerce").fillna(0).astype(int)
    mmt_df = pd.DataFrame(mmt_data).T

    for patient_id, new_val in MMT_COUNT_OVERRIDES.items():
        if patient_id in mmt_df.columns:
            mmt_df.loc[MMT_COUNT_ROW, patient_id] = new_val

    annotations["Rec."] = {
        "annotations": mmt_df,
        "color": RECOMMENDATION_BAR_COLOR,
        "order": 3,
    }

# Diagnosis adapted annotation
if DIAGNOSIS_ADAPTED_ROW in clinical.index:
    dx_annot = pd.DataFrame(
        clinical.loc[DIAGNOSIS_ADAPTED_ROW, sorted_patient_index]).T
    dx_annot.loc[DIAGNOSIS_ADAPTED_ROW] = dx_annot.loc[DIAGNOSIS_ADAPTED_ROW].map(
        DIAGNOSIS_ADAPTED_DISPLAY)
    annotations["Dx adapted"] = {
        "annotations": dx_annot,
        "colors": DIAGNOSIS_ADAPTED_COLORS,
        "order": 4,
    }

print("Annotations built:")
for k, v in annotations.items():
    print(f"  {k}: {v['annotations'].shape}")

# %%
# --- Step 5: Build IHC heatmaps ---

heatmaps = {}

for marker, config in IHC_HEATMAP_CONFIG.items():
    if config is None:
        continue
    if marker not in ihc_heatmap.index:
        print(f"  WARNING: IHC marker '{marker}' not in heatmap data, skipping")
        continue
    hm_row = ihc_heatmap.loc[[marker], sorted_patient_index].astype(float)
    heatmaps[marker] = {
        "heatmap": hm_row,
        "cmap": config["cmap"],
        "vmin": config.get("vmin", hm_row.min().min()),
        "vmax": config.get("vmax", hm_row.max().max()),
    }

print(f"IHC heatmaps: {list(heatmaps.keys())}")

# %%
# --- Step 6: Create OncoPrint and render ---

plot_data = combined_oncoprint.drop(columns=empty_columns)[sorted_patient_index]
op = pop.OncoPrint(plot_data)

fig, axes = op.oncoprint(
    MARKERS,
    annotations=annotations,
    heatmaps=heatmaps,
    topplot=TOP_PLOT,
    rightplot=RIGHT_PLOT,
    figsize=FIGSIZE,
    gap=GAP,
    gene_sort_method="unsorted",
    sample_sort_method="unsorted",
    legend=LEGEND_ENABLED,
    legend_groups=LEGEND_GROUPS if LEGEND_ENABLED else None,
    legend_layout=LEGEND_LAYOUT,
    legend_columns=LEGEND_COLUMNS,
    legend_font_size=LEGEND_FONT_SIZE,
    legend_swatch_size=LEGEND_SWATCH_SIZE,
    legend_scaler_style=LEGEND_SCALER_STYLE,
    strict_markers=True,
)

# Hide the duplicate y-axis labels on the right
axes[1].set_axis_off()

plt.subplots_adjust(bottom=-0.5, top=1.03, hspace=0.5)
plt.show()

# %%
# --- Step 7: Save outputs ---

os.makedirs("output", exist_ok=True)

png_path = os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_oncoprint.png")
svg_path = os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_oncoprint.svg")

fig.savefig(png_path, bbox_inches="tight", transparent=False)
fig.savefig(svg_path, bbox_inches="tight", transparent=True)

print(f"Saved: {png_path}")
print(f"Saved: {svg_path}")
