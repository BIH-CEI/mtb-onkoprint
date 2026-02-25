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
# # MTB Oncoprint — Visualization
#
# This notebook reads the intermediate CSVs produced by `01_data_preparation.ipynb` and generates an oncoprint plot.
#
# **Edit only the cells marked "EDIT FOR YOUR COHORT".**

# %% [markdown]
# ## Configuration — EDIT FOR YOUR COHORT

# %%
# === Cell 1: Study metadata (match notebook 1) ===  EDIT FOR YOUR COHORT

STUDY_DATE = "2025-01-27"
COHORT_NAME = "Uro"
ENTITY_COLUMN = "Entität grob"
RARITY_COLUMN = "Selten ja/nein"

# %%
# === Cell 2: Entity configuration ===  EDIT FOR YOUR COHORT

# Consolidate variant entity names to canonical forms
ENTITY_CONSOLIDATION = {
    "Prostatakarzinom, ADC": "Prostatakarzinom",
    "Prostakarzinom, ADC": "Prostatakarzinom",
    "Prostata, ADC": "Prostatakarzinom",
    "Prostata, NEC/MANEC/neuroendokrine Differenzierungskomponente": "Prostata, NEC/MANEC",
    "Nierenzellkarzinom, papillär": "Nierenzellkarzinom",
}

# Translate German entity names to English for the plot
ENTITY_TRANSLATIONS = {
    "Adenokarzinom der Harnblase": "Bladder cancer, adenocarcinoma",
    "Peniskarzinom": "Penile cancer",
    "Urachuskarzinom": "Urachal carcinoma",
    "andere": "Other tumor",
    "Nierenzellkarzinom": "Non-clear cell renal cell carcinoma",
    "Keimzelltumor": "Germ cell tumor / Testicular stromal tumor",
    "Urothelkarzinom": "Urothelial carcinoma",
    "nicht-invasiv, papillär Urothelkarzinom": "Noninvasive papillary urothelial carcinoma",
    "Blasenkarzinom,PECA": "Bladder cancer, PECA",
    "Blasenkarzinom, mind überwiegend neuroendokrine Differenzierung": "Bladder cancer, NEC or mixed NEC/UC",
    "Prostata, NEC/MANEC": "Prostate cancer, NEC/MANEC",
    "Prostatakarzinom": "Prostate cancer, adenocarcinoma",
    "Urethrakarzinom": "Urethral cancer",
    "Prostata, Basalzell-Adenokarzinom": "Prostate cancer, basal cell adenocarcinoma",
}

# Colors for each translated entity name
ENTITY_COLORS = {
    "Prostate cancer, adenocarcinoma": "red",
    "Prostate cancer, basal cell adenocarcinoma": "coral",
    "Prostate cancer, NEC/MANEC": "salmon",
    "Urothelial carcinoma": "olive",
    "Noninvasive papillary urothelial carcinoma": "mediumseagreen",
    "Penile cancer": "blue",
    "Urachal carcinoma": "violet",
    "Non-clear cell renal cell carcinoma": "purple",
    "Germ cell tumor / Testicular stromal tumor": "darkgreen",
    "Urethral cancer": "yellow",
    "Bladder cancer, adenocarcinoma": "gold",
    "Bladder cancer, NEC or mixed NEC/UC": "yellowgreen",
    "Bladder cancer, PECA": "teal",
    "Other tumor": "darkgrey",
}

# %%
# === Cell 3: Manual corrections ===  EDIT FOR YOUR COHORT

# Override rarity for specific patients (patient_id -> new rarity value)
ENTITY_RARITY_OVERRIDES = {"45": "ja", "47": "ja", "50": "ja", "53": "ja"}

# Override recommendation counts for specific patients
RECOMMENDATION_COUNT_OVERRIDES = {"33": 2}

# Override therapy implementation counts for specific patients
THERAPY_COUNT_OVERRIDES = {"2a": 1, "2b": 1}

# %%
# === Cell 4: Annotation cleaning maps ===  EDIT FOR YOUR COHORT

RECOMMENDATION_VALUE_REPLACEMENTS = {
    "zuvor verstorben": "1", "zuvor verstorben ": "1", "n.e.": "0",
    "2 (insg 7 für Ptn)": "2", "2 (insg 4 für Ptn)": "2",
    "malignes epitheloides Angiomyolipom": "0",
}

THERAPY_IMPLEMENTATION_VALUE_MAP = {"n.e.": "0", "ja": 1, "ja,2": 2}

RARITY_DISPLAY_MAP = {"ja": "yes", "nein": "no"}
RARITY_COLORS = {"yes": "red", "no": "white"}

RECOMMENDATION_BAR_COLOR = "green"
MMT_BAR_COLOR = "red"

GENE_DISPLAY_OVERRIDES = {"ERBB2/Her2": "ERBB2/HER2"}

# %%
# === Cell 5: Plot settings ===  EDIT FOR YOUR COHORT

FIGSIZE = (16, 6)
GAP = 0.2
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
                         [0.4,0.8],[0.4,0.7],[0.25,0.7],[0.25,0.8]]),
        color="red", zindex=5),
    "mut class II": dict(
        marker=Polygon([[0,0.8],[0,1.0],[1,1.0],[1,0.8],
                         [0.75,0.8],[0.75,0.7],[0.6,0.7],[0.6,0.8]]),
        color="red", zindex=5),
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
}

# HOW TO ADD GENE-SPECIFIC VARIANT ANNOTATIONS:
# 1. In the Excel, put the annotation value in the gene's modality column
#    (e.g., "exon 19 del" in EGFR mut column)
# 2. Make sure VALUE_CLEANING_RULES doesn't clean it away
# 3. After transformation it becomes "mut exon 19 del"
# 4. Add a MARKERS entry: "mut exon 19 del": dict(marker=Polygon(...), ...)
# 5. Add it to LEGEND_GROUPS so it appears in the legend
# You can use any Polygon shape within the band boundaries.

# %%
# === Cell 7: Legend configuration ===  EDIT FOR YOUR COHORT

LEGEND_ENABLED = True

# Group legend entries by modality with human-readable labels.
# Each group: (group_title, [(marker_key, display_label), ...])
# Groups are arranged in columns. Within each group, entries are listed vertically.
LEGEND_GROUPS = [
    ("Mutation", [
        ("not assessed", "not assessed"),
        ("mut assessed", "assessed"),
        ("mut VUS", "variant of unclear significance"),
        ("mut pathog/ likely pathog", "pathogenic / likely pathogenic"),
        ("mut class I", "class I (BRAF only)"),
        ("mut class II", "class II (BRAF only)"),
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
        ("CNV dup", "duplication"),
        ("CNV loss", "deletion"),
    ]),
    ("mRNA expression", [
        ("mRNAexp assessed", "assessed"),
        ("mRNAexp overexp", "overexpression"),
        ("mRNAexp underexp", "underexpression"),
    ]),
    ("Immunohistochemistry", [
        ("IHC assessed", "assessed"),
        ("IHC exp", "expression"),
        ("IHC neg.", "negative"),
    ]),
    ("TMB", [
        ("assessed", "assessed"),
        ("high", "high"),
        ("intermediate", "intermediate"),
        ("low", "low"),
    ]),
    ("Other Biomarker", [
        ("assessed", "assessed"),
        ("positive", "mutational signature"),
    ]),
]

# Legend layout
LEGEND_COLUMNS = 4          # max groups per row before wrapping
LEGEND_FONT_SIZE = 9        # text size
LEGEND_SWATCH_SIZE = 0.6    # size of color swatches

# Numeric annotation scaler style in legend (for Rec. and MMT bars)
# "slope"   = triangular fill from min to max
# "stepped" = discrete stepped bars (one per integer value, increasing height)
LEGEND_SCALER_STYLE = "stepped"

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

entities = pd.read_csv(os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_Entities.csv"),
                       index_col=0).T
entities.columns = entities.columns.astype(str)

therapy = pd.read_csv(os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_Therapy.csv"),
                      index_col=0)
therapy.columns = therapy.columns.astype(str)

gene_data_single = pd.read_csv(os.path.join("output", f"{STUDY_DATE}_gene-data-single-rows.csv"),
                               index_col=0)
gene_data_single.columns = gene_data_single.columns.astype(str)

gene_data_aggregated = pd.read_csv(os.path.join("output", f"{STUDY_DATE}_gene-data-aggregated-rows.csv"),
                                   index_col=0)
gene_data_aggregated.columns = gene_data_aggregated.columns.astype(str)

print(f"Loaded: entities {entities.shape}, therapy {therapy.shape}")
print(f"  gene_data_single {gene_data_single.shape}, gene_data_aggregated {gene_data_aggregated.shape}")

# %%
# --- Step 2: Apply manual corrections ---

for patient_id, rarity_val in ENTITY_RARITY_OVERRIDES.items():
    if patient_id in entities.index:
        entities.loc[patient_id, RARITY_COLUMN] = rarity_val
        print(f"  Rarity override: patient {patient_id} -> {rarity_val}")

# Apply gene display overrides
for old_name, new_name in GENE_DISPLAY_OVERRIDES.items():
    gene_data_single.index = pd.Index([s.replace(old_name, new_name) for s in gene_data_single.index])

# %%
# --- Step 2b: Validate data and MARKERS dict ---

from pyoncoprint.vocabulary import validate_marker_values, validate_markers_dict

# Validate that all values in the data have a matching MARKERS entry
combined_data = pd.concat([gene_data_single, gene_data_aggregated])
unknown_values = validate_marker_values(combined_data, strict=True)
if not unknown_values:
    print("All data values are valid.")

# Validate that MARKERS dict covers all base marker values
missing_markers = validate_markers_dict(MARKERS)
if missing_markers:
    raise ValueError(f"MARKERS dict is incomplete — missing: {missing_markers}")
else:
    print("MARKERS dict covers all required base values.")

# %%
# --- Step 3: Consolidate entities, sort patients ---

for old_name, new_name in ENTITY_CONSOLIDATION.items():
    entities[ENTITY_COLUMN] = entities[ENTITY_COLUMN].replace(old_name, new_name)

entity_freq = entities[ENTITY_COLUMN].value_counts()
print(f"Entity frequencies after consolidation:")
print(entity_freq)
print(f"\n{len(entity_freq)} unique entities")

# %%
# --- Step 4: Remove empty patients, sort by entity frequency + rarity ---

combined_oncoprint = pd.concat([gene_data_single, gene_data_aggregated])

# Find patients with no data at all
empty_columns = combined_oncoprint.columns[
    combined_oncoprint.apply(lambda x: all(x.isin(['not assessed', 'normal'])))]
print(f"Removing {len(empty_columns)} empty patients (no data): {list(empty_columns)}")

# Filter entities and data to non-empty patients
needed_patients = [c for c in combined_oncoprint.columns if c not in empty_columns]

# Sort by entity frequency (descending) then rarity
patient_entities = entities.loc[entities.index.isin(needed_patients)].copy()
patient_entity_freq = patient_entities[ENTITY_COLUMN].value_counts().rank(method='first', ascending=True)
patient_entities['_freq_rank'] = patient_entities[ENTITY_COLUMN].map(patient_entity_freq)
sorted_patients = patient_entities.sort_values(
    by=['_freq_rank', RARITY_COLUMN], ascending=False)
sorted_patient_index = sorted_patients.index

print(f"\n{len(sorted_patient_index)} patients with data, sorted by entity frequency + rarity")

# %%
# --- Step 5: Build annotations ---

annotations = {}

# Entity annotation
entity_annot = pd.DataFrame(
    entities.loc[sorted_patient_index, ENTITY_COLUMN]).T
entity_annot.loc[ENTITY_COLUMN] = entity_annot.loc[ENTITY_COLUMN].map(ENTITY_TRANSLATIONS)

annotations['Entity'] = {
    'annotations': entity_annot,
    'colors': ENTITY_COLORS,
    'order': 0,
}

# Rarity annotation
rarity_annot = pd.DataFrame(
    entities.loc[sorted_patient_index, RARITY_COLUMN]).T
rarity_annot.loc[RARITY_COLUMN] = rarity_annot.loc[RARITY_COLUMN].map(RARITY_DISPLAY_MAP)

annotations['Rare Tumor'] = {
    'annotations': rarity_annot,
    'colors': RARITY_COLORS,
    'order': 1,
}

# Recommendation count annotation
rec_data = therapy.loc['Nr Empfehlungen', sorted_patient_index].copy()
for old_val, new_val in RECOMMENDATION_VALUE_REPLACEMENTS.items():
    rec_data = rec_data.replace(old_val, new_val)
rec_df = pd.DataFrame(rec_data.astype(int)).T

for patient_id, new_val in RECOMMENDATION_COUNT_OVERRIDES.items():
    if patient_id in rec_df.columns:
        rec_df.loc['Nr Empfehlungen', patient_id] = new_val

annotations['Rec.'] = {
    'annotations': rec_df,
    'color': RECOMMENDATION_BAR_COLOR,
    'order': 2,
}

# MMT (Molecular Matched Therapy) annotation
mmt_data = therapy.loc['Therapieumsetzung', sorted_patient_index].copy()
for old_val, new_val in THERAPY_IMPLEMENTATION_VALUE_MAP.items():
    mmt_data = mmt_data.replace(old_val, new_val)
mmt_df = pd.DataFrame(mmt_data.astype(int)).T

for patient_id, new_val in THERAPY_COUNT_OVERRIDES.items():
    if patient_id in mmt_df.columns:
        mmt_df.loc['Therapieumsetzung', patient_id] = new_val

annotations['MMT'] = {
    'annotations': mmt_df,
    'color': MMT_BAR_COLOR,
    'order': 3,
}

print("Annotations built:")
for k, v in annotations.items():
    print(f"  {k}: {v['annotations'].shape}")

# %%
# --- Step 6: Create OncoPrint and render ---

plot_data = combined_oncoprint.drop(columns=empty_columns)[sorted_patient_index]
op = pop.OncoPrint(plot_data)

fig, axes = op.oncoprint(
    MARKERS,
    annotations=annotations,
    topplot=TOP_PLOT,
    rightplot=RIGHT_PLOT,
    figsize=FIGSIZE,
    gap=GAP,
    gene_sort_method='unsorted',
    sample_sort_method='unsorted',
    legend=LEGEND_ENABLED,
    legend_groups=LEGEND_GROUPS if LEGEND_ENABLED else None,
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

fig.savefig(png_path, bbox_inches='tight', transparent=False)
fig.savefig(svg_path, bbox_inches='tight', transparent=True)

print(f"Saved: {png_path}")
print(f"Saved: {svg_path}")
