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
# # MTB Oncoprint — Data Preparation
#
# This notebook reads an MTB Excel file and produces intermediate CSVs for the visualization notebook.
#
# **Edit only the cells marked "EDIT FOR YOUR COHORT".**

# %% [markdown]
# ## Configuration — EDIT FOR YOUR COHORT

# %%
# === Cell 1: Study metadata ===  EDIT FOR YOUR COHORT

STUDY_DATE = "2025-01-27"
COHORT_NAME = "Uro"
INPUT_FILENAME = "2025-01-27-Uro_Tabelle_AK.xlsx"
PATIENT_ID_COLUMN = "Nr"
PERSONAL_COLUMNS = ["Geschlecht", "Geb"]

# %%
# === Cell 2: Columns to KEEP ===  EDIT FOR YOUR COHORT
# Everything not listed here and not detected as a gene column will be auto-dropped.

ENTITY_COLUMN = "Entität grob"
RARITY_COLUMN = "Selten ja/nein"

RELEVANT_THERAPY_COLUMNS = [
    "Alteration für Empf", "Nr Empfehlungen",
    "Therapieumsetzung", "erhaltene Therapie (von Empf)",
]

COMPLEX_BIOMARKERS = [
    "TML", "dMMR and/or MSI-h", "HRD Score",
    "AC3 BRCAness Signatur", "AC13 APOBEC Signatur",
]

METADATA_COLUMNS_TO_KEEP = (
    [ENTITY_COLUMN, RARITY_COLUMN]
    + RELEVANT_THERAPY_COLUMNS + COMPLEX_BIOMARKERS
)

# %%
# === Cell 3: Gene column detection mode ===  EDIT FOR YOUR COHORT

# How gene columns are structured in the Excel:
# "one_hot_multimodal" = each gene has up to 5 columns, one per modality
#     (e.g., EGFR mut, EGFR fus, EGFR CNV, EGFR IHC, EGFR mRNAexp)
#     -> auto-detected by modality suffix after column relabeling
# "single_column" = each gene has 1 column with combined/mixed info
#     -> must list gene columns explicitly in GENE_COLUMNS below
#     -> a parsing function splits values into 5 modality rows
GENE_COLUMN_MODE = "one_hot_multimodal"  # "one_hot_multimodal" or "single_column"

MODALITY_SUFFIXES = [" mut", " fus", " CNV", " IHC", " mRNAexp"]

# Only needed for "single_column" mode:
# List gene column names as they appear in the Excel (after whitespace stripping).
# GENE_COLUMNS = ["EGFR", "BRAF", "KRAS", "PTEN", ...]

# Only needed for "single_column" mode:
# Define how to parse a single cell value into modalities.
# This function receives a cell value string and returns a dict
# with keys from {mut, fus, CNV, IHC, mRNAexp} and parsed values.
# def PARSE_GENE_VALUE(value):
#     """Example parser - adapt to your data format."""
#     result = {mod: "n.e." for mod in ["mut", "fus", "CNV", "IHC", "mRNAexp"]}
#     if "amp" in str(value): result["CNV"] = "amp"
#     if "pathog" in str(value): result["mut"] = "pathog/ likely pathog"
#     return result

# %%
# === Cell 4: Column name standardization ===  EDIT FOR YOUR COHORT
# Applied BEFORE gene detection in one_hot_multimodal mode.
# Each rule is (old_substring, new_substring), applied sequentially.

COLUMN_RELABEL_RULES = [
    ("mut", " mut"), ("fus", " fus"), ("PD L1", "PDL1"),
    ("  ", " "), ("NRG 1/2/3", "NRG1/2/3"), (" exp", "exp"),
    ("IHC ", "IHC"), ("ER/ ESR1", "ER/ESR1"),
    ("FGFR2CNV", "FGFR2 CNV"),
    ("PDL1/ CD274 mut", "PDL1/CD274 mut"),
    ("PDL1/ CD274 fus", "PDL1/CD274 fus"),
    ("PDL1 / CD274 CNV", "PDL1/CD274 CNV"),
    ("PDL1/ CD274 IHC", "PDL1/CD274 IHC"),
    ("ERBB2/Her2 mRNA", "ERBB2/Her2 mRNAexp"),
    ("ERBB3/Her3 mRNA", "ERBB3/Her3 mRNAexp"),
    ("BRM1 mRNA", "BRM1 mRNAexp"),
]

# %%
# === Cell 5: Value cleaning rules ===  EDIT FOR YOUR COHORT
# Each rule is (old_value, new_value) applied via df.replace (exact match).

VALUE_CLEANING_RULES = [
    ("e", "e."), ("e. ", "e."), ("n.e", "n.e."), ("n.e.", "n.e."),
    ("Amplifikation", "amp"),
    ("pathog /likely pathog", "pathog/ likely pathog"),
    ("pathog./ likely pathog.", "pathog/ likely pathog"),
    ("pathog./ likely pathog. ", "pathog/ likely pathog"),
    ("exp (Score 1)", "exp"), ("exp (Score 3)", "exp"),
    ("amp ", "amp"), ("VUS ", "VUS"), ("n..e", "n.e."),
    ("n", "n.e."),
]

# %%
# === Cell 6: Complex biomarker display config ===  EDIT FOR YOUR COHORT

from pyoncoprint.vocabulary import COMPLEX_BIOMARKER_VALUE_MAP

COMPLEX_BIOMARKER_DISPLAY_NAMES = {
    "TML": "TMB", "dMMR and/or MSI-h": "dMMR +/- MSI-h",
    "AC3 BRCAness Signatur": "AC3", "AC13 APOBEC Signatur": "AC13",
    "HRD Score": "HRD Score",
}

GENE_FREQUENCY_THRESHOLD = 4

# %% [markdown]
# ---
# ## Processing Logic (do not edit below unless customizing)

# %%
import pandas as pd
import numpy as np
import os

# %%
# --- Step 1: Load Excel + whitespace strip ---

input_path = os.path.join("data", INPUT_FILENAME)
raw_df = pd.read_excel(input_path)

# Strip whitespace from column names
raw_df.columns = [str(c).strip() for c in raw_df.columns]

# Strip whitespace from string cell values
for col in raw_df.select_dtypes(include='object').columns:
    raw_df[col] = raw_df[col].map(lambda x: x.strip() if isinstance(x, str) else x)

print(f"Loaded {len(raw_df)} rows x {len(raw_df.columns)} columns from {INPUT_FILENAME}")

# %%
# --- Step 2: Drop personal columns, set index ---

missing_personal = [c for c in PERSONAL_COLUMNS if c not in raw_df.columns]
if missing_personal:
    print(f"WARNING: Personal columns not found (already removed?): {missing_personal}")

cols_to_drop = [c for c in PERSONAL_COLUMNS if c in raw_df.columns]
df = raw_df.drop(columns=cols_to_drop)

if PATIENT_ID_COLUMN not in df.columns:
    raise ValueError(f"Patient ID column '{PATIENT_ID_COLUMN}' not found. Available: {list(df.columns[:10])}...")

df = df.set_index(PATIENT_ID_COLUMN)
df.index = df.index.astype(str)
print(f"Index set to '{PATIENT_ID_COLUMN}', {len(df)} patients")


# %%
# --- Step 3: Apply column relabeling rules ---

def apply_relabel_rules(columns, rules):
    result = list(columns)
    for i, col in enumerate(result):
        for old, new in rules:
            col = col.replace(old, new)
        result[i] = col
    return result

original_columns = list(df.columns)
df.columns = apply_relabel_rules(df.columns, COLUMN_RELABEL_RULES)

renamed = [(o, n) for o, n in zip(original_columns, df.columns) if o != n]
if renamed:
    print(f"Relabeled {len(renamed)} columns:")
    for o, n in renamed[:20]:
        print(f"  '{o}' -> '{n}'")
    if len(renamed) > 20:
        print(f"  ... and {len(renamed) - 20} more")
else:
    print("No columns were relabeled.")

# %%
# --- Step 4: Detect gene columns ---

if GENE_COLUMN_MODE == "one_hot_multimodal":
    gene_columns = [c for c in df.columns
                    if any(c.endswith(suffix) for suffix in MODALITY_SUFFIXES)]
    print(f"Auto-detected {len(gene_columns)} gene/modality columns by suffix")

elif GENE_COLUMN_MODE == "single_column":
    gene_columns = [c for c in GENE_COLUMNS if c in df.columns]
    missing_genes = [c for c in GENE_COLUMNS if c not in df.columns]
    if missing_genes:
        print(f"WARNING: Gene columns not found: {missing_genes}")
    print(f"Using {len(gene_columns)} explicitly listed gene columns")

else:
    raise ValueError(f"Unknown GENE_COLUMN_MODE: {GENE_COLUMN_MODE}")

# %%
# --- Step 5: Validation report ---

# Check which metadata columns are present
found_metadata = [c for c in METADATA_COLUMNS_TO_KEEP if c in df.columns]
missing_metadata = [c for c in METADATA_COLUMNS_TO_KEEP if c not in df.columns]

if missing_metadata:
    print(f"WARNING: Expected metadata columns not found: {missing_metadata}")
print(f"Found {len(found_metadata)}/{len(METADATA_COLUMNS_TO_KEEP)} metadata columns")
print(f"Found {len(gene_columns)} gene columns")

# Determine what gets dropped
all_kept = set(found_metadata + gene_columns)
dropped = [c for c in df.columns if c not in all_kept]
print(f"\nAuto-dropping {len(dropped)} columns:")
for c in dropped:
    print(f"  - {c}")

# %%
# --- Step 6: Keep only gene + metadata columns ---

columns_to_keep = [c for c in df.columns if c in all_kept]
df = df[columns_to_keep]
print(f"Kept {len(df.columns)} columns ({len(gene_columns)} gene + {len(found_metadata)} metadata)")

# %%
# --- Step 7: Apply value cleaning rules, fill NaN ---

df = df.fillna("n.e.")

for old_val, new_val in VALUE_CLEANING_RULES:
    df = df.replace(old_val, new_val)

print("Value cleaning applied. Sample of unique values in gene columns:")
all_gene_values = set()
for col in gene_columns:
    if col in df.columns:
        all_gene_values.update(df[col].unique())
# Show non-standard values (not e. or n.e.)
non_standard = sorted([v for v in all_gene_values if v not in ('e.', 'n.e.')])
print(f"  Non-standard gene values: {non_standard}")

# %%
# --- Step 8: Export entities, complex biomarkers, therapy CSVs ---

os.makedirs("output", exist_ok=True)

# Entities
entity_cols = [c for c in [ENTITY_COLUMN, RARITY_COLUMN] if c in df.columns]
entities_df = df[entity_cols].T
entities_df.to_csv(os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_Entities.csv"))
print(f"Exported entities: {entities_df.shape}")

# Complex biomarkers
complex_cols = [c for c in COMPLEX_BIOMARKERS if c in df.columns]
complex_df = df[complex_cols].T
complex_df.to_csv(os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_complex.csv"))
print(f"Exported complex biomarkers: {complex_df.shape}")

# Therapy
therapy_cols = [c for c in RELEVANT_THERAPY_COLUMNS if c in df.columns]
therapy_df = df[therapy_cols].T
therapy_df.to_csv(os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_Therapy.csv"))
print(f"Exported therapy: {therapy_df.shape}")

# %%
# --- Step 9: Parse gene/modality, frequency threshold, transform values ---

# Build gene-only dataframe (transposed: genes as rows, patients as columns)
gene_df = df[gene_columns].T

if GENE_COLUMN_MODE == "one_hot_multimodal":
    # Extract gene name and modality from column names
    gene_df['_original_col'] = gene_df.index
    gene_df['_gene'] = gene_df['_original_col'].apply(
        lambda c: next((c[:c.rfind(s)] for s in MODALITY_SUFFIXES if c.endswith(s)), c))
    gene_df['_modality'] = gene_df['_original_col'].apply(
        lambda c: next((s.strip() for s in MODALITY_SUFFIXES if c.endswith(s)), ''))
    gene_df = gene_df.drop(columns='_original_col')

elif GENE_COLUMN_MODE == "single_column":
    # Expand each gene column into 5 modality rows using PARSE_GENE_VALUE
    expanded_rows = []
    modalities = ["mut", "fus", "CNV", "IHC", "mRNAexp"]
    for gene_col in gene_columns:
        for mod in modalities:
            row = {}
            row['_gene'] = gene_col
            row['_modality'] = mod
            for patient in df.index:
                val = df.loc[patient, gene_col]
                parsed = PARSE_GENE_VALUE(str(val))
                row[patient] = parsed.get(mod, 'n.e.')
            expanded_rows.append(row)
    gene_df = pd.DataFrame(expanded_rows)
    gene_df = gene_df.set_index(['_gene', '_modality'])
    gene_df = gene_df.reset_index()

# Count gene frequency (non-trivial entries across all modalities)
patient_cols = [c for c in gene_df.columns if c not in ('_gene', '_modality')]
gene_counts = gene_df.replace('n.e.', np.nan).replace('e.', np.nan)[patient_cols].groupby(
    gene_df['_gene']).count().sum(axis=1).sort_values()

single_genes = gene_counts[gene_counts >= GENE_FREQUENCY_THRESHOLD].index
aggregated_genes = gene_counts[gene_counts < GENE_FREQUENCY_THRESHOLD].index

print(f"Genes above threshold ({GENE_FREQUENCY_THRESHOLD}): {list(single_genes)}")
print(f"Genes below threshold (aggregated as 'Other'): {list(aggregated_genes)}")

# %%
# --- Step 9b: Transform values (prepend modality prefix) ---

from pyoncoprint.transforms import transform_value

for col in patient_cols:
    gene_df[col] = gene_df.apply(
        lambda row: transform_value(row[col], row['_modality']), axis=1)

# Set gene name as index, drop modality helper column
gene_df.index = gene_df['_gene']
gene_df.index.name = 'Gene/Biomarker'
gene_df = gene_df.drop(columns=['_gene', '_modality'])

# %%
# --- Step 10: Post-transform replacements ---

from pyoncoprint.transforms import apply_post_transform_replacements

gene_df = apply_post_transform_replacements(gene_df)

print("Post-transform replacements applied.")

# %%
# --- Step 10b: Filter fusion partners (cohort-specific) ---

from pyoncoprint.transforms import filter_fusion_partners

gene_df = filter_fusion_partners(gene_df, COHORT_NAME)

print(f"Fusion partner filtering applied for cohort '{COHORT_NAME}'.")

# %%
# --- Step 11: Process complex biomarkers ---

from pyoncoprint.transforms import apply_complex_biomarker_map

complex_op = pd.read_csv(os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_complex.csv"),
                         index_col=0, header=0)

# Apply value mapping (shared vocabulary)
complex_op = apply_complex_biomarker_map(complex_op)

# Rename biomarkers
rename_rows = {}
for orig_name, display_name in COMPLEX_BIOMARKER_DISPLAY_NAMES.items():
    if orig_name in complex_op.index and display_name != orig_name:
        rename_rows[orig_name] = display_name

complex_op = complex_op.rename(index=rename_rows)

# Drop any duplicates that arose from renaming (keep the renamed version)
complex_op = complex_op[~complex_op.index.duplicated(keep='last')]

print(f"Complex biomarkers processed: {list(complex_op.index)}")
complex_op

# %%
# --- Step 11b: Validate marker values before export ---

from pyoncoprint.vocabulary import validate_marker_values

gene_single = gene_df.loc[gene_df.index.isin(single_genes)]
combined_for_validation = pd.concat([gene_single, complex_op], axis=0)

unknown = validate_marker_values(combined_for_validation, strict=True)
if not unknown:
    print("All marker values are valid.")

# %%
# --- Step 12: Export gene-data CSVs ---

# Single-row genes (above threshold) + complex biomarkers
gene_single = gene_df.loc[gene_df.index.isin(single_genes)]
gene_data_single = pd.concat([gene_single, complex_op], axis=0)
gene_data_single.to_csv(os.path.join("output", f"{STUDY_DATE}_gene-data-single-rows.csv"))
print(f"Exported gene-data-single-rows: {gene_data_single.shape}")

# Aggregated genes (below threshold) -> labeled as 'Other'
gene_aggregated = gene_df.loc[gene_df.index.isin(aggregated_genes)]
gene_aggregated.index = pd.Index(['Other'] * len(gene_aggregated))
gene_aggregated.to_csv(os.path.join("output", f"{STUDY_DATE}_gene-data-aggregated-rows.csv"))
print(f"Exported gene-data-aggregated-rows: {gene_aggregated.shape}")

print("\nDone! CSVs are in the output/ directory.")
