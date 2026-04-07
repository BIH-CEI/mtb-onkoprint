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
# # MTB Oncoprint — Data Preparation (CUP Cohort)
#
# Reads the pre-transposed CUP CSV and produces intermediate CSVs for the
# visualization notebook.
#
# **Edit only the cells marked "EDIT FOR YOUR COHORT".**

# %% [markdown]
# ## Configuration — EDIT FOR YOUR COHORT

# %%
# === Cell 1: Study metadata ===  EDIT FOR YOUR COHORT

STUDY_DATE = "2026-02-09"
COHORT_NAME = "CUP"
INPUT_FILENAME = "Oncoprin_CUPxMTB_20260209_editHL_transposed.csv"

# %%
# === Cell 2: Clinical metadata rows ===  EDIT FOR YOUR COHORT
# Row labels in the transposed CSV to extract as clinical metadata.

CLINICAL_ROWS = {
    "histology": "Histology",
    "favourable_risk": "Favourable risk (y/n)",
    "mmt_recommended": "MMT recommended (y/n)",
    "mmt_count": "MMT recommended (number; na/not applicable)",
    "diagnosis_adapted": "Diagnosis adapted (y/n)",
    "sex": "Sex (M/F)",
}

# Rows containing free-text summaries (kept as reference, not used in plot)
REFERENCE_ROWS = [
    "IHC assessed (PD-L1, HER2, MMR, TROP2, AR, ER, PR, na)",
    "IHC scores",
    "Fusions detected",
    "Relevant SNVs",
    "All SNVs",
    "Relevant Deletions",
    "Relevant Amplifications",
    "Overexpression/Signatures",
    "Methylation analysis (column not yet final!)",
]

# %%
# === Cell 3: Modality configuration ===  EDIT FOR YOUR COHORT

MODALITY_SUFFIXES = {
    "_SNV": "mut",   # maps to mutation band
    "_Del": "CNV",   # maps to CNV band (loss)
    "_Amp": "CNV",   # maps to CNV band (amp)
    "_Sig": "Sig",   # maps to complex biomarker (positive)
}

# Assessment tracking rows per modality
ASSESSMENT_TRACKING = {
    "_Del": {"assessed": "assessed_Del", "not_assessed": "not_assessed_Del"},
    "_Amp": {"assessed": "assessed_Amp", "not_assessed": "not_assessed_Amp"},
    "_Sig": {"not_assessed": "not_assessed_Sig"},
}

# IHC rows are handled separately as numeric heatmap
IHC_SUFFIX = "_IHC"

# %%
# === Cell 4: Data cleaning rules ===  EDIT FOR YOUR COHORT

# (row_label, action, *args)
# "drop"   = remove row entirely
# "rename" = rename row label
# "merge"  = merge value into target row then drop source
CLEANING_RULES = [
    ("NTRK3,_Sig", "rename", "NTRK3_Sig"),
]

# Manual value fixes: (row_label, patient, new_value)
VALUE_FIXES = [
    # Patient 93 KRAS G12C: only in "All SNVs", not in parsed rows
    ("KRAS_SNV", "Patient 93", "KRAS p.G12C"),
]

# %%
# === Cell 5: Gene-specific SNV rules ===  EDIT FOR YOUR COHORT
# See pyoncoprint/cup_transforms.py for the BRAF class map and KRAS G12C logic.
# Additional BRAF class mappings can be added here if needed.

BRAF_EXTRA_CLASS_MAP = {}  # e.g. {"T241P": "mut class III"}

# %%
# === Cell 6: TMB thresholds ===  EDIT FOR YOUR COHORT

TMB_HIGH_THRESHOLD = 10
TMB_INTERMEDIATE_THRESHOLD = 6

# %%
# === Cell 7: Frequency threshold ===  EDIT FOR YOUR COHORT

GENE_FREQUENCY_THRESHOLD = 3

# %%
# === Cell 8: Fusion partners ===  EDIT FOR YOUR COHORT
# Which fusion partner genes are "interesting" and get their own row.
# Others are collapsed. This is already configured in vocabulary.py
# for CUP but can be overridden here.

# from pyoncoprint.vocabulary import INTERESTING_FUSION_PARTNERS
# FUSION_PARTNERS = INTERESTING_FUSION_PARTNERS.get("CUP", set())

# %% [markdown]
# ---
# ## Processing Logic (do not edit below unless customizing)

# %%
import pandas as pd
import numpy as np
import os

from pyoncoprint.cup_transforms import (
    classify_snv, classify_del, classify_amp, classify_sig,
    classify_tmb, classify_her2, parse_fusions,
)
from pyoncoprint.vocabulary import validate_marker_values

# %%
# --- Step 1: Load and clean ---

input_path = os.path.join("data", INPUT_FILENAME)
df = pd.read_csv(input_path, index_col=0)

# Strip whitespace from all string values
for col in df.columns:
    df[col] = df[col].map(lambda x: x.strip() if isinstance(x, str) else x)

# Strip whitespace from index
df.index = pd.Index([str(idx).strip() for idx in df.index])

patients = df.columns.tolist()
print(f"Loaded {len(df)} rows x {len(patients)} patients from {INPUT_FILENAME}")

# %%
# --- Step 2: Apply cleaning rules ---

for rule in CLEANING_RULES:
    row_label, action = rule[0], rule[1]
    if row_label not in df.index:
        print(f"  WARNING: Row '{row_label}' not found, skipping rule")
        continue
    if action == "drop":
        df = df.drop(index=row_label)
        print(f"  Dropped row '{row_label}'")
    elif action == "rename":
        new_label = rule[2]
        df = df.rename(index={row_label: new_label})
        print(f"  Renamed '{row_label}' -> '{new_label}'")
    elif action == "merge":
        target_row = rule[2]
        if target_row in df.index:
            for col in df.columns:
                src = df.loc[row_label, col]
                tgt = df.loc[target_row, col]
                if pd.notna(src) and str(src).strip() != "":
                    if pd.isna(tgt) or str(tgt).strip() == "":
                        df.loc[target_row, col] = src
        df = df.drop(index=row_label)
        print(f"  Merged '{row_label}' into '{target_row}' and dropped")

# Apply manual value fixes
for row_label, patient, new_value in VALUE_FIXES:
    if row_label in df.index and patient in df.columns:
        df.loc[row_label, patient] = new_value
        print(f"  Fixed {row_label} [{patient}] = '{new_value}'")

# %%
# --- Step 3: Separate clinical, molecular, and IHC rows ---

clinical_labels = list(CLINICAL_ROWS.values())
reference_labels = REFERENCE_ROWS

# Clinical metadata
clinical_df = df.loc[df.index.isin(clinical_labels)].copy()

# TMB row (special handling)
tmb_row = df.loc["TMB"] if "TMB" in df.index else None

# Fusions row
fusions_row = df.loc["Fusions detected"] if "Fusions detected" in df.index else None

# MMR from IHC (extract before IHC processing for complex biomarkers)
mmr_row = df.loc["MMR_IHC"] if "MMR_IHC" in df.index else None

# Determine which rows are molecular features vs metadata
all_non_molecular = set(clinical_labels + reference_labels + ["TMB", "Fusions detected"])
molecular_rows = [r for r in df.index if r not in all_non_molecular]

# Separate IHC rows
ihc_rows = [r for r in molecular_rows if r.endswith(IHC_SUFFIX)]
non_ihc_molecular = [r for r in molecular_rows if not r.endswith(IHC_SUFFIX)]

print(f"Clinical rows: {len(clinical_df)}")
print(f"IHC rows: {len(ihc_rows)}")
print(f"Other molecular rows: {len(non_ihc_molecular)}")

# %%
# --- Step 4: Build assessment status per patient per modality ---

assessment_status = {}  # {suffix: {patient: "assessed" | "not assessed"}}

for suffix, tracking in ASSESSMENT_TRACKING.items():
    status = {}
    not_assessed_row = tracking.get("not_assessed")
    assessed_row = tracking.get("assessed")

    for patient in patients:
        if not_assessed_row and not_assessed_row in df.index:
            val = df.loc[not_assessed_row, patient]
            if pd.notna(val) and str(val).strip() == "not assessed":
                status[patient] = "not assessed"
                continue
        if assessed_row and assessed_row in df.index:
            val = df.loc[assessed_row, patient]
            if pd.notna(val) and str(val).strip() == "assessed":
                status[patient] = "assessed"
                continue
        status[patient] = ""  # unknown
    assessment_status[suffix] = status

print("Assessment status built for:", list(assessment_status.keys()))

# %%
# --- Step 5: Classify SNV values ---

snv_rows = [r for r in non_ihc_molecular
            if r.endswith("_SNV")
            and not r.startswith("not_assessed")
            and not r.startswith("assessed")
            and r != "NA_SNV"]

# Also drop the NA_SNV row if it survived cleaning
if "NA_SNV" in df.index:
    # Merge NA_SNV value for Patient 18 into STK11_SNV if needed
    for patient in patients:
        val = df.loc["NA_SNV", patient]
        if pd.notna(val) and str(val).strip() != "":
            # This is the STK11 p.L137fs that was miscategorized
            if "STK11_SNV" in df.index:
                existing = df.loc["STK11_SNV", patient]
                if pd.isna(existing) or str(existing).strip() == "":
                    df.loc["STK11_SNV", patient] = val
                    print(f"  Merged NA_SNV -> STK11_SNV for {patient}: {val}")

gene_data_rows = []

for row_label in snv_rows:
    gene = row_label.replace("_SNV", "")
    classified = {}
    for patient in patients:
        raw_val = df.loc[row_label, patient]
        classified[patient] = classify_snv(raw_val, gene)
    gene_data_rows.append({"gene": gene, "modality": "SNV", "values": classified})

print(f"Classified {len(snv_rows)} SNV rows")

# %%
# --- Step 6: Classify Del/Amp/Sig values ---

for suffix, vocab_prefix in [("_Del", "Del"), ("_Amp", "Amp"), ("_Sig", "Sig")]:
    feature_rows = [r for r in non_ihc_molecular
                    if r.endswith(suffix)
                    and not r.startswith("not_assessed")
                    and not r.startswith("assessed")]
    classifier = {"_Del": classify_del, "_Amp": classify_amp, "_Sig": classify_sig}[suffix]

    for row_label in feature_rows:
        gene = row_label.replace(suffix, "")
        classified = {}
        status_map = assessment_status.get(suffix, {})
        for patient in patients:
            raw_val = df.loc[row_label, patient]
            patient_status = status_map.get(patient, "")
            classified[patient] = classifier(raw_val, patient_status)
        gene_data_rows.append({"gene": gene, "modality": vocab_prefix, "values": classified})

    print(f"Classified {len(feature_rows)} {suffix} rows")

# %%
# --- Step 7: Parse fusions ---

fusion_entries = {}  # {gene: {patient: "fus det."}}

if fusions_row is not None:
    for patient in patients:
        fusion_val = fusions_row[patient]
        partner_genes = parse_fusions(fusion_val, cohort=COHORT_NAME)
        for gene in partner_genes:
            if gene not in fusion_entries:
                fusion_entries[gene] = {}
            fusion_entries[gene][patient] = "fus det."

    for gene, patient_vals in fusion_entries.items():
        classified = {p: patient_vals.get(p, "") for p in patients}
        gene_data_rows.append({"gene": gene, "modality": "Fus", "values": classified})

    print(f"Parsed fusions: {len(fusion_entries)} partner genes with entries")

# %%
# --- Step 8: Build gene-data matrix ---

# Combine all modality rows per gene
records = []
for entry in gene_data_rows:
    row = {"Gene/Biomarker": entry["gene"]}
    row.update(entry["values"])
    records.append(row)

gene_df = pd.DataFrame(records)
gene_df = gene_df.set_index("Gene/Biomarker")

# Replace empty strings with NaN for frequency counting
print(f"Gene-data matrix: {gene_df.shape}")

# %%
# --- Step 9: Classify TMB + MMR as complex biomarkers ---

complex_rows = {}

if tmb_row is not None:
    tmb_classified = {}
    for patient in patients:
        tmb_classified[patient] = classify_tmb(
            tmb_row[patient], high=TMB_HIGH_THRESHOLD,
            intermediate=TMB_INTERMEDIATE_THRESHOLD)
    complex_rows["TMB"] = tmb_classified

if mmr_row is not None:
    mmr_classified = {}
    for patient in patients:
        val = mmr_row[patient]
        if pd.isna(val) or str(val).strip() == "":
            mmr_classified[patient] = "not assessed"
        elif str(val).strip() == "dMMR":
            mmr_classified[patient] = "positive"
        elif str(val).strip() == "pMMR":
            mmr_classified[patient] = "assessed"
        else:
            mmr_classified[patient] = "not assessed"
    complex_rows["dMMR"] = mmr_classified

# HER2 IHC score as biomarker (4 levels: 0, 1+, 2+, 3+)
her2_row = df.loc["HER2_IHC"] if "HER2_IHC" in df.index else None
if her2_row is not None:
    her2_classified = {}
    for patient in patients:
        her2_classified[patient] = classify_her2(her2_row[patient])
    complex_rows["HER2"] = her2_classified

# Methylation placeholder (future biomarker, column not yet final)
methylation_row = df.loc["Methylation analysis (column not yet final!)"] \
    if "Methylation analysis (column not yet final!)" in df.index else None
if methylation_row is not None:
    methyl_classified = {}
    for patient in patients:
        val = methylation_row[patient]
        if pd.isna(val) or str(val).strip() in ("", "na", "not assessed"):
            methyl_classified[patient] = "not assessed"
        else:
            methyl_classified[patient] = "positive"
    complex_rows["Methylation"] = methyl_classified

complex_df = pd.DataFrame(complex_rows).T
complex_df.index.name = "Gene/Biomarker"
print(f"Complex biomarkers: {list(complex_df.index)}")

# %%
# --- Step 10: Frequency threshold and aggregation ---

# Count truly informative entries per gene (exclude "not assessed", "assessed",
# and vocabulary tokens that just mean "panel was run but nothing found")
_NON_INFORMATIVE = {"", "not assessed", "assessed", "CNV assessed",
                    "mut assessed", "fus assessed", "mRNAexp assessed"}
_informative_mask = gene_df.map(lambda v: v not in _NON_INFORMATIVE)
gene_counts = _informative_mask.sum(axis=1)

# Group by gene name (a gene may have multiple modality rows)
gene_freq = gene_counts.groupby(gene_counts.index).sum().sort_values(ascending=False)

single_genes = gene_freq[gene_freq >= GENE_FREQUENCY_THRESHOLD].index
aggregated_genes = gene_freq[gene_freq < GENE_FREQUENCY_THRESHOLD].index

print(f"Genes above threshold ({GENE_FREQUENCY_THRESHOLD}): {list(single_genes)}")
print(f"Genes below threshold (aggregated): {list(aggregated_genes)}")

# %%
# --- Step 11: Fill empty cells with 'not assessed' ---

# For cells that are still empty (no classification), set to 'not assessed'
gene_df = gene_df.fillna("not assessed")
gene_df = gene_df.replace("", "not assessed")

# %%
# --- Step 12: Validate and export ---

os.makedirs("output", exist_ok=True)

# Gene data: single rows (above threshold) + complex biomarkers
gene_single = gene_df.loc[gene_df.index.isin(single_genes)]
gene_data_single = pd.concat([gene_single, complex_df], axis=0)
gene_data_single.to_csv(os.path.join("output", f"{STUDY_DATE}_gene-data-single-rows.csv"))
print(f"Exported gene-data-single-rows: {gene_data_single.shape}")

# Gene data: aggregated rows (below threshold)
gene_aggregated = gene_df.loc[gene_df.index.isin(aggregated_genes)]
gene_aggregated.index = pd.Index(["Other"] * len(gene_aggregated))
gene_aggregated.to_csv(os.path.join("output", f"{STUDY_DATE}_gene-data-aggregated-rows.csv"))
print(f"Exported gene-data-aggregated-rows: {gene_aggregated.shape}")

# Clinical metadata
clinical_export = clinical_df.copy()
clinical_export.to_csv(os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_clinical.csv"))
print(f"Exported clinical: {clinical_export.shape}")

# Complex biomarkers
complex_df.to_csv(os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_complex.csv"))
print(f"Exported complex: {complex_df.shape}")

# IHC heatmap (numeric values) — exclude rows handled as biomarkers
IHC_BIOMARKER_ROWS = {"HER2_IHC", "MMR_IHC"}  # now in complex_df
ihc_heatmap_rows = [r for r in ihc_rows if r not in IHC_BIOMARKER_ROWS]
ihc_df = df.loc[ihc_heatmap_rows].copy()
# Normalize IHC index: strip _IHC suffix for display
ihc_df.index = pd.Index([r.replace("_IHC", "") for r in ihc_df.index])
# Convert to numeric where possible
for col in ihc_df.columns:
    ihc_df[col] = pd.to_numeric(ihc_df[col], errors="coerce")
ihc_df.to_csv(os.path.join("output", f"{STUDY_DATE}_{COHORT_NAME}_ihc_heatmap.csv"))
print(f"Exported IHC heatmap: {ihc_df.shape}")

# %%
# --- Step 13: Validate marker values ---

combined_for_validation = pd.concat([gene_data_single, gene_aggregated])
# Only validate non-IHC data (IHC is numeric, not vocabulary-mapped)
unknown = validate_marker_values(combined_for_validation, strict=False)
if not unknown:
    print("\nAll marker values are valid.")
else:
    print(f"\nWARNING: Unknown marker values: {unknown}")

print(f"\nDone! CSVs are in the output/ directory.")
