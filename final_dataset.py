import os
import json
import pandas as pd

# ===== CONFIG =====
metadata_path = "metadata.repository.2025-08-13.json"   # Path to GDC metadata file
flat_dir = "tcga_brca_data"                # Folder with patient-renamed expression files (.tsv)
clinical_path = "clinical.tsv"         # Path to clinical file
output_path = "flat_files"     # Output CSV file

# ===== LOAD METADATA =====
with open(metadata_path) as f:
    meta = json.load(f)

uuid_to_case = {}
for item in meta:
    file_id = item["file_id"]  # folder name in original GDC download
    case_ids = [
        e["entity_submitter_id"]
        for e in item.get("associated_entities", [])
        if e.get("entity_type") == "case" and e.get("entity_submitter_id")
    ]
    if case_ids:
        uuid_to_case[file_id] = case_ids[0]

# ===== MERGE EXPRESSION FILES =====
expr_data = []
gene_index = None

for fname in os.listdir(flat_dir):
    if not fname.endswith(".tsv"):
        continue
    patient_id = fname.replace(".tsv", "")
    df = pd.read_csv(os.path.join(flat_dir, fname), sep="\t", comment="#")
    # Keep only genes (drop N_unmapped, etc.)
    df = df[df["gene_id"].str.startswith("ENSG", na=False)].copy()
    # Remove version suffix from gene IDs
    df["gene_id"] = df["gene_id"].str.split(".").str[0]
    # Select FPKM column
    vals = df.set_index("gene_id")["fpkm_unstranded"]
    if gene_index is None:
        gene_index = vals.index
    expr_data.append(vals.rename(patient_id))

expression_df = pd.concat(expr_data, axis=1)
expression_df = expression_df.loc[~expression_df.index.duplicated(keep="first")]

# ===== LOAD CLINICAL DATA =====
clinical_df = pd.read_csv(clinical_path, sep="\t")

# Select key survival-related clinical columns
clinical_cols = [
    "cases.submitter_id",
    "demographic.vital_status",
    "demographic.days_to_death",
    "diagnoses.days_to_last_follow_up",
    "diagnoses.ajcc_pathologic_stage"
]
available_cols = [c for c in clinical_cols if c in clinical_df.columns]
clinical_subset = clinical_df[available_cols].drop_duplicates(subset=["cases.submitter_id"])

# ===== MATCH EXPRESSION + CLINICAL =====
common_patients = expression_df.columns.intersection(clinical_subset["cases.submitter_id"])
expression_df = expression_df[common_patients]
clinical_subset = clinical_subset.set_index("cases.submitter_id").loc[common_patients]

# ===== FINAL DATASET =====
final_df = pd.concat([clinical_subset, expression_df.T], axis=1)
final_df.to_csv(output_path)

print(f"✅ Final dataset saved to {output_path}")
print(f"Shape: {final_df.shape[0]} samples × {final_df.shape[1]} columns")
