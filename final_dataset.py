import os, json, glob, pandas as pd

# ========= EDIT THESE PATHS =========
METADATA_JSON   = "metadata.repository.2025-08-13.json"  
GDC_DOWNLOAD_DIR = "tcga_brca_data"   
CLINICAL_TSV    = "clinical.tsv"       
OUTPUT_CSV      = "tcga_brca_final_dataset.csv"  
EXPR_COL        = "fpkm_unstranded"   


def parse_case_sample(s):
    """From a TCGA submitter ID (aliquot or sample), return (case, sample)."""
    if not s or not s.startswith("TCGA-"):
        return (None, None)
    parts = s.split("-")
    case = "-".join(parts[:3]) if len(parts) >= 3 else None
    sample = "-".join(parts[:4]) if len(parts) >= 4 else None
    return (case, sample)

def best_case(item):
    # 1) try cases[0].submitter_id
    for c in (item.get("cases") or []):
        cs = c.get("submitter_id")
        if cs and cs.startswith("TCGA-"):
            case, _ = parse_case_sample(cs)
            if case: return case
    # 2) try any entity_submitter_id that starts with TCGA-
    for ent in (item.get("associated_entities") or []):
        sid = ent.get("entity_submitter_id") or ent.get("submitter_id")
        if sid and sid.startswith("TCGA-"):
            case, _ = parse_case_sample(sid)
            if case: return case
    return None

def best_sample(item):
    # 1) cases[0].samples[0].submitter_id
    for c in (item.get("cases") or []):
        for s in (c.get("samples") or []):
            sid = s.get("submitter_id")
            if sid and sid.startswith("TCGA-"):
                _, sample = parse_case_sample(sid)
                if sample: return sample
    # 2) any entity_submitter_id with TCGA-
    for ent in (item.get("associated_entities") or []):
        sid = ent.get("entity_submitter_id") or ent.get("submitter_id")
        if sid and sid.startswith("TCGA-"):
            _, sample = parse_case_sample(sid)
            if sample: return sample
    return None

# ---- load metadata ----
with open(METADATA_JSON, "r") as f:
    meta = json.load(f)

# ---- build entries: file path + case/sample ----
entries = []
for item in meta:
    file_id   = item.get("file_id")
    file_name = item.get("file_name")
    # usually here:
    src = os.path.join(GDC_DOWNLOAD_DIR, file_id or "", file_name or "")
    if not os.path.exists(src):
        # fallback: any .tsv inside that UUID folder
        cands = glob.glob(os.path.join(GDC_DOWNLOAD_DIR, file_id or "", "*.tsv"))
        if not cands:
            # nothing to read for this item
            continue
        src = cands[0]

    case   = best_case(item)
    sample = best_sample(item)

    entries.append({"src": src, "case": case, "sample": sample})

# ---- if multiple samples per case, prefer primary tumor (code '01') ----
def sample_code(s):
    # TCGA-XX-YYYY-01A -> '01'
    if not s: return None
    parts = s.split("-")
    return parts[3][:2] if len(parts) >= 4 else None

chosen = {}  # case -> entry
for e in entries:
    case, sample = e["case"], e["sample"]
    if not case:  # cannot use this without a case label
        continue
    if case not in chosen:
        chosen[case] = e
    else:
        # prefer primary tumor ('01')
        old = chosen[case]
        if sample_code(sample) == "01" and sample_code(old["sample"]) != "01":
            chosen[case] = e

# ---- read expression and assemble matrix ----
expr_series = []
for case, e in chosen.items():
    try:
        df = pd.read_csv(e["src"], sep="\t", comment="#")
    except Exception as ex:
        print(f"Skip {e['src']} ({ex})")
        continue
    # keep true genes
    df = df[df["gene_id"].astype(str).str.startswith("ENSG", na=False)].copy()
    # strip version suffix
    df["gene_id"] = df["gene_id"].str.split(".").str[0]
    if EXPR_COL not in df.columns:
        print(f"Skip {e['src']} (no {EXPR_COL})")
        continue
    s = df.set_index("gene_id")[EXPR_COL]
    s.name = case
    expr_series.append(s)

if not expr_series:
    raise SystemExit("No expression files read. Check paths/EXPR_COL/permissions.")

expr = pd.concat(expr_series, axis=1)
expr = expr[~expr.index.duplicated(keep="first")]

# ---- load and prepare clinical ----
clin = pd.read_csv(CLINICAL_TSV, sep="\t")

def pick(colnames, options):
    for o in options:
        if o in colnames: return o
    return None

case_col   = pick(clin.columns, ["cases.submitter_id", "submitter_id", "case_submitter_id"])
vital_col  = pick(clin.columns, ["demographic.vital_status", "vital_status"])
death_col  = pick(clin.columns, ["demographic.days_to_death", "days_to_death"])
lfu_col    = pick(clin.columns, ["diagnoses.days_to_last_follow_up", "days_to_last_follow_up", "days_to_last_followup"])
stage_col  = pick(clin.columns, ["diagnoses.ajcc_pathologic_stage", "ajcc_pathologic_stage"])

need = [c for c in [case_col, vital_col, death_col, lfu_col, stage_col] if c]
clin2 = clin[need].drop_duplicates(subset=[case_col]).copy()

# survival targets
vs = clin2[vital_col].astype(str).str.lower()
d_death = pd.to_numeric(clin2.get(death_col), errors="coerce")
d_lfu   = pd.to_numeric(clin2.get(lfu_col), errors="coerce")
time  = d_death.where(vs.isin(["dead","deceased","1","true"]), d_lfu)
event = vs.isin(["dead","deceased","1","true"]).astype(int)

clin2 = clin2.set_index(case_col)
clin2["time"] = time.values
clin2["event"] = event.values

# ---- align and save ----
common = expr.columns.intersection(clin2.index)
final = pd.concat([clin2.loc[common], expr[common].T], axis=1)
final.to_csv(OUTPUT_CSV)

print(" Saved:", OUTPUT_CSV)
print(" Patients:", final.shape[0], "| Clinical cols:", clin2.shape[1], "| Genes:", expr.shape[0])
