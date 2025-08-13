# flatten_and_rename.py
import os, json, shutil, glob

# ====== CONFIG (EDIT THESE) ======
METADATA_JSON = "metadata.repository.2025-08-13.json"  # Path to GDC metadata file
GDC_DOWNLOAD_DIR = "tcga_brca_data"   # the folder with many UUID subfolders
FLAT_DIR = "flat_files"               # output folder for renamed files
USE_SYMLINKS = True                                    # True = make symlinks instead of copying

os.makedirs(FLAT_DIR, exist_ok=True)

# ====== LOAD METADATA ======
with open(METADATA_JSON) as f:
    meta = json.load(f)

records = []
for item in meta:
    file_id = item.get("file_id")              # matches the UUID folder name
    file_name = item.get("file_name")          # inside that folder
    # pull case barcode and sample submitter id if present
    case_barcode = None
    sample_submitter = None
    for ent in item.get("associated_entities", []):
        if ent.get("entity_type") == "case" and ent.get("entity_submitter_id"):
            case_barcode = ent["entity_submitter_id"]
        if ent.get("entity_type") in {"sample", "aliquot"} and ent.get("entity_submitter_id"):
            sample_submitter = ent["entity_submitter_id"]
    records.append({
        "file_id": file_id,
        "file_name": file_name,
        "case_barcode": case_barcode,
        "sample_submitter": sample_submitter
    })

# ====== SCAN & RENAME ======
seen_case_counts = {}
copied = 0
skipped = 0
for r in records:
    uuid_folder = os.path.join(GDC_DOWNLOAD_DIR, r["file_id"])
    if not os.path.isdir(uuid_folder):
        print(f"⚠️  Missing folder for {r['file_id']}, skipping")
        skipped += 1
        continue

    # find the .tsv file inside (GDC puts exactly one data file in each folder)
    candidates = sorted(glob.glob(os.path.join(uuid_folder, "*.tsv")))
    if not candidates:
        print(f"⚠️  No .tsv inside {uuid_folder}, skipping")
        skipped += 1
        continue
    src = candidates[0]

    case = r["case_barcode"] or "UNKNOWN_CASE"
    sample = r["sample_submitter"]

    # detect duplicates: if a case appears multiple times, append sample id
    seen_case_counts[case] = seen_case_counts.get(case, 0) + 1
    if seen_case_counts[case] > 1 and sample:
        out_name = f"{case}__{sample}.tsv"
    elif seen_case_counts[case] > 1 and not sample:
        out_name = f"{case}__dup{seen_case_counts[case]}.tsv"
    else:
        out_name = f"{case}.tsv"

    dst = os.path.join(FLAT_DIR, out_name)
    if os.path.exists(dst):
        # fallback to non-clobber name
        base, ext = os.path.splitext(out_name)
        dst = os.path.join(FLAT_DIR, f"{base}__{r['file_id']}{ext}")

    if USE_SYMLINKS:
        try:
            os.symlink(src, dst)
        except FileExistsError:
            pass
    else:
        shutil.copy2(src, dst)

    copied += 1
    if copied % 50 == 0:
        print(f"… {copied} files processed")

print(f"✅ Done. Copied/Symlinked: {copied}, Skipped: {skipped}")
print(f"Output folder: {FLAT_DIR}")
