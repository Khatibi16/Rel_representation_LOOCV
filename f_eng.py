import os
import json
import shutil

# ====== CONFIG ======
metadata_path = "/path/to/metadata.cart.json"
download_dir = "/path/to/tcga_brca_data"  # Folder where gdc-client put UUID folders
flat_dir = "/path/to/flat_files"  # Where you want all .tsv files together
os.makedirs(flat_dir, exist_ok=True)

# ====== LOAD METADATA ======
with open(metadata_path) as f:
    meta = json.load(f)

uuid_to_case = {}
uuid_to_filename = {}

for item in meta:
    file_id = item["file_id"]  # UUID that matches folder name
    filename = item["file_name"]
    uuid_to_filename[file_id] = filename
    
    # Find the patient barcode inside 'associated_entities'
    case_ids = [
        e["entity_submitter_id"]
        for e in item.get("associated_entities", [])
        if e.get("entity_type") == "case" and e.get("entity_submitter_id")
    ]
    if case_ids:
        uuid_to_case[file_id] = case_ids[0]  # Take the first case ID

# ====== FLATTEN & RENAME ======
for uuid, filename in uuid_to_filename.items():
    folder_path = os.path.join(download_dir, uuid)
    file_path = os.path.join(folder_path, filename)
    
    if not os.path.exists(file_path):
        print(f"⚠️ Missing file for {uuid}")
        continue
    
    # Rename using patient barcode (if found)
    patient_id = uuid_to_case.get(uuid, "UNKNOWN")
    new_name = f"{patient_id}.tsv"
    dest_path = os.path.join(flat_dir, new_name)
    
    shutil.copy(file_path, dest_path)  # Copy or move: shutil.move(...)
    print(f"✅ Copied {filename} → {new_name}")

print("🎯 Flattening complete. All files are in:", flat_dir)
