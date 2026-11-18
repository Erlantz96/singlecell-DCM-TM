import pysam
import os
import re
import pandas as pd
from collections import defaultdict

# Directory containing the noDup BAM files
#input_dir_nodup = "/mnt/appl/Erlantz/Run2_19102022/Output/sorted/noDup/noDup"
#input_dir_nodup = "/mnt/appl/Erlantz/Run8_scMed_SI_NPC_Kidney/Output/noDup"
#input_dir_nodup =  "/mnt/appl/Erlantz/Run9_RNA-DNA_combi_biotin_U87_ESMEF/MEFs/Output/noDup/"
#input_dir_nodup = "/mnt/appl/Erlantz/Run10_16062023/230610_A00383_0454_AHTC3VDRX2/Erlantz/Output/noDup/"
input_dir_nodup = "/mnt/appl/Erlantz/scNMT/Pool_37/Chase_combi_enterocytes/Output_mouse/noDup/"
#input_dir_nodup = "/mnt/appl/Erlantz/scNMT/pool_32/Combi/Output_mouse/noDup/"
# Directory containing the withDup BAM files
#input_dir_withdup = "/mnt/appl/Erlantz/Run2_19102022/Output/sorted/Dup/"
#input_dir_withdup = "/mnt/appl/Erlantz/Run8_scMed_SI_NPC_Kidney/Output/Dup/"
#input_dir_withdup =  "/mnt/appl/Erlantz/Run9_RNA-DNA_combi_biotin_U87_ESMEF/MEFs/Output/Dup/"
#input_dir_withdup = "/mnt/appl/Erlantz/Run10_16062023/230610_A00383_0454_AHTC3VDRX2/Erlantz/Output/Dup/"
input_dir_withdup = "/mnt/appl/Erlantz/scNMT/Pool_37/Chase_combi_enterocytes/Output_mouse/Dup/"
#input_dir_withdup = "/mnt/appl/Erlantz/scNMT/pool_32/Combi/Output_mouse/Dup/"
# Directory containing the SAM files
#sam_dir = "/mnt/appl/Erlantz/Run2_19102022/Output/Aligned/"
#sam_dir = "/mnt/appl/Erlantz/Run8_scMed_SI_NPC_Kidney/Output/Aligned/"
#sam_dir = "/mnt/appl/Erlantz/Run9_RNA-DNA_combi_biotin_U87_ESMEF/MEFs/Output/Aligned/"
#sam_dir = "/mnt/appl/Erlantz/Run10_16062023/230610_A00383_0454_AHTC3VDRX2/Erlantz/Output/Aligned/"
sam_dir = "/mnt/appl/Erlantz/scNMT/Pool_37/Chase_combi_enterocytes/Output_mouse/Aligned/"
#sam_dir = "/mnt/appl/Erlantz/scNMT/pool_32/Combi/Output_mouse/Aligned/"

#index_dir = "/mnt/appl/Erlantz/Run2_19102022/Output/index/"
#index_dir = "/mnt/appl/Erlantz/Run8_scMed_SI_NPC_Kidney/Output/index/"
#index_dir = "/mnt/appl/Erlantz/Run9_RNA-DNA_combi_biotin_U87_ESMEF/MEFs/Output/index/"
#index_dir = "/mnt/appl/Erlantz/Run10_16062023/230610_A00383_0454_AHTC3VDRX2/Erlantz/Output/index/"
index_dir = "/mnt/appl/Erlantz/scNMT/Pool_37/Chase_combi_enterocytes/Output_mouse/index/"
#index_dir = "/mnt/appl/Erlantz/scNMT/pool_32/Combi/Output_mouse/index/"
#matrix_dir = "/mnt/appl/Erlantz/Run2_19102022/Output/Matrixes/"
#matrix_dir = "/mnt/appl/Erlantz/Run8_scMed_SI_NPC_Kidney/Output/Matrixes/"
#matrix_dir = "/mnt/appl/Erlantz/Run9_RNA-DNA_combi_biotin_U87_ESMEF/MEFs/Output/Matrixes/"
#matrix_dir = "/mnt/appl/Erlantz/Run10_16062023/230610_A00383_0454_AHTC3VDRX2/Erlantz/Output/Matrixes/"
matrix_dir = "/mnt/appl/Erlantz/scNMT/Pool_37/Chase_combi_enterocytes/Output_mouse/Matrixes/"
#matrix_dir = "/mnt/appl/Erlantz/scNMT/pool_32/Combi/Output_mouse/Matrixes/"
# Output directory for CSV files
#output_dir = "/mnt/appl/Erlantz/Run2_19102022/Output/sorted/Counts_barcodes/"
#output_dir = "/mnt/appl/Erlantz/Run8_scMed_SI_NPC_Kidney/Output/Counts_barcodes/"
#output_dir = "/mnt/appl/Erlantz/Run9_RNA-DNA_combi_biotin_U87_ESMEF/MEFs/Output/Counts_barcodes/"
#output_dir = "/mnt/appl/Erlantz/Run10_16062023/230610_A00383_0454_AHTC3VDRX2/Erlantz/Output/Counts_barcodes/"
output_dir = "/mnt/appl/Erlantz/scNMT/Pool_37/Chase_combi_enterocytes/Output/Counts_barcodes/"
#output_dir = "/mnt/appl/Erlantz/scNMT/pool_32/Combi/Output_mouse/Counts_barcodes/"

os.makedirs(output_dir, exist_ok=True)

# Directory containing the original TSV files
#original_dir = "/mnt/appl/Erlantz/Run2_19102022/3_Process_input/"
#original_dir = "/mnt/appl/Erlantz/Run8_scMed_SI_NPC_Kidney/3_Process_input/cell_Stats/"
#original_dir = "/mnt/appl/Erlantz/Run9_RNA-DNA_combi_biotin_U87_ESMEF/MEFs/3_Process_input/"
#original_dir = "/mnt/appl/Erlantz/Run10_16062023/230610_A00383_0454_AHTC3VDRX2/Erlantz/3_Process_input/cell_stats/"
original_dir = "/mnt/appl/Erlantz/scNMT/Pool_37/Chase_combi_enterocytes/3_Process_input_mouse/"
#original_dir = "/mnt/appl/Erlantz/scNMT/pool_32/Combi/3_Process_input/"

# ----------------------------------------
# NAME EXTRACTION FUNCTIONS
# ----------------------------------------

def get_sample_name_bam(filename):
    """Extract sample name from BAM: first three words"""
    match = re.match(r"([^_]+_[^_]+_[^_]+)_", filename)
    return match.group(1) if match else os.path.splitext(filename)[0]

def get_sample_name_tsv(filename):
    """Identify cell_stats files and extract first three words as sample name"""
    if not filename.endswith("_cell_stats.tsv"):
        return None
    match = re.match(r"([^_]+_[^_]+_[^_]+)_", filename)
    return match.group(1) if match else os.path.splitext(filename)[0]

def get_sample_name_index(filename):
    """Identify index files and extract first three words as sample name"""
    if not filename.endswith("_index.tsv"):
        return None
    match = re.match(r"([^_]+_[^_]+_[^_]+)_", filename)
    return match.group(1) if match else os.path.splitext(filename)[0]

def get_sample_name_sam(filename):
    """Extract sample name from SAM: always first three words regardless of suffix"""
    match = re.match(r"([^_]+_[^_]+_[^_]+)_", filename)
    return match.group(1) if match else os.path.splitext(filename)[0]

def get_sample_name_matrix(filename):
    """Extract sample name from matrix file: first three words before _rcMatrix"""
    match = re.match(r"([^_]+_[^_]+_[^_]+)_rcMatrix", filename)
    return match.group(1) if match else os.path.splitext(filename)[0]

def get_assay_type(filename):
    if "DCM" in filename:
        return "DCM"
    elif "CpG" in filename:
        return "CpG"
    return "Unknown"

def get_dup_status(filename):
    return "noDup" if "noDup" in filename else "withDup"

# ----------------------------------------
# HELPER FUNCTIONS
# ----------------------------------------

def save_counts_csv(counts_dict, col_key, sample_name, output_dir):
    """Save counts — key is barcode_number (int), value is count"""
    df = pd.DataFrame(list(counts_dict.items()),
                      columns=["barcode_number", col_key])
    df = df.sort_values("barcode_number").reset_index(drop=True)
    out_path = os.path.join(output_dir, f"{sample_name}_{col_key}_counts.csv")
    df.to_csv(out_path, index=False)
    print(f"  Saved counts: {out_path}  ({len(df)} barcodes)")

def load_existing_counts(sample_name, col_key, output_dir):
    out_path = os.path.join(output_dir, f"{sample_name}_{col_key}_counts.csv")
    if os.path.exists(out_path):
        df = pd.read_csv(out_path)
        if "barcode_number" in df.columns:
            counts = dict(zip(df["barcode_number"].astype(int), df[col_key]))
        elif "combi" in df.columns:
            counts = dict(zip(df["combi"].astype(int), df[col_key]))
        else:
            first_col = df.columns[0]
            counts = dict(zip(df[first_col].astype(int), df[col_key]))
        print(f"  Skipping — loaded existing: {out_path}  ({len(counts)} barcodes)")
        return counts
    return None

def process_bam_directory(bam_dir, sample_counts, sample_to_keys, output_dir):
    """Process all BAM files — barcode extracted based on read name format:
       if '#' in name: split by '#' take index 1
       else:           split by ':' take index 7"""
    for bam_filename in sorted(os.listdir(bam_dir)):
        if not bam_filename.endswith(".bam"):
            continue

        sample_name = get_sample_name_bam(bam_filename)
        assay_type  = get_assay_type(bam_filename)
        dup_status  = get_dup_status(bam_filename)
        bam_path    = os.path.join(bam_dir, bam_filename)
        col_key     = f"{assay_type}_{dup_status}"

        sample_to_keys[sample_name].add(col_key)

        if (sample_name, col_key) in sample_counts and len(sample_counts[(sample_name, col_key)]) > 0:
            print(f"  Skipping BAM (already loaded): {bam_filename}")
            continue

        print(f"Processing BAM: {bam_filename}  →  sample='{sample_name}'  col='{col_key}'")

        input_bam = pysam.AlignmentFile(bam_path, "rb")

        for read in input_bam:
            if read.query_name is None:
                continue

            if "#" in read.query_name:
                parts          = read.query_name.split("#")
                barcode_number = int(parts[1]) if len(parts) > 1 else 0
            else:
                parts          = read.query_name.split(":")
                barcode_number = int(parts[7]) if len(parts) > 7 else 0

            sample_counts[(sample_name, col_key)][barcode_number] += 1

        input_bam.close()

        save_counts_csv(sample_counts[(sample_name, col_key)],
                        col_key, sample_name, output_dir)


def process_sam_directory(sam_dir, sample_counts, sample_to_keys, output_dir):
    """Process all SAM files — barcode extracted from '#' split at index 1
       Skips any SAM containing DCM or CpG in the filename"""
    for sam_filename in sorted(os.listdir(sam_dir)):
        if not sam_filename.endswith(".sam"):
            continue
        if "DCM" in sam_filename or "CpG" in sam_filename:
            continue

        sample_name = get_sample_name_sam(sam_filename)
        sam_path    = os.path.join(sam_dir, sam_filename)
        col_key     = "aligned"

        sample_to_keys[sample_name].add(col_key)

        if (sample_name, col_key) in sample_counts and len(sample_counts[(sample_name, col_key)]) > 0:
            print(f"  Skipping SAM (already loaded): {sam_filename}")
            continue

        print(f"Processing SAM: {sam_filename}  →  sample='{sample_name}'  col='{col_key}'")

        input_sam = pysam.AlignmentFile(sam_path, "r")
        for read in input_sam:
            if read.is_unmapped:
                continue
            if read.query_name is None:
                continue
            parts          = read.query_name.split("#")
            barcode_number = int(parts[1]) if len(parts) > 1 else 0
            sample_counts[(sample_name, col_key)][barcode_number] += 1
        input_sam.close()

        save_counts_csv(sample_counts[(sample_name, col_key)],
                        col_key, sample_name, output_dir)


def get_feature_info(filename):
    """Extract assay and feature from matrix filename"""
    match = re.match(r".+_rcMatrix_(.+?)_(.+?)\.tsv$", filename)
    if match:
        return match.group(1), match.group(2)
    return None, None


def process_matrix_file(filepath):
    """Read matrix in chunks — returns total_counts and features_hit dicts"""
    header       = pd.read_csv(filepath, sep="\t", nrows=0)
    barcode_cols = [c for c in header.columns if c not in ("gene", "region_ID")]
    total_counts = defaultdict(int)
    features_hit = defaultdict(int)

    for chunk in pd.read_csv(filepath, sep="\t", chunksize=1000,
                              dtype={c: "int32" for c in barcode_cols}):
        data = chunk[barcode_cols]
        for col in barcode_cols:
            barcode_number               = int(col)
            col_data                     = data[col]
            total_counts[barcode_number] += int(col_data.sum())
            features_hit[barcode_number] += int((col_data > 0).sum())

    return total_counts, features_hit

# ----------------------------------------
# PRE-SCAN: Register all existing count CSVs
# ----------------------------------------

sample_counts  = defaultdict(lambda: defaultdict(int))
sample_to_keys = defaultdict(set)

count_csv_pattern = re.compile(
    r"^([^_]+_[^_]+_[^_]+)_(aligned|DCM_noDup|DCM_withDup|CpG_noDup|CpG_withDup)_counts\.csv$"
)

print("="*60)
print("PRE-SCAN: checking for existing count CSVs")
print("="*60)

for csv_filename in sorted(os.listdir(output_dir)):
    match = count_csv_pattern.match(csv_filename)
    if not match:
        continue
    sample_name = match.group(1)
    col_key     = match.group(2)
    sample_to_keys[sample_name].add(col_key)
    existing = load_existing_counts(sample_name, col_key, output_dir)
    if existing is not None:
        sample_counts[(sample_name, col_key)] = defaultdict(int, existing)

# ----------------------------------------
# STEP 1a: Count reads per barcode from noDup BAMs
# ----------------------------------------

print("="*60)
print("STEP 1a: Processing noDup BAM files")
print("="*60)

process_bam_directory(input_dir_nodup, sample_counts, sample_to_keys, output_dir)

# ----------------------------------------
# STEP 1b: Count reads per barcode from withDup BAMs
# ----------------------------------------

print("="*60)
print("STEP 1b: Processing withDup BAM files")
print("="*60)

process_bam_directory(input_dir_withdup, sample_counts, sample_to_keys, output_dir)

# ----------------------------------------
# STEP 2: Count aligned reads from SAMs
# ----------------------------------------

print("="*60)
print("STEP 2: Processing SAM files")
print("="*60)

process_sam_directory(sam_dir, sample_counts, sample_to_keys, output_dir)

# ----------------------------------------
# STEP 2b: Process matrix files
# ----------------------------------------

print("="*60)
print("STEP 2b: Processing matrix files")
print("="*60)

matrix_results  = defaultdict(dict)
matrix_col_keys = defaultdict(set)

for matrix_filename in sorted(os.listdir(matrix_dir)):
    if not matrix_filename.endswith(".tsv"):
        continue
    if "rcMatrix" not in matrix_filename:
        continue

    sample_name    = get_sample_name_matrix(matrix_filename)
    assay, feature = get_feature_info(matrix_filename)

    if assay is None:
        print(f"  WARNING: Could not parse assay/feature from {matrix_filename}, skipping.")
        continue

    feature_key = f"{assay}_{feature}"
    cache_path  = os.path.join(output_dir, f"{sample_name}_{feature_key}_matrix_counts.csv")

    if os.path.exists(cache_path):
        print(f"  Skipping — loaded existing: {cache_path}")
        cache_df     = pd.read_csv(cache_path)
        total_counts = dict(zip(cache_df["barcode_number"], cache_df[f"{feature_key}_total"]))
        features_hit = dict(zip(cache_df["barcode_number"], cache_df[f"{feature_key}_n_features"]))
        matrix_results[sample_name][feature_key] = (total_counts, features_hit)
        matrix_col_keys[sample_name].add(feature_key)
        continue

    print(f"Processing matrix: {matrix_filename}  →  sample='{sample_name}'  feature='{feature_key}'")

    matrix_path              = os.path.join(matrix_dir, matrix_filename)
    total_counts, features_hit = process_matrix_file(matrix_path)

    cache_df = pd.DataFrame({
        "barcode_number":            sorted(total_counts.keys()),
        f"{feature_key}_total":      [total_counts[b] for b in sorted(total_counts.keys())],
        f"{feature_key}_n_features": [features_hit[b]  for b in sorted(total_counts.keys())]
    })
    cache_df.to_csv(cache_path, index=False)
    print(f"  Saved matrix counts: {cache_path}  ({len(cache_df)} barcodes)")

    matrix_results[sample_name][feature_key] = (total_counts, features_hit)
    matrix_col_keys[sample_name].add(feature_key)

# ----------------------------------------
# STEP 3: Load index TSV and cell_stats TSV files
# ----------------------------------------

print("="*60)
print("STEP 3: Loading TSV files")
print("="*60)

index_data = {}

for index_filename in sorted(os.listdir(index_dir)):
    sample_name = get_sample_name_index(index_filename)
    if sample_name is None:        # not an index file — skip
        continue

    sample_name = get_sample_name_index(index_filename)
    index_path  = os.path.join(index_dir, index_filename)

    df = pd.read_csv(index_path, sep="\t", header=None,
                     names=["barcode_number", "combi"],
                     dtype={"barcode_number": int, "combi": str})

    index_data[sample_name] = df
    print(f"Loaded index: {index_filename}  →  sample='{sample_name}'  ({len(df)} barcodes)")
    print(f"  Index sample (first 3 rows):\n{df.head(3).to_string()}")

original_data = {}

for tsv_filename in sorted(os.listdir(original_dir)):
    sample_name = get_sample_name_tsv(tsv_filename)
    if sample_name is None:        # not a cell_stats file — skip
        continue

    sample_name = get_sample_name_tsv(tsv_filename)
    tsv_path    = os.path.join(original_dir, tsv_filename)

    df = pd.read_csv(tsv_path, sep="\t", header=None,
                     names=["combi", "readcount", "avg_len"],
                     skiprows=1,
                     dtype={"combi": str})

    df = df[["combi", "readcount"]].copy()

    print(f"  cell_stats sample (first 3 rows):\n{df.head(3).to_string()}")

    if sample_name in index_data:
        index_df       = index_data[sample_name]
        df["combi"]    = df["combi"].str.strip()
        index_df["combi"] = index_df["combi"].str.strip()

        df = df.merge(index_df[["barcode_number", "combi"]], on="combi", how="inner")

        print(f"  After merge: {len(df)} barcodes matched to index")
        df["barcode_number"] = df["barcode_number"].astype(int)
        df = df[["barcode_number", "combi", "readcount"]]

    else:
        print(f"  WARNING: No index file found for '{sample_name}' — generating barcode numbers from row order")
        df.insert(0, "barcode_number", range(len(df)))

    if sample_name in original_data:
        existing = original_data[sample_name]
        original_data[sample_name] = pd.concat([existing, df], ignore_index=True)
    else:
        original_data[sample_name] = df

    print(f"Loaded TSV: {tsv_filename}  →  sample='{sample_name}'  ({len(df)} barcodes)")

# ----------------------------------------
# STEP 4: Merge into one CSV per sample
# ----------------------------------------

print("="*60)
print("STEP 4: Merging into one CSV per sample")
print("="*60)

all_samples = set(sample_to_keys.keys()) | set(original_data.keys())

for sample_name in sorted(all_samples):

    if sample_name not in original_data:
        print(f"WARNING: No TSV found for sample '{sample_name}', skipping.")
        continue

    print(f"\nMerging sample: '{sample_name}'")
    print(f"  Keys found: {sorted(sample_to_keys.get(sample_name, set()))}")

    merged_df = original_data[sample_name].copy()
    merged_df = merged_df.rename(columns={"readcount": "Annotated",
                                           "combi":     "combi_seq"})

    desired_order = ["aligned", "CpG_withDup", "CpG_noDup", "DCM_withDup", "DCM_noDup"]
    for col_key in desired_order:
        count_map          = sample_counts.get((sample_name, col_key), {})
        if not count_map:
            print(f"  WARNING: No counts found for col '{col_key}' — filling with 0")
        merged_df[col_key] = merged_df["barcode_number"].map(count_map).fillna(0).astype(int)

    merged_df["aligned_pct"] = merged_df.apply(
        lambda r: round(r["aligned"] / r["Annotated"] * 100, 2)
        if r["Annotated"] > 0 else 0.0, axis=1)

    merged_df["CpG_unique_pct"] = merged_df.apply(
        lambda r: round(r["CpG_noDup"] / r["CpG_withDup"] * 100, 2)
        if r["CpG_withDup"] > 0 else 0.0, axis=1)

    merged_df["DCM_unique_pct"] = merged_df.apply(
        lambda r: round(r["DCM_noDup"] / r["DCM_withDup"] * 100, 2)
        if r["DCM_withDup"] > 0 else 0.0, axis=1)

    for feature_key in sorted(matrix_col_keys.get(sample_name, set())):
        if feature_key not in matrix_results.get(sample_name, {}):
            continue
        total_counts, features_hit = matrix_results[sample_name][feature_key]
        merged_df[f"{feature_key}_total"]      = merged_df["barcode_number"].map(total_counts).fillna(0).astype(int)
        merged_df[f"{feature_key}_n_features"] = merged_df["barcode_number"].map(features_hit).fillna(0).astype(int)

    base_cols   = ["barcode_number", "combi_seq", "Annotated", "aligned",
                   "CpG_withDup", "CpG_noDup", "DCM_withDup", "DCM_noDup",
                   "aligned_pct", "CpG_unique_pct", "DCM_unique_pct"]
    matrix_cols = [c for c in merged_df.columns if c not in base_cols]
    final_cols  = [c for c in base_cols + sorted(matrix_cols) if c in merged_df.columns]
    merged_df   = merged_df[final_cols]

    output_path = os.path.join(output_dir, f"{sample_name}_merged_counts.csv")
    merged_df.to_csv(output_path, index=False)

    print(f"  Saved merged: {output_path}")
    print(f"  {len(merged_df)} barcodes, columns: {list(merged_df.columns)}")

print("\nDone.")
