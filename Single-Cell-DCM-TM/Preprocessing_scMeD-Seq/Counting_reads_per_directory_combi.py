import pysam
import os
import csv
import re
from collections import defaultdict

# Directory containing the input BAM files
input_dir = "/mnt/appl/Erlantz/Run2_19102022/Output/sorted/noDup/noDup/"

# Output directory for CSV files
output_dir = "/mnt/appl/Erlantz/Run2_19102022/Output/sorted/Counts_barcodes/"
os.makedirs(output_dir, exist_ok=True)

def get_sample_name(filename):
    """Extract sample name: everything before 'CpG' or 'DCM'"""
    match = re.match(r"(.+?)(?=CpG|DCM)", filename)
    return match.group(1).rstrip("_") if match else os.path.splitext(filename)[0]

def get_assay_type(filename):
    """Extract assay type from filename"""
    if "DCM" in filename:
        return "DCM"
    elif "CpG" in filename:
        return "CpG"
    return "Unknown"

def get_dup_status(filename):
    """Check if file is deduplicated or not"""
    return "noDup" if "noDup" in filename else "withDup"

# {sample_name: {barcode: count}}
sample_counts = defaultdict(lambda: defaultdict(int))

# Iterate through all BAM files in the directory
for bam_filename in sorted(os.listdir(input_dir)):
    if not bam_filename.endswith(".bam"):
        continue

    sample_name = get_sample_name(bam_filename)
    assay_type  = get_assay_type(bam_filename)
    dup_status  = get_dup_status(bam_filename)
    bam_path    = os.path.join(input_dir, bam_filename)

    print(f"Processing: {bam_filename}  →  sample='{sample_name}'  assay='{assay_type}'  dup='{dup_status}'")

    input_bam = pysam.AlignmentFile(bam_path, "rb")

    for read in input_bam:
        if read.query_name is None:
            continue

        parts = read.query_name.split("#")

        if len(parts) > 1:
            barcode = parts[1]
        else:
            barcode = read.query_name  # fallback: use full read name

        # Key combines sample + assay + dup status
        key = f"{sample_name}_{assay_type}_{dup_status}"
        sample_counts[key][barcode] += 1

    input_bam.close()

# Write one CSV per sample+assay+dup combination
for sample_key, barcode_data in sample_counts.items():
    output_path = os.path.join(output_dir, f"{sample_key}_counts.csv")

    with open(output_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Barcode", "Count"])
        for barcode, count in sorted(barcode_data.items()):
            writer.writerow([barcode, count])

    print(f"Saved: {output_path}  ({len(barcode_data)} barcodes)")

print("Done.")
