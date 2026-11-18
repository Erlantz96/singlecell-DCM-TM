import pysam
from collections import defaultdict
import subprocess
import os
import re

# === Input/Output Paths ===
#input_dir = "/mnt/appl/Erlantz/scNMT/pool_30_31/DNA/deepseq/fastq/240818_A00383_0709_AHLJ5WDRX2/Creyghton/Output_human/BAM_files/"
#output_dir = "/mnt/appl/Erlantz/scNMT/pool_30_31/DNA/deepseq/fastq/240818_A00383_0709_AHLJ5WDRX2/Creyghton/Output_human/fragment_files/"
#input_dir="/mnt/appl/Erlantz/scNMT/Arina/experiment/Output_mouse/Bam_files/noDup/Plates/Plate_4/"
#output_dir="/mnt/appl/Erlantz/scNMT/Arina/experiment/Output_mouse/Bam_files/noDup/Plates/fragment_files/CpG_files/"
#input_dir="/mnt/appl/Erlantz/scNMT/Pool_40/gDNA/Output_mouse/unique/Bam_files/Plate_3_IGE/"
#output_dir="/mnt/appl/Erlantz/scNMT/Pool_40/gDNA/Output_mouse/unique/Bam_files/fragment_files/"
#input_dir="/mnt/appl/Erlantz/scNMT/Pool_37/Pool_35_resequenced/Output_mouse/Merged_sequence_files/noDup/Plate2_B/DCM/"
#output_dir="/mnt/appl/Erlantz/scNMT/scNMT/Pool_37/Pool_35_resequenced/Output_mouse/Merged_sequence_files/noDup/fragment_files"
#input_dir="/mnt/appl/Erlantz/scNMT/Pool_37/Pool_35_resequenced/Output_mouse/Merged_sequence_files/noDup/Plate1_A/"
#output_dir="/mnt/appl/Erlantz/scNMT/Pool_37/Pool_35_resequenced/Output_mouse/Merged_sequence_files/noDup/fragment_files"
#input_dir = "/mnt/appl/Erlantz/scNMT/Pool_38/CTC_melanoma/Output_human/CpG_only/noDup/"
#output_dir = "/mnt/appl/Erlantz/scNMT/Pool_38/CTC_melanoma/Output_human/CpG_only/noDup/fragment_files/"
#input_dir = "/mnt/appl/Erlantz/Run5_221122/Combi_input/Output_mouse/DCM/"
#output_dir = "/mnt/appl/Erlantz/Run5_221122/Combi_input/Output_mouse/fragment_files/"

input_dir = "/mnt/appl/Erlantz/scNMT/Pool_35/gDNA/Output_mouse/noDup/Plate_2/"
output_dir = "/mnt/appl/Erlantz/scNMT/Pool_35/gDNA/Output_mouse/noDup/Plate_2/fragment_files/"

#index_file_path = "/mnt/appl/Erlantz/Annotations/Sites/MM10_DCM_index.bed" #  MM10_Lpnp1CpG_index.bed  # Replace if needed MM10_DCM_index.bed
index_file_path = "/mnt/appl/Erlantz/Annotations/Sites/MM10_Lpnp1CpG_index.bed" #  MM10_Lpnp1CpG_index.bed  # Replace if needed MM10_DCM_index.bed
# === Create Output Directory ===
os.makedirs(output_dir, exist_ok=True)

# === Load Index File ===
site_to_region = {}
with open(index_file_path, "r") as index_file:
    for line_number, line in enumerate(index_file, start=1):
        parts = line.strip().split()
        if len(parts) < 5:
            print(f"Skipping malformed line {line_number} in index file: {line.strip()}")
            continue
        chrom, start, end, _, site = parts[:5]
        site_to_region[site] = (chrom, int(start), int(end))

# === Get BAM Files with 'DCM' in the Name ===
bam_files = [f for f in os.listdir(input_dir) if f.endswith('.bam') and 'CpG' in f]

for bam_file in bam_files:
    # Extract plate number
    plate_match = re.search(r'^([^_]+)_([^_]+)_', bam_file)
    if plate_match:
        plate_number = plate_match.group(2)
    else:
        print(f"Warning: Could not extract plate number from filename {bam_file}, skipping")
        continue

    base_name = os.path.splitext(bam_file)[0]
    output_tsv_file = os.path.join(output_dir, f"{base_name}_fragments_from_index.tsv")
    output_bgzip_file = output_tsv_file + ".gz"

    print(f"🚀 Processing {bam_file} (Plate {plate_number})...")

    fragment_counts = defaultdict(int)
    bam_path = os.path.join(input_dir, bam_file)
    bam = pysam.AlignmentFile(bam_path, "rb")

    for read in bam:
        try:
            site = read.query_name.split(":")[8]
        except IndexError:
            print(f"⚠️ Malformed read name: {read.query_name}")
            continue

        if site not in site_to_region:
            continue

        chrom, start, end = site_to_region[site]
        barcode = read.query_name.split(":")[7][3:] + f"-1-{plate_number}"
        fragment_counts[(chrom, start, end, barcode)] += 1

    bam.close()

    # Write fragment counts to TSV
    with open(output_tsv_file, "wt") as out:
        for (chrom, start, end, barcode), count in fragment_counts.items():
            out.write(f"{chrom}\t{start}\t{end}\t{barcode}\t{count}\n")

    # === Sort TSV before compression (required for tabix) ===
    sort_result = subprocess.run(
        ["sort", "-k1,1", "-k2,2n", output_tsv_file, "-o", output_tsv_file],
        capture_output=True, text=True
    )
    if sort_result.returncode != 0:
        print(f"❌ Sort failed: {sort_result.stderr}")
        continue

    # === Compress using bgzip ===
    bgzip_result = subprocess.run(["bgzip", "-f", output_tsv_file], capture_output=True, text=True)
    if bgzip_result.returncode != 0:
        print(f"❌ bgzip failed: {bgzip_result.stderr}")
        continue

    # === Index using tabix ===
    tabix_result = subprocess.run(["tabix", "-p", "bed", output_bgzip_file], capture_output=True, text=True)
    if tabix_result.returncode != 0:
        print(f"❌ tabix failed: {tabix_result.stderr}")
        continue

    print(f"✅ Fragment file indexed: {output_bgzip_file}")

print("🎉 All files processed successfully!")

