import glob
import pandas as pd
import subprocess
import os

# Paths
#fragment_dir = "/mnt/appl/Erlantz/scNMT/Arina/experiment/Output_mouse/Bam_files/noDup/Plates/fragment_files/CpG_files/"
#fragment_dir = "/mnt/appl/Erlantz/scNMT/Pool_40/gDNA/Output_mouse/unique/Bam_files/fragment_files/"
#fragment_dir2 = "/mnt/appl/Erlantz/scNMT/Pool_35/gDNA/Output_mouse/noDup/Plate_2/fragment_files/"
fragment_dir2 = "/mnt/appl/Erlantz/scNMT/Pool_40/gDNA/Output_mouse/unique/Bam_files/fragment_files_batch1_2/"
output_prefix = "Pool35_40_merged_CpG_fragments_batch1_2"
output_tsv = os.path.join(fragment_dir2, output_prefix + ".tsv")
output_bgz = output_tsv + ".gz"

# Read all .tsv.gz files
all_files = glob.glob(os.path.join(fragment_dir2, "*.tsv.gz"))
dfs = []

print(f"Found {len(all_files)} fragment files to merge...")

for file in all_files:
    df = pd.read_csv(file, sep="\t", header=None, compression="gzip",
                     names=["chr", "start", "end", "barcode", "count"])
    dfs.append(df)

# Concatenate and sort
merged_df = pd.concat(dfs)
merged_df.sort_values(by=["chr", "start", "end"], inplace=True)

# Write plain TSV file
merged_df.to_csv(output_tsv, sep="\t", header=False, index=False)
print(f"Written merged TSV file: {output_tsv}")

# Compress with bgzip
subprocess.run(["bgzip", "-f", output_tsv])  # -f to overwrite if exists
print(f"Compressed with bgzip: {output_bgz}")

# Index with tabix
subprocess.run(["tabix", "-p", "bed", output_bgz])
print(f"Tabix index created: {output_bgz}.tbi")
