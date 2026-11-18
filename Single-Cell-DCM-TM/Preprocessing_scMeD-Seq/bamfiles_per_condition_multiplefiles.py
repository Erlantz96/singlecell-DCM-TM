import pysam
import os
from multiprocessing import Pool

# ──────────────────────────────────────────────────────────────────────────────
# Define all (input_dir, barcode_info_file) pairs in matching order.
# Each tuple: (input_dir, barcode_info_file)
# ──────────────────────────────────────────────────────────────────────────────
DIR_BARCODE_PAIRS = [
    (
        "/mnt/appl/Erlantz/scNMT/pool_30_31/DNA/deepseq/fastq/240818_A00383_0709_AHLJ5WDRX2/Creyghton/Output_human/BAM_files/Plate_1",
        "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/Pool30_31_plate1_dox.txt",
    ),
    (
        "/mnt/appl/Erlantz/scNMT/pool_30_31/DNA/deepseq/fastq/240818_A00383_0709_AHLJ5WDRX2/Creyghton/Output_human/BAM_files/Plate_2",
        "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/Pool30_31_plate2_dox.txt",
    ),
    (
        "/mnt/appl/Erlantz/scNMT/Pool_34/241204_A00383_0767_BHT7JCDRX2/Creyghton/gDNA/Output_human/noDup/Plate_1",
        "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/Pool34_plate1_dox.txt",
    ),
    (
        "/mnt/appl/Erlantz/scNMT/Pool_34/241204_A00383_0767_BHT7JCDRX2/Creyghton/gDNA/Output_human/noDup/Plate_3",
        "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/Pool34_plate3_dox.txt",
    ),
    (
        "/mnt/appl/Erlantz/scNMT/Pool_36/gDNA/Output_human/noDup/Plate_1",
        "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/U87_plate_1_L1_minus.txt",
    ),
    (
        "/mnt/appl/Erlantz/scNMT/Pool_36/gDNA/Output_human/noDup/Plate_1_L2",
        "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/U87_plate_1_L2_8h.txt",
    ),
    (
        "/mnt/appl/Erlantz/scNMT/Pool_36/gDNA/Output_human/noDup/Plate_2",
        "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/U87_plate_2_L1_12h.txt",
    ),
    (
        "/mnt/appl/Erlantz/scNMT/Pool_36/gDNA/Output_human/noDup/Plate_3",
        "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/U87_plate_3_L2.txt",
    ),
]

# Main output directory
MAIN_OUTPUT_DIR = "/mnt/appl/Erlantz/scNMT/Pool_36/gDNA/Output_human/Ruben_metagene/"

# Number of parallel workers (set to len(DIR_BARCODE_PAIRS) to run all at once,
# or lower if memory/CPU is limited)
N_WORKERS = 4


# ──────────────────────────────────────────────────────────────────────────────
# Core processing function — handles one (input_dir, barcode_info_file) pair
# ──────────────────────────────────────────────────────────────────────────────
def process_pair(args):
    input_dir, barcode_info_file = args

    print(f"[START] {input_dir}\n        barcode file: {barcode_info_file}")

    os.makedirs(MAIN_OUTPUT_DIR, exist_ok=True)

    # Read barcode → cluster mapping
    clusters = {}
    with open(barcode_info_file, "r") as barcode_file:
        for line in barcode_file:
            barcode, cluster = line.strip().split("\t")
            clusters.setdefault(cluster, []).append(barcode)

    # Build a fast barcode → cluster lookup (avoids nested loop per read)
    barcode_to_cluster = {
        barcode: cluster
        for cluster, barcode_list in clusters.items()
        for barcode in barcode_list
    }

    # Derive a clean prefix from the barcode file name, e.g. "Pool30_31_plate1_dox"
    barcode_prefix = os.path.splitext(os.path.basename(barcode_info_file))[0]

    # Create a subdirectory for this barcode file
    plate_output_dir = os.path.join(MAIN_OUTPUT_DIR, barcode_prefix)
    os.makedirs(plate_output_dir, exist_ok=True)

    # Process every BAM file in this input directory
    for input_bam_name in os.listdir(input_dir):
        if not input_bam_name.endswith(".bam"):
            continue

        print(f"  Processing: {input_bam_name}")
        input_bam_path = os.path.join(input_dir, input_bam_name)

        # Determine output suffix
        if "DCM" in input_bam_name:
            output_suffix = "_DCM.bam"
        elif "CpG" in input_bam_name:
            output_suffix = "_CpG.bam"
        else:
            output_suffix = ".bam"

        output_sams = {}

        with pysam.AlignmentFile(input_bam_path, "rb") as input_sam:
            for read in input_sam:
                read_name = read.query_name
                barcode = read_name.split(":")[7][3:]

                cluster = barcode_to_cluster.get(barcode)
                if cluster is None:
                    continue  # barcode not in any condition — skip

                # e.g. Pool30_31_plate1_dox_Plate_condition1.bam
                output_sam_name = f"{barcode_prefix}_Plate_{cluster}{output_suffix}"
                if output_sam_name not in output_sams:
                    output_sam_path = os.path.join(plate_output_dir, output_sam_name)
                    output_sams[output_sam_name] = pysam.AlignmentFile(
                        output_sam_path, "wb", header=input_sam.header
                    )

                output_sams[output_sam_name].write(read)

        for sam_file in output_sams.values():
            sam_file.close()

    print(f"[DONE]  {input_dir}")


# ──────────────────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    with Pool(processes=N_WORKERS) as pool:
        pool.map(process_pair, DIR_BARCODE_PAIRS)

    print("All pairs processed.")
