import pysam
import os

# Directory containing the input SAM files
#input_dir = "/mnt/appl/Erlantz/scNMT/Pool_38/FFPE_Cervix_cancer/4_Readcount_region_input_human/"
#input_dir = "/mnt/appl/Erlantz/Run5_221122/Plate_input/Gastruloids/Output_mouse/noDup/Plate_2"
input_dir = "/mnt/appl/Erlantz/scNMT/pool_30_31/DNA/deepseq/fastq/240818_A00383_0709_AHLJ5WDRX2/Creyghton/Output_human/BAM_files/Plate_1"
input_dir = "/mnt/appl/Erlantz/scNMT/pool_30_31/DNA/deepseq/fastq/240818_A00383_0709_AHLJ5WDRX2/Creyghton/Output_human/BAM_files/Plate_2"
input_dir = "/mnt/appl/Erlantz/scNMT/Pool_34/241204_A00383_0767_BHT7JCDRX2/Creyghton/gDNA/Output_human/noDup/Plate_1"
input_dir = "/mnt/appl/Erlantz/scNMT/Pool_34/241204_A00383_0767_BHT7JCDRX2/Creyghton/gDNA/Output_human/noDup/Plate_3"
input_dir = "/mnt/appl/Erlantz/scNMT/Pool_36/gDNA/Output_human/noDup/Plate_1"
input_dir = "/mnt/appl/Erlantz/scNMT/Pool_36/gDNA/Output_human/noDup/Plate_1_L2"
input_dir = "/mnt/appl/Erlantz/scNMT/Pool_36/gDNA/Output_human/noDup/Plate_2"
input_dir = "/mnt/appl/Erlantz/scNMT/Pool_36/gDNA/Output_human/noDup/Plate_3"


# Path to the input text file containing barcode information and cluster assignments
barcode_info_file = "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/Pool30_31_plate1_dox.txt"
barcode_info_file = "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/Pool30_31_plate2_dox.txt"
barcode_info_file = "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/Pool34_plate1_dox.txt"
barcode_info_file = "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/Pool34_plate3_dox.txt"
barcode_info_file = "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/U87_plate_1_L1_minus.txt"
barcode_info_file = "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/U87_plate_1_L2_8h.txt"
barcode_info_file = "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/U87_plate_2_L1_12h.txt"
barcode_info_file = "/mnt/appl/Erlantz/scNMT/scMedseq_preprocess_Scripts/Plate_Conditions/U87/U87_plate_3_L2.txt"

# Main output directory
main_output_dir = "/mnt/appl/Erlantz/scNMT/Pool_36/gDNA/Output_human/Ruben_metagene/"
os.makedirs(main_output_dir, exist_ok=True)

# Read barcode information and cluster assignments from the text file into a dictionary
clusters = {}
with open(barcode_info_file, "r") as barcode_file:
    for line in barcode_file:
        barcode, cluster = line.strip().split("\t")
        clusters.setdefault(cluster, []).append(barcode)

# Iterate through all SAM files in the directory
for input_sam_file in os.listdir(input_dir):
    if input_sam_file.endswith(".bam"):
        print("Processing:", input_sam_file)
        # Open the input SAM file
        input_sam = pysam.AlignmentFile(os.path.join(input_dir, input_sam_file), "r")

        # Determine the suffix based on the input SAM file name
        if "DCM" in input_sam_file:
            output_suffix = "_DCM.bam"
        elif "CpG" in input_sam_file:
            output_suffix = "_CpG.bam"
        else:
            output_suffix = ".bam"

        # Create a dictionary to store output SAM files for each cluster
        output_sams = {}

        # Iterate through the reads in the input SAM file
        for read in input_sam:
            # Parse the read name
            read_name = read.query_name
            # Split the read name by ":" and get the barcode
            barcode = read_name.split(":")[7][3:]

            # Check if the barcode is in any cluster
            for cluster, barcode_list in clusters.items():
                if barcode in barcode_list:
                    output_sam_name = f"Plate_{cluster}{output_suffix}"

                    # Create or open the SAM file for this cluster
                    output_sam_file = os.path.join(main_output_dir, output_sam_name)
                    if output_sam_name not in output_sams:
                        output_sams[output_sam_name] = pysam.AlignmentFile(output_sam_file, "w", header=input_sam.header)

                    # Write the read to the corresponding SAM file
                    output_sams[output_sam_name].write(read)

        # Close all output SAM files for this input SAM file
        for sam_file in output_sams.values():
            sam_file.close()

        # Close the input SAM file
        input_sam.close()

