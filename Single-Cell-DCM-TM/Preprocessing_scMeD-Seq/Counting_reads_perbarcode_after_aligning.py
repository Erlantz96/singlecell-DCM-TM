import pysam
import os

# Directory containing the input BAM files
#input_dir = "/mnt/appl/Erlantz/Run15_plate_based_CS_RNA-DNA/Deepseq/nova240/Erlantz/Output/new_subset/aligned/NPCEA/"
#input_dir = "/mnt/appl/Erlantz/scNMT/Pool_40/gDNA/Output_mouse/Aligned/Plate_3"
input_dir ="/mnt/appl/Erlantz/Run5_221122/Plate_input/Gastruloids/Output_mouse/Aligned/Plate_1/"
# Read conditions from a text file into a dictionary
conditions = {}
with open("Plate_Conditions/Gastruloids/Gast_plates.txt", "r") as condition_file:
    for line in condition_file:
        columns = line.strip().split()
        if len(columns) == 2:
            condition, second_column = columns
            conditions[int(condition)] = second_column

# Create a dictionary to store counts for each barcode and second column
barcode_counts = {}

# Iterate through all BAM files in the directory
for input_bam_file in os.listdir(input_dir):
    if input_bam_file.endswith(".sam"):
        # Determine the suffix based on the input BAM file name
        if "DCM" in input_bam_file:
            output_suffix = "DCM"
        elif "CpG" in input_bam_file:
            output_suffix = "CpG"
        else:
            output_suffix = "Unknown"

        # Open the input BAM file
        input_bam = pysam.AlignmentFile(os.path.join(input_dir, input_bam_file), "rb")

        # Iterate through the reads in the input BAM file
        for read in input_bam:
            # Check if the read is aligned
            if not read.is_unmapped:
                # Parse the read name
                read_name = read.query_name
                # Split the read name by "#" and get the value at index 8
                read_condition = int(read_name.split("#")[1])

                # Check if the condition matches and retrieve the second column
                if read_condition in conditions:
                    second_column = conditions[read_condition]
                    output_bam_name = f"{second_column}_{output_suffix}.sam"

                    # Create a dictionary key for barcode and second column
                    key = (output_bam_name, second_column)
                    barcode = int(read_name.split("#")[1])

                    # Increment the count for this barcode and second column
                    if key in barcode_counts:
                        if barcode in barcode_counts[key]:
                            barcode_counts[key][barcode] += 1
                        else:
                            barcode_counts[key][barcode] = 1
                    else:
                        barcode_counts[key] = {barcode: 1}

        # Close the input BAM file
        input_bam.close()

# Create count.files.txt for each output_bam_name
for (output_bam_name, second_column), barcode_data in barcode_counts.items():
    output_file_name = f"/mnt/appl/Erlantz/Run5_221122/Plate_input/Gastruloids/Output_mouse/counts/Aligned_plate1_{output_bam_name}_count_barcodes.txt"
    with open(output_file_name, "w") as output_file:
        output_file.write("Barcode\tSecond_Column\tCount\n")
        for barcode, count in barcode_data.items():
            output_file.write(f"{barcode}\t{second_column}\t{count}\n")
