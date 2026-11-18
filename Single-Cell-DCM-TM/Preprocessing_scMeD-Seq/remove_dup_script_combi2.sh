#!/bin/bash

# Define directory and input files
#Directory="/mnt/appl/Erlantz/scNMT/Pool_37/Chase_combi_enterocytes/Output_mouse/"
#Directory="/mnt/appl/Erlantz/Run5_221122/Combi_input/Output_mouse/"
#input_files=("CI_ES_index_DCM.sam" "CI_ES_index_CpG-M.sam" "Exome_capture_index_CpG-M.sam" "Exome_capture_index_DCM.sam")
Directory="/mnt/appl/Erlantz/Run10_16062023/230610_A00383_0454_AHTC3VDRX2/Erlantz/Output/"
input_files=("NPCD2_20h_index_CpG-M.sam" "NPCD2_20h_index_DCM.sam" "NPCD2_minus_index_CpG-M.sam" "NPCD2_minus_index_DCM.sam" "SI_R2_index_CpG-M.sam" "SI_R2_index_DCM.sam")


# Step 0: Convert to BAM files
for file in "${input_files[@]}"; do
    input_path="${Directory}${file}"
    output_file="${file%.sam}.bam"
    output_path="${Directory}${output_file}"
    samtools view -@ 10 -bS "$input_path" > "$output_path"
done

# Step 1: Remove "#" from the first column of each input BAM file
for file in "${input_files[@]}"; do
    input_bam="${Directory}${file%.sam}.bam"
    output_bam="${Directory}${file%.sam}_modified_reads.bam"
    samtools view -@ 10 -h "$input_bam" | awk -F'\t' -v OFS='\t' '{gsub("#", ":", $1); print}' | samtools view -@ 10 -bS - > "$output_bam"
done

# Step 2: Sort BAM files
for file in "${input_files[@]}"; do
    input_bam="${Directory}${file%.sam}_modified_reads.bam"
    output_bam="${Directory}${file%.sam}_sorted.bam"
    samtools sort -@ 10 "$input_bam" -o "$output_bam"
done

# Step 3: Mark duplicates
for file in "${input_files[@]}"; do
    input_bam="${Directory}${file%.sam}_sorted.bam"
    output_bam="${Directory}${file%.sam}_noDup.bam"
    samtools markdup -@ 10 -r --barcode-name "$input_bam" "$output_bam"
done

# Step 4: Convert final BAM to SAM
for file in "${input_files[@]}"; do
    input_bam="${Directory}${file%.sam}_noDup.bam"
    output_sam="${Directory}${file%.sam}_noDup.sam"
    samtools view -@ 10 -h "$input_bam" -o "$output_sam"
done

# Step 5: Remove intermediate files
rm "${Directory}"*_modified_reads*.bam

