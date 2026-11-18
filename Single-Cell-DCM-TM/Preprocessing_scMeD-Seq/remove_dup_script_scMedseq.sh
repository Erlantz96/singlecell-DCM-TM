#!/bin/bash

# Directory containing .sam files
#input_dir="/mnt/appl/Erlantz/scNMT/Pool_40/gDNA/Output_mouse"
input_dir="/mnt/appl/Erlantz/Run5_221122/Plate_input/Gastruloids/Output_mouse"
# Step 0: Convert .sam files to .bam
for file in "$input_dir"/*.sam; do
    output_file="${file%.sam}.bam"
    samtools view -@ 14 -bS "$file" > "$output_file"
done

# Step 1: Remove "-" from the first column
for file in "$input_dir"/*.sam; do
    input_file="${file%.sam}.bam"
    output_file="${file%.sam}_modified_reads.bam"
    samtools view -@ 14 -h "$input_file" | awk -F'\t' -v OFS='\t' '{gsub("#", ":", $1); print}' | samtools view -@ 10 -bS - > "$output_file"
done

# Step 2: Remove "_" from the first column
for file in "$input_dir"/*.sam; do
    input_file="${file%.sam}_modified_reads.bam"
    output_file="${file%.sam}_modified_reads2.bam"
    samtools view -@ 14 -h "$input_file" | awk -F'\t' -v OFS='\t' '{gsub("_", ":", $1); print}' | samtools view -@ 10 -bS - > "$output_file"
done

# Change ":" to "#"
for file in "$input_dir"/*.sam; do
    input_file="${file%.sam}_modified_reads2.bam"
    output_file="${file%.sam}_modified_reads3.bam"
    samtools view -@ 14 -h "$input_file" | sed 's/^\([^\t]*\):\([^\t]*\):\([^\t]*\):\([^\t]*\):\([^\t]*\):\([^\t]*\):\([^\t]*\):\([^\t]*\):\([^\t]*\)\t/\1:\2:\3:\4:\5:\6\7:\8:\9\t/g' | samtools view -@ 10 -b - > "$output_file"
done

# Sort
for file in "$input_dir"/*.sam; do
    input_file="${file%.sam}_modified_reads3.bam"
    output_file="${file%.sam}_sorted.bam"
    samtools sort -@ 14 "$input_file" > "$output_file"
done

# Mark duplicates
for file in "$input_dir"/*.sam; do
    input_file="${file%.sam}_sorted.bam"
    output_file="${file%.sam}_noDup.bam"
    samtools markdup -@ 14 -r --barcode-name "$input_file" "$output_file"
done

# Convert BAM to SAM
for file in "$input_dir"/*.sam; do
    input_file="${file%.sam}_noDup.bam"
    output_file="${file%.sam}_noDup.sam"
    samtools view -@ 14 -h "$input_file" -o "$output_file"
done

# Remove intermediate files
rm "$input_dir"/*_modified_reads*.bam

