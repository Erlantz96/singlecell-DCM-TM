#!/bin/bash
#deeptools environment in the server
fastq=/mnt/appl/Erlantz/scNMT/Pool_40/gDNA/*_R1_001.fastq
for i in $fastq; do
    filename=$(basename "$i")
    trimmed_name="${filename%_R1_TRIM.fastq}_R1_TRIM.fastq"
	cutadapt --cores 5 -m 20 -g AGTTCTACAGTCCGACGATC -o /mnt/appl/Erlantz/scNMT/Pool_40/gDNA/2_Cell_Annotation_input_mouse/${trimmed_name} $i;
done

