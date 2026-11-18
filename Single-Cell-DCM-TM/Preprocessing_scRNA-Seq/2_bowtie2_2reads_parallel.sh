#!/bin/bash
# Step 2: Bowtie2 alignment
# POOL and BOWTIE2_INDEX are read from the environment (set by the orchestrator),
# or default below so the script also runs standalone.

POOL="${POOL:-Pool_44}"
BOWTIE2_INDEX="${BOWTIE2_INDEX:-/mnt/appl/Erlantz/reference/bowtie2/mouse/mm10}"
export POOL BOWTIE2_INDEX

# Make sure the output directory exists
mkdir -p "${POOL}/RNA/extracted/Aligned"

# Define the per-cell alignment function
function run_bowtie2() {
	local INPUT=${1}
	local BASENAME=${INPUT##*/}
	local INPUT_name=${BASENAME%_R1_extracted.fastq.gz}
	local OUTPUT=${INPUT_name}.align.sam
	local INDEX=${2}
	bowtie2 -x "${INDEX}" \
		-1 "${POOL}/RNA/extracted/${INPUT_name}_R1_extracted.fastq.gz" \
		-2 "${POOL}/RNA/extracted/${INPUT_name}_R2_extracted.fastq.gz" \
		-S "${POOL}/RNA/extracted/Aligned/${OUTPUT}" &>> "${POOL}/RNA/extracted/Aligned/fastq_files.log"
}

# Export for GNU Parallel
export -f run_bowtie2

# Run in parallel; the index is passed as the second argument
parallel -j 15 run_bowtie2 {1} {2} ::: $(ls ${POOL}/RNA/extracted/*_R1_extracted.fastq.gz) ::: "${BOWTIE2_INDEX}"
