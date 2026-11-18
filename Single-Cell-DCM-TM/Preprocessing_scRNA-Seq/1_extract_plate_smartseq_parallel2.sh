#!/bin/bash
# Step 1: UMI extraction (SMART-seq plate)
# POOL is read from the environment (set by the orchestrator) or defaults below.

POOL="${POOL:-Pool_44}"
export POOL

# Make sure the output directory exists
mkdir -p "${POOL}/RNA/extracted"

# Define the per-cell extraction function
function run_umiextract() {
	local INPUT=${1}
	local BASENAME=${INPUT##*/}
	local INPUT_name=${BASENAME%_R1_001.fastq.gz}
	echo "${BASENAME}"
	umi_tools extract --bc-pattern=CCCCCCCCCCCNNNNNNNNNNNN \
		--stdin "${POOL}/RNA/${INPUT_name}_R1_001.fastq.gz" \
		--stdout "${POOL}/RNA/extracted/${INPUT_name}_R1_extracted.fastq.gz" \
		--read2-in "${POOL}/RNA/${INPUT_name}_R2_001.fastq.gz" \
		--read2-out="${POOL}/RNA/extracted/${INPUT_name}_R2_extracted.fastq.gz"
}

# Export for GNU Parallel
export -f run_umiextract

# Run in parallel
parallel -j 15 run_umiextract {1} ::: $(ls ${POOL}/RNA/*_R1_001.fastq.gz)
