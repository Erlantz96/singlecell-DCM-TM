#!/bin/bash
# Step 6: Count raw reads per cell
# POOL / POOL_TAG read from environment (set by orchestrator) or default below.

POOL="${POOL:-Pool_44}"
POOL_TAG="${POOL_TAG:-pool_44}"
export POOL POOL_TAG

STATS="${POOL}/RNA/stats_RNA_U87_${POOL_TAG}.txt"

# Start a fresh stats file (overwrites any previous run)
: > "$STATS"

# Define the per-cell function
function run_count_reads() {
	local INPUT=${1}
	local BASENAME=${INPUT##*/}
	local INPUT_name=${BASENAME%R1_extracted.fastq.gz}
	count=$(zcat ${INPUT} | awk 'NR%4==2' | grep -c '')
	echo -e "${BASENAME}\t${count}" >> "${POOL}/RNA/stats_RNA_U87_${POOL_TAG}.txt"
}

# Export for GNU Parallel
export -f run_count_reads

# Run in parallel
parallel -j 15 run_count_reads {1} ::: $(ls ${POOL}/RNA/extracted/*R1_extracted.fastq.gz)
