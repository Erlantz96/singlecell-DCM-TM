#!/bin/bash
# Step 4: Dedup + unique-gene stats
# POOL / POOL_TAG read from environment (set by orchestrator) or default below.

POOL="${POOL:-Pool_44}"
POOL_TAG="${POOL_TAG:-pool_44}"
export POOL POOL_TAG

# Make sure output dirs exist
mkdir -p "${POOL}/RNA/extracted/Aligned/featurecounts/dedup"

# Header for the gene-count summary table
echo -e "File\tUnique_Gene_Count" > "${POOL}/RNA/extracted/Aligned/featurecounts/gene_counts_U87_${POOL_TAG}.tsv"

# Define the per-cell function
function run_genecounts() {
	local INPUT=${1}
	local BASENAME=${INPUT##*/}
	local INPUT_name=${BASENAME%assigned_sorted.bam}
	echo "${BASENAME}"
	# NOTE: original used ${name} in the log filename, which was never defined.
	# Fixed to ${INPUT_name} so each cell gets its own dedup log.
	umi_tools dedup --per-gene --gene-tag=XT \
		--stdin "${POOL}/RNA/extracted/Aligned/featurecounts/${BASENAME}" \
		--log="${POOL}/RNA/extracted/Aligned/featurecounts/dedup/${INPUT_name}.log" \
		> "${POOL}/RNA/extracted/Aligned/featurecounts/dedup/${INPUT_name}.dedup.bam"
	unique_gene_count=$(samtools view -h "${POOL}/RNA/extracted/Aligned/featurecounts/dedup/${INPUT_name}.dedup.bam" | grep -o -P 'XT:Z:\S+' | cut -f3 -d':' | sort -u | wc -l)
	echo -e "${INPUT_name}\t${unique_gene_count}" >> "${POOL}/RNA/extracted/Aligned/featurecounts/gene_counts_U87_${POOL_TAG}.tsv"
}

# Export for GNU Parallel
export -f run_genecounts

# Run in parallel
parallel -j 15 run_genecounts {1} ::: $(ls ${POOL}/RNA/extracted/Aligned/featurecounts/*assigned_sorted.bam)
