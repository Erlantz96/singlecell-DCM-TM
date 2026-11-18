#!/bin/bash
# Step 3: featureCounts + UMI count
# POOL and GTF_FILE are read from the environment (set by the orchestrator),
# or default below so the script also runs standalone.

POOL="${POOL:-Pool_44}"
GTF_FILE="${GTF_FILE:-/mnt/appl/Erlantz/reference/gencode.vM10.annotation.gtf}"
export POOL GTF_FILE

# Make sure the output directories exist
mkdir -p "${POOL}/RNA/extracted/Aligned/featurecounts/umicounts"

# Define the per-cell function
function run_featurecounts() {
	local INPUT=${1}
	local BASENAME=${INPUT##*/}
	local INPUT_name=${BASENAME%.sam}
	local OUTPUT_bam=${INPUT_name}.bam
	local GTF=${2}
	echo "${BASENAME}"
	samtools view -bS "${POOL}/RNA/extracted/Aligned/${BASENAME}" > "${POOL}/RNA/extracted/Aligned/${INPUT_name}.bam"
	samtools sort "${POOL}/RNA/extracted/Aligned/${INPUT_name}.bam" -o "${POOL}/RNA/extracted/Aligned/${INPUT_name}sorted.bam"
	samtools index "${POOL}/RNA/extracted/Aligned/${INPUT_name}sorted.bam"
	featureCounts -a "${GTF}" \
	              -o "${POOL}/RNA/extracted/Aligned/featurecounts/${INPUT_name}gene_Assigned" \
	              -R BAM "${POOL}/RNA/extracted/Aligned/${INPUT_name}sorted.bam" \
	              -T 4 &>> "${POOL}/RNA/extracted/Aligned/featurecounts/featurecounts_files.log"
	samtools sort "${POOL}/RNA/extracted/Aligned/featurecounts/${INPUT_name}sorted.bam.featureCounts.bam" -o "${POOL}/RNA/extracted/Aligned/featurecounts/${INPUT_name}assigned_sorted.bam"
	samtools index "${POOL}/RNA/extracted/Aligned/featurecounts/${INPUT_name}assigned_sorted.bam"
	umi_tools count --per-gene --gene-tag=XT --wide-format-cell-counts \
		-I "${POOL}/RNA/extracted/Aligned/featurecounts/${INPUT_name}assigned_sorted.bam" \
		-S "${POOL}/RNA/extracted/Aligned/featurecounts/umicounts/${INPUT_name}.counts.tsv" &>> "${POOL}/RNA/extracted/Aligned/featurecounts/umicounts/umicounts_files.log"
}

# Export for GNU Parallel
export -f run_featurecounts

# Run in parallel; the GTF is passed as the second argument
parallel -j 15 run_featurecounts {1} {2} ::: $(ls ${POOL}/RNA/extracted/Aligned/*.sam) ::: "${GTF_FILE}"
