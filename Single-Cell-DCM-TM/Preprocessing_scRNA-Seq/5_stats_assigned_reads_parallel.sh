#!/bin/bash
# Step 5: featureCounts assignment summary
# POOL / POOL_TAG read from environment (set by orchestrator) or default below.

POOL="${POOL:-Pool_44}"
POOL_TAG="${POOL_TAG:-pool_44}"
export POOL POOL_TAG

SUMMARY="${POOL}/RNA/combined_summary_U87_${POOL_TAG}.tsv"

# Define the per-file function
function run_assignedsummary() {
	local INPUT=${1}
	local BASENAME=${INPUT##*/}
	local INPUT_name=${BASENAME%_.aligngene_Assigned.summary}
	echo "${BASENAME}"

	assigned=$(grep "Assigned" ${INPUT} | cut -f2)
	unmapped=$(grep "Unassigned_Unmapped" ${INPUT} | cut -f2)
	read_type=$(grep "Unassigned_Read_Type" ${INPUT} | cut -f2)
	singleton=$(grep "Unassigned_Singleton" ${INPUT} | cut -f2)
	mapping_quality=$(grep "Unassigned_MappingQuality" ${INPUT} | cut -f2)
	chimera=$(grep "Unassigned_Chimera" ${INPUT} | cut -f2)
	fragment_length=$(grep "Unassigned_FragmentLength" ${INPUT} | cut -f2)
	duplicate=$(grep "Unassigned_Duplicate" ${INPUT} | cut -f2)
	multi_mapping=$(grep "Unassigned_MultiMapping" ${INPUT} | cut -f2)
	secondary=$(grep "Unassigned_Secondary" ${INPUT} | cut -f2)
	non_split=$(grep "Unassigned_NonSplit" ${INPUT} | cut -f2)
	no_features=$(grep "Unassigned_NoFeatures" ${INPUT} | cut -f2)
	overlapping_length=$(grep "Unassigned_Overlapping_Length" ${INPUT} | cut -f2)
	ambiguity=$(grep "Unassigned_Ambiguity" ${INPUT} | cut -f2)

	echo -e "${INPUT_name}\t${assigned}\t${unmapped}\t${read_type}\t${singleton}\t${mapping_quality}\t${chimera}\t${fragment_length}\t${duplicate}\t${multi_mapping}\t${secondary}\t${non_split}\t${no_features}\t${overlapping_length}\t${ambiguity}" >> "${POOL}/RNA/combined_summary_U87_${POOL_TAG}.tsv"
}

# Write header (overwrites any previous run)
echo -e "File\tAssigned\tUnassigned_Unmapped\tUnassigned_Read_Type\tUnassigned_Singleton\tUnassigned_MappingQuality\tUnassigned_Chimera\tUnassigned_FragmentLength\tUnassigned_Duplicate\tUnassigned_MultiMapping\tUnassigned_Secondary\tUnassigned_NonSplit\tUnassigned_NoFeatures\tUnassigned_Overlapping_Length\tUnassigned_Ambiguity" > "$SUMMARY"

# Export for GNU Parallel
export -f run_assignedsummary

# Run in parallel
parallel -j 8 run_assignedsummary {1} ::: $(ls ${POOL}/RNA/extracted/Aligned/featurecounts/*_Assigned.summary)
