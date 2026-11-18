#!/bin/bash
# Step 7: Build the count matrix (R), with gene annotation from the GTF.
# POOL / POOL_TAG / GTF_FILE read from environment (set by orchestrator) or default below.
set -e

POOL="${POOL:-Pool_44}"
POOL_TAG="${POOL_TAG:-pool_44}"
GTF_FILE="${GTF_FILE:-/mnt/appl/Erlantz/reference/gencode.v45.annotation_human.gtf}"
BASE_DIR="${BASE_DIR:-/mnt/appl/Erlantz/scNMT}"
export POOL POOL_TAG GTF_FILE

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/7_make_matrix.R"

# Paths derived from POOL / POOL_TAG
UMICOUNTS_DIR="${BASE_DIR}/${POOL}/RNA/extracted/Aligned/featurecounts/umicounts/"
METADATA_FILE="${BASE_DIR}/${POOL}/metadata_U87_${POOL_TAG}.csv"
OUTPUT_DIR="${BASE_DIR}/${POOL}/"

# Colors
RED='\033[0;31m'; GREEN='\033[0;32m'; NC='\033[0m'

echo -e "${GREEN}=== U87 UMI Count Preprocessing (R) — ${POOL} ===${NC}\n"

# Validate inputs
[ -f "$R_SCRIPT" ]        || { echo -e "${RED}R script not found: $R_SCRIPT${NC}"; exit 1; }
[ -d "$UMICOUNTS_DIR" ]   || { echo -e "${RED}UMI counts dir not found: $UMICOUNTS_DIR${NC}"; exit 1; }
[ -f "$METADATA_FILE" ]   || { echo -e "${RED}Metadata file not found: $METADATA_FILE${NC}"; exit 1; }
[ -f "$GTF_FILE" ]        || { echo -e "${RED}GTF file not found: $GTF_FILE${NC}"; exit 1; }
command -v Rscript &>/dev/null || { echo -e "${RED}Rscript not found.${NC}"; exit 1; }

echo "  UMI counts: $UMICOUNTS_DIR"
echo "  Metadata:   $METADATA_FILE"
echo "  GTF:        $GTF_FILE"
echo "  Output:     $OUTPUT_DIR"
echo ""

# GTF_FILE is exported, so the R script picks it up via Sys.getenv("GTF_FILE")
Rscript "$R_SCRIPT" \
  --umicounts-dir "$UMICOUNTS_DIR" \
  --metadata "$METADATA_FILE" \
  --output-dir "$OUTPUT_DIR"

if [ $? -eq 0 ]; then
  echo -e "\n${GREEN}Matrix build completed.${NC}"
  echo "  - ${OUTPUT_DIR}U87_counts_merged.csv"
  echo "  - ${OUTPUT_DIR}U87_UMI.csv"
  echo "  - ${OUTPUT_DIR}U87_counts.csv"
else
  echo -e "\n${RED}Matrix build failed.${NC}"
  exit 1
fi
