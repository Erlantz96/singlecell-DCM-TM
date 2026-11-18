#!/bin/bash
# =============================================================================
# Master RNA preprocessing pipeline
# Runs steps 1-7 in order. Each step runs ONLY if the previous one
#   (a) exited with status 0, AND
#   (b) produced a non-empty expected output.
# Stops immediately on the first failure.
#
# Usage:
#   bash run_full_pipeline.sh human          # human reference
#   bash run_full_pipeline.sh mouse          # mouse reference
#   bash run_full_pipeline.sh human Pool_42  # human + a specific pool
# =============================================================================

# --- Configuration ----------------------------------------------------------
BASE_DIR="/mnt/appl/Erlantz/scNMT"        # dir that CONTAINS the Pool_XX folder
GENOME="${1:-human}"                       # 1st arg: human | mouse
POOL="${2:-Pool_44}"                       # 2nd arg: pool directory name
POOL_TAG="$(echo "$POOL" | tr 'A-Z' 'a-z')"  # filename tag, e.g. Pool_44 -> pool_44

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="${BASE_DIR}/pipeline_${POOL}_$(date +%Y%m%d_%H%M%S).log"

export BASE_DIR POOL POOL_TAG

# --- Colors -----------------------------------------------------------------
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'

log() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"; }

# --- Genome / reference selection -------------------------------------------
case "$GENOME" in
  human)
    BOWTIE2_INDEX="/mnt/appl/Erlantz/reference/bowtie2/human/hg38"
    GTF_FILE="/mnt/appl/Erlantz/reference/gencode.v45.annotation_human.gtf"
    ;;
  mouse)
    BOWTIE2_INDEX="/mnt/appl/Erlantz/reference/bowtie2/mouse/mm10"
    GTF_FILE="/mnt/appl/Erlantz/reference/gencode.vM10.annotation.gtf"
    ;;
  *)
    echo -e "${RED}Unknown genome '${GENOME}'. Use 'human' or 'mouse'.${NC}"; exit 1 ;;
esac
export BOWTIE2_INDEX GTF_FILE

# Verify references exist before doing anything expensive
if ! ls "${BOWTIE2_INDEX}".*bt2* >/dev/null 2>&1; then
  echo -e "${RED}Bowtie2 index not found at prefix: ${BOWTIE2_INDEX}${NC}"; exit 1
fi
[ -f "$GTF_FILE" ] || { echo -e "${RED}GTF file not found: ${GTF_FILE}${NC}"; exit 1; }

# --- Output-check helper ----------------------------------------------------
check_output() {
  local pattern="$1" found=0
  for f in $pattern; do [ -s "$f" ] && found=$((found + 1)); done
  if [ "$found" -gt 0 ]; then
    log "  ${GREEN}OK${NC} $found non-empty file(s) matching: $pattern"; return 0
  else
    log "  ${RED}FAIL${NC} no non-empty files matching: $pattern"; return 1
  fi
}

# --- Run a single step ------------------------------------------------------
run_step() {
  local num="$1" script="$2" desc="$3" out_glob="$4"
  local path="${SCRIPT_DIR}/${script}"

  log ""
  log "${YELLOW}=== STEP ${num}: ${desc} ===${NC}"
  log "  Script: ${script}"
  [ -f "$path" ] || { log "  ${RED}Script not found: $path${NC}"; exit 1; }

  local start=$(date +%s)
  ( cd "$BASE_DIR" && bash "$path" ) 2>&1 | tee -a "$LOG_FILE"
  local code=${PIPESTATUS[0]}
  local elapsed=$(( $(date +%s) - start ))

  if [ "$code" -ne 0 ]; then
    log "  ${RED}Step ${num} exited with status ${code} after ${elapsed}s${NC}"
    log "${RED}Pipeline stopped at step ${num}.${NC}"; exit "$code"
  fi
  log "  Exit status 0 (${elapsed}s)"

  if [ -n "$out_glob" ]; then
    check_output "$out_glob" || { log "${RED}Stopped: step ${num} produced no valid output.${NC}"; exit 1; }
  fi
  log "  ${GREEN}Step ${num} complete.${NC}"
}

# =============================================================================
# Pipeline
# =============================================================================
log "${GREEN}########## RNA PREPROCESSING -- ${POOL} (${GENOME}) ##########${NC}"
log "Bowtie2 index: ${BOWTIE2_INDEX}"
log "GTF file:      ${GTF_FILE}"
log "Base dir:      ${BASE_DIR}"
log "Pool tag:      ${POOL_TAG}"
log "Log file:      ${LOG_FILE}"

R="${BASE_DIR}/${POOL}/RNA"

run_step 1 "1_extract_plate_smartseq_parallel2.sh" \
  "UMI extraction" \
  "${R}/extracted/*_R1_extracted.fastq.gz"

run_step 2 "2_bowtie2_2reads_parallel.sh" \
  "Bowtie2 alignment" \
  "${R}/extracted/Aligned/*.align.sam"

run_step 3 "3_feature_counts_parallel.sh" \
  "featureCounts + UMI count" \
  "${R}/extracted/Aligned/featurecounts/umicounts/*.counts.tsv"

run_step 4 "4_dedup_gene_stats_from_bam_parallel.sh" \
  "Dedup + gene stats" \
  "${R}/extracted/Aligned/featurecounts/gene_counts_U87_${POOL_TAG}.tsv"

run_step 5 "5_stats_assigned_reads_parallel.sh" \
  "Assigned-reads summary" \
  "${R}/combined_summary_U87_${POOL_TAG}.tsv"

run_step 6 "6_count_raw_aligned_reads.sh" \
  "Raw read counts" \
  "${R}/stats_RNA_U87_${POOL_TAG}.txt"

run_step 7 "7_run_make_matrix.sh" \
  "Build count matrix (R)" \
  "${BASE_DIR}/${POOL}/U87_counts.csv"

log ""
log "${GREEN}########## PIPELINE COMPLETED (${POOL}, ${GENOME}) ##########${NC}"
log "Full log: ${LOG_FILE}"
