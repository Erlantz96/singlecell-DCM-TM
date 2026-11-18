#!/usr/bin/env Rscript

# U87 UMI Count Preprocessing Pipeline
# Usage: Rscript 7_make_matrix.R [--umicounts-dir path] [--metadata file] [--output-dir path]
# GTF for annotation is read from the GTF_FILE environment variable (set by the
# orchestrator / wrapper), falling back to the human GENCODE default below.

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(dplyr)
}))

# Parse command-line arguments (defaults are placeholders; the wrapper passes real paths)
args <- commandArgs(trailingOnly = TRUE)
umicounts_dir <- "/mnt/appl/Erlantz/scNMT/Pool_44/RNA/extracted/Aligned/featurecounts/umicounts/"
metadata_file <- "/mnt/appl/Erlantz/scNMT/Pool_44/metadata_U87_pool_44.csv"
output_dir    <- "/mnt/appl/Erlantz/scNMT/Pool_44/"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--umicounts-dir" && i < length(args)) {
    umicounts_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--metadata" && i < length(args)) {
    metadata_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output-dir" && i < length(args)) {
    output_dir <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

# GTF comes from the environment (set by the pipeline); fallback = human
gtf_file <- Sys.getenv("GTF_FILE", "/mnt/appl/Erlantz/reference/gencode.v45.annotation_human.gtf")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

log_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

log_msg("Starting U87 UMI count preprocessing pipeline")
log_msg(sprintf("UMI counts directory: %s", umicounts_dir))
log_msg(sprintf("Metadata file: %s", metadata_file))
log_msg(sprintf("Output directory: %s", output_dir))
log_msg(sprintf("GTF file: %s", gtf_file))

# ============================================================================
# STEP 1: Read TSV files from umicounts directory
# ============================================================================
log_msg("Step 1: Reading UMI count files")

if (!dir.exists(umicounts_dir)) {
  log_msg(sprintf("Directory not found: %s", umicounts_dir), "ERROR"); quit(status = 1)
}

file_list_U87 <- list.files(path = umicounts_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(file_list_U87) == 0) {
  log_msg(sprintf("No TSV files found in %s", umicounts_dir), "ERROR"); quit(status = 1)
}
log_msg(sprintf("Found %d TSV files", length(file_list_U87)))

data_frame_U87 <- lapply(file_list_U87, function(file) {
  tryCatch(
    read.table(file, header = TRUE, sep = "\t"),
    error = function(e) {
      log_msg(sprintf("Error reading %s: %s", file, e$message), "WARNING"); NULL
    }
  )
})
data_frame_U87 <- data_frame_U87[!sapply(data_frame_U87, is.null)]
if (length(data_frame_U87) == 0) { log_msg("All files failed to read", "ERROR"); quit(status = 1) }
log_msg(sprintf("Successfully read %d files", length(data_frame_U87)))

# ============================================================================
# STEP 2: Rename second column based on filenames
# ============================================================================
log_msg("Step 2: Renaming count columns")
for (i in seq_along(data_frame_U87)) {
  new_colname <- gsub("\\.tsv$", "", basename(file_list_U87[i]))
  names(data_frame_U87[[i]])[2] <- new_colname
}

# ============================================================================
# STEP 3: Merge data frames
# ============================================================================
log_msg("Step 3: Merging count tables")
merged_data_U87 <- Reduce(
  function(x, y) merge(x, y, by = names(data_frame_U87[[1]])[1], all = TRUE),
  data_frame_U87
)
merged_data_U87[is.na(merged_data_U87)] <- 0
log_msg(sprintf("Merged table dimensions: %d genes x %d samples",
                nrow(merged_data_U87), ncol(merged_data_U87) - 1))

# ============================================================================
# STEP 4: Set gene names as row names
# ============================================================================
rownames(merged_data_U87) <- merged_data_U87[[1]]
merged_data_U87[[1]] <- NULL

# ============================================================================
# STEP 5: Clean column names
# ============================================================================
colnames(merged_data_U87) <- gsub("etc_S[0-9]+_L001\\.align\\.counts$", "", colnames(merged_data_U87))

# ============================================================================
# STEP 6: Load and match metadata
# ============================================================================
log_msg("Step 6: Loading metadata and matching samples")
if (!file.exists(metadata_file)) {
  log_msg(sprintf("Metadata file not found: %s", metadata_file), "ERROR"); quit(status = 1)
}
metadata_U87 <- read.csv(metadata_file, header = TRUE, sep = ",")
log_msg(sprintf("Loaded metadata for %d samples", nrow(metadata_U87)))

merged_data_U87_2 <- merged_data_U87[, match(metadata_U87$Sample, colnames(merged_data_U87))]
if (any(is.na(colnames(merged_data_U87_2))))
  log_msg("Warning: Some metadata samples not found in count matrix", "WARNING")
if (!all(colnames(merged_data_U87_2) == metadata_U87$Sample))
  log_msg("Warning: Column order mismatch after metadata matching", "WARNING")

# ============================================================================
# STEP 6.5: Write merged counts before annotation
# ============================================================================
merged_counts_output <- file.path(output_dir, "U87_counts_merged.csv")
write.csv(merged_data_U87_2, file = merged_counts_output)
log_msg(sprintf("Wrote merged counts to: %s", merged_counts_output))

# ============================================================================
# STEP 7: Gene annotation from the GTF (no biomaRt needed)
# ============================================================================
log_msg(sprintf("Loaded merged counts: %d genes x %d samples",
                nrow(merged_data_U87_2), ncol(merged_data_U87_2)))
cat("First few row names:", head(row.names(merged_data_U87_2)), "\n")
log_msg("Step 7: Annotating genes with GTF")

gtf <- read.delim(gtf_file, header = FALSE, comment.char = "#",
                  sep = "\t", quote = "", stringsAsFactors = FALSE)
gtf <- gtf[gtf$V9 != "", ]

attr      <- gtf$V9
gene_id   <- sub('.*gene_id "([^"]+)".*',   "\\1", attr)
gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", attr)
gene_name[gene_name == attr] <- NA   # rows without a gene_name -> NA

map <- unique(data.frame(ensembl = sub("\\..*", "", gene_id),
                         symbol  = gene_name,
                         stringsAsFactors = FALSE))
map <- map[!is.na(map$symbol) & map$symbol != "", ]

genes  <- row.names(merged_data_U87_2)
genes2 <- sub("\\..*", "", genes)    # strip Ensembl version suffix

sym <- map$symbol[match(genes2, map$ensembl)]
log_msg(sprintf("Mapped %d of %d genes to symbols", sum(!is.na(sym)), length(sym)))
merged_data_U87_2$gene <- sym

# ============================================================================
# STEP 8: Aggregate by gene symbol
# ============================================================================
log_msg("Step 8: Aggregating counts by gene")
merged_data_U87_3 <- merged_data_U87_2 %>%
  group_by(gene) %>%
  summarise(across(everything(), sum), .groups = "drop")
merged_data_U87_3 <- data.frame(merged_data_U87_3)
log_msg(sprintf("After aggregation: %d unique genes", nrow(merged_data_U87_3)))

# ============================================================================
# STEP 9: Remove NA and prepare for output
# ============================================================================
log_msg("Step 9: Cleaning and preparing output")
merged_data_U87_4 <- na.omit(merged_data_U87_3)
log_msg(sprintf("After removing NAs: %d genes", nrow(merged_data_U87_4)))
rownames(merged_data_U87_4) <- merged_data_U87_4$gene
merged_data_U87_4$gene <- NULL

# ============================================================================
# STEP 10: Write output files
# ============================================================================
log_msg("Step 10: Writing output files")
umi_output <- file.path(output_dir, "U87_UMI.csv")
write.csv(colSums(merged_data_U87_4), umi_output)
log_msg(sprintf("Wrote UMI sums to: %s", umi_output))

counts_output <- file.path(output_dir, "U87_counts.csv")
write.csv(merged_data_U87_4, file = counts_output)
log_msg(sprintf("Wrote counts to: %s (%d genes x %d samples)",
                counts_output, nrow(merged_data_U87_4), ncol(merged_data_U87_4)))

log_msg("Pipeline completed successfully")
