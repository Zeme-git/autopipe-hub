#!/usr/bin/env Rscript
# =============================================================================
#  03_normalize.R
#  RNA: SCTransform + LogNormalize + FindVariableFeatures
#  ATAC: FindTopFeatures + RunTFIDF + RunSVD (for both ATAC and macs2 assays)
#  Extracted from: 1_RZRZK_basic_analysis.R (process_expression, process_accessibility)
#  NOTE: Uses Seurat v3 assay version for v4-style object compatibility
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]

# ── Force Seurat v3/v4 assay structure in Seurat v5 ──────────────────────────
options(Seurat.object.assay.version = "v3")

suppressPackageStartupMessages({
  library(yaml)
  library(future)
  library(Seurat)
  library(Signac)
})

cfg <- read_yaml(config_file)

plan("multicore", workers = cfg$n_workers)
options(future.globals.maxSize = cfg$future_max_size_gb * 1024^3)
set.seed(cfg$seed)

cat("=== Step 3: Normalization ===\n")

# ── Helper: Process gene expression (SCTransform + LogNormalize) ──────────────
process_expression <- function(obj, nfeatures = 2000) {
  obj <- SCTransform(obj)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  return(obj)
}

# ── Helper: Process DNA accessibility ────────────────────────────────────────
process_accessibility <- function(obj, min_cutoff = "q0") {
  DefaultAssay(obj) <- "ATAC"
  obj <- FindTopFeatures(obj, min.cutoff = min_cutoff)
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)
  DefaultAssay(obj) <- "macs2"
  obj <- FindTopFeatures(obj, min.cutoff = min_cutoff)
  obj <- RunTFIDF(obj)
  obj <- RunSVD(obj)
  return(obj)
}

# ── Load objects with MACS2 peaks ─────────────────────────────────────────────
wild   <- readRDS("/output/02_peaks/wild_peaks.rds")
mutant <- readRDS("/output/02_peaks/mutant_peaks.rds")

cat("  Processing RNA expression (wild)...\n")
wild <- process_expression(wild, nfeatures = cfg$processing$vst_nfeatures)

cat("  Processing RNA expression (mutant)...\n")
mutant <- process_expression(mutant, nfeatures = cfg$processing$vst_nfeatures)

cat("  Processing ATAC accessibility (wild)...\n")
wild <- process_accessibility(wild, min_cutoff = cfg$processing$top_features_min_cutoff)

cat("  Processing ATAC accessibility (mutant)...\n")
mutant <- process_accessibility(mutant, min_cutoff = cfg$processing$top_features_min_cutoff)

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create("/output/03_norm", showWarnings = FALSE, recursive = TRUE)
saveRDS(wild,   "/output/03_norm/wild_norm.rds")
saveRDS(mutant, "/output/03_norm/mutant_norm.rds")

cat("=== Step 3 complete ===\n")