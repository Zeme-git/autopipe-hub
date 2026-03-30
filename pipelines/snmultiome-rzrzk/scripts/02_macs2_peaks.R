#!/usr/bin/env Rscript
# =============================================================================
#  02_macs2_peaks.R
#  Call peaks with MACS2, remove blacklist regions, create macs2 ChromatinAssay
#  Extracted from: 1_RZRZK_basic_analysis.R (call_macs2_counts, add_peak_to_object)
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
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

cfg <- read_yaml(config_file)

plan("multicore", workers = cfg$n_workers)
options(future.globals.maxSize = cfg$future_max_size_gb * 1024^3)
set.seed(cfg$seed)

cat("=== Step 2: MACS2 peak calling ===\n")

# ── Set annotation ────────────────────────────────────────────────────────────
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

# ── Helper: Call MACS2 peaks, remove blacklist, quantify ──────────────────────
call_macs2_counts <- function(obj, species, macs2_path) {
  peaks <- CallPeaks(obj, macs2.path = macs2_path, verbose = TRUE)
  if (species == "mouse") {
    blacklist <- blacklist_mm10
  } else {
    blacklist <- blacklist_hg38_unified
  }
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist, invert = TRUE)
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(obj),
    features  = peaks,
    cells     = colnames(obj)
  )
  return(macs2_counts)
}

# ── Helper: Add MACS2 ChromatinAssay to object ───────────────────────────────
add_peak_to_object <- function(obj, species, frag_file, macs2_path, annotations) {
  obj[["macs2"]] <- CreateChromatinAssay(
    counts     = call_macs2_counts(obj, species, macs2_path),
    fragments  = frag_file,
    annotation = annotations
  )
  return(obj)
}

# ── Load QC-filtered objects ──────────────────────────────────────────────────
wild   <- readRDS("/output/01_qc/wild_qc.rds")
mutant <- readRDS("/output/01_qc/mutant_qc.rds")

cat("  Calling MACS2 peaks for wild-type...\n")
wild <- add_peak_to_object(wild, cfg$species, cfg$samples$wild$fragment_file,
                            cfg$macs2_path, annotations)

cat("  Calling MACS2 peaks for mutant...\n")
mutant <- add_peak_to_object(mutant, cfg$species, cfg$samples$mutant$fragment_file,
                              cfg$macs2_path, annotations)

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create("/output/02_peaks", showWarnings = FALSE, recursive = TRUE)
saveRDS(wild,   "/output/02_peaks/wild_peaks.rds")
saveRDS(mutant, "/output/02_peaks/mutant_peaks.rds")

cat("=== Step 2 complete ===\n")