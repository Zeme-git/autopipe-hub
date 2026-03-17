#!/usr/bin/env Rscript
# =============================================================================
#  Step 3: Normalization & Dimensional Reduction
# =============================================================================
#  RNA:  SCTransform (for integration) + LogNormalize + FindVariableFeatures
#  ATAC: FindTopFeatures → TF-IDF → SVD  (on both "ATAC" and "macs2" assays)
# =============================================================================

library(future)
library(Seurat)
library(Signac)
library(yaml)

cfg <- yaml.load_file("/pipeline/config.yaml")
plan("multicore", workers = cfg$threads)
options(future.globals.maxSize = cfg$max_memory_gb * 1024^3)
set.seed(1234)

out_dir <- cfg$output_dir
wild   <- readRDS(file.path(out_dir, "rds", "wild_macs2.rds"))
mutant <- readRDS(file.path(out_dir, "rds", "mutant_macs2.rds"))

# ---- RNA processing ---------------------------------------------------------
process_expression <- function(obj) {
    # SCTransform: regularized NB regression → "SCT" assay
    obj <- SCTransform(obj)
    # Standard log-normalization → stored in "RNA" assay (for DE later)
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = cfg$rna$nfeatures)
    return(obj)
}

# ---- ATAC processing --------------------------------------------------------
process_accessibility <- function(obj) {
    for (assay_name in c("ATAC", "macs2")) {
        DefaultAssay(obj) <- assay_name
        obj <- FindTopFeatures(obj, min.cutoff = "q0")
        obj <- RunTFIDF(obj)
        obj <- RunSVD(obj)
    }
    return(obj)
}

cat("Normalizing Wild-type RNA ...\n")
wild   <- process_expression(wild)
cat("Normalizing Mutant RNA ...\n")
mutant <- process_expression(mutant)

cat("Processing Wild-type ATAC ...\n")
wild   <- process_accessibility(wild)
cat("Processing Mutant ATAC ...\n")
mutant <- process_accessibility(mutant)

# ---- Save --------------------------------------------------------------------
saveRDS(wild,   file.path(out_dir, "rds", "wild_normalized.rds"))
saveRDS(mutant, file.path(out_dir, "rds", "mutant_normalized.rds"))
cat("Step 3 complete.\n")