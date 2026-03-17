#!/usr/bin/env Rscript
# =============================================================================
#  Step 2: MACS2 Peak Calling
# =============================================================================
#  Calls peaks using MACS2, removes blacklist regions, quantifies a new
#  FeatureMatrix, and adds a "macs2" ChromatinAssay to each sample.
# =============================================================================

library(future)
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(yaml)

cfg <- yaml.load_file("/pipeline/config.yaml")
plan("multicore", workers = cfg$threads)
options(future.globals.maxSize = cfg$max_memory_gb * 1024^3)
set.seed(1234)

# Rebuild annotation (needed for CreateChromatinAssay)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

# ---- Load QC'd objects -------------------------------------------------------
out_dir <- cfg$output_dir
wild   <- readRDS(file.path(out_dir, "rds", "wild_qc.rds"))
mutant <- readRDS(file.path(out_dir, "rds", "mutant_qc.rds"))

# ---- MACS2 peak calling + new assay -----------------------------------------
add_macs2_assay <- function(obj, species, frag_file) {
    # Call peaks
    peaks <- CallPeaks(obj, macs2.path = cfg$macs2_path, verbose = TRUE)
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

    # Remove blacklist
    bl <- if (species == "mouse") blacklist_mm10 else blacklist_hg38_unified
    peaks <- subsetByOverlaps(x = peaks, ranges = bl, invert = TRUE)

    # Quantify
    macs2_counts <- FeatureMatrix(
        fragments = Fragments(obj),
        features  = peaks,
        cells     = colnames(obj)
    )

    obj[["macs2"]] <- CreateChromatinAssay(
        counts     = macs2_counts,
        fragments  = frag_file,
        annotation = annotations
    )
    return(obj)
}

cat("Calling MACS2 peaks for Wild-type ...\n")
wild   <- add_macs2_assay(wild,   cfg$species, cfg$samples$Wild_type$fragments)

cat("Calling MACS2 peaks for Mutant ...\n")
mutant <- add_macs2_assay(mutant, cfg$species, cfg$samples$Mutant$fragments)

# ---- Save --------------------------------------------------------------------
saveRDS(wild,   file.path(out_dir, "rds", "wild_macs2.rds"))
saveRDS(mutant, file.path(out_dir, "rds", "mutant_macs2.rds"))
cat("Step 2 complete.\n")