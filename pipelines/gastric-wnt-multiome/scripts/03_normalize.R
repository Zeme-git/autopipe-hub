#!/usr/bin/env Rscript
# Step 3: Normalization & Dimensional Reduction

options(Seurat.object.assay.version = "v3")

library(future); library(Seurat); library(Signac); library(yaml)
cfg <- yaml.load_file("/pipeline/config.yaml")
plan("multicore", workers = cfg$threads); options(future.globals.maxSize = cfg$max_memory_gb * 1024^3); set.seed(1234)
out_dir <- cfg$output_dir
wild <- readRDS(file.path(out_dir,"rds","wild_macs2.rds")); mutant <- readRDS(file.path(out_dir,"rds","mutant_macs2.rds"))

process_expression <- function(obj) {
    obj <- SCTransform(obj)
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = cfg$rna$nfeatures)
    return(obj)
}

process_accessibility <- function(obj) {
    for (a in c("ATAC","macs2")) {
        DefaultAssay(obj) <- a
        obj <- FindTopFeatures(obj, min.cutoff = "q0")
        obj <- RunTFIDF(obj)
        obj <- RunSVD(obj)
    }
    return(obj)
}

cat("Normalizing Wild-type ...\n"); wild <- process_expression(wild); wild <- process_accessibility(wild)
cat("Normalizing Mutant ...\n"); mutant <- process_expression(mutant); mutant <- process_accessibility(mutant)
saveRDS(wild, file.path(out_dir,"rds","wild_normalized.rds")); saveRDS(mutant, file.path(out_dir,"rds","mutant_normalized.rds"))
cat("Step 3 complete.\n")