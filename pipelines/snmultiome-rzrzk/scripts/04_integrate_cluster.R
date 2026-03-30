#!/usr/bin/env Rscript
# =============================================================================
#  04_integrate_cluster.R
#  SCT-based integration → PCA → UMAP → Neighbor finding → Clustering
#  Extracted from: 1_RZRZK_basic_analysis.R (integrate_based_on_rna,
#                  generate_umap_integrated)
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

cat("=== Step 4: Integration & Clustering ===\n")

# ── Load normalized objects ───────────────────────────────────────────────────
wild   <- readRDS("/output/03_norm/wild_norm.rds")
mutant <- readRDS("/output/03_norm/mutant_norm.rds")

# ── Integration based on SCT-normalized RNA ───────────────────────────────────
cat("  Finding integration features...\n")
DefaultAssay(wild)   <- "SCT"
DefaultAssay(mutant) <- "SCT"
obj.list <- list(wild, mutant)

features    <- SelectIntegrationFeatures(obj.list)
obj.list    <- PrepSCTIntegration(obj.list, anchor.features = features)
obj.anchors <- FindIntegrationAnchors(
  object.list          = obj.list,
  anchor.features      = features,
  normalization.method = "SCT"
)
obj.combined <- IntegrateData(
  anchorset            = obj.anchors,
  normalization.method = "SCT"
)

cat("  Integration complete. Total cells:", ncol(obj.combined), "\n")

# ── Dimensional reduction & clustering on integrated assay ────────────────────
cat("  Running PCA, UMAP, and clustering...\n")
DefaultAssay(obj.combined) <- "integrated"
obj.combined <- ScaleData(obj.combined, verbose = FALSE)
obj.combined <- RunPCA(obj.combined, npcs = cfg$clustering$npcs, verbose = FALSE)
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:cfg$clustering$dims)
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:cfg$clustering$dims)

# Initial clustering at resolution 0.3 (produces 6 clusters: 0-5)
obj.combined <- FindClusters(
  obj.combined,
  resolution = cfg$clustering$resolution_initial,
  algorithm  = cfg$clustering$algorithm
)

# Create combined cluster + genotype label
obj.combined@meta.data$cl_genotype <- paste0(
  obj.combined@meta.data$integrated_snn_res.0.3, "_",
  obj.combined@meta.data$kras_genotype
)

cat("  Cluster counts at resolution", cfg$clustering$resolution_initial, ":\n")
print(table(obj.combined@meta.data$integrated_snn_res.0.3))

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create("/output/04_cluster", showWarnings = FALSE, recursive = TRUE)
saveRDS(obj.combined, "/output/04_cluster/combined_clustered.rds")

cat("=== Step 4 complete ===\n")