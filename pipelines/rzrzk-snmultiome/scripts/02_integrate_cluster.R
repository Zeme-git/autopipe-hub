#!/usr/bin/env Rscript
###############################################################################
#  Step 6-7: Integration + Clustering
#  SCT-based anchor integration → PCA → UMAP → Clustering (res 0.3 & 0.6)
###############################################################################
suppressPackageStartupMessages({
  library(optparse)
  library(future)
  library(Seurat)
  library(Signac)
  library(dplyr)
})

option_list <- list(
  make_option("--wild",                type="character"),
  make_option("--mutant",              type="character"),
  make_option("--npcs",                type="integer",  default=30),
  make_option("--dims",                type="integer",  default=30),
  make_option("--resolution_initial",  type="double",   default=0.3),
  make_option("--resolution_refine",   type="double",   default=0.6),
  make_option("--algorithm",           type="integer",  default=1),
  make_option("--threads",             type="integer",  default=20),
  make_option("--max_memory_gb",       type="integer",  default=80),
  make_option("--output",              type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

plan("multicore", workers = opt$threads)
options(future.globals.maxSize = opt$max_memory_gb * 1024^3)
future.seed <- TRUE
set.seed(1234)

# ---- Load processed objects --------------------------------------------------
cat(">>> Loading processed objects ...\n")
wild   <- readRDS(opt$wild)
mutant <- readRDS(opt$mutant)

# ---- SCT-based integration --------------------------------------------------
cat(">>> Integrating on SCT assay ...\n")
DefaultAssay(wild)   <- "SCT"
DefaultAssay(mutant) <- "SCT"
obj.list  <- list(wild, mutant)
features  <- SelectIntegrationFeatures(obj.list)
obj.list  <- PrepSCTIntegration(obj.list, anchor.features = features)
obj.anchors <- FindIntegrationAnchors(
  object.list          = obj.list,
  anchor.features      = features,
  normalization.method = "SCT"
)
combined <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT")
rm(wild, mutant, obj.list, obj.anchors); gc()

# ---- PCA + UMAP + initial clustering (resolution 0.3 → 6 clusters) ----------
cat(">>> PCA + UMAP + clustering (res", opt$resolution_initial, ") ...\n")
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = opt$npcs, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:opt$dims)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:opt$dims)
combined <- FindClusters(combined, resolution = opt$resolution_initial,
                         algorithm = opt$algorithm)

# Composite label: cluster_genotype
res_col <- paste0("integrated_snn_res.", opt$resolution_initial)
combined@meta.data$cl_genotype <- paste0(
  combined@meta.data[[res_col]], "_", combined@meta.data$kras_genotype
)

# ---- Refined clustering (resolution 0.6 → splits Pit sub-cluster) -----------
cat(">>> Refined clustering (res", opt$resolution_refine, ") ...\n")
combined <- FindClusters(combined, resolution = opt$resolution_refine)

# Map newly resolved cluster 9 (at res 0.6) → "6" label
combined@meta.data$cl_genotype[
  combined@meta.data$seurat_clusters == "9" &
  combined@meta.data$kras_genotype == "Wild_type"] <- "6_Wild_type"
combined@meta.data$cl_genotype[
  combined@meta.data$seurat_clusters == "9" &
  combined@meta.data$kras_genotype == "Mutant"] <- "6_Mutant"

# ---- Assign cell-type names -------------------------------------------------
# Clu 0 – Neck (Muc6, Cftr)     Clu 1 – Lgr5+ (Lgr5)
# Clu 2 – SPEM (Glipr1)         Clu 3 – Wnt7+ (Wnt7b)
# Clu 4 – Proliferating (Ki67)  Clu 5 – Pre-Pit (Tff1, Gkn1)
# Clu 6 – Pit (Muc5ac, Gkn2)
cluster_name_map <- c(
  "0" = "Neck", "1" = "Lgr5+", "2" = "SPEM", "3" = "Wnt7+",
  "4" = "Proliferating", "5" = "Pre-Pit", "6" = "Pit"
)

for (id in names(cluster_name_map)) {
  wt_label  <- paste0(id, "_Wild_type")
  mut_label <- paste0(id, "_Mutant")
  combined@meta.data$clusters[combined@meta.data$cl_genotype == wt_label]  <- cluster_name_map[id]
  combined@meta.data$clusters[combined@meta.data$cl_genotype == mut_label] <- cluster_name_map[id]
}

cat("   Final cluster counts:\n")
print(table(combined@meta.data$clusters))

# ---- Save --------------------------------------------------------------------
cat(">>> Saving to:", opt$output, "\n")
saveRDS(combined, file = opt$output)
cat(">>> Done.\n")