#!/usr/bin/env Rscript
# Step 4: SCT Integration, Clustering & Cell-Type Annotation

library(future)
library(Seurat)
library(Signac)
library(yaml)

cfg <- yaml.load_file("/pipeline/config.yaml")
plan("multicore", workers = cfg$threads)
options(future.globals.maxSize = cfg$max_memory_gb * 1024^3)
set.seed(1234)

out_dir <- cfg$output_dir
wild   <- readRDS(file.path(out_dir, "rds", "wild_normalized.rds"))
mutant <- readRDS(file.path(out_dir, "rds", "mutant_normalized.rds"))

cat("Finding integration anchors (SCT) ...\n")
DefaultAssay(wild) <- "SCT"; DefaultAssay(mutant) <- "SCT"
obj.list <- list(wild, mutant)
features <- SelectIntegrationFeatures(obj.list)
obj.list <- PrepSCTIntegration(obj.list, anchor.features = features)
anchors  <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, normalization.method = "SCT")
combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
rm(wild, mutant, obj.list, anchors); gc()

cat("Running PCA, UMAP, clustering ...\n")
npcs <- cfg$integration$npcs
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = npcs, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:npcs)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:npcs)
combined <- FindClusters(combined, resolution = cfg$integration$cluster_resolution_coarse, algorithm = cfg$integration$algorithm)
combined@meta.data$cl_genotype <- paste0(combined@meta.data$integrated_snn_res.0.3, "_", combined@meta.data$kras_genotype)

cat("Re-clustering at resolution 0.6 to resolve Pit subtypes ...\n")
combined <- FindClusters(combined, resolution = cfg$integration$cluster_resolution_fine)
combined@meta.data$cl_genotype[combined@meta.data$seurat_clusters == "9" & combined@meta.data$kras_genotype == "Wild_type"] <- "6_Wild_type"
combined@meta.data$cl_genotype[combined@meta.data$seurat_clusters == "9" & combined@meta.data$kras_genotype == "Mutant"] <- "6_Mutant"

cluster_names <- cfg$cluster_names
combined@meta.data$clusters <- cluster_names[sub("_.*", "", combined@meta.data$cl_genotype)]
cat("Cell-type distribution:\n")
print(table(combined@meta.data$clusters))

saveRDS(combined, file.path(out_dir, "rds", "RZ_RZK_annotated.rds"))
cat("Step 4 complete.\n")