#!/usr/bin/env Rscript
# =============================================================================
#  05_visualization.R
#  Generate publication figures:
#    Fig 3A: UMAP colored by cl_genotype
#    Fig 3B: DotPlot of marker genes across clusters
#    Fig 3C: DittoBarPlot of cluster proportions per genotype
#    Fig 3D: Wnt7b FeaturePlot split by genotype
#    Fig S3A: FeaturePlots for all marker genes
#    Fig S3B: Wnt family VlnPlot
#  Also performs refined clustering (res=0.6) to split Pit cells (cluster 9→6)
#  Extracted from: 2_RZRZK_figure.R
#  NOTE: Uses Seurat v3 assay version for v4-style object compatibility
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]

# ── Force Seurat v3/v4 assay structure in Seurat v5 ──────────────────────────
# This ensures [["RNA"]]@data, FeaturePlot()$data column names, etc. work
# identically to the original v4-based analysis scripts.
options(Seurat.object.assay.version = "v3")

suppressPackageStartupMessages({
  library(yaml)
  library(future)
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(dittoSeq)
  library(paletteer)
})

cfg <- read_yaml(config_file)

plan("multicore", workers = cfg$n_workers)
options(future.globals.maxSize = cfg$future_max_size_gb * 1024^3)
set.seed(cfg$seed)

cat("=== Step 5: Visualization ===\n")

# ── Load clustered object ────────────────────────────────────────────────────
ready_to_visualize <- readRDS("/output/04_cluster/combined_clustered.rds")

dir.create("/output/05_visualization", showWarnings = FALSE, recursive = TRUE)

# ── Define color palettes ────────────────────────────────────────────────────
my_col <- c(
  "0_Wild_type" = "#8936EF", "0_Mutant" = "#8936EF",
  "1_Wild_type" = "#F2CA19", "1_Mutant" = "#F2CA19",
  "2_Wild_type" = "#FF00BD", "2_Mutant" = "#FF00BD",
  "3_Wild_type" = "#E11845", "3_Mutant" = "#E11845",
  "4_Wild_type" = "#0057E9", "4_Mutant" = "#0057E9",
  "5_Wild_type" = "#87E911", "5_Mutant" = "#87E911",
  "6_Wild_type" = "#018300", "6_Mutant" = "#018300"
)
cluster_col <- c(
  "Neck" = "#8936EF", "Lgr5+" = "#F2CA19", "SPEM" = "#FF00BD",
  "Wnt7+" = "#E11845", "Proliferating" = "#0057E9",
  "Pre-Pit" = "#87E911", "Pit" = "#018300"
)

# ── Refined clustering (resolution=0.6) to split Pit cells ───────────────────
DefaultAssay(ready_to_visualize) <- "integrated"
RZ.RZK <- FindClusters(ready_to_visualize, resolution = cfg$clustering$resolution_refined)

# Map cluster 9 (from res=0.6) into cluster 6 for cl_genotype
RZ.RZK@meta.data$cl_genotype[RZ.RZK@meta.data$seurat_clusters == "9" &
                               RZ.RZK@meta.data$kras_genotype == "Wild_type"] <- "6_Wild_type"
RZ.RZK@meta.data$cl_genotype[RZ.RZK@meta.data$seurat_clusters == "9" &
                               RZ.RZK@meta.data$kras_genotype == "Mutant"] <- "6_Mutant"

# Assign cell type labels
label_map <- list(
  "0" = "Neck", "1" = "Lgr5+", "2" = "SPEM",
  "3" = "Wnt7+", "4" = "Proliferating", "5" = "Pre-Pit", "6" = "Pit"
)
for (cl_geno in names(label_map)) {
  for (gt in c("Wild_type", "Mutant")) {
    key <- paste0(cl_geno, "_", gt)
    RZ.RZK@meta.data$clusters[RZ.RZK@meta.data$cl_genotype == key] <- label_map[[cl_geno]]
  }
}

# ── Custom FeaturePlot functions (v5-compatible) ──────────────────────────────
# Uses Embeddings() to get UMAP coordinates directly, avoiding $data issues
plot_featureplot <- function(obj, gene) {
  DefaultAssay(obj) <- "RNA"
  # Get UMAP embeddings and expression data directly
  umap_coords <- Embeddings(obj, reduction = "umap")
  expr_data <- FetchData(obj, vars = gene, layer = "data")
  plot_df <- data.frame(
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    gene_expr = expr_data[, 1]
  )
  # Apply cutoffs as in original
  plot_df$gene_expr[plot_df$gene_expr < 0.3] <- 0.3
  plot_df$gene_expr[plot_df$gene_expr > 2] <- 2

  ggplot(data = plot_df, aes(x = UMAP_1, y = UMAP_2, fill = gene_expr)) +
    geom_point(shape = 21, stroke = 0.3, color = "black", alpha = 1, size = 2) +
    scale_fill_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd")) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(fill = gene)
}

plot_featureplot_v3 <- function(obj, gene, min_val, max_val) {
  DefaultAssay(obj) <- "RNA"
  ht_custom_col <- rev(paletteer_c("grDevices::Spectral", 30))
  # Get UMAP embeddings and expression data directly
  umap_coords <- Embeddings(obj, reduction = "umap")
  expr_data <- FetchData(obj, vars = gene, layer = "data")
  plot_df <- data.frame(
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    gene_expr = expr_data[, 1]
  )
  # Apply cutoffs
  plot_df$gene_expr[plot_df$gene_expr < min_val] <- min_val
  plot_df$gene_expr[plot_df$gene_expr > max_val] <- max_val
  # Sort so highest expression is plotted on top (order=TRUE equivalent)
  plot_df <- plot_df[order(plot_df$gene_expr), ]

  ggplot(data = plot_df, aes(x = UMAP_1, y = UMAP_2, fill = gene_expr)) +
    geom_point(shape = 21, stroke = NA, color = "black", alpha = 1, size = 2) +
    scale_fill_gradientn(colours = ht_custom_col, name = "Expression\nLevel") +
    ggtitle(gene) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          plot.subtitle = element_text(hjust = 0.5))
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3A: UMAP colored by cluster + genotype
# ══════════════════════════════════════════════════════════════════════════════
cat("  Generating Fig 3A...\n")
pdf("/output/05_visualization/fig_3A_UMAP.pdf", width = 8, height = 9)
p1 <- DimPlot(RZ.RZK, group.by = "cl_genotype", cols = my_col, pt.size = 1.3) + NoLegend()
print(p1)
dev.off()

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3B: DotPlot of marker genes
# ══════════════════════════════════════════════════════════════════════════════
cat("  Generating Fig 3B...\n")
DefaultAssay(RZ.RZK) <- "RNA"
mylevel <- c("Lgr5+", "SPEM", "Neck", "Proliferating", "Pre-Pit", "Pit", "Wnt7+")
Idents(RZ.RZK) <- RZ.RZK@meta.data$clusters
Idents(RZ.RZK) <- factor(Idents(RZ.RZK), levels = mylevel)

pdf("/output/05_visualization/fig_3B_dotplot.pdf", width = 6, height = 7)
p2 <- DotPlot(RZ.RZK, features = cfg$marker_genes, dot.scale = 10) +
  scale_colour_gradient2(low = "dodgerblue3", mid = "ghostwhite", high = "firebrick",
                         limits = c(-2.5, 2.5)) +
  xlab("Markers") + ylab("Clusters") + coord_flip() +
  theme(axis.text = element_text(size = 15), legend.position = "top")
print(p2)
dev.off()

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3C: DittoBarPlot (cluster proportions per genotype)
# ══════════════════════════════════════════════════════════════════════════════
cat("  Generating Fig 3C...\n")
pdf("/output/05_visualization/fig_3C_dittobarplot.pdf", width = 6, height = 6)
p3 <- dittoBarPlot(RZ.RZK, Idents(RZ.RZK), group.by = "kras_genotype",
                   x.reorder = c(2, 1),
                   var.labels.reorder = c(1, 6, 2, 5, 4, 3, 7),
                   color.panel = cluster_col)
print(p3)
dev.off()

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3D: Wnt7b FeaturePlot split by genotype
# ══════════════════════════════════════════════════════════════════════════════
cat("  Generating Fig 3D...\n")
pdf("/output/05_visualization/fig_3D_Wnt7b_feature.pdf", width = 12, height = 6)
p1_wt  <- plot_featureplot(subset(RZ.RZK, subset = kras_genotype == "Wild_type"), "Wnt7b")
p2_mut <- plot_featureplot(subset(RZ.RZK, subset = kras_genotype == "Mutant"), "Wnt7b")
print(p1_wt + p2_mut)
dev.off()

# ══════════════════════════════════════════════════════════════════════════════
# Fig S3A: FeaturePlots for multiple marker genes
# ══════════════════════════════════════════════════════════════════════════════
cat("  Generating Fig S3A...\n")
pdf("/output/05_visualization/fig_S3A_feature_plots.pdf", width = 10, height = 10)
for (g in cfg$feature_plot_genes) {
  p <- plot_featureplot_v3(RZ.RZK, g$gene, g$min, g$max)
  print(p)
}
dev.off()

# ══════════════════════════════════════════════════════════════════════════════
# Fig S3B: Wnt family VlnPlot by genotype
# ══════════════════════════════════════════════════════════════════════════════
cat("  Generating Fig S3B...\n")
DefaultAssay(RZ.RZK) <- "RNA"
pdf("/output/05_visualization/fig_S3B_wnt_vlnplot.pdf", width = 18, height = 12)
p_wnt <- VlnPlot(RZ.RZK, features = cfg$wnt_family_genes, group.by = "kras_genotype")
print(p_wnt)
dev.off()

# ── Save refined object for downstream motif analysis ─────────────────────────
saveRDS(RZ.RZK, "/output/05_visualization/combined_refined.rds")

cat("=== Step 5 complete ===\n")