#!/usr/bin/env Rscript
# =============================================================================
#  07_motif_enrichment.R
#  Differential accessibility (DA) peaks: Wnt7+ vs Lgr5+ in Mutant & Wild_type
#  Motif enrichment analysis using JASPAR2020
#  Wnt7b locus coverage analysis and statistics
#  Extracted from: 2_RZRZK_figure.R (fig 3H, 3I, S3C)
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
  library(dplyr)
  library(ggplot2)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(motifmatchr)
  library(JASPAR2020)
  library(TFBSTools)
  library(chromVAR)
})

cfg <- read_yaml(config_file)

# Use single core for DA analysis (as in original code)
plan("multicore", workers = 1)
options(future.globals.maxSize = cfg$future_max_size_gb * 1024^3)
set.seed(cfg$seed)

cat("=== Step 7: Motif Enrichment ===\n")

linkage_genome <- BSgenome.Mmusculus.UCSC.mm10

# ── Load refined object ──────────────────────────────────────────────────────
RZ.RZK <- readRDS("/output/05_visualization/combined_refined.rds")

dir.create("/output/07_motif", showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# DA peaks: Wnt7+ (cluster 3) vs Lgr5+ (cluster 1)
# Uses LR test with latent variable = atac_peak_region_fragments
# ══════════════════════════════════════════════════════════════════════════════
DefaultAssay(RZ.RZK) <- "ATAC"

cat("  Finding DA peaks: Wnt7+ vs Lgr5+ (Mutant)...\n")
da_peaks.31.m <- FindMarkers(
  RZ.RZK,
  ident.1    = "3_Mutant",
  ident.2    = "1_Mutant",
  group.by   = "cl_genotype",
  min.pct    = cfg$da_params$min_pct,
  logfc.threshold = cfg$da_params$logfc_threshold,
  test.use   = cfg$da_params$test_use,
  latent.vars = cfg$da_params$latent_vars
)

cat("  Finding DA peaks: Wnt7+ vs Lgr5+ (Wild_type)...\n")
da_peaks.31.w <- FindMarkers(
  RZ.RZK,
  ident.1    = "3_Wild_type",
  ident.2    = "1_Wild_type",
  group.by   = "cl_genotype",
  min.pct    = cfg$da_params$min_pct,
  logfc.threshold = cfg$da_params$logfc_threshold,
  test.use   = cfg$da_params$test_use,
  latent.vars = cfg$da_params$latent_vars
)

# ── Annotate DA peaks with closest genes ──────────────────────────────────────
annotate_da_peaks <- function(obj, da_peaks) {
  closest_genes <- ClosestFeature(obj, regions = rownames(da_peaks))
  rownames(closest_genes) <- closest_genes$query_region
  da_peaks$gene           <- closest_genes$gene_name
  da_peaks$gene_biotype   <- closest_genes$gene_biotype
  da_peaks$closest_region <- closest_genes$closest_region
  da_peaks$distance       <- closest_genes$distance
  return(da_peaks)
}

da_peaks.31.m <- annotate_da_peaks(RZ.RZK, da_peaks.31.m)
da_peaks.31.w <- annotate_da_peaks(RZ.RZK, da_peaks.31.w)

# ══════════════════════════════════════════════════════════════════════════════
# Motif Enrichment using JASPAR2020 vertebrate core collection
# ══════════════════════════════════════════════════════════════════════════════
cat("  Adding motifs (JASPAR2020 vertebrate core)...\n")
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)
RZ.RZK <- AddMotifs(object = RZ.RZK, genome = linkage_genome, pfm = pfm)

# ── Motif enrichment: Mutant (Wnt7+ vs Lgr5+) ───────────────────────────────
cat("  Motif enrichment in Mutant DA peaks...\n")
sig_peaks_m <- rownames(da_peaks.31.m[da_peaks.31.m$p_val_adj < cfg$da_params$motif_p_adj_threshold, ])
motifs_mutant <- FindMotifs(object = RZ.RZK, features = sig_peaks_m)
motifs_mutant$logP <- -log10(motifs_mutant$p.adjust)
write.csv(motifs_mutant, "/output/07_motif/motif_enrichment_mutant.csv", row.names = FALSE)

# ── Motif enrichment: Wild_type (Wnt7+ vs Lgr5+) ─────────────────────────────
cat("  Motif enrichment in Wild_type DA peaks...\n")
sig_peaks_w <- rownames(da_peaks.31.w[da_peaks.31.w$p_val_adj < cfg$da_params$motif_p_adj_threshold, ])
motifs_wt <- FindMotifs(object = RZ.RZK, features = sig_peaks_w)
motifs_wt$logP <- -log10(motifs_wt$p.adjust)
write.csv(motifs_wt, "/output/07_motif/motif_enrichment_wildtype.csv", row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Fig 3I: Enriched motif bar plots
# ══════════════════════════════════════════════════════════════════════════════

plot_motif_barplot <- function(motif_df, title_str, out_file) {
  sig_motifs <- motif_df[motif_df$p.adjust < 0.05, ]
  data <- head(sig_motifs, 15)
  y.orders <- rev(data$motif.name)
  
  p <- ggplot(data, aes(x = logP, y = motif.name, fill = percent.observed)) +
    geom_bar(stat = "identity", width = 0.2) +
    scale_fill_gradientn(colours = c("royalblue", "rosybrown2", "red"), limits = c(0, 60)) +
    scale_y_discrete(limits = y.orders) +
    xlim(0, 30) +
    ggtitle(title_str) +
    xlab("-logP") + ylab("TF motif") +
    theme_linedraw() +
    theme(
      plot.title   = element_text(hjust = 0.5, color = "black", size = 15, face = "bold"),
      axis.title.x = element_text(color = "black", size = 12),
      axis.title.y = element_text(color = "black", size = 12)
    )
  
  ggsave(file = out_file, width = 6, height = 6, dpi = 300, plot = p)
}

cat("  Generating motif bar plots...\n")
plot_motif_barplot(motifs_mutant, "Wnt7b+ vs Lgr5+\n(RZK)",
                   "/output/07_motif/fig_3I_motif_barplot_mutant.pdf")
plot_motif_barplot(motifs_wt, "Wnt7b+ vs Lgr5+\n(RZ)",
                   "/output/07_motif/fig_3I_motif_barplot_wildtype.pdf")

# ══════════════════════════════════════════════════════════════════════════════
# Wnt7b locus: Wilcoxon expression tests + ATAC region enrichment (statistics)
# ══════════════════════════════════════════════════════════════════════════════
cat("  Computing Wnt7b expression and accessibility statistics...\n")

# Wilcoxon tests for Wnt7b expression per cluster
# With v3 assay, [["RNA"]]@data works as expected
wilcox_results <- data.frame()
for (cl in 0:6) {
  wt_cells <- colnames(subset(RZ.RZK, subset = cl_genotype == paste0(cl, "_Wild_type")))
  mt_cells <- colnames(subset(RZ.RZK, subset = cl_genotype == paste0(cl, "_Mutant")))
  
  wt_expr <- RZ.RZK[["RNA"]]@data["Wnt7b", wt_cells]
  mt_expr <- RZ.RZK[["RNA"]]@data["Wnt7b", mt_cells]
  
  res <- wilcox.test(mt_expr, wt_expr)
  wilcox_results <- rbind(wilcox_results, data.frame(
    cluster = cl, p_value = res$p.value,
    mean_mutant = mean(mt_expr), mean_wildtype = mean(wt_expr)
  ))
}
write.csv(wilcox_results, "/output/07_motif/wnt7b_wilcoxon_tests.csv", row.names = FALSE)

# ATAC enrichment tests for Wnt7b regulatory regions
enrichment_results <- data.frame()
for (region_str in cfg$wnt7b_regions) {
  counts_region <- CountsInRegion(
    object  = RZ.RZK,
    assay   = "macs2",
    regions = StringToGRanges(region_str)
  )
  for (cl in 0:6) {
    mt_cells <- colnames(subset(RZ.RZK, subset = cl_genotype == paste0(cl, "_Mutant")))
    wt_cells <- colnames(subset(RZ.RZK, subset = cl_genotype == paste0(cl, "_Wild_type")))
    res <- t.test(counts_region[mt_cells], counts_region[wt_cells])
    enrichment_results <- rbind(enrichment_results, data.frame(
      region = region_str, cluster = cl,
      p_value = res$p.value,
      mean_mutant = res$estimate[1], mean_wildtype = res$estimate[2]
    ))
  }
}
write.csv(enrichment_results, "/output/07_motif/wnt7b_atac_enrichment.csv", row.names = FALSE)

# ── Coverage plot for Wnt7b locus (Fig 3H) ───────────────────────────────────
cat("  Generating Wnt7b coverage plot...\n")
DefaultAssay(RZ.RZK) <- "macs2"
RZ.RZK <- RegionStats(RZ.RZK, genome = linkage_genome)
RZ.RZK <- LinkPeaks(
  object           = RZ.RZK,
  peak.assay       = "macs2",
  expression.assay = "RNA",
  genes.use        = c("Wnt7b")
)

lv.list <- c()
for (i in c(1, 2, 0, 4, 5, 6, 3)) {
  lv.list <- c(lv.list, sprintf("%s_Wild_type", i))
  lv.list <- c(lv.list, sprintf("%s_Mutant", i))
}
RZ.RZK$cl_genotype <- factor(RZ.RZK$cl_genotype, levels = lv.list)

pdf("/output/07_motif/fig_3H_coverage_Wnt7b.pdf", width = 8, height = 8)
p_cov <- CoveragePlot(
  object             = RZ.RZK,
  region             = "Wnt7b",
  features           = "Wnt7b",
  expression.assay   = "RNA",
  group.by           = "cl_genotype",
  extend.upstream    = 0,
  extend.downstream  = 0
)
print(p_cov)
dev.off()

cat("=== Step 7 complete ===\n")