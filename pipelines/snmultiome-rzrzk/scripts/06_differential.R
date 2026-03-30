#!/usr/bin/env Rscript
# =============================================================================
#  06_differential.R
#  Differential gene expression (MAST) between Mutant vs Wild_type
#  KRAS signaling pathway overlap analysis (MSigDB Hallmark)
#  Extracted from: 1_RZRZK_basic_analysis.R (DEG section)
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
  library(msigdbr)
  library(gprofiler2)
  library(MAST)
})

cfg <- read_yaml(config_file)

plan("multicore", workers = cfg$n_workers)
options(future.globals.maxSize = cfg$future_max_size_gb * 1024^3)
set.seed(cfg$seed)

cat("=== Step 6: Differential Expression ===\n")

# ── Load clustered object ────────────────────────────────────────────────────
obj <- readRDS("/output/04_cluster/combined_clustered.rds")

dir.create("/output/06_differential", showWarnings = FALSE, recursive = TRUE)

# ── DEG: Mutant vs Wild_type using MAST ──────────────────────────────────────
DefaultAssay(obj) <- "RNA"
cat("  Running FindMarkers (MAST): Mutant vs Wild_type...\n")
DEG <- FindMarkers(
  obj,
  group.by = "kras_genotype",
  ident.1  = "Mutant",
  test.use = cfg$deg$test_use
)
write.csv(DEG, "/output/06_differential/DEG_scMAST.csv")
cat("  Total DEGs found:", nrow(DEG), "\n")

# ── Filter significant up-regulated DEGs ──────────────────────────────────────
sig_up <- DEG[DEG$p_val_adj < cfg$deg$p_val_adj_threshold &
              DEG$avg_log2FC > cfg$deg$avg_log2FC_threshold, ]
cat("  Significant up-regulated genes (p_adj <", cfg$deg$p_val_adj_threshold,
    ", log2FC >", cfg$deg$avg_log2FC_threshold, "):", nrow(sig_up), "\n")

# ── KRAS signaling pathway overlap (Hallmark) ────────────────────────────────
cat("  Performing KRAS signaling pathway overlap...\n")
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

kras_signaling <- hallmark[hallmark$gs_name == "HALLMARK_KRAS_SIGNALING_UP", ]
kras_genes     <- kras_signaling$gene_symbol
cat("  Hallmark KRAS_SIGNALING_UP genes:", length(kras_genes), "\n")

# Convert mouse gene symbols to human orthologs via gprofiler2
converted <- gconvert(query = rownames(sig_up), organism = "hsapiens", target = "ENSG")
cat("  Converted DEGs to human orthologs:", length(converted$name), "\n")

# Find intersection
overlap <- intersect(converted$name, kras_genes)
cat("  KRAS signaling overlap genes:", length(overlap), "\n")

# Save overlap results
overlap_df <- data.frame(
  gene             = overlap,
  source           = "HALLMARK_KRAS_SIGNALING_UP",
  n_kras_genes     = length(kras_genes),
  n_deg_orthologs  = length(converted$name),
  n_overlap        = length(overlap)
)
write.csv(overlap_df, "/output/06_differential/DEG_KRAS_overlap.csv", row.names = FALSE)

cat("=== Step 6 complete ===\n")