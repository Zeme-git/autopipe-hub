#!/usr/bin/env Rscript
###############################################################################
#  Step 8: Differential Expression (MAST) + KRAS Hallmark overlap
###############################################################################
suppressPackageStartupMessages({
  library(optparse)
  library(future)
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(msigdbr)
  library(gprofiler2)
})

option_list <- list(
  make_option("--input",            type="character"),
  make_option("--test_use",         type="character", default="MAST"),
  make_option("--padj_threshold",   type="double",    default=0.01),
  make_option("--log2fc_threshold", type="double",    default=1),
  make_option("--output",           type="character"),
  make_option("--threads",          type="integer",   default=20),
  make_option("--max_memory_gb",    type="integer",   default=80)
)
opt <- parse_args(OptionParser(option_list=option_list))

plan("multicore", workers = opt$threads)
options(future.globals.maxSize = opt$max_memory_gb * 1024^3)
set.seed(1234)

# ---- Load object -------------------------------------------------------------
cat(">>> Loading clustered object ...\n")
obj <- readRDS(opt$input)
DefaultAssay(obj) <- "RNA"

# ---- DEG: Mutant vs Wild-type (MAST) ----------------------------------------
cat(">>> Running FindMarkers (MAST): Mutant vs Wild_type ...\n")
DEG <- FindMarkers(obj, group.by = "kras_genotype",
                   ident.1 = "Mutant", test.use = opt$test_use)

write.csv(DEG, file = opt$output)
cat("   Total DEGs tested:", nrow(DEG), "\n")

# ---- KRAS Hallmark overlap ---------------------------------------------------
cat(">>> HALLMARK_KRAS_SIGNALING_UP overlap ...\n")
sig_deg <- DEG[DEG$p_val_adj < opt$padj_threshold &
               DEG$avg_log2FC > opt$log2fc_threshold, ]
cat("   DEGs (padj <", opt$padj_threshold, ", LFC >", opt$log2fc_threshold, "):",
    nrow(sig_deg), "\n")

hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
kras_up <- hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_KRAS_SIGNALING_UP"]

# Convert mouse genes to human orthologs
ortho <- gconvert(query = rownames(sig_deg), organism = "hsapiens", target = "ENSG")
overlap <- intersect(ortho$name, kras_up)

cat("   KRAS_SIGNALING_UP genes  :", length(kras_up), "\n")
cat("   Converted DEG orthologs  :", length(ortho$name), "\n")
cat("   Overlap                  :", length(overlap), "\n")
cat("   Overlap genes            :", paste(overlap, collapse=", "), "\n")

# Save overlap summary
summary_df <- data.frame(
  metric = c("total_deg", "sig_deg", "kras_up_genes", "orthologs", "overlap"),
  value  = c(nrow(DEG), nrow(sig_deg), length(kras_up),
             length(ortho$name), length(overlap))
)
write.csv(summary_df,
          file = sub("\\.csv$", "_kras_summary.csv", opt$output),
          row.names = FALSE)

cat(">>> Done.\n")