#!/usr/bin/env Rscript
# =============================================================================
#  Step 7: Motif Enrichment Analysis (Figures 3I, S3C)
# =============================================================================
#  Uses JASPAR2020 vertebrate motifs on DA peaks from Wnt7+ vs Lgr5+
#  comparisons in both Mutant and Wild-type conditions.
# =============================================================================

library(future)
library(Seurat)
library(Signac)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(chromVAR)
library(yaml)

cfg <- yaml.load_file("/pipeline/config.yaml")
plan("multicore", workers = cfg$threads)
options(future.globals.maxSize = cfg$max_memory_gb * 1024^3)
set.seed(1234)

out_dir <- cfg$output_dir
linkage_genome <- BSgenome.Mmusculus.UCSC.mm10

RZ.RZK <- readRDS(file.path(out_dir, "rds", "RZ_RZK_annotated.rds"))
da_mut <- read.csv(file.path(out_dir, "tables", "DA_peaks_Wnt7_vs_Lgr5_mutant.csv"),
                   row.names = 1)
da_wt  <- read.csv(file.path(out_dir, "tables", "DA_peaks_Wnt7_vs_Lgr5_wildtype.csv"),
                   row.names = 1)

# ---- Add JASPAR2020 motifs ---------------------------------------------------
cat("Adding JASPAR2020 vertebrate motifs ...\n")
DefaultAssay(RZ.RZK) <- "ATAC"
pfm <- getMatrixSet(x = JASPAR2020,
                    opts = list(collection = "CORE",
                                tax_group  = "vertebrates",
                                all_versions = FALSE))
RZ.RZK <- AddMotifs(object = RZ.RZK, genome = linkage_genome, pfm = pfm)

# ---- Helper: motif enrichment + bar plot -------------------------------------
run_motif_enrichment <- function(obj, da_peaks, title_str, out_file) {
    sig_peaks <- rownames(da_peaks[da_peaks$p_val_adj < cfg$da$padj_threshold_motif, ])
    enriched  <- FindMotifs(object = obj, features = sig_peaks)

    enriched$logP <- -log10(enriched$p.adjust)
    data     <- enriched[1:15, ]
    y.orders <- rev(data$motif.name)

    p <- ggplot(data, aes(x = logP, y = motif.name, fill = percent.observed)) +
        geom_bar(stat = "identity", width = 0.2) +
        scale_fill_gradientn(colours = c("royalblue", "rosybrown2", "red"),
                             limits = c(0, 60)) +
        scale_y_discrete(limits = y.orders) +
        xlim(0, 30) +
        ggtitle(title_str) + xlab("-log10(p.adjust)") + ylab("TF motif") +
        theme_linedraw() +
        theme(
            plot.title   = element_text(hjust = 0.5, size = 15, face = "bold"),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12)
        )

    ggsave(p, filename = out_file, width = 6, height = 6, dpi = 300)
    return(enriched)
}

fig_dir <- file.path(out_dir, "figures")

# ---- Figure 3I: Enriched motifs — Wnt7+ vs Lgr5+ in Mutant ------------------
cat("  Fig 3I: Motif enrichment (RZK / Mutant) ...\n")
motifs_mut <- run_motif_enrichment(
    RZ.RZK, da_mut,
    "Wnt7b+ vs Lgr5+\n(RZK)",
    file.path(fig_dir, "fig_3I_motifs_RZK.pdf")
)

# ---- Figure S3C: Enriched motifs — Wnt7+ vs Lgr5+ in Wild-type --------------
cat("  Fig S3C: Motif enrichment (RZ / Wild-type) ...\n")
motifs_wt <- run_motif_enrichment(
    RZ.RZK, da_wt,
    "Wnt7b+ vs Lgr5+\n(RZ)",
    file.path(fig_dir, "fig_S3C_motifs_RZ.pdf")
)

cat("Step 7 complete.\n")