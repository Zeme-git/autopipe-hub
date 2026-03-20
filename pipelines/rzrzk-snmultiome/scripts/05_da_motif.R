#!/usr/bin/env Rscript
###############################################################################
#  Step 10: Differential Accessibility + Motif Enrichment
#  DA peaks (Wnt7+ vs Lgr5+) in Mutant & WT + JASPAR2020 motif enrichment
###############################################################################
suppressPackageStartupMessages({
  library(optparse)
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

option_list <- list(
  make_option("--input",            type="character"),
  make_option("--species",          type="character", default="mouse"),
  make_option("--min_pct",          type="double",    default=0.05),
  make_option("--logfc_threshold",  type="double",    default=0.1),
  make_option("--test_use",         type="character", default="LR"),
  make_option("--latent_vars",      type="character", default="atac_peak_region_fragments"),
  make_option("--padj_motif",       type="double",    default=0.005),
  make_option("--outdir",           type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Single thread for DA analysis (memory-safe, same as original script)
plan("multicore", workers = 1)
set.seed(1234)

linkage_genome <- BSgenome.Mmusculus.UCSC.mm10

# ---- Load object -------------------------------------------------------------
cat(">>> Loading clustered object ...\n")
RZ.RZK <- readRDS(opt$input)
DefaultAssay(RZ.RZK) <- "ATAC"

dir.create(file.path(opt$outdir, "figures"), showWarnings=FALSE, recursive=TRUE)

# ---- DA peaks: Wnt7+ (cluster 3) vs Lgr5+ (cluster 1) ----------------------

# -- Mutant
cat(">>> DA peaks: Wnt7+ vs Lgr5+ in Mutant ...\n")
da_m <- FindMarkers(
  RZ.RZK, ident.1="3_Mutant", ident.2="1_Mutant",
  group.by="cl_genotype", min.pct=opt$min_pct,
  logfc.threshold=opt$logfc_threshold,
  test.use=opt$test_use, latent.vars=opt$latent_vars
)
cg_m <- ClosestFeature(RZ.RZK, regions=rownames(da_m))
rownames(cg_m) <- cg_m$query_region
da_m$gene         <- cg_m$gene_name
da_m$gene_biotype <- cg_m$gene_biotype
da_m$distance     <- cg_m$distance
write.csv(da_m, file.path(opt$outdir, "da_peaks_mutant.csv"))

# -- Wild-type
cat(">>> DA peaks: Wnt7+ vs Lgr5+ in Wild-type ...\n")
da_w <- FindMarkers(
  RZ.RZK, ident.1="3_Wild_type", ident.2="1_Wild_type",
  group.by="cl_genotype", min.pct=opt$min_pct,
  logfc.threshold=opt$logfc_threshold,
  test.use=opt$test_use, latent.vars=opt$latent_vars
)
cg_w <- ClosestFeature(RZ.RZK, regions=rownames(da_w))
rownames(cg_w) <- cg_w$query_region
da_w$gene         <- cg_w$gene_name
da_w$gene_biotype <- cg_w$gene_biotype
da_w$distance     <- cg_w$distance
write.csv(da_w, file.path(opt$outdir, "da_peaks_wildtype.csv"))

# ---- Motif enrichment (JASPAR2020 vertebrate core) ---------------------------
cat(">>> Adding motifs (JASPAR2020) ...\n")
pfm <- getMatrixSet(x=JASPAR2020,
                    opts=list(collection="CORE", tax_group="vertebrates",
                              all_versions=FALSE))
RZ.RZK <- AddMotifs(object=RZ.RZK, genome=linkage_genome, pfm=pfm)

# -- Mutant enriched motifs
cat(">>> Motif enrichment: Mutant (Wnt7+ vs Lgr5+) ...\n")
motifs_m <- FindMotifs(object=RZ.RZK,
  features=rownames(da_m[da_m$p_val_adj < opt$padj_motif, ]))
motifs_m$logP <- -log10(motifs_m$p.adjust)

p_3I <- ggplot(motifs_m[1:15,], aes(x=logP, y=motif.name, fill=percent.observed)) +
  geom_bar(stat="identity", width=0.2) +
  scale_fill_gradientn(colours=c("royalblue","rosybrown2","red"), limits=c(0,60)) +
  scale_y_discrete(limits=rev(motifs_m$motif.name[1:15])) + xlim(0,30) +
  ggtitle("Wnt7b+ vs Lgr5+\n(RZK)") + xlab("-logP") + ylab("TF motif") +
  theme_linedraw() +
  theme(plot.title=element_text(hjust=0.5, size=15, face="bold"),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12))
ggsave(p_3I, filename=file.path(opt$outdir,"figures","fig_3I_motifs_Wnt_vs_Lgr5_RZK.pdf"),
       width=6, height=6, dpi=300)

# -- Wild-type enriched motifs
cat(">>> Motif enrichment: Wild-type (Wnt7+ vs Lgr5+) ...\n")
motifs_w <- FindMotifs(object=RZ.RZK,
  features=rownames(da_w[da_w$p_val_adj < opt$padj_motif, ]))
motifs_w$logP <- -log10(motifs_w$p.adjust)

p_S3C <- ggplot(motifs_w[1:15,], aes(x=logP, y=motif.name, fill=percent.observed)) +
  geom_bar(stat="identity", width=0.2) +
  scale_fill_gradientn(colours=c("royalblue","rosybrown2","red"), limits=c(0,60)) +
  scale_y_discrete(limits=rev(motifs_w$motif.name[1:15])) + xlim(0,30) +
  ggtitle("Wnt7b+ vs Lgr5+\n(RZ)") + xlab("-logP") + ylab("TF motif") +
  theme_linedraw() +
  theme(plot.title=element_text(hjust=0.5, size=15, face="bold"),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12))
ggsave(p_S3C, filename=file.path(opt$outdir,"figures","fig_S3C_motifs_Wnt_vs_Lgr5_RZ.pdf"),
       width=6, height=6, dpi=300)

# Save motif tables
write.csv(motifs_m, file.path(opt$outdir, "motifs_enriched_mutant.csv"), row.names=FALSE)
write.csv(motifs_w, file.path(opt$outdir, "motifs_enriched_wildtype.csv"), row.names=FALSE)

cat(">>> Done.\n")