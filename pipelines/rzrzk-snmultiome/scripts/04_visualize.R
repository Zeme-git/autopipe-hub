#!/usr/bin/env Rscript
###############################################################################
#  Step 9: Visualisation — Figures 3A-3E, 3H, S3A, S3B
###############################################################################
suppressPackageStartupMessages({
  library(optparse)
  library(future)
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(dittoSeq)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(paletteer)
})

option_list <- list(
  make_option("--input",          type="character"),
  make_option("--species",        type="character", default="mouse"),
  make_option("--outdir",         type="character"),
  make_option("--threads",        type="integer",   default=20),
  make_option("--max_memory_gb",  type="integer",   default=80)
)
opt <- parse_args(OptionParser(option_list=option_list))

plan("multicore", workers = opt$threads)
options(future.globals.maxSize = opt$max_memory_gb * 1024^3)
set.seed(1234)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
linkage_genome <- BSgenome.Mmusculus.UCSC.mm10

# ---- Colour palettes ---------------------------------------------------------
cl_geno_cols <- c(
  "0_Wild_type"="#8936EF","0_Mutant"="#8936EF",
  "1_Wild_type"="#F2CA19","1_Mutant"="#F2CA19",
  "2_Wild_type"="#FF00BD","2_Mutant"="#FF00BD",
  "3_Wild_type"="#E11845","3_Mutant"="#E11845",
  "4_Wild_type"="#0057E9","4_Mutant"="#0057E9",
  "5_Wild_type"="#87E911","5_Mutant"="#87E911",
  "6_Wild_type"="#018300","6_Mutant"="#018300"
)
cluster_cols <- c(
  "Neck"="#8936EF","Lgr5+"="#F2CA19","SPEM"="#FF00BD","Wnt7+"="#E11845",
  "Proliferating"="#0057E9","Pre-Pit"="#87E911","Pit"="#018300"
)
mylevel <- c("Lgr5+","SPEM","Neck","Proliferating","Pre-Pit","Pit","Wnt7+")

# ---- Load object -------------------------------------------------------------
cat(">>> Loading clustered object ...\n")
RZ.RZK <- readRDS(opt$input)

# ---- Helper functions --------------------------------------------------------
plot_featureplot <- function(obj, gene) {
  DefaultAssay(obj) <- "RNA"
  p <- FeaturePlot(obj, features=gene, pt.size=1, min.cutoff=0.3, max.cutoff=2)
  ggplot(data=p$data, aes_string(x="UMAP_1", y="UMAP_2", fill=gene)) +
    geom_point(shape=21, stroke=0.3, color="black", alpha=1, size=2) +
    scale_fill_gradientn(colours=brewer.pal(9,"YlOrRd")) +
    theme_void() +
    theme(plot.title=element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5))
}

plot_featureplot_v3 <- function(obj, gene, min_cut, max_cut) {
  DefaultAssay(obj) <- "RNA"
  ht_custom_col <- rev(paletteer_c("grDevices::Spectral", 30))
  p <- FeaturePlot(obj, features=gene, pt.size=1,
                   min.cutoff=min_cut, max.cutoff=max_cut, order=TRUE)
  ggplot(data=p$data, aes_string(x="UMAP_1", y="UMAP_2", fill=gene)) +
    geom_point(shape=21, stroke=NA, color="black", alpha=1, size=2) +
    scale_fill_gradientn(colours=ht_custom_col, name="Expression\nLevel") +
    ggtitle(gene) +
    theme_void() +
    theme(plot.title=element_text(hjust=0.5, size=16))
}

plot_linkage_cl_genotype <- function(obj, assay, target_gene) {
  DefaultAssay(obj) <- assay
  obj <- RegionStats(obj, genome=linkage_genome)
  obj <- LinkPeaks(object=obj, peak.assay=assay,
                   expression.assay="RNA", genes.use=target_gene)
  lv.list <- c()
  for (i in c(1,2,0,4,5,6,3)) {
    lv.list <- c(lv.list, sprintf("%s_Wild_type",i), sprintf("%s_Mutant",i))
  }
  obj$cl_genotype <- factor(obj$cl_genotype, levels=lv.list)
  CoveragePlot(object=obj, region=target_gene, features=target_gene,
               expression.assay="RNA", group.by="cl_genotype",
               extend.upstream=0, extend.downstream=0)
}

# ==== FIGURE 3A — UMAP by cluster × genotype =================================
cat(">>> Fig 3A\n")
p3A <- DimPlot(RZ.RZK, group.by="cl_genotype", cols=cl_geno_cols, pt.size=1.3) + NoLegend()
ggsave(p3A, filename=file.path(opt$outdir,"fig_3A_UMAP.pdf"), width=8, height=9)

# ==== FIGURE 3B — DotPlot of markers =========================================
cat(">>> Fig 3B\n")
DefaultAssay(RZ.RZK) <- "RNA"
Idents(RZ.RZK) <- factor(RZ.RZK@meta.data$clusters, levels=mylevel)
marker_genes <- c("Lgr5","Glipr1","Cd44","Cftr","Muc6",
                  "Mki67","Foxm1","Hmgb2","Top2a","Smc2",
                  "Gkn1","Tff1","Gkn2","Muc5ac","Wnt7b")
p3B <- DotPlot(RZ.RZK, features=marker_genes, dot.scale=10) +
  scale_colour_gradient2(low="dodgerblue3", mid="ghostwhite",
                         high="firebrick", limits=c(-2.5,2.5)) +
  xlab("Markers") + ylab("Clusters") + coord_flip() +
  theme(axis.text=element_text(size=15), legend.position="top")
ggsave(p3B, filename=file.path(opt$outdir,"fig_3B_dotplot.pdf"), width=6, height=7)

# ==== FIGURE 3C — dittoBarPlot ================================================
cat(">>> Fig 3C\n")
p3C <- dittoBarPlot(RZ.RZK, Idents(RZ.RZK), group.by="kras_genotype",
                    x.reorder=c(2,1), var.labels.reorder=c(1,6,2,5,4,3,7),
                    color.panel=cluster_cols)
ggsave(p3C, filename=file.path(opt$outdir,"fig_3C_dittobarplot.pdf"), width=6, height=6)

# ==== FIGURE 3D — Wnt7b FeaturePlot split by genotype ========================
cat(">>> Fig 3D\n")
p3D_wt  <- plot_featureplot(subset(RZ.RZK, subset=kras_genotype=="Wild_type"), "Wnt7b")
p3D_mut <- plot_featureplot(subset(RZ.RZK, subset=kras_genotype=="Mutant"),    "Wnt7b")
p3D <- p3D_wt + p3D_mut
ggsave(p3D, filename=file.path(opt$outdir,"fig_3D_Wnt7b_feature.pdf"), width=12, height=6, dpi=300)

# ==== FIGURE 3E — Wnt7b-expressing cell % per cluster × genotype =============
cat(">>> Fig 3E\n")
cl_labels <- c("0_Wild_type","0_Mutant","1_Wild_type","1_Mutant",
               "2_Wild_type","2_Mutant","3_Wild_type","3_Mutant",
               "4_Wild_type","4_Mutant","5_Wild_type","5_Mutant",
               "6_Wild_type","6_Mutant")
pct_vec <- c()
for (cl in cl_labels) {
  sub <- subset(RZ.RZK, subset=cl_genotype==cl)
  counts <- sub[["RNA"]]@data["Wnt7b",]
  pct_vec <- c(pct_vec, sum(counts!=0)/length(counts)*100)
}
df_3E <- data.frame(Clusters=cl_labels, Percentage=pct_vec)
x_order <- c("1_Wild_type","1_Mutant","2_Wild_type","2_Mutant",
             "0_Wild_type","0_Mutant","4_Wild_type","4_Mutant",
             "5_Wild_type","5_Mutant","6_Wild_type","6_Mutant",
             "3_Wild_type","3_Mutant")
p3E <- ggplot(df_3E, aes(Clusters, Percentage, fill=Clusters)) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=x_order) + ylim(0,100) +
  theme_classic() + scale_fill_manual(values=cl_geno_cols) +
  ggtitle("Wnt7b-expressing cells (%)") +
  theme(plot.title=element_text(hjust=0.5)) + NoLegend()
ggsave(p3E, filename=file.path(opt$outdir,"fig_3E_Wnt7b_pct.pdf"), width=10, height=6, dpi=300)

# Fisher tests + DotPlot for Wnt7b split by genotype
cl_pairs <- list(c("0_Wild_type","0_Mutant"),c("1_Wild_type","1_Mutant"),
                 c("2_Wild_type","2_Mutant"),c("3_Wild_type","3_Mutant"),
                 c("4_Wild_type","4_Mutant"),c("5_Wild_type","5_Mutant"),
                 c("6_Wild_type","6_Mutant"))
fisher_results <- data.frame(cluster=character(), p_value=numeric(), odds_ratio=numeric(),
                             stringsAsFactors=FALSE)
for (pair in cl_pairs) {
  s1 <- subset(RZ.RZK, subset=cl_genotype==pair[1])
  s2 <- subset(RZ.RZK, subset=cl_genotype==pair[2])
  c1 <- s1[["RNA"]]@data["Wnt7b",]; c2 <- s2[["RNA"]]@data["Wnt7b",]
  wnt <- c(sum(c1!=0), sum(c2!=0)); notwnt <- c(sum(c1==0), sum(c2==0))
  ft <- fisher.test(matrix(c(wnt, notwnt), nrow=2), alternative="less")
  fisher_results <- rbind(fisher_results,
    data.frame(cluster=paste(pair, collapse=" vs "),
               p_value=ft$p.value, odds_ratio=ft$estimate))
}
write.csv(fisher_results, file.path(opt$outdir,"fig_3E_fisher_tests.csv"), row.names=FALSE)

# Wnt7b DotPlot split by genotype
wnt7b.RZ  <- RZ.RZK[["RNA"]]@data["Wnt7b",] * (RZ.RZK$kras_genotype=="Wild_type")
wnt7b.RZK <- RZ.RZK[["RNA"]]@data["Wnt7b",] * (RZ.RZK$kras_genotype=="Mutant")
RZ.RZK[["NEW"]] <- CreateAssayObject(data=rbind(wnt7b.RZ, wnt7b.RZK))
DefaultAssay(RZ.RZK) <- "NEW"
Idents(RZ.RZK) <- factor(RZ.RZK@meta.data$clusters, levels=mylevel)
p3E_dot <- DotPlot(RZ.RZK, features=c("wnt7b.RZK","wnt7b.RZ"), dot.scale=10) +
  scale_colour_gradient2(low="dodgerblue3", mid="ghostwhite",
                         high="firebrick", limits=c(-2.5,2.5)) +
  xlab("Markers") + ylab("Clusters") + coord_flip() +
  theme(axis.text=element_text(size=15), legend.position="top")
ggsave(p3E_dot, filename=file.path(opt$outdir,"fig_3E_dotplot_wnt7b.pdf"), width=15, height=3)

# ==== FIGURE 3H — Coverage / Linkage plot for Wnt7b ==========================
cat(">>> Fig 3H\n")
p3H <- plot_linkage_cl_genotype(RZ.RZK, "macs2", "Wnt7b")
ggsave(p3H, filename=file.path(opt$outdir,"fig_3H_linkage_coverage.pdf"), width=8, height=8, dpi=300)

# Wilcoxon tests
wilcox_results <- data.frame(cluster=character(), p_value=numeric(), stringsAsFactors=FALSE)
for (i in 0:6) {
  id1 <- paste0(i,"_Wild_type"); id2 <- paste0(i,"_Mutant")
  c1 <- colnames(subset(RZ.RZK, subset=cl_genotype==id1))
  c2 <- colnames(subset(RZ.RZK, subset=cl_genotype==id2))
  DefaultAssay(RZ.RZK) <- "RNA"
  v1 <- RZ.RZK[["RNA"]]@data["Wnt7b",c1]; v2 <- RZ.RZK[["RNA"]]@data["Wnt7b",c2]
  wt <- wilcox.test(v1, v2)
  wilcox_results <- rbind(wilcox_results,
    data.frame(cluster=paste(id1,"vs",id2), p_value=wt$p.value))
}
write.csv(wilcox_results, file.path(opt$outdir,"fig_3H_wilcoxon_tests.csv"), row.names=FALSE)

# ATAC count enrichment tests for Wnt7b regulatory regions
regions <- c("chr15-85568400-85569900","chr15-85574150-85575200")
for (r in seq_along(regions)) {
  count_region <- CountsInRegion(object=RZ.RZK, assay="macs2",
                                 regions=StringToGRanges(regions[r]))
  enrich_res <- data.frame(cluster=character(), p_value=numeric(),
                           mean_mut=numeric(), mean_wt=numeric(), stringsAsFactors=FALSE)
  for (i in 0:6) {
    id_m <- paste0(i,"_Mutant"); id_w <- paste0(i,"_Wild_type")
    cells_m <- colnames(subset(RZ.RZK, subset=cl_genotype==id_m))
    cells_w <- colnames(subset(RZ.RZK, subset=cl_genotype==id_w))
    tt <- t.test(count_region[cells_m], count_region[cells_w])
    enrich_res <- rbind(enrich_res,
      data.frame(cluster=paste(id_m,"vs",id_w), p_value=tt$p.value,
                 mean_mut=tt$estimate[1], mean_wt=tt$estimate[2]))
  }
  write.csv(enrich_res, file.path(opt$outdir,
    paste0("fig_3H_atac_enrichment_region",r,".csv")), row.names=FALSE)
}

# ==== FIGURE S3A — FeaturePlots for all marker genes =========================
cat(">>> Fig S3A\n")
gene_cutoffs <- list(
  c("Muc5ac",0,3),c("Gkn2",0,3),c("Tff1",0.5,6),c("Gkn1",0.5,5),
  c("Smc2",0.5,3),c("Top2a",0.5,4),c("Hmgb2",1,4),c("Foxm1",0,2),
  c("Mki67",0.3,3),c("Muc6",0.3,4.5),c("Cftr",0.5,3.5),c("Cd44",2,4.2),
  c("Glipr1",1.5,4.5),c("Lgr5",0.5,3.5)
)
for (gc in gene_cutoffs) {
  p <- plot_featureplot_v3(RZ.RZK, gc[[1]], as.numeric(gc[[2]]), as.numeric(gc[[3]]))
  ggsave(p, filename=file.path(opt$outdir, paste0("fig_S3A_",gc[[1]],".pdf")),
         width=10, height=10)
}

# ==== FIGURE S3B — Wnt family VlnPlot ========================================
cat(">>> Fig S3B\n")
DefaultAssay(RZ.RZK) <- "RNA"
wnt_genes <- c("Wnt1","Wnt2","Wnt2b","Wnt3","Wnt3a","Wnt4","Wnt5a","Wnt5b",
               "Wnt6","Wnt7a","Wnt7b","Wnt8a","Wnt8b","Wnt9a","Wnt9b",
               "Wnt10a","Wnt10b","Wnt11","Wnt16")
pS3B <- VlnPlot(object=RZ.RZK, features=wnt_genes, group.by="kras_genotype")
ggsave(pS3B, filename=file.path(opt$outdir,"fig_S3B_Wnt_family.pdf"), width=16, height=16)

cat(">>> All figures saved to:", opt$outdir, "\n")