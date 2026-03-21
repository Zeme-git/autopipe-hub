#!/usr/bin/env Rscript
# Step 5: Visualization — Figures 3A-E, 3H, S3A-B + Statistics

options(Seurat.object.assay.version = "v3")

library(future); library(Seurat); library(Signac); library(ggplot2)
library(RColorBrewer); library(dittoSeq); library(BSgenome.Mmusculus.UCSC.mm10); library(yaml)

cfg <- yaml.load_file("/pipeline/config.yaml")
plan("multicore", workers = cfg$threads); options(future.globals.maxSize = cfg$max_memory_gb * 1024^3); set.seed(1234)
out_dir <- cfg$output_dir
dir.create(file.path(out_dir,"figures"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(out_dir,"tables"), showWarnings=FALSE, recursive=TRUE)
RZ.RZK <- readRDS(file.path(out_dir,"rds","RZ_RZK_annotated.rds"))
linkage_genome <- BSgenome.Mmusculus.UCSC.mm10

my_col <- c("0_Wild_type"="#8936EF","0_Mutant"="#8936EF","1_Wild_type"="#F2CA19","1_Mutant"="#F2CA19",
            "2_Wild_type"="#FF00BD","2_Mutant"="#FF00BD","3_Wild_type"="#E11845","3_Mutant"="#E11845",
            "4_Wild_type"="#0057E9","4_Mutant"="#0057E9","5_Wild_type"="#87E911","5_Mutant"="#87E911",
            "6_Wild_type"="#018300","6_Mutant"="#018300")
cluster_col <- c("Neck"="#8936EF","Lgr5+"="#F2CA19","SPEM"="#FF00BD","Wnt7+"="#E11845",
                 "Proliferating"="#0057E9","Pre-Pit"="#87E911","Pit"="#018300")
mylevel <- c("Lgr5+","SPEM","Neck","Proliferating","Pre-Pit","Pit","Wnt7+")

get_gene_expr <- function(obj, gene, assay_name = "RNA") {
    DefaultAssay(obj) <- assay_name
    tryCatch({
        return(GetAssayData(obj, slot = "data")[gene, ])
    }, error = function(e) {
        tryCatch({
            return(obj[[assay_name]]@data[gene, ])
        }, error = function(e2) {
            return(LayerData(obj, assay = assay_name, layer = "data")[gene, ])
        })
    })
}

plot_featureplot <- function(obj, gene) {
    DefaultAssay(obj) <- "RNA"
    p <- FeaturePlot(obj, features=gene, pt.size=1, min.cutoff=0.3, max.cutoff=2)
    ggplot(data=p$data, aes_string(x="UMAP_1", y="UMAP_2", fill=gene)) +
        geom_point(shape=21, stroke=0.3, color="black", alpha=1, size=2) +
        scale_fill_gradientn(colours=brewer.pal(9,"YlOrRd")) +
        theme_void() + theme(plot.title=element_text(hjust=0.5))
}

fig_dir <- file.path(out_dir,"figures")
stats_rows <- list()

cat("  Fig 3A ...\n")
p_3A <- DimPlot(RZ.RZK, group.by="cl_genotype", cols=my_col, pt.size=1.3) + NoLegend()
ggsave(p_3A, filename=file.path(fig_dir,"fig_3A_UMAP.pdf"), width=8, height=9)

cat("  Fig 3B ...\n")
DefaultAssay(RZ.RZK) <- "RNA"
existing_clusters <- unique(RZ.RZK@meta.data$clusters); existing_clusters <- existing_clusters[!is.na(existing_clusters)]
mylevel_filtered <- mylevel[mylevel %in% existing_clusters]
Idents(RZ.RZK) <- factor(RZ.RZK@meta.data$clusters, levels=mylevel_filtered)
existing_markers <- cfg$marker_genes[cfg$marker_genes %in% rownames(RZ.RZK)]
p_3B <- DotPlot(RZ.RZK, features=existing_markers, dot.scale=10) +
    scale_colour_gradient2(low="dodgerblue3", mid="ghostwhite", high="firebrick", limits=c(-2.5,2.5)) +
    xlab("Markers") + ylab("Clusters") + coord_flip() +
    theme(axis.text=element_text(size=15), legend.position="top")
ggsave(p_3B, filename=file.path(fig_dir,"fig_3B_dotplot.pdf"), width=6, height=7)

cat("  Fig 3C ...\n")
existing_cols <- cluster_col[names(cluster_col) %in% existing_clusters]
p_3C <- dittoBarPlot(RZ.RZK, var="clusters", group.by="kras_genotype", color.panel=existing_cols)
ggsave(p_3C, filename=file.path(fig_dir,"fig_3C_composition.pdf"), width=6, height=6)

cat("  Fig 3D ...\n")
tryCatch({
    p_wt <- plot_featureplot(subset(RZ.RZK, subset=kras_genotype=="Wild_type"), "Wnt7b")
    p_mut <- plot_featureplot(subset(RZ.RZK, subset=kras_genotype=="Mutant"), "Wnt7b")
    p_3D <- p_wt + p_mut
    ggsave(p_3D, filename=file.path(fig_dir,"fig_3D_Wnt7b_feature.pdf"), width=12, height=6, dpi=300)
}, error=function(e) cat(sprintf("    Fig 3D skipped: %s\n", e$message)))

cat("  Fig 3E ...\n")
all_cl_genotypes <- unique(RZ.RZK@meta.data$cl_genotype)
cl_list <- list(c("0_Wild_type","0_Mutant"), c("1_Wild_type","1_Mutant"),
                c("2_Wild_type","2_Mutant"), c("3_Wild_type","3_Mutant"),
                c("4_Wild_type","4_Mutant"), c("5_Wild_type","5_Mutant"),
                c("6_Wild_type","6_Mutant"))
pp <- c(); cl_vec <- c()
for(pair in cl_list) {
    wnt_counts <- c(); non_wnt <- c()
    for(j in pair) {
        if(!(j %in% all_cl_genotypes)) { pp<-c(pp,0); cl_vec<-c(cl_vec,j); wnt_counts<-c(wnt_counts,0); non_wnt<-c(non_wnt,0); next }
        sub_obj <- subset(RZ.RZK, subset=cl_genotype==j)
        counts <- get_gene_expr(sub_obj, "Wnt7b")
        gc_n <- sum(counts != 0); total <- length(counts); pct <- gc_n/total*100
        wnt_counts <- c(wnt_counts, gc_n); non_wnt <- c(non_wnt, total-gc_n)
        pp <- c(pp, pct); cl_vec <- c(cl_vec, j)
    }
    if(sum(wnt_counts)>0 && sum(non_wnt)>0) {
        tryCatch({
            fisher.df <- matrix(c(wnt_counts, non_wnt), nrow=2)
            fisher_res <- fisher.test(fisher.df, alternative="less")
            stats_rows[[length(stats_rows)+1]] <- data.frame(
                comparison=paste(pair, collapse=" vs "), test="Fisher_exact",
                p_value=fisher_res$p.value, statistic=fisher_res$estimate, stringsAsFactors=FALSE)
        }, error=function(e) cat(sprintf("    Fisher skipped: %s\n", e$message)))
    }
}
df_3E <- data.frame(Clusters=cl_vec, Percentage=pp)
df_3E <- df_3E[df_3E$Clusters %in% all_cl_genotypes, ]
p_3E <- ggplot(df_3E, aes(Clusters, Percentage, fill=Clusters)) +
    geom_bar(stat="identity") + ylim(0,100) + theme_classic() +
    scale_fill_manual(values=my_col) + ggtitle("Wnt7b-expressing cells (%)") +
    theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1)) + NoLegend()
ggsave(p_3E, filename=file.path(fig_dir,"fig_3E_Wnt7b_percent.pdf"), width=10, height=6, dpi=300)

cat("  Fig 3H ...\n")
tryCatch({
    DefaultAssay(RZ.RZK) <- "macs2"
    RZ.RZK <- RegionStats(RZ.RZK, genome=linkage_genome)
    RZ.RZK <- LinkPeaks(object=RZ.RZK, peak.assay="macs2", expression.assay="RNA", genes.use="Wnt7b")
    lv <- c(); for(i in c(1,2,0,4,5,6,3)) { lv <- c(lv, sprintf("%s_Wild_type",i), sprintf("%s_Mutant",i)) }
    lv <- lv[lv %in% all_cl_genotypes]; RZ.RZK$cl_genotype <- factor(RZ.RZK$cl_genotype, levels=lv)
    p_3H <- CoveragePlot(object=RZ.RZK, region="Wnt7b", features="Wnt7b",
                          expression.assay="RNA", group.by="cl_genotype", extend.upstream=0, extend.downstream=0)
    ggsave(p_3H, filename=file.path(fig_dir,"fig_3H_Wnt7b_coverage.pdf"), width=8, height=8, dpi=300)
}, error=function(e) {
    cat(sprintf("    Fig 3H skipped: %s\n", e$message))
    pdf(file.path(fig_dir,"fig_3H_Wnt7b_coverage.pdf")); plot.new(); title("Fig 3H - see log"); dev.off()
})

cat("  Fig 3H statistics ...\n")
for(cl in 0:6) {
    wt_id <- paste0(cl,"_Wild_type"); mut_id <- paste0(cl,"_Mutant")
    if(!(wt_id %in% all_cl_genotypes) || !(mut_id %in% all_cl_genotypes)) next
    tryCatch({
        wt_cells <- colnames(subset(RZ.RZK, subset=cl_genotype==wt_id))
        mut_cells <- colnames(subset(RZ.RZK, subset=cl_genotype==mut_id))
        if(length(wt_cells)<2 || length(mut_cells)<2) next
        wt_expr <- get_gene_expr(RZ.RZK, "Wnt7b")[wt_cells]; mut_expr <- get_gene_expr(RZ.RZK, "Wnt7b")[mut_cells]
        res <- wilcox.test(wt_expr, mut_expr)
        stats_rows[[length(stats_rows)+1]] <- data.frame(comparison=sprintf("Cluster%d_WT_vs_Mut_Wnt7b_expr",cl), test="Wilcoxon", p_value=res$p.value, statistic=NA, stringsAsFactors=FALSE)
    }, error=function(e) cat(sprintf("    Wilcoxon skipped cl %d: %s\n",cl,e$message)))
}
for(region_str in cfg$wnt7b_regions) {
    tryCatch({
        region_counts <- CountsInRegion(object=RZ.RZK, assay="macs2", regions=StringToGRanges(region_str))
        for(cl in 0:6) {
            wt_id <- paste0(cl,"_Wild_type"); mut_id <- paste0(cl,"_Mutant")
            if(!(wt_id %in% all_cl_genotypes) || !(mut_id %in% all_cl_genotypes)) next
            cells_m <- colnames(subset(RZ.RZK, subset=cl_genotype==mut_id)); cells_w <- colnames(subset(RZ.RZK, subset=cl_genotype==wt_id))
            if(length(cells_m)<2 || length(cells_w)<2) next
            res <- t.test(region_counts[cells_m], region_counts[cells_w])
            stats_rows[[length(stats_rows)+1]] <- data.frame(comparison=sprintf("Cluster%d_Mut_vs_WT_%s",cl,region_str), test="t_test", p_value=res$p.value, statistic=res$statistic, stringsAsFactors=FALSE)
        }
    }, error=function(e) cat(sprintf("    Region test skipped %s: %s\n",region_str,e$message)))
}

cat("  Fig S3B ...\n")
DefaultAssay(RZ.RZK) <- "RNA"
existing_wnt <- cfg$wnt_family_genes[cfg$wnt_family_genes %in% rownames(RZ.RZK)]
if(length(existing_wnt) > 0) {
    p_S3B <- VlnPlot(object=RZ.RZK, features=existing_wnt, group.by="kras_genotype")
    ggsave(p_S3B, filename=file.path(fig_dir,"fig_S3B_Wnt_family_violin.pdf"), width=20, height=12)
} else { pdf(file.path(fig_dir,"fig_S3B_Wnt_family_violin.pdf")); plot.new(); title("No Wnt genes found"); dev.off() }

if(length(stats_rows) > 0) { stats_df <- do.call(rbind, stats_rows) } else { stats_df <- data.frame(comparison=character(), test=character(), p_value=numeric(), statistic=numeric()) }
write.csv(stats_df, file.path(out_dir,"tables","statistics_summary.csv"), row.names=FALSE)
cat("Step 5 complete.\n")