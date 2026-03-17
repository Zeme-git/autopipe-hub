#!/usr/bin/env Rscript
# Step 6: Differential Expression & Differential Accessibility

library(future)
library(Seurat)
library(Signac)
library(dplyr)
library(MAST)
library(msigdbr)
library(gprofiler2)
library(yaml)

cfg <- yaml.load_file("/pipeline/config.yaml")
plan("multicore", workers = cfg$threads)
options(future.globals.maxSize = cfg$max_memory_gb * 1024^3)
set.seed(1234)

out_dir <- cfg$output_dir
RZ.RZK <- readRDS(file.path(out_dir, "rds", "RZ_RZK_annotated.rds"))

cat("Running MAST DE: Mutant vs Wild_type ...\n")
DefaultAssay(RZ.RZK) <- "RNA"
DEG <- FindMarkers(RZ.RZK, group.by="kras_genotype", ident.1="Mutant", test.use=cfg$de$test_use)
write.csv(DEG, file.path(out_dir,"tables","RZRZK_DEG_scMAST.csv"))
DEG_sig <- DEG[DEG$p_val_adj < cfg$de$padj_threshold & DEG$avg_log2FC > cfg$de$log2fc_threshold, ]
cat(sprintf("  Significant up-regulated DEGs: %d\n", nrow(DEG_sig)))

cat("Computing KRAS Hallmark overlap ...\n")
hallmark <- msigdbr(species="Homo sapiens",category="H") %>% dplyr::select(gs_name,gene_symbol)
kras_up <- hallmark$gene_symbol[hallmark$gs_name=="HALLMARK_KRAS_SIGNALING_UP"]
human_orth <- gconvert(query=rownames(DEG_sig),organism="hsapiens",target="ENSG")
overlap <- intersect(human_orth$name,kras_up)
cat(sprintf("  KRAS-UP: %d | DEGs: %d | Overlap: %d\n",length(kras_up),length(human_orth$name),length(overlap)))

cat("Running DA peaks (LR test) ...\n")
DefaultAssay(RZ.RZK) <- "ATAC"
plan("multicore", workers=1)
da_args <- list(min.pct=cfg$da$min_pct, logfc.threshold=cfg$da$logfc_threshold, test.use=cfg$da$test_use, latent.vars=cfg$da$latent_vars)

da_mut <- do.call(FindMarkers, c(list(object=RZ.RZK,ident.1="3_Mutant",ident.2="1_Mutant",group.by="cl_genotype"),da_args))
da_wt  <- do.call(FindMarkers, c(list(object=RZ.RZK,ident.1="3_Wild_type",ident.2="1_Wild_type",group.by="cl_genotype"),da_args))

annotate_da <- function(obj, da_peaks) {
    closest <- ClosestFeature(obj, regions=rownames(da_peaks))
    rownames(closest) <- closest$query_region
    da_peaks$gene <- closest$gene_name; da_peaks$gene_biotype <- closest$gene_biotype; da_peaks$distance <- closest$distance
    return(da_peaks)
}
da_mut <- annotate_da(RZ.RZK, da_mut)
da_wt  <- annotate_da(RZ.RZK, da_wt)

write.csv(da_mut, file.path(out_dir,"tables","DA_peaks_Wnt7_vs_Lgr5_mutant.csv"))
write.csv(da_wt,  file.path(out_dir,"tables","DA_peaks_Wnt7_vs_Lgr5_wildtype.csv"))
cat("Step 6 complete.\n")