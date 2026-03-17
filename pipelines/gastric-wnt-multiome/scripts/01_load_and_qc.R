#!/usr/bin/env Rscript
# =============================================================================
#  Step 1: Load 10x Multiome Data & Quality Control
# =============================================================================

library(future)
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(yaml)

cfg <- yaml.load_file("/pipeline/config.yaml")
plan("multicore", workers = cfg$threads)
options(future.globals.maxSize = cfg$max_memory_gb * 1024^3)
set.seed(1234)

cat("Building gene annotation for mm10 ...\n")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

generate_object <- function(sample_cfg) {
    cat(sprintf("  Loading: %s\n", sample_cfg$h5))
    metadata   <- read.csv(sample_cfg$metadata, header = TRUE, row.names = 1)
    input.10x  <- Read10X_h5(sample_cfg$h5)
    rna.counts <- input.10x$`Gene Expression`
    obj <- CreateSeuratObject(counts = rna.counts, assay = "RNA", meta.data = metadata)
    mt.pattern <- ifelse(cfg$species == "mouse", "^mt-", "^MT-")
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt.pattern)
    atac.counts   <- input.10x$Peaks
    grange.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
    grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac.counts   <- atac.counts[as.vector(grange.use), ]
    obj[["ATAC"]] <- CreateChromatinAssay(
        counts = atac.counts, sep = c(":", "-"), genome = genome(annotations),
        fragments = sample_cfg$fragments, min.cells = 10, annotation = annotations)
    obj@meta.data$kras_genotype <- sample_cfg$genotype
    return(obj)
}

qc_and_filter <- function(obj) {
    DefaultAssay(obj) <- "ATAC"
    obj <- NucleosomeSignal(obj)
    obj <- TSSEnrichment(obj)
    qc <- cfg$qc
    obj <- subset(x = obj,
        subset = nCount_ATAC < qc$nCount_ATAC_max & nCount_RNA < qc$nCount_RNA_max &
                 nCount_ATAC > qc$nCount_ATAC_min & nCount_RNA > qc$nCount_RNA_min &
                 nucleosome_signal < qc$nucleosome_signal_max &
                 TSS.enrichment > qc$TSS_enrichment_min & percent.mt < qc$percent_mt_max)
    return(obj)
}

cat("Creating Wild-type object ...\n")
wild   <- generate_object(cfg$samples$Wild_type)
cat("Creating Mutant object ...\n")
mutant <- generate_object(cfg$samples$Mutant)
cat("Running QC and filtering ...\n")
wild   <- qc_and_filter(wild)
mutant <- qc_and_filter(mutant)
cat(sprintf("  Wild-type cells after QC: %d\n", ncol(wild)))
cat(sprintf("  Mutant cells after QC:    %d\n", ncol(mutant)))

out_dir <- cfg$output_dir
dir.create(file.path(out_dir, "rds"), showWarnings = FALSE, recursive = TRUE)
saveRDS(wild,   file.path(out_dir, "rds", "wild_qc.rds"))
saveRDS(mutant, file.path(out_dir, "rds", "mutant_qc.rds"))
cat("Step 1 complete.\n")