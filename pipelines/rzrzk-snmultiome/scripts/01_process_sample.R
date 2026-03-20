#!/usr/bin/env Rscript
###############################################################################
#  Step 1-5: Per-sample processing
#  Load 10x h5 → QC → MACS2 peak calling → RNA + ATAC processing
###############################################################################
suppressPackageStartupMessages({
  library(optparse)
  library(future)
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

# ---- CLI arguments -----------------------------------------------------------
option_list <- list(
  make_option("--h5",        type="character"),
  make_option("--fragments", type="character"),
  make_option("--metadata",  type="character"),
  make_option("--species",   type="character", default="mouse"),
  make_option("--genotype",  type="character"),
  make_option("--nCount_ATAC_min",        type="integer",  default=1000),
  make_option("--nCount_ATAC_max",        type="integer",  default=100000),
  make_option("--nCount_RNA_min",         type="integer",  default=1000),
  make_option("--nCount_RNA_max",         type="integer",  default=25000),
  make_option("--nucleosome_signal_max",  type="double",   default=2),
  make_option("--TSS_enrichment_min",     type="double",   default=1),
  make_option("--percent_mt_max",         type="double",   default=15),
  make_option("--nfeatures",              type="integer",  default=2000),
  make_option("--macs2_path",             type="character", default="/opt/conda/bin/macs2"),
  make_option("--threads",                type="integer",  default=20),
  make_option("--max_memory_gb",          type="integer",  default=80),
  make_option("--output",                 type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---- Setup -------------------------------------------------------------------
plan("multicore", workers = opt$threads)
options(future.globals.maxSize = opt$max_memory_gb * 1024^3)
future.seed <- TRUE
set.seed(1234)

# ---- Gene annotation ---------------------------------------------------------
cat(">>> Setting up gene annotation ...\n")
if (opt$species == "mouse") {
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "mm10"
  mt_pattern <- "^mt-"
  blacklist <- blacklist_mm10
} else {
  stop("Only mouse (mm10) is currently supported.")
}

# ---- Create Seurat object ----------------------------------------------------
cat(">>> Loading 10x data:", opt$h5, "\n")
metadata   <- read.csv(file = opt$metadata, header = TRUE, row.names = 1)
input.10x  <- Read10X_h5(opt$h5)
rna.counts <- input.10x$`Gene Expression`
obj <- CreateSeuratObject(counts = rna.counts, assay = "RNA", meta.data = metadata)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)

# ATAC assay (standard chromosomes only)
atac.counts   <- input.10x$Peaks
grange.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac.counts   <- atac.counts[as.vector(grange.use), ]

obj[["ATAC"]] <- CreateChromatinAssay(
  counts     = atac.counts,
  sep        = c(":", "-"),
  genome     = genome(annotations),
  fragments  = opt$fragments,
  min.cells  = 10,
  annotation = annotations
)

# ---- Add genotype metadata ---------------------------------------------------
obj@meta.data$kras_genotype <- opt$genotype
cat("   Cells loaded:", ncol(obj), " genotype:", opt$genotype, "\n")

# ---- Quality control & filtering --------------------------------------------
cat(">>> QC filtering ...\n")
DefaultAssay(obj) <- "ATAC"
obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)

obj <- subset(
  x = obj,
  subset = nCount_ATAC       < opt$nCount_ATAC_max &
           nCount_RNA        < opt$nCount_RNA_max  &
           nCount_ATAC       > opt$nCount_ATAC_min &
           nCount_RNA        > opt$nCount_RNA_min  &
           nucleosome_signal < opt$nucleosome_signal_max &
           TSS.enrichment    > opt$TSS_enrichment_min    &
           percent.mt        < opt$percent_mt_max
)
cat("   Cells after QC:", ncol(obj), "\n")

# ---- MACS2 peak calling -----------------------------------------------------
cat(">>> MACS2 peak calling ...\n")
peaks <- CallPeaks(obj, macs2.path = opt$macs2_path, verbose = TRUE)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist, invert = TRUE)

macs2_counts <- FeatureMatrix(
  fragments = Fragments(obj),
  features  = peaks,
  cells     = colnames(obj)
)

obj[["macs2"]] <- CreateChromatinAssay(
  counts     = macs2_counts,
  fragments  = opt$fragments,
  annotation = annotations
)

# ---- RNA processing (SCTransform + LogNormalize) ----------------------------
cat(">>> RNA processing ...\n")
obj <- SCTransform(obj)
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = opt$nfeatures)

# ---- ATAC processing (TF-IDF + SVD on both peak sets) ----------------------
cat(">>> ATAC processing ...\n")
DefaultAssay(obj) <- "ATAC"
obj <- FindTopFeatures(obj, min.cutoff = "q0")
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)

DefaultAssay(obj) <- "macs2"
obj <- FindTopFeatures(obj, min.cutoff = "q0")
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)

# ---- Save processed object ---------------------------------------------------
cat(">>> Saving to:", opt$output, "\n")
saveRDS(obj, file = opt$output)
cat(">>> Done.\n")