#!/usr/bin/env Rscript
# =============================================================================
#  01_load_and_qc.R
#  Load 10X snMultiome data (RNA + ATAC), compute QC metrics, filter cells
#  Extracted from: 1_RZRZK_basic_analysis.R
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
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

# ── Load configuration ───────────────────────────────────────────────────────
cfg <- read_yaml(config_file)

# ── Set up parallelism ────────────────────────────────────────────────────────
plan("multicore", workers = cfg$n_workers)
options(future.globals.maxSize = cfg$future_max_size_gb * 1024^3)
future.seed <- TRUE
set.seed(cfg$seed)

cat("=== Step 1: Load data and QC ===\n")

# ── Helper: Set genome annotations (mouse mm10) ──────────────────────────────
set_annotation <- function(species) {
  if (species == "mouse") {
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotation) <- "UCSC"
    genome(annotation) <- "mm10"
  }
  return(annotation)
}

# ── Helper: Create Seurat object from 10X h5 + ATAC fragments ────────────────
generate_object <- function(h5_file, frag_file, meta_file, species, annotations) {
  cat("  Loading:", h5_file, "\n")
  
  # Read per-barcode metadata from Cell Ranger
  metadata <- read.csv(file = meta_file, header = TRUE, row.names = 1)
  
  # Read 10X h5 (contains both Gene Expression and Peaks)
  input.10x <- Read10X_h5(h5_file)
  
  # --- RNA assay ---
  rna.counts <- input.10x$`Gene Expression`
  obj <- CreateSeuratObject(counts = rna.counts, assay = "RNA", meta.data = metadata)
  
  # Mitochondrial percentage (mouse: "^mt-")
  if (species == "mouse") {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  } else {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  }
  
  # --- ATAC assay (standard chromosomes only) ---
  atac.counts <- input.10x$Peaks
  grange.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac.counts <- atac.counts[as.vector(grange.use), ]
  
  obj[["ATAC"]] <- CreateChromatinAssay(
    counts    = atac.counts,
    sep       = c(":", "-"),
    genome    = genome(annotations),
    fragments = frag_file,
    min.cells = 10,
    annotation = annotations
  )
  
  return(obj)
}

# ── Helper: QC metrics (nucleosome signal + TSS enrichment) ───────────────────
qc_function <- function(obj) {
  DefaultAssay(obj) <- "ATAC"
  obj <- NucleosomeSignal(obj)
  obj <- TSSEnrichment(obj)
  return(obj)
}

# ── Helper: Subset cells based on QC thresholds ──────────────────────────────
subsetting_objects_QC <- function(obj, cfg) {
  obj <- qc_function(obj)
  obj <- subset(
    x = obj,
    subset = nCount_ATAC < cfg$qc$nCount_ATAC_max &
             nCount_RNA  < cfg$qc$nCount_RNA_max &
             nCount_ATAC > cfg$qc$nCount_ATAC_min &
             nCount_RNA  > cfg$qc$nCount_RNA_min &
             nucleosome_signal < cfg$qc$nucleosome_signal_max &
             TSS.enrichment    > cfg$qc$TSS_enrichment_min &
             percent.mt        < cfg$qc$percent_mt_max
  )
  return(obj)
}

# ── Main execution ────────────────────────────────────────────────────────────

# 1. Set genome annotation
annotations <- set_annotation(cfg$species)

# 2. Create Seurat objects for each sample
wild <- generate_object(
  cfg$samples$wild$h5_file,
  cfg$samples$wild$fragment_file,
  cfg$samples$wild$meta_file,
  cfg$species,
  annotations
)
mutant <- generate_object(
  cfg$samples$mutant$h5_file,
  cfg$samples$mutant$fragment_file,
  cfg$samples$mutant$meta_file,
  cfg$species,
  annotations
)

# 3. Add genotype metadata
wild@meta.data$kras_genotype   <- cfg$samples$wild$genotype
mutant@meta.data$kras_genotype <- cfg$samples$mutant$genotype

cat("  Wild cells before QC:", ncol(wild), "\n")
cat("  Mutant cells before QC:", ncol(mutant), "\n")

# 4. QC violin plots (before filtering)
dir.create("/output/01_qc", showWarnings = FALSE, recursive = TRUE)

pdf("/output/01_qc/qc_vlnplot_wild.pdf", width = 12, height = 6)
DefaultAssay(wild) <- "ATAC"
wild_qc <- qc_function(wild)
print(VlnPlot(wild_qc, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 4, pt.size = 0))
dev.off()

pdf("/output/01_qc/qc_vlnplot_mutant.pdf", width = 12, height = 6)
DefaultAssay(mutant) <- "ATAC"
mutant_qc <- qc_function(mutant)
print(VlnPlot(mutant_qc, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 4, pt.size = 0))
dev.off()

# 5. Filter cells
wild   <- subsetting_objects_QC(wild, cfg)
mutant <- subsetting_objects_QC(mutant, cfg)

cat("  Wild cells after QC:", ncol(wild), "\n")
cat("  Mutant cells after QC:", ncol(mutant), "\n")

# 6. Save filtered objects
saveRDS(wild,   "/output/01_qc/wild_qc.rds")
saveRDS(mutant, "/output/01_qc/mutant_qc.rds")

cat("=== Step 1 complete ===\n")