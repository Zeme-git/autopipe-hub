# snMultiome RZ vs RZK Analysis Pipeline

Reproduces the snMultiome (joint scRNA-seq + scATAC-seq) analysis from the paper:
**"Epithelial WNT secretion drives niche escape of developing gastric cancer"**

## What it does

Analyzes 10X Genomics snMultiome data from two mouse gastric tissue samples (RZ wild-type and RZK mutant) through a 7-step workflow: QC filtering → MACS2 peak calling → normalization (SCTransform + TF-IDF/SVD) → SCT-based integration and clustering → publication-quality visualization → differential gene expression with KRAS pathway overlap → motif enrichment analysis at the Wnt7b locus.

## Required Inputs

Place the following 10X Cell Ranger ARC output directories in the input folder:

```
input/
├── Mouse_RZ_WEN/
│   ├── filtered_feature_bc_matrix.h5
│   ├── atac_fragments.tsv.gz
│   ├── atac_fragments.tsv.gz.tbi
│   └── per_barcode_metrics.csv
└── Mouse_RZK_WEN/
    ├── filtered_feature_bc_matrix.h5
    ├── atac_fragments.tsv.gz
    ├── atac_fragments.tsv.gz.tbi
    └── per_barcode_metrics.csv
```

## Expected Outputs

| Directory | Files | Description |
|-----------|-------|-------------|
| `01_qc/` | `*_qc.rds`, QC VlnPlots | Filtered Seurat objects + QC violin plots |
| `02_peaks/` | `*_peaks.rds` | Objects with MACS2-called peaks |
| `03_norm/` | `*_norm.rds` | Normalized objects (SCT + TF-IDF) |
| `04_cluster/` | `combined_clustered.rds` | Integrated + clustered combined object |
| `05_visualization/` | Fig 3A-D, S3A-B PDFs | Publication figures (UMAP, DotPlot, etc.) |
| `06_differential/` | DEG CSVs | MAST DEGs + KRAS pathway overlap |
| `07_motif/` | Motif CSVs + barplots | DA peaks, motif enrichment, Wnt7b stats |

## Configuration

All parameters are in `config.yaml`, including QC thresholds, clustering resolutions, marker genes, and analysis parameters. Key QC thresholds:
- ATAC counts: 1,000–100,000
- RNA counts: 1,000–25,000
- Nucleosome signal < 2
- TSS enrichment > 1
- Mitochondrial % < 15

## How to Run

```bash
# Build Docker image
docker build -t snmultiome-rzrzk .

# Run pipeline
docker run --rm \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  snmultiome-rzrzk \
  conda run -n pipeline snakemake --cores 8 -s /pipeline/Snakefile
```