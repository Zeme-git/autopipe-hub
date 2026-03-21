# Gastric WNT snMultiome Pipeline

Reproduces the single-nucleus multiome (joint snRNA-seq + scATAC-seq) analysis from **"Epithelial WNT secretion drives niche escape of developing gastric cancer"**, comparing wild-type (RZ) and KRAS-mutant (RZK) mouse gastric tissue.

## Pipeline Steps

| Step | Rule | Script | Description |
|------|------|--------|-------------|
| 1 | `load_and_qc` | 01_load_and_qc.R | Load 10x .h5, QC filter |
| 2 | `macs2_peaks` | 02_macs2_peaks.R | MACS2 peak calling + mm10 blacklist removal |
| 3 | `normalize` | 03_normalize.R | SCTransform + LogNorm (RNA), TF-IDF + SVD (ATAC & MACS2) |
| 4 | `integrate_cluster` | 04_integrate_cluster.R | SCT integration, PCA → UMAP, clustering, cell-type annotation |
| 5 | `visualization` | 05_visualization.R | Figures 3A-E, 3H, S3B + statistics |
| 6 | `differential` | 06_differential.R | DEG (MAST) + KRAS hallmark overlap + DA peaks |
| 7 | `motif_enrichment` | 07_motif_enrichment.R | JASPAR2020 motif enrichment → Figures 3I, S3C |

All scripts use `options(Seurat.object.assay.version = "v3")` for Seurat v5 → v4 compatibility.

## How to Run

```bash
docker build -t gastric-wnt-multiome .

docker run --rm \
    -v /path/to/data:/input:ro \
    -v /path/to/output:/output \
    gastric-wnt-multiome \
    snakemake --snakefile /pipeline/Snakefile --cores 20 -p
```

## Required Inputs

Place under a single directory with `Mouse_RZ_WEN/` and `Mouse_RZK_WEN/` subdirectories, each containing 10x CellRanger ARC outputs: `filtered_feature_bc_matrix.h5`, `atac_fragments.tsv.gz`, `atac_fragments.tsv.gz.tbi`, `per_barcode_metrics.csv`.

## Cell Types

| Cluster | Cell Type | Key Markers |
|---------|-----------|-------------|
| 0 | Neck | Muc6, Cftr |
| 1 | Lgr5+ | Lgr5 |
| 2 | SPEM | Glipr1 |
| 3 | Wnt7+ | Wnt7b |
| 4 | Proliferating | Mki67, Foxm1 |
| 5 | Pre-Pit | Tff1, Gkn1 |
| 6 | Pit | Muc5ac, Gkn2 |