# gastric-wnt-multiome

Reproduces the single-nucleus multiome (joint snRNA-seq + scATAC-seq) analysis from **"Epithelial WNT secretion drives niche escape of developing gastric cancer"**, generating Figures 3A-I and S3A-C.

## v1.0.1 Changes
- **Dockerfile**: 2-stage conda install with bioconda pre-built binaries for Bioconductor packages (Rhtslib, Rsamtools, rtracklayer, BSgenome, CNEr, TFBSTools) — resolves GCC 15 compilation failures
- **Seurat v5 compatibility**: All R scripts use `options(Seurat.object.assay.version = "v3")` for stable v4-style integration workflow
- **Robust visualization**: Dynamic cluster handling, removed hardcoded reorder indices, tryCatch error handling
- **Validated**: Successfully tested end-to-end with real 10x snMultiome data

## Required Inputs

Mount your data directory at `/input` with two sample directories:

```
/input/
├── Mouse_RZ_WEN/              # Wild-type (RZ)
│   ├── filtered_feature_bc_matrix.h5
│   ├── atac_fragments.tsv.gz
│   ├── atac_fragments.tsv.gz.tbi
│   └── per_barcode_metrics.csv
└── Mouse_RZK_WEN/             # KRAS-mutant (RZK)
    ├── filtered_feature_bc_matrix.h5
    ├── atac_fragments.tsv.gz
    ├── atac_fragments.tsv.gz.tbi
    └── per_barcode_metrics.csv
```

## Expected Outputs

| Directory | Contents |
|-----------|----------|
| `rds/` | Seurat objects at each stage (QC, MACS2, normalized, annotated) |
| `figures/` | PDFs: UMAP (3A), DotPlot (3B), composition (3C), Wnt7b feature (3D), Wnt7b percentage (3E), coverage (3H), motifs (3I, S3C), Wnt family violin (S3B) |
| `tables/` | CSVs: DEG list (MAST), DA peaks (mutant & WT), statistics summary |
| `logs/` | Per-step log files |

## Pipeline Steps

1. **load_and_qc** — Load 10x h5, compute ATAC QC, filter cells
2. **macs2_peaks** — Call peaks with MACS2, add macs2 assay
3. **normalize** — SCTransform + LogNorm (RNA); TF-IDF + SVD (ATAC)
4. **integrate_cluster** — SCT integration, PCA, UMAP, clustering, annotation
5. **visualization** — Figures 3A-E, 3H, S3B + statistical tests
6. **differential** — DEG (MAST), DA peaks (LR), KRAS hallmark overlap
7. **motif_enrichment** — JASPAR2020 motifs, Figures 3I, S3C

## How to Run

```bash
docker build --no-cache -t gastric-wnt-multiome .
docker run --rm \
    -v /path/to/data:/input:ro \
    -v /path/to/output:/output \
    gastric-wnt-multiome \
    conda run --no-capture-output -n pipeline \
    snakemake --cores 20 --snakefile /pipeline/Snakefile
```

## Configuration

Edit `config.yaml` to adjust QC thresholds, clustering resolution, DE/DA parameters, gene lists, thread count, and memory limits.