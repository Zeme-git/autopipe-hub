# RZRZK snMultiome Pipeline

Reproduces the single-nucleus multiome (joint snRNA-seq + scATAC-seq) analysis from **"Epithelial WNT secretion drives niche escape of developing gastric cancer"**, generating Figures 3A–I and S3A–C.

## Required Inputs

Mount your data directory at `/input`. It must contain two sample directories:

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

These are standard 10x Genomics Cell Ranger ARC outputs.

## Expected Outputs

Written to `/output`:

| Directory | Contents |
|-----------|----------|
| `rds/` | Seurat objects at each stage (QC, MACS2, normalized, final annotated) |
| `figures/` | PDFs: UMAP (3A), DotPlot (3B), composition (3C), Wnt7b feature (3D), Wnt7b percentage (3E), coverage (3H), motifs (3I, S3C), violins (S3B), individual markers (S3A/) |
| `tables/` | CSVs: DEG list (MAST), DA peaks (Wnt7+ vs Lgr5+ in mutant & WT), statistics summary |
| `logs/` | Per-step log files |

## Pipeline Steps

1. **load_and_qc** — Load 10x h5, compute ATAC QC, filter cells
2. **macs2_peaks** — Call peaks with MACS2, add macs2 assay
3. **normalize** — SCTransform + LogNorm (RNA); TF-IDF + SVD (ATAC)
4. **integrate_cluster** — SCT integration → PCA → UMAP → clustering → annotation
5. **visualization** — Figures 3A–E, 3H, S3A–B + statistical tests
6. **differential** — DEG (MAST), DA peaks (LR), KRAS hallmark overlap
7. **motif_enrichment** — JASPAR2020 motifs → Figures 3I, S3C

## How to Run

```bash
# Build the Docker image
docker build -t rzrzk-multiome .

# Run the pipeline
docker run --rm \
    -v /path/to/your/data:/input:ro \
    -v /path/to/output:/output \
    rzrzk-multiome \
    snakemake --cores 20 --snakefile /pipeline/Snakefile
```

## Configuration

Edit `config.yaml` to adjust:

- **QC thresholds**: ATAC/RNA count bounds, nucleosome signal, TSS enrichment, %MT
- **Clustering**: PCA dimensions, resolution (coarse 0.3, fine 0.6)
- **DE parameters**: MAST test, p-value and log2FC cutoffs
- **DA parameters**: LR test, min.pct, logfc threshold
- **Parallelization**: thread count, max memory
- **Gene lists**: marker genes, Wnt family, feature plot specifications