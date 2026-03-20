# RZ vs RZK snMultiome Analysis Pipeline

Reproduces the snMultiome (joint snRNA-seq + scATAC-seq) analysis from **"Epithelial WNT secretion drives niche escape of developing gastric cancer"**, comparing wild-type (RZ) and KRAS-mutant (RZK) mouse gastric tissue.

## Pipeline Steps

| Step | Rule | Description |
|------|------|-------------|
| 1-5 | `process_sample` | Load 10x .h5, QC filter, MACS2 peak calling, SCTransform + LogNorm (RNA), TF-IDF + SVD (ATAC) |
| 6-7 | `integrate_and_cluster` | SCT-based integration, PCA → UMAP, clustering at res 0.3 & 0.6 |
| 8 | `deg_analysis` | DEG (MAST) Mutant vs WT; HALLMARK_KRAS_SIGNALING_UP overlap |
| 9 | `visualize` | Figures 3A-3E, 3H, S3A, S3B (UMAP, DotPlot, dittoBar, FeaturePlot, Coverage) |
| 10 | `da_motif` | Differential accessibility (Wnt7+ vs Lgr5+) + JASPAR2020 motif enrichment |

## Required Inputs

Place under a single input directory with this structure:

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

All files are standard 10x Genomics CellRanger ARC outputs.

## Expected Outputs

```
output/
├── processed_wild.rds              # Per-sample processed objects
├── processed_mutant.rds
├── clustered_object.rds            # Integrated + clustered object
├── RZRZK_DEG_scMAST.csv           # DEG results
├── RZRZK_DEG_scMAST_kras_summary.csv
├── da_peaks_mutant.csv             # DA peaks (Wnt7+ vs Lgr5+)
├── da_peaks_wildtype.csv
├── motifs_enriched_mutant.csv      # JASPAR motif enrichment
├── motifs_enriched_wildtype.csv
├── figures/
│   ├── fig_3A_UMAP.pdf
│   ├── fig_3B_dotplot.pdf
│   ├── fig_3C_dittobarplot.pdf
│   ├── fig_3D_Wnt7b_feature.pdf
│   ├── fig_3E_Wnt7b_pct.pdf
│   ├── fig_3H_linkage_coverage.pdf
│   ├── fig_3I_motifs_Wnt_vs_Lgr5_RZK.pdf
│   ├── fig_S3A_*.pdf               # Per-gene FeaturePlots
│   ├── fig_S3B_Wnt_family.pdf
│   └── fig_S3C_motifs_Wnt_vs_Lgr5_RZ.pdf
└── logs/
```

## How to Run

```bash
# Build the Docker image
docker build -t rzrzk-snmultiome .

# Run the pipeline
docker run --rm \
    -v /path/to/your/data:/input:ro \
    -v /path/to/output:/output \
    rzrzk-snmultiome \
    conda run --no-capture-output -n base \
    snakemake --cores 20 --snakefile /pipeline/Snakefile
```

## Configuration

Edit `config.yaml` to adjust:

- **QC thresholds**: ATAC/RNA count limits, nucleosome signal, TSS enrichment, %MT
- **Clustering**: PCA dimensions, resolutions (0.3 initial, 0.6 refined)
- **DE parameters**: test method (MAST), adjusted p-value / LFC cutoffs
- **DA parameters**: min.pct, LFC threshold, test (LR), motif p-value cutoff
- **Threads / memory**: parallelisation settings

## Cell Type Annotations

| Cluster | Cell Type | Key Markers |
|---------|-----------|-------------|
| 0 | Neck | Muc6, Cftr |
| 1 | Lgr5+ | Lgr5 |
| 2 | SPEM | Glipr1 |
| 3 | Wnt7+ | Wnt7b |
| 4 | Proliferating | Mki67, Foxm1, Stmn1 |
| 5 | Pre-Pit | Tff1, Gkn1 |
| 6 | Pit | Muc5ac, Gkn2 |