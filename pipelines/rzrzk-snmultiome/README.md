# RZ vs RZK snMultiome Analysis Pipeline

Reproduces the snMultiome (joint snRNA-seq + scATAC-seq) analysis from **"Epithelial WNT secretion drives niche escape of developing gastric cancer"**, comparing wild-type (RZ) and KRAS-mutant (RZK) mouse gastric tissue.

## Pipeline Steps

| Step | Rule | Script | Description |
|------|------|--------|-------------|
| 1 | `load_and_qc` | 01_load_and_qc.R | Load 10x .h5, QC filter (ATAC/RNA counts, nucleosome signal, TSS enrichment, %MT) |
| 2 | `macs2_peaks` | 02_macs2_peaks.R | MACS2 peak calling + mm10 blacklist removal |
| 3 | `normalize` | 03_normalize.R | SCTransform + LogNorm (RNA), TF-IDF + SVD (ATAC & MACS2) |
| 4 | `integrate_cluster` | 04_integrate_cluster.R | SCT integration, PCA → UMAP, clustering (res 0.3 & 0.6), cell-type annotation |
| 5 | `visualization` | 05_visualization.R | Figures 3A-E, 3H, S3B + statistics (Fisher, Wilcoxon, t-tests) |
| 6 | `differential` | 06_differential.R | DEG (MAST) + KRAS hallmark overlap + DA peaks (Wnt7+ vs Lgr5+) |
| 7 | `motif_enrichment` | 07_motif_enrichment.R | JASPAR2020 motif enrichment → Figures 3I, S3C |

**Important:** All scripts use `options(Seurat.object.assay.version = "v3")` for Seurat v5 → v4 compatibility mode.

## Required Inputs

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

```
output/
├── rds/                            # Seurat objects (intermediate + final)
├── figures/
│   ├── fig_3A_UMAP.pdf
│   ├── fig_3B_dotplot.pdf
│   ├── fig_3C_composition.pdf
│   ├── fig_3D_Wnt7b_feature.pdf
│   ├── fig_3E_Wnt7b_percent.pdf
│   ├── fig_3H_Wnt7b_coverage.pdf
│   ├── fig_3I_motifs_RZK.pdf
│   ├── fig_S3B_Wnt_family_violin.pdf
│   └── fig_S3C_motifs_RZ.pdf
├── tables/
│   ├── RZRZK_DEG_scMAST.csv
│   ├── DA_peaks_Wnt7_vs_Lgr5_mutant.csv
│   ├── DA_peaks_Wnt7_vs_Lgr5_wildtype.csv
│   └── statistics_summary.csv
└── logs/
```

## How to Run

```bash
docker build -t rzrzk-snmultiome .

docker run --rm \
    -v /path/to/data:/input:ro \
    -v /path/to/output:/output \
    rzrzk-snmultiome \
    snakemake --snakefile /pipeline/Snakefile --cores 20 -p
```

## Cell Type Annotations

| Cluster | Cell Type | Key Markers |
|---------|-----------|-------------|
| 0 | Neck | Muc6, Cftr |
| 1 | Lgr5+ | Lgr5 |
| 2 | SPEM | Glipr1 |
| 3 | Wnt7+ | Wnt7b |
| 4 | Proliferating | Mki67, Foxm1 |
| 5 | Pre-Pit | Tff1, Gkn1 |
| 6 | Pit | Muc5ac, Gkn2 |