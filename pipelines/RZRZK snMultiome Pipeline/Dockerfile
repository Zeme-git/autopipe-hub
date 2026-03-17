# =============================================================================
#  Dockerfile — RZ vs RZK snMultiome Pipeline (AutoPipe)
# =============================================================================
FROM condaforge/miniforge3:latest

# ---- System dependencies for R graphics & HDF5 ------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libpng-dev libbz2-dev liblzma-dev libreadline-dev \
        libx11-dev libxt-dev zlib1g-dev libgl1-mesa-glx \
        libhdf5-dev libv8-dev \
        libharfbuzz-dev libfribidi-dev libfreetype6-dev \
        libjpeg-dev libcairo2-dev libtiff5-dev libicu-dev \
        build-essential gfortran procps \
    && apt-get clean

# ---- Conda: R 4.1, snakemake, bash, core R packages -------------------------
RUN conda install -y -c conda-forge -c bioconda \
        snakemake-minimal \
        bash \
        r-base=4.1 \
        r-yaml r-rcpp r-rcpparmadillo r-rcppeigen \
        r-dplyr r-tidyr r-purrr r-tibble r-rlang r-vctrs \
        r-ggplot2 r-scales r-ggrepel r-gridextra \
        r-ggforce r-ggraph r-tidygraph r-igraph r-graphlayouts \
        r-colorspace r-viridis r-viridislite r-rcolorbrewer \
        r-hdf5r r-matrix r-mgcv r-mass r-lattice r-nlme \
        r-glue r-cli r-crayon r-digest r-fansi r-utf8 \
        r-lifecycle r-pillar r-pkgconfig r-generics r-ellipsis \
        r-systemfonts r-bit r-bit64 r-isoband \
        r-labeling r-farver r-tweenr r-polyclip r-munsell \
        r-r6 r-withr r-cpp11 r-gtable \
        r-curl r-httr r-hmisc \
        macs2 \
    && conda clean -afy

# Replace system bash with conda bash to prevent GLIBC mismatch
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# ---- CRAN packages (Seurat v4 + extras) -------------------------------------
RUN Rscript -e " \
    options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2023-01-01')); \
    install.packages(c('Seurat', 'clustree', 'hash', 'paletteer', \
                        'readxl', 'ggnewscale', 'cowplot', 'dendextend'))"

# ---- Bioconductor packages --------------------------------------------------
RUN Rscript -e " \
    if (!requireNamespace('BiocManager', quietly = TRUE)) \
        install.packages('BiocManager'); \
    BiocManager::install(version = '3.14', ask = FALSE); \
    BiocManager::install(c( \
        'Signac', \
        'biovizBase', \
        'EnsDb.Mmusculus.v79', \
        'BSgenome.Mmusculus.UCSC.mm10', \
        'dittoSeq', \
        'motifmatchr', \
        'JASPAR2020', \
        'TFBSTools', \
        'chromVAR', \
        'msigdbr', \
        'gprofiler2', \
        'MAST' \
    ), ask = FALSE)"

# ---- Copy pipeline files -----------------------------------------------------
COPY Snakefile config.yaml /pipeline/
COPY scripts/ /pipeline/scripts/

WORKDIR /pipeline