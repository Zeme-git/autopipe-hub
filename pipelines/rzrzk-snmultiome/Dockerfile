FROM condaforge/miniforge3:latest

# ---------- System libraries for R compilation --------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gfortran \
        libcurl4-openssl-dev libssl-dev libxml2-dev libpng-dev \
        libbz2-dev liblzma-dev libreadline-dev zlib1g-dev \
        libx11-dev libxt-dev libgl1-mesa-glx libhdf5-dev \
        libharfbuzz-dev libfribidi-dev libfreetype6-dev libjpeg-dev \
        libcairo2-dev libtiff5-dev libicu-dev libgsl-dev libv8-dev \
    && rm -rf /var/lib/apt/lists/*

# ---------- Conda: R 4.1 + core compiled packages + snakemake + macs2 --------
RUN conda install -y -c conda-forge -c bioconda \
        snakemake-minimal \
        bash \
        r-base=4.1 \
        r-rcpp r-rcpparmadillo r-rcppeigen \
        r-matrix r-lattice r-nlme r-mgcv r-mass \
        r-ggplot2 r-dplyr r-tidyr r-purrr r-tibble r-scales \
        r-hdf5r r-igraph r-rcolorbrewer \
        r-ggrepel r-gridextra r-cowplot \
        r-httr r-curl r-hmisc \
        r-cli r-crayon r-glue r-rlang r-vctrs r-lifecycle \
        r-digest r-systemfonts \
        r-readxl r-readr r-dendextend r-ggnewscale \
        r-optparse \
        macs2 \
    && conda clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# ---------- CRAN packages (Seurat v4, Signac, supporting pkgs) ---------------
RUN Rscript -e " \
    options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2023-01-01')); \
    install.packages(c('Seurat', 'SeuratObject', 'Signac', \
                       'clustree', 'hash', 'paletteer', 'future'))"

# ---------- BiocManager + Bioconductor 3.14 packages -------------------------
RUN Rscript -e " \
    if (!requireNamespace('BiocManager', quietly=TRUE)) \
        install.packages('BiocManager', repos='https://cran.r-project.org'); \
    BiocManager::install(version = '3.14', ask = FALSE); \
    BiocManager::install(c( \
        'EnsDb.Mmusculus.v79', \
        'BSgenome.Mmusculus.UCSC.mm10', \
        'biovizBase', \
        'motifmatchr', \
        'JASPAR2020', \
        'TFBSTools', \
        'chromVAR', \
        'dittoSeq', \
        'msigdbr', \
        'gprofiler2', \
        'MAST' \
    ), ask = FALSE, update = FALSE)"

# ---------- Pipeline files ----------------------------------------------------
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/