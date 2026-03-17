# ============================================================
# Whole Brain Perturb-seq — Analysis environment
# ============================================================
# Base: scverse rapids_singlecell (GPU single-cell stack already installed).
# Adds R (for Figure 2, 3, 4 Rmd) and repo-specific Python deps.
# Build: docker build -t wholebrainperturbseq .
# Run:   docker run -it --rm -v $(pwd):/workspace/wholebrainperturbseq wholebrainperturbseq bash
# ============================================================

FROM ghcr.io/scverse/rapids_singlecell:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV PYTHONPATH=/workspace/wholebrainperturbseq
ENV WORKDIR=/workspace/wholebrainperturbseq

# ---- 1. System and R (CRAN) ---------------------------------
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    ca-certificates \
    gnupg \
    wget \
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && echo "deb [arch=amd64] https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" > /etc/apt/sources.list.d/cran_r.list \
    && apt-get update -y \
    && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# ---- 2. R packages (Figure 2, 3, 4) -------------------------
COPY environments/install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R \
    && rm /tmp/install_r_packages.R

# ---- 3. Extra Python deps (base has scanpy, anndata, pertpy, rapids_singlecell) ----
RUN pip install --no-cache-dir huggingface_hub

# ---- 4. Working directory and entrypoint ---------------------
WORKDIR ${WORKDIR}

# Copy repo into image (override with bind mount at run time if preferred)
COPY . ${WORKDIR}

CMD ["/bin/bash"]
