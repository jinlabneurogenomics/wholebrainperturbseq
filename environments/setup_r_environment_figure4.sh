#!/usr/bin/env bash
# ============================================================
# R Environment Setup for Figure 4 Analysis
# ============================================================
# Run this script once to install R and all required packages.
# Tested on Ubuntu/Debian. Adjust for macOS (brew) if needed.
# ============================================================

set -e

# ---- 1. Install R (Ubuntu/Debian) -------------------------
if ! command -v R &>/dev/null; then
  echo ">>> Installing R..."
  sudo apt-get update -y
  sudo apt-get install -y --no-install-recommends \
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
    libgit2-dev
  echo ">>> R installed: $(R --version | head -1)"
else
  echo ">>> R already installed: $(R --version | head -1)"
fi

# ---- 2. Install R packages --------------------------------
echo ">>> Installing R packages..."
Rscript install_r_packages.R

echo ""
echo "============================================================"
echo " Setup complete. Open Figure4_notebook.Rmd in RStudio or"
echo " render with: Rscript -e \"rmarkdown::render('Figure4_notebook.Rmd')\""
echo "============================================================"
