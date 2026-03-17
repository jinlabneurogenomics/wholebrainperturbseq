#!/usr/bin/env Rscript
# ============================================================
# Install R packages required for Whole Brain Perturb-seq
# (Figure 2, 3, 4 and related notebooks)
# ============================================================
# Run from repo root or environments/: Rscript install_r_packages.R
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Packages used across Figure 2, 3, 4 Rmd/notebooks
packages <- c(
  # tidyverse and core
  "tidyverse",   # dplyr, ggplot2, tidyr, forcats, readr, purrr
  "data.table",
  "readxl",
  "patchwork",
  "cowplot",
  "ggbeeswarm",
  "ggpubr",
  "ggrepel",
  "pheatmap",
  "viridis",
  "scales",
  # Figure 3 and general
  "arrow",
  "RColorBrewer",
  # Figure 1 (optional; some notebooks use pak)
  "Matrix",
  "reshape2",
  # Rendering
  "rmarkdown",
  "knitr"
)

# biomaRt is from Bioconductor (not CRAN) — install first so it's available for all R versions
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  cat("Installing BiocManager and biomaRt (Bioconductor)...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  }
  BiocManager::install("biomaRt", update = FALSE, ask = FALSE)
}

cat("Installing", length(packages), "R packages from CRAN...\n")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing", pkg, "...\n")
    tryCatch(
      install.packages(pkg, dependencies = TRUE, quiet = TRUE),
      error = function(e) {
        cat("  ERROR installing", pkg, ":", conditionMessage(e), "\n")
        quit(save = "no", status = 1)
      }
    )
  } else {
    cat("  ", pkg, "already installed\n")
  }
}

# Verify critical packages load (fail script if any missing so Docker build fails)
cat("\nVerifying packages load:\n")
required <- c("tidyverse", "data.table", "readxl", "patchwork", "pheatmap", "ggpubr", "rmarkdown")
failed <- character(0)
for (pkg in required) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat("  ", pkg, ":", if (ok) "OK" else "FAILED", "\n")
  if (!ok) failed <- c(failed, pkg)
}
if (length(failed) > 0) {
  cat("\nERROR: Required packages failed to load:", paste(failed, collapse = ", "), "\n")
  quit(save = "no", status = 1)
}

cat("\nDone.\n")
