# Figure 4 – NDD Gene × Perturbation Analysis

Run the notebooks in this folder **in order**:

## 1. HomologConvert.Rmd (run first)

**Homolog conversion: human gene annotations to mouse**

- Loads the gene annotation database and BioMart human–mouse homolog table.
- Converts human gene IDs to mouse and builds mouse gene sets (NDD, neuropsychiatry, neurodegenerative, control groups).
- Writes `ref_ndd_list_edit.RDS` under `DATA_DIR/ndd_perturbation/` for downstream use.

**Run this notebook first** so that the RDS file exists before running the figure panels.

## 2. figure4_panels.Rmd (run second)

**Figure 4 – NDD Gene × Perturbation Analysis**

- Uses the gene database (and any outputs from HomologConvert) to generate Figure 4 panels.
- Expects `DATA_DIR` to point to the analysis directory (e.g. `/data/wholebrain_crispr_atlas/analysis/`).
- Saves outputs under `./Figure4_output`.

---

## How to run

From the repo root (or from this folder), run in order:

```bash
cd notebooks/Figure_4

# 1. Homolog conversion (creates ref_ndd_list_edit.RDS)
Rscript -e 'rmarkdown::render("HomologConvert.Rmd")'

# 2. Figure 4 panels (writes to ./Figure4_output)
Rscript -e 'rmarkdown::render("figure4_panels.Rmd")'
```

From R or RStudio, you can instead run each `.Rmd` interactively, or render with:

```r
setwd("notebooks/Figure_4")  # or your path to this folder
rmarkdown::render("HomologConvert.Rmd")
rmarkdown::render("figure4_panels.Rmd")
```

