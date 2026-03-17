# Figure 1 — Whole-Brain CRISPR Atlas Overview

Reproduces panels B, G, H, I, K from Figure 1

## Requirements

```
pip install huggingface_hub anndata zarr numpy scipy tqdm scanpy
```

An R environment with `Seurat`, `ggplot2`, and `patchwork` is needed for the `.rmd` panels.

## Steps

### 1. Download and preprocess

**Step 1a — Download:** Run `scripts/download.py` to fetch 14 h5ad shards from HuggingFace (`perturbai/wholebrain_crispr_atlas`) into `output-dir/h5ads/`.

```bash
python download.py --output-dir /data/wholebrain
```

**Step 1b — Assemble zarr and preprocess:** Run `scripts/umap_preprocess.py` to assemble `raw_wholebrain.zarr` from the h5ads, create `preprocess_config.yaml`, and run filtering and log-normalisation.

```bash
python umap_preprocess.py \
  --output-dir /data/wholebrain \
  --use-genes-path /data/use_genes.txt
```

Use `--skip-zarr` if the raw zarr is already built; use `--skip-preprocess` to only write the config. Outputs: `raw_wholebrain.zarr` and `filtered_lognorm.zarr` under `--output-dir`.

**Step 1c — Impute and visualize:** Run `impute_visualize.py` with `umap_config.yaml` to run hierarchical KNN imputation, UMAP, and category plots. The config sets `data_path`, `out_dir`, and `use_genes_path`; ensure they match the paths from steps 1a–1b (e.g. `filtered_lognorm.zarr` and your use-genes file). This script will create `results.h5ad` file needed for the next step. Outputs go to `out_dir` (e.g. `umap_results/`).

```bash
python impute_visualize.py umap_config.yaml
```

### 2. UMAP embedding — panel B

Open and run `01-UMAP_plot_Figure_1b_6M.ipynb`. Point `DATA_PATH` at `results.h5ad` created from Step 1c. Generates the 6M-cell UMAP coloured by cell type.

### 3. Validate paper numbers

Run `02-validate_paper_numbers.ipynb` to reproduce cell counts and guide-assignment statistics cited in the manuscript.

### 4. Cell-type distribution — panel H

Run `03-cells_disitribution_Figure_1H.ipynb`. Produces the stacked bar chart of cell-type proportions across mice.

### 5. Perturbation distribution — panels G, I, K

The following R markdown notebook covers the remaining panels:

```bash
Rscript -e "rmarkdown::render('04-perturbaion_distribution_Figure_1_GIK.rmd')"
```

These plot per-perturbation cell counts, guide efficiency, and target-gene coverage across brain regions.

## Output files

All figures are saved to `../Figures/Figure_1/` by default. Adjust the `out_dir` variable at the top of each notebook/script to change the destination.
