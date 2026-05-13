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
python download.py --output-dir /data/wholebrain_crispr_atlas
```


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
