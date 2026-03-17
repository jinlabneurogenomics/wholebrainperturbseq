# Whole Brain Perturb-seq

Analysis code and notebooks for a **whole-brain in vivo CRISPR perturbation screen**. This repository reproduces manuscript figures and supplementary data tables.

## Docker

Build and run the analysis environment in a container (includes R for Figure 2–4 Rmd and Python with rapids-singlecell for notebooks and scripts).

**Build the image:**

```bash
docker build -t wholebrainperturbseq .
```

**Run a shell with the repo and data mounted (CPU):**

```bash
docker run -it --rm \
  -v $(pwd):/workspace/wholebrainperturbseq \
  -v /path/to/your/data:/data \
  wholebrainperturbseq bash
```

**Run with GPU support** (for Figure 1 GPU-accelerated steps). Requires [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) (nvidia-docker) installed on the host.

```bash
docker run -it --rm --gpus all \
  -v $(pwd):/workspace/wholebrainperturbseq \
  -v /path/to/your/data:/data \
  wholebrainperturbseq bash
```

- **Repo mount:** `-v $(pwd):/workspace/wholebrainperturbseq` — your local repo is the working directory; scripts and notebooks use this path.
- **Data mount:** `-v /path/to/your/data:/data` — directory for data download, outputs, figures, etc. Replace `/path/to/your/data` with the host path to your data directory.

## Data download

To download the whole-brain CRISPR atlas h5ad shards from HuggingFace, use [`scripts/download.py`](scripts/download.py):

```bash
python scripts/download.py --output-dir /path/to/output
```

- **`--output-dir`** (required) — root directory where data will be written; the script creates an `h5ads/` subdirectory with the shards.
- **`--hf-token`** (optional) — HuggingFace token; only needed for gated or private datasets (the default repo is public).

Example with the data mount at `/data`:

```bash
python scripts/download.py --output-dir /data/wholebrain_crispr_atlas
```

## Supplementary data

Supplementary data tables and figure-level inputs (e.g. for Figures 2–4 and Data S4–S7) can be downloaded here:

- [Supplementary data (link to be added)](PLACEHOLDER_SUPPLEMENTARY_DATA_URL)

Place the downloaded files in the paths expected by each notebook (e.g. `Figure4_Codes/Figure4_inputfiles/`, or the `data/` paths referenced in the Rmd files).

## Preprocessing for downstream analysis

[`scripts/processing_script.sh`](scripts/processing_script.sh) prepares the downloaded h5ad shards for downstream steps such as **differential gene expression (DEGs)** and **E-distance**. It (1) processes each batch and extracts guide/expression data, (2) applies QC filters per batch, (3) merges batches into one AnnData, and (4) runs final QC on the merged object.

**Prerequisite:** Downloaded data at `/data/wholebrain_crispr_atlas/h5ads/` (e.g. from the [Data download](#data-download) step with `--output-dir /data/wholebrain_crispr_atlas`).

**Run from the repository root** (e.g. inside the Docker container):

```bash
cd scripts && bash processing_script.sh
```

The script reads from `BASE_DATA_PATH` (set to `/data/wholebrain_crispr_atlas` in the script) and writes to `processed/` under that path. The final merged, filtered object is `processed/WB8588_screen_gex_filtered.h5ad`, which you can use for DEG notebooks and E-distance analysis.

## Notebooks

Notebooks are under `notebooks/`, grouped by figure. Rendered **HTML** outputs are available for some notebooks. Detailed steps for **Figure 1** are in [`notebooks/Figure_1/README.md`](notebooks/Figure_1/README.md).

| Notebook | Outputs | Description | Rendered HTML |
|----------|---------|-------------|---------------|
| **Figure 1** ([`notebooks/Figure_1/`](notebooks/Figure_1/)) | Figure 1 B, G, H, I, K | Python notebooks (UMAP 1b, validation, cell-type distribution 1H) and R **Rmd** ([`04-perturbaion_distribution_Figure_1_GIK.rmd`](notebooks/Figure_1/04-perturbaion_distribution_Figure_1_GIK.rmd)) for panels G, I, K. See linked README for scripts and config. | [`04-perturbaion_distribution_Figure_1_GIK.nb.html`](notebooks/Figure_1/04-perturbaion_distribution_Figure_1_GIK.nb.html) |
| [`cell_type_differential_abundance_figure2_s5.ipynb`](notebooks/Figure_2/Figure_2_celltype_depletion_analysis_a-b-c/cell_type_differential_abundance_figure2_s5.ipynb) | Figure 2 A/B/C, Figure S5 A/B/C, Data S3 | Fisher’s exact test (differential abundance), MAGeCK2, fitness loss / log-odds heatmaps. | — |
| [`DEGs_computation.ipynb`](notebooks/DEGs_computation.ipynb) | Figure S6 B, Data S4 A/B | Wilcoxon DE per cell type; potency/magnitude heatmaps; exports supplementary DE table. | — |
| [`figure2_panels.ipynb`](notebooks/Figure_2/Figure_2_e-f-g-h-i/figure2_panels.ipynb) | Figure 2 E, F, G, H, I | Panels from supplementary data: top-50 DEG barplot, Jaccard heatmap, shared-DEG LFC heatmap, Tsc1/Tsc2 scatter, SWI/SNF heatmap. | — |
| [`pairwise_edistance_figure3_and_data_s5.ipynb`](notebooks/Figure_3/pairwise_edistance_figure3_and_data_s5.ipynb) | Figure 3A, Figure S8, Data S5 | Pairwise E-distance per cell type (and whole brain); clustermaps. | — |
| [`figure3_panels.Rmd`](notebooks/Figure_3/figure3_panels.Rmd) | Figure 3 D, E, F | Grin/Gria analysis: effect-size dotplot, Grin2a/Grin2b correlation, shared-DEG heatmap (L2-3 IT CTX Glut). | [`figure3_panels.html`](notebooks/Figure_3/figure3_panels.html), [`figure3_panels_old.html`](notebooks/Figure_3/figure3_panels_old.html) |
| [`figure4_panels.Rmd`](notebooks/Figure_4/figure4_panels.Rmd) | Figure 4, Data S6/S7 | NDD gene × perturbation analysis; disorder-stratified DEG burden and risk-gene enrichment per cell type. | — |
