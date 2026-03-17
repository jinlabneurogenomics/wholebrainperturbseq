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
python scripts/download.py --output-dir /data/WholeBrainPerturbSeq
```

## Supplementary data

Supplementary data tables and figure-level inputs (e.g. for Figures 2–4 and Data S4–S7) can be downloaded here:

- [Supplementary data (link to be added)](PLACEHOLDER_SUPPLEMENTARY_DATA_URL)

Place the downloaded files in the paths expected by each notebook (e.g. `Figure4_Codes/Figure4_inputfiles/`, or the `data/` paths referenced in the Rmd files).

## Preprocessing for downstream analysis

[`scripts/processing_script.sh`](scripts/processing_script.sh) prepares the downloaded h5ad shards for downstream steps such as **differential gene expression (DEGs)** and **E-distance**. It (1) processes each batch and extracts guide/expression data, (2) applies QC filters per batch, (3) merges batches into one AnnData, and (4) runs final QC on the merged object.

**Prerequisite:** Downloaded data at `/data/WholeBrainPerturbSeq/h5ads/` (e.g. from the [Data download](#data-download) step with `--output-dir /data/WholeBrainPerturbSeq`).

**Run from the repository root** (e.g. inside the Docker container):

```bash
cd scripts && bash processing_script.sh
```

The script reads from `BASE_DATA_PATH` (set to `/data/WholeBrainPerturbSeq` in the script) and writes to `processed/` under that path. The final merged, filtered object is `processed/WB8588_screen_gex_filtered.h5ad`, which you can use for DEG notebooks and E-distance analysis.

## Notebooks

| Notebook | Outputs | Description |
|----------|---------|-------------|
| **cell_type_differential_abundance_figure2_s5** | Figure 2 A/B/C, Figure S5 A/B/C, Data 1 | Fisher’s exact test (differential abundance), MAGeCK2, fitness loss / log-odds heatmaps |
| **deg_analysis_figure_s6_and_data_s4** | Figure S6 B, Data S4 A/B | Wilcoxon DE per cell type, potency/magnitude heatmaps, supplementary DE table |
| **pairwise_edistance_figure3_and_data_s5** | Figure 3A, Figure S8, Data S5 | Pairwise E-distance per cell type (and whole brain), clustermaps |
| **NDD_gene_perturbation_analysis_figure4** | Figure 4, Data S6/S7 | Disorder-stratified DEG burden and downstream risk gene enrichment per cell type |
