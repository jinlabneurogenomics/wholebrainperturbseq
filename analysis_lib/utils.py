from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Literal, Tuple

import scanpy as sc
from anndata import AnnData
import numpy as np
import pandas as pd
from tqdm import tqdm

import re


def deduplicate_gene_symbols(
    adata, gene_symbol_col="gene_symbols", gene_id_col="gene_ids"
):
    """
    Deduplicate features in adata so that each gene_symbol maps to only one gene_id.
    Keeps the row with the smallest numeric portion of the gene_id for each gene_symbol.
    Operates on columns 'gene_symbols' and 'gene_ids' in adata.var.

    Parameters
    ----------
    adata : AnnData
        The AnnData object to deduplicate.

    Returns
    -------
    AnnData
        Deduplicated AnnData object (new copy if deduplication was needed, else original object).
    """
    # Preconditions & quick checks
    assert (
        gene_symbol_col in adata.var.columns
    ), f"adata.var must contain '{gene_symbol_col}'"
    assert gene_id_col in adata.var.columns, f"adata.var must contain '{gene_id_col}'"

    gene_symbols = adata.var[gene_symbol_col]
    print(f"{gene_symbol_col} is unique:", gene_symbols.is_unique)

    if not gene_symbols.is_unique:
        var_df = adata.var.copy()

        # Robust numeric extraction from Ensembl gene ids (e.g., ENSMUSG00000123456 -> 123456)
        def numeric_id(gid):
            s = str(gid) if pd.notna(gid) else ""
            m = re.search(r"(\d+)$", s)
            return int(m.group(1)) if m else np.inf

        var_df["numeric_id"] = var_df[gene_id_col].map(numeric_id)

        # Show only the duplicated gene_symbol rows (for logging)
        dup_mask = var_df[gene_symbol_col].duplicated(keep=False)
        print(f"Non-unique {gene_symbol_col} values (prior to drop):")
        print(var_df.loc[dup_mask, [gene_symbol_col, gene_id_col]])

        # For each gene_symbol group, keep the row with the smallest numeric_id
        grp = var_df.groupby(gene_symbol_col, sort=False, group_keys=False)
        idx_keep = grp["numeric_id"].idxmin()

        # Preserve the ORIGINAL order of features in adata.var_names
        keep_set = set(idx_keep.tolist())
        features_to_keep = [name for name in adata.var_names if name in keep_set]

        # Optional: inspect which features will be dropped
        dropped = [name for name in adata.var_names if name not in keep_set]
        if dropped:
            print(
                f"Dropping {len(dropped)} duplicate feature(s). Example:", dropped[:10]
            )

        # Subset adata to keep only the selected features (public API; updates X, var, layers, obsm/varm shapes)
        adata = adata[:, features_to_keep].copy()

        # No need to overwrite adata.var; the helper column lived only in var_df
        print("After de-duplication, n_vars:", adata.n_vars)
        adata.var_names = adata.var[gene_symbol_col]
    return adata


def get_processed_data_for_deg(
    data_path: str,
    *,
    gene_symbol_col: str = "gene_symbols",
    gene_id_col: str = "gene_ids",
    mt_prefix: str = "mt-",
    min_cells_per_gene: int = 0,
    min_genes_per_cell: int = 0,
    max_pct_mt: Optional[float] = None,
    min_total_counts: Optional[int] = None,
    max_total_counts: Optional[int] = None,
    subclass_col: str = "subclass_name",
    target_gene_col: str = "gene_target",
    min_cells_per_subclass: int = 50,
    remove_multiple_guides: bool = True,
    remove_unassigned: bool = True,
    unassigned_label: str = "Negative",
    remove_doublets: bool = True,
    doublet_column: str = "scDblFinder.class",
    log_ambient_mse_norm_min: Optional[float] = None,
    verbose: bool = True,
    filter_clusters: Optional[Sequence[str]] = ["1", "2", "3", "6", "17", "57", "83", "NA"],
) -> AnnData:
    """
    Load and preprocess Perturbai adata for downstream DEG analysis, following the steps:
        - Load `.h5ad`
        - Filter genes by cells expressed
        - Compute QC metrics and filter cells by mitochondrial percentage and total counts
        - Remove rows with unassigned target gene
        - Filter out rare subclasses and optionally subset to top-N or specified subclasses

    Parameters
    ----------
    data_path : str
        Path to the `.h5ad` file.
    mt_prefix : str, default "mt-"
        Prefix used to flag mitochondrial genes.
    min_cells_per_gene : int, default 50
        Keep genes expressed in at least this many cells.
    min_genes_per_cell : int, default 0
        Keep cells with at least this many genes detected.
    max_pct_mt : Optional[float], default 10.0
        Maximum allowed percent mitochondrial counts per cell.
    min_total_counts : Optional[int], default 1000
        Minimum allowed total counts per cell.
    max_total_counts : Optional[int], default 10000
        Maximum allowed total counts per cell.
    subclass_col : str, default "subclass_name"
        Column indicating cell type/subclass.
    target_gene_col : str, default "target_gene"
        Column indicating target gene; rows with `unassigned_label` are removed if present.
    channel_col : str, default "channel"
        Column indicating channel; used for optional mouse mapping.
    min_cells_per_subclass : int, default 50
        Keep subclasses with at least this many cells.
    top_n_cell_types : Optional[int], default 10
        If set, subset to the top-N most abundant subclasses.
    restrict_cell_types : Optional[Sequence[str]]
        If provided, subset to this explicit list of subclasses.
    log_ambient_mse_norm_min : Optional[float], default None
        If set, keep only cells with obs["log_ambient_mse_norm"] > this value.
    remove_multiple_guides : bool, default True
        If set, remove cells with multiple guides.
    remove_doublets : bool, default True
        If set, remove doublets.
    verbose : bool, default True
        Print progress and filtering summaries.

    Returns
    -------
    AnnData
        The processed and optionally subset AnnData object suitable for DEG.
    """

    if verbose:
        print("Loading data...")

    with tqdm(
        total=1, desc="Reading h5ad metadata", unit="file", disable=not verbose
    ) as pbar:
        adata_backed = sc.read_h5ad(data_path, backed="r")
        pbar.update(1)

    if verbose:
        print(f"Dataset: {adata_backed.n_obs:,} cells, {adata_backed.n_vars:,} genes")

    # Identify cells to keep based on metadata filters before loading full data
    cell_mask = np.ones(adata_backed.n_obs, dtype=bool)
    if remove_unassigned and target_gene_col in adata_backed.obs.columns:
        unassigned_mask = adata_backed.obs[target_gene_col] != unassigned_label
        cell_mask &= unassigned_mask
        if verbose:
            print(f"  Pre-filtering unassigned: {(~unassigned_mask).sum():,} cells")
    # Early subclass filter
    if subclass_col in adata_backed.obs.columns and min_cells_per_subclass > 0:
        subclass_counts = adata_backed.obs[subclass_col].value_counts()
        valid_subclasses = subclass_counts[
            subclass_counts >= min_cells_per_subclass
        ].index
        subclass_mask = adata_backed.obs[subclass_col].isin(valid_subclasses)
        cell_mask &= subclass_mask
        if verbose:
            print(f"  Pre-filtering rare subclasses: {(~subclass_mask).sum():,} cells")
    if remove_doublets:
        cell_mask &= adata_backed.obs[doublet_column] == "singlet"
        if verbose:
            print(f"  Pre-filtering doublets: {(~cell_mask).sum():,} cells")
    if log_ambient_mse_norm_min is not None and "log_ambient_mse_norm" in adata_backed.obs.columns:
        ambient_mask = adata_backed.obs["log_ambient_mse_norm"] > log_ambient_mse_norm_min
        cell_mask &= ambient_mask
        if verbose:
            print(f"  Pre-filtering log_ambient_mse_norm <= {log_ambient_mse_norm_min}: {(~ambient_mask).sum():,} cells")
    if remove_multiple_guides and "num_guides" in adata_backed.obs.columns:
        guides_mask = adata_backed.obs["num_guides"] == 1
        cell_mask &= guides_mask
        if verbose:
            print(f"  Pre-filtering num_guides > 1: {(~guides_mask).sum():,} cells")
    if min_genes_per_cell is not None and min_genes_per_cell > 0:
        cell_mask &= adata_backed.obs["num_genes"] >= min_genes_per_cell
        if verbose:
            print(f"  Pre-filtering num_genes < {min_genes_per_cell}: {(~cell_mask).sum():,} cells")
    if filter_clusters is not None:
        cluster_mask = ~adata_backed.obs["cluster"].isin(filter_clusters)
        cell_mask &= cluster_mask
        if verbose:
            print(f"  Pre-filtering clusters: {(~cluster_mask).sum():,} cells")
            
    cells_to_keep = cell_mask.sum()
    if verbose:
        print(
            f"  Loading {cells_to_keep:,} cells after pre-filtering ({cells_to_keep/adata_backed.n_obs*100:.1f}% of total)"
        )

    # Load only the filtered cells into memory
    with tqdm(
        total=1, desc="Loading filtered data", unit="step", disable=not verbose
    ) as pbar:
        adata = adata_backed[cell_mask].to_memory()
        pbar.update(1)

    if verbose:
        print(f"Loaded data: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    # Filter genes and cells
    if min_cells_per_gene and min_cells_per_gene > 0:
        with tqdm(
            total=1, desc="Filtering genes", unit="step", disable=not verbose
        ) as pbar:
            sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
            pbar.update(1)
        if verbose:
            print(
                f"After gene filter (min_cells={min_cells_per_gene}): {adata.n_vars:,} genes"
            )

    adata = deduplicate_gene_symbols(
        adata, gene_symbol_col=gene_symbol_col, gene_id_col=gene_id_col
    )

    # Identify mitochondrial genes using robust symbol detection
    symbols = adata.var[gene_symbol_col].astype(str).str.lower()
    adata.var["mt"] = symbols.str.startswith(mt_prefix.lower())

    # QC metrics
    if ("total_counts" not in adata.obs.columns) or ("pct_counts_mt" not in adata.obs.columns) or ("n_genes_by_counts" not in adata.obs.columns):
        with tqdm(
            total=1, desc="Computing QC metrics", unit="step", disable=not verbose
        ) as pbar:
            sc.pp.calculate_qc_metrics(
                adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
            )
            pbar.update(1)

    # Mitochondrial filter
    if "pct_counts_mt" in adata.obs.columns and max_pct_mt is not None:
        before = adata.n_obs
        adata = adata[adata.obs["pct_counts_mt"] < float(max_pct_mt)].copy()
        if verbose:
            removed = before - adata.n_obs
            print(
                f"Filtered by pct_counts_mt < {max_pct_mt}: removed {removed:,} cells; kept {adata.n_obs:,}"
            )

    # Total counts filters
    if "total_counts" not in adata.obs.columns:
        adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).flatten()
    if min_total_counts is not None:
        before = adata.n_obs
        adata = adata[adata.obs["total_counts"] >= int(min_total_counts)].copy()
        if verbose:
            removed = before - adata.n_obs
            print(
                f"Filtered by total_counts >= {min_total_counts}: removed {removed:,} cells; kept {adata.n_obs:,}"
            )
    if max_total_counts is not None:
        before = adata.n_obs
        adata = adata[adata.obs["total_counts"] < int(max_total_counts)].copy()
        if verbose:
            removed = before - adata.n_obs
            print(
                f"Filtered by total_counts < {max_total_counts}: removed {removed:,} cells; kept {adata.n_obs:,}"
            )

    return adata
