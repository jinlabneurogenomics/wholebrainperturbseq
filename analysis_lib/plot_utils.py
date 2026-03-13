import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def summarize_de_genes_per_condition(
    df_sig: pd.DataFrame,
    min_de_genes: int,
    lfc_threshold: float = 1.0,
    plot: bool = True,
    verbose: bool = True,
    group_name_col: str = "group_name",
    target_gene_col: str = "target_gene",
    savepath: str | None = None,
    dpi: int = 150,
    format: str | None = None,
) -> tuple:
    """
    Summarize number of DE genes and cells per condition (group_name x target_gene).
    Adds mean absolute log fold change (LFC) for each condition.

    Args:
        df_sig: Filtered dataframe of DE results (columns include names, logfoldchanges,
            n_pert_matched, and the group_name_col / target_gene_col).
        min_de_genes: Minimum number of DE genes per condition to highlight/report.
        lfc_threshold: Threshold for |logfoldchange| for a gene to count as DE (default=1.0).
        plot: Whether to plot n_de_genes vs n_cells scatter (default True).
        verbose: Whether to print counts (default True).
        group_name_col: Column name for cell type / group (default "group_name").
        target_gene_col: Column name for target gene (default "target_gene").
        savepath: If provided and plot=True, save figure to this path.
        dpi: DPI for saved figure (default 150).
        format: Image format when saving (default None, inferred from savepath).

    Returns:
        Tuple of (per_condition_df, n_conditions, n_target_genes, n_group_names).
        per_condition_df has columns condition, n_de_genes, n_cells, mean_abs_lfc,
        plus group_name_col and target_gene_col.
    """
    df_sig = df_sig.copy()
    df_sig["condition"] = df_sig[group_name_col] + " x " + df_sig[target_gene_col].astype(str)
    df_sig_filtered = df_sig[df_sig["logfoldchanges"].abs() >= lfc_threshold]

    de_genes_per_condition = df_sig_filtered.groupby("condition")["names"].nunique()
    n_cells_per_condition = df_sig_filtered.groupby("condition")["n_pert_matched"].max()
    mean_abs_lfc_per_condition = df_sig_filtered.groupby("condition")["logfoldchanges"].apply(
        lambda x: x.abs().mean()
    )

    per_condition_df = pd.DataFrame({
        "n_de_genes": de_genes_per_condition,
        "n_cells": n_cells_per_condition,
        "mean_abs_lfc": mean_abs_lfc_per_condition,
    }).reset_index()

    per_condition_df[[group_name_col, target_gene_col]] = per_condition_df["condition"].str.split(
        " x ", n=1, expand=True
    )

    to_show = per_condition_df[per_condition_df["n_de_genes"] >= min_de_genes]
    n_conditions = len(to_show)
    if verbose:
        print("Number of conditions with at least ", min_de_genes, " DE genes: ", n_conditions)

    n_target_genes = len(to_show[target_gene_col].value_counts())
    n_group_names = len(to_show[group_name_col].value_counts())
    if verbose:
        print("Number of unique target genes: ", n_target_genes)
        print("Number of unique group names: ", n_group_names)

    if plot:
        from scipy.stats import gaussian_kde

        x = per_condition_df["n_de_genes"].values
        y = per_condition_df["n_cells"].values
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

        plt.scatter(x, y, c=z, cmap="viridis")
        plt.xlabel("n_de_genes")
        plt.ylabel("n_cells")
        plt.title("DE genes per condition vs n_cells")
        plt.colorbar(label="Density")
        if savepath:
            plt.savefig(savepath, bbox_inches="tight", dpi=dpi, format=format)

    return per_condition_df, n_conditions, n_target_genes, n_group_names
