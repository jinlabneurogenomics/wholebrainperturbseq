import os
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed
from typing import Optional, Union
from anndata import AnnData
from tqdm import tqdm


def get_matched_adata(
    adata_group: AnnData,
    pert_cells: pd.DataFrame,
    ctrl_cells: pd.DataFrame,
    do_normalization: bool = True,
) -> AnnData:
    """
    Prepare a combined adata from perturbation and control cells for DEG analysis.

    Optionally performs normalization and log-transform.
    """
    sel_idx = list(pert_cells.index) + list(ctrl_cells.index)
    adata_sub = adata_group[sel_idx].copy()
    # Use object dtype to avoid AnnData auto-converting to categorical (and logging about it)
    deg_group = pd.Series("control", index=adata_sub.obs_names, dtype=object)
    deg_group.loc[pert_cells.index] = "treatment"
    adata_sub.obs["deg_group"] = pd.Categorical(deg_group, categories=["control", "treatment"])

    if do_normalization:
        # target norm and log transform for DEG
        sc.pp.normalize_total(adata_sub, target_sum=1e4)
        sc.pp.log1p(adata_sub)

    adata_sub.obs_names_make_unique()
    return adata_sub


def run_wilcoxon_de_per_cell_type(
    adata: sc.AnnData,
    *,
    subclass_col: str = "group_name",
    labels: str = "gene_target",
    control_labels: list[str] = ["Non_target"],
    min_cells_per_group: int = 20,
    out_dir_root: Optional[str] = None,
    out_dir: Optional[str] = None,
    method: str = "wilcoxon",
    n_jobs: int = 1,
    match_libsize: bool = False,
) -> Union[pd.DataFrame, dict, None]:
    """Run Wilcoxon DE per cell type (perturbation vs control) and return combined results.

    For each value of ``adata.obs[subclass_col]``, filters targets and controls by
    ``min_cells_per_group``, builds a combined AnnData of perturbation and control
    cells (no matching), and runs ``scanpy.tl.rank_genes_groups(method=method)``.
    Results are concatenated across all cell types.

    If multiple ``control_labels`` are provided, returns a dict mapping control
    label → combined DataFrame; if a single control is provided, returns the
    combined DataFrame.

    Args:
        adata: AnnData object with counts and metadata in ``obs``.
        subclass_col: ``obs`` column containing the cell-type labels to iterate over.
        labels: ``obs`` column containing target/control labels (e.g. ``"gene_target"``).
        control_labels: List of labels considered controls.
        min_cells_per_group: Minimum number of cells required per label within a cell type.
        out_dir_root: Optional directory to create (no outputs are written).
        out_dir: Alias for out_dir_root; if provided, overrides out_dir_root.
        method: Scanpy rank_genes_groups method (default ``"wilcoxon"``).
        n_jobs: Parallel workers for target contrasts (1 = sequential).
        match_libsize: Ignored; retained for API compatibility (no matching is performed).

    Returns:
        - If ``len(control_labels) == 1``: a single combined ``pd.DataFrame``.
        - Else: a ``dict[str, pd.DataFrame]`` mapping each control label to its results.
    """
    if out_dir is not None:
        out_dir_root = out_dir
    if out_dir_root is not None:
        os.makedirs(out_dir_root, exist_ok=True)

    # Initialize accumulator per control label
    results_by_control: dict[str, list[pd.DataFrame]] = {
        ctrl: [] for ctrl in control_labels
    }

    # Iterate all cell types
    subclasses = [c for c in pd.unique(adata.obs[subclass_col]) if pd.notna(c)]
    for ct in tqdm(sorted(subclasses), desc="Processing cell types"):
        ad_sub = adata[adata.obs[subclass_col] == ct].copy()
        # Filter genes expressed in less than 10 cells
        gene_counts = (ad_sub.X > 0).sum(axis=0)
        # If ad_sub.X is sparse, convert result to flat array
        try:
            gene_counts = gene_counts.A1
        except AttributeError:
            gene_counts = gene_counts.ravel()
        genes_to_keep = ad_sub.var_names[gene_counts >= 10]
        ad_sub = ad_sub[:, genes_to_keep].copy()
        counts = ad_sub.obs[labels].value_counts()

        valid_targets = [
            t
            for t, c in counts.items()
            if t not in control_labels and c >= min_cells_per_group
        ]
        present_controls = [
            c for c in control_labels if counts.get(c, 0) >= min_cells_per_group
        ]
        if len(valid_targets) == 0 or len(present_controls) == 0:
            continue

        for ctrl in present_controls:
            per_ctrl_dfs: list[pd.DataFrame] = []
            tasks = (
                delayed(compute_wilcoxon_de_genes)(
                    ad_sub,
                    target=tgt,
                    control_label=ctrl,
                    labels=labels,
                    subclass_col=subclass_col,
                    method=method,
                )
                for tgt in valid_targets
            )
            run_in_parallel = n_jobs is not None and n_jobs != 1
            dfs = (
                Parallel(n_jobs=n_jobs, prefer="threads")(tasks)
                if run_in_parallel
                else [task() for task in tasks]
            )
            for df in dfs:
                if df is None or df.empty:
                    continue
                per_ctrl_dfs.append(df)
            if per_ctrl_dfs:
                results_by_control[ctrl].append(
                    pd.concat(per_ctrl_dfs, ignore_index=True)
                )

    # Combine per-control across all cell types
    combined_by_control: dict[str, pd.DataFrame] = {}
    for ctrl, df_list in results_by_control.items():
        if df_list:
            combined_by_control[ctrl] = pd.concat(df_list, ignore_index=True)

    # Return similar to pydeseq2 function
    if len(control_labels) == 1:
        return combined_by_control.get(control_labels[0], None)
    return combined_by_control


def compute_wilcoxon_de_genes(
    adata_sub: sc.AnnData,
    *,
    target: str,
    control_label: str,
    labels: str = "gene_target",
    subclass_col: str = "group_name",
    method: str = "wilcoxon",
    save_path: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """Compute DE for one target vs one control label (all cells, no matching).

    Returns a DataFrame with standard Scanpy DE output plus metadata columns,
    or None if not enough cells.
    """
    ctrl_cells = adata_sub.obs[adata_sub.obs[labels] == control_label]
    pert_cells = adata_sub.obs[adata_sub.obs[labels] == target]
    if len(ctrl_cells) == 0 or len(pert_cells) == 0:
        return None

    ad_pair = get_matched_adata(adata_sub, pert_cells, ctrl_cells)

    # Compute DE
    sc.tl.rank_genes_groups(
        ad_pair,
        groupby="deg_group",
        groups=["treatment"],
        reference="control",
        method=method,
        use_raw=False,
    )
    df = sc.get.rank_genes_groups_df(ad_pair, group="treatment")

    # Build results table
    df[subclass_col] = (
        adata_sub.obs[subclass_col].iloc[0] if subclass_col in adata_sub.obs else None
    )
    df[labels] = target
    df["control_label"] = control_label
    df["n_pert_matched"] = len(pert_cells)
    df["n_ctrl_matched"] = len(ctrl_cells)

    if save_path is not None:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)) or ".", exist_ok=True)
        df.to_csv(save_path, index=False)

    return df



def compute_gene_stats(
    df, adata, target_gene_col, cell_type_col, control_label: str = "Non_target", verbose: bool = True
):
    """
    Compute and add columns 'frac_expr_perturb', 'frac_expr_control', 'mean_expr_perturb',
    'mean_expr_control', 'mean_nnz_perturb', and 'mean_nnz_control' to df (fraction of cells
    expressing each gene, mean expression, and mean expression among expressing cells in
    perturbation vs control).

    Args:
        df: DataFrame of significant DEGs with 'names', target_gene_col, cell_type_col columns.
        adata: Annotated data object with .obs having columns target_gene_col, cell_type_col.
        target_gene_col: Name of the condition column (e.g., "gene_target").
        cell_type_col: Name of the group/celltype column (e.g., "group_name").
        control_label: Control condition label used for control stats (default "Non_target").
        verbose: Whether to print progress.

    Returns:
        df with new columns frac_expr_perturb, frac_expr_control, mean_expr_perturb,
        mean_expr_control, mean_nnz_perturb, mean_nnz_control.
    """
    import numpy as np
    import pandas as pd
    from collections import defaultdict
    from tqdm import tqdm

    def compute_fraction_expressed_vec(X):
        # X: 2D (n_cells, n_genes) or 1D (n_cells,) array for single gene
        if hasattr(X, "toarray"):
            X = X.toarray()
        return (X > 0).mean(axis=0) if X.size > 0 else np.nan

    def compute_mean_expr_vec(X):
        # X: 2D (n_cells, n_genes) or 1D (n_cells,) array for single gene
        if hasattr(X, "toarray"):
            X = X.toarray()
        return np.asarray(X).mean(axis=0) if X.size > 0 else np.nan

    def compute_mean_nnz_vec(X):
        # X: 2D (n_cells, n_genes) or 1D (n_cells,); per-gene mean of values > 0 (nan if none)
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X)
        if X.size == 0:
            return np.nan
        if X.ndim == 1:
            X = X.reshape(-1, 1)
        result = np.full(X.shape[1], np.nan, dtype=float)
        for j in range(X.shape[1]):
            col = X[:, j]
            pos = col > 0
            if pos.any():
                result[j] = col[pos].mean()
        return result[0] if result.size == 1 else result

    # Get all unique (target_gene, group) and (control_label, group) pairs in df
    perturb_conditions = df[[cell_type_col, target_gene_col]].drop_duplicates()
    control_conditions = df[[cell_type_col]].drop_duplicates().copy()
    control_conditions[target_gene_col] = control_label
    all_conditions = pd.concat([perturb_conditions, control_conditions], ignore_index=True)

    # (target_gene, group) -> set of gene names needed from df for this pair
    condition2genes = defaultdict(set)
    for idx, row in df.iterrows():
        condition2genes[(row[target_gene_col], row[cell_type_col])].add(row["names"])
        condition2genes[(control_label, row[cell_type_col])].add(row["names"])

    # For each unique (condition, cell type), compute gene-wise fraction expressed, mean expression, and mean nnz
    fractions = dict()  # (condition, cell type, gene) -> value
    means = dict()      # (condition, cell type, gene) -> value
    mean_nnz = dict()   # (condition, cell type, gene) -> value

    iterator = tqdm(all_conditions.iterrows(), total=len(all_conditions), desc="Computing fract. expressed & mean expr.") if verbose else all_conditions.iterrows()
    for _, row in iterator:
        target_gene = row[target_gene_col]
        group = row[cell_type_col]
        genes_this = [g for g in condition2genes[(target_gene, group)] if g in adata.var_names]
        if not genes_this:
            continue
        cond_mask = (adata.obs[target_gene_col] == target_gene) & (adata.obs[cell_type_col] == group)
        # If mask empty, skip
        if cond_mask.sum() == 0:
            for gene in genes_this:
                fractions[(target_gene, group, gene)] = np.nan
                means[(target_gene, group, gene)] = np.nan
                mean_nnz[(target_gene, group, gene)] = np.nan
            continue
        X = adata[cond_mask, genes_this].X  # shape (n_cells, len(genes_this))
        vec_frac = compute_fraction_expressed_vec(X)
        vec_mean = compute_mean_expr_vec(X)
        vec_mean_nnz = compute_mean_nnz_vec(X)
        # vec may be scalar if one gene, else array
        if np.isscalar(vec_frac) or (getattr(vec_frac, "shape", ()) == ()):
            vec_frac = [vec_frac]
        if np.isscalar(vec_mean) or (getattr(vec_mean, "shape", ()) == ()):
            vec_mean = [vec_mean]
        if np.isscalar(vec_mean_nnz) or (getattr(vec_mean_nnz, "shape", ()) == ()):
            vec_mean_nnz = [vec_mean_nnz]
        for i, gene in enumerate(genes_this):
            fractions[(target_gene, group, gene)] = vec_frac[i]
            means[(target_gene, group, gene)] = vec_mean[i]
            mean_nnz[(target_gene, group, gene)] = vec_mean_nnz[i]

    # Assign fractions, means, and mean_nnz back to df
    frac_expr_perturb = []
    frac_expr_control = []
    mean_expr_perturb = []
    mean_expr_control = []
    mean_nnz_perturb = []
    mean_nnz_control = []

    for idx, row in df.iterrows():
        group = row[cell_type_col]
        gene = row['names']
        pert_label = row[target_gene_col]
        frac_perturb = fractions.get((pert_label, group, gene), np.nan)
        frac_control = fractions.get((control_label, group, gene), np.nan)
        mean_perturb = means.get((pert_label, group, gene), np.nan)
        mean_control = means.get((control_label, group, gene), np.nan)
        nnz_perturb = mean_nnz.get((pert_label, group, gene), np.nan)
        nnz_control = mean_nnz.get((control_label, group, gene), np.nan)
        frac_expr_perturb.append(frac_perturb)
        frac_expr_control.append(frac_control)
        mean_expr_perturb.append(mean_perturb)
        mean_expr_control.append(mean_control)
        mean_nnz_perturb.append(nnz_perturb)
        mean_nnz_control.append(nnz_control)

    df = df.copy()
    df["frac_expr_perturb"] = frac_expr_perturb
    df["frac_expr_control"] = frac_expr_control
    df["mean_expr_perturb"] = mean_expr_perturb
    df["mean_expr_control"] = mean_expr_control
    df["mean_nnz_perturb"] = mean_nnz_perturb
    df["mean_nnz_control"] = mean_nnz_control
    return df