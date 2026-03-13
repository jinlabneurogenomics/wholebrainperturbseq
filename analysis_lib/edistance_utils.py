from typing import Dict
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import pertpy as pt
from tqdm import tqdm

def compute_edistance_per_celltype(
    adata,
    *,
    perturbation_col: str = "gene_target",
    cell_type_col: str = "group_name",
    n_pcs: int = 50,
    n_top_genes: int = 2000,
    out_dir: str = "/data/perturbai/WB8588_screen/analysis/edistance_results_pertpy",
    n_jobs: int = -1, 
) -> Dict[str, pd.DataFrame]:
    os.makedirs(out_dir, exist_ok=True)
    cell_types = (
        adata.obs[cell_type_col]
        .value_counts()
        .sort_values(ascending=True)
        .index
        .tolist()
    )
    distance = pt.tl.Distance("edistance", obsm_key="X_pca")
    results: Dict[str, pd.DataFrame] = {}
    for ct in tqdm(cell_types):
        print(f"Computing E-distance for {ct}")
        df = distance.pairwise(adata, groupby=perturbation_col, n_jobs=n_jobs)
        results[ct] = df
        df.to_csv(f"{out_dir}/{ct.split(' ')[0]}.csv")

    return results


def plot_and_save_clustermaps(results_dict,  out_dir="edistance_heatmaps", save=False):
    import os
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    plt.close('all')
    for cell_type, sub_df in results_dict.items():
        if len(sub_df) == 1:
            continue

        # 3. Plot clustermap (show all x/y labels)
        n_rows, n_cols = sub_df.shape
        # Scale figure so labels fit; use small font when many labels
        max_labels = max(n_rows, n_cols)
        fontsize = max(8, 8 - max_labels // 20)
        fig_h = max(10, min(30, 8 + n_rows * 0.15))
        fig_w = max(10, min(30, 8 + n_cols * 0.15))
        sub_df[sub_df < 0] = 0

        # --- Here, capture the normed matrix for colorbar labeling
        normed_sub_df = np.log1p(sub_df)

        g = sns.clustermap(
            normed_sub_df,
            figsize=(fig_w, fig_h),
            method="weighted",
            xticklabels=True,
            yticklabels=True,
            cmap="crest"
        )

        # Add a title to the clustermap
        print(f"Clustermap: {cell_type}")

        # Ensure all tick labels are visible: rotate x, set font size
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xticklabels(),
            rotation=90,
            ha="right",
            fontsize=fontsize,
        )
        g.ax_heatmap.set_yticklabels(
            g.ax_heatmap.get_yticklabels(),
            rotation=0,
            fontsize=fontsize,
        )

        # 5. Set color legend label for log1p(E-distance)
        colorbar = g.cax
        colorbar.set_ylabel('log1p(E-distance)', fontsize=fontsize+2)

        if save:
            # Save PDF
            pdf_path = os.path.join(out_dir, f"{cell_type}.pdf")
            g.savefig(pdf_path, dpi=300, format="pdf")
            plt.close('all')
        else:
            plt.show()
