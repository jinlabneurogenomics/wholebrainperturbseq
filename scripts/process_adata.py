#!/usr/bin/env python
"""
Process raw h5ad files: concatenate, separate guide/gene data, extract top guides, annotate target gene.

python process_adata.py  /data/WholeBrainPerturbSeq/h5ads --min_guide_UMI 3
"""

import argparse
import glob
import os

import numpy as np
import scanpy as sc
from scipy import sparse
import sys
sys.path.append("/workspace/wholebrainperturbseq")
from analysis_lib.utils import deduplicate_gene_symbols

def process_adata(
    data_path: str,
    output_dir: str,
    output_name: str,
    guide_output_name: str = "WB8588_screen_crispr_guides.h5ad",
    min_guide_UMI: int = 3,
) -> None:
    """
    Load and process h5ad files from data_path and save processed adata to output_dir.
    Saves two outputs:
      - Gene expression anndata with guide annotations in obs.
      - CRISPR guide UMI anndata with the same obs metadata.
    data_path can be a single .h5ad file or a directory containing *.h5ad files.
    """
    # Load all h5ad files
    if os.path.isfile(data_path):
        print(f"Processing single h5ad file: {data_path}")
        adata_files = [data_path]
    elif "*.h5ad" in data_path:
        print(f"Processing multiple h5ad files: {data_path}")
        adata_files = sorted(glob.glob(data_path))
    else:
        print(f"Processing directory: {data_path}")
        adata_files = sorted(glob.glob(os.path.join(data_path, "*.h5ad")))
        if not adata_files:
            raise FileNotFoundError(f"No .h5ad files found in {data_path}")
    print(f"Found {len(adata_files)} h5ad files to concatenate:")
    for f in adata_files:
        print(f"  - {os.path.basename(f)}")

    tables = []
    import re

    batch_pattern = re.compile(r"WB8588_([a-zA-Z0-9]+)_[a-zA-Z0-9]+\.h5ad$")

    for adata_file in adata_files:
        filename = os.path.basename(adata_file)
        match = batch_pattern.match(filename)
        if not match:
            raise ValueError(
                f"Filename '{filename}' does not match expected format WB8588_{{batch}}_{{channel}}.h5ad"
            )
        adata_temp = sc.read_h5ad(adata_file)
        print(
            f"  Loaded {os.path.basename(adata_file)}: {adata_temp.shape[0]:,} cells × {adata_temp.shape[1]:,} genes"
        )
        tables.append(adata_temp)

    adata = sc.concat(tables, join="outer", merge="unique", uns_merge="unique")

    # Make observation names unique
    adata.obs_names_make_unique()
    # Add gene symbols and IDs
    adata.var["gene_symbols"] = adata.var.index.tolist()
     # Identify mitochondrial genes using robust symbol detection
    adata = deduplicate_gene_symbols(
        adata, gene_symbol_col="gene_symbols", gene_id_col="gene_ids"
    )
    # Identify mitochondrial genes using robust symbol detection
    symbols = adata.var["gene_symbols"].astype(str).str.lower()
    adata.var["mt"] = symbols.str.startswith("mt-").astype(bool)
    # Compute QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata.obs.rename(columns={"predicted_subclass": "subclass_name"}, inplace=True)
    adata.obs.rename(columns={"predicted_group": "group_name"}, inplace=True)
    adata.obs.rename(columns={"predicted_class": "class_name"}, inplace=True)
    adata.obs.rename(columns={"predicted_supertype": "supertype_name"}, inplace=True)
    adata.obs.rename(columns={"predicted_cluster": "cluster_name"}, inplace=True)
    adata.obs.rename(columns={"guide_call": "target_guide"}, inplace=True)
    adata.obs.rename(columns={"source": "mouse"}, inplace=True)
    # Extract CRISPR Guide Capture features and save in separate layer
    # Separate guide and gene counts
    guide_mask = adata.var["feature_types"] == "CRISPR Guide Capture"
    gene_mask = adata.var["feature_types"] == "Gene Expression"
    adata_genes = adata[:, gene_mask].copy()
    adata_guides = adata[:, guide_mask].copy()
    print(f"\nGuide data shape: {adata_guides.shape}")
    print(f"Gene data shape: {adata_genes.shape}")

    # Extract top 2 guides per cell from guide count matrix
    guide_matrix = adata_guides.X
    guide_names = np.array(adata_guides.var.index)
    if sparse.issparse(guide_matrix):
        guide_matrix = guide_matrix.toarray()
    n_cells = guide_matrix.shape[0]
    top_2_indices = np.argsort(-guide_matrix, axis=1)[:, :2]

    # Get UMI counts for top 2 guides
    top_guide_1_umi = guide_matrix[np.arange(n_cells), top_2_indices[:, 0]]
    top_guide_2_umi = guide_matrix[np.arange(n_cells), top_2_indices[:, 1]]

    # Check that the values are the same as the ones reporte in the original obs
    assert np.all(adata_genes.obs["guide_umi_second"] == top_guide_2_umi)
    assert np.all(adata_genes.obs["guide_umi_top"] == top_guide_1_umi)

    # Get guide labels for top 2
    top_guide_1_names = guide_names[top_2_indices[:, 0]]
    # Apply Negative label if top guide UMI is less than min_guide_UMI
    top_guide_1_umi = guide_matrix[np.arange(n_cells), top_2_indices[:, 0]]
    top_guide_1_names = np.where(
        top_guide_1_umi < min_guide_UMI, "unassigned", top_guide_1_names
    )
    top_guide_2_names = guide_names[top_2_indices[:, 1]]
    top_guide_2_names = np.where(
        top_guide_2_umi < min_guide_UMI, "unassigned", top_guide_2_names
    )

    # Store in obs
    adata_genes.obs["top_guide_1"] = top_guide_1_names
    adata_genes.obs["top_guide_2"] = top_guide_2_names

    # Save processed adata
    os.makedirs(output_dir, exist_ok=True)
    gene_output_path = os.path.join(output_dir, output_name)
    guide_output_path = os.path.join(output_dir, guide_output_name)

    adata_genes.write_h5ad(gene_output_path)
    adata_guides.write_h5ad(guide_output_path)

    print(f"\nSaved gene expression adata to: {gene_output_path}")
    print(f"Saved CRISPR guide adata to: {guide_output_path}")
    print(
        f"Gene expression shape: {adata_genes.shape[0]:,} cells × {adata_genes.shape[1]:,} genes"
    )
    print(
        f"CRISPR guide shape: {adata_guides.shape[0]:,} cells × {adata_guides.shape[1]:,} guides"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Process raw h5ad files: concatenate, separate guide/gene data, extract top guides."
    )
    parser.add_argument(
        "--data_path",
        type=str,
        required=True,
        help="Path to a single .h5ad file or a directory containing raw h5ad files (*.h5ad)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Path to directory where processed adata will be saved",
    )
    parser.add_argument(
        "--output_name",
        type=str,
        required=False,
        help="Name of the output file",
    )
    parser.add_argument(
        "--guide_output_name",
        type=str,
        required=False,
        help="Name of the CRISPR guide output file",
    )
    parser.add_argument(
        "--min_guide_UMI",
        type=int,
        required=False,
        default=3,
        help="Minimum guide UMI count required to assign a guide; otherwise label as unassigned",
    )

    args = parser.parse_args()

    process_adata(
        data_path=args.data_path,
        output_dir=args.output_dir,
        output_name=args.output_name,
        guide_output_name=args.guide_output_name,
        min_guide_UMI=args.min_guide_UMI,
    )


if __name__ == "__main__":
    main()
