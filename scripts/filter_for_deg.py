#!/usr/bin/env python
"""
Filter processed WholeBrainPerturbSeq data.

This script:
1. Loads data processed by process_adata.py using get_processed_data_for_deg
2. Applies QC filters (mt%, total counts, subclass, remove cells with no guides assigned)
3. Saves the filtered data to output directory

Usage example:
    python filter_for_deg.py \
        --data_path /data/WholeBrainPerturbSeq/processed/WB8588_screen_gex.h5ad \
        --output_dir /data/WholeBrainPerturbSeq/filtered \
        --output_name WB8588_screen_gex_filtered.h5ad \
        --min_cells_per_gene 10 \
        --min_cells_per_subclass 50 \
        --remove_unassigned \
        --remove_multiple_guides \
        --remove_doublets
"""

from __future__ import annotations

import argparse
import os
import sys

from tqdm import tqdm

sys.path.append("/workspace/wholebrainperturbseq")

from analysis_lib.utils import get_processed_data_for_deg


def main():
    parser = argparse.ArgumentParser(
        description="Filter processed WholeBrainPerturbSeq data."
    )

    # Input/Output arguments
    parser.add_argument(
        "--data_path",
        type=str,
        required=True,
        help="Path to processed h5ad file (output of process_adata.py)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory to save filtered data",
    )
    parser.add_argument(
        "--output_name",
        type=str,
        default="filtered.h5ad",
        help="Name of output file",
    )

    # Filter arguments (passed to get_processed_data_for_deg)
    parser.add_argument(
        "--min_cells_per_gene",
        type=int,
        default=10,
        help="Keep genes expressed in at least this many cells (default: 10)",
    )
    parser.add_argument(
        "--min_genes_per_cell",
        type=int,
        default=0,
        help="Keep cells with at least this many genes detected (default: 0)",
    )
    parser.add_argument(
        "--max_pct_mt",
        type=float,
        default=None,
        help="Filter cells above this mitochondrial percentage (default: None)",
    )
    parser.add_argument(
        "--min_total_counts",
        type=float,
        default=0,
        help="Filter cells below this total counts (default: 2000)",
    )
    parser.add_argument(
        "--min_cells_per_subclass",
        type=int,
        default=0,
        help="Keep subclasses with at least this many cells (default: 50)",
    )
    parser.add_argument(
        "--remove_unassigned",
        action="store_true",
        help="Remove cells with unassigned target gene",
    )
    parser.add_argument(
        "--unassigned_label",
        type=str,
        default="Negative",
        help="Label used for unassigned cells (default: 'Negative')",
    )
    parser.add_argument(
        "--subclass_col",
        type=str,
        default="subclass_name",
        help="Column name for subclass/cell type (default: 'subclass_name')",
    )
    parser.add_argument(
        "--target_gene_col",
        type=str,
        default="gene_target",
        help="Column name for target gene (default: 'gene_target')",
    )
    parser.add_argument(
        "--log_ambient_mse_norm_min",
        type=float,
        default=None,
        help="Keep cells with log_ambient_mse_norm above this value (default: None)",
    )
    parser.add_argument(
        "--remove_multiple_guides",
        action="store_true",
        help="Remove cells with multiple guides (default: False)",
    )
    parser.add_argument(
        "--remove_doublets",
        action="store_true",
        help="Remove doublets (default: False)",
    )

    parser.add_argument(
        "--write_zarr",
        action="store_true",
        help="Write data to zarr file (default: False)",
    )
    parser.add_argument(
        "--write_h5ad",
        action="store_true",
        help="Write data to h5ad file (default: False)",
    )

    args = parser.parse_args()

    print("Filter Pipeline")
    print(f"\nInput: {args.data_path}")
    print(f"Output: {os.path.join(args.output_dir, args.output_name)}")

    # Step 1: Load and filter data using get_processed_data_for_deg
    print("\n" + "-" * 40)
    print("Step 1: Loading and filtering data")
    print("-" * 40)

    with tqdm(total=1, desc="Loading and filtering", unit="step") as pbar:
        adata = get_processed_data_for_deg(
            args.data_path,
            min_cells_per_gene=args.min_cells_per_gene,
            min_genes_per_cell=args.min_genes_per_cell,
            min_cells_per_subclass=args.min_cells_per_subclass,
            remove_unassigned=args.remove_unassigned,
            unassigned_label=args.unassigned_label,
            subclass_col=args.subclass_col,
            target_gene_col=args.target_gene_col,
            min_total_counts=args.min_total_counts,
            max_total_counts=None,
            max_pct_mt=args.max_pct_mt,
            log_ambient_mse_norm_min=args.log_ambient_mse_norm_min,
            remove_multiple_guides=args.remove_multiple_guides,
            verbose=True,
        )
        pbar.update(1)

    print(f"\nAfter filtering: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    # Step 2: Save output
    print("Step 2: Saving output")

    os.makedirs(args.output_dir, exist_ok=True)
    output_path = os.path.join(args.output_dir, args.output_name)

    if args.write_zarr:
        with tqdm(total=1, desc="Writing zarr file", unit="file") as pbar:
            adata.write_zarr(output_path)
            pbar.update(1)
    
    if args.write_h5ad:
        with tqdm(total=1, desc="Writing h5ad file", unit="file") as pbar:
            adata.write_h5ad(output_path)
            pbar.update(1)

    print(f"\nSaved to: {output_path}")
    print(f"Final shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Print summary
    print("Summary")
    print(f"Target gene distribution (top 20):")
    print(adata.obs[args.target_gene_col].value_counts().head(20).to_string())


if __name__ == "__main__":
    main()
