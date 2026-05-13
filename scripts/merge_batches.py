#!/usr/bin/env python
"""
Merge multiple batch h5ad files into a single concatenated file.

python merge_batches.py --data_path batch1.h5ad --data_path batch2.h5ad --output_dir /path/to/output --output_name merged.h5ad
"""

import argparse
import os

import scanpy as sc

def merge_batches(
    data_paths: list,
    output_dir: str,
    output_name: str,
    write_zarr: bool = False,
    write_h5ad: bool = True,
) -> None:
    """
    Load and merge multiple h5ad files into a single concatenated file.
    
    Parameters
    ----------
    data_paths : list
        List of paths to h5ad files to merge.
    output_dir : str
        Directory where the merged file will be saved.
    output_name : str
        Name of the output file.
    """
    if not data_paths:
        raise ValueError("At least one --data_path argument is required")
    
    print(f"Found {len(data_paths)} h5ad files to merge:")
    for path in data_paths:
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")
        print(f"  - {os.path.basename(path)}")
    
    # Load all h5ad files
    adatas = []
    for data_path in data_paths:
        adata = sc.read_h5ad(data_path)
        print(
            f"  Loaded {os.path.basename(data_path)}: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes"
        )
        adatas.append(adata)
    
    # Concatenate all batches
    print("\nConcatenating batches...")
    merged_adata = sc.concat(adatas, join="outer", merge="unique", uns_merge="unique")
    
    # Make observation names unique
    merged_adata.obs_names_make_unique()
    
    print(
        f"Merged dataset: {merged_adata.shape[0]:,} cells × {merged_adata.shape[1]:,} genes"
    )
    
    # Save merged adata
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_name)
    
    if write_zarr:
        merged_adata.write_zarr(output_path)
    if write_h5ad:
        merged_adata.write_h5ad(output_path)
    
    print(f"\nSaved merged adata to: {output_path}")
    print(
        f"Final shape: {merged_adata.shape[0]:,} cells × {merged_adata.shape[1]:,} genes"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Merge multiple batch h5ad files into a single concatenated file."
    )
    parser.add_argument(
        "--data_path",
        type=str,
        action="append",
        required=True,
        help="Path to an h5ad file to merge (can be specified multiple times)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Directory where the merged adata will be saved",
    )
    parser.add_argument(
        "--output_name",
        type=str,
        required=True,
        help="Name of the output file",
    )
    parser.add_argument(
        "--write_zarr",
        type=bool,
        default=False,
        help="Write the merged adata to a zarr file",
    )
    parser.add_argument(
        "--write_h5ad",
        type=bool,
        default=True,
        help="Write the merged adata to a h5ad file",
    )

    args = parser.parse_args()

    merge_batches(
        data_paths=args.data_path,
        output_dir=args.output_dir,
        output_name=args.output_name,
        write_zarr=args.write_zarr,
        write_h5ad=args.write_h5ad,
    )


if __name__ == "__main__":
    main()