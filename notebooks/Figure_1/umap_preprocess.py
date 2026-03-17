#!/usr/bin/env python
"""
umap_preprocess.py

Assembles raw zarr from downloaded h5ad shards, creates preprocess_config.yaml,
and runs preprocess.py. Assumes h5ads have been downloaded with download.py
(output-dir/h5ads/).

Usage:
    python umap_preprocess.py \\
        --output-dir /path/to/output \\
        --use-genes-path /path/to/use_genes.txt \\
        [--chunk-size 20000] \\
        [--skip-zarr] \\
        [--skip-preprocess]
"""

from __future__ import annotations

import argparse
import sys
import subprocess
import textwrap
import time
import yaml
from pathlib import Path

import numpy as np

PREPROCESS_SCRIPT = Path(__file__).parent / "preprocess.py"


# ---------------------------------------------------------------------------
# Assemble raw zarr from h5ad shards
# ---------------------------------------------------------------------------
def build_raw_zarr(h5ad_paths: list[Path], out_zarr_path: Path, chunk_size: int = 20_000) -> None:
    """
    Streams per-batch h5ad shards into a single zarr store without loading all X
    into memory at once.

    Pipeline:
      1. Pass 1: read obs/var metadata from each shard to compute total
         shape and build concatenated obs; write skeleton AnnData (no X) to zarr
         to persist obs/var/uns; pre-allocate the full X array.
      2. Pass 2: stream each shard's X into the pre-allocated zarr rows.
    """
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    import zarr

    if out_zarr_path.exists():
        print(f"\n[build_raw_zarr] zarr already exists at {out_zarr_path}, skipping.")
        return

    print(f"\n[build_raw_zarr] Pass 1 — reading metadata from {len(h5ad_paths)} shards ...")
    t0 = time.time()

    obs_frames = []
    n_vars = None
    var_df = None
    uns = {}
    for p in h5ad_paths:
        a = ad.read_h5ad(str(p), backed="r")
        obs_frames.append(a.obs.copy())
        if n_vars is None:
            n_vars = a.n_vars
            var_df = a.var.copy()
            uns = dict(a.uns)
        a.file.close()

    obs_all = pd.concat(obs_frames, join="outer")
    n_obs = len(obs_all)
    print(f"  Total shape: {n_obs:,} × {n_vars:,}  ({(time.time() - t0) / 60:.1f} min)")

    # Write obs/var/uns skeleton (no X) to zarr
    skeleton = ad.AnnData(obs=obs_all, var=var_df, uns=uns)
    skeleton.X = None
    print(f"  Writing metadata to {out_zarr_path} ...")
    skeleton.write_zarr(str(out_zarr_path))
    del skeleton, obs_all, obs_frames

    # Cap chunk_size so each chunk stays under the 2 GB blosc limit
    max_chunk_rows = int((2**31 - 1) / (n_vars * 4))
    if chunk_size > max_chunk_rows:
        print(f"  chunk_size {chunk_size} → capped to {max_chunk_rows} (blosc 2 GB limit)")
        chunk_size = max_chunk_rows

    # use_consolidated=False avoids stale .zmetadata that doesn't include X yet.
    Z = zarr.open_group(str(out_zarr_path), mode="a", zarr_format=2, use_consolidated=False)
    Z_X = Z.create_array(
        "X",
        shape=(n_obs, n_vars),
        chunks=(chunk_size, n_vars),
        dtype="float32",
        overwrite=True,
    )

    # Pass 2: stream X from each shard into the pre-allocated array
    print(f"\n[build_raw_zarr] Pass 2 — streaming X from each shard ...")
    t1 = time.time()
    row = 0
    for p in h5ad_paths:
        a = ad.read_h5ad(str(p))  # full load — one shard at a time
        X = a.X
        if sp.issparse(X):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float32)
        Z_X[row : row + a.n_obs] = X
        row += a.n_obs
        elapsed = time.time() - t1
        print(f"  {p.name}  rows …{row:,}  ({elapsed / 60:.1f} min)", flush=True)
        del a, X

    Z_X.attrs["encoding-type"] = "array"
    Z_X.attrs["encoding-version"] = "0.2.0"
    Z_X.attrs["shape"] = [n_obs, n_vars]

    zarr.consolidate_metadata(str(out_zarr_path))
    print(f"  Done in {(time.time() - t0) / 60:.1f} min total", flush=True)


def write_preprocess_config(
    config_path: Path,
    raw_data_path: Path,
    out_zarr_path: Path,
    use_genes_path: Path,
) -> None:
    cfg = {
        "raw_data_path": str(raw_data_path),
        "out_zarr_path": str(out_zarr_path),
        "use_genes_path": str(use_genes_path),
        "sparse_chunk_size": 20000,
        "filter": True,
        "normalize": True,
    }
    with open(config_path, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
    print(f"\n[config] Written to {config_path}")
    print(textwrap.indent(yaml.dump(cfg, default_flow_style=False, sort_keys=False), "  "))


def run_preprocess(config_path: Path) -> None:
    print(f"\n[preprocess] Running: python {PREPROCESS_SCRIPT} {config_path}")
    result = subprocess.run(
        [sys.executable, str(PREPROCESS_SCRIPT), str(config_path)],
        check=False,
    )
    if result.returncode != 0:
        sys.exit(f"preprocess.py exited with code {result.returncode}")
    print("\n[preprocess] Done.")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Assemble raw zarr from h5ads, create preprocess config, and run preprocess.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
            Example:
              python umap_preprocess.py \\
                --output-dir /data/wholebrain \\
                --use-genes-path /data/use_genes.txt
        """),
    )
    p.add_argument("--output-dir", required=True,
                   help="Root directory containing h5ads/ (from download.py) and where raw zarr will be written")
    p.add_argument("--use-genes-path", required=True,
                   help="Path to use_genes.txt (marker gene list for normalization)")
    p.add_argument("--chunk-size", type=int, default=20_000,
                   help="Row chunk size for zarr storage (default: 20000)")
    p.add_argument("--skip-zarr", action="store_true",
                   help="Skip zarr assembly (assumes raw_wholebrain.zarr already exists)")
    p.add_argument("--skip-preprocess", action="store_true",
                   help="Only write preprocess_config.yaml, do not run preprocess.py")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    output_dir = Path(args.output_dir).resolve()
    h5ad_dir = output_dir / "h5ads"
    raw_zarr_path = output_dir / "raw_wholebrain.zarr"
    filtered_zarr_path = output_dir / "filtered_lognorm.zarr"
    config_path = output_dir / "preprocess_config.yaml"
    use_genes_path = Path(args.use_genes_path).resolve()

    # Step 1: Assemble raw zarr from h5ads (unless skipped)
    print("=" * 60)
    print("Step 1/3  Assemble raw zarr from h5ad shards")
    print("=" * 60)
    if args.skip_zarr:
        print("  [skipped via --skip-zarr]")
        if not raw_zarr_path.exists():
            sys.exit(f"Raw zarr not found at {raw_zarr_path}. Run without --skip-zarr or run download.py first.")
    else:
        h5ad_paths = sorted(h5ad_dir.glob("*.h5ad"))
        if not h5ad_paths:
            sys.exit(f"No h5ad files found in {h5ad_dir}. Run download.py first.")
        build_raw_zarr(h5ad_paths, raw_zarr_path, chunk_size=args.chunk_size)

    # Step 2: Write preprocess config
    print("\n" + "=" * 60)
    print("Step 2/3  Write preprocess config")
    print("=" * 60)
    write_preprocess_config(config_path, raw_zarr_path, filtered_zarr_path, use_genes_path)

    # Step 3: Run preprocess (unless skipped)
    if args.skip_preprocess:
        print(f"\n[skip-preprocess] When ready, run:\n  python {PREPROCESS_SCRIPT} {config_path}")
    else:
        print("\n" + "=" * 60)
        print("Step 3/3  Run preprocess.py")
        print("=" * 60)
        run_preprocess(config_path)


if __name__ == "__main__":
    main()
