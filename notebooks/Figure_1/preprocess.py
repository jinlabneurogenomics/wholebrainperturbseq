from __future__ import annotations

import gc
import sys
import os
import time
import numpy as np
import yaml
import dask
from dask_cuda import LocalCUDACluster
from dask.distributed import Client
import dask.array as da
import rmm
import cupy as cp

from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(pool_allocator=True, initial_pool_size=2**33)  # 8 GB seed; pool grows dynamically up to free VRAM
cp.cuda.set_allocator(rmm_cupy_allocator)

import anndata as ad
import zarr
import rapids_singlecell as rsc
import scipy.sparse as sp
from packaging.version import parse as parse_version

OUT_DIR = "/data/WholeBrainPerturbSeq/for_umap"
# Handle version-specific import for anndata experimental
if parse_version(ad.__version__) < parse_version("0.12.0rc1"):
    print("anndata version is less than 0.12.0rc1")
    from anndata.experimental import read_elem_as_dask as read_dask
else:
    print("anndata version is greater than 0.12.0rc1")
    from anndata.experimental import read_elem_lazy as read_dask

# ---------------------------------------------------------------------------
# Constants / config
# ---------------------------------------------------------------------------
def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)

if __name__ == "__main__":
    config = load_config(sys.argv[1])  # python script.py config.yaml
    print(sys.argv[1])

    SPARSE_CHUNK_SIZE = config["sparse_chunk_size"]

    cluster = LocalCUDACluster(
        CUDA_VISIBLE_DEVICES="1",
        threads_per_worker=2,
        protocol="tcp",
        rmm_allocator_external_lib_list="cupy",
        memory_limit="auto",
    )
    client = Client(cluster)
    print(client)

    data_pth = config["raw_data_path"]  # 7.7M cells minus leiden clusters

    f = zarr.open(data_pth, mode="r", use_consolidated=False)
    X = f["X"]
    shape = X.attrs.get("shape", None) or list(X.shape)

    adata = ad.AnnData(
        X = read_dask(X, (SPARSE_CHUNK_SIZE, shape[1])),
        obs = ad.io.read_elem(f["obs"]),
        var = ad.io.read_elem(f["var"]),
        uns = ad.io.read_elem(f["uns"]),
    )

    # apply filter
    if config["filter"]:
        adata = adata[(adata.obs["scDblFinder.class"]=="singlet") 
                & (adata.obs["num_genes"]>=2000)
                & (adata.obs["log_ambient_mse_norm"]>0.09)].copy()
        is_gex = ~adata.var_names.str.contains("_")
        adata = adata[:,is_gex].copy()

    # reorganize chunks for downstream (like old code: applied before GPU pipeline)
    adata = adata[
        np.lexsort((
            adata.obs["predicted_subclass"].values,
            adata.obs["predicted_class"].values,
            adata.obs["neighborhood"].values,
        ))
    ].copy()

    # load marker genes
    use_genes_path = config["use_genes_path"]
    with open(use_genes_path, 'r') as f:
        use_genes = [line.strip() for line in f if line.strip()]

    out_zarr_path = config.get(
        "out_zarr_path",
        os.path.join(OUT_DIR, "260306_subgenome_filter_lognorm_6M.zarr"),
    )

    # normalize
    if config["normalize"] and "log1p" not in adata.uns:
        print("Log normalizing and subsetting to markers")
        adata.var["marker_genes"] = adata.var_names.isin(use_genes)
        rsc.get.anndata_to_GPU(adata)
        rsc.pp.normalize_total(adata, target_sum=1e4)
        rsc.pp.log1p(adata)
        adata = adata[:, adata.var.marker_genes].copy()
        rsc.get.anndata_to_CPU(adata)
        # Fix 1: anndata_to_CPU now returns a lazy dask array in the current rapids_singlecell
        # version. Materialize with synchronous scheduler (one zarr-chunk at a time on GPU,
        # ~1.5 GB peak) so the dask shuffle from the sort above executes without OOM.
        if isinstance(adata.X, da.Array):
            print("  materializing dask X (synchronous scheduler) ...", flush=True)
            with dask.config.set(scheduler="synchronous"):
                adata.X = adata.X.compute()
            if hasattr(adata.X, "get"):
                adata.X = adata.X.get()
        gc.collect()
    elif "log1p" in adata.uns:
        print("Data already normalized, skipping")

    # All GPU work done — shut down the cluster before writing.
    client.close()
    cluster.close()
    gc.collect()

    print("writing zarr")
    if not os.path.exists(out_zarr_path):
        print(f"Writing filtered + sorted adata to {out_zarr_path}")
        t0 = time.time()

        # Fix 2: adata.write_zarr SEGFAULTs with an active CUDA context in the current
        # anndata/zarr versions. Write obs/var/uns via write_elem instead.
        X_data = adata.X
        adata.X = None
        Z = zarr.open_group(out_zarr_path, mode="w", zarr_format=2, use_consolidated=False)
        Z.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
        ad.io.write_elem(Z, "obs", adata.obs)
        ad.io.write_elem(Z, "var", adata.var)
        ad.io.write_elem(Z, "uns", dict(adata.uns))
        ad.io.write_elem(Z, "obsm", dict(adata.obsm))
        adata.X = X_data

        # Blosc codec limit: chunks must be < 2 GiB (2^31 - 1 bytes)
        MAX_BLOSC_BYTES = 2**31 - 1
        WRITE_CHUNK = min(SPARSE_CHUNK_SIZE, MAX_BLOSC_BYTES // (adata.n_vars * 4))
        Z_X = Z.require_array(
            "X",
            shape=(adata.n_obs, adata.n_vars),
            chunks=(WRITE_CHUNK, adata.n_vars),
            dtype="float32",
            overwrite=True,
        )
        for start in range(0, adata.n_obs, WRITE_CHUNK):
            end = min(start + WRITE_CHUNK, adata.n_obs)
            block = adata.X[start:end]
            if sp.issparse(block):
                block = block.toarray()
            Z_X[start:end] = np.asarray(block, dtype=np.float32)
        Z_X.attrs["encoding-type"] = "array"
        Z_X.attrs["encoding-version"] = "0.2.0"
        zarr.consolidate_metadata(out_zarr_path)
        elapsed = time.time() - t0
        print(f"Done writing zarr in {elapsed/60:.2f}min ({elapsed:.1f}s)", flush=True)
    else:
        print(f"Found existing filtered zarr at {out_zarr_path}, skipping write")
