from __future__ import annotations

import os
import time
import shutil
from typing import Optional
import importlib
import psutil
import yaml
import sys
import gc
from datetime import datetime
sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)

import numpy as np
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import zarr

import matplotlib as mpl
import matplotlib.pyplot as plt 
mpl.rcParams['pdf.fonttype'] = 42
sc.settings.n_jobs = 16
sc.set_figure_params(figsize=(6,6), vector_friendly = True)

import dask.array as da
from dask_cuda import LocalCUDACluster
from dask.distributed import Client, as_completed, get_client, wait

from packaging.version import parse as parse_version

import rmm
import cupy as cp
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(pool_allocator=True, initial_pool_size=2**33)  # 8 GB initial pool
cp.cuda.set_allocator(rmm_cupy_allocator)
import rapids_singlecell as rsc
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed as thread_as_completed
import threading

import logging

if parse_version(importlib.metadata.version('anndata')) < parse_version("0.12.0rc1"):
    from anndata.experimental import read_elem_as_dask as read_dask
else:
    from anndata.experimental import read_elem_lazy as read_dask

class TaggedLogger:
    def __init__(self, logger, tag):
        self.logger = logger
        self.tag    = tag

    def info(self, msg, *args, **kwargs):
        self.logger.info(f"[{self.tag}] {msg}", *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        self.logger.warning(f"[{self.tag}] {msg}", *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self.logger.error(f"[{self.tag}] {msg}", *args, **kwargs)

    def debug(self, msg, *args, **kwargs):
        self.logger.debug(f"[{self.tag}] {msg}", *args, **kwargs)

# config / constants
def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)

# Module-level stubs — populated from config in __main__
config: dict = {}
SPARSE_CHUNK_SIZE: int = 20_000
PARENT_OBS_KEY: dict = {}
CACHE_KEY_HIER3: str = "impute_knn_cache_hier3"
HIER_LEVELS: list = []
SEED: int = 0
CELL_BATCH_SIZE: int = 400_000
OUT_DIR: str = "."
GPU_KNN_BATCH: int = 150_000
VRAM_LIMIT_GB: float = 18.0  # GPU memory threshold; groups above this fall back to CPU KNN

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("impute")

# helpers
def _level_to_out_layer(level):
    return f"imputed_{level}" if level is not None else "imputed_global"


def _ensure_cache_hier3(adata, cache_key=CACHE_KEY_HIER3):
    if cache_key not in adata.uns or not isinstance(adata.uns[cache_key], dict):
        adata.uns[cache_key] = {}
    return adata.uns[cache_key]


def _get_source_matrix_for_cells(adata, cell_indices, source_layer=None):
    """Materialise only the requested rows; works for dask, sparse, or dense X."""
    X = adata.X if source_layer is None else adata.layers[source_layer]
    try:
        import dask.array as da
        if isinstance(X, da.Array):
            block = X[cell_indices].compute(scheduler="threads")
            return block.toarray().astype(np.float32) if sp.issparse(block) \
                   else np.asarray(block, dtype=np.float32)
    except ImportError:
        pass
    if sp.issparse(X):
        return X[cell_indices].toarray().astype(np.float32)
    return np.asarray(X[cell_indices], dtype=np.float32)

def _validate_hier3_inputs(adata, levels):
    missing = [k for k in levels if k is not None and k not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing required obs keys: {missing}")
    if adata.n_obs < 2:
        raise ValueError("Need at least 2 cells.")

# marker gene helpers
def _get_marker_genes_for_group(adata, rgg_key, group, n_top=20, n_knn=10):
    try:
        df  = sc.get.rank_genes_groups_df(adata, group=group, key=rgg_key)
        pos = df[df["logfoldchanges"] > 0].sort_values("scores", ascending=False)
        imp = pos.head(n_top)["names"].tolist()
        return imp[:n_knn], imp
    except Exception as e:
        print(f"  [marker genes] rgg_key={rgg_key} group={group}: {e}. All-gene fallback.")
        return None, None


def _resolve_parent_key_for_group(adata, group_key, grp):
    parent_col = PARENT_OBS_KEY[group_key]
    mask       = adata.obs[group_key].astype(str) == str(grp)
    parents    = adata.obs.loc[mask, parent_col].unique()
    if len(parents) == 0:
        raise ValueError(f"No cells for group_key={group_key}, grp={grp}")
    if len(parents) > 1:
        raise ValueError(f"grp='{grp}' maps to multiple parents: {parents}")
    return str(parents[0])


def _build_gene_dicts_for_level(adata, group_key, n_top=20, n_knn=10):
    log    = logging.getLogger("gene dicts")
    groups = adata.obs[group_key].astype(str).unique()
    knn_g, imp_g, rgg_map = {}, {}, {}
    all_imp_genes = set()
    all_knn_genes = set()

    for grp in groups:
        try:
            pv      = _resolve_parent_key_for_group(adata, group_key, grp)
            rgg_key = f"{pv}_rank_genes_groups"
        except Exception as e:
            log.info(f"Skipping grp={grp}: {e}")
            continue
        knn_genes, imp_genes = _get_marker_genes_for_group(
            adata, rgg_key=rgg_key, group=grp, n_top=n_top, n_knn=n_knn
        )
        if imp_genes: all_imp_genes.update(imp_genes)
        if knn_genes: all_knn_genes.update(knn_genes)
        rgg_map[grp] = rgg_key
    
    # Convert to sorted lists for consistent ordering
    all_imp_genes = sorted(all_imp_genes)
    all_knn_genes = sorted(all_knn_genes)

    # Second pass — assign the full union gene set to every group
    for grp in groups:
        knn_g[grp] = all_knn_genes
        imp_g[grp] = all_imp_genes
        # if grp in rgg_map:  # only groups that succeeded above
        #     knn_g[grp] = all_knn_genes
        #     imp_g[grp] = all_imp_genes
    



    log.info(f"group_key={group_key} "
          f"groups_with_markers={len(knn_g)}/{len(groups)} "
          f"unique_knn_genes={len(all_knn_genes)} "
          f"unique_imp_genes={len(all_imp_genes)}")
    return knn_g, imp_g, rgg_map

# gpu worker
def _knn_worker(
    X_dense: np.ndarray,   # float32 (n_grp, n_features) — already sliced on driver
    k: int,
    metric: str,
    exclude_self: bool,
    gpu_knn_batch: int = GPU_KNN_BATCH,
) -> np.ndarray:
    """
    Fits cuML NearestNeighbors on the full group matrix, then queries in
    sub-batches of `gpu_knn_batch` rows to bound VRAM usage.
    Returns LOCAL int64 indices of shape (n_grp, k_eff).
    Falls back to sklearn on any CUDA error.
    """
    t0 = time.time()
    n_grp  = X_dense.shape[0]
    k_eff  = min(k, n_grp - 1) if exclude_self else min(k, n_grp)
    n_q    = k_eff + (1 if exclude_self else 0)

    vram_needed_gb = n_grp * X_dense.shape[1] * 4 / 1e9
    use_gpu = vram_needed_gb < VRAM_LIMIT_GB

    print(f"  [_knn_worker] n={n_grp} vram_needed={vram_needed_gb:.1f}GB "
          f"use_gpu={use_gpu}", flush=True)

    if n_grp <= 1:
        return np.zeros((n_grp, max(k_eff, 1)), dtype=np.int64)

    if use_gpu:
        try:
            import cupy as cp
            from cuml.neighbors import NearestNeighbors as cuNN

            nn = cuNN(n_neighbors=n_q, metric=metric)
            nn.fit(X_dense)

            chunks = []
            for s in range(0, n_grp, gpu_knn_batch):
                e = min(s + gpu_knn_batch, n_grp)
                chunks.append(cp.asnumpy(nn.kneighbors(X_dense[s:e], return_distance=False)))
            idx_local = np.vstack(chunks)
        except Exception as cuda_err:
            print(f"  [_knn_worker] CUDA fallback: {cuda_err}", flush=True)
            use_gpu = False
    
    if not use_gpu:
        print(f"  [_knn_worker] using sklearn CPU n={n_grp}", flush=True)
        from sklearn.neighbors import NearestNeighbors as skNN
        nn = skNN(n_neighbors=n_q, metric=metric, n_jobs=-1)
        nn.fit(X_dense)
        idx_local = nn.kneighbors(X_dense, return_distance=False)

    if exclude_self:
        idx_local = idx_local[:, 1:]

    print(f"  [_knn_worker] done n={n_grp} k={k_eff} total={time.time()-t0:.1f}s", flush=True)
    return idx_local  # LOCAL indices

# cpu worker (non-deterministic)
def _knn_worker_cpu(X_dense, k, metric, exclude_self):
    import hnswlib
    import logging, time
    log    = logging.getLogger("impute")
    t0     = time.time()
    n_grp  = X_dense.shape[0]
    k_eff  = min(k, n_grp - 1) if exclude_self else min(k, n_grp)
    n_q    = k_eff + (1 if exclude_self else 0)
    n_feat = X_dense.shape[1]
    space  = "cosine" if metric == "cosine" else "l2"

    log.info(f"[knn_worker_cpu] hnswlib build n={n_grp} space={space} n_feat={n_feat}")
    index = hnswlib.Index(space=space, dim=n_feat)
    index.init_index(
        max_elements=n_grp,
        ef_construction=100,  # higher = more accurate but slower build
        M=16,                 # higher = more accurate but more RAM
    )
    index.add_items(X_dense, num_threads=-1)
    log.info(f"[knn_worker_cpu] index built elapsed={time.time()-t0:.1f}s")

    index.set_ef(50)  # query accuracy, must be >= k
    labels, _ = index.knn_query(X_dense, k=n_q, num_threads=-1)
    log.info(f"[knn_worker_cpu] query done elapsed={time.time()-t0:.1f}s")

    # Self-removal
    if exclude_self:
        labels = labels[:, 1:]

    log.info(f"[knn_worker_cpu] done total={time.time()-t0:.1f}s")
    return labels

def _build_worker_gpu_map(client: Client) -> dict[str, int]:
    """
    Returns {worker_address: gpu_device_id} by querying each worker's
    CUDA_VISIBLE_DEVICES env.  Works with LocalCUDACluster which sets
    CUDA_VISIBLE_DEVICES=<single_id> per worker.
    """
    import os

    def _get_device():
        cvd = os.environ.get("CUDA_VISIBLE_DEVICES", "")
        try:
            return int(cvd.split(",")[0])
        except (ValueError, IndexError):
            return -1  # CPU worker

    futures = {w: client.submit(_get_device, workers=[w], pure=False)
               for w in client.scheduler_info()["workers"]}
    return {w: f.result() for w, f in futures.items()}

def _greedy_bin_pack(
    groups_by_size: list[tuple[str, int]],   # [(grp, n_cells), ...] sorted DESC
    n_bins: int,
) -> list[list[str]]:
    """
    Assign groups to `n_bins` bins (GPUs) greedily: always assign the next
    group to the least-loaded bin.  Returns list of n_bins lists of group names.
    """
    bins  = [[] for _ in range(n_bins)]
    loads = [0] * n_bins
    for grp, sz in groups_by_size:
        i = int(np.argmin(loads))
        bins[i].append(grp)
        loads[i] += sz
    return bins

# RAM weighted for parallel computing
class MemoryManager:
    def __init__(self, total_bytes):
        self.total = total_bytes
        self.available = total_bytes
        self.cond = threading.Condition()

    def acquire(self, amount):
        with self.cond:
            while self.available < amount:
                self.cond.wait()
            self.available -= amount

    def release(self, amount):
        with self.cond:
            self.available += amount
            self.cond.notify_all()

def compute_knn_indices_parallel(
    adata,
    group_key: Optional[str],
    source_layer=None,
    k: int = 15,
    metric: str = "cosine",
    exclude_self: bool = True,
    knn_genes_per_group=None,
    prefetch_workers: int = 1,  # per-worker prefetch slots
    min_free_gb: float = 20.0,
) -> tuple[np.ndarray, np.ndarray, dict]:
    from concurrent.futures import ThreadPoolExecutor

    log = TaggedLogger(logger, "knn_parallel")
    t0      = time.time()
    n_cells = adata.n_obs
    client  = get_client()

    # --- Discover workers -------------------------------------------------
    worker_gpu_map = _build_worker_gpu_map(client)
    gpu_workers    = {w: d for w, d in worker_gpu_map.items() if d >= 0}
    worker_list    = list(gpu_workers.keys())
    n_gpus         = len(worker_list)
    if n_gpus == 0:
        raise RuntimeError("No GPU workers found. Check LocalCUDACluster setup.")
    log.info(f"Found {n_gpus} GPU worker(s): "
             f"{[(w[-12:], d) for w, d in gpu_workers.items()]}")

    # --- Group setup ------------------------------------------------------
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names)} \
                  if knn_genes_per_group else None

    if group_key is None:
        groups = np.full(n_cells, "__global__", dtype=object)
    else:
        groups = np.asarray(adata.obs[group_key].astype("string").fillna("NA"))

    unique_groups = np.unique(groups)

    group_sizes = [(g, int(np.sum(groups == g))) for g in unique_groups]
    group_sizes.sort(key=lambda x: -x[1])

    # Bin-pack onto GPUs
    gpu_bins      = _greedy_bin_pack(group_sizes, n_bins=n_gpus)
    grp_to_worker = {}
    for gpu_i, grp_list in enumerate(gpu_bins):
        for grp in grp_list:
            grp_to_worker[grp] = worker_list[gpu_i]

    # --- Preallocate outputs ----------------------------------------------
    k_pad          = min(k, max(1, n_cells - 1)) if exclude_self else min(k, n_cells)
    neighbor_idx   = np.full((n_cells, k_pad), -1, dtype=np.int64)
    k_eff_per_cell = np.zeros(n_cells, dtype=np.int64)

    # --- Handle trivial groups --------------------------------------------
    non_trivial = []
    for grp, n_grp in group_sizes:
        if n_grp <= 1:
            i = np.where(groups == grp)[0][0]
            neighbor_idx[i, :] = i
            k_eff_per_cell[i]  = 0 if exclude_self else 1
        else:
            non_trivial.append((grp, n_grp))

    log.info(f"{len(non_trivial)} non-trivial groups, "
             f"prefetch_per_worker={prefetch_workers} n_gpus={n_gpus}")

    ### Weighted Semaphore For RAM ###
    RAM_LIMIT = psutil.virtual_memory().available / 1e9
    log.info(f"Creating memory manager with RAM limit {RAM_LIMIT:.1f} GB")
    mem = MemoryManager(RAM_LIMIT)

    ### Producer ###
    def _fetch_and_submit(grp, n_grp):
        knn_gene_list = knn_genes_per_group.get(str(grp)) if knn_genes_per_group else None
        valid_cols    = ([gene_to_idx[g] for g in knn_gene_list if g in gene_to_idx]
                        if knn_gene_list and gene_to_idx else None)
        n_feat      = len(valid_cols) if valid_cols else adata.n_vars

        required_gb = max(min_free_gb, n_grp * n_feat * 4 / 1e9 * 2.5) # double overhead
        target_worker = grp_to_worker.get(str(grp), worker_list[0])

        mem.acquire(required_gb)
        try:
            grp_cells = np.where(groups == grp)[0]
            k_eff     = min(k, n_grp - 1) if exclude_self else min(k, n_grp)

            log.info(f"fetch starting grp={grp} n_cells={n_grp} "
                    f"required memory={required_gb:.1f} "
                    f"free={psutil.virtual_memory().available/1e9:.1f}GB")
            X_grp = _get_source_matrix_for_cells(adata, grp_cells, source_layer)
            log.info(f"fetch done grp={grp}")

            if valid_cols:
                X_grp = X_grp[:, valid_cols]
            
            # Check VRAM — run CPU sklearn on driver if too large to scatter
            vram_needed_gb = n_grp * X_grp.shape[1] * 4 / 1e9

            def _check_vram():
                import subprocess
                r = subprocess.run(
                    ['nvidia-smi', '--query-gpu=index,memory.used,memory.free',
                    '--format=csv,noheader'],
                    capture_output=True, text=True
                )
                return r.stdout.strip()
            
            vram_status = client.submit(_check_vram, workers=[target_worker], pure=False).result()
            log.info(f"grp={grp} target_worker={target_worker} vram_needed={vram_needed_gb:.1f}GB VRAM status: {vram_status}")

            if vram_needed_gb >= VRAM_LIMIT_GB:
                log.info(f"grp={grp} too large for GPU ({vram_needed_gb:.1f} GB), "
                        f"running sklearn on driver")
                idx_local = _knn_worker_cpu(X_grp, k, metric, exclude_self)
                del X_grp
                fut = client.submit(lambda x: x, idx_local, pure=False)
                return fut, grp_cells, k_eff

            X_fut = client.scatter(X_grp, workers=[target_worker], broadcast=False)
            del X_grp
            fut = client.submit(
                _knn_worker, X_fut, k, metric, exclude_self,
                workers=[target_worker], pure=False,
            )
            X_fut.release()
            return fut, grp_cells, k_eff
        finally:
            mem.release(required_gb)

    # --- Per-worker group queues ------------------------------------------
    # Split non_trivial into per-worker lists, preserving size-sorted order
    worker_queues = {w: [] for w in worker_list}
    for grp, n_grp in non_trivial:
        target = grp_to_worker.get(str(grp), worker_list[0])
        worker_queues[target].append((grp, n_grp))
    worker_iters = {w: iter(q) for w, q in worker_queues.items()}

    # Per-worker prefetch slots: {worker -> {thread_future -> (grp, n_grp)}}
    prefetch_per_worker = {w: {} for w in worker_list}

    futures_map = {}   # dask future.key -> (grp, grp_cells, k_eff)
    active      = []
    n_submitted = 0
    n_done      = 0

    with ThreadPoolExecutor(max_workers=n_gpus * prefetch_workers) as pool:

        def _refill_worker(worker):
            """Top up prefetch slots for a specific worker."""
            pq = prefetch_per_worker[worker]
            while len(pq) < prefetch_workers:
                try:
                    grp, n_grp = next(worker_iters[worker])
                    pf = pool.submit(_fetch_and_submit, grp, n_grp)
                    pq[pf] = (grp, n_grp)
                    log.debug(f"prefetch queued grp={grp} worker={worker[-12:]}")
                except StopIteration:
                    break

        def _refill_all():
            for w in worker_list:
                _refill_worker(w)

        def _collect_prefetch(worker):
            """Move any completed prefetch futures for worker into active."""
            nonlocal n_submitted
            pq = prefetch_per_worker[worker]
            done_pf = [pf for pf in pq if pf.done()]
            for pf in done_pf:
                grp_r, n_grp_r = pq.pop(pf)
                try:
                    fut, grp_cells_r, k_eff_r = pf.result()
                    futures_map[fut.key] = (grp_r, grp_cells_r, k_eff_r)
                    active.append(fut)
                    n_submitted += 1
                    log.info(f"submitted {n_submitted}/{len(non_trivial)} "
                             f"grp={grp_r} elapsed={time.time()-t0:.1f}s")
                except Exception as e:
                    log.error(f"submit failed grp={grp_r}: {e}")
                    nonlocal n_done
                    n_done += 1

        # Seed prefetch for all workers
        _refill_all()

        while n_done < len(non_trivial):

            # Try to move any ready prefetches into active
            _refill_all()
            for w in worker_list:
                _collect_prefetch(w)

            if not active:
                log.debug(f"active empty, waiting for prefetch — "
                          f"queued={sum(len(v) for v in prefetch_per_worker.values())}")
                time.sleep(1.0)
                continue

            done_futures, still_running = wait(active, return_when="FIRST_COMPLETED")
            active = list(still_running)

            for future in done_futures:
                grp, grp_cells, k_eff = futures_map[future.key]
                try:
                    local_idx = future.result()
                    mapped = grp_cells[local_idx]
                    neighbor_idx[np.ix_(grp_cells, np.arange(k_eff))] = mapped
                    k_eff_per_cell[grp_cells] = k_eff
                except Exception as e:
                    log.error(f"KNN failed grp={grp}: {e}")
                future.release()
                n_done += 1
                log.info(f"collected {n_done}/{len(non_trivial)} "
                         f"grp={grp} elapsed={time.time()-t0:.1f}s")

                # Refill the worker that just freed up
                target_worker = grp_to_worker.get(str(grp), worker_list[0])
                _collect_prefetch(target_worker)
                _refill_worker(target_worker)

    # Force worker GC
    client.run(lambda: gc.collect())

    # --- Pad -1 entries ---------------------------------------------------
    bad_mask = neighbor_idx < 0
    if bad_mask.any():
        rows, cols = np.where(bad_mask)
        neighbor_idx[rows, cols] = rows

    if np.any(neighbor_idx < 0):
        raise RuntimeError("Negative neighbor indices remain after padding.")

    elapsed = time.time() - t0
    log.info(f"group_key={group_key} n_gpus={n_gpus} "
             f"k_eff_min={k_eff_per_cell.min()} k_eff_max={k_eff_per_cell.max()} "
             f"elapsed={elapsed:.2f}s")
    return neighbor_idx, k_eff_per_cell, {
        "group_key": group_key, "source_layer": source_layer,
        "k": int(k), "metric": metric,
        "n_gpus": n_gpus, "elapsed_sec": float(elapsed),
        "k_eff_min": int(k_eff_per_cell.min()),
        "k_eff_max": int(k_eff_per_cell.max()),
    }

# for each row, impute the feature space expression from the mean of nearest 15 neighbors
def impute_streaming(
    adata,
    knn_idx: np.ndarray,
    k_eff_per_cell: np.ndarray,
    source_layer=None,
    out_layer: str = "imputed_counts",
    imp_genes_per_group=None,
    group_key: Optional[str] = None,
    cell_batch_size: Optional[int] = 20_000,
):
    log = TaggedLogger(logger, "impute_streaming")
    log.info(f"Starting imputation")
    t0      = time.time()
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    k_pad   = knn_idx.shape[1]

    if cell_batch_size is None:
        cell_batch_size = n_cells

    # use zarr to save intermediate layers to save RAM
    out_path = os.path.join(OUT_DIR,f"{out_layer}.zarr")
    Z_out = zarr.open(
        out_path, mode="w",
        shape=(n_cells, n_genes),
        chunks=(SPARSE_CHUNK_SIZE, n_genes),
        dtype="float32",
    )

    # Pre-build selective imputation structures
    if imp_genes_per_group is not None:
        if group_key is None:
            raise ValueError("group_key required with imp_genes_per_group")
        gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}
        groups      = np.asarray(adata.obs[group_key].astype("string").fillna("NA"))
        grp_to_cols = {
            grp: np.asarray([gene_to_idx[g] for g in gl if g in gene_to_idx], dtype=np.int64)
            for grp, gl in imp_genes_per_group.items()
        }
        grp_to_cols = {k: v for k, v in grp_to_cols.items() if len(v) > 0}
    else:
        groups = grp_to_cols = None

    for start in range(0, n_cells, cell_batch_size):
        end       = min(start + cell_batch_size, n_cells)
        batch_sz  = end - start

        batch_knn  = knn_idx[start:end]          # (batch, k_pad)
        batch_keff = k_eff_per_cell[start:end]   # (batch,)
        batch_rows = np.arange(start, end)

        # All global row indices needed: batch + their neighbours
        all_needed = np.unique(np.concatenate([batch_rows, batch_knn.ravel()]))

        # Materialise only those rows
        X_needed = _get_source_matrix_for_cells(adata, all_needed, source_layer)
        # Map global → local
        g2l = np.full(n_cells, -1, dtype=np.int64)
        g2l[all_needed] = np.arange(len(all_needed), dtype=np.int64)

        # --- Vectorised KNN mean ------------------------------------------
        # Remap neighbour global indices to local
        local_knn = g2l[batch_knn]               # (batch, k_pad)

        # Replace the X_nbr block entirely
        X_mean = np.zeros((batch_sz, n_genes), dtype=np.float32)
        k_eff_safe = np.maximum(batch_keff, 1).astype(np.float32)

        for ki in range(k_pad):
            mask = batch_keff > ki           # (batch,) — cells that have a ki-th neighbour
            if not mask.any():
                break
            rows = local_knn[mask, ki]       # local indices of ki-th neighbour
            X_mean[mask] += X_needed[rows]   # accumulate — pure numpy, no 3D tensor

        X_mean /= k_eff_safe[:, None]

        # Cells with k_eff==0: fall back to source
        zero_keff = batch_keff == 0
        if np.any(zero_keff):
            src_local = g2l[batch_rows[zero_keff]]
            X_mean[zero_keff] = X_needed[src_local]

        if grp_to_cols is not None:
            batch_groups = groups[start:end]
            src_local    = g2l[batch_rows]
            X_src_batch  = X_needed[src_local]  # (batch, genes) — already in memory

            # Start from source, selectively overwrite marker cols with KNN mean
            X_out_batch = X_src_batch.copy()  # base is source

            for grp, cols in grp_to_cols.items():
                mask = batch_groups == grp
                if not np.any(mask):
                    continue
                X_out_batch[np.ix_(mask, cols)] = X_mean[np.ix_(mask, cols)]

            Z_out[start:end] = X_out_batch.astype(np.float32)
        else:
            Z_out[start:end] = X_mean.astype(np.float32)

        if (start // cell_batch_size) % 2 == 0:
            pct = 100 * end / n_cells
            log.info(f"{end:,}/{n_cells:,} rows ({pct:.1f}%) done …")
    
    # With zarr and Z_out from loop
    adata.layers[out_layer] = da.from_zarr(out_path)

    elapsed = time.time() - t0
    log.info(f"out_layer={out_layer} shape={Z_out.shape} "
          f"elapsed={elapsed:.2f}s")
    return {"out_layer": out_layer, "elapsed_sec": float(elapsed),
            "shape": tuple(Z_out.shape)}

# orchestrator
def run_hierarchical_impute_parallel(
    adata,
    levels=HIER_LEVELS,
    k: int = 15,
    metric: str = "cosine",
    exclude_self: bool = True,
    cache_key: str = CACHE_KEY_HIER3,
    n_marker_knn: int = 10,
    n_marker_imp: int = 20,
    cell_batch_size: int = 20_000,
):
    _validate_hier3_inputs(adata, levels)
    cache = _ensure_cache_hier3(adata, cache_key=cache_key)

    round_outs    = [_level_to_out_layer(lv) for lv in levels]
    round_sources = [None] + round_outs[:-1]
    round_meta    = []
    t_all         = time.time()

    prefetch_workers_per_level = {
        None: 1,                  # global — one huge group
        "neighborhood": 1,        # few very large groups
        "predicted_class": 2,     # medium groups
        "predicted_subclass": 3,  # many small groups
    }

    for i, level in enumerate(levels):
        source_layer = round_sources[i]
        out_layer    = round_outs[i]

        print(f"\n[round {i+1}/{len(levels)}] level={level} "
              f"source_layer={source_layer} -> out_layer={out_layer}")

        # Check if imputation output already exists — skip entire round
        out_path = os.path.join(OUT_DIR, f"{out_layer}.zarr")
        if os.path.exists(out_path):
            print(f"  [round {i+1}] Found existing {out_path}, loading and skipping round")
            adata.layers[out_layer] = da.from_zarr(out_path)
            if source_layer is not None and source_layer in adata.layers:
                del adata.layers[source_layer]
            continue

        # Check if KNN already computed for this round
        use_gene_subset = level in PARENT_OBS_KEY  # initialize before if/else
        imp_genes = None
        rgg_map   = {}
        knn_path = os.path.join(OUT_DIR, f"knn_{out_layer}.npz")
        if os.path.exists(knn_path):
            print(f"  [knn] Loading cached KNN from {knn_path}")
            npz = np.load(knn_path)
            knn_idx = npz["knn_idx"]
            k_eff   = npz["k_eff"]
            knn_meta = {"group_key": level, "source_layer": source_layer,
                        "k": int(k), "metric": metric, "cached": True}
        else:
            if use_gene_subset:
                knn_genes, imp_genes, rgg_map = _build_gene_dicts_for_level(
                    adata, group_key=level, n_top=n_marker_imp, n_knn=n_marker_knn
                )
            else:
                knn_genes = imp_genes = None
                rgg_map   = {}

            prefetch_w = prefetch_workers_per_level.get(level, 1)

            knn_idx, k_eff, knn_meta = compute_knn_indices_parallel(
                adata,
                group_key=level,
                source_layer=source_layer,
                k=k,
                metric=metric,
                exclude_self=exclude_self,
                knn_genes_per_group=knn_genes,
                prefetch_workers=prefetch_w,
            )
            # Save knn results for this round
            np.savez(knn_path, knn_idx=knn_idx, k_eff=k_eff)
            print(f"  [knn] saved to {knn_path}", flush=True)

        imp_meta = impute_streaming(
            adata,
            knn_idx=knn_idx,
            k_eff_per_cell=k_eff,
            source_layer=source_layer,
            out_layer=out_layer,
            imp_genes_per_group=imp_genes if use_gene_subset else None,
            group_key=level if use_gene_subset else None,
            cell_batch_size=cell_batch_size,
        )

        # Free source layer memory
        if source_layer is not None and source_layer in adata.layers:
            del adata.layers[source_layer]

        # Memory snapshot / troubleshooting
        import psutil
        rss = psutil.Process().memory_info().rss / 1e9
        print(f"\n[memory] round {i+1} complete", flush=True)
        print(f"  RSS (driver process): {rss:.1f} GB", flush=True)
        print(f"  adata.layers keys: {list(adata.layers.keys())}", flush=True)
        for name, layer in adata.layers.items():
            if isinstance(layer, np.ndarray):
                print(f"    {name}: numpy {layer.shape} {layer.nbytes/1e9:.1f} GB", flush=True)
            else:
                print(f"    {name}: {type(layer).__name__} {layer.shape}", flush=True)
        print(f"  knn_idx: {knn_idx.nbytes/1e9:.1f} GB", flush=True)
        print(f"  k_eff:   {k_eff.nbytes/1e6:.1f} MB", flush=True)

        round_meta.append({
            "level": level, "source_layer": source_layer, "out_layer": out_layer,
            "knn": knn_meta, "impute": imp_meta, "rgg_key_map": rgg_map,
        })

    total = time.time() - t_all
    cache.update({
        "levels": [str(lv) if lv is not None else None for lv in levels],
        "rounds": round_meta,
        "params": {
            "k": int(k), "metric": metric, "exclude_self": bool(exclude_self),
            "n_marker_knn": int(n_marker_knn), "n_marker_imp": int(n_marker_imp),
            "cell_batch_size": cell_batch_size if cell_batch_size is not None else adata.n_obs,
        },
        "elapsed_sec": float(total),
    })
    print(f"\n[run_hierarchical_impute_parallel] total elapsed={total:.2f}s")
    return round_meta

if __name__ == "__main__":
    config = load_config(sys.argv[1])
    print(sys.argv[1])
    SPARSE_CHUNK_SIZE = config["sparse_chunk_size"]
    PARENT_OBS_KEY    = config.get("parent_obs_key") or {}
    CACHE_KEY_HIER3   = config["cache_key_hier3"]
    HIER_LEVELS       = config["hier_levels"]
    SEED              = config["seed"]
    CELL_BATCH_SIZE   = config["cell_batch_size"]
    OUT_DIR           = config["out_dir"]
    os.makedirs(OUT_DIR, exist_ok=True)
    logging.getLogger().addHandler(logging.FileHandler(os.path.join(OUT_DIR, "run.log")))

    cluster = LocalCUDACluster(
        CUDA_VISIBLE_DEVICES="0,1",
        threads_per_worker=1,
        protocol="tcp",
        rmm_allocator_external_lib_list="cupy",
        memory_limit="auto",
    )
    client = Client(cluster)
    print(client)

    if config["use_previous_results"]:
        print("Using previous results")
        results_h5ad_path = config["previous_results_path"]
        adata = ad.read_h5ad(results_h5ad_path)
    else:
        data_pth = config["data_path"]
        f = zarr.open(data_pth, mode="r", use_consolidated=False)
        X_raw     = f["X"]
        shape = getattr(X_raw, "shape", None)
        if shape is None:
            shape = X_raw.attrs["shape"]
        print(f"Dataset shape: {shape}")

        adata = ad.AnnData(
            X    = read_dask(X_raw, (SPARSE_CHUNK_SIZE, shape[1])),  # stays dask
            obs  = ad.io.read_elem(f["obs"]),
            var  = ad.io.read_elem(f["var"]),
            obsm = ad.io.read_elem(f["obsm"]),
            uns  = ad.io.read_elem(f["uns"]),
        )
        # No adata.X.compute() — lazy throughout

        # Load use_genes from the provided file path
        use_genes_path = config["use_genes_path"]
        with open(use_genes_path, 'r') as f:
            # Split, strip, and ignore any blank lines
            use_genes = [line.strip() for line in f if line.strip()]

        adata_rank = ad.read_h5ad(config["adata_rank_path"])
        for key, val in adata_rank.uns.items():
            if key not in adata.uns:
                adata.uns[key] = val
        del adata_rank

        if config["filter"]:
            adata = adata[(adata.obs["scDblFinder.class"]=="singlet") 
                & (adata.obs["num_genes"]>=2000)
                & (adata.obs["log_ambient_mse_norm"]>0.09)].copy()
            is_gex = ~adata.var_names.str.contains("_")
            adata = adata[:,is_gex].copy()
            adata = adata[
                np.lexsort((
                    adata.obs["predicted_subclass"].values,
                    adata.obs["predicted_class"].values,
                    adata.obs["neighborhood"].values,
                ))
            ].copy()

        # # run normalization and          
        if config.get("normalize", True) and "log1p" not in adata.uns:
            out_zarr_path = os.path.join(OUT_DIR, "lognorm.zarr")
            NORMALIZE_PATH = config["normalize_path"] if "normalize_path" in config else None
            if NORMALIZE_PATH is not None:
                print("Using normalized path as X")
                adata.var["marker_genes"] = adata.var_names.isin(use_genes)
                adata = adata[:,adata.var.marker_genes].copy()
                adata.X = da.from_zarr(NORMALIZE_PATH)
            else:
                print("Log normalizing and subsetting to markers")
                adata.var["marker_genes"] = adata.var_names.isin(use_genes)
                rsc.get.anndata_to_GPU(adata)
                rsc.pp.normalize_total(adata,target_sum = 1e4)
                rsc.pp.log1p(adata)
                adata = adata[:,adata.var.marker_genes].copy()
                rsc.get.anndata_to_CPU(adata)
                adata.X = adata.X.persist()
                if config.get("overwrite", False):
                    if os.path.exists(out_zarr_path):
                        print("Overwrite is true, removing log normalized subdirectory")
                        shutil.rmtree(out_zarr_path)
                if not os.path.exists(out_zarr_path):
                    print(f"Writing filtered+sorted adata to {out_zarr_path}...")
                    t0 = time.time()
                    # Open the store and create X dataset
                    Z_X = zarr.open(
                        out_zarr_path, mode="w",
                        shape=(adata.n_obs, adata.n_vars),
                        chunks=(SPARSE_CHUNK_SIZE, adata.n_vars),
                        dtype="float32",
                    )
                    if isinstance(adata.X, da.Array):
                        print("  X is dask — streaming to zarr via da.store")
                        X_dense = adata.X.map_blocks(
                            lambda block: block.toarray().astype(np.float32) if sp.issparse(block) else np.asarray(block, dtype=np.float32),
                            dtype=np.float32,
                        )
                        da.store(X_dense, Z_X)
                        # Add anndata encoding metadata so read_elem_lazy can read it
                        Z_X.attrs["encoding-type"] = "array"
                        Z_X.attrs["encoding-version"] = "0.2.0"
                    elif sp.issparse(adata.X):
                        print("  X is sparse — writing chunk by chunk")
                        X_csr = adata.X.tocsr()
                        for start in range(0, adata.n_obs, SPARSE_CHUNK_SIZE):
                            end   = min(start + SPARSE_CHUNK_SIZE, adata.n_obs)
                            chunk = X_csr[start:end].toarray().astype(np.float32)
                            Z_X[start:end] = chunk
                            if (start // SPARSE_CHUNK_SIZE) % 20 == 0:
                                pct = 100 * end / adata.n_obs
                                print(f"  [write X] {end:,}/{adata.n_obs:,} ({pct:.1f}%)", flush=True)
                    else:
                        print("  X is dense numpy — writing chunk by chunk")
                        X_np = np.asarray(adata.X, dtype=np.float32)
                        for start in range(0, adata.n_obs, SPARSE_CHUNK_SIZE):
                            end   = min(start + SPARSE_CHUNK_SIZE, adata.n_obs)
                            Z_X[start:end] = X_np[start:end]
                            if (start // SPARSE_CHUNK_SIZE) % 20 == 0:
                                pct = 100 * end / adata.n_obs
                                print(f"  [write X] {end:,}/{adata.n_obs:,} ({pct:.1f}%)", flush=True)
                    elapsed = time.time() - t0
                    print(f"  Done writing lognorm zarr in {elapsed/60:.2f}min ({elapsed:.1f}s)", flush=True)
                    adata.X = da.from_zarr(out_zarr_path)
                else:
                    print(f"Found existing filtered zarr at {out_zarr_path}, skipping write")
                    adata.X = da.from_zarr(out_zarr_path)
                    print(adata.X.shape)
                    print(adata.X.chunks)
        else:
            adata.var["marker_genes"] = adata.var_names.isin(use_genes)
            adata = adata[:,adata.var.marker_genes].copy()

    round_meta = run_hierarchical_impute_parallel(
        adata,
        levels=HIER_LEVELS,
        k=config["knn"]["k"],
        metric=config["knn"]["metric"],
        exclude_self=False,
        n_marker_knn=config["knn"]["n_marker_knn"],
        n_marker_imp=config["knn"]["n_marker_imp"],
        cell_batch_size=CELL_BATCH_SIZE,
    )
    
    if config["downstream"]:
        log = TaggedLogger(logger, "downstream")
        today = datetime.now().strftime("%y%m%d_%H%M")
        downstream_out_dir = os.path.join(OUT_DIR, today)
        os.makedirs(downstream_out_dir, exist_ok=True)
        sc.settings.figdir = downstream_out_dir
        shutil.copy(sys.argv[1], os.path.join(downstream_out_dir, os.path.basename(sys.argv[1])))

        final_out_layer = _level_to_out_layer(HIER_LEVELS[-1])
        adata.X = adata.layers[final_out_layer]

        # If it's a numpy array, wrap it as dask first so rsc can chunk it
        if not isinstance(adata.X, da.Array):
            adata.X = da.from_array(adata.X, chunks=(SPARSE_CHUNK_SIZE, adata.n_vars))

        rng = np.random.default_rng(SEED)
        cl = adata.obs["predicted_cluster"].astype("string").fillna("NA")
        nbhd = adata.obs["neighborhood"].astype("string").fillna("NA")
        cluster_to_cells = cl.groupby(cl).indices

        print(list(cluster_to_cells.keys())[:40])

        pca_idx  = []
        umap_idx = []

        def sample_size(n, max_cells):
            min_cells = 50 # change from 50 or so 03/10 12pm
            return min(n, int(np.clip(80 * n**0.6, min_cells, max_cells)))
                
        SUBSAMPLING = config.get("subsampling_scheme", "yao")

        print(SUBSAMPLING)

        if SUBSAMPLING == "yao":
            for idx in cluster_to_cells.values():
                perm = rng.permutation(idx)
                pca_idx.append(perm[:100])
                umap_idx.append(perm[:1000])
        elif SUBSAMPLING == "yao_custom":
            for idx in cluster_to_cells.values():
                perm = rng.permutation(idx)
                pca_idx.append(perm[:100])
                umap_idx.append(perm[:4000])
                #umap_idx.append(perm[:6000])
        elif SUBSAMPLING == "proportional_v1":
            for idx in cluster_to_cells.values():   
                perm = rng.permutation(idx)
                pca_idx.append(perm[:100])
                umap_sample = sample_size(len(idx), 6000)
                umap_idx.append(perm[:umap_sample])
        elif SUBSAMPLING == "proportional":
            for neighborhood, group in adata.obs.groupby(nbhd):
                max_cells = 6_000
                for cluster, sub in group.groupby(cl.loc[group.index]):
                    idx = adata.obs_names.get_indexer(sub.index)
                    perm = rng.permutation(idx)
                    pca_idx.append(perm[:100])
                    umap_idx.append(perm[:sample_size(len(idx), max_cells)])
        elif SUBSAMPLING == "rare":
            for neighborhood, group in adata.obs.groupby(nbhd):
                if neighborhood in ["MB-HB-CB-GABA", "HY-EA-Glut-GABA", "Subpallium-GABA; HY-EA-Glut-GABA"]:
                    max_cells = 12_000
                else:
                    max_cells = 6_000
                for cluster, sub in group.groupby(cl.loc[group.index]):
                    idx = adata.obs_names.get_indexer(sub.index)
                    perm = rng.permutation(idx)
                    pca_idx.append(perm[:100])
                    umap_idx.append(perm[:sample_size(len(idx), max_cells)])


        pca_idx  = np.sort(np.concatenate(pca_idx))
        umap_idx = np.sort(np.concatenate(umap_idx))

        log.info(f"PCA Subsample Size: {len(pca_idx)}, Slice: {pca_idx}")
        log.info(f"UMAP Subsample Size: {len(umap_idx)}, Slice: {umap_idx}")

        if "X_pca" not in adata.obsm:
            if config.get("subsample_pca", True):
                log.info("Calcluating PCA with stratified subsample")
                adata_pca = adata[pca_idx].copy()
                # PCA + neighbors + UMAP
                rsc.get.anndata_to_GPU(adata_pca)
                start_time = time.time()
                rsc.pp.pca(adata_pca, n_comps=100, mask_var=None, random_state=SEED)
                rsc.get.anndata_to_CPU(adata_pca)
                adata_pca.obsm["X_pca"]=adata_pca.obsm["X_pca"].compute()
                elapsed_time = time.time() - start_time
                log.info(f"PCA calculation completed in {elapsed_time:.2f} seconds")

                log.info("Projecting PCA")
                start_time = time.time()
                PCs = adata_pca.varm["PCs"]        # genes × PCs
                mean = adata.X.mean(axis=0).compute()  # full-dataset mean for correct centering

                PCs_da = da.from_array(PCs, chunks=(PCs.shape[0], -1))
                Xp = adata.X @ PCs_da
                adata.obsm["X_pca"] = Xp - (mean @ PCs)
                adata.obsm["X_pca"] = adata.obsm["X_pca"].compute()
                elapsed_time = time.time() - start_time
                log.info(f"PCA projection completed in {elapsed_time:.2f} seconds")

                adata.X = None
                adata.uns.pop("impute_knn_cache_hier3", None)
                adata.layers.pop(final_out_layer, None)
                adata.write_h5ad(os.path.join(downstream_out_dir, "pca_results.h5ad"))
            elif config.get("subsample_pca_nbhd", False):
                # Cap contribution per neighborhood regardless of cluster count
                n_per_neighborhood = 50_000
                rng = np.random.default_rng(SEED)
                sampled_idx = []

                for nbhd in adata.obs["neighborhood"].unique():
                    mask = np.where(adata.obs["neighborhood"] == nbhd)[0]
                    n_take = min(len(mask), n_per_neighborhood)
                    sampled_idx.append(rng.choice(mask, n_take, replace=False))

                sampled_idx = np.sort(np.concatenate(sampled_idx))

                adata_pca = adata[sampled_idx].copy()
                adata_pca.X = adata_pca.X.persist()

                # PCA + neighbors + UMAP
                rsc.get.anndata_to_GPU(adata_pca)

                print("Calcluating PCA with nbhd subsample")
                start_time = time.time()
                rsc.pp.pca(adata_pca, n_comps=100, mask_var=None, random_state=SEED)
                rsc.get.anndata_to_CPU(adata_pca)
                adata_pca.obsm["X_pca"]=adata_pca.obsm["X_pca"].compute()
                elapsed_time = time.time() - start_time
                print(f"PCA calculation completed in {elapsed_time:.2f} seconds")

                # Project all 1.7M cells
                PCs = adata_pca.varm["PCs"]
                mean = adata_pca.X.mean(axis=0)
                PCs_da = da.from_array(PCs, chunks=(PCs.shape[0], -1))
                adata.obsm["X_pca"] = (adata.X @ PCs_da - (mean @ PCs)).compute()
                del adata_pca
            else:
                print("Calcluating PCA, no subsample")
                start_time = time.time()
                rsc.get.anndata_to_GPU(adata)
                rsc.pp.pca(adata, n_comps=100, mask_var=None, random_state=SEED)
                rsc.get.anndata_to_CPU(adata)
                adata.obsm["X_pca"]=adata.obsm["X_pca"].compute()
                elapsed_time = time.time() - start_time
                print(f"PCA calculation completed in {elapsed_time:.2f} seconds")
            # X_centered = adata.X - mean
            # adata.obsm["X_pca"] = X_centered @ PCs
            # adata.obsm["X_pca"]=adata.obsm["X_pca"].compute()
        else:
            log.info("Using PCA from previous results")

        UMAP_METRIC = config["umap_metric"]
        UMAP_NEIGHBORS = config["umap_neighbors"]
        UMAP_MIN_DIST = config.get("umap_min_dist", 0.4)

        if config.get("subsample_umap", True):
            adata_umap = adata[umap_idx].copy()

            # --- Neighbors + UMAP on subsample -----------------------------------
            log.info("Calculating neighbors on stratified subsample")
            t0 = time.time()
            rsc.pp.neighbors(adata_umap, use_rep="X_pca", n_neighbors=UMAP_NEIGHBORS,
                            metric=UMAP_METRIC, random_state=SEED)
            log.info(f"Neighbors done elapsed={time.time()-t0:.1f}s")

            log.info("Calculating UMAP on subsample")
            t0 = time.time()
            rsc.tl.umap(adata_umap, min_dist=UMAP_MIN_DIST, random_state=SEED, init_pos="spectral")
            log.info(f"UMAP calculation completed in {time.time()-t0:.1f}s")

            plt_categories = config["plt_categories"]
            for cat in plt_categories:
                sc.pl.umap(adata_umap, color=cat, save=f"_{cat}_sub.pdf", size=0.15, alpha=0.6)

            umap      = adata_umap.obsm["X_umap"]
            center    = umap.mean(axis=0)
            dist      = np.linalg.norm(umap - center, axis=1)
            cutoff    = np.percentile(dist, 99.9)
            keep_mask = dist <= cutoff
            adata_umap = adata_umap[keep_mask].copy()

            for cat in plt_categories:
                sc.pl.umap(adata_umap, color=cat, save=f"_{cat}_flt.pdf", size=0.15, alpha=0.6)

            # --- Re-derive umap_idx from obs_names after all filtering -----------
            obs_name_to_idx = {name: i for i, name in enumerate(adata.obs_names)}
            umap_idx        = np.array([obs_name_to_idx[name]
                                        for name in adata_umap.obs_names])

            # --- Project full dataset into subsample UMAP ------------------------
            pca_full   = adata.obsm["X_pca"]
            pca_sub    = adata_umap.obsm["X_pca"]
            umap_sub   = adata_umap.obsm["X_umap"]

            non_sub_mask           = np.ones(adata.n_obs, dtype=bool)
            non_sub_mask[umap_idx] = False
            non_sub_global         = np.where(non_sub_mask)[0]

            log.info(f"pca_sub shape: {pca_sub.shape}, dtype: {pca_sub.dtype}")
            log.info(f"pca_sub nbytes: {pca_sub.nbytes / 1024 / 1024:.2f} MB")

            import cupy as cp
            from cuml.neighbors import NearestNeighbors as cuNN

            nn = cuNN(n_neighbors=UMAP_NEIGHBORS, metric="euclidean", algorithm="brute")
            nn.fit(pca_sub.astype(np.float32))

            pca_nonsub = pca_full[non_sub_global].astype(np.float32)
            n_nonsub   = len(non_sub_global)
            all_idx    = np.empty((n_nonsub, UMAP_NEIGHBORS), dtype=np.int64)

            t0 = time.time()
            batch_sz = 100_000
            for s in range(0, n_nonsub, batch_sz):
                e            = min(s + batch_sz, n_nonsub)
                query_batch  = cp.asarray(pca_nonsub[s:e])
                I            = nn.kneighbors(query_batch, n_neighbors=UMAP_NEIGHBORS,
                                            return_distance=False)
                all_idx[s:e] = cp.asnumpy(I)
                log.info(f"KNN projection {e}/{n_nonsub} elapsed={time.time()-t0:.1f}s")

            umap_projected = umap_sub[all_idx].mean(axis=1).astype(np.float32)

            umap_full                 = np.empty((adata.n_obs, 2), dtype=np.float32)
            umap_full[umap_idx]       = umap_sub
            umap_full[non_sub_global] = umap_projected
            adata.obsm["X_umap"]      = umap_full

            for cat in plt_categories:
                sc.pl.umap(adata, color=cat, save=f"_{cat}_full.pdf", size=0.15, alpha=0.6)

            # --- Filter outliers on full projected UMAP --------------------------
            center    = umap_full.mean(axis=0)
            dist      = np.linalg.norm(umap_full - center, axis=1)
            cutoff    = np.percentile(dist, 99.9)
            keep_mask = dist <= cutoff
            adata     = adata[keep_mask].copy()
            log.info(f"UMAP filter removed {(~keep_mask).sum()} cells, "
                    f"{keep_mask.sum()} remaining")

            for cat in plt_categories:
                sc.pl.umap(adata, color=cat, save=f"_{cat}_full_flt.pdf", size=0.15, alpha=0.6)
        else:
            #adata.X = adata.X.persist()
            print("Calcluating Neighbors")
            start_time = time.time()
            rsc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=UMAP_NEIGHBORS, metric=UMAP_METRIC, random_state=SEED)
            elapsed_time = time.time() - start_time
            print(f"Neighbor calculation completed in {elapsed_time:.2f} seconds")

            print("Calcluating UMAP")
            start_time = time.time()
            rsc.tl.umap(adata, min_dist=UMAP_MIN_DIST, random_state=SEED, init_pos="spectral")
            elapsed_time = time.time() - start_time
            print(f"UMAP calculation completed in {elapsed_time:.2f} seconds")

            plt_categories = config["plt_categories"]

            print("Starting Plots")
            for cat in plt_categories:
                sc.pl.umap(adata,
                    color=cat,
                    save=f"_{cat}.pdf")

            UMAP_FILTER = config.get("umap_filter", False)
            if UMAP_FILTER:
                umap = adata.obsm["X_umap"]
                center = umap.mean(axis=0)
                dist = np.linalg.norm(umap - center, axis=1)
                cutoff = np.percentile(dist, 99.9)
                keep_mask = dist <= cutoff
                adata = adata[keep_mask].copy()
                print("Starting Plots")
                for cat in plt_categories:
                    sc.pl.umap(adata,
                        color=cat,
                        save=f"_{cat}_flt.pdf")

        # save only metadata needed to recreate UMAP
        adata.X = None
        adata.uns.pop("impute_knn_cache_hier3", None)
        adata.layers.pop(final_out_layer, None)
        adata.write_h5ad(os.path.join(downstream_out_dir, "results.h5ad"))
        print("Finished Saving")
        client.shutdown() 
        cluster.close()