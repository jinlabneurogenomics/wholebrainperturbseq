#!/usr/bin/env python
"""
download.py

Downloads the wholebrain CRISPR atlas from HuggingFace: metadata, analysis folder,
and h5ad shards. Assemble raw zarr and run preprocess via umap_preprocess.py.

Usage:
    python download.py --output-dir /data/wholebrain_crispr_atlas [--hf-token TOKEN]
    python download.py --output-dir /data/wholebrain_crispr_atlas --skip-h5ads   # metadata + analysis only

Requirements:
    pip install huggingface_hub
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import textwrap

# ---------------------------------------------------------------------------
# HuggingFace dataset info
# ---------------------------------------------------------------------------
HF_REPO_ID = "perturbai/wholebrain_crispr_atlas"
HF_REPO_TYPE = "dataset"

METADATA_FILES = [
    "metadata/all_obs.parquet",
]

H5AD_FILES = [
    "h5ads/WB8588_1_1.h5ad",
    "h5ads/WB8588_1_2.h5ad",
    "h5ads/WB8588_1_3.h5ad",
    "h5ads/WB8588_2_1.h5ad",
    "h5ads/WB8588_2_2.h5ad",
    "h5ads/WB8588_2_3.h5ad",
    "h5ads/WB8588_3_1.h5ad",
    "h5ads/WB8588_3_2.h5ad",
    "h5ads/WB8588_3_3.h5ad",
    "h5ads/WB8588_3_4.h5ad",
    "h5ads/WB8588_4_1.h5ad",
    "h5ads/WB8588_4_2.h5ad",
    "h5ads/WB8588_4_3.h5ad",
    "h5ads/WB8588_4_4.h5ad",
]


# ---------------------------------------------------------------------------
# Step 1: Download metadata from HuggingFace
# ---------------------------------------------------------------------------
def download_metadata(output_dir: Path, hf_token: str | None = None) -> list[Path]:
    try:
        from huggingface_hub import hf_hub_download
    except ImportError:
        sys.exit("huggingface_hub not installed. Run: pip install huggingface_hub")

    meta_dir = output_dir / "metadata"
    meta_dir.mkdir(parents=True, exist_ok=True)

    local_paths = []
    for hf_path in METADATA_FILES:
        local_path = output_dir / hf_path
        if local_path.exists():
            print(f"  [skip] {local_path.name} already present")
        else:
            print(f"  Downloading {hf_path} ...")
            downloaded = hf_hub_download(
                repo_id=HF_REPO_ID,
                filename=hf_path,
                repo_type=HF_REPO_TYPE,
                local_dir=str(output_dir),
                token=hf_token,
            )
            print(f"  -> {downloaded}")
        local_paths.append(local_path)

    return local_paths


# ---------------------------------------------------------------------------
# Step 2: Download analysis folder from HuggingFace
# ---------------------------------------------------------------------------
def download_analysis(output_dir: Path, hf_token: str | None = None) -> Path:
    try:
        from huggingface_hub import snapshot_download
    except ImportError:
        sys.exit("huggingface_hub not installed. Run: pip install huggingface_hub")

    print("  Downloading analysis folder ...")
    snapshot_download(
        repo_id=HF_REPO_ID,
        repo_type=HF_REPO_TYPE,
        local_dir=str(output_dir),
        allow_patterns=["analysis/*"],
        token=hf_token,
    )
    return output_dir / "analysis"


# ---------------------------------------------------------------------------
# Step 3: Download h5ad shards from HuggingFace
# ---------------------------------------------------------------------------
def download_h5ads(output_dir: Path, hf_token: str | None = None) -> list[Path]:
    try:
        from huggingface_hub import hf_hub_download
    except ImportError:
        sys.exit("huggingface_hub not installed. Run: pip install huggingface_hub")

    h5ad_dir = output_dir / "h5ads"
    h5ad_dir.mkdir(parents=True, exist_ok=True)

    local_paths = []
    for hf_path in H5AD_FILES:
        local_path = h5ad_dir / Path(hf_path).name
        if local_path.exists():
            print(f"  [skip] {local_path.name} already present")
        else:
            print(f"  Downloading {hf_path} ...")
            downloaded = hf_hub_download(
                repo_id=HF_REPO_ID,
                filename=hf_path,
                repo_type=HF_REPO_TYPE,
                local_dir=str(output_dir),
                token=hf_token,
            )
            print(f"  -> {downloaded}")
        local_paths.append(local_path)

    return local_paths


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Download wholebrain CRISPR atlas h5ad shards from HuggingFace",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
            Example:
              python download.py --output-dir /data/wholebrain

            Then run umap_preprocess.py to assemble zarr and run preprocess.
        """),
    )
    p.add_argument("--output-dir", required=True,
                   help="Root directory for downloaded data and outputs")
    p.add_argument("--hf-token", default=None,
                   help="HuggingFace token (not needed for public datasets)")
    p.add_argument("--skip-h5ads", action="store_true",
                   help="Skip downloading h5ad shards (e.g. to fetch only metadata and analysis)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Download metadata from HuggingFace")
    print("=" * 60)
    download_metadata(output_dir, hf_token=args.hf_token)

    print("\n" + "=" * 60)
    print("Download analysis folder from HuggingFace")
    print("=" * 60)
    download_analysis(output_dir, hf_token=args.hf_token)

    if not args.skip_h5ads:
        print("\n" + "=" * 60)
        print("Download h5ad shards from HuggingFace")
        print("=" * 60)
        download_h5ads(output_dir, hf_token=args.hf_token)
        print("\nDone. Run umap_preprocess.py to assemble raw zarr and run preprocess.")
    else:
        print("\nSkipped h5ad shards (--skip-h5ads). Done.")


if __name__ == "__main__":
    main()
