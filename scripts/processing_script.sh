#!/bin/bash
# Process multiple batches of wholebrain_crispr_atlas data: process, QC, merge, and final QC

set -e  # Exit on any error
set -u  # Exit on undefined variable

# Run from script directory so process_adata.py, filter_for_deg.py, merge_batches.py are found
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# Configuration
NUM_BATCHES=${1:-4}  # Number of batches
BASE_DATA_PATH="/data/wholebrain_crispr_atlas"
SCREEN_NAME="WB8588_screen"

ORIGINAL_DIR="${BASE_DATA_PATH}/h5ads"
PROCESSED_DIR="${BASE_DATA_PATH}/processed"

# QC parameters for individual batches
MIN_CELLS_PER_GENE=0
MIN_TOTAL_COUNTS=0
MAX_PCT_MT=100
MIN_CELLS_PER_SUBCLASS=0
MIN_GENES_PER_CELL=2000
LOG_AMBIENT_MSE_NORM_MIN=0.09
UNASSIGNED_LABEL="Negative"

# QC parameters for final merged file
FINAL_MIN_CELLS_PER_GENE=50
FINAL_MIN_CELLS_PER_SUBCLASS=0

echo "=========================================="
echo "Processing ${NUM_BATCHES} batches"
echo "Base path: ${BASE_DATA_PATH}"
echo "=========================================="
echo ""

# Step 1: Process each batch
echo "Step 1: Processing individual batches..."
for i in $(seq 1 ${NUM_BATCHES}); do
    echo "Processing batch ${i}..."
    python process_adata.py \
        --data_path "${ORIGINAL_DIR}/WB8588_${i}_*.h5ad" \
        --output_dir "${PROCESSED_DIR}/batch${i}" \
        --output_name "${SCREEN_NAME}_batch${i}_gex.h5ad" \
        --guide_output_name "${SCREEN_NAME}_batch${i}_crispr_guides.h5ad"
    echo "✓ Batch ${i} processed"
done
echo ""

# Step 2: QC for each batch
echo "Step 2: Filter for DEG for individual batches..."
for i in $(seq 1 ${NUM_BATCHES}); do
    echo "QC for batch ${i}..."
    python filter_for_deg.py \
        --data_path "${PROCESSED_DIR}/batch${i}/${SCREEN_NAME}_batch${i}_gex.h5ad" \
        --output_dir "${PROCESSED_DIR}/batch${i}" \
        --output_name "${SCREEN_NAME}_batch${i}_gex_filtered.h5ad" \
        --remove_multiple_guides \
        --remove_doublets \
        --remove_unassigned \
        --min_cells_per_gene ${MIN_CELLS_PER_GENE} \
        --min_genes_per_cell ${MIN_GENES_PER_CELL} \
        --min_total_counts ${MIN_TOTAL_COUNTS} \
        --max_pct_mt ${MAX_PCT_MT} \
        --min_cells_per_subclass ${MIN_CELLS_PER_SUBCLASS} \
        --log_ambient_mse_norm_min ${LOG_AMBIENT_MSE_NORM_MIN} \
        --unassigned_label "${UNASSIGNED_LABEL}" \
        --write_h5ad
    echo "✓ Batch ${i} QC completed"
done
echo ""

# Step 3: Merge all batches
echo "Step 3: Merging all batches..."
MERGE_CMD="python merge_batches.py"
for i in $(seq 1 ${NUM_BATCHES}); do
    MERGE_CMD="${MERGE_CMD} --data_path ${PROCESSED_DIR}/batch${i}/${SCREEN_NAME}_batch${i}_gex_filtered.h5ad"
done
MERGE_CMD="${MERGE_CMD} --output_dir ${PROCESSED_DIR} --output_name ${SCREEN_NAME}_gex_filtered.h5ad" --write_zarr True --write_h5ad True

echo "Executing: ${MERGE_CMD}"
eval ${MERGE_CMD}
echo "✓ Batches merged"
echo ""

# Step 4: Final QC on merged file
echo "Step 4: Final QC on merged file..."
python filter_for_deg.py \
    --data_path "${PROCESSED_DIR}/${SCREEN_NAME}_gex_filtered.h5ad" \
    --output_dir "${PROCESSED_DIR}/" \
    --output_name "${SCREEN_NAME}_gex_filtered.h5ad" \
    --min_cells_per_gene ${FINAL_MIN_CELLS_PER_GENE} \
    --min_cells_per_subclass ${FINAL_MIN_CELLS_PER_SUBCLASS} \
    --remove_unassigned \
    --unassigned_label "${UNASSIGNED_LABEL}" \
    --write_zarr
echo "✓ Final QC completed"
echo ""

echo "=========================================="
echo "All processing steps completed successfully!"
echo "=========================================="