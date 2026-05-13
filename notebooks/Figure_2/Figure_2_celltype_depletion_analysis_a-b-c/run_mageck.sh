#!/usr/bin/env bash
# Run MAGeCK2 test for each cell type (count table: observed vs expected).
# Prepends RRA to PATH if the binary exists (edit RRA_BIN if your install path differs).
set -e

MAGECK_DIR="/workspace/mageck2"

# Ensure mageck2 is available; if not, add MAGECK_DIR to PATH or install
if ! command -v mageck2 &>/dev/null; then
  if [[ -x "${MAGECK_DIR}/bin/mageck2" ]]; then
    export PATH="${MAGECK_DIR}/bin:${PATH}"
  fi
fi
if ! command -v mageck2 &>/dev/null; then
  echo "MAGeCK2 not found. Attempting installation into ${MAGECK_DIR} ..."
  if [[ -f "${MAGECK_DIR}/setup.py" ]]; then
    pip install -e "${MAGECK_DIR}"
  else
    mkdir -p "$(dirname "${MAGECK_DIR}")"
    git clone --depth 1 https://github.com/davidliwei/mageck2.git "${MAGECK_DIR}"
    pip install -e "${MAGECK_DIR}"
    (cd "${MAGECK_DIR}/rra" && make) 2>/dev/null || true
  fi
  export PATH="${MAGECK_DIR}/bin:${PATH}"
fi
if ! command -v mageck2 &>/dev/null; then
  echo "ERROR: MAGeCK2 could not be installed. Install manually and ensure mageck2 is on PATH." >&2
  exit 1
fi
RRA_BIN="/workspace/mageck2/rra/bin/RRA"
if [[ -x "$RRA_BIN" ]]; then
  export PATH="$(dirname "$RRA_BIN"):$PATH"
fi

if [ $# -lt 1 ]; then
  echo "Usage: $0 <OUT_DIR>" >&2
  exit 1
fi
OUT_DIR="$1"
BASE_INPUT="${OUT_DIR}/mageck_input"
BASE_OUTPUT="${OUT_DIR}/mageck_results"

CELL_IDS=(Pallium-Glut TH-EPI-Glut Subpallium-GABA MB-HB-Glut-Sero-Dopa MB-HB-CB-GABA HY-EA-Glut-GABA)

for id in "${CELL_IDS[@]}"; do
  count_table="${BASE_INPUT}/${id}.txt"
  out_prefix="${BASE_OUTPUT}/${id}/mageck_results"
  out_dir="${BASE_OUTPUT}/${id}"
  if [[ ! -f "$count_table" ]]; then
    echo "SKIP $id: count table not found: $count_table"
    continue
  fi
  mkdir -p "$out_dir"
  echo "Running MAGeCK2 test for cell $id ..."
  mageck2 test -k "$count_table" -t observed -c expected -n "$out_prefix"
  echo "Done $id -> $out_prefix.gene_summary.txt"
done

echo "All done."