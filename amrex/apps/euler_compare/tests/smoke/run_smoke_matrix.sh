#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
cd "$ROOT"

cd amrex/apps/euler_compare
make -j2 DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE
cd ../../..

OUTPUT="${SMOKE_OUTPUT_DIR:-${1:-$ROOT/../euler_compare_smoke_$(date +%Y%m%d_%H%M%S)}}"
mkdir -p "$OUTPUT"
python3 scripts/run_project_smoke_matrix.py --output "$OUTPUT"
echo "smoke_output=$OUTPUT"
