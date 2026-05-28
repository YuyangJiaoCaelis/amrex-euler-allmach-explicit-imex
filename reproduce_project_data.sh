#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

build() {
  (cd amrex/apps/euler_compare && make -j2 DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE)
}

riemann() {
  python3 scripts/run_riemann_1d_convergence.py
  python3 scripts/plot_riemann_1d_efficiency.py
  python3 scripts/plot_riemann_1d_overlays.py
  python3 scripts/plot_riemann_1d_refinement_zoom.py
}

gresho() {
  python3 scripts/run_gresho_report_rows.py --row-timeout-sec 7200 --force
  python3 scripts/combine_gresho_reproduction_data.py
  python3 scripts/plot_gresho_context_figures.py
  python3 scripts/plot_gresho_density_error.py
}

advection() {
  python3 scripts/run_advection_blob_efficiency.py
  python3 scripts/plot_advection_blob_efficiency.py
  python3 scripts/plot_advection_blob_timesnaps.py
}

shock_density_bubble() {
  EXE=./amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex
  INPUT=amrex/apps/euler_compare/inputs-ci
  run_shock_explicit() {
    local nx=$1
    local ny=$2
    local root=$3
    local riemann=$4
    local prefix=$5
    local row_id="${prefix}_${nx}x${ny}"
    local final_csv="$root/${row_id}_summary.csv"
    local snapshot_dir="$root/${row_id}_snapshots"
    mkdir -p "$root"
    python3 scripts/run_manifest.py run \
      --root "$ROOT" \
      --row-id "$row_id" \
      --output-dir "$root" \
      --log "$root/${row_id}.log" \
      --command-file "$root/commands/${row_id}.txt" \
      --input-file "$INPUT" \
      --output-file "$final_csv" \
      --output-glob "$snapshot_dir/*" \
      --output-class "${RUN_OUTPUT_CLASS:-candidate}" \
      --notes "same-gamma shock-density-bubble explicit row" \
      -- "$EXE" "$INPUT" \
      max_step=200000 stop_time=0.3 amr.n_cell="$nx $ny" amr.max_grid_size=32 amr.plot_int=-1 \
      geometry.prob_lo="0 0" geometry.prob_hi="2 0.5" geometry.is_periodic="0 0" \
      euler.problem=shock_density_bubble_2d euler.method=explicit euler.spatial_order=2 \
      euler.slope_limiter=minmod euler.riemann="$riemann" euler.cfl=0.45 \
      euler.final_csv="$final_csv" \
      euler.shock_density_bubble_snapshot_dir="$snapshot_dir"
  }

  run_shock_imex() {
    local nx=$1
    local ny=$2
    local root=$3
    local prefix=imex_t1s2_bdltv20
    local row_id="${prefix}_${nx}x${ny}"
    local final_csv="$root/${row_id}_summary.csv"
    local snapshot_dir="$root/${row_id}_snapshots"
    mkdir -p "$root"
    python3 scripts/run_manifest.py run \
      --root "$ROOT" \
      --row-id "$row_id" \
      --output-dir "$root" \
      --log "$root/${row_id}.log" \
      --command-file "$root/commands/${row_id}.txt" \
      --input-file "$INPUT" \
      --output-file "$final_csv" \
      --output-glob "$snapshot_dir/*" \
      --output-class "${RUN_OUTPUT_CLASS:-candidate}" \
      --notes "same-gamma shock-density-bubble IMEX row" \
      -- "$EXE" "$INPUT" \
      max_step=200000 stop_time=0.3 amr.n_cell="$nx $ny" amr.max_grid_size=32 amr.plot_int=-1 \
      geometry.prob_lo="0 0" geometry.prob_hi="2 0.5" geometry.is_periodic="0 0" \
      euler.problem=shock_density_bubble_2d euler.method=imex euler.spatial_order=2 \
      euler.imex_form=bdltv20_t1_s2_source_map_picard euler.imex_cfl=0.45 \
      euler.imex_picard_iterations=4 euler.imex_solver_tol=1e-8 euler.imex_solver_max_iter=200 \
      euler.imex_acoustic_startup=0 euler.imex_acoustic_cfl_cap=0.0 \
      euler.imex_pressure_stabilization=off euler.imex_predictor_dissipation=material \
      euler.bdltv20_paper_pressure_solver=gmres \
      euler.final_csv="$final_csv" \
      euler.shock_density_bubble_snapshot_dir="$snapshot_dir"
  }

  BASE_EFF=results/amrex/same_gamma_shock_density_bubble_efficiency_divisible_2026-05-20
  BASE_SMOKE=results/amrex/project_same_gamma_shock_density_bubble_smoke_2026-05-20
  BASE_REPORT=results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20
  BASE_HIRES=results/amrex/project_same_gamma_shock_density_bubble_hires_reference_2026-05-20

  for root_n in "$BASE_EFF 128 32" "$BASE_SMOKE 160 40" "$BASE_REPORT 320 80"; do
    set -- $root_n
    root=$1
    nx=$2
    ny=$3
    run_shock_explicit "$nx" "$ny" "$root" hllc hllc_o2
    run_shock_explicit "$nx" "$ny" "$root" xie_am_hllc_p lowmach_corrected_hllc_p_o2
    run_shock_imex "$nx" "$ny" "$root"
  done

  # Long high-resolution reference row used by the shock-density-bubble checks.
  run_shock_explicit 640 160 "$BASE_HIRES" hllc hllc_o2

  python3 scripts/render_shock_density_bubble.py
  python3 scripts/summarise_shock_density_bubble_grid_consistency.py
  python3 scripts/compare_shock_density_bubble_hllc_reference_ladder.py
  python3 scripts/compare_shock_density_bubble_hllc640_reference.py
  python3 scripts/compare_shock_density_bubble_efficiency.py
}

case "${1:-help}" in
  build) build ;;
  riemann) riemann ;;
  gresho) gresho ;;
  advection) advection ;;
  shock_density_bubble) shock_density_bubble ;;
  all) build; riemann; gresho; advection; shock_density_bubble ;;
  *)
    echo "Usage: $0 {build|riemann|gresho|advection|shock_density_bubble|all}"
    ;;
esac
