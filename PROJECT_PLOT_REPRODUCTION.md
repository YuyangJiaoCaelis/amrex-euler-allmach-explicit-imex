# Report Plot Reproduction

Build the executable before running the AMReX rows:

```bash
cd amrex/apps/euler_compare
make -j2 DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE
cd ../../..
```

The scripts create fresh outputs under `results/amrex/`,
`project_outputs/test_packages/`, and `project_outputs/project_evidence/`.

Several figures use measured serial wall time on the machine that runs the
scripts. Those figures should reproduce the same schemes, grids, error trends,
and ordering, but their exact pixel values can change between machines or
between repeated runs because the horizontal axis is newly timed. Field
snapshots, exact-solution overlays, and reference-error plots are deterministic
up to the usual plotting-library rendering details.

## Exact Riemann Figures

```bash
python3 scripts/run_riemann_1d_convergence.py
python3 scripts/plot_riemann_1d_overlays.py
python3 scripts/plot_riemann_1d_efficiency.py
python3 scripts/plot_riemann_1d_refinement_zoom.py
```

| Report figure | Generated file |
|---|---|
| Fig. 4.1 exact Riemann overlays | `project_outputs/test_packages/riemann_1d/figures/riemann_1d_primary_overlays_full_domain.png` |
| Fig. 4.2 IMEX Toro 3/4 refinement | `project_outputs/test_packages/riemann_1d/figures/riemann_1d_toro34_imex_refinement_wave_zoom.png` |
| Fig. 4.3 density error versus wall time | `project_outputs/test_packages/riemann_1d_convergence/figures/riemann_1d_density_error_vs_walltime_with_grid_labels.png` |

These rows use exact Dirichlet data in both coordinate directions for the one-dimensional Riemann profile embedded in a 2D AMReX box.
The Riemann runner also generates 800-cell rows for the Toro 3/4 refinement
zoom in Fig. 4.2. The report-facing efficiency plots use grids up to 400 cells
so the additional refinement row does not alter Fig. 4.3.

## Gresho Vortex Figures

The Gresho evidence uses exact Dirichlet field boundaries in both directions,
non-periodic geometry, and final time \(0.4\pi\) for every row. The low-Mach
explicit rows need a high step budget because the acoustic CFL is severe at
Mach \(0.001\).

```bash
python3 scripts/run_gresho_report_rows.py --row-timeout-sec 7200
```

Then combine the row summaries and regenerate the report figures with:

```bash
python3 scripts/combine_gresho_reproduction_data.py
python3 scripts/plot_gresho_context_figures.py
python3 scripts/plot_gresho_density_error.py
```

| Report figure | Generated file |
|---|---|
| Fig. 4.4 Gresho field snapshot | `project_outputs/test_packages/gresho_vortex/figures/amrex_gresho_field_snapshot_mach0p01_n48_t0p4pi.png` |
| Fig. 4.5 target-error cost bar chart | `project_outputs/test_packages/gresho_vortex/figures/amrex_gresho_target_error_cost_bar_primary.png` |
| Fig. 4.6 pressure error versus wall time | `project_outputs/report_figures_cued_template/amrex_gresho_pressure_error_vs_walltime_with_grid_labels_polished.png` |
| Fig. 4.7 velocity error versus wall time | `project_outputs/report_figures_cued_template/amrex_gresho_velocity_error_vs_walltime_with_grid_labels_polished.png` |
| Fig. 4.8 density error versus wall time | `project_outputs/report_figures_cued_template/amrex_gresho_density_error_vs_walltime_with_grid_labels_polished.png` |

The pressure and velocity polished copies are produced from the same combined Gresho CSV and stored with the LaTeX report figures.

## Periodic Advection-Blob Figures

```bash
python3 scripts/run_advection_blob_efficiency.py
python3 scripts/plot_advection_blob_efficiency.py
python3 scripts/plot_advection_blob_timesnaps.py
```

| Report figure | Generated file |
|---|---|
| Fig. 4.9 density snapshots | `project_outputs/test_packages/advection_blob/figures/amrex_advection_blob_density_timesnaps_current.png` |
| Fig. 4.10 density error versus time | `project_outputs/test_packages/advection_blob/figures/amrex_advection_blob_density_error_vs_time_current.png` |
| Fig. 4.11 density error versus wall time | `project_outputs/report_figures_cued_template/amrex_advection_blob_density_error_vs_walltime_with_grid_labels_polished.png` |
| Fig. 4.12 density error versus grid size | `project_outputs/test_packages/advection_blob/figures/amrex_advection_blob_density_error_vs_grid_primary.png` |

The advection-blob test uses periodic boundaries, matching the test description in the project brief.
The IMEX density-grid evidence uses the report-facing
`euler.bdltv20_paper_t1_s2=advection_blob_periodic_2d` path with
CFL-controlled time stepping, four Picard iterations, and solver tolerance
`1e-8`.

## Same-Gamma Shock-Density-Bubble Figures

The case is a single-material Cartesian ideal-gas high-speed test with \(\gamma=1.4\), fixed post-shock inflow on the left, outflow on the right/top, and reflective symmetry on the lower boundary.

The wrapper `./reproduce_project_data.sh shock_density_bubble` runs the full set of rows needed for Figs. 4.13-4.16: 128x32, 160x40, and 320x80 rows for all three schemes, plus the 640x160 HLLC reference. The 640x160 row is the long reference run.

Representative 320x80 report-grid commands:

```bash
./amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex amrex/apps/euler_compare/inputs-ci \
  max_step=200000 stop_time=0.3 amr.n_cell="320 80" amr.max_grid_size=32 amr.plot_int=-1 \
  geometry.prob_lo="0 0" geometry.prob_hi="2 0.5" geometry.is_periodic="0 0" \
  euler.problem=shock_density_bubble_2d euler.method=explicit euler.spatial_order=2 \
  euler.slope_limiter=minmod euler.riemann=hllc euler.cfl=0.45 \
  euler.final_csv=results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20/hllc_o2_320x80_summary.csv \
  euler.shock_density_bubble_snapshot_dir=results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20/hllc_o2_320x80_snapshots

./amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex amrex/apps/euler_compare/inputs-ci \
  max_step=200000 stop_time=0.3 amr.n_cell="320 80" amr.max_grid_size=32 amr.plot_int=-1 \
  geometry.prob_lo="0 0" geometry.prob_hi="2 0.5" geometry.is_periodic="0 0" \
  euler.problem=shock_density_bubble_2d euler.method=explicit euler.spatial_order=2 \
  euler.slope_limiter=minmod euler.riemann=xie_am_hllc_p euler.cfl=0.45 \
  euler.final_csv=results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20/lowmach_corrected_hllc_p_o2_320x80_summary.csv \
  euler.shock_density_bubble_snapshot_dir=results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20/lowmach_corrected_hllc_p_o2_320x80_snapshots

./amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex amrex/apps/euler_compare/inputs-ci \
  max_step=200000 stop_time=0.3 amr.n_cell="320 80" amr.max_grid_size=32 amr.plot_int=-1 \
  geometry.prob_lo="0 0" geometry.prob_hi="2 0.5" geometry.is_periodic="0 0" \
  euler.problem=shock_density_bubble_2d euler.method=imex euler.spatial_order=2 \
  euler.imex_form=bdltv20_t1_s2_source_map_picard euler.imex_cfl=0.45 \
  euler.imex_picard_iterations=4 euler.imex_solver_tol=1e-8 euler.imex_solver_max_iter=200 \
  euler.imex_acoustic_startup=0 euler.imex_acoustic_cfl_cap=0.0 \
  euler.final_csv=results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20/imex_t1s2_bdltv20_320x80_summary.csv \
  euler.shock_density_bubble_snapshot_dir=results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20/imex_t1s2_bdltv20_320x80_snapshots
```

Then run:

```bash
python3 scripts/render_shock_density_bubble.py
python3 scripts/summarise_shock_density_bubble_grid_consistency.py
python3 scripts/compare_shock_density_bubble_hllc_reference_ladder.py
python3 scripts/compare_shock_density_bubble_hllc640_reference.py
python3 scripts/compare_shock_density_bubble_efficiency.py
```

| Report figure | Generated file |
|---|---|
| Fig. 4.13 density snapshots | `project_outputs/test_packages/shock_density_bubble_cartesian_same_gamma/figures/shock_density_bubble_cartesian_three_scheme_density_snapshots_320x80.png` |
| Fig. 4.14 HLLC grid-reference ladder | `project_outputs/test_packages/shock_density_bubble_cartesian_same_gamma/figures/shock_density_bubble_cartesian_hllc_grid_reference_ladder.png` |
| Fig. 4.15 errors against 640x160 HLLC reference | `project_outputs/test_packages/shock_density_bubble_cartesian_same_gamma/figures/shock_density_bubble_cartesian_hllc640_reference_errors.png` |
| Fig. 4.16 shock-density-bubble error versus wall time | `project_outputs/test_packages/shock_density_bubble_cartesian_same_gamma/figures/shock_density_bubble_cartesian_hllc640_efficiency.png` |
