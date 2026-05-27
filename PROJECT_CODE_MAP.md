# MPhil project Code Map

The AMReX results in the report are produced by the application entry point:

`amrex/apps/euler_compare/src/main.cpp`

The entry point calls the AMReX application driver declared in:

`amrex/apps/euler_compare/include/euler_compare/application.hpp`

The implementation is compiled from:

`amrex/apps/euler_compare/src/euler_compare_app.cpp`

with internal numerical modules under:

`amrex/apps/euler_compare/src/euler_compare/`

For the high-level data flow, see `ARCHITECTURE.md`.

## Source Layout

| File | Main contents |
|---|---|
| `src/main.cpp` | Thin executable entry point. |
| `include/euler_compare/application.hpp` | Public declaration of the AMReX application driver. |
| `src/euler_compare_app.cpp` | AMReX/Eigen includes, conserved-component indices, and ordered inclusion of internal solver modules. |
| `src/euler_compare/core/` | Runtime configuration, state algebra, diagnostic records, option parsing, scheme selectors, and report-run consistency checks. |
| `src/euler_compare/problems/` | Exact Riemann states, Gresho/advection/shock-density-bubble initial conditions, and physical ghost-cell fills. |
| `src/euler_compare/numerics/` | Reconstruction, HLLC, Low-Mach Corrected HLLC-P, explicit finite-volume update, timestep rules, and error diagnostics. |
| `src/euler_compare/imex/` | BDLTV20/Toro-Vazquez source-map states, sparse pressure-row support, host GMRES solve, momentum correction, and energy closeout. |
| `src/euler_compare/io/` | Plotfile, CSV, error, snapshot, and MPhil project output writers. |
| `src/euler_compare/app/` | AMReX setup, time loop, timestep selection, final status output. |

## Report-Facing Runtime Selectors

| Report scheme name | Runtime selector |
|---|---|
| Explicit O2 HLLC | `euler.method=explicit`, `euler.spatial_order=2`, `euler.riemann=hllc` |
| Explicit O2 Low-Mach Corrected HLLC-P | `euler.method=explicit`, `euler.spatial_order=2`, `euler.riemann=xie_am_hllc_p` |
| IMEX T1/S2 BDLTV20 | `euler.method=imex`, `euler.spatial_order=2`, `euler.imex_form=bdltv20_t1_s2_source_map_picard` |

## Report-Facing Problem Selectors

| Report evidence | Runtime selector |
|---|---|
| Exact Riemann problems | `euler.problem=toro1` with `prob.test=1` (Sod), `93` (Lax/RP2), `3`, or `4` |
| Gresho vortex | `euler.problem=gresho_vortex` |
| Periodic advection blob | `euler.problem=advection_blob` |
| Same-gamma shock-density-bubble | `euler.problem=shock_density_bubble_2d` |

## Figure-Generation Scripts

| Evidence group | Scripts |
|---|---|
| Exact Riemann | `scripts/run_riemann_1d_convergence.py`, `scripts/plot_riemann_1d_overlays.py`, `scripts/plot_riemann_1d_efficiency.py`, `scripts/plot_riemann_1d_refinement_zoom.py` |
| Gresho vortex | `scripts/run_gresho_reference_matrix.py`, `scripts/combine_gresho_reproduction_data.py`, `scripts/plot_gresho_context_figures.py`, `scripts/plot_gresho_density_error.py` |
| Periodic advection blob | `scripts/run_advection_blob_efficiency.py`, `scripts/plot_advection_blob_efficiency.py`, `scripts/plot_advection_blob_timesnaps.py` |
| Same-gamma shock-density-bubble | `scripts/render_shock_density_bubble.py`, `scripts/summarise_shock_density_bubble_grid_consistency.py`, `scripts/compare_shock_density_bubble_hllc_reference_ladder.py`, `scripts/compare_shock_density_bubble_hllc640_reference.py`, `scripts/compare_shock_density_bubble_efficiency.py` |

The reproduction commands in `PROJECT_PLOT_REPRODUCTION.md` and `reproduce_project_data.sh` call the project-facing paths listed above.

For a compact post-build check before longer runs, use:

```bash
amrex/apps/euler_compare/tests/smoke/run_smoke_matrix.sh ../euler_compare_smoke
```
