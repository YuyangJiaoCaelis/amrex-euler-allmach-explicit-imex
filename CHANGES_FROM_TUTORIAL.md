# Changes From AMReX Tutorial Layout

This code package follows the single-level AMReX application style used by the online AMReX tutorials: GNUmake build files, parameter-file inputs, `MultiFab` state storage, box/tile loops, ghost-cell fills, and face-based finite-volume updates.

## AMReX Tutorial Elements Used

The AMReX framework source is supplied under `external/AMReX/`. The project application keeps the tutorial-style structure:

| Area | Role in this package |
|---|---|
| GNUmake application build | `amrex/apps/euler_compare/GNUmakefile`, `Make.package` |
| Runtime input pattern | `amrex/apps/euler_compare/inputs-ci` plus command-line parameter overrides |
| AMReX data layout | cell-centred `MultiFab` conserved state and ghost cells |
| Loop structure | AMReX box/tile iteration over cell and face ranges |
| Output pattern | CSV summaries and optional AMReX plotfiles/snapshots |

## Project Code Written For MPhil project

The MPhil project solver and analysis additions are project-specific:

| File or directory | Project-specific contents |
|---|---|
| `amrex/apps/euler_compare/src/main.cpp` | Thin executable entry point. |
| `amrex/apps/euler_compare/src/euler_compare_app.cpp` | AMReX application build unit that includes the internal numerical modules. |
| `amrex/apps/euler_compare/include/euler_compare/application.hpp` | Public declaration of the application driver. |
| `amrex/apps/euler_compare/src/euler_compare/core/` | Runtime configuration, state algebra, diagnostics, scheme selectors, and report-run consistency checks. |
| `amrex/apps/euler_compare/src/euler_compare/problems/` | Exact Riemann, Gresho vortex, periodic density-blob advection, and same-gamma shock-density-bubble setup and boundary handling. |
| `amrex/apps/euler_compare/src/euler_compare/numerics/` | Explicit finite-volume update, MUSCL/minmod reconstruction, HLLC, and Low-Mach Corrected HLLC-P. |
| `amrex/apps/euler_compare/src/euler_compare/imex/` | IMEX T1/S2 BDLTV20 pressure-split path and pressure solve. |
| `amrex/apps/euler_compare/src/euler_compare/io/` | CSV, snapshot, and project-case output routines. |
| `amrex/apps/euler_compare/src/euler_compare/app/` | AMReX driver loop and final status output. |
| `scripts/riemann_exact.py` | Exact Riemann reference solver used for one-dimensional overlays and error norms. |
| `scripts/run_riemann_1d_convergence.py` and refresh scripts | Riemann refinement rows and report figures. |
| `scripts/run_gresho_reference_matrix.py` and Gresho plotting scripts | Gresho vortex rows, source-data combination, and report figures. |
| `scripts/run_advection_blob_efficiency.py` and advection plotting scripts | Periodic density-blob runs and report figures. |
| `scripts/render_shock_density_bubble.py` and shock-density-bubble comparison scripts | Same-gamma Cartesian high-speed test figures, reference checks, and efficiency plots. |

## Report Scheme Selectors

Report figures use these selectors:

| Report scheme name | Runtime selector |
|---|---|
| Explicit O2 HLLC | `euler.method=explicit`, `euler.spatial_order=2`, `euler.riemann=hllc` |
| Explicit O2 Low-Mach Corrected HLLC-P | `euler.method=explicit`, `euler.spatial_order=2`, `euler.riemann=xie_am_hllc_p` |
| IMEX T1/S2 BDLTV20 | `euler.method=imex`, `euler.spatial_order=2`, `euler.imex_form=bdltv20_t1_s2_source_map_picard` |

The figure-generation commands in `PROJECT_PLOT_REPRODUCTION.md` identify the paths used for the submitted MPhil project results. `PROJECT_CODE_MAP.md` gives a source-level guide to the corresponding AMReX driver regions.
