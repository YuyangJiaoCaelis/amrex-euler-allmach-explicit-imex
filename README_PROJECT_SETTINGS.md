# Report Settings

The report settings are defined by the scripts and commands in `PROJECT_PLOT_REPRODUCTION.md` and `reproduce_project_data.sh`. The `inputs-ci` file is a compact smoke-test template; the report rows set grid size, boundary condition, final time, scheme, solver tolerance, and output path through command-line overrides.

## Scheme Names

| Report name | Main runtime settings |
|---|---|
| Explicit O2 HLLC | `euler.method=explicit`, `euler.spatial_order=2`, `euler.riemann=hllc` |
| Explicit O2 Low-Mach Corrected HLLC-P | `euler.method=explicit`, `euler.spatial_order=2`, `euler.riemann=xie_am_hllc_p` |
| IMEX T1/S2 BDLTV20 | `euler.method=imex`, `euler.spatial_order=2`, `euler.imex_form=bdltv20_t1_s2_source_map_picard` |

## Test Settings

| Test group | Boundary and time settings |
|---|---|
| Riemann problems | Exact Dirichlet boundary data for the one-dimensional profile embedded in a 2D AMReX box; final times follow the Toro/Lax cases used by the Riemann runner. |
| Gresho vortex | Exact Dirichlet field boundary, Mach numbers `0.001`, `0.01`, and `0.1`, final time `0.4*pi`. |
| Periodic density-blob advection | Periodic boundaries, final time `0.5`, density blob transported through the domain. |
| Same-gamma shock-density-bubble | Cartesian ideal-gas case with `gamma=1.4`, fixed post-shock inflow on the left, outflow on the right/top, reflective lower boundary, final time `0.3`, snapshots at `0,0.06,0.12,0.18,0.24,0.30`. |

## Output Locations

Fresh reproduction runs write under:

| Output root | Contents |
|---|---|
| `results/amrex/` | Row summaries, logs, and snapshots from direct runs. |
| `project_outputs/test_packages/` | Regenerated report-package figures and source CSV files. |
| `project_outputs/report_figures_cued_template/` | Polished copies sized for the CUED LaTeX report template. |
| `project_outputs/project_evidence/` | Convenience copies of selected figure and CSV outputs. |

Large rows, especially high-resolution Gresho and shock-density-bubble reference rows, may require substantially longer wall time than the smoke rows.
