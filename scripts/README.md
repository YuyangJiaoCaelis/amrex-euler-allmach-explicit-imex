# MPhil project Scripts

These scripts run the AMReX rows and regenerate the figures used in MPhil project.

| Script | Purpose |
|---|---|
| `run_manifest.py` | Shared helper and command wrapper for writing JSON provenance manifests for AMReX run rows. |
| `run_mpi_explicit_agreement.py` | Runs serial-vs-MPI agreement smoke or candidate checks for retained explicit Report 2 rows. |
| `run_riemann_1d_convergence.py` | Runs the Riemann refinement rows for the three report schemes. |
| `riemann_exact.py` | Exact Riemann solver imported by the Riemann runner and refresh scripts. |
| `plot_riemann_1d_overlays.py` | Builds the exact-Riemann overlay and IMEX Toro 3/4 refinement figures from the Riemann package. |
| `plot_riemann_1d_efficiency.py` | Builds Riemann error-versus-wall-time figures with grid labels. |
| `plot_riemann_1d_refinement_zoom.py` | Builds the Toro 3/4 IMEX refinement and pointwise-error figure. |
| `run_gresho_reference_matrix.py` | Runs the Gresho vortex rows with exact Dirichlet field boundaries and \(t=0.4\pi\). |
| `combine_gresho_reproduction_data.py` | Combines the three Gresho summary CSV files into the source CSV used by the Gresho plots. |
| `plot_gresho_context_figures.py` | Builds the Gresho field snapshot and target-error cost bar chart. |
| `plot_gresho_density_error.py` | Regenerates the polished Gresho pressure, velocity, and density error figures from the Gresho source CSV. |
| `run_advection_blob_efficiency.py` | Runs the periodic advection-blob rows and builds error-versus-wall-time plots. |
| `plot_advection_blob_efficiency.py` | Builds the advection-blob density error-versus-wall-time and error-versus-grid report figures. |
| `plot_advection_blob_timesnaps.py` | Builds advection-blob density snapshots and density-error history. |
| `render_shock_density_bubble.py` | Builds the same-gamma shock-density-bubble density contact sheet and animation from AMReX snapshots. |
| `summarise_shock_density_bubble_grid_consistency.py` | Summarises grid/domain consistency for the same-gamma shock-density-bubble snapshots. |
| `compare_shock_density_bubble_hllc_reference_ladder.py` | Compares HLLC grid-reference behaviour for the same-gamma shock-density-bubble. |
| `compare_shock_density_bubble_hllc640_reference.py` | Compares the 320x80 rows with the averaged 640x160 HLLC numerical reference. |
| `compare_shock_density_bubble_efficiency.py` | Builds final-time shock-density-bubble error-versus-wall-time curves. |
