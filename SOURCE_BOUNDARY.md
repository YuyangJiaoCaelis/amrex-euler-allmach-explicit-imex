# Source Boundary

This package is the report-facing AMReX Euler comparison source. It keeps the source needed to build and reproduce the submitted figures for:

- Explicit O2 HLLC.
- Explicit O2 Low-Mach Corrected HLLC-P.
- IMEX T1/S2 BDLTV20.

The retained AMReX tests are:

- one-dimensional exact Riemann problems embedded in the 2D AMReX grid;
- Gresho vortex;
- periodic density-blob advection;
- same-gamma Cartesian shock-density-bubble.

The public header `amrex/apps/euler_compare/include/euler_compare/application.hpp`
declares the application driver. The implementation modules under
`amrex/apps/euler_compare/src/euler_compare/` are grouped into `core/`,
`problems/`, `numerics/`, `imex/`, `io/`, and `app/`. This separation keeps
device-oriented explicit kernels, host-side IMEX pressure solves, physical
problem setup, and output routines visible as different parts of the codebase.

The enum types retain a few additional slope limiter, Riemann solver, and problem identifiers used while checking solver behaviour. The report reproduction scripts use only the three scheme selectors and four test families listed above.

Earlier two-material shock-bubble transfer files and standalone development checks have been removed from this preserved copy. The remaining code paths are the paths used by the report reproduction scripts in `PROJECT_PLOT_REPRODUCTION.md`.
