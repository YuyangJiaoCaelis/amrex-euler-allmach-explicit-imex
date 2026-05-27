# Source Boundary

This repository is the AMReX Euler comparison source. It keeps the source needed to build and reproduce the project figures for:

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

The enum types retain a few additional slope-limiter, Riemann-solver, and problem identifiers used by the solver configuration system. The reproduction scripts use the three scheme selectors and four test families listed above.

The main reproducible code paths are listed in `PROJECT_PLOT_REPRODUCTION.md`.
