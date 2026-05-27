# AMReX Euler Comparison Architecture

This package keeps the AMReX tutorial-style GNUmake application layout, while
separating the project implementation into a thin executable entry point, one
AMReX application build unit, internal numerical modules, smoke tests, and
longer reproduction benchmarks.

The public executable entry point is
`amrex/apps/euler_compare/src/main.cpp`. It only calls the application driver
declared in `amrex/apps/euler_compare/include/euler_compare/application.hpp`.
The AMReX driver and solver implementation are compiled from
`amrex/apps/euler_compare/src/euler_compare_app.cpp`, which includes internal
modules from `amrex/apps/euler_compare/src/euler_compare/`. This keeps the
supplied GNUmake build simple while making the source tree readable as a normal
scientific-computing application. The internal modules are grouped by role:

1. `core/solver_types.hpp`, `core/diagnostics.hpp`,
   `core/basic_diagnostics.hpp`,
   `imex_diagnostics.hpp`, `problem_metrics.hpp`, `run_config.hpp`,
   `config_parse_helpers.hpp`, and `config_reader.hpp` define the state types,
   diagnostics, runtime options, parsers, and consistency checks.
2. `problems/` contains analytic problem states, exact Riemann sampling,
   initial conditions, boundary fills, and the problem setup include hub.
3. `numerics/` contains reconstruction, explicit finite-volume fluxes,
   positivity/accuracy diagnostics, and timestep calculations.
4. `imex/` contains the BDLTV20/Toro-Vazquez pressure-split IMEX T1/S2 update,
   including host sparse pressure-system assembly and the Eigen GMRES solve
   path.
5. `io/` writes CSV rows, profile data, snapshots, and plotfiles.
6. `app/driver_main.hpp` owns AMReX initialisation, timestep selection,
   dispatch, and final status reporting.

## Main Data Flow

For explicit rows, cell-centred conserved variables are filled with the selected
test state, ghost cells are filled from the matching physical boundary rule,
face states are reconstructed with the selected slope limiter, HLLC or
Low-Mach Corrected HLLC-P fluxes are evaluated on faces, and flux differences
advance the cell state.

For IMEX T1/S2 rows, the explicit Toro-Vazquez material/source-map stage first
forms the star state. The pressure-split stage assembles a sparse pressure row
from the star state and lagged coefficients, solves it with Eigen GMRES, applies
the pressure-gradient momentum correction, and closes energy with the ideal-gas
EOS.

## CPU/GPU Transition Boundary

The project description requires later CPU/GPU-oriented efficiency comparison.
This package is therefore organised so the hardware-sensitive parts are visible:

- explicit finite-volume work is concentrated in `numerics/` and uses AMReX
  `ParallelFor`/reduction patterns that can be compiled for device execution;
- problem setup and boundary fills are kept in `problems/`, where device-safe
  ghost filling can be checked separately from numerical flux changes;
- the current IMEX pressure solve is deliberately isolated in `imex/` and uses a
  host Eigen GMRES path. This is suitable for the present serial CPU evidence
  but is not yet the final MPI/CUDA pressure-solve implementation;
- output and report plotting stay outside the solver kernels, so timing and
  profiling can later distinguish numerical work from I/O.

Before claiming Report-2 CPU/GPU efficiency, the project still needs separate
MPI and CUDA verification: explicit kernels should be built and profiled with
`USE_CUDA=TRUE`, and the IMEX pressure solve needs either a documented CPU
baseline comparison or a reviewed AMReX/MLMG or other parallel pressure-solve
path. The current layout is intended to make those next steps clean rather than
to claim they are already complete.

## Verification Entry Points

Build:

```bash
cd amrex/apps/euler_compare
make -j2 DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE
```

Figure reproduction:

```bash
./reproduce_project_data.sh riemann
./reproduce_project_data.sh gresho
./reproduce_project_data.sh advection
./reproduce_project_data.sh shock_density_bubble
```

Compact solver smoke matrix:

```bash
amrex/apps/euler_compare/tests/smoke/run_smoke_matrix.sh ../euler_compare_smoke
```

Concise reproduction records are stored under `verification/`. Fresh smoke logs
are written to the output path supplied to
`amrex/apps/euler_compare/tests/smoke/run_smoke_matrix.sh`.

For long reproduction runs, prefer an external output root:

```bash
amrex/apps/euler_compare/benchmarks/report_reproduction/run_report_reproduction.sh ../euler_compare_full_reproduction
```
