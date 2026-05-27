# Final Architecture Pass - 2026-05-26

## Scope

This pass prepared the future-development AMReX Euler package for long-term
research use before any further full reproduction run. The frozen Report-1
snapshot was left untouched.

## Structural Changes

- Replaced the former application-sized entry point with a thin executable
  entry point:
  `amrex/apps/euler_compare/src/main.cpp`.
- Added a public application declaration:
  `amrex/apps/euler_compare/include/euler_compare/application.hpp`.
- Moved the AMReX application implementation into:
  `amrex/apps/euler_compare/src/euler_compare_app.cpp`.
- Moved internal numerical modules into:
  `amrex/apps/euler_compare/src/euler_compare/`.
- Grouped internal modules into hardware-relevant domains:
  `core/`, `problems/`, `numerics/`, `imex/`, `io/`, and `app/`.
- Updated `Make.package` so AMReX GNUmake compiles both `main.cpp` and
  `euler_compare_app.cpp`.
- Updated documentation to describe the new `src/`, `include/`, `tests/`, and
  `benchmarks/` layout.
- Added `REPORT2_HARDWARE_TRANSITION.md` to make clear that the project
  ultimately requires CPU/MPI/CUDA comparison and that the current host Eigen
  IMEX pressure solve is not yet a GPU/MPI pressure-solve implementation.
- Updated the smoke wrapper so smoke outputs go to a caller-supplied or
  timestamped external directory, including AMReX plotfiles generated during
  smoke verification.
- Updated the full reproduction wrapper so a full run can write `results/` and
  `project_outputs/` through temporary links into an external output root.
- Archived the interrupted partial reproduction outputs outside this source
  package at:
  `/Users/yuyangjiao/Desktop/MPhilresearch/research_codebase_interrupted_outputs_2026-05-26`.

## Verification

Build:

```bash
cd amrex/apps/euler_compare
make -j2 DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE
```

Result: passed.

Compact smoke:

```bash
amrex/apps/euler_compare/tests/smoke/run_smoke_matrix.sh \
  /Users/yuyangjiao/Desktop/MPhilresearch/research_codebase_verification_2026-05-26/hardware_architecture_smoke
```

Result: all 12 rows passed.

Rows covered:

- Riemann Sod: Explicit O2 HLLC, Explicit O2 Low-Mach Corrected HLLC-P,
  IMEX T1/S2 BDLTV20.
- Gresho Mach 0.01 at 32x32: the same three schemes.
- Periodic advection blob at 32x32: the same three schemes.
- Same-gamma shock-density-bubble at 160x40: the same three schemes.

The smoke summary is:

`/Users/yuyangjiao/Desktop/MPhilresearch/research_codebase_verification_2026-05-26/hardware_architecture_smoke/summary.csv`

This is a serial CPU verification of the reorganised source tree. It does not
claim MPI, CUDA, GPU, or multi-core performance readiness.

Final source-tree hygiene after verification:

- no `results/`;
- no `project_outputs/`;
- no AMReX plotfile directories;
- no `tmp_build_dir`;
- no executable build product;
- no Python cache files.

## Freeze Decision

The future-development package is now suitable for the next full reproduction
run. Do not request another layout polish before that run unless one of these
specific problems appears:

- build failure;
- smoke-matrix failure;
- a reported scheme is missing or renamed incorrectly;
- a generated figure cannot be traced to a command and data file;
- earlier development-only shock-bubble paths re-enter the report-facing source
  path.

The next full reproduction should use an external output root, for example:

```bash
amrex/apps/euler_compare/benchmarks/report_reproduction/run_report_reproduction.sh \
  /Users/yuyangjiao/Desktop/MPhilresearch/research_codebase_full_reproduction_2026-05-26
```
