# AMReX Euler All-Mach Explicit/IMEX

This repository contains a two-dimensional AMReX finite-volume code for comparing explicit and IMEX schemes for the compressible Euler equations. It is organised as a reproducible scientific-computing codebase: the solver sources, test cases, plotting scripts, and figure-regeneration commands are kept together.

## Contents

| Path | Contents |
|---|---|
| `amrex/apps/euler_compare/src/main.cpp`, `src/euler_compare_app.cpp`, `src/euler_compare/`, and `include/euler_compare/application.hpp` | Project-written AMReX Euler comparison driver. The public `main.cpp` is a thin entry point, `euler_compare_app.cpp` owns the AMReX application build unit, and the internal modules separate configuration, test problems, explicit HLLC and Low-Mach Corrected HLLC-P updates, IMEX T1/S2 BDLTV20 pressure-split updates, outputs, and the main time loop. |
| `amrex/apps/euler_compare/GNUmakefile`, `Make.package`, `inputs-ci` | Project build and smoke-test files following the AMReX online tutorial/single-level application layout. |
| `amrex/apps/euler_compare/tests/` | Smoke-test entry points for checking the reported schemes before longer runs. |
| `amrex/apps/euler_compare/benchmarks/` | Benchmark/reproduction entry points for longer report-data regeneration. |
| `scripts/` | Project-written Python scripts that run the MPhil project evidence rows and regenerate the report figures. |
| `external/AMReX/` | Third-party AMReX framework source needed for the local GNUmake build. |
| `external/eigen3/` | Third-party Eigen headers used by the host sparse linear algebra and GMRES pressure solve. |
| `PROJECT_PLOT_REPRODUCTION.md` | Command map from each project figure to the scripts and AMReX rows that produce it. |
| `PROJECT_CODE_MAP.md` | Short guide to the project-facing regions of the AMReX driver. |
| `ARCHITECTURE.md` | Higher-level source layout and data-flow guide for the AMReX implementation. |
| `REPORT2_HARDWARE_TRANSITION.md` | CPU/MPI/CUDA transition plan and current hardware-readiness boundary. |
| `SOURCE_BOUNDARY.md` | Short statement of the schemes, tests, and helper modules retained in this source tree. |
| `scripts/run_project_smoke_matrix.py` | Compact build-check matrix for the three reported schemes on the four retained test families. |
| `THIRD_PARTY_CODE.md` | Authorship and external-code notes. |
| `CHANGES_FROM_TUTORIAL.md` | Short provenance note separating AMReX tutorial layout conventions from project-specific solver and analysis code. |
| `README_PROJECT_SETTINGS.md` | Compact list of report scheme selectors, boundary conditions, final times, and output roots. |
| `verification/` | Short record of the latest completed build and reproduction checks. Long logs and generated figures are excluded from the source tree. |

The project figures use the three scheme selectors below. For code reading, start with `ARCHITECTURE.md` and `PROJECT_CODE_MAP.md`; they point to the source files used by the reproduced figures.

| Report name | Code selector |
|---|---|
| Explicit O2 HLLC | `euler.method=explicit`, `euler.spatial_order=2`, `euler.riemann=hllc` |
| Explicit O2 Low-Mach Corrected HLLC-P | `euler.method=explicit`, `euler.spatial_order=2`, `euler.riemann=xie_am_hllc_p` |
| IMEX T1/S2 BDLTV20 | `euler.method=imex`, `euler.imex_form=bdltv20_t1_s2_source_map_picard`, `euler.spatial_order=2` |

The shock-density-bubble evidence uses a single-material same-gamma Cartesian case. The AMReX driver is compiled as one application, so the supplied GNUmake build is the same build used for all reproduced runs. The selectable IMEX path is the BDLTV20 T1/S2 source-map pressure-split scheme.

## Build

From the package root:

```bash
cd amrex/apps/euler_compare
make -j2 DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE
```

The executable is:

```bash
amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex
```

Compact solver smoke matrix:

```bash
amrex/apps/euler_compare/tests/smoke/run_smoke_matrix.sh ../euler_compare_smoke
```

This command writes fresh smoke logs outside the source tree unless an explicit
output path inside the package is supplied.

## Python Requirements

The plotting scripts use Python 3 with:

```bash
python3 -m pip install -r requirements.txt
```

## Reproducing Figures

The full figure-generation commands are listed in `PROJECT_PLOT_REPRODUCTION.md`. The wrapper below provides the same command groups:

```bash
./reproduce_project_data.sh riemann
./reproduce_project_data.sh gresho
./reproduce_project_data.sh advection
./reproduce_project_data.sh shock_density_bubble
```

For a source-tree-clean full reproduction run, use the benchmark wrapper with an
external output directory:

```bash
amrex/apps/euler_compare/benchmarks/report_reproduction/run_report_reproduction.sh ../euler_compare_full_reproduction
```

Some rows, especially the 640x160 shock-density-bubble reference and high-resolution Gresho sweeps, can take a long time on a laptop.
These targets rerun AMReX and regenerate fresh outputs from the supplied source.
Measured serial wall times may change with hardware and current load.
