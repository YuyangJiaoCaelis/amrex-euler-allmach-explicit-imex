# Authorship And External Code

## Project-Written Code

- `amrex/apps/euler_compare/src/main.cpp`
- `amrex/apps/euler_compare/src/euler_compare_app.cpp`
- `amrex/apps/euler_compare/src/euler_compare/**/*.hpp`
- `amrex/apps/euler_compare/include/euler_compare/application.hpp`
- `amrex/apps/euler_compare/GNUmakefile`
- `amrex/apps/euler_compare/Make.package`
- `amrex/apps/euler_compare/inputs-ci`
- `scripts/*.py`
- `reproduce_project_data.sh`

The AMReX application layout and makefile structure follow the public AMReX tutorial style for a single-level application. The Euler schemes, MPhil project test cases, metrics, and plotting scripts were written for this project.

The preserved source copy excludes earlier development-only shock-bubble transfer files. The retained source is the AMReX application and reproduction scripts used for the reported schemes and tests.

## External Code

| Path | Role |
|---|---|
| `external/AMReX/` | AMReX framework source used for `MultiFab`, geometry, tiling, boundary handling, reductions, and the GNUmake build system. |
| `external/eigen3/` | Eigen header library used for sparse matrices and host linear solves in the IMEX pressure step. |

The external directories are included so the marker can build the supplied AMReX app without reconstructing the dependency layout.
