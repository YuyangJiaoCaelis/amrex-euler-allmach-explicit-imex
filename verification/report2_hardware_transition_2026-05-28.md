# Report 2 Hardware Transition Check - 2026-05-28

Status: current local hardware-transition check; no CUDA hardware was available
locally.

Scope: local build and MPI smoke check on macOS ARM64. This is not performance
evidence.

## Commands Checked

- `make -j2 DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE`
- `make -j2 DIM=2 USE_MPI=TRUE USE_OMP=FALSE USE_CUDA=FALSE`
- `make -n DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=TRUE`
- serial-vs-`mpirun -np 2` periodic advection-blob explicit HLLC smoke row

## Result

- Serial CPU build passed.
- MPI CPU build passed with Homebrew Open MPI.
- The first MPI smoke exposed incomplete custom `euler.final_csv` output:
  512 MPI rows versus 1024 serial rows for a 32x32 grid.
- After changing per-cell report CSV/snapshot output to gather rank-local rows
  to the AMReX I/O processor, the same smoke produced 1024 rows in both serial
  and MPI outputs.
- The post-fix serial/MPI comparison had maximum absolute numeric difference 0
  across common numeric CSV columns.
- CUDA was not built locally because `nvcc`, `nvidia-smi`, and Nsight tooling
  were unavailable.

## Remaining Boundary

The code is not yet ready for Report 2 hardware-efficiency claims. The next
required checks are serial/MPI agreement across the retained explicit test
families, a real CUDA build on GPU hardware, and an explicit decision on whether
IMEX remains CPU-only or receives a parallel pressure-solve implementation.
