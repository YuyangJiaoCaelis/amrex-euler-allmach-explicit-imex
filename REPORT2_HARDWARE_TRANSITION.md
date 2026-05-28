# Report 2 Hardware-Comparison Transition

The project description requires accuracy-versus-time comparison on different
hardware or parallelisation approaches, with MPI expected where possible and
CUDA/GPU comparison as an important target. This package is therefore a
numerical-method baseline and a hardware-transition base, not a completed
CPU/GPU performance package.

## Current Status

| Area | Current state | Report-2 implication |
|---|---|---|
| AMReX data layout | Cell-centred `MultiFab` state, AMReX geometry, ghost cells, reductions, and plotfile support. | Suitable base for MPI/CUDA-oriented AMReX work. |
| Explicit schemes | Reconstruction, boundary fills, flux updates, timesteps, and diagnostics use AMReX `MFIter`, `ParallelFor`, and reduction patterns in the report-facing path. A local 2-rank MPI advection-blob smoke check passed after making report CSV output gather rank-local rows. | Extend serial/MPI agreement checks across the retained explicit test families before using MPI timing claims. |
| IMEX pressure solve | T1/S2 BDLTV20 pressure assembly and GMRES solve are isolated under `src/euler_compare/imex/`, using a host Eigen sparse solve for the current serial CPU evidence. | Needs a reviewed MPI/GPU pressure-solve path, or a clearly scoped CPU-only IMEX comparison, before hardware-efficiency claims. |
| Output/analysis | CSV, snapshots, and plotting are separated from solver kernels. Per-cell report CSV/snapshot output gathers MPI rank-local rows to the I/O processor. | Timing can later exclude or separately report output cost. |
| Test families | Riemann, Gresho, periodic advection blob, and same-gamma shock-density-bubble are retained. | These are the starting test set for CPU/GPU scaling and accuracy/cost comparisons. |

## Required Steps Before CPU/GPU Claims

1. Confirm serial/MPI numerical agreement for the retained explicit test
   families, beyond the current advection-blob smoke row.
2. Build explicit rows with `USE_CUDA=TRUE` on available GPU hardware and
   confirm numerical agreement against the serial CPU baseline.
3. Profile explicit kernels with AMReX TinyProfiler and, on GPU, NVIDIA Nsight
   or an equivalent CUDA profiler.
4. Decide the IMEX hardware route:
   - keep IMEX as a serial or multi-core CPU baseline and compare against
     explicit CPU/GPU results with that limitation clearly stated; or
   - replace the host Eigen pressure solve with a reviewed AMReX/MLMG or other
     parallel sparse pressure-solve path and validate it before timing claims.
5. Only after numerical agreement is established, produce accuracy-versus-time
   plots separated by hardware backend.

## Latest Local Check

On 2026-05-28, the local macOS ARM64 environment built:

- serial CPU: `make -j2 DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE`;
- MPI CPU: `make -j2 DIM=2 USE_MPI=TRUE USE_OMP=FALSE USE_CUDA=FALSE`.

A 32x32 periodic advection-blob explicit HLLC smoke row was run both serially
and with `mpirun -np 2`. The serial and MPI `euler.final_csv` outputs both had
1024 rows and matched with maximum absolute numeric difference 0. This is a
launch and output-consistency smoke check, not performance evidence.

CUDA was not tested locally because `nvcc`, `nvidia-smi`, and Nsight tools were
not available on this machine. A dry-run make parse reached the expected CUDA
toolchain calls and reported missing CUDA tools.

## Architecture Rule For Future Work

Keep hardware-sensitive code in the existing domains:

- `src/euler_compare/numerics/` for device-safe explicit kernels;
- `src/euler_compare/problems/` for boundary and initial-condition kernels;
- `src/euler_compare/imex/` for pressure-solve implementation choices;
- `src/euler_compare/io/` for output kept out of kernel timing.

Do not mix plotting, report packaging, or exploratory solver experiments into
the report-facing numerical kernels.
