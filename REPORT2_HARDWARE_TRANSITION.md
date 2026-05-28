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
| Explicit schemes | Reconstruction, boundary fills, flux updates, timesteps, and diagnostics use AMReX `MFIter`, `ParallelFor`, and reduction patterns in the report-facing path. A local 2-rank MPI smoke matrix now passes for the retained explicit Riemann, Gresho, advection-blob, and shock-density-bubble rows after making report CSV and shock snapshot output gather rank-local rows. | Promote the smoke gate to candidate/frozen report grids before using MPI timing claims. |
| IMEX pressure solve | T1/S2 BDLTV20 pressure assembly and GMRES solve are isolated under `src/euler_compare/imex/`, using a host Eigen sparse solve for the current serial CPU evidence. | Needs a reviewed MPI/GPU pressure-solve path, or a clearly scoped CPU-only IMEX comparison, before hardware-efficiency claims. |
| Output/analysis | CSV, snapshots, and plotting are separated from solver kernels. Per-cell report CSV/snapshot output gathers MPI rank-local rows to the I/O processor. | Timing can later exclude or separately report output cost. |
| Test families | Riemann, Gresho, periodic advection blob, and same-gamma shock-density-bubble are retained. | These are the starting test set for CPU/GPU scaling and accuracy/cost comparisons. |

## Required Steps Before CPU/GPU Claims

1. Confirm serial/MPI numerical agreement on the report-scale explicit test
   rows, using the smoke matrix as the pre-timing gate.
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

The local explicit agreement smoke matrix was run with
`scripts/run_mpi_explicit_agreement.py` and `mpirun -np 2`. It covered Riemann
Sod, Gresho Mach 0.01, periodic advection blob, and shock-density-bubble rows
for both HLLC and low-Mach HLLC-P. All eight serial/MPI CSV comparisons matched
row counts and had maximum absolute numeric difference 0 at tolerance `1e-10`.
This is a launch, output-consistency, and small-grid numerical-agreement smoke
check, not performance evidence.

The first candidate pre-timing matrix was also run on 2026-05-28 from clean
commit `188d1dff4c36e0f5779e15a36fda82ad0ee3d868`, using
`scripts/run_mpi_explicit_agreement.py --case-set candidate`. It covered
Riemann Sod `400x5`, Gresho Mach 0.01 `32x32`, periodic advection blob
`128x128`, and shock-density-bubble `160x40` rows for HLLC and low-Mach HLLC-P.
All eight serial/MPI CSV comparisons matched row counts and had maximum
absolute numeric difference 0 at tolerance `1e-10`; shock rows also matched
gathered snapshot CSVs. This is candidate numerical-agreement evidence, not MPI
scaling or performance evidence.

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
