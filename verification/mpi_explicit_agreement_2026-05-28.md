# MPI Explicit Agreement Check 2026-05-28

Branch: `exp/mpi-explicit-agreement`

Purpose: establish a compact local serial-vs-2-rank-MPI smoke gate for the
retained explicit Report 2 families before moving toward hardware timing.

## Change Under Test

- Added `scripts/run_mpi_explicit_agreement.py`.
- Fixed `write_shock_density_bubble_snapshot` so all MPI ranks participate in
  the gathered snapshot CSV write. Previously non-I/O ranks returned before the
  collective gather, causing two-rank shock snapshot rows to stall.

## Commands

```sh
cd /Users/yuyangjiao/Desktop/MPhilresearch/amrex_euler_allmach_explicit_imex/amrex/apps/euler_compare
make -j2 DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE
make -j2 DIM=2 USE_MPI=TRUE USE_OMP=FALSE USE_CUDA=FALSE

cd /Users/yuyangjiao/Desktop/MPhilresearch/amrex_euler_allmach_explicit_imex
python3 scripts/run_mpi_explicit_agreement.py \
  --output /tmp/report2_mpi_explicit_agreement_2026-05-28 \
  --mpi-ranks 2 \
  --tolerance 1e-10 \
  --row-timeout-sec 20 \
  --force
```

## Result

The runner completed successfully with summary:

- Riemann Sod explicit HLLC and low-Mach HLLC-P: pass.
- Gresho Mach 0.01 explicit HLLC and low-Mach HLLC-P: pass.
- Periodic advection blob explicit HLLC and low-Mach HLLC-P: pass.
- Shock-density-bubble explicit HLLC and low-Mach HLLC-P: pass.

For all eight rows, serial and 2-rank MPI outputs had matching row counts and
maximum numeric CSV difference `0` at tolerance `1e-10`.

Summary artifact:
`/tmp/report2_mpi_explicit_agreement_2026-05-28/summary.csv`

## Claim Boundary

This is a launch/output/numerical-agreement smoke check on small local grids,
not performance evidence and not a CUDA/GPU validation. MPI timing claims still
need candidate or frozen runs with clean provenance, larger report-relevant
grids, and timing methodology that separates solver cost from output cost.
