# MPI Candidate Agreement Check 2026-05-28

Branch: `codex/report2-next-research-step`

Commit: `188d1dff4c36e0f5779e15a36fda82ad0ee3d868`

Purpose: run the first Report 2 candidate pre-timing agreement matrix for
retained explicit rows before producing MPI timing evidence.

## Commands

```sh
cd /Users/yuyangjiao/Desktop/MPhilresearch/amrex_euler_allmach_explicit_imex/amrex/apps/euler_compare
make -j2 DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE
make -j2 DIM=2 USE_MPI=TRUE USE_OMP=FALSE USE_CUDA=FALSE

cd /Users/yuyangjiao/Desktop/MPhilresearch/amrex_euler_allmach_explicit_imex
python3 scripts/run_mpi_explicit_agreement.py \
  --case-set candidate \
  --output /tmp/report2_candidate_mpi_explicit_agreement_2026-05-28 \
  --mpi-ranks 2 \
  --tolerance 1e-10 \
  --row-timeout-sec 600 \
  --force
```

## Result

All eight candidate rows passed serial-vs-2-rank-MPI agreement:

- Riemann Sod, `400 x 5`, HLLC and low-Mach HLLC-P.
- Gresho vortex, Mach `0.01`, `32 x 32`, HLLC and low-Mach HLLC-P.
- Periodic advection blob, `128 x 128`, HLLC and low-Mach HLLC-P.
- Shock-density-bubble, `160 x 40`, HLLC and low-Mach HLLC-P.

For every row, serial and MPI status values were `ok`, row counts matched, and
maximum numeric CSV difference was `0` at tolerance `1e-10`. Shock rows also
compared two gathered snapshot CSV files per scheme with maximum numeric
difference `0`.

The first manifest checked recorded:

- branch: `codex/report2-next-research-step`;
- dirty: `false`;
- output class: `candidate`.

Summary artifact:
`/tmp/report2_candidate_mpi_explicit_agreement_2026-05-28/summary.csv`

## Claim Boundary

This is candidate numerical agreement for explicit MPI rows. It is not MPI
scaling or performance evidence. Timing claims still require repeated rank
scans, machine metadata, and frozen evidence rows.
