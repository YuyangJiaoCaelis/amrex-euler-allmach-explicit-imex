# Report 2 Candidate Evidence Matrix

This matrix defines the first candidate evidence rows for the Report 2
hardware-comparison stage. It is intentionally narrow: it turns the existing
small-grid MPI smoke gate into a controlled pre-timing agreement gate without
adding new schemes or exploratory diagnostics to the live solver.

Candidate outputs must be generated outside the source tree and must include
run manifests. A row may be promoted to frozen report evidence only after the
working tree is clean, the command is recorded, numerical agreement has passed,
and the claim boundary is written into the report text.

## Current Branch

- Active branch: `codex/report2-next-research-step`
- Local gate: `scripts/run_mpi_explicit_agreement.py --case-set candidate`
- Output class: `candidate`

## Candidate Agreement Rows

| Family | Grid and final time | Schemes | Current purpose |
|---|---|---|---|
| Riemann Sod | `400 x 5`, `t = 0.2` | Explicit O2 HLLC; Explicit O2 Low-Mach HLLC-P | Check serial/MPI agreement on a report-scale discontinuous row. |
| Gresho vortex | `32 x 32`, Mach `0.01`, `t = 0.4*pi` | Explicit O2 HLLC; Explicit O2 Low-Mach HLLC-P | Check serial/MPI agreement on the low-Mach reference controls. |
| Periodic advection blob | `128 x 128`, `t = 0.5` | Explicit O2 HLLC; Explicit O2 Low-Mach HLLC-P | Check serial/MPI agreement on the advection efficiency endpoint. |
| Shock-density-bubble | `160 x 40`, `t = 0.03` | Explicit O2 HLLC; Explicit O2 Low-Mach HLLC-P | Check serial/MPI agreement for summary and gathered snapshot CSV output. |

The candidate gate compares serial CPU output against a 2-rank MPI run with
matching row counts and numeric CSV differences below the requested tolerance.
Shock-density-bubble rows also compare gathered snapshot CSV files.

## Command

```sh
cd /Users/yuyangjiao/Desktop/MPhilresearch/amrex_euler_allmach_explicit_imex
python3 scripts/run_mpi_explicit_agreement.py \
  --case-set candidate \
  --output /tmp/report2_candidate_mpi_explicit_agreement_2026-05-28 \
  --mpi-ranks 2 \
  --tolerance 1e-10 \
  --row-timeout-sec 600 \
  --force
```

Before running candidate rows, rebuild the local serial and MPI executables:

```sh
cd /Users/yuyangjiao/Desktop/MPhilresearch/amrex_euler_allmach_explicit_imex/amrex/apps/euler_compare
make -j2 DIM=2 USE_MPI=FALSE USE_OMP=FALSE USE_CUDA=FALSE
make -j2 DIM=2 USE_MPI=TRUE USE_OMP=FALSE USE_CUDA=FALSE
```

## Promotion Rules

A candidate row can support later MPI timing only if:

- serial and MPI status values are `ok` or `passed`;
- failure categories are `none` or `ok`;
- final CSV row counts match;
- maximum numeric difference is within tolerance;
- shock snapshot CSV comparisons pass where snapshots are requested;
- manifests record a clean branch, build flags, command line, rank count, and
  output hashes.

This agreement matrix is not itself performance evidence. It is the gate before
repeated timing runs and accuracy-versus-time plots.

## Candidate Rank-Scan Timing

After the candidate agreement matrix passes, use the rank-scan runner to collect
first-pass MPI timing rows without per-cell final CSV output:

```sh
cd /Users/yuyangjiao/Desktop/MPhilresearch/amrex_euler_allmach_explicit_imex
python3 scripts/run_mpi_explicit_rank_scan.py \
  --case-set candidate \
  --output /tmp/report2_candidate_mpi_rank_scan_2026-05-28 \
  --ranks 1,2,4 \
  --repeats 3 \
  --row-timeout-sec 600 \
  --force
```

This runner records both the Python driver wall time and the application
reported wall time. It disables final CSV writing by default so the first
rank-scan is closer to solver/log timing rather than output-throughput timing.
Use the agreement matrix, not this timing runner, to validate gathered CSV
contents.

## Latest Local Result

On 2026-05-28, the candidate matrix passed from clean commit
`188d1dff4c36e0f5779e15a36fda82ad0ee3d868`. All eight rows matched serial and
2-rank MPI CSV row counts with maximum numeric difference `0` at tolerance
`1e-10`; shock rows also matched gathered snapshot CSVs. The summary artifact is
`/tmp/report2_candidate_mpi_explicit_agreement_2026-05-28/summary.csv`.

The first candidate rank scan also passed on 2026-05-28 from clean commit
`1b43ab09025a4b7ccd4d4d33d59f4b767d424389`:

- command: `scripts/run_mpi_explicit_rank_scan.py --case-set candidate --ranks 1,2,4 --repeats 3`;
- result: all 72 launches passed;
- summary artifacts:
  `/tmp/report2_candidate_mpi_rank_scan_2026-05-28/summary.csv` and
  `/tmp/report2_candidate_mpi_rank_scan_2026-05-28/rank_summary.csv`.

This rank scan is local candidate trend evidence only. It should guide the
next frozen timing design, not be treated as final Report 2 performance data.

## Not Yet Supported

- CUDA/GPU claims: no local CUDA build or GPU run has passed.
- IMEX MPI/GPU claims: the current IMEX pressure solve remains a host
  single-rank Eigen GMRES path.
- MPI scaling claims: rank scans and repeated timing runs are still required.
- Cross-machine wall-time claims: machine metadata and frozen evidence rows are
  still required.
