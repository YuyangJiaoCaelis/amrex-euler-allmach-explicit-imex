# MPI Rank-Scan Candidate Timing 2026-05-28

Branch: `codex/report2-next-research-step`

Commit used for timing manifests: `1b43ab09025a4b7ccd4d4d33d59f4b767d424389`

Purpose: collect a first controlled MPI rank-scan timing dataset for retained
explicit Report 2 candidate rows after serial-vs-MPI agreement had passed.

## Command

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

The runner used the MPI executable for all rows and did not write per-cell
final CSV files, so these timings are closer to solver/log timings than output
throughput timings. Output CSV correctness is covered separately by
`verification/mpi_candidate_agreement_2026-05-28.md`.

## Result

All 72 launches passed: 8 candidate cases x 3 ranks x 3 repeats.

The first manifest checked recorded:

- branch: `codex/report2-next-research-step`;
- dirty: `false`;
- output class: `candidate`;
- commit: `1b43ab09025a4b7ccd4d4d33d59f4b767d424389`.

Median Python-driver wall times and speedups:

| Case | MPI 1 median s | MPI 2 speedup | MPI 4 speedup |
|---|---:|---:|---:|
| Riemann Sod HLLC | 0.6766 | 1.53 | 1.89 |
| Riemann Sod low-Mach HLLC-P | 0.8990 | 1.62 | 2.28 |
| Gresho HLLC | 4.4190 | 1.60 | 1.32 |
| Gresho low-Mach HLLC-P | 17.7957 | 1.99 | 2.56 |
| Advection blob HLLC | 2.3160 | 1.79 | 2.15 |
| Advection blob low-Mach HLLC-P | 4.7859 | 1.84 | 2.22 |
| Shock-density-bubble HLLC | 0.2102 | 1.06 | 1.10 |
| Shock-density-bubble low-Mach HLLC-P | 0.2564 | 1.15 | 1.26 |

Artifacts:

- `/tmp/report2_candidate_mpi_rank_scan_2026-05-28/summary.csv`
- `/tmp/report2_candidate_mpi_rank_scan_2026-05-28/rank_summary.csv`

## Claim Boundary

This is candidate local MPI trend evidence on a macOS ARM64 workstation, not
frozen report evidence and not a cross-machine benchmark. The shock rows are
too small for meaningful scaling interpretation. Before final Report 2 claims,
repeat timing on the target hardware, freeze the evidence branch, and generate
accuracy-versus-time plots from frozen outputs.
