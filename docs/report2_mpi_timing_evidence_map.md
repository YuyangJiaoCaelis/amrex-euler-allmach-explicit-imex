# Report 2 MPI Timing Evidence Map

This document maps the current candidate MPI timing work to the evidence that
can later appear in Report 2. It separates usable scaling rows from rows that
are only launch/overhead diagnostics.

## Candidate Inputs

- Agreement gate:
  `/tmp/report2_candidate_mpi_explicit_agreement_2026-05-28/summary.csv`
- Rank scan:
  `/tmp/report2_candidate_mpi_rank_scan_2026-05-28/rank_summary.csv`
- Candidate plots:
  `/tmp/report2_candidate_mpi_rank_scan_figures_2026-05-28/`

The candidate rank scan used MPI ranks `1,2,4` and 3 repeats per row. It did
not write per-cell final CSVs during timing; numerical CSV agreement is covered
by the separate agreement gate.

## Candidate Reading

| Family | Current reading | Frozen evidence decision |
|---|---|---|
| Riemann Sod | Moderate-to-useful speedup for both explicit schemes. | Keep in frozen rank scan as a discontinuous 1D-in-2D reference row. |
| Gresho vortex | Low-Mach HLLC-P scales usefully; HLLC regresses at 4 ranks on the small `32 x 32` row. | Keep, but treat the HLLC 4-rank behaviour as overhead/domain-split evidence unless a larger grid is added. |
| Periodic advection blob | Useful speedup for both explicit schemes. | Keep as the cleanest MPI scaling row for the explicit solver. |
| Shock-density-bubble | Very weak speedup at `160 x 40`; timings are launch/overhead dominated. | Do not use `160 x 40` for a scaling claim. Either increase the grid or label it as a launch/output-path sanity row. |

## Recommended Frozen Timing Set

For a Report 2 frozen evidence run:

- branch/tag: `report2/frozen-mpi-timing-YYYY-MM-DD`;
- ranks: `1,2,4`;
- repeats: at least `5`, preferably `7` if wall time allows;
- output: external `runs/frozen/...` root, not inside this source tree;
- final CSV output: disabled during timing rows;
- agreement: rerun `scripts/run_mpi_explicit_agreement.py --case-set candidate`
  first from the same clean commit;
- plots: regenerate with `scripts/plot_mpi_rank_scan.py`.

Suggested frozen rows:

- Riemann Sod `400 x 5`: HLLC and low-Mach HLLC-P;
- Gresho Mach `0.01`: HLLC and low-Mach HLLC-P, with possible larger-grid
  follow-up if the HLLC 4-rank regression remains important;
- Advection blob `128 x 128`: HLLC and low-Mach HLLC-P;
- Shock-density-bubble: increase beyond `160 x 40` for scaling evidence, or
  keep `160 x 40` only as a sanity/overhead row.

## Report Claim Boundary

The current candidate data support the statement that the explicit AMReX rows
now have a reproducible MPI timing workflow and show promising local scaling on
selected rows. They do not yet support final MPI performance claims, GPU claims,
or IMEX parallel-performance claims.
