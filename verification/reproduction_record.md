# Full Reproduction Record

## Scope

This run regenerated the full AMReX evidence set from the structured source
package:

the repository root.

The run used the external-output wrapper:

```bash
amrex/apps/euler_compare/benchmarks/report_reproduction/run_report_reproduction.sh \
  ../euler_compare_full_reproduction
```

## Result

The wrapper completed successfully with exit code 0.

The rerun regenerated the four report test families:

- exact Riemann tests;
- Gresho vortex;
- periodic density-blob advection;
- same-gamma Cartesian shock-density-bubble.

The retained report schemes were:

- Explicit O2 HLLC;
- Explicit O2 Low-Mach Corrected HLLC-P;
- IMEX T1/S2 BDLTV20.

## Output Check

The generated output root contained:

- 100 PNG figures;
- 2 GIF animations;
- 358 CSV files.

A status scan over the relevant rerun summary and metrics CSV files checked 46
CSV files and found 0 non-ok entries.

The generated figures were compared against the final project figure set. After
refreshing the project source with the regenerated figures, all 16 result-figure
files matched the rerun outputs by file hash.

## Timing Note

Figures with a wall-clock-time axis are empirical serial CPU timings. Repeated
runs should preserve scheme ordering and trends, while exact wall times can
change with machine load, thermal state, compiler/runtime state, and hardware.
