# Full Report-Data Reproduction Pass - 2026-05-26

## Scope

This pass reran the full Report-1 AMReX evidence set from the clean structured
source package:

`/Users/yuyangjiao/Desktop/MPhilresearch/amrex_euler_research_codebase_2026-05-26`

The run used the external-output wrapper:

```bash
amrex/apps/euler_compare/benchmarks/report_reproduction/run_report_reproduction.sh \
  /Users/yuyangjiao/Desktop/MPhilresearch/research_codebase_full_reproduction_2026-05-26_clean_structured_190802
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

All Chapter-4 report figures were compared against the final Overleaf source.
After refreshing the final Overleaf source with the regenerated figures, all 16
Chapter-4 figure files in the source and in the Overleaf zip matched the rerun
outputs by file hash.

## Timing Note

Figures with a wall-clock-time axis are empirical serial CPU timings. Repeated
runs should preserve scheme ordering and trends, while exact wall times can
change with machine load, thermal state, compiler/runtime state, and hardware.
