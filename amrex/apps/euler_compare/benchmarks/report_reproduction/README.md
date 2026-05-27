# Report-Data Reproduction Benchmark

This directory contains the long reproduction entry point for the AMReX Euler
comparison results. It rebuilds the executable and then runs the same command
groups listed in `PROJECT_PLOT_REPRODUCTION.md`.

Run all report-data groups from the package root with an external output
directory:

```bash
amrex/apps/euler_compare/benchmarks/report_reproduction/run_report_reproduction.sh ../euler_compare_full_reproduction
```

The wrapper creates temporary `results` and `project_outputs` links during the
run and removes the links afterward, so generated data are kept outside the
source package. Running the script without an output directory uses the legacy
in-package output locations from `reproduce_project_data.sh`.

The full run can take a long time on a laptop, especially the high-resolution
Gresho rows and the 640x160 shock-density-bubble reference row.
