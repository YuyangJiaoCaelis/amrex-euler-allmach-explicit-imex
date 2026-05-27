# Smoke Tests

This directory contains the short pre-flight check for the AMReX Euler
comparison app. The smoke matrix builds the executable and runs the three
reported schemes on the four retained test families:

- exact Riemann;
- Gresho vortex;
- periodic advection blob;
- same-gamma shock-density-bubble.

Run from the package root with:

```bash
amrex/apps/euler_compare/tests/smoke/run_smoke_matrix.sh ../euler_compare_smoke
```

Fresh logs are written to the supplied output directory. If no directory is
given, the wrapper creates a timestamped sibling directory next to the source
package.
