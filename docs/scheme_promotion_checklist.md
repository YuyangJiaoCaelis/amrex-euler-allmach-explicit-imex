# Scheme Promotion Checklist

Use this checklist before any new or modified scheme becomes report-facing.

## Identity

- Scheme has a stable report name.
- Runtime selector is unique and not reused from an exploratory variant.
- Selector is documented in `README.md`, `PROJECT_CODE_MAP.md`, and
  `README_PROJECT_SETTINGS.md`.
- Report text and plotting scripts use the same scheme name.

## Numerical Scope

- Governing equations, flux split, reconstruction, timestep policy, and
  pressure-solve treatment are written down.
- Boundary conditions are documented for every report test family.
- Any adaptation from a paper method is labelled as an adaptation, not as a
  full reproduction, unless fully validated.
- Unsupported claim boundaries are explicit, especially for MPI/CUDA/GPU,
  shock-density-bubble, all-Mach behaviour, and pressure-solver scalability.

## Test Coverage

- Smoke row exists for each retained test family where the scheme is claimed:
  Riemann, Gresho vortex, periodic advection blob, and same-gamma
  shock-density-bubble.
- Reference comparison exists where appropriate:
  exact Riemann solution, analytic Gresho/advection references, or a higher
  resolution shock-density-bubble reference.
- Positivity/nonfinite status is checked for the relevant rows.
- Python plotting/running scripts pass syntax checks.

## Reproducibility

- Full command is documented in the reproduction map.
- Output root is outside the source tree or explicitly ignored.
- Candidate run records git commit, clean/dirty state, build flags, machine,
  command line, exit code, and output hashes.
- A frozen rerun reproduces the candidate evidence before the figure is used in
  the report.

## Code Hygiene

- No old prototype, candidate, scaffold, audit, or diary-style path remains in
  the default build.
- Diagnostic fields are either report-facing and documented, or removed from
  live code.
- Third-party code remains under `external/`; project code remains under the
  AMReX app and `scripts/`.
- `main.cpp` remains a thin entry point.

## Promotion Decision

- Numerical evidence supports the exact claim being made.
- Documentation, source, scripts, and verification records are updated in the
  same change.
- The final report-facing figure or table can be traced to one clean commit and
  one reproduction command.
