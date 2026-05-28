# Report 2 Research Workflow

This document records the operating rules for Report 2 development in this
clean AMReX repository. Longer reviews, audits, and planning notes belong in
the parent-level `Review/` folder, not in this source tree.

The goal is to prevent exploratory work, one-off diagnostics, and stale output
from becoming confused with report-facing schemes or figures.

## Repository Roles

- `main`: polished public baseline. Update only after a successful frozen
  reproduction.
- `codex/report2-next-research-step`: active Report 2 development branch.
- `exp/<short-name>`: short-lived experiment branches for CUDA, MPI,
  pressure-solver, or numerical-method trials.
- `report2/frozen-YYYY-MM-DD`: frozen evidence branches for report-facing
  figures. These branches should be tagged and not rewritten.

No Report 2 figure should be promoted unless it can be traced to a frozen branch
or tag, a clean working tree, and a recorded command.

## Output Policy

Generated outputs must stay outside the nested source repo whenever practical.

Use external run roots with three classes:

- `runs/exploratory/`: learning runs. Do not cite in the report.
- `runs/candidate/`: possible evidence rows awaiting checks.
- `runs/frozen/`: report-facing reruns from frozen branches or tags.

The source repo should not contain `results/`, `project_outputs/`, plotfile
directories, executable build products, object files, Python bytecode, or crash
backtraces.

## Scheme Boundaries

Report-facing schemes are the documented selectors in `README.md`,
`PROJECT_CODE_MAP.md`, and `README_PROJECT_SETTINGS.md`.

Exploratory schemes or diagnostic variants must not share a report-facing
selector. If an exploratory implementation is needed, keep it on an experiment
branch or behind an explicit build/runtime guard. Do not leave failed prototypes
or alternate candidate paths reachable from the default build.

When a scheme is promoted, update the scheme selector documentation in the same
change as the source and runner updates.

## Run Metadata

Every candidate or frozen run should record enough metadata to reconstruct the
evidence:

- git branch, commit, and dirty/clean status;
- build flags, especially `USE_MPI`, `USE_OMP`, `USE_CUDA`, `DEBUG`,
  `PROFILE`, and compiler;
- AMReX/Eigen provenance or package version;
- machine tag, CPU/GPU details, MPI ranks, OpenMP threads, and CUDA device;
- full command line and input-file hashes;
- start/end timestamps, wall time, exit code, and output root;
- hashes of emitted CSV, PNG, GIF, or table artifacts used in the report.

Use `scripts/run_manifest.py` for new AMReX runner work whenever possible. It
can be imported by Python runners or used as a command wrapper from shell
scripts.

A wall-time figure is not hardware evidence unless the metadata records the
backend and machine used to generate it.

## Figure Promotion Rule

A generated figure may become report evidence only when:

1. its source command is documented in `PROJECT_PLOT_REPRODUCTION.md` or a
   Report 2 reproduction map;
2. its source data are from a candidate or frozen run, not an exploratory run;
3. the run records clean git state and relevant build flags;
4. any required numerical checks have passed;
5. the figure/csv hashes are recorded under `verification/`;
6. captions and report text state the supported claim boundary.

Serial CPU wall time must be labelled as serial CPU timing. MPI, CUDA, GPU, or
multi-core claims require separate backend runs and agreement checks.

## Diagnostics And Audits

Diagnostics that are part of live report evidence may stay in the solver, but
their output columns must be stable and documented.

One-off audits should live as scripts or parent-level review notes. They should
read generated data and report pass/fail summaries; they should not add bulky
temporary counters or trial logic to report-facing solver paths.

## Documentation Cadence

Treat documentation as part of the code:

- source or selector changes require matching updates to `PROJECT_CODE_MAP.md`
  and `README_PROJECT_SETTINGS.md`;
- reproduction-script changes require updates to `PROJECT_PLOT_REPRODUCTION.md`
  or its Report 2 replacement;
- MPI/CUDA readiness changes require updates to
  `REPORT2_HARDWARE_TRANSITION.md`;
- new verification claims require a concise record under `verification/`.

Verification records should describe the current validated state, not become a
long diary of every failed attempt.

## Weekly Pruning

Review exploratory and candidate outputs regularly. Each candidate should be
promoted, demoted, or discarded with a short rationale outside this source repo.

Missing pruning is an early warning that the workspace is becoming overgrown.
