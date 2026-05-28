# Report 2 Validation Ladder

Date: 2026-05-28

Purpose: define the workflow that should be established before the next research
stage, so new results are credible, reproducible, and not contaminated by trial
runs or overgrown diagnostics.

## Principle

Report 2 should separate three things:

1. development probes that answer a local question;
2. candidate verification runs that may become report evidence;
3. frozen report runs produced from a clean commit with a recorded manifest.

Only frozen report runs should be used for figures or headline claims.

## Ladder

| Level | Gate | Required output | Current status |
|---|---|---|---|
| 0 | Build and hygiene | Clean build, no stale generated files inside source tree | Needs repeated before each frozen run |
| 1 | Selector trace | Command selects the intended scheme route | Partly documented in `docs/report2_scheme_credibility_matrix.md` |
| 2 | Formula-to-code trace | Numerical formula maps to code anchors | Explicit trace and direct IMEX AP boundary documented |
| 3 | Smoke reproduction | Small rows run to completion with manifests | Present from candidate probes |
| 4 | Reference reproduction | Report 1/Report 2 reference scripts reproduce expected rows | Recorded in `verification/reproduction_record.md` |
| 5 | Parallel agreement | Serial/MPI agreement for explicit rows | Recorded in `verification/mpi_candidate_agreement_2026-05-28.md` |
| 6 | Rank-scan timing | Rank-count timing with machine/build labels | Candidate data in `verification/mpi_rank_scan_2026-05-28.md` |
| 7 | AP behaviour | Mach sweep tests timestep and pressure perturbation behaviour, with epsilon scaling labelled explicitly | Candidate high-background-pressure probe in `verification/imex_ap_probe_2026-05-28.md` |
| 8 | Robustness stress | Shock-density-bubble positivity, conservation, residual table | Still needed if shock row is central |
| 9 | Frozen report evidence | Clean commit, clean output root, manifest, hashes, captions | Still needed before final Report 2 figures |

Level 7 is not passed by a tiny smoke run. For an AP-style claim, the minimum
candidate scope is:

- at least one meaningful Gresho characteristic time, not just two to four
  steps;
- at least two grid resolutions, preferably `{16, 32, 64}` for a readable
  trend;
- fixed Picard count and solver tolerance across the sweep;
- explicit labelling of high-background-pressure Mach sweeps versus literal
  `euler.bdltv20_paper_epsilon` sweeps;
- positivity and pressure-residual columns in the summary table.

## Scheme-Specific Gates

| Scheme or route | Must pass before use in Report 2 |
|---|---|
| Explicit O2 HLLC | HLLC formula trace, exact Riemann/Gresho/shock rows, serial/MPI agreement, acoustic-CFL cost label |
| Explicit O2 Low-Mach HLLC-P | Xie flux trace, Gresho low-Mach comparison against HLLC, shock robustness check, acoustic-CFL cost label |
| Direct BDLTV20 T1/S2 IMEX | Direct selector trace, material timestep record, AP derivation boundary, clean Mach sweep, pressure residual record |
| Shock-adapted source-map IMEX | Source-map selector trace, shock positivity/conservation/residual table, clear label that it is not the formal AP proof route |

## Workflow To Avoid Contamination

Use this directory discipline:

| Run type | Location | Rule |
|---|---|---|
| Scratch probes | `/Users/yuyangjiao/Desktop/MPhilresearch/runs/scratch/` | May be deleted; never cite directly |
| Candidate checks | `/Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/` | Keep manifest and verification note; can inform decisions |
| Frozen report runs | `/Users/yuyangjiao/Desktop/MPhilresearch/runs/frozen/report2/` | Only run from a clean commit or tagged commit; use in report figures |
| Review notes | `/Users/yuyangjiao/Desktop/MPhilresearch/Review/` | Human/LLM reviews, not raw numerical evidence |
| Progress notes | `/Users/yuyangjiao/Desktop/MPhilresearch/YuyangProgressNote/` | Project memory and claim-boundary discoveries |

Run outputs should stay outside the code repo unless the file is a deliberate
small verification note. The repo should contain scripts, source code, docs, and
selected verification records, not large trial artifacts.

## Promotion Checklist

Before a candidate run becomes report evidence, record:

- exact command and working directory;
- git branch, commit hash, and `git status --short`;
- build type, compiler, MPI/CUDA setting, and AMReX options;
- output directory and manifest path;
- scheme selector values;
- expected claim category: accuracy, low-Mach correction, AP behaviour,
  robustness, or timing;
- failure/guard counters, pressure-solver residuals, and positivity status where
  relevant.

If any of these are missing, the result remains a candidate check, not a report
result.

## Grandfathered Evidence

Report-1 reproduction rows may remain useful context, but they should be labelled
as grandfathered evidence unless rerun under this ladder with manifests and a
clean output root. Report 2 figures should not silently mix frozen ladder runs
with older unlabelled trial outputs.

## Solver-Code Cool-Off

For report-facing solver files under
`amrex/apps/euler_compare/src/euler_compare/numerics/` or
`amrex/apps/euler_compare/src/euler_compare/imex/`, use a 24-hour second-look
rule before promotion or merge whenever feasible. Documentation and runner
scripts can move faster, but solver changes should get one calm reread before
they become evidence-producing code.

## Immediate Solid Next Step

The next research step should be:

1. keep the emitted `imex_route_tag` visible in IMEX commands, manifests, CSVs,
   and captions;
2. freeze the direct BDLTV20 AP/all-Mach probe from a clean commit into
   `/Users/yuyangjiao/Desktop/MPhilresearch/runs/frozen/report2/imex_ap_probe/`;
3. keep explicit HLLC and Low-Mach HLLC-P in the same Mach sweep as acoustic-CFL
   cost baselines;
4. label whether the run is a high-background-pressure Mach sweep or a literal
   `euler.bdltv20_paper_epsilon` sweep;
5. generate one compact table with step count, pressure-perturbation relative
   error, velocity error, pressure residual, and status;
6. keep the pressure-route decision in force: IMEX remains serial CPU
   numerical-method evidence until a separate pressure-solver experiment is
   validated.

This gives Report 2 a defensible path: explicit schemes as credible baselines,
direct BDLTV20 as the all-Mach route and AP-proof candidate, and shock
source-map IMEX as a clearly labelled adaptation rather than an overclaimed
proof.
