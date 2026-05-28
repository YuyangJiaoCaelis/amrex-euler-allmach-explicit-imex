# IMEX AP-Behaviour Probe, 2026-05-28

Status: current candidate AP-like behaviour note; not frozen report evidence.

This is a numerical AP-behaviour probe, not a formal AP proof.

## Commands

Primary mixed IMEX/explicit run:

```text
python3 scripts/run_imex_ap_probe.py --output-dir /Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/imex_ap_probe_2026-05-28 --cases imex_source_map,imex_direct_bdltv20,explicit_hllc,explicit_xie_am_hllc_p --machs 0.1,0.01,0.001 --grid 16 --target-time 0.05 --row-timeout-sec 180
```

Larger IMEX-only check:

```text
python3 scripts/run_imex_ap_probe.py --output-dir /Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/imex_ap_probe_2026-05-28_n32_imex --cases imex_source_map,imex_direct_bdltv20 --machs 0.1,0.01,0.001 --grid 32 --target-time 0.05 --row-timeout-sec 240
```

Both runs used:

- branch: `codex/report2-next-research-step`;
- base commit: `96e7f1f50ef9c138b322ed9bf8d5538f0a14c1fb`;
- executable: `amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex`;
- problem: Gresho vortex with exact Dirichlet field boundary;
- low-Mach construction: Report-1-style high thermodynamic background pressure;
- direct-driver epsilon: default `euler.bdltv20_paper_epsilon=1.0` unless
  overridden in a future run;
- IMEX Picard iterations: `4`;
- IMEX solver tolerance: `1e-8`;
- CFL values: explicit `0.4`, IMEX `0.8`;
- output roots outside the source repo under
  `/Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/`.

The AP probe script classifies IMEX rows with step-count ratio above `1.5` as
failing the material-timestep probe, and pressure-perturbation relative-error
ratio above `100` as material-timestep-only evidence. These thresholds are
screening flags for candidate review, not mathematical pass/fail criteria.

The manifests record a dirty tree because this AP note and probe runner were
new uncommitted files during the candidate run. Treat these rows as candidate
evidence, not frozen report evidence.

## Result Summary

Grid `16 x 16`, target time `0.05`:

| Case | Mach 0.1 steps | Mach 0.01 steps | Mach 0.001 steps | Pressure-perturbation relative L1 trend |
|---|---:|---:|---:|---|
| Source-map IMEX | 2 | 2 | 2 | `1.27e-1` to `1.17e3` |
| Direct BDLTV20 IMEX | 2 | 2 | 2 | `1.96e-2` to `1.71e-2` |
| Explicit HLLC | 43 | 402 | 4001 | acoustic-CFL step growth |
| Explicit Low-Mach HLLC-P | 43 | 403 | 4003 | acoustic-CFL step growth |

Grid `32 x 32`, target time `0.05`:

| Case | Mach 0.1 steps | Mach 0.01 steps | Mach 0.001 steps | Pressure-perturbation relative L1 trend |
|---|---:|---:|---:|---|
| Source-map IMEX | 4 | 4 | 4 | `1.92e-2` to `1.60e2` |
| Direct BDLTV20 IMEX | 4 | 4 | 4 | `4.59e-3` to `5.72e-3` |

## Interpretation

The explicit HLLC and Low-Mach HLLC-P rows are not AP schemes in the
time-integration sense. Their step counts increase by about two orders of
magnitude as Mach decreases from `0.1` to `0.001`, consistent with acoustic-CFL
control.

The source-map IMEX route passes the material-timestep part of the AP-behaviour
probe: the step count is independent of Mach, and the printed material rate is
constant while the acoustic rate scales with sound speed. However, its pressure
perturbation relative error grows strongly as Mach decreases. This supports the
existing code label `epsilon_1_dimensional_no_formal_ap_claim`; do not claim a
formal AP theorem for this route.

The direct BDLTV20 route is much closer to an AP candidate in this probe: its
step count is Mach-independent and its pressure-perturbation relative error
stays comparable across Mach `0.1`, `0.01`, and `0.001` on both tested grids.
This is promising AP-like low-Mach behaviour evidence, but still not a proof of
the formal `euler.bdltv20_paper_epsilon -> 0` AP limit.

## Report Claim Boundary

Safe claim:

```text
The direct BDLTV20 IMEX route shows candidate AP-like low-Mach behaviour in
short Gresho Mach sweeps, with material-CFL step counts and stable
pressure-perturbation relative errors. The source-map IMEX route currently
supports material-CFL timestepping but not a formal AP claim. Explicit HLLC and
Low-Mach HLLC-P are acoustic-CFL baselines, not AP schemes.
```

Required before a stronger claim:

- promote the discrete epsilon-to-zero derivation for the direct BDLTV20
  pressure row and momentum correction into report prose, with assumptions;
- rerun from a clean/frozen commit if the results are promoted to report
  evidence;
- extend the probe to longer target time and at least one additional grid if it
  becomes a report-facing AP claim;
- add a direct-driver `euler.bdltv20_paper_epsilon` sweep before using the data
  as numerical support for the formal epsilon-limit argument.
