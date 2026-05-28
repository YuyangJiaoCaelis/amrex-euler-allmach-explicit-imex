# Report 2 Scheme Credibility Matrix

Date: 2026-05-28

Purpose: define what must be shown, and what must not be claimed, for each
report-facing numerical scheme before moving into the next research stage.

This document is a credibility guide, not a figure-generation record. It should
be updated when a claim is promoted, weakened, or replaced by new evidence.

## Scheme Claim Summary

| Scheme or route | Current role | Credible claim today | Do not claim |
|---|---|---|---|
| Explicit O2 HLLC | Standard explicit finite-volume baseline | Conservative shock-capturing Euler baseline with MUSCL reconstruction and HLLC flux | AP, all-Mach timestep, low-Mach pressure robustness |
| Explicit O2 Low-Mach HLLC-P | Explicit low-Mach-corrected baseline | Same explicit FV baseline plus a low-Mach pressure/flux correction that improves low-Mach behaviour in tests | AP in the timestep/cost sense |
| Direct BDLTV20 T1/S2 IMEX | Main all-Mach IMEX evidence route | BDLTV20 T1/S2, first order in time and second order in space, with material-CFL pressure-split IMEX behaviour and a documented formal epsilon-limit derivation boundary | Complete rigorous AP theorem beyond the stated assumptions |
| Shock-adapted BDLTV20 source-map IMEX | Robustness/stress-test IMEX route | BDLTV20/Toro-Vazquez source-map IMEX adaptation for same-gamma shock-density-bubble | Direct paper-driver equivalence or separately proven AP behaviour |

## Code Traceability Matrix

| Component | Code anchor | Credibility check |
|---|---|---|
| HLLC flux | `amrex/apps/euler_compare/src/euler_compare/numerics/reconstruction_fluxes.hpp`, `hllc_flux` | Check wave speeds, star speed, star states, physical-flux upwinding, and fallback behaviour. |
| Low-Mach HLLC-P | `reconstruction_fluxes.hpp`, `xie_am_hllc_p_flux` | Check Xie-style sensor, pressure correction, fallback to HLLC, and low-Mach limiting behaviour. |
| Explicit FV update | `amrex/apps/euler_compare/src/euler_compare/numerics/explicit_update.hpp`, `euler_fv_one_step` | Check flux-difference conservation, ghost fills, reconstruction, and diagnostics. |
| Explicit timestep | `amrex/apps/euler_compare/src/euler_compare/numerics/diagnostic_metrics_timestep.hpp`, `compute_euler_dt` | Records acoustic CFL control through `|u| + c`; this is why explicit schemes are not AP. |
| Direct BDLTV20 driver | `amrex/apps/euler_compare/src/euler_compare/imex/imex_bdltv20_direct_driver.hpp` | Check BDLTV20 T1/S2 formula path, material timestep, pressure row, momentum correction, and energy closeout. |
| Source-map IMEX driver | `amrex/apps/euler_compare/src/euler_compare/imex/imex_source_map_step.hpp` | Check shock-adapted pressure row, host GMRES, source-map diagnostics, positivity, and stated limitations. |
| IMEX route dispatch | `amrex/apps/euler_compare/src/euler_compare/app/driver_main.hpp` | Check that non-`off` `euler.bdltv20_paper_t1_s2` dispatches to the direct driver. |
| IMEX route disambiguation | Runtime stdout/CSV stamp | Report-facing IMEX rows now emit `imex_route_tag=bdltv20_direct_paper_driver` or `imex_route_tag=bdltv20_source_map`; captions should still show `euler.bdltv20_paper_t1_s2`. |
| Runtime consistency checks | `amrex/apps/euler_compare/src/euler_compare/core/config_reader.hpp` | Check that report-facing selector combinations are accepted or rejected intentionally. |

## Evidence Family Matrix

| Evidence family | HLLC | Low-Mach HLLC-P | Direct BDLTV20 T1/S2 IMEX | Shock-adapted source-map IMEX |
|---|---|---|---|---|
| Riemann | Shock/contact baseline | Shock/contact baseline with low-Mach correction inactive or limited where appropriate | Direct BDLTV20 paper-driver evidence | Not the main route |
| Gresho Mach sweep | Explicit acoustic-CFL cost baseline | Main low-Mach explicit correction evidence | Main all-Mach/AP-behaviour evidence | Useful only as secondary diagnostic |
| Periodic advection blob | Smooth transport baseline | Smooth transport low-Mach corrected baseline | Smooth IMEX accuracy/efficiency evidence | Not the main route |
| Shock-density-bubble | Shock robustness baseline | Shock robustness baseline with low-Mach correction | Not currently supported by direct paper driver | Main IMEX shock robustness route |
| MPI rank scan | Current explicit parallel evidence | Current explicit parallel evidence | Not supported by host direct-driver pressure solve | Not supported by host source-map pressure solve |
| CUDA/GPU | Future explicit backend evidence | Future explicit backend evidence | Future work only after pressure solver route decision | Future work only after pressure solver route decision |

## Required Credibility Work By Scheme

### Explicit O2 HLLC

What makes it credible:

- Conservative finite-volume update from shared face fluxes.
- Standard HLLC wave-speed and star-state construction.
- Riemann benchmark agreement against exact solutions.
- Shock-density-bubble comparison against the high-resolution HLLC reference.
- Serial/MPI agreement and rank-scan timing for explicit rows.

Remaining checks before stronger Report 2 claims:

- Add a short equation-to-code note mapping HLLC formulas to `hllc_flux`.
- Keep explicit timing claims separated by backend: serial CPU, MPI CPU, CUDA GPU.
- Do not describe HLLC as all-Mach or AP; its timestep is acoustic-CFL.

### Explicit O2 Low-Mach HLLC-P

What makes it credible:

- It is not just relabelled HLLC: the Xie low-Mach pressure correction path is
  implemented separately from standard HLLC.
- Gresho low-Mach rows can show whether pressure/velocity errors improve over
  standard HLLC.
- Riemann and shock-density-bubble rows check that the correction does not break
  discontinuous-flow robustness.

Remaining checks before stronger Report 2 claims:

- Write a formula-to-code trace for the Xie sensor, modified pressure star, and
  pressure correction term.
- Present it as a low-Mach corrected explicit flux, not an AP scheme.
- Keep the explicit acoustic step-count growth as part of the evidence, because
  it clarifies the difference between low-Mach correction and AP timestepping.

### Direct BDLTV20 T1/S2 IMEX

What makes it credible:

- It is the cleaner BDLTV20 T1/S2 paper-driver route used by Riemann, Gresho, and
  advection evidence.
- It uses material-CFL timestepping in the configured report rows.
- The focused AP-behaviour probe from 2026-05-28 showed stable step counts and
  pressure-perturbation relative errors across Mach `0.1`, `0.01`, and `0.001`
  on short Report-1-style high-background-pressure Gresho probes.
- The formal epsilon-limit derivation boundary is documented separately in
  `docs/report2_imex_ap_claim_boundary.md`.

Remaining checks before a formal AP/all-Mach claim:

- Write the discrete epsilon-to-zero derivation for the direct-driver pressure
  row and pressure-gradient momentum correction into report-ready prose.
- State assumptions: well-prepared data, fixed mesh, material-CFL timestep,
  sufficiently converged pressure solve, and no fallback guards changing the
  limiting equation.
- Do not conflate the existing high-background-pressure Mach sweep with a
  literal `euler.bdltv20_paper_epsilon -> 0` proof.
- Rerun AP-behaviour evidence from a clean or frozen commit if used in Report 2.
- Extend the AP probe beyond the short smoke horizon before using it as a
  report-facing result.

### Shock-Adapted BDLTV20 Source-Map IMEX

What makes it credible:

- It preserves the BDLTV20/Toro-Vazquez pressure-split structure: star state,
  pressure solve, pressure-gradient momentum correction, and energy closeout.
- It adapts the route to the same-gamma shock-density-bubble case and records
  positivity, residual, pressure-solver, and boundary diagnostics.
- It provides stress-test evidence that the IMEX idea can be exercised beyond
  smooth benchmark rows.

Remaining checks before stronger claims:

- Keep it labelled as a shock-adapted source-map IMEX route.
- Do not use it as the formal AP proof route.
- Add a small positivity/conservation/residual table for shock-density-bubble
  if it becomes a central Report 2 result.
- Decide whether future Report 2 work will keep it as a CPU-only robustness row
  or replace the pressure solve with an AMReX-native parallel path.

## Minimum Next-Stage Credibility Gate

Before promoting new Report 2 figures or claims, each retained row should pass:

1. Selector check: the command uses the intended report-facing scheme selector.
2. Code-trace check: the selector maps to the documented flux or IMEX route.
3. Numerical-status check: status is `ok` or explicitly explained.
4. Positivity check: no nonfinite state, nonpositive density, or nonpositive
   pressure unless the row is explicitly a failed stress test.
5. Claim-boundary check: the caption states whether the row supports accuracy,
   low-Mach behaviour, AP-like behaviour, robustness, or hardware timing.
6. Reproducibility check: command, git state, build flags, output root, and
   hashes are recorded by the run manifest.

The broader workflow gate is recorded in `docs/report2_validation_ladder.md`.

## Immediate Follow-Up Tasks

| Priority | Task | Output |
|---|---|---|
| P0 | Maintain the two-route IMEX stamp in all report-facing outputs | emitted `imex_route_tag`; optional future direct `imex_form` |
| P1 | Maintain HLLC and Low-Mach HLLC-P formula-to-code trace | `docs/report2_explicit_flux_trace.md` |
| P1 | Maintain the direct BDLTV20 AP derivation boundary | expanded `docs/report2_imex_ap_claim_boundary.md` |
| P1 | Use the validation ladder before promoting new figures | `docs/report2_validation_ladder.md` |
| P1 | Rerun AP probe from a clean commit if promoted | frozen/candidate run directory plus verification record |
| P2 | Add shock-source-map robustness table | verification note or report table source |
| P2 | Maintain IMEX pressure-solver route boundary for Report 2 hardware claims | `docs/imex_pressure_route_decision.md`; `REPORT2_HARDWARE_TRANSITION.md` |

## Safe Report Wording

Use:

```text
The explicit HLLC and Low-Mach HLLC-P methods are acoustic-CFL explicit
baselines. Low-Mach HLLC-P is assessed as a low-Mach flux correction, not as an
AP time integrator. The direct BDLTV20 T1/S2 route is the main AP-proof
candidate and all-Mach IMEX evidence route. The shock-density-bubble row uses a
related BDLTV20 source-map adaptation and is reported as robustness/stress
evidence.
```

Avoid:

```text
All schemes are AP.
```

Avoid:

```text
All IMEX rows use the same direct BDLTV20 paper-driver implementation.
```

Avoid:

```text
The shock-density-bubble IMEX adaptation is formally AP-proven.
```
