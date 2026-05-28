# Report 2 IMEX AP Claim Boundary

This note separates three different statements that are easy to blur:

1. an asymptotic-preserving mathematical proof;
2. numerical AP-like evidence from Mach sweeps;
3. a low-Mach correction for an otherwise explicit acoustic-CFL method.

## AP Meaning Used Here

For the low-Mach Euler limit, an AP scheme should satisfy both properties below
for well-prepared data as the nondimensional Mach parameter tends to zero:

- the stable timestep is controlled by material velocity, not sound speed;
- at fixed mesh and timestep, the discrete update has a finite incompressible
  limit, usually a pressure/velocity projection enforcing a discrete divergence
  constraint.

Numerical Mach sweeps can support this claim, but they do not prove it by
themselves. A proof must take the discrete equations and pass to the
epsilon-to-zero limit.

## Current IMEX Paths

There are two relevant IMEX routes in this codebase.

The normal source-map driver is selected by:

```text
euler.method=imex
euler.imex_form=bdltv20_t1_s2_source_map_picard
euler.bdltv20_paper_t1_s2=off
```

This route is used by the general AMReX time loop and by the shock-density-bubble
IMEX rows. Its diagnostics deliberately state:

```text
source_scaling = epsilon_1_dimensional_no_formal_ap_claim
```

Because this path has no active nondimensional epsilon parameter in the update,
it cannot currently receive a formal epsilon-to-zero AP proof as implemented.
It can only be tested for AP-like low-Mach behaviour through dimensional Mach
sweeps.

The direct BDLTV20 T1/S2 driver is selected by:

```text
euler.method=imex
euler.imex_form=bdltv20_t1_s2_source_map_picard
euler.bdltv20_paper_t1_s2=<case>
```

This route is used by several smooth/Riemann report-generation scripts. It has
an explicit `euler.bdltv20_paper_epsilon` parameter, stores

```text
E = p/(gamma - 1) + epsilon * rho |u|^2 / 2
```

and applies the implicit pressure correction with a `1/epsilon` pressure
gradient. This is the route for which the formal AP consistency argument below
has direct code support, but only for the epsilon-scaled direct-driver scheme
and under well-prepared low-Mach assumptions.

## Selector Route Stamp

Both IMEX routes share the legacy public selector

```text
euler.imex_form=bdltv20_t1_s2_source_map_picard
```

and are separated by the secondary value of `euler.bdltv20_paper_t1_s2`. To make
this visible in evidence, the code now emits an `imex_route_tag` in report-facing
stdout and CSV outputs.

Safe route labels are:

| Route label | Runtime stamp | Required selector evidence |
|---|---|
| Direct BDLTV20 paper driver | `imex_route_tag=bdltv20_direct_paper_driver` | `euler.bdltv20_paper_t1_s2=<case other than off>` |
| Shock/source-map IMEX adaptation | `imex_route_tag=bdltv20_source_map` | `euler.bdltv20_paper_t1_s2=off` |

Future cleanup may still add a distinct direct-driver `imex_form` value, but the
route stamp is the current guard against ambiguous report labels.

## Formal AP Derivation Boundary For The Direct Driver

The statement below is the strongest claim that currently has direct code
support. It is still a formal consistency argument, not a completed journal-level
proof.

For the direct BDLTV20 T1/S2 driver, assume:

- well-prepared low-Mach data with pressure split
  `p^n = p0^n + epsilon * pi^n + O(epsilon^2)`;
- `grad_h p0^n = 0`, or boundary data compatible with a spatially constant
  thermodynamic pressure background;
- density and velocity remain `O(1)` as `epsilon -> 0`;
- the mesh and material-CFL timestep are fixed independently of sound speed;
- the pressure solve is converged accurately enough for the residual to be
  smaller than the leading-order low-Mach balance;
- no fallback or positivity guard changes the limiting equation.

The direct driver stores scaled energy as

```text
E = p/(gamma - 1) + epsilon * rho |u|^2 / 2.
```

Its explicit material step forms a star state `U*` with material Lax-Friedrichs
transport. In compact finite-volume notation,

```text
U* = U^n - dt * D_h F_material(U^n),
```

where the code uses material wave speeds, not acoustic wave speeds, in this
star update. The pressure correction then has the form

```text
(rho u)^(n+1) = (rho u)^* - dt / epsilon * grad_h p^(n+1),
```

implemented in code by centred pressure differences divided by `epsilon`.

The pressure row has the structure

```text
epsilon/(gamma - 1) * p^(n+1)
  - dt^2 * L_h(h_hat, p^(n+1))
  = epsilon * R^* + boundary terms,
```

where `L_h` is the BDLTV20 face-enthalpy pressure operator assembled from the
lagged Picard state. This is the same structure recorded by the direct driver:
`eps / (gamma - 1)` on the equation-of-state diagonal, pressure-operator
coefficients proportional to `dt^2`, and a right-hand side multiplied by
`epsilon`.

Insert the well-prepared expansion

```text
p^(n+1) = p0^(n+1) + epsilon * pi^(n+1) + O(epsilon^2).
```

Since `grad_h p0 = 0` and the pressure operator annihilates a compatible
constant background pressure, the pressure correction becomes

```text
(rho u)^(n+1)
  = (rho u)^* - dt * grad_h pi^(n+1) + O(epsilon).
```

Thus the `1/epsilon` factor does not make the limiting momentum update singular;
it converts the small pressure perturbation into the finite incompressible
pressure correction.

For the pressure equation, subtract the constant thermodynamic-pressure balance
and divide by `epsilon`. The leading-order equation has the form

```text
- dt^2 * L_h(h0, pi^(n+1))
  = R0^* - p0^(n+1)/(gamma - 1) + lower-order known terms,
```

with the precise lower-order terms determined by the material star update,
enthalpy lagging, and boundary treatment. This is a finite elliptic equation for
the dynamic pressure `pi^(n+1)`. It supplies the pressure projection needed to
obtain the incompressible-limit velocity constraint at the discrete level.

The formal AP mechanism is therefore:

1. the timestep remains material-CFL controlled;
2. the pressure perturbation scaling cancels the apparent `1/epsilon`
   singularity in the momentum correction;
3. the pressure row tends to a finite elliptic dynamic-pressure solve;
4. the energy closeout returns
   `E^(n+1) = p^(n+1)/(gamma - 1) + epsilon * kinetic`, so kinetic energy does
   not reintroduce acoustic stiffness in the limit.

This derivation should only be attached to the direct epsilon-scaled BDLTV20
driver. It does not automatically apply to the source-map shock route, because
that route currently records dimensional `epsilon_1` source scaling and has no
active epsilon parameter to send to zero.

## Numerical AP Evidence Status

The focused probe recorded in
`verification/imex_ap_probe_2026-05-28.md` supports the distinction above, but
it must be labelled carefully. That probe used the Report-1-style low-Mach setup
with a large thermodynamic background pressure and the direct-driver epsilon at
its recorded run value. It is numerical all-Mach/AP-like evidence, not the same
thing as taking the formal `euler.bdltv20_paper_epsilon -> 0` limit.

- Direct BDLTV20 T1/S2 kept constant step counts as Mach decreased from `0.1`
  to `0.001`, and pressure-perturbation relative errors stayed bounded on the
  short Gresho probe.
- The shock-adapted source-map IMEX route also kept material-timestep step
  counts, but its pressure-perturbation relative error degraded strongly at low
  Mach in the same probe.
- Explicit HLLC and Low-Mach HLLC-P showed the expected acoustic-CFL step-count
  growth.

That makes the direct BDLTV20 route the credible Report 2 all-Mach candidate and
the only plausible route for an AP proof. The source-map route remains useful
robustness evidence, not the proof route.

## Explicit HLLC And Low-Mach HLLC-P

The explicit HLLC and Low-Mach HLLC-P rows are not AP schemes in the
time-integration/cost sense. Their timestep is computed from the acoustic
spectral radius:

```text
(|u| + c) / dx + (|v| + c) / dy.
```

Therefore their step count grows as sound speed grows. Low-Mach HLLC-P may
reduce low-Mach flux dissipation and improve accuracy, but it remains an
explicit acoustic-CFL baseline.

## Safe Report Claim Today

The safe wording after this review is:

```text
The direct BDLTV20 T1/S2 IMEX path has a documented formal epsilon-limit
consistency argument under well-prepared assumptions, and the current candidate
Mach sweeps show material-CFL all-Mach behaviour. The candidate Mach sweeps are
not, by themselves, a rigorous AP theorem. Explicit HLLC and Low-Mach HLLC-P are
acoustic-CFL explicit baselines, not AP schemes.
```

If the epsilon-scaled derivation is promoted into the report, use wording like:

```text
For well-prepared epsilon-scaled data, the direct BDLTV20 T1/S2 discrete pressure
row formally tends to a finite incompressible pressure equation as epsilon tends
to zero, while the timestep remains material-CFL controlled.
```
