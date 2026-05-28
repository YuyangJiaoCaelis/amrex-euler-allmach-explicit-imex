# Report 2 Explicit Flux Formula-To-Code Trace

Date: 2026-05-28

Scope: Report-facing explicit schemes:

- Explicit O2 HLLC;
- Explicit O2 Low-Mach Corrected HLLC-P, selected by `euler.riemann=xie_am_hllc_p`.

This trace is intended to make the explicit baselines credible for Report 2.
It does not claim that either explicit method is AP. Both remain acoustic-CFL
schemes.

## Shared Finite-Volume Update

The explicit update has the conservative finite-volume form

```text
U_ij^(n+1) = U_ij^n
  - dt/dx * (F_(i+1/2,j) - F_(i-1/2,j))
  - dt/dy * (G_(i,j+1/2) - G_(i,j-1/2)).
```

Code trace:

| Formula role | Code |
|---|---|
| Reconstruct left/right primitive face states | `reconstruct_state` in `reconstruction_fluxes.hpp` |
| Select interface flux | `reconstructed_interface_flux` in `reconstruction_fluxes.hpp` |
| Update conservative variables by flux differences | `euler_fv_one_step` in `explicit_update.hpp` |
| Fill problem-specific ghost cells before update | `fill_problem_ghosts` called inside `euler_fv_one_step` |

The implementation updates `rho`, `mx`, `my`, and `E` only through face-flux
differences. Therefore, away from physical boundary fluxes, the discrete update
has the expected conservative finite-volume structure.

## MUSCL Reconstruction

For second-order explicit rows, the code reconstructs primitive states using
limited slopes:

```text
q_(i+1/2)^L = q_i + 0.5 * limiter(q_i - q_(i-1), q_(i+1) - q_i)
q_(i+1/2)^R = q_(i+1) - 0.5 * limiter(q_(i+2) - q_(i+1), q_(i+1) - q_i)
```

Code trace:

| Formula role | Code |
|---|---|
| Compute one-sided primitive differences | `primitive_slope` |
| Apply configured limiter | `limited_slope` through `primitive_slope` |
| Return first-order state when `spatial_order <= 1` | `reconstruct_state` |
| Return second-order face state otherwise | `reconstruct_state` |

Report-facing rows use `euler.spatial_order=2` and `euler.slope_limiter=minmod`.

## HLLC Flux

The HLLC implementation uses the standard left, right, and contact-wave
structure.

The left and right signal estimates are:

```text
S_L = min(u_n,L - c_L, u_n,R - c_R)
S_R = max(u_n,L + c_L, u_n,R + c_R)
```

The contact speed is:

```text
S_* =
(p_R - p_L + rho_L u_L (S_L - u_L) - rho_R u_R (S_R - u_R))
/
(rho_L (S_L - u_L) - rho_R (S_R - u_R)).
```

The flux returns:

```text
F_L,                         if S_L >= 0
F_L + S_L (U_*L - U_L),       if S_L < 0 <= S_*
F_R + S_R (U_*R - U_R),       if S_* < 0 <= S_R
F_R,                         if S_R <= 0
```

Code trace:

| Formula role | Code |
|---|---|
| Physical Euler flux | `physical_flux` |
| Sound speed | `sound_speed` |
| Normal velocity | `normal_velocity` |
| `S_L`, `S_R` | `hllc_flux` |
| Contact speed `S_*` | `hllc_flux` |
| HLLC star states | `build_hllc_star_state` |
| Degeneracy fallback | `hllc_flux` returns `rusanov_flux` |
| Degeneracy diagnostic count | `count_hllc_degenerate_star_states` |

Credibility claim:

```text
The HLLC row is a conservative explicit finite-volume Euler baseline using
MUSCL/minmod reconstruction and an HLLC approximate Riemann solver, with
Rusanov fallback only for degenerate star-state cases.
```

The fallback rate is measurable through the existing
`count_hllc_degenerate_star_states` diagnostic.

Do not claim:

```text
HLLC is AP or all-Mach in the timestep sense.
```

The timestep uses the acoustic spectral rate:

```text
(|u| + c) / dx + (|v| + c) / dy.
```

That is the expected explicit compressible Euler CFL, and it intentionally makes
HLLC a cost baseline against IMEX at low Mach.

## Low-Mach Corrected HLLC-P

The report-facing low-Mach corrected explicit row uses:

```text
euler.riemann=xie_am_hllc_p
```

Here `HLLC-P` means a pressure-corrected HLLC-family flux, following the
Xie-style all-Mach/low-Mach pressure correction. This is not a relabelled HLLC
path. It enters `xie_am_hllc_p_flux`, which:

1. computes HLLC-like wave speeds and star states;
2. forms Roe-like averaged quantities;
3. computes a local Mach-based pressure blending factor `theta`;
4. computes a pressure-ratio shock sensor `f`;
5. modifies the star pressure and adds the pressure correction vector `phi_p`;
6. falls back to standard HLLC if the correction is nonfinite or degenerate.

Code trace:

| Formula role | Code |
|---|---|
| Select Xie low-Mach path | `interface_flux`, case `RiemannKind::XieAmHllcP` |
| Local pressure-ratio sensor | `xie_pressure_ratio_from_states` |
| Multi-face shock sensor | `xie_multiface_pressure_sensor` |
| Cubed sensor sharpening | `xie_cubed_sensor` |
| Local Mach pressure blending | `theta` in `xie_am_hllc_p_flux` |
| Modified star pressure | `p_star_star`, `p_star_star_star` in `xie_am_hllc_p_flux` |
| Pressure correction vector | `phi_p` in `xie_am_hllc_p_flux` |
| Fallback to HLLC | guarded returns inside `xie_am_hllc_p_flux` |

The pressure blending in code has the form:

```text
p** = theta * p_* + (1 - theta) * 0.5 * (p_L + p_R)
p*** = f * p** + (1 - f) * p_*
```

The pressure correction vector has the form:

```text
phi_p = beta_p * (1, u_hat, v_hat, 0.5 * |u_hat|^2)
```

with `beta_p` depending on the pressure jump, wave-speed product, Roe sound
speed, and the pressure sensor.

On strong shocks or degenerate correction states, the HLLC-P path falls back to
standard HLLC by design. The low-Mach correction should therefore be interpreted
as acting mainly in smooth low-Mach regions, while the shock rows check that the
correction path does not damage discontinuous-flow robustness.

Credibility claim:

```text
The Low-Mach HLLC-P row is a conservative explicit finite-volume baseline using
the same MUSCL update as HLLC, but with a separate low-Mach pressure-correction
flux path selected by `xie_am_hllc_p`.
```

Do not claim:

```text
Low-Mach HLLC-P is AP in the timestep/cost sense.
```

It remains explicit and is stepped by the same acoustic-CFL `compute_euler_dt`
rule as HLLC.

## Evidence Already Available

| Evidence | Status | Claim supported |
|---|---|---|
| Exact Riemann rows | Present from Report 1 reproduction | Shock/contact baseline behaviour |
| Gresho Mach sweep | Present from Report 1 reproduction | Low-Mach error/cost comparison |
| Advection blob | Present from Report 1 reproduction | Smooth transport behaviour |
| Shock-density-bubble | Present from Report 1 reproduction | Stress/robustness comparison |
| Serial-vs-2-rank MPI agreement | Passed on 2026-05-28 | Explicit rows are MPI-consistent on candidate rows |
| MPI rank scan | Candidate timing collected on 2026-05-28 | Local candidate MPI trend evidence, not frozen benchmark |

## Remaining Explicit-Scheme Credibility Tasks

Before final Report 2 promotion:

1. Add equation references in report prose for standard HLLC and Xie HLLC-P.
2. Include a small table comparing HLLC and Low-Mach HLLC-P step counts at low
   Mach to show both remain acoustic-CFL explicit schemes.
3. Use Gresho error-vs-cost plots to show whether the low-Mach correction buys
   accuracy per wall time despite acoustic-CFL cost.
4. Keep MPI timing labels specific to machine, rank count, build flags, and
   output class.

## Safe Report Wording

```text
The explicit HLLC method is used as a standard conservative shock-capturing
finite-volume baseline. The Low-Mach HLLC-P method uses the same explicit
finite-volume update but selects a distinct Xie-style low-Mach pressure
correction flux. Both explicit schemes remain acoustic-CFL limited and are not
claimed to be AP schemes.
```
