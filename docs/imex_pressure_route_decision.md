# IMEX Pressure Route Decision

Date: 2026-05-28

Decision: for the next Report 2 research stage, keep the current BDLTV20 IMEX
pressure solve as a serial CPU host-GMRES route. Do not claim MPI/GPU hardware
efficiency for IMEX until a separate reviewed pressure-solver route exists.

## Why

The project description points toward accuracy-versus-time comparison across
hardware or parallelisation approaches. The explicit HLLC and Low-Mach HLLC-P
paths already use AMReX data structures and device-oriented loops, so they are
the appropriate first target for MPI/CUDA timing work.

The current IMEX routes use host sparse matrix assembly and an Eigen GMRES
pressure solve. That is useful for numerical-method evidence, AP/all-Mach
candidate checks, and shock robustness, but it is not yet a credible MPI/GPU
pressure-solve implementation.

## Report 2 Boundary

Use these labels:

| Route | Report 2 role | Hardware claim allowed now |
|---|---|---|
| Explicit HLLC | Explicit shock-capturing baseline | Serial/MPI/CUDA after agreement checks |
| Explicit Low-Mach HLLC-P | Explicit low-Mach-corrected baseline | Serial/MPI/CUDA after agreement checks |
| Direct BDLTV20 paper-driver IMEX | AP/all-Mach candidate and numerical-method evidence | Serial CPU only |
| Shock/source-map IMEX | Robustness/stress evidence | Serial CPU only |

Do not place IMEX wall times on the same hardware-scaling claim line as explicit
MPI/CUDA timings unless the caption states that IMEX is using the serial host
pressure solve.

## Future Experimental Route

If IMEX hardware timing becomes a research target, use a separate experiment
branch, for example:

```text
exp/mlmg-pressure-solve
```

The experiment must pass these gates before promotion:

1. same pressure equation and boundary policy documented against the current
   host route;
2. serial agreement against the host Eigen GMRES route on Gresho/Riemann rows;
3. MPI agreement and rank-count checks;
4. CUDA/GPU build and agreement checks if GPU timing is claimed;
5. updated AP claim boundary if the pressure operator changes.

Until then, Report 2 should proceed with explicit-kernel hardware work and
IMEX numerical-method validation as separate evidence tracks.
