# Report 2 Next-Research Readiness, 2026-05-28

Status: current readiness gate for moving from credibility cleanup into the next
research stage.

## Commit Under Test

- branch: `codex/report2-next-research-step`
- commit: `2585cf6a04ea78696236f9f7dddc41ba04c01779`
- tree state during candidate AP manifests: clean (`git_dirty=false`)

## What Changed

This gate follows the Report 2 credibility cleanup:

- explicit HLLC and Low-Mach HLLC-P formula-to-code trace;
- direct BDLTV20 AP claim-boundary note;
- scheme credibility matrix;
- validation ladder and run-output discipline;
- IMEX pressure-route decision;
- runtime `imex_route_tag` stamp for distinguishing direct BDLTV20 paper-driver
  rows from shock/source-map IMEX rows.

## Clean Smoke Gate

Command:

```text
SMOKE_OUTPUT_DIR=/Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/smoke_route_tag_clean_2585cf6_2026-05-28 amrex/apps/euler_compare/tests/smoke/run_smoke_matrix.sh
```

Result: all 12 smoke rows passed.

The IMEX logs now show:

- Riemann/Gresho/advection direct-driver rows:
  `imex_route_tag=bdltv20_direct_paper_driver`;
- shock-density-bubble source-map row:
  `imex_route_tag=bdltv20_source_map`.

Summary CSV SHA-256:

```text
f75cf5128a1b76c2fba8543217d72fef229af1b10a7d829a4f61af71b7b54ff6  /Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/smoke_route_tag_clean_2585cf6_2026-05-28/summary.csv
```

## Clean AP Candidate Probe

Command:

```text
PYTHONDONTWRITEBYTECODE=1 python3 scripts/run_imex_ap_probe.py --output-dir /Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/imex_ap_probe_clean_2585cf6_2026-05-28 --cases imex_source_map,imex_direct_bdltv20,explicit_hllc,explicit_xie_am_hllc_p --machs 0.1,0.01,0.001 --grid 16 --target-time 0.05 --row-timeout-sec 240 --direct-epsilons 1,0.1,0.01
```

Result: all 18 rows completed with status `ok`; manifests record
`git_dirty=false`.

Key interpretation:

- explicit HLLC and Low-Mach HLLC-P remain acoustic-CFL baselines;
- source-map IMEX keeps material-timestep step counts but pressure perturbation
  error grows strongly at low Mach;
- direct BDLTV20 paper-driver rows keep constant step counts across Mach and
  remain the only credible AP-proof candidate route.

Summary CSV SHA-256:

```text
464dfe1d4cfdb490388f2c5affc8585165940e7792fba20f8d5013646f8d3611  /Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/imex_ap_probe_clean_2585cf6_2026-05-28/summary.csv
```

## Longer Direct-IMEX Epsilon Sweep

Commands:

```text
PYTHONDONTWRITEBYTECODE=1 python3 scripts/run_imex_ap_probe.py --output-dir /Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/imex_direct_epsilon_sweep_clean_2585cf6_n16_2026-05-28 --cases imex_direct_bdltv20 --machs 0.01 --grid 16 --target-time 0.4 --row-timeout-sec 300 --direct-epsilons 1,0.1,0.01

PYTHONDONTWRITEBYTECODE=1 python3 scripts/run_imex_ap_probe.py --output-dir /Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/imex_direct_epsilon_sweep_clean_2585cf6_n32_2026-05-28 --cases imex_direct_bdltv20 --machs 0.01 --grid 32 --target-time 0.4 --row-timeout-sec 300 --direct-epsilons 1,0.1,0.01
```

Result: all six direct BDLTV20 rows completed with status `ok`; manifests record
`git_dirty=false`.

Compact results:

| Grid | Epsilon | Steps | Pressure pert. rel. L1 | Velocity L1 | Pressure residual rel. Linf |
|---:|---:|---:|---:|---:|---:|
| 16 | 1 | 14 | `3.91e-2` | `4.31e-2` | `9.88e-9` |
| 16 | 0.1 | 14 | `1.12e-1` | `4.31e-2` | `9.80e-9` |
| 16 | 0.01 | 14 | `1.21e-1` | `4.30e-2` | `9.99e-9` |
| 32 | 1 | 30 | `1.09e-2` | `1.47e-2` | `1.00e-8` |
| 32 | 0.1 | 30 | `1.12e-1` | `1.47e-2` | `9.99e-9` |
| 32 | 0.01 | 30 | `1.22e-1` | `1.47e-2` | `9.88e-9` |

Summary CSV SHA-256:

```text
8470f59c1bbc8d41938eaf65e464f8c52e887059ede0abf2ec656d9227c0571b  /Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/imex_direct_epsilon_sweep_clean_2585cf6_n16_2026-05-28/summary.csv
d3503850fd6cc08e0a1c8b01dfce9baec3c8a9442cfa35796182582b9a9f06ca  /Users/yuyangjiao/Desktop/MPhilresearch/runs/candidate/imex_direct_epsilon_sweep_clean_2585cf6_n32_2026-05-28/summary.csv
```

## Readiness Conclusion

The branch is ready to move into the next research target with these boundaries:

- explicit schemes are the hardware/MPI/CUDA timing target;
- direct BDLTV20 paper-driver IMEX is the AP/all-Mach numerical-method target;
- shock/source-map IMEX is robustness evidence, not the AP-proof route;
- IMEX hardware timing is out of scope until a separate pressure-solver
  experiment is created and validated.

These are candidate readiness checks, not final frozen report figures. Final
Report 2 evidence still requires frozen run roots and figure/table hashes.
