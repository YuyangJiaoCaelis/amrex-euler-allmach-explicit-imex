#!/usr/bin/env python3
"""Run a compact solver smoke matrix for the reported Euler schemes."""

from __future__ import annotations

import argparse
import csv
import subprocess
import time
from pathlib import Path


def parse_key_values(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in text.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output directory for logs and summary.csv.",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    exe = root / "amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex"
    if not exe.exists():
        raise SystemExit(f"Executable not found: {exe}")

    output = args.output or root.parent / f"{root.name}_smoke_matrix"
    output.mkdir(parents=True, exist_ok=True)

    base = ["max_step=100000"]
    imex_base = [
        "euler.imex_acoustic_startup=0",
        "euler.imex_acoustic_cfl_cap=0.0",
        "euler.imex_pressure_stabilization=off",
    ]
    shock_snapshots = "euler.shock_density_bubble_snapshot_times=0,0.03"

    rows: list[tuple[str, list[str]]] = [
        (
            "riemann_sod_explicit_hllc",
            base
            + [
                "amr.n_cell=128 4",
                "geometry.prob_lo=0 0",
                "geometry.prob_hi=1 0.03125",
                "geometry.is_periodic=0 0",
                "euler.problem=toro1",
                "euler.method=explicit",
                "euler.spatial_order=2",
                "euler.riemann=hllc",
                "euler.field_boundary=exact_dirichlet",
                "stop_time=0.02",
                "euler.cfl=0.45",
            ],
        ),
        (
            "riemann_sod_explicit_lowmach_hllcp",
            base
            + [
                "amr.n_cell=128 4",
                "geometry.prob_lo=0 0",
                "geometry.prob_hi=1 0.03125",
                "geometry.is_periodic=0 0",
                "euler.problem=toro1",
                "euler.method=explicit",
                "euler.spatial_order=2",
                "euler.riemann=xie_am_hllc_p",
                "euler.field_boundary=exact_dirichlet",
                "stop_time=0.02",
                "euler.cfl=0.45",
            ],
        ),
        (
            "riemann_sod_imex_t1s2_bdltv20",
            base
            + imex_base
            + [
                "amr.n_cell=128 4",
                "geometry.prob_lo=0 0",
                "geometry.prob_hi=1 0.03125",
                "geometry.is_periodic=0 0",
                "euler.problem=toro1",
                "euler.method=imex",
                "euler.spatial_order=2",
                "euler.slope_limiter=minmod",
                "euler.imex_form=bdltv20_t1_s2_source_map_picard",
                "euler.bdltv20_paper_t1_s2=toro_xy_exact_dirichlet_2d",
                "euler.bdltv20_paper_pressure_solver=gmres",
                "euler.bdltv20_paper_t1_s2_dt=0.002",
                "euler.field_boundary=exact_dirichlet",
                "stop_time=0.02",
                "euler.imex_picard_iterations=4",
                "euler.imex_solver_tol=1e-10",
                "euler.imex_solver_max_iter=200",
            ],
        ),
        (
            "gresho_m0p01_explicit_hllc_n32",
            base
            + [
                "amr.n_cell=32 32",
                "geometry.prob_lo=0 0",
                "geometry.prob_hi=1 1",
                "geometry.is_periodic=0 0",
                "euler.problem=gresho_vortex",
                "euler.method=explicit",
                "euler.spatial_order=2",
                "euler.riemann=hllc",
                "euler.field_boundary=exact_dirichlet",
                "euler.mach=0.01",
                "stop_time=0.05",
                "euler.cfl=0.45",
            ],
        ),
        (
            "gresho_m0p01_explicit_lowmach_hllcp_n32",
            base
            + [
                "amr.n_cell=32 32",
                "geometry.prob_lo=0 0",
                "geometry.prob_hi=1 1",
                "geometry.is_periodic=0 0",
                "euler.problem=gresho_vortex",
                "euler.method=explicit",
                "euler.spatial_order=2",
                "euler.riemann=xie_am_hllc_p",
                "euler.field_boundary=exact_dirichlet",
                "euler.mach=0.01",
                "stop_time=0.05",
                "euler.cfl=0.45",
            ],
        ),
        (
            "gresho_m0p01_imex_t1s2_bdltv20_n32",
            base
            + imex_base
            + [
                "amr.n_cell=32 32",
                "geometry.prob_lo=0 0",
                "geometry.prob_hi=1 1",
                "geometry.is_periodic=0 0",
                "euler.problem=gresho_vortex",
                "euler.method=imex",
                "euler.spatial_order=2",
                "euler.slope_limiter=minmod",
                "euler.imex_form=bdltv20_t1_s2_source_map_picard",
                "euler.bdltv20_paper_t1_s2=gresho_exact_dirichlet_2d",
                "euler.bdltv20_paper_pressure_solver=gmres",
                "euler.bdltv20_paper_t1_s2_dt=0.002",
                "euler.field_boundary=exact_dirichlet",
                "euler.mach=0.01",
                "stop_time=0.05",
                "euler.imex_picard_iterations=4",
                "euler.imex_solver_tol=1e-10",
                "euler.imex_solver_max_iter=200",
            ],
        ),
        (
            "advection_blob_explicit_hllc_n32",
            base
            + [
                "amr.n_cell=32 32",
                "geometry.prob_lo=0 0",
                "geometry.prob_hi=1 1",
                "geometry.is_periodic=1 1",
                "euler.problem=advection_blob",
                "euler.method=explicit",
                "euler.spatial_order=2",
                "euler.riemann=hllc",
                "stop_time=0.05",
                "euler.cfl=0.45",
            ],
        ),
        (
            "advection_blob_explicit_lowmach_hllcp_n32",
            base
            + [
                "amr.n_cell=32 32",
                "geometry.prob_lo=0 0",
                "geometry.prob_hi=1 1",
                "geometry.is_periodic=1 1",
                "euler.problem=advection_blob",
                "euler.method=explicit",
                "euler.spatial_order=2",
                "euler.riemann=xie_am_hllc_p",
                "stop_time=0.05",
                "euler.cfl=0.45",
            ],
        ),
        (
            "advection_blob_imex_t1s2_bdltv20_n32",
            base
            + imex_base
            + [
                "amr.n_cell=32 32",
                "geometry.prob_lo=0 0",
                "geometry.prob_hi=1 1",
                "geometry.is_periodic=1 1",
                "euler.problem=advection_blob",
                "euler.method=imex",
                "euler.spatial_order=2",
                "euler.slope_limiter=minmod",
                "euler.imex_form=bdltv20_t1_s2_source_map_picard",
                "euler.bdltv20_paper_t1_s2=advection_blob_periodic_2d",
                "euler.bdltv20_paper_pressure_solver=gmres",
                "euler.bdltv20_paper_t1_s2_dt=0.002",
                "stop_time=0.05",
                "euler.imex_picard_iterations=4",
                "euler.imex_solver_tol=1e-10",
                "euler.imex_solver_max_iter=200",
            ],
        ),
        (
            "shock_density_bubble_explicit_hllc_160x40",
            base
            + [
                shock_snapshots,
                "euler.problem=shock_density_bubble_2d",
                "euler.method=explicit",
                "euler.spatial_order=2",
                "euler.riemann=hllc",
                "amr.n_cell=160 40",
                "stop_time=0.03",
                "euler.cfl=0.45",
                f"euler.shock_density_bubble_snapshot_dir={output / 'snap_hllc'}",
            ],
        ),
        (
            "shock_density_bubble_explicit_lowmach_hllcp_160x40",
            base
            + [
                shock_snapshots,
                "euler.problem=shock_density_bubble_2d",
                "euler.method=explicit",
                "euler.spatial_order=2",
                "euler.riemann=xie_am_hllc_p",
                "amr.n_cell=160 40",
                "stop_time=0.03",
                "euler.cfl=0.45",
                f"euler.shock_density_bubble_snapshot_dir={output / 'snap_lmhllcp'}",
            ],
        ),
        (
            "shock_density_bubble_imex_t1s2_bdltv20_160x40",
            base
            + imex_base
            + [
                shock_snapshots,
                "euler.problem=shock_density_bubble_2d",
                "euler.method=imex",
                "euler.spatial_order=2",
                "euler.slope_limiter=minmod",
                "euler.imex_form=bdltv20_t1_s2_source_map_picard",
                "euler.bdltv20_paper_t1_s2=off",
                "euler.bdltv20_paper_pressure_solver=gmres",
                "euler.bdltv20_paper_t1_s2_dt=0.002",
                "amr.n_cell=160 40",
                "stop_time=0.03",
                "euler.imex_picard_iterations=4",
                "euler.imex_solver_tol=1e-10",
                "euler.imex_solver_max_iter=200",
                f"euler.shock_density_bubble_snapshot_dir={output / 'snap_imex'}",
            ],
        ),
    ]

    summary = output / "summary.csv"
    fields = [
        "row",
        "returncode",
        "status",
        "failure_category",
        "rho_min",
        "pressure_min",
        "completed_steps",
        "completed_time",
        "wall_time_sec",
        "log",
    ]
    with summary.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for name, row_args in rows:
            log = output / f"{name}.log"
            plot_dir = output / "plotfiles" / name
            plot_dir.mkdir(parents=True, exist_ok=True)
            start = time.time()
            proc = subprocess.run(
                [str(exe)] + row_args + [f"amr.plot_file={plot_dir / 'plt'}"],
                cwd=exe.parent,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            elapsed = time.time() - start
            log.write_text(proc.stdout)
            values = parse_key_values(proc.stdout)
            status = (
                values.get("project_amrex_euler_compare_status")
                or values.get("bdltv20_paper_t1_s2_status")
                or ""
            )
            failure = values.get(
                "failure_category", "none" if status in {"ok", "passed"} else ""
            )
            rho_min = values.get("rho_min") or values.get("bdltv20_paper_rho_min", "")
            pressure_min = values.get("pressure_min") or values.get(
                "bdltv20_paper_pressure_min", ""
            )
            completed_steps = values.get("completed_steps") or values.get(
                "bdltv20_paper_steps", ""
            )
            completed_time = values.get("completed_time") or values.get(
                "bdltv20_paper_final_time", ""
            )
            writer.writerow(
                {
                    "row": name,
                    "returncode": proc.returncode,
                    "status": status,
                    "failure_category": failure,
                    "rho_min": rho_min,
                    "pressure_min": pressure_min,
                    "completed_steps": completed_steps,
                    "completed_time": completed_time,
                    "wall_time_sec": f"{elapsed:.6g}",
                    "log": str(log),
                }
            )
            print(name, proc.returncode, status, failure, rho_min, pressure_min)
            if proc.returncode != 0 or status not in {"ok", "passed"} or failure not in {
                "ok",
                "none",
            }:
                return proc.returncode or 2
    print(f"summary {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
