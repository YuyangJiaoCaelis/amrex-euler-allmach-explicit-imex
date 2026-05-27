#!/usr/bin/env python3
"""Run the AMReX Gresho reference matrix.

The reference controls are uniform across all rows: exact Dirichlet field
boundaries in x and y, non-periodic geometry, and final time 0.4*pi.

- Explicit O2 HLLC
- Explicit O2 LowMachCorrected HLLC-P
- AMReX IMEX T1/S2
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex"
DEFAULT_OUT = ROOT / "results/amrex/project_gresho_reference_matrix"


@dataclass(frozen=True)
class Scheme:
    key: str
    label: str
    role: str
    color: str
    marker: str


SCHEMES = (
    Scheme("explicit_o2_hllc", "AMReX Explicit O2 HLLC", "primary_explicit", "#4C78A8", "o"),
    Scheme(
        "explicit_o2_xie_am_hllc_p",
        "AMReX Explicit O2 LowMachCorrected HLLC-P",
        "literature_explicit_low_mach",
        "#F28E2B",
        "P",
    ),
    Scheme(
        "imex_t1s2_bdltv20",
        "IMEX T1/S2 BDLTV20",
        "primary_imex_target",
        "#2ca02c",
        "^",
    ),
)


def split_csv_floats(text: str) -> list[float]:
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def split_csv_ints(text: str) -> list[int]:
    return [int(item.strip()) for item in text.split(",") if item.strip()]


def split_csv_strings(text: str) -> list[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def pressure_from_mach(mach: float, gamma: float, density: float) -> float:
    return density / (gamma * mach * mach)


def parse_kv(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        if "=" not in line or line.startswith("$ "):
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if re.fullmatch(r"[A-Za-z0-9_]+", key):
            out[key] = value.strip()
    return out


def finite_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        out = float(value)
    except ValueError:
        return None
    return out if math.isfinite(out) else None


def profile_metrics(path: Path, mach: float, gamma: float, density: float) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"empty final CSV: {path}")

    def get(row: dict[str, str], *keys: str) -> float:
        for key in keys:
            value = row.get(key)
            if value not in (None, ""):
                return float(value)
        raise KeyError(keys)

    parsed = []
    for row in rows:
        rho = get(row, "rho")
        u = get(row, "u")
        v = get(row, "v")
        p = get(row, "p", "pressure")
        exact_rho = get(row, "exact_rho")
        exact_u = get(row, "exact_u")
        exact_v = get(row, "exact_v")
        exact_p = get(row, "exact_p", "exact_pressure")
        parsed.append(
            {
                "rho": rho,
                "u": u,
                "v": v,
                "p": p,
                "exact_rho": exact_rho,
                "exact_u": exact_u,
                "exact_v": exact_v,
                "exact_p": exact_p,
                "rho_error": (
                    get(row, "rho_error", "density_error")
                    if row.get("rho_error") not in (None, "")
                    or row.get("density_error") not in (None, "")
                    else rho - exact_rho
                ),
                "u_error": get(row, "u_error") if row.get("u_error") not in (None, "") else u - exact_u,
                "v_error": get(row, "v_error") if row.get("v_error") not in (None, "") else v - exact_v,
                "velocity_error": (
                    get(row, "velocity_error")
                    if row.get("velocity_error") not in (None, "")
                    else math.hypot(u - exact_u, v - exact_v)
                ),
                "p_error": (
                    get(row, "p_error", "pressure_error")
                    if row.get("p_error") not in (None, "")
                    or row.get("pressure_error") not in (None, "")
                    else p - exact_p
                ),
                "pressure_perturbation_error": (
                    get(row, "pressure_perturbation_error")
                    if row.get("pressure_perturbation_error") not in (None, "")
                    else math.nan
                ),
            }
        )

    def l1(key: str) -> float:
        return sum(abs(row[key]) for row in parsed) / len(parsed)

    pressure_background = pressure_from_mach(mach, gamma, density)
    p_pert_errors = [
        abs(row["pressure_perturbation_error"])
        if math.isfinite(row.get("pressure_perturbation_error", math.nan))
        else abs((row["p"] - pressure_background) - (row["exact_p"] - pressure_background))
        for row in parsed
    ]
    exact_pert_scale = (
        sum(abs(row["exact_p"] - pressure_background) for row in parsed) / len(parsed)
    )
    velocity_errors = []
    for row in parsed:
        if math.isfinite(row["velocity_error"]):
            velocity_errors.append(abs(row["velocity_error"]))
        else:
            velocity_errors.append(math.hypot(row["u_error"], row["v_error"]))
    exact_ke = sum(
        0.5 * row["exact_rho"] * (row["exact_u"] ** 2 + row["exact_v"] ** 2)
        for row in parsed
    )
    numerical_ke = sum(0.5 * row["rho"] * (row["u"] ** 2 + row["v"] ** 2) for row in parsed)
    p_l1 = sum(abs(row["p_error"]) for row in parsed) / len(parsed)
    return {
        "density_l1_error": f"{l1('rho_error'):.12e}",
        "velocity_l1_error": f"{sum(velocity_errors) / len(velocity_errors):.12e}",
        "pressure_l1_error": f"{p_l1:.12e}",
        "pressure_perturbation_l1_error": f"{sum(p_pert_errors) / len(p_pert_errors):.12e}",
        "pressure_perturbation_l1_relative_error": (
            f"{(sum(p_pert_errors) / len(p_pert_errors)) / exact_pert_scale:.12e}"
            if exact_pert_scale > 0.0
            else "nan"
        ),
        "kinetic_energy_ratio": f"{numerical_ke / exact_ke:.12e}" if exact_ke > 0.0 else "nan",
        "rho_min": f"{min(row['rho'] for row in parsed):.12e}",
        "pressure_min": f"{min(row['p'] for row in parsed):.12e}",
    }


def shell_join(command: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in command)


def command_for(
    scheme: Scheme,
    mach: float,
    n: int,
    args: argparse.Namespace,
    final_csv: Path,
) -> list[str]:
    pressure = pressure_from_mach(mach, args.gamma, args.density)
    exact_dirichlet = args.field_boundary == "exact_dirichlet"
    euler_field_boundary = "exact_dirichlet" if exact_dirichlet else "outflow"
    base = [
        str(APP),
        "inputs-ci",
        "euler.problem=gresho_vortex",
        f"euler.field_boundary={euler_field_boundary}",
        f"euler.mach={mach:.17g}",
        f"euler.pressure={pressure:.17g}",
        f"euler.density_outer={args.density:.17g}",
        "euler.vortex_center=0.0 0.0",
        f"euler.gamma={args.gamma:.17g}",
        f"euler.final_csv={final_csv}",
        f"amr.n_cell={n} {n}",
        f"amr.max_grid_size={args.max_grid_size}",
        "amr.plot_int=-1",
        "geometry.prob_lo=-0.5 -0.5",
        "geometry.prob_hi=0.5 0.5",
        "geometry.is_periodic=0 0" if exact_dirichlet else "geometry.is_periodic=1 1",
    ]

    if scheme.key == "explicit_o2_hllc":
        return base + [
            "euler.method=explicit",
            "euler.riemann=hllc",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            f"euler.cfl={args.explicit_cfl:.17g}",
            f"stop_time={args.target_time:.17g}",
            f"max_step={args.explicit_max_step}",
        ]
    if scheme.key == "explicit_o2_xie_am_hllc_p":
        return base + [
            "euler.method=explicit",
            "euler.riemann=xie_am_hllc_p",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            f"euler.cfl={args.explicit_cfl:.17g}",
            f"stop_time={args.target_time:.17g}",
            f"max_step={args.explicit_max_step}",
        ]
    if scheme.key == "imex_t1s2_bdltv20":
        return base + [
            (
                "euler.bdltv20_paper_t1_s2=gresho_exact_dirichlet_2d"
                if exact_dirichlet
                else "euler.bdltv20_paper_t1_s2=gresho_periodic_2d"
            ),
            "euler.method=imex",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            "euler.imex_picard_iterations=4",
            "euler.imex_acoustic_startup=0",
            "euler.imex_acoustic_cfl_cap=0.0",
            "euler.imex_pressure_stabilization=off",
            "euler.imex_predictor_dissipation=material",
            "euler.imex_solver_tol=1e-8",
            "euler.imex_solver_max_iter=1000",
            "euler.bdltv20_paper_pressure_solver=gmres",
            "euler.bdltv20_paper_t1_s2_dt=-1",
            f"euler.bdltv20_paper_t1_s2_max_steps={args.imex_max_step}",
            f"euler.imex_cfl={args.clean_imex_cfl:.17g}",
            f"stop_time={args.target_time:.17g}",
        ]
    raise ValueError(scheme.key)


def run_command(command: list[str], log_path: Path, timeout: float) -> tuple[int, str, str, float]:
    start = time.perf_counter()
    proc = subprocess.run(
        command,
        cwd=str(APP.parent),
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    elapsed = time.perf_counter() - start
    log_path.write_text(
        "$ "
        + shell_join(command)
        + f"\n# measured_driver_wall_time_sec={elapsed:.12g}\n\n[stdout]\n"
        + proc.stdout
        + "\n[stderr]\n"
        + proc.stderr
    )
    return proc.returncode, proc.stdout, proc.stderr, elapsed


def make_plots(rows: list[dict[str, str]], output_dir: Path) -> None:
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)
    ok_rows = [
        row
        for row in rows
        if row["status"] in {"ok", "failed_quality_gate"}
        and finite_float(row["driver_wall_time_sec"]) is not None
        and finite_float(row["pressure_perturbation_l1_relative_error"]) is not None
    ]
    if not ok_rows:
        return

    # Low-to-high Mach left-to-right, so the lower-Mach panels sit on the left.
    machs = sorted({float(row["mach"]) for row in ok_rows})
    fig, axes = plt.subplots(1, len(machs), figsize=(7.2 * len(machs), 9.5), sharey=True)
    if len(machs) == 1:
        axes = [axes]
    scheme_by_key = {scheme.key: scheme for scheme in SCHEMES}
    for ax, mach in zip(axes, machs):
        for scheme in SCHEMES:
            sub = sorted(
                [
                    row
                    for row in ok_rows
                    if row["scheme_key"] == scheme.key and abs(float(row["mach"]) - mach) < 1e-15
                ],
                key=lambda row: int(row["n"]),
            )
            if not sub:
                continue
            has_gate_fail = any(row["status"] == "failed_quality_gate" for row in sub)
            ax.scatter(
                [float(row["driver_wall_time_sec"]) for row in sub],
                [float(row["pressure_perturbation_l1_relative_error"]) for row in sub],
                color=scheme.color if not has_gate_fail else "white",
                edgecolors=scheme.color,
                marker=scheme.marker,
                linewidths=1.3 if has_gate_fail else 0.8,
                s=38,
                zorder=3,
                label=scheme.label + (" (gate fail)" if has_gate_fail else ""),
            )
            if len(sub) >= 3:
                xs = [math.log10(float(row["driver_wall_time_sec"])) for row in sub]
                ys = [
                    math.log10(float(row["pressure_perturbation_l1_relative_error"]))
                    for row in sub
                ]
                x_mean = sum(xs) / len(xs)
                y_mean = sum(ys) / len(ys)
                denom = sum((x - x_mean) ** 2 for x in xs)
                if denom > 0.0:
                    slope = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys)) / denom
                    intercept = y_mean - slope * x_mean
                    x0, x1 = min(xs), max(xs)
                    fit_x = [10.0**x0, 10.0**x1]
                    fit_y = [10.0 ** (intercept + slope * x0), 10.0 ** (intercept + slope * x1)]
                    ax.plot(
                        fit_x,
                        fit_y,
                        color=scheme.color,
                        linestyle="-",
                        linewidth=2.0,
                        alpha=0.9,
                        zorder=2,
                    )
                    ax.annotate(
                        f"fit {slope:.2f}",
                        (fit_x[-1], fit_y[-1]),
                        xytext=(5, -10),
                        textcoords="offset points",
                        fontsize=7,
                        color=scheme.color,
                    )
            for row in sub:
                ax.annotate(
                    f"N={row['n']}" + (" gate" if row["status"] == "failed_quality_gate" else ""),
                    (
                        float(row["driver_wall_time_sec"]),
                        float(row["pressure_perturbation_l1_relative_error"]),
                    ),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=7.5,
                )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(True, which="both", color="#dddddd", linewidth=0.7, alpha=0.8)
        ax.set_title(f"Gresho Mach {mach:g}")
        ax.set_xlabel("Measured driver wall time (s)")
    axes[0].set_ylabel("Pressure perturbation relative L1 error")
    axes[-1].legend(loc="best", fontsize=7.5)
    fig.suptitle("AMReX Gresho reference matrix", fontsize=12.5)
    fig.text(
        0.5,
        0.02,
        "Mach increases left-to-right. Shared vertical axis; markers are rows; solid lines are least-squares fits. Hollow/gate markers, when present, mark quality-gate failures.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    fig.savefig(fig_dir / "amrex_gresho_reference_pressure_efficiency.png", dpi=220)
    fig.savefig(
        fig_dir / "amrex_gresho_reference_pressure_efficiency_fit_only.png",
        dpi=220,
    )
    plt.close(fig)

    for mach in machs:
        fig, ax = plt.subplots(figsize=(10.5, 4.8))
        sub = [row for row in ok_rows if abs(float(row["mach"]) - mach) < 1e-15]
        labels = []
        xs = []
        values = []
        colors = []
        for n in sorted({int(row["n"]) for row in sub}):
            for scheme in SCHEMES:
                match = [row for row in sub if int(row["n"]) == n and row["scheme_key"] == scheme.key]
                if not match:
                    continue
                labels.append(
                    f"{scheme.label}\nN={n}"
                    + ("\n(gate)" if match[0]["status"] == "failed_quality_gate" else "")
                )
                xs.append(len(xs))
                values.append(float(match[0]["driver_wall_time_sec"]))
                colors.append(scheme.color)
        ax.bar(xs, values, color=colors)
        ax.set_yscale("log")
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=7.5)
        ax.set_ylabel("Measured driver wall time (s)")
        ax.set_title(f"AMReX Gresho reference wall time, Mach {mach:g}")
        ax.grid(axis="y", which="both", color="#dddddd", linewidth=0.7, alpha=0.8)
        fig.tight_layout()
        fig.savefig(fig_dir / f"amrex_gresho_reference_walltime_mach_{mach:g}.png", dpi=220)
        plt.close(fig)


def write_fit_summary(rows: list[dict[str, str]], output_dir: Path) -> None:
    fit_path = output_dir / "fit_summary.csv"
    valid_rows = [
        row
        for row in rows
        if row["status"] in {"ok", "failed_quality_gate"}
        and finite_float(row["driver_wall_time_sec"]) is not None
        and finite_float(row["pressure_perturbation_l1_relative_error"]) is not None
    ]
    out_rows: list[dict[str, str]] = []
    for mach in sorted({float(row["mach"]) for row in valid_rows}):
        for scheme in SCHEMES:
            sub = sorted(
                [
                    row
                    for row in valid_rows
                    if row["scheme_key"] == scheme.key and abs(float(row["mach"]) - mach) < 1e-15
                ],
                key=lambda row: int(row["n"]),
            )
            row_out = {
                "mach": f"{mach:.12g}",
                "scheme_key": scheme.key,
                "scheme_label": scheme.label,
                "valid_points": str(len(sub)),
                "n_values": ";".join(row["n"] for row in sub),
                "statuses": ";".join(row["status"] for row in sub),
                "fit_space": "log10(pressure_perturbation_l1_relative_error)_vs_log10(driver_wall_time_sec)",
                "fit_status": "insufficient_points",
                "slope": "",
                "intercept": "",
            }
            if len(sub) >= 3:
                xs = [math.log10(float(row["driver_wall_time_sec"])) for row in sub]
                ys = [
                    math.log10(float(row["pressure_perturbation_l1_relative_error"]))
                    for row in sub
                ]
                x_mean = sum(xs) / len(xs)
                y_mean = sum(ys) / len(ys)
                denom = sum((x - x_mean) ** 2 for x in xs)
                if denom > 0.0:
                    slope = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys)) / denom
                    intercept = y_mean - slope * x_mean
                    row_out.update(
                        {
                            "fit_status": "least_squares",
                            "slope": f"{slope:.12e}",
                            "intercept": f"{intercept:.12e}",
                        }
                    )
            out_rows.append(row_out)
    with fit_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mach",
                "scheme_key",
                "scheme_label",
                "valid_points",
                "n_values",
                "statuses",
                "fit_space",
                "fit_status",
                "slope",
                "intercept",
            ],
        )
        writer.writeheader()
        writer.writerows(out_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--machs", default="0.1,0.01,0.001")
    parser.add_argument("--grids", default="32,64")
    parser.add_argument("--target-time", type=float, default=0.4 * math.pi)
    parser.add_argument("--field-boundary", choices=["exact_dirichlet", "periodic"], default="exact_dirichlet")
    parser.add_argument(
        "--allow-non-reference-controls",
        action="store_true",
        help=(
            "Opt out of the project Gresho reference controls. By default this "
            "runner requires exact Dirichlet x/y boundaries and final_time=0.4*pi."
        ),
    )
    parser.add_argument("--clean-imex-dt", type=float, default=0.0)
    parser.add_argument("--clean-imex-steps", type=int, default=0)
    parser.add_argument("--clean-imex-cfl", type=float, default=0.15)
    parser.add_argument("--clean-imex-picard-iterations", type=int, default=2)
    parser.add_argument("--explicit-cfl", type=float, default=0.4)
    parser.add_argument("--explicit-max-step", type=int, default=200000)
    parser.add_argument("--imex-max-step", type=int, default=10000)
    parser.add_argument("--row-timeout-sec", type=float, default=300.0)
    parser.add_argument("--max-grid-size", type=int, default=64)
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument(
        "--schemes",
        default=",".join(scheme.key for scheme in SCHEMES),
        help="Comma-separated scheme keys to run. Default: all project Gresho reference classifications.",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if not args.allow_non_reference_controls:
        expected_target_time = 0.4 * math.pi
        if args.field_boundary != "exact_dirichlet":
            raise SystemExit(
                "Refusing Gresho run: use --field-boundary exact_dirichlet for "
                "project Gresho reference comparison rows, or pass "
                "--allow-non-reference-controls for a clearly labelled comparison run."
            )
        if abs(args.target_time - expected_target_time) > 1.0e-12:
            raise SystemExit(
                "Refusing Gresho run: use --target-time 0.4*pi "
                f"({expected_target_time:.17g}) for project Gresho reference comparison rows, "
                "or pass --allow-non-reference-controls for a clearly labelled "
                "comparison run."
            )
    selected_schemes = split_csv_strings(args.schemes)
    valid_schemes = {scheme.key for scheme in SCHEMES}
    unknown_schemes = [scheme for scheme in selected_schemes if scheme not in valid_schemes]
    if unknown_schemes:
        raise SystemExit(
            f"unknown --schemes entries: {','.join(unknown_schemes)}; "
            f"valid: {','.join(sorted(valid_schemes))}"
        )
    schemes = tuple(scheme for scheme in SCHEMES if scheme.key in set(selected_schemes))

    output_dir = args.output_dir
    logs_dir = output_dir / "logs"
    fields_dir = output_dir / "fields"
    command_dir = output_dir / "commands"
    for path in (logs_dir, fields_dir, command_dir):
        path.mkdir(parents=True, exist_ok=True)
    output_dir = output_dir.resolve()
    logs_dir = logs_dir.resolve()
    fields_dir = fields_dir.resolve()
    command_dir = command_dir.resolve()

    rows: list[dict[str, str]] = []
    for mach in split_csv_floats(args.machs):
        for n in split_csv_ints(args.grids):
            for scheme in schemes:
                row_id = f"gresho_m{mach:g}_{scheme.key}_n{n}".replace(".", "p")
                final_csv = fields_dir / f"{row_id}_final.csv"
                log_path = logs_dir / f"{row_id}.log"
                command_path = command_dir / f"{row_id}.txt"
                command = command_for(scheme, mach, n, args, final_csv)
                command_path.write_text(shell_join(command) + "\n")
                row = {
                    "row_id": row_id,
                    "scheme_key": scheme.key,
                    "scheme_label": scheme.label,
                    "scheme_role": scheme.role,
                    "problem": "gresho_vortex",
                    "mach": f"{mach:.12g}",
                    "n": str(n),
                    "target_time": f"{args.target_time:.12g}",
                    "clean_imex_dt": f"{args.clean_imex_dt:.12g}",
                    "clean_imex_steps": str(args.clean_imex_steps),
                    "status": "dry_run" if args.dry_run else "not_run",
                    "returncode": "",
                    "driver_wall_time_sec": "",
                    "app_wall_time_sec": "",
                    "completed_steps": "",
                    "final_time": "",
                    "time_order": "1" if scheme.key.startswith("imex_") else "explicit",
                    "space_order": "2",
                    "boundary_condition": (
                        "xy_exact_dirichlet_gresho_field" if args.field_boundary == "exact_dirichlet"
                        else "periodic_2d"
                    ),
                    "field_boundary": args.field_boundary,
                    "final_csv": str(final_csv),
                    "run_log": str(log_path),
                    "command_file": str(command_path),
                    "claim_scope": (
                        "AMReX Gresho reference matrix; "
                        "exact Dirichlet and 0.4*pi final time when field_boundary=exact_dirichlet; "
                        "not CPU/GPU evidence"
                    ),
                }
                if not args.dry_run:
                    print(f"[run] {row_id}", flush=True)
                    try:
                        returncode, stdout, stderr, elapsed = run_command(
                            command, log_path, args.row_timeout_sec
                        )
                    except subprocess.TimeoutExpired as exc:
                        log_path.write_text(
                            "$ "
                            + shell_join(command)
                            + f"\n# timed_out_after_sec={args.row_timeout_sec:.12g}\n"
                            + (exc.stdout or "")
                            + (exc.stderr or "")
                        )
                        row.update({"status": "timed_out", "returncode": "timeout"})
                        rows.append(row)
                        continue
                    text = stdout + "\n" + stderr
                    kv = parse_kv(text)
                    row.update(
                        {
                            "returncode": str(returncode),
                            "driver_wall_time_sec": f"{elapsed:.12g}",
                            "app_wall_time_sec": kv.get("wall_time_sec", ""),
                            "completed_steps": kv.get(
                                "completed_steps",
                                kv.get("bdltv20_paper_steps", ""),
                            ),
                            "final_time": kv.get(
                                "completed_time",
                                kv.get("bdltv20_paper_final_time", ""),
                            ),
                        }
                    )
                    clean_status = kv.get("bdltv20_paper_t1_s2_status", "")
                    project_status = kv.get("project_amrex_euler_compare_status", "")
                    if returncode != 0 and not (
                        final_csv.exists()
                        and clean_status == "failed"
                    ):
                        row["status"] = "failed_run"
                    elif not final_csv.exists():
                        row["status"] = "failed_missing_final_csv"
                    else:
                        if scheme.key == "imex_t1s2_bdltv20":
                            row["status"] = "ok" if clean_status == "passed" else "failed_quality_gate"
                        else:
                            row["status"] = (
                                "ok"
                                if clean_status in {"", "passed"}
                                and project_status in {
                                    "",
                                    "ok",
                                    "bdltv20_paper_t1_s2_direct_exit",
                                }
                                else "failed_quality_gate"
                            )
                        row.update(profile_metrics(final_csv, mach, args.gamma, args.density))
                        row["rho_min"] = kv.get("rho_min", row.get("rho_min", ""))
                        row["pressure_min"] = kv.get("pressure_min", row.get("pressure_min", ""))
                        row["nonfinite_count"] = kv.get(
                            "nonfinite_count", kv.get("imex_clean_contains_nan", "")
                        )
                rows.append(row)

    summary_path = output_dir / "summary.csv"
    fieldnames = sorted({key for row in rows for key in row})
    preferred = [
        "row_id",
        "status",
        "scheme_key",
        "scheme_label",
        "scheme_role",
        "problem",
        "mach",
        "n",
        "target_time",
        "final_time",
        "completed_steps",
        "driver_wall_time_sec",
        "pressure_perturbation_l1_relative_error",
        "density_l1_error",
        "velocity_l1_error",
        "pressure_l1_error",
        "kinetic_energy_ratio",
        "time_order",
        "space_order",
        "claim_scope",
        "final_csv",
        "run_log",
        "command_file",
    ]
    ordered = preferred + [key for key in fieldnames if key not in preferred]
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered)
        writer.writeheader()
        writer.writerows(rows)
    if not args.dry_run:
        write_fit_summary(rows, output_dir)
        make_plots(rows, output_dir)
    print(summary_path)


if __name__ == "__main__":
    main()
