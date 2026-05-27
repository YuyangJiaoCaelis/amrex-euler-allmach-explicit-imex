#!/usr/bin/env python3
"""Run and plot AMReX advection-blob efficiency rows with fit-only lines."""

from __future__ import annotations

import csv
import math
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex"
INPUT = ROOT / "amrex/apps/euler_compare/inputs-ci"
OUT = ROOT / "results/amrex/advection_blob_explicit_imex_efficiency_refresh_2026-05-17"
GRIDS = (32, 40, 48, 56, 64, 80, 96, 112, 128)


@dataclass(frozen=True)
class Scheme:
    name: str
    selector: str
    color: str
    marker: str
    imex: bool
    args: tuple[str, ...]


SCHEMES = (
    Scheme(
        "AMReX Explicit O2 HLLC",
        "explicit_o2_hllc",
        "#1f77b4",
        "o",
        False,
        (
            "euler.method=explicit",
            "euler.problem=advection_blob",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            "euler.riemann=hllc",
            "euler.cfl=0.4",
            "max_step=200000",
            "stop_time=0.5",
        ),
    ),
    Scheme(
        "AMReX Explicit O2 LowMachCorrected HLLC-P",
        "explicit_o2_lowmachcorrected_hllcp",
        "#9467bd",
        "s",
        False,
        (
            "euler.method=explicit",
            "euler.problem=advection_blob",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            "euler.riemann=xie_am_hllc_p",
            "euler.cfl=0.4",
            "max_step=200000",
            "stop_time=0.5",
        ),
    ),
    Scheme(
        "IMEX T1/S2 BDLTV20",
        "imex_bdltv20_paper_t1s2",
        "#2ca02c",
        "^",
        True,
        (
            "euler.method=imex",
            "euler.problem=advection_blob",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            "euler.bdltv20_paper_t1_s2=advection_blob_periodic_2d",
            "euler.imex_picard_iterations=4",
            "euler.imex_acoustic_startup=0",
            "euler.imex_acoustic_cfl_cap=0.0",
            "euler.imex_pressure_stabilization=off",
            "euler.imex_predictor_dissipation=material",
            "euler.imex_solver_tol=1e-8",
            "euler.imex_solver_max_iter=1000",
            "euler.bdltv20_paper_pressure_solver=gmres",
            "euler.bdltv20_paper_t1_s2_dt=-1",
            "euler.bdltv20_paper_t1_s2_max_steps=200000",
            "euler.imex_cfl=0.4",
            "stop_time=0.5",
        ),
    ),
)


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    for sub in ("commands", "logs", "fields", "figures"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    for n in GRIDS:
        for scheme in SCHEMES:
            rows.append(run_or_read(n, scheme))

    write_csv(OUT / "advection_blob_efficiency_summary.csv", rows)
    plot_efficiency(
        rows,
        "velocity_l1_error",
        "velocity L1 error at t=0.5",
        OUT / "figures/advection_blob_velocity_error_vs_walltime_fit_only.png",
    )
    plot_efficiency(
        rows,
        "pressure_l1_error",
        "pressure L1 error at t=0.5",
        OUT / "figures/advection_blob_pressure_error_vs_walltime_fit_only.png",
    )
    plot_velocity_pressure_combined(rows)
    plot_density_reference(rows)
    write_readme(rows)
    print(OUT)


def run_or_read(n: int, scheme: Scheme) -> dict[str, str]:
    stem = f"advection_blob_{scheme.selector}_n{n}_t0p5"
    csv_path = OUT / "fields" / f"{stem}.csv"
    log_path = OUT / "logs" / f"{stem}.log"
    cmd_path = OUT / "commands" / f"{stem}.txt"
    cmd = [
        str(EXE),
        str(INPUT),
        f"amr.n_cell={n} {n}",
        f"amr.max_grid_size={n}",
        "amr.plot_int=-1",
        "geometry.is_periodic=1 1",
        "geometry.prob_lo=0.0 0.0",
        "geometry.prob_hi=1.0 1.0",
        *scheme.args,
        f"euler.final_csv={csv_path}",
    ]
    if not csv_path.exists() or not log_path.exists():
        cmd_path.write_text(" ".join(subprocess.list2cmdline([part]) for part in cmd) + "\n")
        start = time.perf_counter()
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        wall = time.perf_counter() - start
        log_path.write_text(proc.stdout)
        (OUT / "logs" / f"{stem}.driver_wall_time").write_text(f"{wall:.12g}\n")
        print(f"{stem} rc={proc.returncode} wall={wall:.3f}s")
    else:
        wall = float((OUT / "logs" / f"{stem}.driver_wall_time").read_text().strip())
        print(f"{stem} cached wall={wall:.3f}s")

    log = parse_log(log_path)
    errors = summarize_field(csv_path, scheme.imex)
    if scheme.imex:
        status = log.get("bdltv20_paper_t1_s2_status", "")
        completed_steps = log.get("bdltv20_paper_steps", log.get("imex_clean_evolving_steps", ""))
        completed_time = log.get("bdltv20_paper_final_time", log.get("imex_clean_evolving_final_time", ""))
    else:
        status = log.get("project_amrex_euler_compare_status", "")
        completed_steps = log.get("completed_steps", "")
        completed_time = log.get("completed_time", "")
    return {
        "scheme": scheme.name,
        "selector": scheme.selector,
        "n": str(n),
        "status": status,
        "driver_wall_time_sec": f"{wall:.12g}",
        "emitted_wall_time_sec": log.get("wall_time_sec", ""),
        "completed_steps": completed_steps,
        "completed_time": completed_time,
        "density_l1_error": f"{errors['density_l1_error']:.12g}",
        "density_linf_error": f"{errors['density_linf_error']:.12g}",
        "velocity_l1_error": f"{errors['velocity_l1_error']:.12g}",
        "velocity_linf_error": f"{errors['velocity_linf_error']:.12g}",
        "pressure_l1_error": f"{errors['pressure_l1_error']:.12g}",
        "pressure_linf_error": f"{errors['pressure_linf_error']:.12g}",
        "source_csv": str(csv_path),
        "source_log": str(log_path),
        "source_command": str(cmd_path),
    }


def parse_log(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(errors="replace") as f:
        for line in f:
            text = line.strip()
            if "=" in text:
                key, value = text.split("=", 1)
                values[key.strip()] = value.strip()
    return values


def summarize_field(path: Path, imex: bool) -> dict[str, float]:
    density: list[float] = []
    velocity: list[float] = []
    pressure: list[float] = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if "density_error" in row:
                density.append(abs(float(row["density_error"])))
                velocity.append(abs(float(row["velocity_error"])))
                pressure.append(abs(float(row["pressure_error"])))
            elif "rho_error" in row:
                density.append(abs(float(row["rho_error"])))
                velocity.append(math.hypot(float(row["u_error"]), float(row["v_error"])))
                pressure.append(abs(float(row["p_error"])))
            else:
                density.append(abs(float(row["rho"]) - float(row["exact_rho"])))
                velocity.append(
                    math.hypot(float(row["u"]) - float(row["exact_u"]), float(row["v"]) - float(row["exact_v"]))
                )
                pressure.append(abs(float(row["p"]) - float(row["exact_p"])))
    return {
        "density_l1_error": float(np.mean(density)),
        "density_linf_error": float(np.max(density)),
        "velocity_l1_error": float(np.mean(velocity)),
        "velocity_linf_error": float(np.max(velocity)),
        "pressure_l1_error": float(np.mean(pressure)),
        "pressure_linf_error": float(np.max(pressure)),
    }


def plot_efficiency(rows: list[dict[str, str]], metric: str, ylabel: str, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.8), constrained_layout=True)
    fit_rows: list[dict[str, str]] = []
    for scheme in SCHEMES:
        data = scheme_rows(rows, scheme.name)
        x = np.array([float(row["driver_wall_time_sec"]) for row in data])
        y = np.array([float(row[metric]) for row in data])
        n = np.array([int(row["n"]) for row in data])
        valid = (x > 0.0) & (y > 0.0) & np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[valid], y[valid], color=scheme.color, marker=scheme.marker, s=48, label=scheme.name)
        for xx, yy, nn in zip(x[valid], y[valid], n[valid]):
            ax.annotate(str(nn), (xx, yy), xytext=(4, 3), textcoords="offset points", fontsize=7)
        if valid.sum() >= 3:
            coeff = np.polyfit(np.log10(x[valid]), np.log10(y[valid]), 1)
            fit_x = np.linspace(x[valid].min(), x[valid].max(), 240)
            fit_y = 10 ** (coeff[1] + coeff[0] * np.log10(fit_x))
            ax.plot(
                fit_x,
                fit_y,
                color=scheme.color,
                linewidth=2.0,
                label=f"{scheme.name} least-squares slope {coeff[0]:.2f}",
            )
            fit_rows.append(
                {
                    "metric": metric,
                    "scheme": scheme.name,
                    "fit_space": f"log10({metric})_vs_log10(driver_wall_time_sec)",
                    "point_count": str(int(valid.sum())),
                    "slope": f"{coeff[0]:.12g}",
                    "intercept": f"{coeff[1]:.12g}",
                }
            )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("driver wall time (s)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"Advection blob {ylabel} vs wall time")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7)
    fig.savefig(output, dpi=220)
    plt.close(fig)
    if fit_rows:
        fit_path = OUT / f"{metric}_fit_summary.csv"
        write_csv(fit_path, fit_rows)


def plot_velocity_pressure_combined(rows: list[dict[str, str]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(15.0, 5.8), constrained_layout=True)
    for ax, metric, ylabel in (
        (axes[0], "velocity_l1_error", "velocity L1 error at t=0.5"),
        (axes[1], "pressure_l1_error", "pressure L1 error at t=0.5"),
    ):
        for scheme in SCHEMES:
            data = scheme_rows(rows, scheme.name)
            x = np.array([float(row["driver_wall_time_sec"]) for row in data])
            y = np.array([float(row[metric]) for row in data])
            n = np.array([int(row["n"]) for row in data])
            valid = (x > 0.0) & (y > 0.0) & np.isfinite(x) & np.isfinite(y)
            ax.scatter(x[valid], y[valid], color=scheme.color, marker=scheme.marker, s=45, label=scheme.name)
            for xx, yy, nn in zip(x[valid], y[valid], n[valid]):
                ax.annotate(str(nn), (xx, yy), xytext=(4, 3), textcoords="offset points", fontsize=7)
            if valid.sum() >= 3:
                coeff = np.polyfit(np.log10(x[valid]), np.log10(y[valid]), 1)
                fit_x = np.linspace(x[valid].min(), x[valid].max(), 240)
                fit_y = 10 ** (coeff[1] + coeff[0] * np.log10(fit_x))
                ax.plot(fit_x, fit_y, color=scheme.color, linewidth=2.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("driver wall time (s)")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.25)
    axes[0].legend(fontsize=7)
    fig.suptitle("Advection blob velocity/pressure efficiency: markers plus least-squares fits only")
    fig.savefig(OUT / "figures/advection_blob_velocity_pressure_error_vs_walltime_fit_only.png", dpi=220)
    plt.close(fig)


def plot_density_reference(rows: list[dict[str, str]]) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.8), constrained_layout=True)
    for scheme in SCHEMES:
        data = scheme_rows(rows, scheme.name)
        x = np.array([float(row["driver_wall_time_sec"]) for row in data])
        y = np.array([float(row["density_l1_error"]) for row in data])
        n = np.array([int(row["n"]) for row in data])
        valid = (x > 0.0) & (y > 0.0) & np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[valid], y[valid], color=scheme.color, marker=scheme.marker, s=48, label=scheme.name)
        for xx, yy, nn in zip(x[valid], y[valid], n[valid]):
            ax.annotate(str(nn), (xx, yy), xytext=(4, 3), textcoords="offset points", fontsize=7)
        if valid.sum() >= 3:
            coeff = np.polyfit(np.log10(x[valid]), np.log10(y[valid]), 1)
            fit_x = np.linspace(x[valid].min(), x[valid].max(), 240)
            fit_y = 10 ** (coeff[1] + coeff[0] * np.log10(fit_x))
            ax.plot(fit_x, fit_y, color=scheme.color, linewidth=2.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("driver wall time (s)")
    ax.set_ylabel("density L1 error at t=0.5")
    ax.set_title("Advection blob density efficiency reference")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7)
    fig.savefig(OUT / "figures/advection_blob_density_error_vs_walltime_fit_only.png", dpi=220)
    plt.close(fig)


def scheme_rows(rows: list[dict[str, str]], scheme_name: str) -> list[dict[str, str]]:
    return sorted((row for row in rows if row["scheme"] == scheme_name), key=lambda row: int(row["n"]))


def write_readme(rows: list[dict[str, str]]) -> None:
    lines = [
        "# AMReX Advection-Blob Efficiency Refresh",
        "",
        "Grid ladder: `32,40,48,56,64,80,96,112,128`.",
        "",
        "Figures use scatter markers for data and solid least-squares fit lines only. They do not draw dot-to-dot connecting lines.",
        "",
        "## Figures",
        "",
        "- `figures/advection_blob_velocity_pressure_error_vs_walltime_fit_only.png`",
        "- `figures/advection_blob_velocity_error_vs_walltime_fit_only.png`",
        "- `figures/advection_blob_pressure_error_vs_walltime_fit_only.png`",
        "- `figures/advection_blob_density_error_vs_walltime_fit_only.png`",
        "",
        "## Caveat",
        "",
        "For the advection blob, exact velocity and pressure are constant. Velocity and pressure errors therefore mostly measure preservation/roundoff, not the interface-transport accuracy. Density L1 remains the physically meaningful blob-shape error.",
        "",
        "| scheme | n | wall time (s) | density L1 | velocity L1 | pressure L1 | status |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in sorted(rows, key=lambda item: (int(item["n"]), item["scheme"])):
        lines.append(
            f"| {row['scheme']} | {row['n']} | {float(row['driver_wall_time_sec']):.4g} | "
            f"{float(row['density_l1_error']):.4g} | {float(row['velocity_l1_error']):.4g} | "
            f"{float(row['pressure_l1_error']):.4g} | {row['status']} |"
        )
    (OUT / "README.md").write_text("\n".join(lines) + "\n")


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
