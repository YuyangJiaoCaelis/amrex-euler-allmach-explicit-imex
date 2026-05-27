#!/usr/bin/env python3
"""Refresh 1D Riemann figures from the current AMReX runs."""

from __future__ import annotations

import csv
import shutil
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt

from riemann_exact import RiemannState, sample, solve_star_region


ROOT = Path(__file__).resolve().parents[1]
CONV = ROOT / "project_outputs/test_packages/riemann_1d_convergence"
RUNS = ROOT / "results/amrex/project_riemann_1d_convergence_current_t1s2_2026-05-17"
PACKAGE = ROOT / "project_outputs/test_packages/riemann_1d"
EVIDENCE = ROOT / "project_outputs/project_evidence"

GAMMA = 1.4


@dataclass(frozen=True)
class Case:
    key: str
    label: str
    toro_test: int
    final_time: float
    left: RiemannState
    right: RiemannState


@dataclass(frozen=True)
class Scheme:
    key: str
    label: str
    color: str


CASES = [
    Case("sod_toro1", "Sod / Toro 1", 1, 0.2, RiemannState(1.0, 0.0, 1.0), RiemannState(0.125, 0.0, 0.1)),
    Case("lax_rp2", "Lax / BDLTV20 RP2", 93, 0.14, RiemannState(0.445, 1.698, 3.528), RiemannState(0.5, 0.0, 0.571)),
    Case("toro3_strong_shock", "Toro 3", 3, 0.012, RiemannState(1.0, 0.0, 1000.0), RiemannState(1.0, 0.0, 0.01)),
    Case("toro4_strong_shock", "Toro 4", 4, 0.035, RiemannState(1.0, 0.0, 0.01), RiemannState(1.0, 0.0, 100.0)),
]

SCHEMES = [
    Scheme("explicit_o2_hllc", "Explicit O2 HLLC", "#1f77b4"),
    Scheme("explicit_o2_lowmach_hllcp", "Explicit O2 Low-Mach Corrected HLLC-P", "#9467bd"),
    Scheme("imex_t1s2_bdltv20", "IMEX T1/S2 BDLTV20", "#2ca02c"),
]


def main() -> None:
    for path in [
        PACKAGE / "figures",
        PACKAGE / "source_data/fields",
        PACKAGE / "source_data/logs",
        PACKAGE / "source_data/commands",
        EVIDENCE / "figures",
        EVIDENCE / "data",
    ]:
        path.mkdir(parents=True, exist_ok=True)

    copy_current_sources()
    plot_primary_overlay(PACKAGE / "figures/riemann_1d_primary_overlays_full_domain.png")
    plot_toro34_imex_refinement(PACKAGE / "figures/riemann_1d_toro34_imex_refinement_wave_zoom.png")
    write_primary_error_summary()
    refresh_readme()

    for name in [
        "riemann_1d_primary_overlays_full_domain.png",
        "riemann_1d_toro34_imex_refinement_wave_zoom.png",
    ]:
        shutil.copy2(PACKAGE / "figures" / name, EVIDENCE / "figures" / name)
    for name in [
        "riemann_1d_density_error_vs_walltime_with_grid_labels.png",
        "riemann_1d_velocity_error_vs_walltime_with_grid_labels.png",
        "riemann_1d_pressure_error_vs_walltime_with_grid_labels.png",
    ]:
        shutil.copy2(CONV / "figures" / name, EVIDENCE / "figures" / name)
    for name in [
        "riemann_1d_convergence_metrics.csv",
        "riemann_1d_observed_order_summary.csv",
        "run_summary.csv",
    ]:
        shutil.copy2(CONV / "source_data" / name, EVIDENCE / "data" / name)

    print(PACKAGE / "figures/riemann_1d_primary_overlays_full_domain.png")
    print(PACKAGE / "figures/riemann_1d_toro34_imex_refinement_wave_zoom.png")


def copy_current_sources() -> None:
    for case in CASES:
        for scheme in SCHEMES:
            stem = f"{case.key}_{scheme.key}_n400"
            copy_if_exists(CONV / "source_data/fields" / f"{stem}.csv", PACKAGE / "source_data/fields" / f"{stem}.csv")
            copy_if_exists(CONV / "source_data/logs" / f"{stem}.log", PACKAGE / "source_data/logs" / f"{stem}.log")
            copy_if_exists(CONV / "source_data/commands" / f"{stem}.command.txt", PACKAGE / "source_data/commands" / f"{stem}.command.txt")
    for case_key in ["toro3_strong_shock", "toro4_strong_shock"]:
        for grid in [100, 200]:
            stem = f"{case_key}_imex_t1s2_bdltv20_n{grid}"
            copy_if_exists(CONV / "source_data/fields" / f"{stem}.csv", PACKAGE / "source_data/fields" / f"{stem}.csv")
            copy_if_exists(CONV / "source_data/logs" / f"{stem}.log", PACKAGE / "source_data/logs" / f"{stem}.log")
            copy_if_exists(CONV / "source_data/commands" / f"{stem}.command.txt", PACKAGE / "source_data/commands" / f"{stem}.command.txt")
        stem = f"{case_key}_imex_t1s2_bdltv20_n800"
        copy_if_exists(RUNS / "fields" / f"{stem}.csv", PACKAGE / "source_data/fields" / f"{stem}.csv")
        copy_if_exists(RUNS / "logs" / f"{stem}.log", PACKAGE / "source_data/logs" / f"{stem}.log")
        copy_if_exists(RUNS / "commands" / f"{stem}.command.txt", PACKAGE / "source_data/commands" / f"{stem}.command.txt")


def copy_if_exists(src: Path, dst: Path) -> None:
    src = restore_canonical_artifact(src)
    shutil.copy2(src, dst)


def restore_canonical_artifact(path: Path) -> Path:
    if path.exists():
        return path

    matching_paths = sorted(path.parent.glob(f"{path.stem} *{path.suffix}"))
    if len(matching_paths) == 1:
        matching_paths[0].replace(path)
        return path

    raise SystemExit(f"Missing required source: {path}")


def profile_path(case: Case, scheme: Scheme, n: int = 400) -> Path:
    return PACKAGE / "source_data/fields" / f"{case.key}_{scheme.key}_n{n}.csv"


def read_profile(path: Path) -> dict[str, list[float]]:
    with path.open(newline="") as f:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(f)]
    y_values = sorted({row["y"] for row in rows})
    y_mid = 0.5 * (y_values[0] + y_values[-1])
    y_pick = min(y_values, key=lambda value: abs(value - y_mid))
    line = sorted((row for row in rows if row["y"] == y_pick), key=lambda row: row["x"])
    return {
        "x": [row["x"] for row in line],
        "rho": [row["rho"] for row in line],
        "u": [row["u"] for row in line],
        "p": [row["p"] for row in line],
    }


def exact_profile(case: Case, xs: list[float]) -> dict[str, list[float]]:
    star = solve_star_region(case.left, case.right, GAMMA)
    states = [sample(x, case.final_time, case.left, case.right, GAMMA, 0.5, star) for x in xs]
    return {
        "rho": [state.rho for state in states],
        "u": [state.u for state in states],
        "p": [state.p for state in states],
    }


def plot_primary_overlay(output: Path) -> None:
    quantities = [("rho", "Density"), ("u", "Velocity u"), ("p", "Pressure")]
    fig, axes = plt.subplots(len(CASES), len(quantities), figsize=(16.5, 13.4), constrained_layout=True)
    for row_index, case in enumerate(CASES):
        reference = read_profile(profile_path(case, SCHEMES[0]))
        exact = exact_profile(case, reference["x"])
        for col_index, (quantity, label) in enumerate(quantities):
            ax = axes[row_index][col_index]
            ax.plot(reference["x"], exact[quantity], color="black", linestyle="--", linewidth=1.7, label="Exact Riemann")
            for scheme in SCHEMES:
                profile = read_profile(profile_path(case, scheme))
                ax.plot(profile["x"], profile[quantity], color=scheme.color, linewidth=1.35, label=scheme.label)
            ax.set_title(label if row_index == 0 else "")
            ax.set_xlabel("x")
            ax.set_ylabel(case.label if col_index == 0 else label)
            ax.grid(True, alpha=0.25)
            if row_index == 0 and col_index == 0:
                ax.legend(fontsize=8, frameon=True)
    fig.suptitle("AMReX exact Riemann overlays for current report schemes", fontsize=16)
    fig.savefig(output, dpi=240)
    plt.close(fig)


def plot_toro34_imex_refinement(output: Path) -> None:
    cases = [CASES[2], CASES[3]]
    quantities = [("rho", "Density"), ("u", "Velocity u"), ("p", "Pressure")]
    grids = [100, 200, 400, 800]
    colors = ["#9ecae1", "#6baed6", "#3182bd", "#08519c"]
    fig, axes = plt.subplots(4, 3, figsize=(18.5, 15.0), constrained_layout=True)
    for case_block, case in enumerate(cases):
        row_top = case_block * 2
        row_err = row_top + 1
        exact_ref = read_profile(profile_path(case, SCHEMES[2], 400))
        exact = exact_profile(case, exact_ref["x"])
        for col_index, (quantity, label) in enumerate(quantities):
            ax = axes[row_top][col_index]
            ax_err = axes[row_err][col_index]
            ax.plot(exact_ref["x"], exact[quantity], color="black", linestyle="--", linewidth=1.9, label="Exact Riemann")
            ax_err.axhline(0.0, color="0.25", linewidth=0.9)
            for grid, color in zip(grids, colors):
                profile = read_profile(profile_path(case, SCHEMES[2], grid))
                exact_grid = exact_profile(case, profile["x"])
                ax.plot(profile["x"], profile[quantity], color=color, linewidth=1.5, label=f"n={grid}")
                err = [a - b for a, b in zip(profile[quantity], exact_grid[quantity])]
                ax_err.plot(profile["x"], err, color=color, linewidth=1.35, label=f"n={grid}")
            ax.set_title(f"{case.label} {label}")
            ax.set_xlabel("x")
            ax.set_ylabel(label)
            ax.grid(True, alpha=0.24)
            ax_err.set_title(f"{case.label} pointwise error")
            ax_err.set_xlabel("x")
            ax_err.set_ylabel(f"{label}: numerical - exact")
            ax_err.grid(True, alpha=0.24)
            if case_block == 0 and col_index == 0:
                ax.legend(fontsize=8)
    fig.suptitle("IMEX T1/S2 BDLTV20: Toro 3/4 refinement against exact Riemann", fontsize=18)
    fig.savefig(output, dpi=240)
    plt.close(fig)


def write_primary_error_summary() -> None:
    rows: list[dict[str, str]] = []
    for case in CASES:
        for scheme in SCHEMES:
            profile = read_profile(profile_path(case, scheme))
            exact = exact_profile(case, profile["x"])
            for quantity in ["rho", "u", "p"]:
                errors = [abs(a - b) for a, b in zip(profile[quantity], exact[quantity])]
                rows.append({
                    "case": case.key,
                    "case_label": case.label,
                    "scheme": scheme.label,
                    "ncell": "400",
                    "quantity": quantity,
                    "l1_error": f"{sum(errors) / len(errors):.12g}",
                    "linf_error": f"{max(errors):.12g}",
                    "boundary": "xy_exact_dirichlet",
                    "exact_solution_source": "scripts/riemann_exact.py",
                })
    with (PACKAGE / "riemann_error_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def refresh_readme() -> None:
    text = """# MPhil project AMReX 1D Riemann Package

This package has been regenerated from the current AMReX `euler_compare` rows.

## Schemes

- `Explicit O2 HLLC`
- `Explicit O2 Low-Mach Corrected HLLC-P` (`euler.riemann=xie_am_hllc_p`)
- `IMEX T1/S2 BDLTV20`

## Controls

- Domain: `x in [0,1]`, thin 2D AMReX strip with `amr.n_cell=N 5`.
- Boundaries: `geometry.is_periodic=0 0`, `euler.field_boundary=exact_dirichlet`.
- Riemann cases: Sod/Toro 1 at `t=0.2`, Lax/BDLTV20 RP2 at `t=0.14`, Toro 3 at `t=0.012`, and Toro 4 at `t=0.035`.
- IMEX rows use `euler.bdltv20_paper_t1_s2=toro_xy_exact_dirichlet_2d`, `spatial_order=2`, `slope_limiter=minmod`, and four Picard iterations.

## Report Figures

- `figures/riemann_1d_primary_overlays_full_domain.png`
- `figures/riemann_1d_toro34_imex_refinement_wave_zoom.png`

The low-Mach corrected explicit rows use the source-backed Xie AM-HLLC-P path.
"""
    (PACKAGE / "README.md").write_text(text)


if __name__ == "__main__":
    main()
