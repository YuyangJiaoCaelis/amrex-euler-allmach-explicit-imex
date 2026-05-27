#!/usr/bin/env python3
"""Refresh Riemann efficiency plots with grid-size labels."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "project_outputs/test_packages/riemann_1d_convergence"
METRICS = PACKAGE / "source_data/riemann_1d_convergence_metrics.csv"
FIG_DIR = PACKAGE / "figures"
REPORT_EFFICIENCY_MAX_GRID = 400

CASES = [
    ("sod_toro1", "Sod / Toro 1"),
    ("lax_rp2", "Lax / BDLTV20 RP2"),
    ("toro3_strong_shock", "Toro 3"),
    ("toro4_strong_shock", "Toro 4"),
]

SCHEMES = [
    ("Explicit O2 HLLC", "Explicit O2 HLLC", "#1f77b4"),
    ("Explicit O2 Low-Mach HLLC-P", "Explicit O2 LowMachCorrected HLLC-P", "#9467bd"),
    ("IMEX T1/S2 BDLTV20", "IMEX T1/S2 BDLTV20", "#2ca02c"),
]

QUANTITIES = [
    ("rho", "Density", "riemann_1d_density_error_vs_walltime_primary.png", "riemann_1d_density_error_vs_walltime_with_grid_labels.png"),
    ("u", "Velocity", "riemann_1d_velocity_error_vs_walltime_primary.png", "riemann_1d_velocity_error_vs_walltime_with_grid_labels.png"),
    ("p", "Pressure", "riemann_1d_pressure_error_vs_walltime_primary.png", "riemann_1d_pressure_error_vs_walltime_with_grid_labels.png"),
]


def main() -> None:
    rows = list(csv.DictReader(METRICS.open()))
    verify_current_scheme_sources(rows)
    report_rows = [row for row in rows if int(float(row["ncell"])) <= REPORT_EFFICIENCY_MAX_GRID]
    for quantity, label, filename, labelled_filename in QUANTITIES:
        plot_quantity(report_rows, quantity, label, FIG_DIR / filename)
        plot_quantity(report_rows, quantity, label, FIG_DIR / labelled_filename)
    refresh_manifest()


def verify_current_scheme_sources(rows: list[dict[str, str]]) -> None:
    schemes = sorted({row["scheme"] for row in rows})
    expected = sorted(scheme[0] for scheme in SCHEMES)
    if schemes != expected:
        raise SystemExit(f"Unexpected Riemann schemes: {schemes}; expected {expected}")
    if {row["boundary"] for row in rows} != {"xy_exact_dirichlet"}:
        raise SystemExit("Riemann efficiency rows are not all xy_exact_dirichlet.")
    if any(row["scheme"] == "IMEX T1/S2 BDLTV20" and row["bdltv20_paper_status"] != "passed" for row in rows):
        raise SystemExit("At least one IMEX row is not bdltv20_paper_t1_s2_status=passed.")


def plot_quantity(rows: list[dict[str, str]], quantity: str, ylabel: str, output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.2), constrained_layout=True)
    for ax, (case_key, case_label) in zip(axes.ravel(), CASES):
        case_rows = [row for row in rows if row["case"] == case_key and row["quantity"] == quantity]
        for source_label, plot_label, color in SCHEMES:
            data = sorted(
                [row for row in case_rows if row["scheme"] == source_label],
                key=lambda row: float(row["wall_time_sec"]),
            )
            if not data:
                continue
            xs = [float(row["wall_time_sec"]) for row in data]
            ys = [float(row["l1_error"]) for row in data]
            grids = [row["ncell"] for row in data]
            ax.scatter(xs, ys, color=color, s=38, label=plot_label, zorder=3)
            annotate_grid_sizes(ax, xs, ys, grids, color)
            positive = [
                (x, y)
                for x, y in zip(xs, ys)
                if x > 0.0 and y > 0.0 and math.isfinite(x) and math.isfinite(y)
            ]
            if len(positive) >= 2:
                log_x = [math.log(x) for x, _ in positive]
                log_y = [math.log(y) for _, y in positive]
                slope, intercept = least_squares(log_x, log_y)
                fit_x = [min(x for x, _ in positive), max(x for x, _ in positive)]
                fit_y = [math.exp(intercept + slope * math.log(x)) for x in fit_x]
                ax.plot(fit_x, fit_y, color=color, linewidth=1.6)
        ax.set_title(case_label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Wall time (s)")
        ax.set_ylabel(f"{ylabel} L1 error")
        ax.grid(True, which="both", alpha=0.22)
    axes.ravel()[0].legend(fontsize=8, frameon=True)
    fig.suptitle(
        f"AMReX 1D Riemann efficiency: {ylabel} error vs wall time\n"
        "markers are grid runs; solid lines are log-log least-square fits",
        fontsize=14,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240)
    plt.close(fig)
    print(output)


def annotate_grid_sizes(ax, xs: list[float], ys: list[float], grids: list[str], color: str) -> None:
    for index, (x, y, grid) in enumerate(zip(xs, ys, grids)):
        dx = 8 if index % 2 == 0 else -16
        dy = 9 if index % 2 == 0 else -13
        ax.annotate(
            grid,
            xy=(x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            color=color,
            fontsize=9.0,
            fontweight="bold",
            ha="left" if dx > 0 else "right",
            va="bottom" if dy > 0 else "top",
            bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": color, "linewidth": 0.6, "alpha": 0.86},
            zorder=6,
            clip_on=False,
        )


def refresh_manifest() -> None:
    manifest = PACKAGE / "MANIFEST.csv"
    rows: list[dict[str, str]] = []
    for path in sorted(PACKAGE.rglob("*")):
        if path.is_file() and path.name != "MANIFEST.csv":
            rows.append(
                {
                    "path": str(path),
                    "role": classify_path(path),
                    "source": "run_riemann_1d_convergence.py; plot_riemann_1d_efficiency.py",
                }
            )
    with manifest.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["path", "role", "source"])
        writer.writeheader()
        writer.writerows(rows)


def classify_path(path: Path) -> str:
    if path.suffix == ".png":
        return "report figure"
    if path.name.endswith(".csv"):
        return "source data"
    if path.suffix == ".log":
        return "run log"
    if path.name.endswith(".command.txt"):
        return "reproduction command"
    return "documentation"


def least_squares(xs: list[float], ys: list[float]) -> tuple[float, float]:
    xbar = sum(xs) / len(xs)
    ybar = sum(ys) / len(ys)
    denom = sum((x - xbar) ** 2 for x in xs)
    if denom == 0.0:
        return 0.0, ybar
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / denom
    return slope, ybar - slope * xbar


if __name__ == "__main__":
    main()
