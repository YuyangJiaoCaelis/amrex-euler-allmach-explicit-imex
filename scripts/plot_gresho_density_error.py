#!/usr/bin/env python3
"""Create Gresho error figures used in the report.

This script reads the evidence AMReX Gresho CSV only; it does not run solvers.
"""

from __future__ import annotations

import csv
import math
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "project_outputs/test_packages/gresho_vortex"
REPORT_FIGS = ROOT / "project_outputs/report_figures_cued_template"
EVIDENCE_FIGURES = ROOT / "project_outputs/project_evidence/figures"
SOURCE = PACKAGE / "source_data/amrex_gresho_efficiency_combined_data.csv"

SCHEME_ORDER = (
    "Explicit O2 HLLC",
    "Explicit O2 Low-Mach Corrected HLLC-P",
    "IMEX T1/S2 BDLTV20",
)

SCHEME_ALIASES = {
    "AMReX Explicit O2 HLLC": "Explicit O2 HLLC",
    "AMReX Explicit O2 LowMachCorrected HLLC-P": "Explicit O2 Low-Mach Corrected HLLC-P",
    "Explicit O2 HLLC": "Explicit O2 HLLC",
    "Explicit O2 LowMachCorrected HLLC-P": "Explicit O2 Low-Mach Corrected HLLC-P",
    "Explicit O2 Low-Mach Corrected HLLC-P": "Explicit O2 Low-Mach Corrected HLLC-P",
    "AMReX IMEX T1/S2": "IMEX T1/S2 BDLTV20",
    "IMEX T1/S2 BDLTV20": "IMEX T1/S2 BDLTV20",
}

COLORS = {
    "Explicit O2 HLLC": "#1f77b4",
    "Explicit O2 Low-Mach Corrected HLLC-P": "#9467bd",
    "IMEX T1/S2 BDLTV20": "#2ca02c",
}

MARKERS = {
    "Explicit O2 HLLC": "o",
    "Explicit O2 Low-Mach Corrected HLLC-P": "s",
    "IMEX T1/S2 BDLTV20": "o",
}


def read_rows() -> list[dict[str, str]]:
    with SOURCE.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row["plot_scheme"] = SCHEME_ALIASES.get(row.get("scheme_label", ""), row.get("scheme_label", ""))
    return [row for row in rows if row.get("plot_scheme") in SCHEME_ORDER]


def row_float(row: dict[str, str], *keys: str) -> float:
    for key in keys:
        text = row.get(key, "")
        if text:
            return float(text)
    return math.nan


METRICS = (
    (
        "pressure_perturbation_l1_relative_error",
        r"Pressure perturbation relative $L^1$ error",
        "amrex_gresho_pressure_error_vs_walltime_with_grid_labels_polished.png",
    ),
    (
        "velocity_l1_error",
        r"Velocity $L^1$ error",
        "amrex_gresho_velocity_error_vs_walltime_with_grid_labels_polished.png",
    ),
    (
        "density_l1_error",
        r"Density $L^1$ error",
        "amrex_gresho_density_error_vs_walltime_with_grid_labels_polished.png",
    ),
)


def make_plot(metric: str, ylabel: str, filename: str) -> Path:
    rows = read_rows()
    machs = sorted({float(row["mach"]) for row in rows})
    use_report_full_canvas = filename in {
        "amrex_gresho_pressure_error_vs_walltime_with_grid_labels_polished.png",
        "amrex_gresho_velocity_error_vs_walltime_with_grid_labels_polished.png",
    }
    figsize = (15.8, 5.35) if use_report_full_canvas else (15.4, 4.7)
    fig, axes = plt.subplots(1, len(machs), figsize=figsize, sharey=True, constrained_layout=True)
    if len(machs) == 1:
        axes = [axes]

    for ax, mach in zip(axes, machs):
        for scheme in SCHEME_ORDER:
            data = [
                row
                for row in rows
                if row["plot_scheme"] == scheme and abs(float(row["mach"]) - mach) < 1.0e-14
            ]
            if not data:
                continue
            x = np.array([row_float(row, "driver_wall_time_sec", "wall_time_sec") for row in data])
            y = np.array([row_float(row, metric) for row in data])
            n = np.array([int(row["n"]) for row in data])
            valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
            if not valid.any():
                continue
            ax.scatter(
                x[valid],
                y[valid],
                s=42,
                color=COLORS[scheme],
                marker=MARKERS[scheme],
                label=scheme,
                zorder=3,
            )
            for xx, yy, nn in zip(x[valid], y[valid], n[valid]):
                ax.annotate(
                    str(nn),
                    (xx, yy),
                    xytext=(4, 3),
                    textcoords="offset points",
                    fontsize=7,
                    color=COLORS[scheme],
                    weight="bold",
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec=COLORS[scheme], lw=0.6, alpha=0.88),
                )
            if valid.sum() >= 3:
                coeff = np.polyfit(np.log10(x[valid]), np.log10(y[valid]), 1)
                xs = np.logspace(np.log10(x[valid].min()), np.log10(x[valid].max()), 100)
                ax.plot(xs, 10.0 ** np.polyval(coeff, np.log10(xs)), color=COLORS[scheme], lw=1.7)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"Mach {mach:g}")
        ax.set_xlabel("CPU wall time (s)")
        ax.grid(True, which="both", alpha=0.25)

    axes[0].set_ylabel(ylabel)
    axes[0].legend(loc="lower left", fontsize=8, frameon=True)

    out = PACKAGE / "figures" / filename
    out.parent.mkdir(parents=True, exist_ok=True)
    if use_report_full_canvas:
        fig.savefig(out, dpi=240)
    else:
        fig.savefig(out, dpi=240, bbox_inches="tight")
    plt.close(fig)

    REPORT_FIGS.mkdir(parents=True, exist_ok=True)
    EVIDENCE_FIGURES.mkdir(parents=True, exist_ok=True)
    shutil.copy2(out, REPORT_FIGS / out.name)
    shutil.copy2(out, EVIDENCE_FIGURES / out.name)
    return out


def update_manifest(path: Path) -> None:
    manifest = PACKAGE / "MANIFEST.csv"
    if not manifest.exists():
        return
    text = manifest.read_text()
    if str(path) in text:
        return
    with manifest.open("a") as handle:
        handle.write(f"{path},report figure,polished Gresho density error plot\n")


def main() -> None:
    for metric, ylabel, filename in METRICS:
        out = make_plot(metric, ylabel, filename)
        update_manifest(out)
        print(out)


if __name__ == "__main__":
    main()
