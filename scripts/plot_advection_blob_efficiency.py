#!/usr/bin/env python3
"""Create advection-blob density efficiency figures used in the report."""

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
SOURCE = ROOT / "results/amrex/advection_blob_explicit_imex_efficiency_refresh_2026-05-17/advection_blob_efficiency_summary.csv"
PACKAGE = ROOT / "project_outputs/test_packages/advection_blob"
REPORT_FIGS = ROOT / "project_outputs/report_figures_cued_template"
EVIDENCE_FIGS = ROOT / "project_outputs/project_evidence/figures"
PACKAGE_SOURCE = PACKAGE / "source_data/advection_blob_periodic_efficiency_summary.csv"

SCHEME_ORDER = (
    "Explicit O2 HLLC",
    "Explicit O2 Low-Mach HLLC-P",
    "IMEX T1/S2 BDLTV20",
)

SCHEME_ALIASES = {
    "AMReX Explicit O2 HLLC": "Explicit O2 HLLC",
    "AMReX Explicit O2 LowMachCorrected HLLC-P": "Explicit O2 Low-Mach HLLC-P",
    "Explicit O2 HLLC": "Explicit O2 HLLC",
    "Explicit O2 LowMachCorrected HLLC-P": "Explicit O2 Low-Mach HLLC-P",
    "Explicit O2 Low-Mach Corrected HLLC-P": "Explicit O2 Low-Mach HLLC-P",
    "Explicit O2 Low-Mach HLLC-P": "Explicit O2 Low-Mach HLLC-P",
    "AMReX IMEX T1/S2": "IMEX T1/S2 BDLTV20",
    "IMEX T1/S2 BDLTV20": "IMEX T1/S2 BDLTV20",
}

COLORS = {
    "Explicit O2 HLLC": "#1f77b4",
    "Explicit O2 Low-Mach HLLC-P": "#9467bd",
    "IMEX T1/S2 BDLTV20": "#2ca02c",
}

MARKERS = {
    "Explicit O2 HLLC": "o",
    "Explicit O2 Low-Mach HLLC-P": "s",
    "IMEX T1/S2 BDLTV20": "^",
}


def read_rows() -> list[dict[str, str]]:
    with SOURCE.open(newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]
    for row in rows:
        row["plot_scheme"] = SCHEME_ALIASES.get(row.get("scheme", ""), row.get("scheme", ""))
    return [row for row in rows if row["plot_scheme"] in SCHEME_ORDER]


def row_float(row: dict[str, str], key: str) -> float:
    text = row.get(key, "")
    return float(text) if text else math.nan


def save_to_report_locations(path: Path, report_copy: bool) -> None:
    EVIDENCE_FIGS.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, EVIDENCE_FIGS / path.name)
    if report_copy:
        REPORT_FIGS.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, REPORT_FIGS / path.name)


def plot_density_walltime(rows: list[dict[str, str]]) -> Path:
    fig, ax = plt.subplots(figsize=(9.71, 6.93), constrained_layout=True)
    for scheme in SCHEME_ORDER:
        data = sorted(
            [row for row in rows if row["plot_scheme"] == scheme],
            key=lambda row: int(float(row["n"])),
        )
        x = np.array([row_float(row, "driver_wall_time_sec") for row in data])
        y = np.array([row_float(row, "density_l1_error") for row in data])
        n = np.array([int(float(row["n"])) for row in data])
        valid = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
        ax.scatter(x[valid], y[valid], color=COLORS[scheme], marker=MARKERS[scheme], s=48, label=scheme)
        for xx, yy, nn in zip(x[valid], y[valid], n[valid]):
            ax.annotate(str(nn), (xx, yy), xytext=(4, 3), textcoords="offset points", fontsize=7, color=COLORS[scheme])
        if np.count_nonzero(valid) >= 3:
            coeff = np.polyfit(np.log10(x[valid]), np.log10(y[valid]), 1)
            fit_x = np.geomspace(float(np.min(x[valid])), float(np.max(x[valid])), 160)
            fit_y = 10.0 ** np.polyval(coeff, np.log10(fit_x))
            ax.plot(fit_x, fit_y, color=COLORS[scheme], linewidth=1.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Wall time (s)")
    ax.set_ylabel(r"Density $L^1$ error")
    ax.set_title(
        "markers are grid runs; solid lines are log-log least-square fits\n"
        "AMReX advection-blob efficiency: Density error vs wall time"
    )
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    out = REPORT_FIGS / "amrex_advection_blob_density_error_vs_walltime_with_grid_labels_polished.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)
    save_to_report_locations(out, report_copy=False)
    return out


def plot_density_grid(rows: list[dict[str, str]]) -> Path:
    fig, ax = plt.subplots(figsize=(9.6, 6.11), constrained_layout=True)
    for scheme in SCHEME_ORDER:
        data = sorted(
            [row for row in rows if row["plot_scheme"] == scheme],
            key=lambda row: int(float(row["n"])),
        )
        x = np.array([row_float(row, "n") for row in data])
        y = np.array([row_float(row, "density_l1_error") for row in data])
        valid = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
        ax.scatter(x[valid], y[valid], color=COLORS[scheme], marker=MARKERS[scheme], s=48, label=scheme)
        n = np.array([int(float(row["n"])) for row in data])
        for xx, yy, nn in zip(x[valid], y[valid], n[valid]):
            ax.annotate(str(nn), (xx, yy), xytext=(4, 3), textcoords="offset points", fontsize=7)
        if np.count_nonzero(valid) >= 3:
            coeff = np.polyfit(np.log10(x[valid]), np.log10(y[valid]), 1)
            fit_x = np.geomspace(float(np.min(x[valid])), float(np.max(x[valid])), 160)
            fit_y = 10.0 ** np.polyval(coeff, np.log10(fit_x))
            ax.plot(fit_x, fit_y, color=COLORS[scheme], linewidth=1.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("grid cells per direction")
    ax.set_ylabel("density L1 error at t=0.5")
    ax.set_title("AMReX advection blob: explicit schemes vs current IMEX")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    out = PACKAGE / "figures/amrex_advection_blob_density_error_vs_grid_primary.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)
    save_to_report_locations(out, report_copy=True)
    return out


def copy_source(rows: list[dict[str, str]]) -> None:
    PACKAGE_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    with PACKAGE_SOURCE.open("w", newline="") as handle:
        fieldnames = [key for key in rows[0].keys() if key != "plot_scheme"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def main() -> None:
    rows = read_rows()
    copy_source(rows)
    print(plot_density_walltime(rows))
    print(plot_density_grid(rows))


if __name__ == "__main__":
    main()
