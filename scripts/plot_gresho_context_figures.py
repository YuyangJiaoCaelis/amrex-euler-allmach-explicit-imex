#!/usr/bin/env python3
"""Create Gresho field-snapshot and target-error cost figures."""

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
EVIDENCE_FIGS = ROOT / "project_outputs/project_evidence/figures"
SOURCE = PACKAGE / "source_data/amrex_gresho_efficiency_combined_data.csv"

SCHEME_ORDER = (
    "Explicit O2 HLLC",
    "Explicit O2 Low-Mach Corrected HLLC-P",
    "IMEX T1/S2 BDLTV20",
)

SCHEME_KEYS = {
    "explicit_o2_hllc": "Explicit O2 HLLC",
    "explicit_o2_xie_am_hllc_p": "Explicit O2 Low-Mach Corrected HLLC-P",
    "imex_t1s2_bdltv20": "IMEX T1/S2 BDLTV20",
    "imex_bdltv20_paper_t1s2": "IMEX T1/S2 BDLTV20",
}

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


def read_rows() -> list[dict[str, str]]:
    with SOURCE.open(newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]
    for row in rows:
        label = SCHEME_KEYS.get(row.get("scheme_key", ""))
        if not label:
            label = SCHEME_ALIASES.get(row.get("scheme_label", ""), row.get("scheme_label", ""))
        row["plot_scheme"] = label
    return [row for row in rows if row["plot_scheme"] in SCHEME_ORDER]


def field_value(row: dict[str, str], *keys: str) -> float:
    for key in keys:
        text = row.get(key, "")
        if text:
            return float(text)
    return math.nan


def read_field(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]
    if not rows:
        raise ValueError(f"empty field CSV: {path}")
    xs = np.array(sorted({float(row["x"]) for row in rows}))
    ys = np.array(sorted({float(row["y"]) for row in rows}))
    xi = {x: i for i, x in enumerate(xs)}
    yi = {y: j for j, y in enumerate(ys)}
    fields = {
        "x": xs,
        "y": ys,
        "u": np.empty((len(ys), len(xs))),
        "v": np.empty((len(ys), len(xs))),
        "p": np.empty((len(ys), len(xs))),
        "exact_u": np.empty((len(ys), len(xs))),
        "exact_v": np.empty((len(ys), len(xs))),
        "exact_p": np.empty((len(ys), len(xs))),
    }
    for row in rows:
        i = xi[float(row["x"])]
        j = yi[float(row["y"])]
        fields["u"][j, i] = field_value(row, "u")
        fields["v"][j, i] = field_value(row, "v")
        fields["p"][j, i] = field_value(row, "p", "pressure")
        fields["exact_u"][j, i] = field_value(row, "exact_u")
        fields["exact_v"][j, i] = field_value(row, "exact_v")
        fields["exact_p"][j, i] = field_value(row, "exact_p", "exact_pressure")
    return fields


def choose_snapshot_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected = []
    for scheme in SCHEME_ORDER:
        matches = [
            row
            for row in rows
            if row["plot_scheme"] == scheme
            and abs(float(row["mach"]) - 0.01) < 1.0e-14
            and int(float(row["n"])) == 48
            and row.get("final_csv")
        ]
        if not matches:
            raise RuntimeError(f"missing Gresho Mach 0.01 n=48 field row for {scheme}")
        selected.append(matches[0])
    return selected


def save_to_report_locations(path: Path) -> None:
    REPORT_FIGS.mkdir(parents=True, exist_ok=True)
    EVIDENCE_FIGS.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, REPORT_FIGS / path.name)
    shutil.copy2(path, EVIDENCE_FIGS / path.name)


def make_field_snapshot(rows: list[dict[str, str]]) -> Path:
    selected = choose_snapshot_rows(rows)
    scheme_fields = [(row["plot_scheme"], read_field(Path(row["final_csv"]))) for row in selected]
    reference = scheme_fields[0][1]
    fields = [
        (
            "Exact",
            {
                "x": reference["x"],
                "y": reference["y"],
                "u": reference["exact_u"],
                "v": reference["exact_v"],
                "p": reference["exact_p"],
            },
        ),
        *scheme_fields,
    ]
    speeds = [np.hypot(item["u"], item["v"]) for _, item in fields]
    pressure_background = 1.0 / (1.4 * 0.01 * 0.01)
    pressure_pert = [item["p"] - pressure_background for _, item in fields]
    speed_min = min(float(np.min(s)) for s in speeds)
    speed_max = max(float(np.max(s)) for s in speeds)
    pressure_min = min(float(np.nanmin(p)) for p in pressure_pert)
    pressure_max = max(float(np.nanmax(p)) for p in pressure_pert)

    fig, axes = plt.subplots(2, len(fields), figsize=(16.5, 7.6), constrained_layout=True)
    im_speed = None
    im_pressure = None
    for col, (label, item) in enumerate(fields):
        extent = (item["x"].min(), item["x"].max(), item["y"].min(), item["y"].max())
        im_speed = axes[0, col].imshow(
            speeds[col],
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap="viridis",
            vmin=speed_min,
            vmax=speed_max,
        )
        axes[0, col].set_title(label, fontsize=10)
        axes[0, col].set_xlabel("x")
        if col == 0:
            axes[0, col].set_ylabel(r"Speed $|\mathbf{u}|$")

        im_pressure = axes[1, col].imshow(
            pressure_pert[col],
            origin="lower",
            extent=extent,
            aspect="equal",
            cmap="coolwarm",
            vmin=pressure_min,
            vmax=pressure_max,
        )
        axes[1, col].set_xlabel("x")
        if col == 0:
            axes[1, col].set_ylabel(r"Pressure perturbation $p-p_0$")

    fig.colorbar(im_speed, ax=axes[0, :].tolist(), fraction=0.025, pad=0.02)
    fig.colorbar(im_pressure, ax=axes[1, :].tolist(), fraction=0.025, pad=0.02)
    fig.suptitle(
        r"AMReX Gresho field snapshot, Mach 0.01, $n=48$, $t=0.4\pi$"
        "\nExact Dirichlet boundary; local CPU evidence only",
        fontsize=16,
    )
    out = PACKAGE / "figures/amrex_gresho_field_snapshot_mach0p01_n48_t0p4pi.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)
    save_to_report_locations(out)
    return out


def estimate_target_time(rows: list[dict[str, str]], metric: str, target: float) -> dict[tuple[str, float], dict[str, str]]:
    estimates: dict[tuple[str, float], dict[str, str]] = {}
    for scheme in SCHEME_ORDER:
        for mach in sorted({float(row["mach"]) for row in rows}):
            data = [
                row
                for row in rows
                if row["plot_scheme"] == scheme
                and abs(float(row["mach"]) - mach) < 1.0e-14
                and row.get(metric, "")
            ]
            x = np.array([field_value(row, "driver_wall_time_sec", "wall_time_sec") for row in data])
            y = np.array([field_value(row, metric) for row in data])
            valid = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
            status = "no_decreasing_fit"
            value = math.nan
            if np.count_nonzero(valid) >= 3:
                coeff = np.polyfit(np.log10(x[valid]), np.log10(y[valid]), 1)
                if coeff[0] < 0.0:
                    value = 10.0 ** ((math.log10(target) - coeff[1]) / coeff[0])
                    status = "interpolated" if x[valid].min() <= value <= x[valid].max() else "extrapolated"
            estimates[(scheme, mach)] = {
                "scheme": scheme,
                "mach": f"{mach:g}",
                "metric": metric,
                "target": f"{target:g}",
                "target_wall_time_sec": "" if not math.isfinite(value) else f"{value:.12g}",
                "estimate_status": status,
            }
    return estimates


def format_seconds(seconds: float) -> str:
    if not math.isfinite(seconds):
        return "no fit"
    if seconds < 60.0:
        return f"{seconds:.2g}s"
    minutes = seconds / 60.0
    if minutes < 60.0:
        return f"{minutes:.2g}m"
    hours = minutes / 60.0
    return f"{hours:.2g}h"


def make_target_cost(rows: list[dict[str, str]]) -> Path:
    metrics = (
        ("pressure_perturbation_l1_relative_error", "Pressure perturbation relative L1"),
        ("velocity_l1_error", "Velocity L1"),
    )
    target = 5.0e-3
    machs = sorted({float(row["mach"]) for row in rows})
    estimates = {metric: estimate_target_time(rows, metric, target) for metric, _ in metrics}
    fig, axes = plt.subplots(len(metrics), len(machs), figsize=(18.0, 8.8), sharey=False)
    for row_idx, (metric, metric_title) in enumerate(metrics):
        finite_values = [
            float(item["target_wall_time_sec"])
            for item in estimates[metric].values()
            if item["target_wall_time_sec"]
        ]
        ymin = max(min(finite_values) * 0.2, 1.0e-1) if finite_values else 1.0e-1
        ymax = max(finite_values) * 8.0 if finite_values else 1.0
        for col_idx, mach in enumerate(machs):
            ax = axes[row_idx, col_idx]
            xs = np.arange(len(SCHEME_ORDER))
            heights = []
            labels = []
            colors = []
            edgecolors = []
            hatches = []
            for scheme in SCHEME_ORDER:
                item = estimates[metric][(scheme, mach)]
                value_text = item["target_wall_time_sec"]
                if value_text:
                    value = float(value_text)
                    heights.append(value)
                    labels.append(format_seconds(value))
                    colors.append(COLORS[scheme])
                    edgecolors.append("black" if item["estimate_status"] == "extrapolated" else COLORS[scheme])
                    hatches.append("")
                else:
                    heights.append(ymin)
                    labels.append("no fit")
                    colors.append("white")
                    edgecolors.append(COLORS[scheme])
                    hatches.append("xx")
            bars = ax.bar(xs, heights, color=colors, edgecolor=edgecolors, linewidth=1.5)
            for bar, hatch, label, color in zip(bars, hatches, labels, colors):
                bar.set_hatch(hatch)
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    bar.get_height() * 1.22,
                    label,
                    ha="center",
                    va="bottom",
                    rotation=22,
                    fontsize=8.5,
                    color="0.15" if color == "white" else color,
                    fontweight="bold",
                )
            ax.set_yscale("log")
            ax.set_ylim(ymin, ymax)
            ax.set_xticks(xs)
            ax.set_xticklabels(["HLLC", "Low-Mach\nHLLC-P", "IMEX\nBDLTV20"])
            ax.set_title(f"Mach {mach:g}")
            ax.grid(True, axis="y", which="both", alpha=0.25)
            if col_idx == 0:
                ax.set_ylabel(f"Time to {metric_title}\n<= {target:g} (s)")
    fig.suptitle(
        "AMReX Gresho target-error cost comparison\n"
        "bars show local CPU time inferred from grid sweeps; black outline means extrapolated fit",
        fontsize=16,
    )
    out = PACKAGE / "figures/amrex_gresho_target_error_cost_bar_primary.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)
    save_to_report_locations(out)

    source_out = PACKAGE / "source_data/amrex_gresho_target_time_colour_bar_current.csv"
    with source_out.open("w", newline="") as handle:
        fieldnames = ["scheme", "mach", "metric", "target", "target_wall_time_sec", "estimate_status"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for metric, _ in metrics:
            writer.writerows(estimates[metric].values())
    return out


def main() -> None:
    rows = read_rows()
    print(make_field_snapshot(rows))
    print(make_target_cost(rows))


if __name__ == "__main__":
    main()
