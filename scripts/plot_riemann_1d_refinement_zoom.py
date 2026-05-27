#!/usr/bin/env python3
"""Refresh the 1D Riemann shock/contact zoom figure."""

from __future__ import annotations

import csv
import shutil
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt

from riemann_exact import RiemannState, sample, solve_star_region


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "project_outputs/test_packages/riemann_1d"
EVIDENCE = ROOT / "project_outputs/project_evidence"
FIELDS = PACKAGE / "source_data/fields"
FIGURE_NAME = "riemann_1d_primary_overlays_wave_zoom.png"


@dataclass(frozen=True)
class Series:
    case: str
    row_label: str
    toro_test: int
    final_time: float
    scheme: str
    path: Path
    color: str


STATES = {
    1: (RiemannState(1.0, 0.0, 1.0), RiemannState(0.125, 0.0, 0.1)),
    3: (RiemannState(1.0, 0.0, 1000.0), RiemannState(1.0, 0.0, 0.01)),
    4: (RiemannState(1.0, 0.0, 0.01), RiemannState(1.0, 0.0, 100.0)),
    93: (RiemannState(0.445, 1.698, 3.528), RiemannState(0.5, 0.0, 0.571)),
}

COLORS = {
    "Explicit O2 HLLC": "#1f77b4",
    "Explicit O2 Low-Mach HLLC-P": "#9467bd",
    "IMEX T1/S2 BDLTV20": "#2ca02c",
}

MIN_PANEL_WIDTH = 0.035
MAX_PANEL_WIDTH = 0.120
ROI_THRESHOLD_FRACTION = 0.20
ROI_EXPAND_POINTS = 8


def main() -> None:
    series = [
        Series("sod_toro1", "Sod / Toro 1", 1, 0.2, "Explicit O2 HLLC",
               FIELDS / "sod_toro1_explicit_o2_hllc_n400.csv", COLORS["Explicit O2 HLLC"]),
        Series("sod_toro1", "Sod / Toro 1", 1, 0.2, "Explicit O2 Low-Mach HLLC-P",
               FIELDS / "sod_toro1_explicit_o2_lowmach_hllcp_n400.csv", COLORS["Explicit O2 Low-Mach HLLC-P"]),
        Series("sod_toro1", "Sod / Toro 1", 1, 0.2, "IMEX T1/S2 BDLTV20",
               FIELDS / "sod_toro1_imex_t1s2_bdltv20_n400.csv", COLORS["IMEX T1/S2 BDLTV20"]),
        Series("lax_rp2", "Lax / BDLTV20 RP2", 93, 0.14, "Explicit O2 HLLC",
               FIELDS / "lax_rp2_explicit_o2_hllc_n400.csv", COLORS["Explicit O2 HLLC"]),
        Series("lax_rp2", "Lax / BDLTV20 RP2", 93, 0.14, "Explicit O2 Low-Mach HLLC-P",
               FIELDS / "lax_rp2_explicit_o2_lowmach_hllcp_n400.csv", COLORS["Explicit O2 Low-Mach HLLC-P"]),
        Series("lax_rp2", "Lax / BDLTV20 RP2", 93, 0.14, "IMEX T1/S2 BDLTV20",
               FIELDS / "lax_rp2_imex_t1s2_bdltv20_n400.csv", COLORS["IMEX T1/S2 BDLTV20"]),
        Series("toro3_strong_shock", "Toro 3", 3, 0.012, "Explicit O2 HLLC",
               FIELDS / "toro3_strong_shock_explicit_o2_hllc_n400.csv", COLORS["Explicit O2 HLLC"]),
        Series("toro3_strong_shock", "Toro 3", 3, 0.012, "Explicit O2 Low-Mach HLLC-P",
               FIELDS / "toro3_strong_shock_explicit_o2_lowmach_hllcp_n400.csv", COLORS["Explicit O2 Low-Mach HLLC-P"]),
        Series("toro3_strong_shock", "Toro 3", 3, 0.012, "IMEX T1/S2 BDLTV20",
               FIELDS / "toro3_strong_shock_imex_t1s2_bdltv20_n400.csv", COLORS["IMEX T1/S2 BDLTV20"]),
        Series("toro4_strong_shock", "Toro 4", 4, 0.035, "Explicit O2 HLLC",
               FIELDS / "toro4_strong_shock_explicit_o2_hllc_n400.csv", COLORS["Explicit O2 HLLC"]),
        Series("toro4_strong_shock", "Toro 4", 4, 0.035, "Explicit O2 Low-Mach HLLC-P",
               FIELDS / "toro4_strong_shock_explicit_o2_lowmach_hllcp_n400.csv", COLORS["Explicit O2 Low-Mach HLLC-P"]),
        Series("toro4_strong_shock", "Toro 4", 4, 0.035, "IMEX T1/S2 BDLTV20",
               FIELDS / "toro4_strong_shock_imex_t1s2_bdltv20_n400.csv", COLORS["IMEX T1/S2 BDLTV20"]),
    ]
    for item in series:
        if not item.path.exists():
            raise SystemExit(f"Missing source profile: {item.path}")

    output = PACKAGE / "figures" / FIGURE_NAME
    plot_zoom(series, output)
    (EVIDENCE / "figures").mkdir(parents=True, exist_ok=True)
    shutil.copy2(output, EVIDENCE / "figures" / FIGURE_NAME)
    print(output)
    print(EVIDENCE / "figures" / FIGURE_NAME)


def read_profile(path: Path) -> dict[str, list[float]]:
    with path.open(newline="") as f:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(f)]
    y_values = sorted({row["y"] for row in rows})
    y_pick = min(y_values, key=lambda value: abs(value - 0.5 * (y_values[0] + y_values[-1])))
    line = sorted((row for row in rows if row["y"] == y_pick), key=lambda row: row["x"])
    return {
        "x": [row["x"] for row in line],
        "rho": [row["rho"] for row in line],
        "u": [row["u"] for row in line],
        "p": [row["p"] for row in line],
    }


def exact_profile(item: Series, xs: list[float]) -> dict[str, list[float]]:
    left, right = STATES[item.toro_test]
    star = solve_star_region(left, right, 1.4)
    values = [sample(x, item.final_time, left, right, 1.4, 0.5, star) for x in xs]
    return {
        "rho": [state.rho for state in values],
        "u": [state.u for state in values],
        "p": [state.p for state in values],
    }


def plot_zoom(series: list[Series], output: Path) -> None:
    cases = list(dict.fromkeys(item.case for item in series))
    quantities = [("rho", "Density"), ("u", "Velocity"), ("p", "Pressure")]
    fig, axes = plt.subplots(
        len(cases),
        len(quantities),
        figsize=(17.5, 13.2),
        squeeze=False,
        constrained_layout=True,
    )
    for row_index, case in enumerate(cases):
        case_series = [item for item in series if item.case == case]
        reference = read_profile(case_series[0].path)
        exact = exact_profile(case_series[0], reference["x"])
        for col_index, (key, title) in enumerate(quantities):
            ax = axes[row_index][col_index]
            profiles = [(item, read_profile(item.path)) for item in case_series]
            roi = find_roi(reference["x"], [exact[key]] + [profile[key] for _, profile in profiles])
            exact_x, exact_y = slice_xy(reference["x"], exact[key], roi)
            ax.plot(exact_x, exact_y, color="black", linestyle="--", linewidth=2.0, label="Exact")
            for item in case_series:
                profile = dict(profiles)[item]
                xs, ys = slice_xy(profile["x"], profile[key], roi)
                ax.plot(
                    xs,
                    ys,
                    color=item.color,
                    linewidth=1.7,
                    marker="o",
                    markersize=2.0,
                    label=item.scheme,
                )
            ax.set_xlim(roi[0], roi[1])
            set_local_ylim(ax, [exact_y] + [slice_xy(profile["x"], profile[key], roi)[1] for _, profile in profiles])
            ax.grid(True, alpha=0.25)
            ax.set_title(title if row_index == 0 else "")
            if col_index == 0:
                ax.set_ylabel(case_series[0].row_label)
                ax.text(
                    0.03,
                    0.93,
                    f"t = {case_series[0].final_time:g}",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=10,
                    bbox={"facecolor": "white", "edgecolor": "0.8", "alpha": 0.72},
                )
            if row_index == len(cases) - 1:
                ax.set_xlabel("x")
            if row_index == 0 and col_index == 0:
                ax.legend(loc="best", fontsize=9)
    fig.suptitle("AMReX 1D Riemann tests: data-selected wave-feature zooms", fontsize=18)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def find_roi(xs: list[float], curves: list[list[float]]) -> tuple[float, float]:
    if len(xs) < 4:
        return min(xs), max(xs)
    scores = [0.0 for _ in range(len(xs) - 1)]
    for values in curves:
        spread = max(values) - min(values)
        scale = spread if spread > 1.0e-14 else 1.0
        for i in range(len(xs) - 1):
            scores[i] += abs(values[i + 1] - values[i]) / scale
    smooth = []
    for i in range(len(scores)):
        lo = max(0, i - 2)
        hi = min(len(scores), i + 3)
        smooth.append(sum(scores[lo:hi]) / (hi - lo))
    peak = max(range(len(smooth)), key=lambda i: smooth[i])
    cutoff = ROI_THRESHOLD_FRACTION * smooth[peak]
    left = peak
    while left > 0 and smooth[left - 1] >= cutoff:
        left -= 1
    right = peak
    while right + 1 < len(smooth) and smooth[right + 1] >= cutoff:
        right += 1
    left = max(0, left - ROI_EXPAND_POINTS)
    right = min(len(xs) - 2, right + ROI_EXPAND_POINTS)
    center = 0.5 * (xs[peak] + xs[peak + 1])
    xmin = xs[left]
    xmax = xs[right + 1]
    if xmax - xmin < MIN_PANEL_WIDTH:
        half = 0.5 * MIN_PANEL_WIDTH
        xmin, xmax = center - half, center + half
    if xmax - xmin > MAX_PANEL_WIDTH:
        half = 0.5 * MAX_PANEL_WIDTH
        xmin, xmax = center - half, center + half
    xmin = max(min(xs), xmin)
    xmax = min(max(xs), xmax)
    return xmin, xmax


def slice_xy(xs: list[float], ys: list[float], roi: tuple[float, float]) -> tuple[list[float], list[float]]:
    selected = [(x, y) for x, y in zip(xs, ys) if roi[0] <= x <= roi[1]]
    if len(selected) < 2:
        center = 0.5 * (roi[0] + roi[1])
        index = min(range(len(xs)), key=lambda i: abs(xs[i] - center))
        lo = max(0, index - 1)
        hi = min(len(xs), index + 2)
        selected = list(zip(xs[lo:hi], ys[lo:hi]))
    return [x for x, _ in selected], [y for _, y in selected]


def set_local_ylim(ax, curve_values: list[list[float]]) -> None:
    values = [value for curve in curve_values for value in curve]
    ymin = min(values)
    ymax = max(values)
    spread = ymax - ymin
    pad = 0.08 * spread if spread > 1.0e-14 else 0.1
    ax.set_ylim(ymin - pad, ymax + pad)


if __name__ == "__main__":
    main()
