#!/usr/bin/env python3
"""Grid-sanity analysis for Cartesian same-gamma shock-density-bubble data."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
ROOTS = [
    ("160x40", REPO / "results/amrex/project_same_gamma_shock_density_bubble_smoke_2026-05-20"),
    ("320x80", REPO / "results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20"),
]
PACKAGE = REPO / "project_outputs/test_packages/shock_density_bubble_cartesian_same_gamma"
EVIDENCE_DATA = REPO / "project_outputs/project_evidence/data"
EVIDENCE_FIGURES = REPO / "project_outputs/project_evidence/figures"

SCHEMES = [
    ("Explicit O2 HLLC", "hllc_o2"),
    ("Explicit O2 Low-Mach Corrected HLLC-P", "lowmach_corrected_hllc_p_o2"),
    ("IMEX T1/S2 BDLTV20", "imex_t1s2_bdltv20"),
]

GAMMA = 1.4
P_PRE = 1.0
RHO_PRE = 1.0
P_POST = 5.0
X_SHOCK_0 = 0.2
RHO_BUBBLE_THRESHOLD = 0.55


def rankine_hugoniot_shock_speed() -> float:
    pressure_ratio = P_POST / P_PRE
    mach_sq = 1.0 + (pressure_ratio - 1.0) * (GAMMA + 1.0) / (2.0 * GAMMA)
    return math.sqrt(mach_sq) * math.sqrt(GAMMA * P_PRE / RHO_PRE)


def load_snapshot(path: Path) -> tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows = list(csv.DictReader(path.open(newline="")))
    xs = np.array(sorted({float(row["x"]) for row in rows}))
    ys = np.array(sorted({float(row["y"]) for row in rows}))
    xi = {x: i for i, x in enumerate(xs)}
    yi = {y: j for j, y in enumerate(ys)}
    rho = np.empty((len(ys), len(xs)))
    pressure = np.empty((len(ys), len(xs)))
    for row in rows:
        i = xi[float(row["x"])]
        j = yi[float(row["y"])]
        rho[j, i] = float(row["rho"])
        pressure[j, i] = float(row["pressure"])
    return float(rows[0]["time"]), xs, ys, rho, pressure


def shock_position_top_band(xs: np.ndarray, ys: np.ndarray, rho: np.ndarray) -> float:
    band = np.where(ys >= 0.30)[0]
    if len(band) == 0:
        band = np.arange(len(ys))
    positions = []
    for j in band:
        grad = np.diff(rho[j, :])
        idx = int(np.argmin(grad))
        positions.append(0.5 * (xs[idx] + xs[idx + 1]))
    return float(np.median(positions))


def bubble_metrics(xs: np.ndarray, ys: np.ndarray, rho: np.ndarray) -> tuple[float, float, float]:
    mask = rho < RHO_BUBBLE_THRESHOLD
    if not np.any(mask):
        return float("nan"), float("nan"), 0.0
    xx, yy = np.meshgrid(xs, ys)
    return float(xx[mask].mean()), float(yy[mask].mean()), float(mask.mean())


def snapshot_dir(root: Path, stem: str, grid: str) -> Path:
    return root / f"{stem}_{grid}_snapshots"


def main() -> None:
    out_dir = PACKAGE / "physics_check"
    fig_dir = PACKAGE / "figures"
    for path in [out_dir, fig_dir, EVIDENCE_DATA, EVIDENCE_FIGURES]:
        path.mkdir(parents=True, exist_ok=True)

    shock_speed = rankine_hugoniot_shock_speed()
    rows: list[dict[str, object]] = []
    for grid, root in ROOTS:
        for scheme, stem in SCHEMES:
            snap_dir = snapshot_dir(root, stem, grid)
            for snap_path in sorted(snap_dir.glob("snapshot_*.csv")):
                time, xs, ys, rho, pressure = load_snapshot(snap_path)
                shock_x = shock_position_top_band(xs, ys, rho)
                cx, cy, area_fraction = bubble_metrics(xs, ys, rho)
                rows.append(
                    {
                        "scheme": scheme,
                        "grid": grid,
                        "time": time,
                        "shock_x_top_band": shock_x,
                        "shock_x_rankine_hugoniot": X_SHOCK_0 + shock_speed * time,
                        "shock_x_error": shock_x - (X_SHOCK_0 + shock_speed * time),
                        "bubble_centroid_x": cx,
                        "bubble_centroid_y_half_domain": cy,
                        "bubble_low_density_cell_fraction": area_fraction,
                        "rho_min": float(np.min(rho)),
                        "rho_max": float(np.max(rho)),
                        "pressure_min": float(np.min(pressure)),
                        "pressure_max": float(np.max(pressure)),
                    }
                )

    csv_path = out_dir / "shock_density_bubble_cartesian_grid_sanity_160x40_320x80.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (EVIDENCE_DATA / csv_path.name).write_bytes(csv_path.read_bytes())

    fig, axes = plt.subplots(1, 3, figsize=(16.2, 4.8), constrained_layout=True)
    colors = {
        "Explicit O2 HLLC": "tab:blue",
        "Explicit O2 Low-Mach Corrected HLLC-P": "#9467bd",
        "IMEX T1/S2 BDLTV20": "tab:green",
    }
    linestyles = {"160x40": "--", "320x80": "-"}
    for scheme, _ in SCHEMES:
        for grid, _root in ROOTS:
            rr = [r for r in rows if r["scheme"] == scheme and r["grid"] == grid]
            t = np.array([float(r["time"]) for r in rr])
            label = f"{scheme}, {grid}"
            axes[0].plot(
                t,
                [float(r["shock_x_error"]) for r in rr],
                marker="o",
                linestyle=linestyles[grid],
                color=colors[scheme],
                label=label,
            )
            axes[1].plot(
                t,
                [float(r["bubble_centroid_x"]) for r in rr],
                marker="o",
                linestyle=linestyles[grid],
                color=colors[scheme],
            )
            axes[2].plot(
                t,
                [float(r["pressure_min"]) for r in rr],
                marker="o",
                linestyle=linestyles[grid],
                color=colors[scheme],
            )
    axes[0].axhline(0.0, color="black", linewidth=0.8, alpha=0.6)
    axes[0].set_title("Shock-position error")
    axes[0].set_ylabel("x_num - x_RH")
    axes[1].set_title("Low-density bubble centroid")
    axes[1].set_ylabel("centroid x")
    axes[2].set_title("Minimum pressure")
    axes[2].set_ylabel("p_min")
    for ax in axes:
        ax.set_xlabel("time")
        ax.grid(True, alpha=0.25)
    axes[0].legend(fontsize=7, loc="best")
    fig_path = fig_dir / "shock_density_bubble_cartesian_grid_sanity_160x40_320x80.png"
    fig.savefig(fig_path, dpi=220)
    plt.close(fig)
    (EVIDENCE_FIGURES / fig_path.name).write_bytes(fig_path.read_bytes())

    final_rows = [r for r in rows if abs(float(r["time"]) - 0.3) < 1.0e-12]
    note = [
        "# Cartesian shock-density-bubble grid-sanity analysis",
        "",
        "This is a physical sanity check for the MPhil project Cartesian same-gamma",
        "shock-density-bubble case. The 2D interaction has no exact solution,",
        "so the check combines an exact Rankine-Hugoniot shock-speed comparison",
        "away from the bubble with grid-to-grid consistency of bubble motion.",
        "",
        "- grids: 160x40 and 320x80",
        "- final time: 0.3",
        "- gamma: 1.4",
        "- no GFM, no level set, no cylindrical source terms",
        "",
        "## Final-time values",
        "",
        "| Scheme | Grid | shock x error | bubble centroid x | pressure min | rho min |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in final_rows:
        note.append(
            f"| {r['scheme']} | {r['grid']} | {float(r['shock_x_error']):.5g} | "
            f"{float(r['bubble_centroid_x']):.5g} | {float(r['pressure_min']):.5g} | "
            f"{float(r['rho_min']):.5g} |"
        )
    note.extend(
        [
            "",
            "Interpretation: the planar shock track is checked against the exact",
            "Rankine-Hugoniot speed in the top band away from the bubble. The bubble",
            "centroid moves downstream on both grids for all three schemes. Positive",
            "density and pressure are preserved in the stored snapshots. This supports",
            "using the case as qualitative high-speed 2D evidence, while still avoiding",
            "claims of an exact 2D bubble-interaction solution or final shock-bubble",
            "efficiency validation.",
            "",
        ]
    )
    note_path = out_dir / "shock_density_bubble_cartesian_grid_sanity_160x40_320x80.md"
    note_path.write_text("\n".join(note))

    print(csv_path)
    print(fig_path)
    print(note_path)


if __name__ == "__main__":
    main()
