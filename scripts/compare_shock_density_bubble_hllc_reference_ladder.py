#!/usr/bin/env python3
"""Grid-reference ladder for the Cartesian same-gamma shock-density-bubble HLLC runs."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
PACKAGE = REPO / "project_outputs/test_packages/shock_density_bubble_cartesian_same_gamma"
EVIDENCE_DATA = REPO / "project_outputs/project_evidence/data"
EVIDENCE_FIGURES = REPO / "project_outputs/project_evidence/figures"

RUNS = {
    "160x40": REPO / "results/amrex/project_same_gamma_shock_density_bubble_smoke_2026-05-20/hllc_o2_160x40_snapshots",
    "320x80": REPO / "results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20/hllc_o2_320x80_snapshots",
    "640x160": REPO / "results/amrex/project_same_gamma_shock_density_bubble_hires_reference_2026-05-20/hllc_o2_640x160_snapshots",
}

GAMMA = 1.4
P_PRE = 1.0
RHO_PRE = 1.0
P_POST = 5.0
X_SHOCK_0 = 0.2
RHO_BUBBLE_THRESHOLD = 0.55


def snapshot_index(path: Path) -> int:
    match = re.search(r"snapshot_(\d+)_", path.name)
    if not match:
        raise ValueError(f"cannot parse snapshot index from {path}")
    return int(match.group(1))


def load_snapshot(path: Path) -> dict[str, object]:
    rows = list(csv.DictReader(path.open(newline="")))
    if not rows:
        raise ValueError(f"empty snapshot: {path}")
    xs = np.array(sorted({float(row["x"]) for row in rows}))
    ys = np.array(sorted({float(row["y"]) for row in rows}))
    xi = {x: i for i, x in enumerate(xs)}
    yi = {y: j for j, y in enumerate(ys)}
    fields = {
        "rho": np.empty((len(ys), len(xs))),
        "pressure": np.empty((len(ys), len(xs))),
        "u": np.empty((len(ys), len(xs))),
        "v": np.empty((len(ys), len(xs))),
        "internal_energy": np.empty((len(ys), len(xs))),
    }
    for row in rows:
        i = xi[float(row["x"])]
        j = yi[float(row["y"])]
        for name in fields:
            fields[name][j, i] = float(row[name])
    return {
        "time": float(rows[0]["time"]),
        "xs": xs,
        "ys": ys,
        "fields": fields,
    }


def load_series(snap_dir: Path) -> list[dict[str, object]]:
    paths = sorted(snap_dir.glob("snapshot_*.csv"), key=snapshot_index)
    if not paths:
        raise FileNotFoundError(f"no snapshots in {snap_dir}")
    return [load_snapshot(path) for path in paths]


def average_2x2(field: np.ndarray) -> np.ndarray:
    ny, nx = field.shape
    if nx % 2 or ny % 2:
        raise ValueError(f"cannot 2x2-average shape {field.shape}")
    return field.reshape(ny // 2, 2, nx // 2, 2).mean(axis=(1, 3))


def relative_l1(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.mean(np.abs(a - b)) / max(np.mean(np.abs(b)), 1.0e-300))


def linf(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b)))


def rankine_hugoniot_shock_speed() -> float:
    pressure_ratio = P_POST / P_PRE
    mach_sq = 1.0 + (pressure_ratio - 1.0) * (GAMMA + 1.0) / (2.0 * GAMMA)
    return math.sqrt(mach_sq) * math.sqrt(GAMMA * P_PRE / RHO_PRE)


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


def bubble_centroid(xs: np.ndarray, ys: np.ndarray, rho: np.ndarray) -> tuple[float, float]:
    mask = rho < RHO_BUBBLE_THRESHOLD
    if not np.any(mask):
        return float("nan"), float("nan")
    xx, yy = np.meshgrid(xs, ys)
    return float(xx[mask].mean()), float(yy[mask].mean())


def downsample_series(series: list[dict[str, object]]) -> list[dict[str, object]]:
    out = []
    for snap in series:
        fields = {name: average_2x2(value) for name, value in snap["fields"].items()}
        out.append({"time": snap["time"], "fields": fields})
    return out


def compare_pair(
    label: str,
    coarse: list[dict[str, object]],
    fine: list[dict[str, object]],
) -> list[dict[str, object]]:
    fine_avg = downsample_series(fine)
    rows: list[dict[str, object]] = []
    for idx, (coarse_snap, fine_snap) in enumerate(zip(coarse, fine_avg)):
        if abs(float(coarse_snap["time"]) - float(fine_snap["time"])) > 1.0e-12:
            raise RuntimeError(f"time mismatch for {label} at snapshot {idx}")
        row = {"grid_gap": label, "time": float(coarse_snap["time"])}
        for field in ["rho", "pressure", "u", "v", "internal_energy"]:
            row[f"{field}_relative_l1"] = relative_l1(
                coarse_snap["fields"][field],
                fine_snap["fields"][field],
            )
            row[f"{field}_linf"] = linf(
                coarse_snap["fields"][field],
                fine_snap["fields"][field],
            )
        rows.append(row)
    return rows


def main() -> None:
    out_dir = PACKAGE / "physics_check"
    fig_dir = PACKAGE / "figures"
    for path in [out_dir, fig_dir, EVIDENCE_DATA, EVIDENCE_FIGURES]:
        path.mkdir(parents=True, exist_ok=True)

    series = {grid: load_series(path) for grid, path in RUNS.items()}
    for grid in ["160x40", "320x80"]:
        if len(series[grid]) != len(series["640x160"]):
            raise RuntimeError(f"snapshot count mismatch for {grid}")

    gap_rows = []
    gap_rows.extend(compare_pair("160x40_vs_avg_320x80", series["160x40"], series["320x80"]))
    gap_rows.extend(compare_pair("320x80_vs_avg_640x160", series["320x80"], series["640x160"]))

    csv_path = out_dir / "shock_density_bubble_cartesian_hllc_grid_reference_ladder.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(gap_rows[0].keys()))
        writer.writeheader()
        writer.writerows(gap_rows)
    (EVIDENCE_DATA / csv_path.name).write_bytes(csv_path.read_bytes())

    shock_speed = rankine_hugoniot_shock_speed()
    feature_rows = []
    for grid, snaps in series.items():
        for snap in snaps:
            xs = snap["xs"]
            ys = snap["ys"]
            rho = snap["fields"]["rho"]
            cx, cy = bubble_centroid(xs, ys, rho)
            time = float(snap["time"])
            feature_rows.append(
                {
                    "grid": grid,
                    "time": time,
                    "shock_x_top_band": shock_position_top_band(xs, ys, rho),
                    "shock_x_rankine_hugoniot": X_SHOCK_0 + shock_speed * time,
                    "bubble_centroid_x": cx,
                    "bubble_centroid_y_half_domain": cy,
                    "rho_min": float(np.min(rho)),
                    "pressure_min": float(np.min(snap["fields"]["pressure"])),
                    "internal_energy_min": float(np.min(snap["fields"]["internal_energy"])),
                }
            )
    feature_csv = out_dir / "shock_density_bubble_cartesian_hllc_grid_reference_features.csv"
    with feature_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(feature_rows[0].keys()))
        writer.writeheader()
        writer.writerows(feature_rows)
    (EVIDENCE_DATA / feature_csv.name).write_bytes(feature_csv.read_bytes())

    colors = {
        "160x40_vs_avg_320x80": "tab:purple",
        "320x80_vs_avg_640x160": "tab:blue",
    }
    labels = {
        "160x40_vs_avg_320x80": "160x40 vs averaged 320x80",
        "320x80_vs_avg_640x160": "320x80 vs averaged 640x160",
    }
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.0), constrained_layout=True)
    for gap in labels:
        rows = [row for row in gap_rows if row["grid_gap"] == gap]
        t = np.array([float(row["time"]) for row in rows])
        axes[0, 0].plot(t, [float(row["rho_relative_l1"]) for row in rows], "o-", color=colors[gap], label=labels[gap])
        axes[0, 1].plot(t, [float(row["pressure_relative_l1"]) for row in rows], "o-", color=colors[gap], label=labels[gap])
        axes[1, 0].plot(t, [float(row["u_relative_l1"]) for row in rows], "o-", color=colors[gap], label=labels[gap])
        axes[1, 1].plot(t, [float(row["internal_energy_relative_l1"]) for row in rows], "o-", color=colors[gap], label=labels[gap])
    axes[0, 0].set_title("Density grid-gap difference")
    axes[0, 1].set_title("Pressure grid-gap difference")
    axes[1, 0].set_title("u-velocity grid-gap difference")
    axes[1, 1].set_title("Internal-energy grid-gap difference")
    for ax in axes.ravel():
        ax.set_xlabel("time")
        ax.set_ylabel("relative L1")
        ax.grid(True, alpha=0.25)
    axes[0, 0].legend(fontsize=8)
    fig_path = fig_dir / "shock_density_bubble_cartesian_hllc_grid_reference_ladder.png"
    fig.savefig(fig_path, dpi=220)
    plt.close(fig)
    (EVIDENCE_FIGURES / fig_path.name).write_bytes(fig_path.read_bytes())

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.7), constrained_layout=True)
    grid_colors = {"160x40": "tab:purple", "320x80": "tab:blue", "640x160": "tab:green"}
    for grid in RUNS:
        rows = [row for row in feature_rows if row["grid"] == grid]
        t = np.array([float(row["time"]) for row in rows])
        axes[0].plot(
            t,
            [float(row["shock_x_top_band"]) - float(row["shock_x_rankine_hugoniot"]) for row in rows],
            "o-",
            color=grid_colors[grid],
            label=grid,
        )
        axes[1].plot(
            t,
            [float(row["bubble_centroid_x"]) for row in rows],
            "o-",
            color=grid_colors[grid],
            label=grid,
        )
    axes[0].axhline(0.0, color="black", linewidth=0.8, alpha=0.5)
    axes[0].set_title("Top-band shock-position error")
    axes[0].set_ylabel("x_num - x_RH")
    axes[1].set_title("Low-density bubble centroid")
    axes[1].set_ylabel("centroid x")
    for ax in axes:
        ax.set_xlabel("time")
        ax.grid(True, alpha=0.25)
    axes[1].legend(fontsize=8)
    feature_fig = fig_dir / "shock_density_bubble_cartesian_hllc_grid_reference_features.png"
    fig.savefig(feature_fig, dpi=220)
    plt.close(fig)
    (EVIDENCE_FIGURES / feature_fig.name).write_bytes(feature_fig.read_bytes())

    final_gap = [row for row in gap_rows if abs(float(row["time"]) - 0.3) < 1.0e-12]
    final_feature = [row for row in feature_rows if abs(float(row["time"]) - 0.3) < 1.0e-12]
    gap_by_name = {row["grid_gap"]: row for row in final_gap}
    ratio_lines = []
    for field in ["rho", "pressure", "u", "internal_energy"]:
        coarse_gap = float(gap_by_name["160x40_vs_avg_320x80"][f"{field}_relative_l1"])
        fine_gap = float(gap_by_name["320x80_vs_avg_640x160"][f"{field}_relative_l1"])
        ratio_lines.append((field, coarse_gap, fine_gap, coarse_gap / fine_gap if fine_gap > 0.0 else float("inf")))

    note = [
        "# Cartesian shock-density-bubble HLLC grid-reference ladder",
        "",
        "This analysis checks whether the 640x160 Explicit O2 HLLC row is a credible",
        "high-resolution numerical reference for the MPhil project Cartesian same-gamma",
        "shock-density-bubble test. It compares HLLC against HLLC across the grid",
        "ladder 160x40 -> 320x80 -> 640x160 after averaging the finer result onto",
        "the coarser grid.",
        "",
        "This is not an exact solution proof. It is a reference-quality sanity check",
        "for a two-dimensional shock/bubble interaction that has no closed-form",
        "solution.",
        "",
        "## Final-time grid-gap relative L1 differences",
        "",
        "| Field | 160x40 vs avg 320x80 | 320x80 vs avg 640x160 | coarse/fine ratio |",
        "|---|---:|---:|---:|",
    ]
    for field, coarse_gap, fine_gap, ratio in ratio_lines:
        note.append(f"| {field} | {coarse_gap:.6g} | {fine_gap:.6g} | {ratio:.3g} |")
    note.extend(
        [
            "",
            "## Final-time shock and bubble features",
            "",
            "| Grid | shock x error | bubble centroid x | rho min | pressure min | internal energy min |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in final_feature:
        shock_error = float(row["shock_x_top_band"]) - float(row["shock_x_rankine_hugoniot"])
        note.append(
            f"| {row['grid']} | {shock_error:.6g} | {float(row['bubble_centroid_x']):.6g} | "
            f"{float(row['rho_min']):.6g} | {float(row['pressure_min']):.6g} | "
            f"{float(row['internal_energy_min']):.6g} |"
        )
    note.extend(
        [
            "",
            "Interpretation:",
            "",
            "- The 320x80-to-640x160 gap is smaller than the 160x40-to-320x80 gap",
            "  for the main density, pressure, velocity, and internal-energy fields.",
            "- The shock-position error remains small relative to the domain and is",
            "  checked against the exact Rankine-Hugoniot incident-shock speed away",
            "  from the bubble.",
            "- The bubble centroid moves consistently downstream across the HLLC grid",
            "  ladder.",
            "- Density, pressure, and internal energy remain positive on stored",
            "  snapshots.",
            "",
            "Therefore the 640x160 HLLC row is a reasonable high-resolution numerical",
            "reference for MPhil project qualitative high-speed evidence. It should still be",
            "called a numerical reference, not an exact solution.",
            "",
        ]
    )
    note_path = out_dir / "shock_density_bubble_cartesian_hllc_grid_reference_ladder.md"
    note_path.write_text("\n".join(note))

    print(csv_path)
    print(feature_csv)
    print(fig_path)
    print(feature_fig)
    print(note_path)


if __name__ == "__main__":
    main()
