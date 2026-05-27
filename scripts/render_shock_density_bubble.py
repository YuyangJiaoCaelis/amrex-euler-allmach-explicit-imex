#!/usr/bin/env python3
"""Package the MPhil project Cartesian same-gamma shock-density-bubble AMReX runs.

The input snapshots are CSV files emitted by euler_compare for the
single-material Cartesian shock-density-bubble case. The plotted view mirrors
the half-domain about y=0 to show the full symmetry pair.
"""

from __future__ import annotations

import csv
import math
import re
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


REPO = Path(__file__).resolve().parents[1]
RESULT_ROOT = REPO / "results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20"
PACKAGE_ROOT = REPO / "project_outputs/test_packages/shock_density_bubble_cartesian_same_gamma"
EVIDENCE_FIGURES = REPO / "project_outputs/project_evidence/figures"
EVIDENCE_DATA = REPO / "project_outputs/project_evidence/data"

GAMMA = 1.4
P_PRE = 1.0
RHO_PRE = 1.0
P_POST = 5.0
X_SHOCK_0 = 0.2
RHO_BUBBLE_THRESHOLD = 0.55

SCHEMES = [
    {
        "key": "hllc",
        "label": "Explicit O2 HLLC",
        "snapshots": RESULT_ROOT / "hllc_o2_320x80_snapshots",
        "summary": RESULT_ROOT / "hllc_o2_320x80_summary.csv",
    },
    {
        "key": "lowmach",
        "label": "Explicit O2 Low-Mach Corrected HLLC-P",
        "snapshots": RESULT_ROOT / "lowmach_corrected_hllc_p_o2_320x80_snapshots",
        "summary": RESULT_ROOT / "lowmach_corrected_hllc_p_o2_320x80_summary.csv",
    },
    {
        "key": "imex",
        "label": "IMEX T1/S2 BDLTV20",
        "snapshots": RESULT_ROOT / "imex_t1s2_bdltv20_320x80_snapshots",
        "summary": RESULT_ROOT / "imex_t1s2_bdltv20_320x80_summary.csv",
    },
]


def snapshot_index(path: Path) -> int:
    match = re.search(r"snapshot_(\d+)_", path.name)
    if not match:
        raise ValueError(f"Cannot parse snapshot index from {path}")
    return int(match.group(1))


def rankine_hugoniot_state() -> dict[str, float]:
    pressure_ratio = P_POST / P_PRE
    mach_sq = 1.0 + (pressure_ratio - 1.0) * (GAMMA + 1.0) / (2.0 * GAMMA)
    mach = math.sqrt(mach_sq)
    sound_pre = math.sqrt(GAMMA * P_PRE / RHO_PRE)
    shock_speed = mach * sound_pre
    density_ratio = ((GAMMA + 1.0) * mach_sq) / ((GAMMA - 1.0) * mach_sq + 2.0)
    rho_post = RHO_PRE * density_ratio
    u_post = shock_speed * (1.0 - RHO_PRE / rho_post)
    return {
        "mach": mach,
        "shock_speed": shock_speed,
        "rho_post": rho_post,
        "u_post": u_post,
        "p_post": P_POST,
    }


def load_snapshot(path: Path) -> dict[str, object]:
    rows = list(csv.DictReader(path.open(newline="")))
    if not rows:
        raise ValueError(f"Empty snapshot {path}")

    xs = np.array(sorted({float(row["x"]) for row in rows}))
    ys = np.array(sorted({float(row["y"]) for row in rows}))
    x_to_i = {x: i for i, x in enumerate(xs)}
    y_to_j = {y: j for j, y in enumerate(ys)}

    fields = {
        "rho": np.empty((len(ys), len(xs))),
        "pressure": np.empty((len(ys), len(xs))),
        "schlieren": np.empty((len(ys), len(xs))),
    }
    for row in rows:
        i = x_to_i[float(row["x"])]
        j = y_to_j[float(row["y"])]
        fields["rho"][j, i] = float(row["rho"])
        fields["pressure"][j, i] = float(row["pressure"])
        fields["schlieren"][j, i] = float(row["schlieren"])

    return {
        "time": float(rows[0]["time"]),
        "step": int(rows[0]["step"]),
        "xs": xs,
        "ys": ys,
        "fields": fields,
    }


def mirrored(field: np.ndarray) -> np.ndarray:
    return np.vstack([field[::-1, :], field])


def collect_snapshots() -> list[list[dict[str, object]]]:
    all_series = []
    for scheme in SCHEMES:
        paths = sorted(scheme["snapshots"].glob("snapshot_*.csv"), key=snapshot_index)
        if len(paths) != 6:
            raise RuntimeError(f"{scheme['label']} expected 6 snapshots, found {len(paths)}")
        all_series.append([load_snapshot(path) for path in paths])
    return all_series


def render_frame(series: list[list[dict[str, object]]], frame_index: int, out_path: Path) -> None:
    time = float(series[0][frame_index]["time"])
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4), constrained_layout=True)
    fig.suptitle(
        f"AMReX Cartesian same-gamma shock-density-bubble, t={time:.2f} "
        "(320x80 half-domain, mirrored about y=0)",
        fontsize=13,
    )
    im = None
    for ax, scheme, snapshots in zip(axes, SCHEMES, series):
        rho = mirrored(snapshots[frame_index]["fields"]["rho"])
        im = ax.imshow(
            rho,
            origin="lower",
            extent=(0.0, 2.0, -0.5, 0.5),
            cmap="viridis",
            vmin=0.1,
            vmax=2.82,
            interpolation="nearest",
            aspect="equal",
        )
        ax.set_title(scheme["label"], fontsize=10)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_xlim(0.0, 1.35)
        ax.set_ylim(-0.32, 0.32)
        ax.axhline(0.0, color="white", linewidth=0.5, alpha=0.45)
    if im is not None:
        fig.colorbar(im, ax=axes, shrink=0.82, label="Density")
    fig.savefig(out_path, dpi=170)
    plt.close(fig)


def render_contact_sheet(series: list[list[dict[str, object]]], out_path: Path) -> None:
    selected = [0, 2, 4, 5]
    fig, axes = plt.subplots(len(selected), 3, figsize=(14.5, 10.4), constrained_layout=True)
    fig.suptitle(
        "AMReX Cartesian same-gamma shock-density-bubble density snapshots "
        "(320x80, mirrored full-bubble view)",
        fontsize=13,
    )
    im = None
    for row_id, frame_index in enumerate(selected):
        time = float(series[0][frame_index]["time"])
        for col_id, (scheme, snapshots) in enumerate(zip(SCHEMES, series)):
            ax = axes[row_id, col_id]
            rho = mirrored(snapshots[frame_index]["fields"]["rho"])
            im = ax.imshow(
                rho,
                origin="lower",
                extent=(0.0, 2.0, -0.5, 0.5),
                cmap="viridis",
                vmin=0.1,
                vmax=2.82,
                interpolation="nearest",
                aspect="equal",
            )
            if row_id == 0:
                ax.set_title(scheme["label"], fontsize=10)
            if col_id == 0:
                ax.set_ylabel(f"t={time:.2f}\ny")
            else:
                ax.set_ylabel("y")
            ax.set_xlabel("x")
            ax.set_xlim(0.0, 1.35)
            ax.set_ylim(-0.32, 0.32)
            ax.axhline(0.0, color="white", linewidth=0.5, alpha=0.45)
    if im is not None:
        fig.colorbar(im, ax=axes, shrink=0.75, label="Density")
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


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


def bubble_centroid(xs: np.ndarray, ys: np.ndarray, rho: np.ndarray) -> tuple[float, float, float]:
    mask = rho < RHO_BUBBLE_THRESHOLD
    if not np.any(mask):
        return (float("nan"), float("nan"), 0.0)
    xx, yy = np.meshgrid(xs, ys)
    return (float(xx[mask].mean()), float(yy[mask].mean()), float(mask.mean()))


def render_physics_check(series: list[list[dict[str, object]]], out_csv: Path, out_fig: Path, out_note: Path) -> None:
    rh = rankine_hugoniot_state()
    rows = []
    for scheme, snapshots in zip(SCHEMES, series):
        for snap in snapshots:
            xs = snap["xs"]
            ys = snap["ys"]
            rho = snap["fields"]["rho"]
            x_num = shock_position_top_band(xs, ys, rho)
            x_ref = X_SHOCK_0 + rh["shock_speed"] * float(snap["time"])
            cx, cy, area_fraction = bubble_centroid(xs, ys, rho)
            rows.append(
                {
                    "scheme": scheme["label"],
                    "time": float(snap["time"]),
                    "shock_x_top_band": x_num,
                    "shock_x_rankine_hugoniot": x_ref,
                    "shock_x_error": x_num - x_ref,
                    "bubble_centroid_x": cx,
                    "bubble_centroid_y_half_domain": cy,
                    "bubble_low_density_cell_fraction": area_fraction,
                }
            )

    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), constrained_layout=True)
    for scheme in SCHEMES:
        scheme_rows = [row for row in rows if row["scheme"] == scheme["label"]]
        t = np.array([row["time"] for row in scheme_rows])
        shock_x = np.array([row["shock_x_top_band"] for row in scheme_rows])
        centroid_x = np.array([row["bubble_centroid_x"] for row in scheme_rows])
        axes[0].plot(t, shock_x, "o-", label=scheme["label"])
        axes[1].plot(t, centroid_x, "o-", label=scheme["label"])
    tref = np.array(sorted({row["time"] for row in rows}))
    axes[0].plot(tref, X_SHOCK_0 + rh["shock_speed"] * tref, "k--", label="Rankine-Hugoniot shock speed")
    axes[0].set_xlabel("time")
    axes[0].set_ylabel("top-band shock position")
    axes[0].set_title("Planar shock speed away from the bubble")
    axes[0].grid(True, alpha=0.25)
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("low-density bubble centroid x")
    axes[1].set_title("Downstream motion of the density bubble")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(fontsize=8, loc="best")
    fig.savefig(out_fig, dpi=220)
    plt.close(fig)

    max_err_by_scheme = {}
    for scheme in SCHEMES:
        vals = [
            abs(row["shock_x_error"])
            for row in rows
            if row["scheme"] == scheme["label"] and row["time"] > 0.0
        ]
        max_err_by_scheme[scheme["label"]] = max(vals)

    lines = [
        "# Cartesian same-gamma shock-density-bubble physics check",
        "",
        "This note compares the AMReX Cartesian same-gamma shock-density-bubble",
        "setup with the Rankine-Hugoniot post-shock state. It is a qualitative",
        "high-speed 2D robustness check, not a two-material air-helium or GFM",
        "shock-bubble validation.",
        "",
        "## Rankine-Hugoniot check",
        "",
        f"- gamma = {GAMMA}",
        f"- pressure ratio p_post/p_pre = {P_POST / P_PRE}",
        f"- incident Mach number = {rh['mach']:.8f}",
        f"- shock speed = {rh['shock_speed']:.8f}",
        f"- computed post-shock density = {rh['rho_post']:.12f}",
        f"- computed post-shock lab-frame velocity = {rh['u_post']:.12f}",
        "",
        "Maximum top-band shock-position error over t>0:",
    ]
    for scheme_label, value in max_err_by_scheme.items():
        lines.append(f"- {scheme_label}: {value:.5f}")
    lines.extend(
        [
            "",
            "Interpretation: the imposed shock state and boundary are coherent if the",
            "top-band shock position follows the Rankine-Hugoniot line. Bubble",
            "centroid motion should be monotone downstream after interaction.",
            "",
        ]
    )
    out_note.write_text("\n".join(lines))


def write_readme(out_dir: Path) -> None:
    summaries = []
    for scheme in SCHEMES:
        with scheme["summary"].open(newline="") as handle:
            row = next(csv.DictReader(handle))
        summaries.append(
            f"| {scheme['label']} | {row['status']} | {row['steps_accepted']} | "
            f"{float(row['wall_time_sec']):.3f} | {float(row['rho_min']):.4f} | "
            f"{float(row['pressure_min']):.4f} | {float(row['internal_energy_min']):.4f} |"
        )

    text = "\n".join(
        [
            "# Cartesian same-gamma shock-density-bubble MPhil project package",
            "",
            "This package uses a single-material ideal-gas shock-density-bubble test",
            "in Cartesian coordinates. It keeps the same Rankine-Hugoniot shock",
            "state and low-density bubble geometry used by the Clawpack-style setup,",
            "but it does not include cylindrical `1/r` geometric source terms.",
            "",
            "## Claim scope",
            "",
            "- Qualitative high-speed 2D robustness/validation evidence for AMReX Euler schemes.",
            "- Single material, gamma=1.4.",
            "- Cartesian coordinates; no cylindrical geometric source terms.",
            "- No GFM, no level set, no fractional material fields, no variable gamma.",
            "- Not an air-helium/GFM validation and not cylindrical/axisymmetric IMEX validation.",
            "",
            "## Figures",
            "",
            "- `figures/shock_density_bubble_cartesian_three_scheme_density_320x80.gif`: density animation.",
            "- `figures/shock_density_bubble_cartesian_three_scheme_density_snapshots_320x80.png`: report-use contact sheet.",
            "- `physics_check/shock_density_bubble_cartesian_physics_check_320x80.md`: Rankine-Hugoniot shock-speed sanity check.",
            "",
            "## Report-grid summary",
            "",
            "| Scheme | Status | Steps | Wall time (s) | rho_min | pressure_min | internal_energy_min |",
            "|---|---:|---:|---:|---:|---:|---:|",
            *summaries,
            "",
            "The IMEX row is a Cartesian ideal-gas high-speed stress result. It should",
            "not be used as evidence for cylindrical/axisymmetric IMEX shock-bubble",
            "closure, because that source-aware pressure split has not been derived.",
            "",
        ]
    )
    (out_dir / "README.md").write_text(text)


def main() -> None:
    figures_dir = PACKAGE_ROOT / "figures"
    data_dir = PACKAGE_ROOT / "data"
    log_dir = PACKAGE_ROOT / "logs"
    physics_dir = PACKAGE_ROOT / "physics_check"
    frame_dir = figures_dir / "frames"
    for path in [figures_dir, data_dir, log_dir, physics_dir, frame_dir, EVIDENCE_FIGURES, EVIDENCE_DATA]:
        path.mkdir(parents=True, exist_ok=True)

    series = collect_snapshots()
    frame_paths = []
    for idx in range(6):
        frame_path = frame_dir / f"shock_density_bubble_cartesian_density_frame_{idx:02d}.png"
        render_frame(series, idx, frame_path)
        frame_paths.append(frame_path)

    gif_path = figures_dir / "shock_density_bubble_cartesian_three_scheme_density_320x80.gif"
    images = [Image.open(path).convert("P", palette=Image.Palette.ADAPTIVE) for path in frame_paths]
    images[0].save(gif_path, save_all=True, append_images=images[1:], duration=650, loop=0)

    sheet_path = figures_dir / "shock_density_bubble_cartesian_three_scheme_density_snapshots_320x80.png"
    render_contact_sheet(series, sheet_path)

    check_csv = physics_dir / "shock_density_bubble_cartesian_physics_check_320x80.csv"
    check_fig = figures_dir / "shock_density_bubble_cartesian_physics_check_320x80.png"
    check_note = physics_dir / "shock_density_bubble_cartesian_physics_check_320x80.md"
    render_physics_check(series, check_csv, check_fig, check_note)

    for scheme in SCHEMES:
        shutil.copy2(scheme["summary"], data_dir / scheme["summary"].name)
        log_path = RESULT_ROOT / scheme["summary"].name.replace("_summary.csv", ".log")
        if log_path.exists():
            shutil.copy2(log_path, log_dir / log_path.name)
        shutil.copy2(scheme["summary"], EVIDENCE_DATA / scheme["summary"].name)

    for path in [gif_path, sheet_path, check_fig]:
        shutil.copy2(path, EVIDENCE_FIGURES / path.name)
    shutil.copy2(check_csv, EVIDENCE_DATA / check_csv.name)

    write_readme(PACKAGE_ROOT)
    print(f"wrote {gif_path}")
    print(f"wrote {sheet_path}")
    print(f"wrote {check_fig}")
    print(f"wrote {PACKAGE_ROOT / 'README.md'}")


if __name__ == "__main__":
    main()
