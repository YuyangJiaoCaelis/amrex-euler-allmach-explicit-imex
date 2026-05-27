#!/usr/bin/env python3
"""Shock-density-bubble error-vs-wall-time analysis against HLLC 640x160.

The grids are restricted to sizes that coarsen exactly from the 640x160
reference, so the comparison uses conservative block averages rather than
interpolation.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
PACKAGE = REPO / "project_outputs/test_packages/shock_density_bubble_cartesian_same_gamma"
EVIDENCE_DATA = REPO / "project_outputs/project_evidence/data"
EVIDENCE_FIGURES = REPO / "project_outputs/project_evidence/figures"
REFERENCE_DIR = (
    REPO
    / "results/amrex/project_same_gamma_shock_density_bubble_hires_reference_2026-05-20"
    / "hllc_o2_640x160_snapshots"
)

SCHEMES = [
    {
        "label": "Explicit O2 HLLC",
        "color": "tab:blue",
        "rows": [
            (
                128,
                32,
                REPO / "results/amrex/same_gamma_shock_density_bubble_efficiency_divisible_2026-05-20",
                "hllc_o2",
            ),
            (
                160,
                40,
                REPO / "results/amrex/project_same_gamma_shock_density_bubble_smoke_2026-05-20",
                "hllc_o2",
            ),
            (
                320,
                80,
                REPO / "results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20",
                "hllc_o2",
            ),
        ],
    },
    {
        "label": "Explicit O2 Low-Mach Corrected HLLC-P",
        "color": "#9467bd",
        "rows": [
            (
                128,
                32,
                REPO / "results/amrex/same_gamma_shock_density_bubble_efficiency_divisible_2026-05-20",
                "lowmach_corrected_hllc_p_o2",
            ),
            (
                160,
                40,
                REPO / "results/amrex/project_same_gamma_shock_density_bubble_smoke_2026-05-20",
                "lowmach_corrected_hllc_p_o2",
            ),
            (
                320,
                80,
                REPO / "results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20",
                "lowmach_corrected_hllc_p_o2",
            ),
        ],
    },
    {
        "label": "IMEX T1/S2 BDLTV20",
        "color": "tab:green",
        "rows": [
            (
                128,
                32,
                REPO / "results/amrex/same_gamma_shock_density_bubble_efficiency_divisible_2026-05-20",
                "imex_t1s2_bdltv20",
            ),
            (
                160,
                40,
                REPO / "results/amrex/project_same_gamma_shock_density_bubble_smoke_2026-05-20",
                "imex_t1s2_bdltv20",
            ),
            (
                320,
                80,
                REPO / "results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20",
                "imex_t1s2_bdltv20",
            ),
        ],
    },
]


def snapshot_index(path: Path) -> int:
    match = re.search(r"snapshot_(\d+)_", path.name)
    if not match:
        raise ValueError(f"cannot parse snapshot index from {path}")
    return int(match.group(1))


def load_snapshot(path: Path) -> dict[str, object]:
    rows = list(csv.DictReader(path.open(newline="")))
    if not rows:
        raise ValueError(f"empty snapshot {path}")
    xs = np.array(sorted({float(row["x"]) for row in rows}))
    ys = np.array(sorted({float(row["y"]) for row in rows}))
    xi = {x: i for i, x in enumerate(xs)}
    yi = {y: j for j, y in enumerate(ys)}
    fields = {
        "rho": np.empty((len(ys), len(xs))),
        "pressure": np.empty((len(ys), len(xs))),
    }
    for row in rows:
        i = xi[float(row["x"])]
        j = yi[float(row["y"])]
        fields["rho"][j, i] = float(row["rho"])
        fields["pressure"][j, i] = float(row["pressure"])
    return {"time": float(rows[0]["time"]), "fields": fields}


def average_to_shape(field: np.ndarray, ny: int, nx: int) -> np.ndarray:
    ref_ny, ref_nx = field.shape
    if ref_ny % ny != 0 or ref_nx % nx != 0:
        raise ValueError(f"reference shape {field.shape} does not coarsen to {(ny, nx)}")
    fy = ref_ny // ny
    fx = ref_nx // nx
    return field.reshape(ny, fy, nx, fx).mean(axis=(1, 3))


def relative_l1(field: np.ndarray, ref: np.ndarray) -> float:
    return float(np.mean(np.abs(field - ref)) / max(np.mean(np.abs(ref)), 1.0e-300))


def load_summary(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return next(csv.DictReader(handle))


def log_fit(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    if np.count_nonzero(mask) < 2:
        return None
    coeff = np.polyfit(np.log10(x[mask]), np.log10(y[mask]), 1)
    xx = np.geomspace(float(np.min(x[mask])), float(np.max(x[mask])), 80)
    yy = 10.0 ** np.polyval(coeff, np.log10(xx))
    return xx, yy


def main() -> None:
    out_dir = PACKAGE / "physics_check"
    fig_dir = PACKAGE / "figures"
    for path in [out_dir, fig_dir, EVIDENCE_DATA, EVIDENCE_FIGURES]:
        path.mkdir(parents=True, exist_ok=True)

    ref_final = load_snapshot(sorted(REFERENCE_DIR.glob("snapshot_*.csv"), key=snapshot_index)[-1])
    rows: list[dict[str, object]] = []
    for scheme in SCHEMES:
        for nx, ny, root, prefix in scheme["rows"]:
            snap_dir = root / f"{prefix}_{nx}x{ny}_snapshots"
            summary_path = root / f"{prefix}_{nx}x{ny}_summary.csv"
            final = load_snapshot(sorted(snap_dir.glob("snapshot_*.csv"), key=snapshot_index)[-1])
            summary = load_summary(summary_path)
            rho_ref = average_to_shape(ref_final["fields"]["rho"], ny, nx)
            p_ref = average_to_shape(ref_final["fields"]["pressure"], ny, nx)
            rows.append(
                {
                    "scheme": scheme["label"],
                    "grid_nx": nx,
                    "grid_ny": ny,
                    "wall_time_sec": float(summary["wall_time_sec"]),
                    "steps_accepted": int(summary["steps_accepted"]),
                    "rho_relative_l1_vs_hllc640": relative_l1(final["fields"]["rho"], rho_ref),
                    "pressure_relative_l1_vs_hllc640": relative_l1(
                        final["fields"]["pressure"], p_ref
                    ),
                    "status": summary["status"],
                    "failure_category": summary["failure_category"],
                    "rho_min": float(summary["rho_min"]),
                    "pressure_min": float(summary["pressure_min"]),
                    "internal_energy_min": float(summary["internal_energy_min"]),
                }
            )

    csv_path = out_dir / "shock_density_bubble_cartesian_hllc640_efficiency_comparison.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (EVIDENCE_DATA / csv_path.name).write_bytes(csv_path.read_bytes())

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.8), constrained_layout=True)
    metrics = [
        ("rho_relative_l1_vs_hllc640", "Density relative L1"),
        ("pressure_relative_l1_vs_hllc640", "Pressure relative L1"),
    ]
    for ax, (metric, title) in zip(axes, metrics):
        for scheme in SCHEMES:
            sr = [row for row in rows if row["scheme"] == scheme["label"]]
            x = np.array([float(row["wall_time_sec"]) for row in sr])
            y = np.array([float(row[metric]) for row in sr])
            ax.scatter(x, y, color=scheme["color"], s=42, label=scheme["label"])
            fit = log_fit(x, y)
            if fit is not None:
                ax.plot(fit[0], fit[1], color=scheme["color"], linewidth=1.8)
            for row, xx, yy in zip(sr, x, y):
                ax.annotate(
                    f"{row['grid_nx']}",
                    (xx, yy),
                    textcoords="offset points",
                    xytext=(5, 4),
                    fontsize=8,
                    color=scheme["color"],
                )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Wall time (s)")
        ax.set_ylabel(title)
        ax.set_title(title + " vs wall time")
        ax.grid(True, which="both", alpha=0.25)
    axes[1].legend(fontsize=8, loc="best")
    fig.suptitle(
        "AMReX same-gamma shock-density-bubble: final-time error vs wall time\n"
        "markers are grids 128, 160, 320; solid lines are log-log least-square fits",
        fontsize=12,
    )
    fig_path = fig_dir / "shock_density_bubble_cartesian_hllc640_efficiency.png"
    fig.savefig(fig_path, dpi=220)
    plt.close(fig)
    (EVIDENCE_FIGURES / fig_path.name).write_bytes(fig_path.read_bytes())

    final_320 = [row for row in rows if int(row["grid_nx"]) == 320]
    note = [
        "# Cartesian shock-density-bubble error-vs-wall-time analysis",
        "",
        "This analysis uses the same-gamma Cartesian shock-density-bubble rows and",
        "compares final-time density and pressure against the averaged 640x160",
        "Explicit O2 HLLC numerical reference. Only grids that coarsen exactly",
        "from 640x160 are used: 128x32, 160x40, and 320x80.",
        "",
        "This is a preliminary cost/error trend, not a final shock-bubble IMEX",
        "efficiency claim. The reference is HLLC-family and the interaction has",
        "no closed-form exact solution.",
        "",
        "## 320x80 headline values",
        "",
        "| Scheme | wall time (s) | density rel L1 | pressure rel L1 |",
        "|---|---:|---:|---:|",
    ]
    for row in final_320:
        note.append(
            f"| {row['scheme']} | {float(row['wall_time_sec']):.3f} | "
            f"{float(row['rho_relative_l1_vs_hllc640']):.6g} | "
            f"{float(row['pressure_relative_l1_vs_hllc640']):.6g} |"
        )
    note.extend(
        [
            "",
            "Interpretation: the IMEX row is close to the explicit rows",
            "in final-time relative L1 for this same-gamma high-speed stress case,",
            "but it is much more expensive on the current serial CPU path because",
            "each time step assembles and solves the pressure system. This supports",
            "using the case as qualitative high-speed evidence, while reserving",
            "shock-bubble efficiency conclusions for later validated CPU/GPU runs.",
            "",
        ]
    )
    note_path = out_dir / "shock_density_bubble_cartesian_hllc640_efficiency_comparison.md"
    note_path.write_text("\n".join(note))

    print(csv_path)
    print(fig_path)
    print(note_path)


if __name__ == "__main__":
    main()
