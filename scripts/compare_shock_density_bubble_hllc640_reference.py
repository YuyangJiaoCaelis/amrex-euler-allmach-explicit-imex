#!/usr/bin/env python3
"""Compare 320x80 Cartesian shock-density-bubble runs to a 640x160 HLLC reference."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
COARSE_ROOT = REPO / "results/amrex/project_same_gamma_shock_density_bubble_report_grid_2026-05-20"
REFERENCE_ROOT = REPO / "results/amrex/project_same_gamma_shock_density_bubble_hires_reference_2026-05-20"
PACKAGE = REPO / "project_outputs/test_packages/shock_density_bubble_cartesian_same_gamma"
EVIDENCE_DATA = REPO / "project_outputs/project_evidence/data"
EVIDENCE_FIGURES = REPO / "project_outputs/project_evidence/figures"

REFERENCE = {
    "label": "Explicit O2 HLLC 640x160 reference",
    "snapshots": REFERENCE_ROOT / "hllc_o2_640x160_snapshots",
}

SCHEMES = [
    ("Explicit O2 HLLC", COARSE_ROOT / "hllc_o2_320x80_snapshots"),
    ("Explicit O2 Low-Mach Corrected HLLC-P", COARSE_ROOT / "lowmach_corrected_hllc_p_o2_320x80_snapshots"),
    ("IMEX T1/S2 BDLTV20", COARSE_ROOT / "imex_t1s2_bdltv20_320x80_snapshots"),
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
    return {
        "time": float(rows[0]["time"]),
        "xs": xs,
        "ys": ys,
        "fields": fields,
    }


def average_2x2_to_coarse(field: np.ndarray) -> np.ndarray:
    ny, nx = field.shape
    if ny % 2 != 0 or nx % 2 != 0:
        raise ValueError(f"expected even shape, got {field.shape}")
    return field.reshape(ny // 2, 2, nx // 2, 2).mean(axis=(1, 3))


def mirrored(field: np.ndarray) -> np.ndarray:
    return np.vstack([field[::-1, :], field])


def load_series(snap_dir: Path) -> list[dict[str, object]]:
    return [load_snapshot(path) for path in sorted(snap_dir.glob("snapshot_*.csv"), key=snapshot_index)]


def relative_l1(diff: np.ndarray, ref: np.ndarray) -> float:
    return float(np.mean(np.abs(diff)) / max(np.mean(np.abs(ref)), 1.0e-300))


def linf(diff: np.ndarray) -> float:
    return float(np.max(np.abs(diff)))


def main() -> None:
    out_dir = PACKAGE / "physics_check"
    fig_dir = PACKAGE / "figures"
    for path in [out_dir, fig_dir, EVIDENCE_DATA, EVIDENCE_FIGURES]:
        path.mkdir(parents=True, exist_ok=True)

    ref_series = load_series(REFERENCE["snapshots"])
    coarse_ref = []
    for snap in ref_series:
        coarse_ref.append(
            {
                "time": snap["time"],
                "rho": average_2x2_to_coarse(snap["fields"]["rho"]),
                "pressure": average_2x2_to_coarse(snap["fields"]["pressure"]),
            }
        )

    rows: list[dict[str, object]] = []
    final_fields: dict[str, dict[str, np.ndarray]] = {}
    for label, snap_dir in SCHEMES:
        series = load_series(snap_dir)
        if len(series) != len(coarse_ref):
            raise RuntimeError(f"{label} snapshot count mismatch")
        for idx, (snap, ref) in enumerate(zip(series, coarse_ref)):
            if abs(float(snap["time"]) - float(ref["time"])) > 1.0e-12:
                raise RuntimeError(f"{label} time mismatch at snapshot {idx}")
            rho = snap["fields"]["rho"]
            pressure = snap["fields"]["pressure"]
            rho_diff = rho - ref["rho"]
            p_diff = pressure - ref["pressure"]
            rows.append(
                {
                    "scheme": label,
                    "time": float(snap["time"]),
                    "rho_relative_l1_vs_hllc640": relative_l1(rho_diff, ref["rho"]),
                    "pressure_relative_l1_vs_hllc640": relative_l1(p_diff, ref["pressure"]),
                    "rho_linf_vs_hllc640": linf(rho_diff),
                    "pressure_linf_vs_hllc640": linf(p_diff),
                }
            )
            if idx == len(series) - 1:
                final_fields[label] = {
                    "rho": rho,
                    "pressure": pressure,
                    "rho_diff": rho_diff,
                    "pressure_diff": p_diff,
                }
    final_fields[REFERENCE["label"]] = {
        "rho": coarse_ref[-1]["rho"],
        "pressure": coarse_ref[-1]["pressure"],
        "rho_diff": np.zeros_like(coarse_ref[-1]["rho"]),
        "pressure_diff": np.zeros_like(coarse_ref[-1]["pressure"]),
    }

    csv_path = out_dir / "shock_density_bubble_cartesian_hllc640_reference_comparison.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (EVIDENCE_DATA / csv_path.name).write_bytes(csv_path.read_bytes())

    colors = {
        "Explicit O2 HLLC": "tab:blue",
        "Explicit O2 Low-Mach Corrected HLLC-P": "#9467bd",
        "IMEX T1/S2 BDLTV20": "tab:green",
    }
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    for label, _ in SCHEMES:
        rr = [row for row in rows if row["scheme"] == label]
        t = np.array([float(row["time"]) for row in rr])
        axes[0].plot(
            t,
            [float(row["rho_relative_l1_vs_hllc640"]) for row in rr],
            "o-",
            color=colors[label],
            label=label,
        )
        axes[1].plot(
            t,
            [float(row["pressure_relative_l1_vs_hllc640"]) for row in rr],
            "o-",
            color=colors[label],
            label=label,
        )
    axes[0].set_title("Density difference from averaged HLLC 640x160")
    axes[0].set_ylabel("relative L1")
    axes[1].set_title("Pressure difference from averaged HLLC 640x160")
    axes[1].set_ylabel("relative L1")
    for ax in axes:
        ax.set_xlabel("time")
        ax.grid(True, alpha=0.25)
    axes[1].legend(fontsize=8)
    metric_fig = fig_dir / "shock_density_bubble_cartesian_hllc640_reference_errors.png"
    fig.savefig(metric_fig, dpi=220)
    plt.close(fig)
    (EVIDENCE_FIGURES / metric_fig.name).write_bytes(metric_fig.read_bytes())

    # Final-time field comparison: show reference, coarse fields, and density differences.
    order = [REFERENCE["label"], *(label for label, _ in SCHEMES)]
    fig, axes = plt.subplots(2, 4, figsize=(17.0, 6.7), constrained_layout=False)
    fig.subplots_adjust(left=0.045, right=0.90, top=0.88, bottom=0.10, wspace=0.25, hspace=0.38)
    vmax = 2.82
    vmin = 0.1
    max_abs_diff = max(float(np.max(np.abs(final_fields[label]["rho_diff"]))) for label in order[1:])
    for col, label in enumerate(order):
        rho = mirrored(final_fields[label]["rho"])
        im = axes[0, col].imshow(
            rho,
            origin="lower",
            extent=(0.0, 2.0, -0.5, 0.5),
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
            aspect="equal",
        )
        axes[0, col].set_title(label, fontsize=9)
        axes[0, col].set_xlim(0.0, 1.35)
        axes[0, col].set_ylim(-0.32, 0.32)
        axes[0, col].set_xlabel("x")
        axes[0, col].set_ylabel("y")
        axes[0, col].axhline(0.0, color="white", linewidth=0.5, alpha=0.45)
        diff = mirrored(final_fields[label]["rho_diff"])
        im_diff = axes[1, col].imshow(
            diff,
            origin="lower",
            extent=(0.0, 2.0, -0.5, 0.5),
            cmap="coolwarm",
            vmin=-max_abs_diff,
            vmax=max_abs_diff,
            interpolation="nearest",
            aspect="equal",
        )
        axes[1, col].set_title("density - HLLC 640 avg", fontsize=9)
        axes[1, col].set_xlim(0.0, 1.35)
        axes[1, col].set_ylim(-0.32, 0.32)
        axes[1, col].set_xlabel("x")
        axes[1, col].set_ylabel("y")
        axes[1, col].axhline(0.0, color="black", linewidth=0.5, alpha=0.35)
    cax_density = fig.add_axes([0.92, 0.56, 0.015, 0.28])
    cax_diff = fig.add_axes([0.92, 0.15, 0.015, 0.28])
    fig.colorbar(im, cax=cax_density, label="Density")
    fig.colorbar(im_diff, cax=cax_diff, label="Density difference")
    field_fig = fig_dir / "shock_density_bubble_cartesian_hllc640_reference_final_density.png"
    fig.savefig(field_fig, dpi=220)
    plt.close(fig)
    (EVIDENCE_FIGURES / field_fig.name).write_bytes(field_fig.read_bytes())

    final_rows = [row for row in rows if abs(float(row["time"]) - 0.3) < 1.0e-12]
    note = [
        "# Cartesian shock-density-bubble HLLC 640x160 reference comparison",
        "",
        "This analysis compares the existing 320x80 AMReX Cartesian same-gamma",
        "shock-density-bubble rows against a newly generated 640x160 Explicit",
        "O2 HLLC reference, averaged back to 320x80. It is a grid-reference",
        "sanity check for the genuinely 2D interaction, not an exact solution.",
        "",
        "## Reference run",
        "",
        "- scheme: Explicit O2 HLLC",
        "- grid: 640x160",
        "- status: ok",
        "- final time: 0.3",
        "- pressure_min: 1.0",
        "- rho_min: 0.25113",
        "",
        "## Final-time relative L1 difference from averaged HLLC 640x160",
        "",
        "| Scheme | density rel L1 | pressure rel L1 | density Linf | pressure Linf |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in final_rows:
        note.append(
            f"| {row['scheme']} | {float(row['rho_relative_l1_vs_hllc640']):.6g} | "
            f"{float(row['pressure_relative_l1_vs_hllc640']):.6g} | "
            f"{float(row['rho_linf_vs_hllc640']):.6g} | "
            f"{float(row['pressure_linf_vs_hllc640']):.6g} |"
        )
    note.extend(
        [
            "",
            "Interpretation: the explicit 320x80 HLLC row is the closest coarse row",
            "to the high-resolution HLLC reference, as expected. The low-Mach",
            "corrected explicit row is very similar in this supersonic case. The",
            "IMEX T1/S2 row remains positive and captures the same large-scale",
            "shock/bubble motion, but has larger density/pressure differences near",
            "the compressed bubble and transmitted wave. This is acceptable as a",
            "MPhil project qualitative high-speed stress check, but it should not be",
            "sold as final shock-bubble accuracy or efficiency validation.",
            "",
        ]
    )
    note_path = out_dir / "shock_density_bubble_cartesian_hllc640_reference_comparison.md"
    note_path.write_text("\n".join(note))

    print(csv_path)
    print(metric_fig)
    print(field_fig)
    print(note_path)


if __name__ == "__main__":
    main()
