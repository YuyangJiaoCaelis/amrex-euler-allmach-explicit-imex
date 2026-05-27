#!/usr/bin/env python3
"""Create a report-ready AMReX advection-blob density time-snapshot figure.

The figure uses the schemes reported in MPhil project:

- Explicit O2 HLLC
- Explicit O2 Low-Mach Corrected HLLC-P
- IMEX T1/S2 BDLTV20

It runs only the missing n=128 periodic-blob intermediate snapshots and reuses
existing final-time fields when available.
"""

from __future__ import annotations

import csv
import math
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex"
INPUT = ROOT / "amrex/apps/euler_compare/inputs-ci"
OUT = ROOT / "results/amrex/project_advection_blob_timesnaps_current_2026-05-19"
PACKAGE = ROOT / "project_outputs/test_packages/advection_blob"
REPORT_FIGS = ROOT / "project_outputs/report_figures_cued_template"
N = 128
SNAPSHOT_TIMES = (0.0, 0.25, 0.5)
ERROR_TIMES = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
BLOB_CX = 0.35
BLOB_CY = 0.35
BLOB_RADIUS = 0.15
BLOB_U = 1.0
BLOB_V = 0.5
GAMMA = 1.4


@dataclass(frozen=True)
class Scheme:
    name: str
    short: str
    args: tuple[str, ...]
    existing_final_csv: Path | None = None


SCHEMES = (
    Scheme(
        "Explicit O2 HLLC",
        "explicit_o2_hllc",
        (
            "euler.method=explicit",
            "euler.problem=advection_blob",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            "euler.riemann=hllc",
            "euler.cfl=0.4",
            "max_step=200000",
        ),
        ROOT
        / "results/amrex/advection_blob_explicit_imex_efficiency_refresh_2026-05-17/fields/advection_blob_explicit_o2_hllc_n128_t0p5.csv",
    ),
    Scheme(
        "Explicit O2 Low-Mach Corrected HLLC-P",
        "explicit_o2_lowmachcorrected_hllcp",
        (
            "euler.method=explicit",
            "euler.problem=advection_blob",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            "euler.riemann=xie_am_hllc_p",
            "euler.cfl=0.4",
            "max_step=200000",
        ),
        ROOT
        / "results/amrex/advection_blob_explicit_imex_efficiency_refresh_2026-05-17/fields/advection_blob_explicit_o2_lowmachcorrected_hllcp_n128_t0p5.csv",
    ),
    Scheme(
        "IMEX T1/S2 BDLTV20",
        "imex_t1s2_bdltv20",
        (
            "euler.method=imex",
            "euler.problem=advection_blob",
            "euler.bdltv20_paper_t1_s2=advection_blob_periodic_2d",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            "euler.imex_picard_iterations=4",
            "euler.imex_acoustic_startup=0",
            "euler.imex_acoustic_cfl_cap=0.0",
            "euler.imex_pressure_stabilization=off",
            "euler.imex_predictor_dissipation=material",
            "euler.imex_solver_tol=1e-8",
            "euler.imex_solver_max_iter=1000",
            "euler.bdltv20_paper_pressure_solver=gmres",
            "euler.bdltv20_paper_t1_s2_dt=-1",
            "euler.bdltv20_paper_t1_s2_max_steps=200000",
            "euler.imex_cfl=0.4",
        ),
        ROOT
        / "results/amrex/current_bdltv20_t1s2_project_refresh_tol1e8_2026-05-17/advection/fields/advection_blob_imex_bdltv20_paper_t1s2_n128_t0p5.csv",
    ),
)


def time_label(t: float) -> str:
    return f"{t:g}".replace(".", "p")


def exact_density(x: np.ndarray, y: np.ndarray, t: float) -> np.ndarray:
    cx = (BLOB_CX + BLOB_U * t) % 1.0
    cy = (BLOB_CY + BLOB_V * t) % 1.0
    dx = np.minimum(np.abs(x - cx), 1.0 - np.abs(x - cx))
    dy = np.minimum(np.abs(y - cy), 1.0 - np.abs(y - cy))
    r2 = dx * dx + dy * dy
    return np.where(r2 <= BLOB_RADIUS * BLOB_RADIUS, 2.0, 1.0)


def command_for(scheme: Scheme, t: float, csv_path: Path) -> list[str]:
    cmd = [
        str(EXE),
        str(INPUT),
        f"amr.n_cell={N} {N}",
        f"amr.max_grid_size={N}",
        "amr.plot_int=-1",
        "geometry.is_periodic=1 1",
        "geometry.prob_lo=0.0 0.0",
        "geometry.prob_hi=1.0 1.0",
        f"euler.blob_center={BLOB_CX} {BLOB_CY}",
        f"euler.blob_radius={BLOB_RADIUS}",
        f"euler.velocity={BLOB_U} {BLOB_V}",
        *scheme.args,
        f"stop_time={t:.12g}",
        f"euler.final_csv={csv_path}",
    ]
    return cmd


def run_snapshot(scheme: Scheme, t: float) -> Path:
    if t == 0.5 and scheme.existing_final_csv and scheme.existing_final_csv.exists():
        return scheme.existing_final_csv

    stem = f"advection_blob_{scheme.short}_n{N}_t{time_label(t)}"
    csv_path = OUT / "fields" / f"{stem}.csv"
    log_path = OUT / "logs" / f"{stem}.log"
    cmd_path = OUT / "commands" / f"{stem}.txt"
    if csv_path.exists():
        return csv_path

    cmd = command_for(scheme, t, csv_path)
    cmd_path.write_text(" ".join(subprocess.list2cmdline([p]) for p in cmd) + "\n")
    start = time.perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    elapsed = time.perf_counter() - start
    log_path.write_text(proc.stdout)
    (OUT / "logs" / f"{stem}.driver_wall_time").write_text(f"{elapsed:.12g}\n")
    if proc.returncode != 0 or not csv_path.exists():
        raise RuntimeError(f"{scheme.name} t={t} failed; see {log_path}")
    return csv_path


def read_field(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    rho: list[float] = []
    exact: list[float] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))
            rho.append(float(row["rho"]))
            exact_key = "exact_rho" if "exact_rho" in row else None
            exact.append(float(row[exact_key]) if exact_key else math.nan)

    xvals = np.array(sorted(set(xs)))
    yvals = np.array(sorted(set(ys)))
    nx = len(xvals)
    ny = len(yvals)
    return (
        xvals,
        yvals,
        np.array(rho).reshape((ny, nx)),
        np.array(exact).reshape((ny, nx)),
    )


def write_summary(rows: list[dict[str, str]]) -> None:
    summary = OUT / "advection_blob_current_timesnap_summary.csv"
    with summary.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    package_copy = PACKAGE / "source_data" / summary.name
    package_copy.write_text(summary.read_text())


def make_figure(rows: list[dict[str, str]]) -> None:
    fields: dict[tuple[str, float], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    exact_fields: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for row in rows:
        t = float(row["time"])
        x, y, rho, exact = read_field(Path(row["source_csv"]))
        fields[(row["scheme"], t)] = (x, y, rho)
        if t not in exact_fields:
            if np.isfinite(exact).all():
                exact_fields[t] = (x, y, exact)
            else:
                xx, yy = np.meshgrid(x, y)
                exact_fields[t] = (x, y, exact_density(xx, yy, t))

    row_labels = ["Exact", *(s.name for s in SCHEMES)]
    field_keys = ["Exact", *(s.name for s in SCHEMES)]
    fig, axes = plt.subplots(
        len(row_labels),
        len(SNAPSHOT_TIMES),
        figsize=(9.4, 8.4),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    cmap = "viridis"
    vmin, vmax = 1.0, 2.0
    for c, t in enumerate(SNAPSHOT_TIMES):
        for r, label in enumerate(row_labels):
            ax = axes[r, c]
            if r == 0:
                x, y, data = exact_fields[t]
            else:
                x, y, data = fields[(field_keys[r], t)]
            im = ax.imshow(
                data,
                origin="lower",
                extent=(x.min(), x.max(), y.min(), y.max()),
                vmin=vmin,
                vmax=vmax,
                cmap=cmap,
                interpolation="nearest",
                aspect="equal",
            )
            ax.contour(
                x,
                y,
                data,
                levels=[1.5],
                colors="white",
                linewidths=0.8,
                alpha=0.9,
            )
            if r == 0:
                ax.set_title(f"$t={t:g}$")
            if c == 0:
                ax.set_ylabel(label, fontsize=9)
            ax.set_xticks([0.0, 0.5, 1.0])
            ax.set_yticks([0.0, 0.5, 1.0])
    fig.supxlabel("$x$")
    fig.supylabel("$y$")
    cbar = fig.colorbar(im, ax=axes, location="right", shrink=0.82)
    cbar.set_label("Density")

    out = PACKAGE / "figures" / "amrex_advection_blob_density_timesnaps_current.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    report_out = REPORT_FIGS / out.name
    report_out.write_bytes(out.read_bytes())
    evidence = ROOT / "project_outputs/project_evidence/figures" / out.name
    evidence.write_bytes(out.read_bytes())
    plt.close(fig)


def make_error_time_figure(rows: list[dict[str, str]]) -> None:
    colors = {
        "Explicit O2 HLLC": "#1f77b4",
        "Explicit O2 Low-Mach Corrected HLLC-P": "#9467bd",
        "IMEX T1/S2 BDLTV20": "#2ca02c",
    }
    labels = {
        "Explicit O2 HLLC": "Explicit O2 HLLC",
        "Explicit O2 Low-Mach Corrected HLLC-P": "Explicit O2 Low-Mach HLLC-P",
        "IMEX T1/S2 BDLTV20": "IMEX T1/S2 BDLTV20",
    }
    by_scheme: dict[str, list[tuple[float, float]]] = {s.name: [] for s in SCHEMES}
    for row in rows:
        t = float(row["time"])
        if t not in ERROR_TIMES:
            continue
        _, _, rho, exact = read_field(Path(row["source_csv"]))
        if not np.isfinite(exact).all():
            x, y, _, _ = read_field(Path(row["source_csv"]))
            xx, yy = np.meshgrid(x, y)
            exact = exact_density(xx, yy, t)
        density_l1 = float(np.mean(np.abs(rho - exact)))
        by_scheme[row["scheme"]].append((t, density_l1))

    fig, ax = plt.subplots(figsize=(7.1, 4.25), constrained_layout=True)
    draw_order = (
        "Explicit O2 Low-Mach Corrected HLLC-P",
        "Explicit O2 HLLC",
        "IMEX T1/S2 BDLTV20",
    )
    styles = {
        "Explicit O2 HLLC": dict(marker="o", linestyle="-", linewidth=1.4, markersize=4.2, zorder=3),
        "Explicit O2 Low-Mach Corrected HLLC-P": dict(marker="s", linestyle="--", linewidth=1.5, markersize=4.0, zorder=2),
        "IMEX T1/S2 BDLTV20": dict(marker="o", linestyle="-", linewidth=1.6, markersize=4.5, zorder=4),
    }
    for scheme_name in draw_order:
        values = sorted(by_scheme[scheme_name])
        times = [v[0] for v in values]
        errors = [v[1] for v in values]
        ax.plot(
            times,
            errors,
            color=colors[scheme_name],
            label=labels[scheme_name],
            **styles[scheme_name],
        )
    ax.set_xlabel("Physical time")
    ax.set_ylabel(r"Density $L^1$ error")
    ax.grid(True, alpha=0.28)
    ax.legend(frameon=True, fontsize=8)

    out = PACKAGE / "figures" / "amrex_advection_blob_density_error_vs_time_current.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    report_out = REPORT_FIGS / out.name
    report_out.write_bytes(out.read_bytes())
    evidence = ROOT / "project_outputs/project_evidence/figures" / out.name
    evidence.write_bytes(out.read_bytes())
    plt.close(fig)


def update_manifest() -> None:
    manifest = PACKAGE / "MANIFEST.csv"
    existing = manifest.read_text() if manifest.exists() else "path,role,source\n"
    additions = [
        (
            PACKAGE / "figures/amrex_advection_blob_density_timesnaps_current.png",
            "report figure",
            "current AMReX advection blob density time snapshots",
        ),
        (
            PACKAGE / "source_data/advection_blob_current_timesnap_summary.csv",
            "source data",
            "current AMReX advection blob density time snapshots",
        ),
        (
            PACKAGE / "figures/amrex_advection_blob_density_error_vs_time_current.png",
            "report figure",
            "current AMReX advection blob density L1 error against physical time",
        ),
    ]
    with manifest.open("a") as f:
        for path, role, source in additions:
            if str(path) not in existing:
                f.write(f"{path},{role},{source}\n")


def main() -> None:
    for sub in ("fields", "logs", "commands"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    (PACKAGE / "figures").mkdir(parents=True, exist_ok=True)
    (PACKAGE / "source_data").mkdir(parents=True, exist_ok=True)
    REPORT_FIGS.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    times = tuple(sorted(set(SNAPSHOT_TIMES).union(ERROR_TIMES)))
    for scheme in SCHEMES:
        for t in times:
            csv_path = run_snapshot(scheme, t)
            rows.append(
                {
                    "scheme": scheme.name,
                    "selector": scheme.short,
                    "n": str(N),
                    "time": f"{t:g}",
                    "source_csv": str(csv_path),
                }
            )
    write_summary(rows)
    make_figure(rows)
    make_error_time_figure(rows)
    update_manifest()
    print(PACKAGE / "figures/amrex_advection_blob_density_timesnaps_current.png")


if __name__ == "__main__":
    main()
