#!/usr/bin/env python3
"""Run and package project 1D Riemann convergence evidence under AMReX.

The package is intentionally narrow:
- exact Dirichlet boundaries in x/y;
- current explicit O2 HLLC and Low-Mach HLLC-P paths;
- current BDLTV20 T1/S2 paper IMEX path;
- discontinuous Riemann refinement evidence only, not a smooth formal order proof.
"""

from __future__ import annotations

import csv
import math
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt

from riemann_exact import RiemannState, sample, solve_star_region
from run_manifest import environment_build_flags, utc_now, write_manifest


ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex"
INPUTS = ROOT / "amrex/apps/euler_compare/inputs-ci"
RESULTS = ROOT / "results/amrex/project_riemann_1d_convergence_current_t1s2_2026-05-17"
PACKAGE = ROOT / "project_outputs/test_packages/riemann_1d_convergence"

GRIDS = [50, 100, 200, 400, 800]
GAMMA = 1.4
Y_CELLS = 5


@dataclass(frozen=True)
class Case:
    key: str
    label: str
    toro_test: int
    final_time: float
    left: RiemannState
    right: RiemannState
    note: str


@dataclass(frozen=True)
class Scheme:
    key: str
    label: str
    method: str
    riemann: str | None
    color: str


CASES = [
    Case(
        "sod_toro1",
        "Sod / Toro 1",
        1,
        0.2,
        RiemannState(1.0, 0.0, 1.0),
        RiemannState(0.125, 0.0, 0.1),
        "standard shock tube",
    ),
    Case(
        "lax_rp2",
        "Lax / BDLTV20 RP2",
        93,
        0.14,
        RiemannState(0.445, 1.698, 3.528),
        RiemannState(0.5, 0.0, 0.571),
        "BDLTV20-style RP2/Lax comparison case",
    ),
    Case(
        "toro3_strong_shock",
        "Toro 3",
        3,
        0.012,
        RiemannState(1.0, 0.0, 1000.0),
        RiemannState(1.0, 0.0, 0.01),
        "strong shock/contact sanity case",
    ),
    Case(
        "toro4_strong_shock",
        "Toro 4",
        4,
        0.035,
        RiemannState(1.0, 0.0, 0.01),
        RiemannState(1.0, 0.0, 100.0),
        "strong shock/contact sanity case",
    ),
]

SCHEMES = [
    Scheme("explicit_o2_hllc", "Explicit O2 HLLC", "explicit", "hllc", "#1f77b4"),
    Scheme(
        "explicit_o2_lowmach_hllcp",
        "Explicit O2 Low-Mach HLLC-P",
        "explicit",
        "xie_am_hllc_p",
        "#9467bd",
    ),
    Scheme("imex_t1s2_bdltv20", "IMEX T1/S2 BDLTV20", "imex", None, "#2ca02c"),
]


def main() -> None:
    require_inputs()
    reset_dir(RESULTS)
    reset_dir(PACKAGE)
    for subdir in [
        RESULTS / "fields",
        RESULTS / "logs",
        RESULTS / "commands",
        RESULTS / "manifests",
        PACKAGE / "figures",
        PACKAGE / "source_data/fields",
        PACKAGE / "source_data/logs",
        PACKAGE / "source_data/commands",
    ]:
        subdir.mkdir(parents=True, exist_ok=True)

    run_rows: list[dict[str, str]] = []
    for case in CASES:
        for scheme in SCHEMES:
            for n in GRIDS:
                run_rows.append(run_case(case, scheme, n))

    write_csv(RESULTS / "run_summary.csv", run_rows)
    metric_rows = build_metric_rows(run_rows)
    write_csv(RESULTS / "riemann_1d_convergence_metrics.csv", metric_rows)
    order_rows = build_order_rows(metric_rows)
    write_csv(RESULTS / "riemann_1d_observed_order_summary.csv", order_rows)

    copy_source_artifacts(run_rows)
    write_csv(PACKAGE / "source_data/riemann_1d_convergence_metrics.csv", metric_rows)
    write_csv(PACKAGE / "source_data/riemann_1d_observed_order_summary.csv", order_rows)
    write_csv(PACKAGE / "source_data/run_summary.csv", run_rows)
    make_plots(metric_rows)
    write_readme(run_rows, order_rows)
    write_readiness_check(run_rows, metric_rows, order_rows)
    write_manifest()
    print(f"Wrote {PACKAGE}")


def require_inputs() -> None:
    if not EXE.exists():
        raise SystemExit(f"Missing executable: {EXE}")
    if not INPUTS.exists():
        raise SystemExit(f"Missing inputs file: {INPUTS}")


def reset_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)


def run_case(case: Case, scheme: Scheme, n: int) -> dict[str, str]:
    stem = f"{case.key}_{scheme.key}_n{n}"
    field_csv = RESULTS / "fields" / f"{stem}.csv"
    log_path = RESULTS / "logs" / f"{stem}.log"
    command_path = RESULTS / "commands" / f"{stem}.command.txt"
    manifest_path = RESULTS / "manifests" / f"{stem}.manifest.json"
    cmd = base_command(case, n, field_csv)
    if scheme.method == "explicit":
        cmd += [
            "max_step=200000",
            "euler.method=explicit",
            f"euler.riemann={scheme.riemann}",
            "euler.cfl=0.4",
        ]
    else:
        imex_steps = imex_step_count(case, n)
        imex_dt = case.final_time / imex_steps
        cmd += [
            f"max_step={imex_steps}",
            "euler.method=imex",
            "euler.imex_cfl=0.95",
            "euler.imex_acoustic_startup=0",
            "euler.imex_acoustic_cfl_cap=0.0",
            "euler.imex_pressure_stabilization=off",
            "euler.imex_predictor_dissipation=material",
            "euler.imex_solver_tol=1e-10",
            "euler.imex_solver_max_iter=1000",
            "euler.imex_picard_iterations=4",
            "euler.bdltv20_paper_t1_s2=toro_xy_exact_dirichlet_2d",
            "euler.bdltv20_paper_pressure_solver=gmres",
            f"euler.bdltv20_paper_t1_s2_dt={imex_dt:.17g}",
            f"euler.bdltv20_paper_t1_s2_max_steps={imex_steps}",
        ]

    command_path.write_text(" ".join(shell_quote(part) for part in cmd) + "\n")
    start_utc = utc_now()
    started = time.perf_counter()
    completed = subprocess.run(
        cmd,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    elapsed = time.perf_counter() - started
    end_utc = utc_now()
    log_path.write_text(completed.stdout)
    write_manifest(
        manifest_path,
        root=ROOT,
        row_id=stem,
        command=cmd,
        start_utc=start_utc,
        end_utc=end_utc,
        wall_time_s=elapsed,
        exit_code=completed.returncode,
        output_root=RESULTS,
        output_class="candidate",
        input_files=[INPUTS],
        output_files=[field_csv, log_path, command_path],
        build_flags=environment_build_flags({"DIM": "2"}),
        notes="1D Riemann refinement row embedded in a 2D AMReX grid",
    )
    log_keys = parse_log(log_path)
    row = {
        "case": case.key,
        "case_label": case.label,
        "toro_test": str(case.toro_test),
        "case_note": case.note,
        "scheme": scheme.label,
        "scheme_key": scheme.key,
        "ncell": str(n),
        "ny": str(Y_CELLS),
        "boundary": "xy_exact_dirichlet",
        "final_time_requested": fmt(case.final_time),
        "return_code": str(completed.returncode),
        "subprocess_wall_time_sec": fmt(elapsed),
        "project_status": log_keys.get("project_amrex_euler_compare_status", ""),
        "bdltv20_paper_status": log_keys.get("bdltv20_paper_t1_s2_status", ""),
        "bdltv20_paper_boundary_policy": log_keys.get("bdltv20_paper_boundary_policy", ""),
        "field_boundary": log_keys.get("field_boundary", ""),
        "geometry_is_periodic": log_keys.get("geometry_is_periodic", ""),
        "completed_steps": log_keys.get("completed_steps", log_keys.get("bdltv20_paper_steps", "")),
        "completed_time": log_keys.get("completed_time", log_keys.get("bdltv20_paper_final_time", "")),
        "wall_time_sec": log_keys.get("wall_time_sec", fmt(elapsed)),
        "rho_min": log_keys.get("rho_min", log_keys.get("bdltv20_paper_rho_min", "")),
        "pressure_min": log_keys.get("pressure_min", log_keys.get("bdltv20_paper_pressure_min", "")),
        "bdltv20_paper_pressure_solver_failure_count": log_keys.get(
            "bdltv20_paper_pressure_solver_failure_count", ""
        ),
        "bdltv20_paper_pressure_relative_residual_linf_max": log_keys.get(
            "bdltv20_paper_pressure_relative_residual_linf_max", ""
        ),
        "field_csv": str(field_csv),
        "log": str(log_path),
        "command": str(command_path),
        "manifest": str(manifest_path),
    }
    ok = completed.returncode == 0 and field_csv.exists()
    print(f"{'OK' if ok else 'ERROR'} {case.key} {scheme.key} n={n}")
    if not ok:
        print(completed.stdout[-2000:])
    return row


def base_command(case: Case, n: int, field_csv: Path) -> list[str]:
    return [
        str(EXE),
        str(INPUTS),
        "amr.plot_int=-1",
        f"amr.n_cell={n} {Y_CELLS}",
        f"amr.max_grid_size={n}",
        "geometry.prob_lo=0.0 0.0",
        "geometry.prob_hi=1.0 0.1",
        "geometry.is_periodic=0 0",
        f"prob.gamma={GAMMA}",
        "euler.problem=toro1",
        f"prob.test={case.toro_test}",
        "prob.x0=0.5",
        "euler.field_boundary=exact_dirichlet",
        "euler.spatial_order=2",
        "euler.slope_limiter=minmod",
        f"stop_time={case.final_time:.17g}",
        f"euler.final_csv={field_csv}",
    ]


def imex_step_count(case: Case, n: int) -> int:
    """Scale accepted n=400 Riemann IMEX timesteps to the active grid."""
    reference_steps = {
        "sod_toro1": 194,
        "lax_rp2": 220,
        "toro3_strong_shock": 240,
        "toro4_strong_shock": 240,
    }[case.key]
    return max(1, int(round(reference_steps * n / 400.0)))


def build_metric_rows(run_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    cases_by_key = {case.key: case for case in CASES}
    for run in run_rows:
        if run["return_code"] != "0" or not Path(run["field_csv"]).exists():
            continue
        case = cases_by_key[run["case"]]
        x_values, numerical = read_midline_profile(Path(run["field_csv"]))
        exact = exact_profile(case, x_values)
        for quantity, values in numerical.items():
            errors = [abs(a - b) for a, b in zip(values, exact[quantity])]
            rows.append(
                {
                    **{key: run[key] for key in [
                        "case",
                        "case_label",
                        "toro_test",
                        "case_note",
                        "scheme",
                        "scheme_key",
                        "ncell",
                        "ny",
                        "boundary",
                        "final_time_requested",
                        "project_status",
                        "bdltv20_paper_status",
                        "bdltv20_paper_boundary_policy",
                        "field_boundary",
                        "geometry_is_periodic",
                        "completed_steps",
                        "completed_time",
                        "wall_time_sec",
                        "rho_min",
                        "pressure_min",
                        "bdltv20_paper_pressure_solver_failure_count",
                        "bdltv20_paper_pressure_relative_residual_linf_max",
                        "field_csv",
                        "log",
                        "command",
                    ]},
                    "quantity": quantity,
                    "l1_error": fmt(sum(errors) / len(errors)),
                    "linf_error": fmt(max(errors)),
                    "exact_solution_source": "scripts/riemann_exact.py",
                }
            )
    return rows


def build_order_rows(metric_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in metric_rows:
        grouped.setdefault((row["case"], row["scheme"], row["quantity"]), []).append(row)

    out: list[dict[str, str]] = []
    for (case_key, scheme, quantity), rows in sorted(grouped.items()):
        rows = sorted(rows, key=lambda row: int(row["ncell"]))
        if len(rows) < 2:
            continue
        log_n = [math.log(float(row["ncell"])) for row in rows if float(row["l1_error"]) > 0.0]
        log_e = [math.log(float(row["l1_error"])) for row in rows if float(row["l1_error"]) > 0.0]
        slope, intercept = least_squares(log_n, log_e)
        out.append(
            {
                "case": case_key,
                "case_label": rows[0]["case_label"],
                "scheme": scheme,
                "quantity": quantity,
                "grids": " ".join(row["ncell"] for row in rows),
                "l1_errors": " ".join(row["l1_error"] for row in rows),
                "least_squares_slope_log_error_vs_log_n": fmt(slope),
                "observed_l1_order_minus_slope": fmt(-slope),
                "fit_intercept": fmt(intercept),
                "claim_limit": "Riemann discontinuity refinement sanity check; not smooth formal order proof",
            }
        )
    return out


def read_midline_profile(path: Path) -> tuple[list[float], dict[str, list[float]]]:
    with path.open(newline="") as f:
        all_rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(f)]
    if not all_rows:
        raise SystemExit(f"No data rows in {path}")
    y_values = sorted({row["y"] for row in all_rows})
    y_mid = 0.5 * (y_values[0] + y_values[-1])
    selected_y = min(y_values, key=lambda y: abs(y - y_mid))
    rows = sorted((row for row in all_rows if row["y"] == selected_y), key=lambda row: row["x"])
    return [row["x"] for row in rows], {
        "rho": [row["rho"] for row in rows],
        "u": [row["u"] for row in rows],
        "p": [row["p"] for row in rows],
    }


def exact_profile(case: Case, x_values: list[float]) -> dict[str, list[float]]:
    star = solve_star_region(case.left, case.right, GAMMA)
    states = [sample(x, case.final_time, case.left, case.right, GAMMA, 0.5, star) for x in x_values]
    return {
        "rho": [state.rho for state in states],
        "u": [state.u for state in states],
        "p": [state.p for state in states],
    }


def make_plots(metric_rows: list[dict[str, str]]) -> None:
    for quantity, label in [("rho", "Density"), ("u", "Velocity"), ("p", "Pressure")]:
        plot_error(metric_rows, quantity, label, "ncell", "Grid cells in x", PACKAGE / "figures" / f"riemann_1d_{quantity}_l1_error_vs_grid.png")
        plot_error(metric_rows, quantity, label, "wall_time_sec", "Wall time (s)", PACKAGE / "figures" / f"riemann_1d_{quantity}_l1_error_vs_walltime.png")


def plot_error(rows: list[dict[str, str]], quantity: str, ylabel: str, x_key: str, xlabel: str, output: Path) -> None:
    case_keys = [case.key for case in CASES]
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.8), sharey=True, constrained_layout=True)
    axes_flat = axes.ravel()
    color_by_scheme = {scheme.label: scheme.color for scheme in SCHEMES}
    for ax, case_key in zip(axes_flat, case_keys):
        case_rows = [row for row in rows if row["case"] == case_key and row["quantity"] == quantity]
        for scheme in SCHEMES:
            data = sorted(
                [row for row in case_rows if row["scheme"] == scheme.label],
                key=lambda row: float(row[x_key]),
            )
            if not data:
                continue
            xs = [float(row[x_key]) for row in data]
            ys = [float(row["l1_error"]) for row in data]
            ax.scatter(xs, ys, s=32, color=color_by_scheme[scheme.label], label=scheme.label, zorder=3)
            positive = [(x, y) for x, y in zip(xs, ys) if x > 0.0 and y > 0.0 and math.isfinite(x) and math.isfinite(y)]
            if len(positive) >= 2:
                log_x = [math.log(item[0]) for item in positive]
                log_y = [math.log(item[1]) for item in positive]
                slope, intercept = least_squares(log_x, log_y)
                fit_x = [min(xs), max(xs)]
                fit_y = [math.exp(intercept + slope * math.log(x)) for x in fit_x]
                ax.plot(fit_x, fit_y, linewidth=1.4, color=color_by_scheme[scheme.label])
        label = next(case.label for case in CASES if case.key == case_key)
        ax.set_title(label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.22)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f"{ylabel} L1 error")
    axes_flat[0].legend(fontsize=8)
    fig.suptitle(f"AMReX 1D Riemann refinement: {ylabel} L1 error")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def copy_source_artifacts(run_rows: list[dict[str, str]]) -> None:
    for run in run_rows:
        for key, target_subdir in [
            ("field_csv", "fields"),
            ("log", "logs"),
            ("command", "commands"),
            ("manifest", "manifests"),
        ]:
            src = restore_canonical_artifact(Path(run[key]))
            dst = PACKAGE / "source_data" / target_subdir / src.name
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)


def restore_canonical_artifact(path: Path) -> Path:
    if path.exists():
        return path

    matching_paths = sorted(path.parent.glob(f"{path.stem} *{path.suffix}"))
    if len(matching_paths) == 1:
        matching_paths[0].replace(path)
        return path

    raise SystemExit(f"Missing required source: {path}")


def parse_log(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open(errors="replace") as f:
        for line in f:
            text = line.strip()
            if "=" not in text:
                continue
            key, value = text.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def write_readme(run_rows: list[dict[str, str]], order_rows: list[dict[str, str]]) -> None:
    passes = sum(row["return_code"] == "0" for row in run_rows)
    imex_rows = [row for row in run_rows if row["scheme"] == "IMEX T1/S2 BDLTV20"]
    imex_passes = sum(row["return_code"] == "0" and row["bdltv20_paper_status"] == "passed" for row in imex_rows)
    text = f"""# MPhil project AMReX 1D Riemann Accuracy/Convergence Package

This package contains a compact one-dimensional Riemann refinement sanity check for the
current AMReX schemes:

- `Explicit O2 HLLC`
- `Explicit O2 Low-Mach HLLC-P`
- `IMEX T1/S2 BDLTV20`

## Cases

- Sod / Toro 1, `t = 0.2`
- Lax / BDLTV20 RP2, `t = 0.14`
- Toro 3, `t = 0.012`
- Toro 4, `t = 0.035`

All rows use `geometry.is_periodic=0 0` and `euler.field_boundary=exact_dirichlet`.
The IMEX rows use `euler.bdltv20_paper_t1_s2=toro_xy_exact_dirichlet_2d`,
`spatial_order=2`, `slope_limiter=minmod`, and four Picard iterations.

## Figures

- `figures/riemann_1d_rho_l1_error_vs_grid.png`
- `figures/riemann_1d_u_l1_error_vs_grid.png`
- `figures/riemann_1d_p_l1_error_vs_grid.png`
- `figures/riemann_1d_rho_l1_error_vs_walltime.png`
- `figures/riemann_1d_u_l1_error_vs_walltime.png`
- `figures/riemann_1d_p_l1_error_vs_walltime.png`

The plotted lines are least-square log-log fits, not dot-to-dot connecting lines.

## Status

- Run rows completed: {passes}/{len(run_rows)}.
- IMEX rows with `bdltv20_paper_t1_s2_status=passed`: {imex_passes}/{len(imex_rows)}.
- Observed-order rows: {len(order_rows)}.

## Claim Limits

This is project preliminary validation evidence: exact-Riemann comparison and grid
refinement sanity checks for discontinuous problems. Since shocks and contacts dominate
the global error, these plots should not be presented as smooth formal second-order
accuracy proofs.

This package does not claim shock-bubble readiness, AP proof, production readiness,
full/3D BDLTV20, second-order time IMEX, AMR/MLMG/MPI/CUDA/GPU readiness, or CPU/GPU
efficiency.
"""
    (PACKAGE / "README.md").write_text(text)


def write_readiness_check(
    run_rows: list[dict[str, str]],
    metric_rows: list[dict[str, str]],
    order_rows: list[dict[str, str]],
) -> None:
    failed = [row for row in run_rows if row["return_code"] != "0"]
    non_exact = [
        row
        for row in run_rows
        if row.get("geometry_is_periodic") not in ("0 0", "")
        or row.get("field_boundary") not in ("exact_dirichlet", "")
    ]
    finite_metrics = sum(
        math.isfinite(float(row["l1_error"])) and math.isfinite(float(row["linf_error"]))
        for row in metric_rows
    )
    text = f"""# Riemann 1D Convergence Readiness Check

## Checks

- Failed executable rows: {len(failed)}.
- Rows with non-exact/periodic boundary evidence: {len(non_exact)}.
- Finite metric rows: {finite_metrics}/{len(metric_rows)}.
- Observed-order fit rows: {len(order_rows)}.

## Interpretation

Use this as a compact project validation table/figure set for 1D Riemann tests.
For discontinuous tests, error reduction with refinement is the meaningful sanity check;
do not describe the fitted slopes as clean formal second-order convergence.
"""
    (PACKAGE / "REPORT_READINESS_CHECK.md").write_text(text)


def write_manifest() -> None:
    rows: list[dict[str, str]] = []
    for path in sorted(PACKAGE.rglob("*")):
        if path.is_file():
            rows.append(
                {
                    "path": str(path),
                    "role": classify_path(path),
                    "source": "run_riemann_1d_convergence.py",
                }
            )
    write_csv(PACKAGE / "MANIFEST.csv", rows)


def classify_path(path: Path) -> str:
    if path.suffix == ".png":
        return "report figure"
    if path.name.endswith(".csv"):
        return "source data"
    if path.suffix == ".log":
        return "run log"
    if path.name.endswith(".command.txt"):
        return "reproduction command"
    return "documentation"


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise SystemExit(f"No rows to write for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def least_squares(xs: list[float], ys: list[float]) -> tuple[float, float]:
    if len(xs) != len(ys) or len(xs) < 2:
        return 0.0, ys[0] if ys else 0.0
    xbar = sum(xs) / len(xs)
    ybar = sum(ys) / len(ys)
    denom = sum((x - xbar) ** 2 for x in xs)
    if denom == 0.0:
        return 0.0, ybar
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / denom
    intercept = ybar - slope * xbar
    return slope, intercept


def shell_quote(text: str) -> str:
    if all(char.isalnum() or char in "/._=-" for char in text):
        return text
    return '"' + text.replace('"', '\\"') + '"'


def fmt(value: float) -> str:
    return f"{value:.12g}"


if __name__ == "__main__":
    main()
