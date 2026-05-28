#!/usr/bin/env python3
"""Run a focused low-Mach AP-behaviour probe.

The probe is deliberately small and claim-scoped. It writes outputs outside the
source tree by default:

    <parent>/runs/candidate/imex_ap_probe

It is evidence for AP-like behaviour, not a mathematical proof.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import time
from datetime import datetime
from pathlib import Path

from run_manifest import environment_build_flags, shell_join, utc_now, write_manifest


ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex"
INPUTS = ROOT / "amrex/apps/euler_compare/inputs-ci"
DEFAULT_OUT = ROOT.parent / "runs/candidate/imex_ap_probe"
PROBE_VERSION = "2026-05-28.v2"
DEFAULT_STEP_RATIO_THRESHOLD = 1.5
DEFAULT_PRESSURE_RATIO_THRESHOLD = 100.0


def pressure_from_mach(mach: float, gamma: float, density: float) -> float:
    return density / (gamma * mach * mach)


def split_floats(text: str) -> list[float]:
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def split_strings(text: str) -> list[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def parse_kv(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key:
            out[key] = value
    return out


def finite_float(value: str | None) -> float | None:
    if value in (None, ""):
        return None
    try:
        out = float(value)
    except ValueError:
        return None
    return out if math.isfinite(out) else None


def profile_metrics(path: Path, mach: float, gamma: float, density: float) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return {}

    def get(row: dict[str, str], *keys: str) -> float:
        for key in keys:
            value = row.get(key)
            if value not in (None, ""):
                return float(value)
        raise KeyError(keys)

    parsed: list[dict[str, float]] = []
    for row in rows:
        rho = get(row, "rho")
        u = get(row, "u")
        v = get(row, "v")
        p = get(row, "p", "pressure")
        exact_rho = get(row, "exact_rho")
        exact_u = get(row, "exact_u")
        exact_v = get(row, "exact_v")
        exact_p = get(row, "exact_p", "exact_pressure")
        rho_error = (
            get(row, "rho_error", "density_error")
            if row.get("rho_error") not in (None, "")
            or row.get("density_error") not in (None, "")
            else rho - exact_rho
        )
        velocity_error = (
            get(row, "velocity_error")
            if row.get("velocity_error") not in (None, "")
            else math.hypot(u - exact_u, v - exact_v)
        )
        pressure_error = (
            get(row, "p_error", "pressure_error")
            if row.get("p_error") not in (None, "")
            or row.get("pressure_error") not in (None, "")
            else p - exact_p
        )
        parsed.append(
            {
                "rho": rho,
                "u": u,
                "v": v,
                "p": p,
                "exact_rho": exact_rho,
                "exact_u": exact_u,
                "exact_v": exact_v,
                "exact_p": exact_p,
                "rho_error": rho_error,
                "velocity_error": velocity_error,
                "pressure_error": pressure_error,
            }
        )

    pressure_background = pressure_from_mach(mach, gamma, density)
    p_pert_errors = [
        abs((row["p"] - pressure_background) - (row["exact_p"] - pressure_background))
        for row in parsed
    ]
    exact_pert_scale = (
        sum(abs(row["exact_p"] - pressure_background) for row in parsed) / len(parsed)
    )
    exact_ke = sum(
        0.5 * row["exact_rho"] * (row["exact_u"] ** 2 + row["exact_v"] ** 2)
        for row in parsed
    )
    numerical_ke = sum(0.5 * row["rho"] * (row["u"] ** 2 + row["v"] ** 2) for row in parsed)
    return {
        "density_l1_error": f"{sum(abs(row['rho_error']) for row in parsed) / len(parsed):.12e}",
        "velocity_l1_error": f"{sum(abs(row['velocity_error']) for row in parsed) / len(parsed):.12e}",
        "pressure_l1_error": f"{sum(abs(row['pressure_error']) for row in parsed) / len(parsed):.12e}",
        "pressure_perturbation_l1_error": f"{sum(p_pert_errors) / len(p_pert_errors):.12e}",
        "pressure_perturbation_l1_relative_error": (
            f"{(sum(p_pert_errors) / len(p_pert_errors)) / exact_pert_scale:.12e}"
            if exact_pert_scale > 0.0
            else "nan"
        ),
        "kinetic_energy_ratio_from_csv": f"{numerical_ke / exact_ke:.12e}" if exact_ke > 0.0 else "nan",
    }


def common_args(args: argparse.Namespace, mach: float, final_csv: Path) -> list[str]:
    pressure = pressure_from_mach(mach, args.gamma, args.density)
    return [
        str(APP),
        "inputs-ci",
        "euler.problem=gresho_vortex",
        "euler.field_boundary=exact_dirichlet",
        f"euler.mach={mach:.17g}",
        f"euler.pressure={pressure:.17g}",
        f"euler.density_outer={args.density:.17g}",
        "euler.vortex_center=0.0 0.0",
        f"euler.gamma={args.gamma:.17g}",
        f"euler.final_csv={final_csv}",
        f"amr.n_cell={args.grid} {args.grid}",
        f"amr.max_grid_size={args.grid}",
        "amr.plot_int=-1",
        "geometry.prob_lo=-0.5 -0.5",
        "geometry.prob_hi=0.5 0.5",
        "geometry.is_periodic=0 0",
    ]


def command_for(
    case: str,
    mach: float,
    args: argparse.Namespace,
    final_csv: Path,
    direct_epsilon: float,
) -> list[str]:
    base = common_args(args, mach, final_csv)
    if case == "imex_source_map":
        return base + [
            "euler.method=imex",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            "euler.imex_form=bdltv20_t1_s2_source_map_picard",
            f"euler.imex_picard_iterations={args.imex_picard_iterations}",
            "euler.imex_acoustic_startup=0",
            "euler.imex_acoustic_cfl_cap=0.0",
            "euler.imex_pressure_stabilization=off",
            "euler.imex_predictor_dissipation=material",
            f"euler.imex_solver_tol={args.imex_solver_tol:.17g}",
            f"euler.imex_solver_max_iter={args.imex_solver_max_iter}",
            "euler.bdltv20_paper_t1_s2=off",
            f"euler.imex_cfl={args.imex_cfl:.17g}",
            f"stop_time={args.target_time:.17g}",
            f"max_step={args.imex_max_step}",
        ]
    if case == "imex_direct_bdltv20":
        return base + [
            "euler.method=imex",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            "euler.imex_form=bdltv20_t1_s2_source_map_picard",
            "euler.bdltv20_paper_t1_s2=gresho_exact_dirichlet_2d",
            f"euler.bdltv20_paper_epsilon={direct_epsilon:.17g}",
            f"euler.imex_picard_iterations={args.imex_picard_iterations}",
            "euler.imex_acoustic_startup=0",
            "euler.imex_acoustic_cfl_cap=0.0",
            "euler.imex_pressure_stabilization=off",
            "euler.imex_predictor_dissipation=material",
            f"euler.imex_solver_tol={args.imex_solver_tol:.17g}",
            f"euler.imex_solver_max_iter={args.imex_solver_max_iter}",
            "euler.bdltv20_paper_pressure_solver=gmres",
            "euler.bdltv20_paper_t1_s2_dt=-1",
            f"euler.bdltv20_paper_t1_s2_max_steps={args.imex_max_step}",
            f"euler.imex_cfl={args.imex_cfl:.17g}",
            f"stop_time={args.target_time:.17g}",
        ]
    if case == "explicit_hllc":
        return base + [
            "euler.method=explicit",
            "euler.riemann=hllc",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            f"euler.cfl={args.explicit_cfl:.17g}",
            f"stop_time={args.target_time:.17g}",
            f"max_step={args.explicit_max_step}",
        ]
    if case == "explicit_xie_am_hllc_p":
        return base + [
            "euler.method=explicit",
            "euler.riemann=xie_am_hllc_p",
            "euler.spatial_order=2",
            "euler.slope_limiter=minmod",
            f"euler.cfl={args.explicit_cfl:.17g}",
            f"stop_time={args.target_time:.17g}",
            f"max_step={args.explicit_max_step}",
        ]
    raise ValueError(case)


def run_row(command: list[str], log_path: Path, timeout: float) -> tuple[int | str, str, float, str, str]:
    start_utc = utc_now()
    start = time.perf_counter()
    try:
        proc = subprocess.run(
            command,
            cwd=str(APP.parent),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        wall = time.perf_counter() - start
        end_utc = utc_now()
        text = (exc.stdout or "") + "\n" + (exc.stderr or "")
        log_path.write_text(
            "$ " + shell_join(command) + f"\n# timed_out_after_sec={timeout:.12g}\n\n" + text
        )
        return "timeout", text, wall, start_utc, end_utc
    wall = time.perf_counter() - start
    end_utc = utc_now()
    log_path.write_text(
        "$ " + shell_join(command) + f"\n# measured_driver_wall_time_sec={wall:.12g}\n\n" + proc.stdout
    )
    return proc.returncode, proc.stdout, wall, start_utc, end_utc


def verdict_for_case(rows: list[dict[str, str]], args: argparse.Namespace) -> tuple[str, str]:
    ok_rows = [row for row in rows if row["status"] == "ok"]
    if len(ok_rows) != len(rows):
        return "failed_rows", "At least one row did not complete successfully."
    steps = [finite_float(row.get("completed_steps")) for row in ok_rows]
    steps = [value for value in steps if value is not None and value > 0.0]
    step_ratio = max(steps) / min(steps) if steps else math.inf
    pressure_errors = [
        finite_float(row.get("pressure_perturbation_l1_relative_error")) for row in ok_rows
    ]
    pressure_errors = [value for value in pressure_errors if value is not None and value > 0.0]
    pressure_ratio = max(pressure_errors) / min(pressure_errors) if pressure_errors else math.inf

    case = ok_rows[0]["case"]
    if case.startswith("explicit_"):
        return (
            "not_ap_by_acoustic_cfl",
            f"Explicit acoustic-CFL baseline: step ratio across Mach sweep is {step_ratio:.3g}.",
        )
    if step_ratio > args.verdict_step_ratio_threshold:
        return (
            "fails_material_timestep_probe",
            f"IMEX step count changes by factor {step_ratio:.3g}; expected near-constant material-CFL count.",
        )
    if pressure_ratio > args.verdict_pressure_ratio_threshold:
        return (
            "material_timestep_only_pressure_error_grows",
            f"Step count is Mach-independent, but pressure-perturbation relative error changes by factor {pressure_ratio:.3g}.",
        )
    return (
        "passes_small_ap_behaviour_probe",
        f"Step count ratio {step_ratio:.3g}; pressure-error ratio {pressure_ratio:.3g}.",
    )


def write_interpretation(rows: list[dict[str, str]], output_dir: Path, args: argparse.Namespace) -> None:
    by_case: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        group = row["case"]
        if row.get("direct_epsilon"):
            group = f"{row['case']} eps={row['direct_epsilon']}"
        by_case.setdefault(group, []).append(row)

    lines = [
        "# IMEX AP Behaviour Probe",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Probe version: {args.probe_version}",
        f"Grid: {args.grid} x {args.grid}",
        f"Target time: {args.target_time:.12g}",
        f"Mach values: {args.machs}",
        f"Direct epsilon values: {args.direct_epsilons or args.direct_epsilon}",
        f"IMEX Picard iterations: {args.imex_picard_iterations}",
        f"IMEX solver tolerance: {args.imex_solver_tol:.12g}",
        f"CFL values: explicit {args.explicit_cfl:.12g}, IMEX {args.imex_cfl:.12g}",
        (
            "Verdict thresholds: step ratio "
            f"{args.verdict_step_ratio_threshold:.12g}, pressure-error ratio "
            f"{args.verdict_pressure_ratio_threshold:.12g}"
        ),
        "",
        "This is numerical AP-like evidence only; it is not a mathematical proof.",
        "",
        "## Verdicts",
        "",
        "| Case | Verdict | Reason |",
        "|---|---|---|",
    ]
    for case in sorted(by_case):
        verdict, reason = verdict_for_case(by_case[case], args)
        lines.append(f"| `{case}` | `{verdict}` | {reason} |")

    lines.extend(
        [
            "",
            "## Key Rows",
            "",
            "| Case | Mach | Direct epsilon | Status | Steps | Sound speed | Material rate | Acoustic rate | Pressure pert. rel. L1 |",
            "|---|---:|---:|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        lines.append(
            "| {case} | {mach} | {eps} | {status} | {steps} | {sound} | {material} | {acoustic} | {perr} |".format(
                case=f"`{row['case']}`",
                mach=row["mach"],
                eps=row.get("direct_epsilon", ""),
                status=row["status"],
                steps=row.get("completed_steps", ""),
                sound=row.get("reference_sound_speed", ""),
                material=row.get("imex_timestep_material_rate_max", ""),
                acoustic=row.get("imex_timestep_acoustic_rate_max", ""),
                perr=row.get("pressure_perturbation_l1_relative_error", ""),
            )
        )

    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "- Explicit HLLC and Low-Mach HLLC-P are acoustic-CFL baselines, not AP schemes.",
            "- The source-map IMEX route can show Mach-independent timesteps, but the code labels it as dimensional epsilon=1 with no formal AP claim.",
            "- The direct BDLTV20 route is the only route here with an epsilon-scaled AP derivation boundary.",
            "- Check `bdltv20_paper_epsilon` in the summary before using this as evidence for a literal epsilon-limit claim.",
            "- A report AP theorem still needs the written assumptions and discrete limit argument; this probe only checks behaviour.",
            "- Very short runs are smoke evidence only; report-facing AP evidence needs longer horizon and multiple grids.",
            "",
            f"Summary CSV: `{output_dir / 'summary.csv'}`",
        ]
    )
    (output_dir / "interpretation.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--machs", default="0.1,0.01,0.001")
    parser.add_argument(
        "--cases",
        default="imex_source_map,imex_direct_bdltv20,explicit_hllc,explicit_xie_am_hllc_p",
    )
    parser.add_argument("--grid", type=int, default=16)
    parser.add_argument("--target-time", type=float, default=0.05)
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--imex-cfl", type=float, default=0.8)
    parser.add_argument("--explicit-cfl", type=float, default=0.4)
    parser.add_argument("--imex-picard-iterations", type=int, default=4)
    parser.add_argument("--imex-solver-tol", type=float, default=1.0e-8)
    parser.add_argument("--imex-solver-max-iter", type=int, default=1000)
    parser.add_argument("--imex-max-step", type=int, default=10000)
    parser.add_argument("--explicit-max-step", type=int, default=200000)
    parser.add_argument("--row-timeout-sec", type=float, default=180.0)
    parser.add_argument("--direct-epsilon", type=float, default=1.0)
    parser.add_argument(
        "--direct-epsilons",
        default="",
        help="Optional comma-separated epsilon sweep for imex_direct_bdltv20 rows.",
    )
    parser.add_argument("--probe-version", default=PROBE_VERSION)
    parser.add_argument(
        "--verdict-step-ratio-threshold",
        type=float,
        default=DEFAULT_STEP_RATIO_THRESHOLD,
    )
    parser.add_argument(
        "--verdict-pressure-ratio-threshold",
        type=float,
        default=DEFAULT_PRESSURE_RATIO_THRESHOLD,
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    logs_dir = output_dir / "logs"
    fields_dir = output_dir / "fields"
    commands_dir = output_dir / "commands"
    manifests_dir = output_dir / "manifests"
    for path in (logs_dir, fields_dir, commands_dir, manifests_dir):
        path.mkdir(parents=True, exist_ok=True)

    cases = split_strings(args.cases)
    direct_epsilon_values = (
        split_floats(args.direct_epsilons) if args.direct_epsilons else [args.direct_epsilon]
    )
    valid_cases = {
        "imex_source_map",
        "imex_direct_bdltv20",
        "explicit_hllc",
        "explicit_xie_am_hllc_p",
    }
    unknown = [case for case in cases if case not in valid_cases]
    if unknown:
        raise SystemExit(f"unknown cases: {','.join(unknown)}")

    rows: list[dict[str, str]] = []
    for mach in split_floats(args.machs):
        for case in cases:
            eps_values = direct_epsilon_values if case == "imex_direct_bdltv20" else [args.direct_epsilon]
            for direct_epsilon in eps_values:
                eps_part = f"_eps{direct_epsilon:g}" if case == "imex_direct_bdltv20" else ""
                row_id = f"ap_probe_{case}{eps_part}_m{mach:g}_n{args.grid}".replace(".", "p")
                final_csv = fields_dir / f"{row_id}_final.csv"
                log_path = logs_dir / f"{row_id}.log"
                command_path = commands_dir / f"{row_id}.txt"
                manifest_path = manifests_dir / f"{row_id}.manifest.json"
                command = command_for(case, mach, args, final_csv, direct_epsilon)
                command_path.write_text(shell_join(command) + "\n")
                row: dict[str, str] = {
                    "row_id": row_id,
                    "probe_version": args.probe_version,
                    "case": case,
                    "mach": f"{mach:.12g}",
                    "direct_epsilon": f"{direct_epsilon:.12g}" if case == "imex_direct_bdltv20" else "",
                    "n": str(args.grid),
                    "target_time": f"{args.target_time:.12g}",
                    "imex_picard_iterations": str(args.imex_picard_iterations),
                    "imex_solver_tol": f"{args.imex_solver_tol:.12g}",
                    "imex_cfl": f"{args.imex_cfl:.12g}",
                    "explicit_cfl": f"{args.explicit_cfl:.12g}",
                    "verdict_step_ratio_threshold": f"{args.verdict_step_ratio_threshold:.12g}",
                    "verdict_pressure_ratio_threshold": f"{args.verdict_pressure_ratio_threshold:.12g}",
                    "final_csv": str(final_csv),
                    "run_log": str(log_path),
                    "command_file": str(command_path),
                    "manifest": str(manifest_path),
                    "status": "dry_run" if args.dry_run else "not_run",
                }
                if args.dry_run:
                    rows.append(row)
                    continue
                print(f"[run] {row_id}", flush=True)
                returncode, stdout, wall, start_utc, end_utc = run_row(
                    command, log_path, args.row_timeout_sec
                )
                kv = parse_kv(stdout)
                write_manifest(
                    manifest_path,
                    root=ROOT,
                    row_id=row_id,
                    command=command,
                    start_utc=start_utc,
                    end_utc=end_utc,
                    wall_time_s=wall,
                    exit_code=returncode,
                    output_root=output_dir,
                    output_class="candidate",
                    input_files=[INPUTS],
                    output_files=[final_csv, log_path, command_path],
                    build_flags=environment_build_flags({"DIM": "2"}),
                    notes=(
                        f"Focused low-Mach AP-behaviour probe {args.probe_version}; "
                        "not a formal proof."
                    ),
                )
                row.update(
                    {
                        "returncode": str(returncode),
                        "driver_wall_time_sec": f"{wall:.12g}",
                        "app_wall_time_sec": kv.get("wall_time_sec", ""),
                        "completed_steps": kv.get("completed_steps", kv.get("bdltv20_paper_steps", "")),
                        "final_time": kv.get("completed_time", kv.get("bdltv20_paper_final_time", "")),
                    "project_status": kv.get("project_amrex_euler_compare_status", ""),
                    "imex_route_tag": kv.get("imex_route_tag", ""),
                    "bdltv20_paper_status": kv.get("bdltv20_paper_t1_s2_status", ""),
                        "reference_sound_speed": kv.get("reference_sound_speed", ""),
                        "velocity_scale": kv.get("velocity_scale", ""),
                        "imex_timestep_material_rate_max": kv.get("imex_timestep_material_rate_max", ""),
                        "imex_timestep_acoustic_rate_max": kv.get("imex_timestep_acoustic_rate_max", ""),
                        "imex_acoustic_startup_steps": kv.get("imex_acoustic_startup_steps", ""),
                        "imex_acoustic_cap_steps": kv.get("imex_acoustic_cap_steps", ""),
                        "imex_solver_status": kv.get("imex_solver_status", ""),
                        "imex_solver_iterations_max": kv.get("imex_solver_iterations_max", ""),
                        "bdltv20_paper_epsilon": kv.get("bdltv20_paper_epsilon", ""),
                        "bdltv20_paper_timestep_policy": kv.get("bdltv20_paper_timestep_policy", ""),
                        "bdltv20_paper_gmres_iterations_max": kv.get(
                            "bdltv20_paper_gmres_iterations_max", ""
                        ),
                        "bdltv20_paper_pressure_relative_residual_linf_max": kv.get(
                            "bdltv20_paper_pressure_relative_residual_linf_max", ""
                        ),
                        "rho_min": kv.get("rho_min", kv.get("bdltv20_paper_rho_min", "")),
                        "pressure_min": kv.get("pressure_min", kv.get("bdltv20_paper_pressure_min", "")),
                        "nonfinite_count": kv.get("nonfinite_count", kv.get("bdltv20_paper_nonfinite_count", "")),
                    }
                )
                row.update(profile_metrics(final_csv, mach, args.gamma, args.density))
                ok_project = row["project_status"] in {"", "ok"}
                ok_direct = row["bdltv20_paper_status"] in {"", "passed"}
                row["status"] = (
                    "ok"
                    if returncode == 0 and final_csv.exists() and ok_project and ok_direct
                    else "failed"
                )
                rows.append(row)

    summary_path = output_dir / "summary.csv"
    fieldnames = sorted({key for row in rows for key in row})
    preferred = [
        "row_id",
        "probe_version",
        "status",
        "case",
        "imex_route_tag",
        "mach",
        "direct_epsilon",
        "n",
        "target_time",
        "imex_picard_iterations",
        "imex_solver_tol",
        "imex_cfl",
        "explicit_cfl",
        "completed_steps",
        "final_time",
        "reference_sound_speed",
        "velocity_scale",
        "imex_timestep_material_rate_max",
        "imex_timestep_acoustic_rate_max",
        "bdltv20_paper_epsilon",
        "bdltv20_paper_timestep_policy",
        "pressure_perturbation_l1_relative_error",
        "density_l1_error",
        "velocity_l1_error",
        "driver_wall_time_sec",
        "run_log",
        "manifest",
    ]
    ordered = preferred + [key for key in fieldnames if key not in preferred]
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered)
        writer.writeheader()
        writer.writerows(rows)
    if not args.dry_run:
        write_interpretation(rows, output_dir, args)
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
