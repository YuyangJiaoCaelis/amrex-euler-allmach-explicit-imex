#!/usr/bin/env python3
"""Run candidate MPI rank-scan timing rows for retained explicit schemes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import statistics
from pathlib import Path

from run_mpi_explicit_agreement import (
    explicit_cases,
    failure_from_log,
    parse_key_values,
    parse_float,
    status_from_log,
    with_build_env,
)


def split_csv_ints(text: str) -> list[int]:
    values = [int(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise ValueError("expected at least one integer rank")
    if any(value <= 0 for value in values):
        raise ValueError("MPI ranks must be positive")
    return values


def selected_cases(case_set: str, names: str) -> list[object]:
    cases = explicit_cases(case_set)
    if not names:
        return cases
    requested = {name.strip() for name in names.split(",") if name.strip()}
    available = {case.row_id for case in cases}
    unknown = sorted(requested - available)
    if unknown:
        raise SystemExit(
            "Unknown --cases entries: "
            + ",".join(unknown)
            + "; valid entries: "
            + ",".join(sorted(available))
        )
    return [case for case in cases if case.row_id in requested]


def finite_float(text: str) -> float | None:
    value = parse_float(text)
    if value is None or not math.isfinite(value):
        return None
    return value


def median_text(values: list[float]) -> str:
    if not values:
        return ""
    return f"{statistics.median(values):.12g}"


def mean_text(values: list[float]) -> str:
    if not values:
        return ""
    return f"{statistics.mean(values):.12g}"


def best_text(values: list[float]) -> str:
    if not values:
        return ""
    return f"{min(values):.12g}"


def run_rank_row(
    *,
    root: Path,
    output: Path,
    mpi_exe: Path,
    mpirun: str,
    case: object,
    case_set: str,
    ranks: int,
    repeat: int,
    output_class: str,
    row_timeout_sec: float,
    write_final_csv: bool,
) -> dict[str, str]:
    logs_dir = output / "logs"
    commands_dir = output / "commands"
    manifests_dir = output / "manifests"
    fields_dir = output / "fields"
    for directory in (logs_dir, commands_dir, manifests_dir, fields_dir):
        directory.mkdir(parents=True, exist_ok=True)

    row_id = f"{case.row_id}_mpi{ranks}_rep{repeat:02d}"
    log_path = logs_dir / f"{row_id}.log"
    command_path = commands_dir / f"{row_id}.txt"
    manifest_path = manifests_dir / f"{row_id}.manifest.json"
    final_csv = fields_dir / f"{row_id}.csv"
    output_files = [final_csv] if write_final_csv else []
    command = [
        mpirun,
        "-np",
        str(ranks),
        str(mpi_exe),
        *case.args,
    ]
    if write_final_csv:
        command.append(f"euler.final_csv={final_csv}")

    returncode = with_build_env(
        True,
        command,
        root=root,
        row_id=row_id,
        output_dir=output,
        log_path=log_path,
        command_path=command_path,
        manifest_path=manifest_path,
        output_files=output_files,
        output_globs=[],
        output_class=output_class,
        notes=(
            f"{case_set} explicit MPI rank-scan timing row; "
            "per-cell final CSV output disabled unless --write-final-csv is set."
        ),
        cwd=root,
        timeout_sec=row_timeout_sec,
    )

    values = parse_key_values(log_path)
    manifest = json.loads(manifest_path.read_text())
    status = status_from_log(values)
    failure = failure_from_log(values, status)
    final_csv_written = final_csv.exists() if write_final_csv else False
    ok = (
        returncode == 0
        and status in {"ok", "passed"}
        and failure in {"ok", "none"}
        and (not write_final_csv or final_csv_written)
    )

    return {
        "row_id": row_id,
        "case_id": case.row_id,
        "case_set": case_set,
        "output_class": output_class,
        "mpi_ranks": str(ranks),
        "repeat": str(repeat),
        "returncode": str(returncode),
        "status": status,
        "failure_category": failure,
        "passed": "yes" if ok else "no",
        "driver_wall_time_sec": f"{float(manifest['wall_time_s']):.12g}",
        "app_wall_time_sec": values.get("wall_time_sec", ""),
        "completed_steps": values.get("completed_steps", ""),
        "completed_time": values.get("completed_time", ""),
        "rho_min": values.get("rho_min", ""),
        "pressure_min": values.get("pressure_min", ""),
        "nonfinite_count": values.get("nonfinite_count", ""),
        "initial_mass": values.get("initial_mass", ""),
        "final_mass": values.get("final_mass", ""),
        "mass_drift": values.get("mass_drift", ""),
        "initial_energy": values.get("initial_energy", ""),
        "final_energy": values.get("final_energy", ""),
        "energy_drift": values.get("energy_drift", ""),
        "write_final_csv": "yes" if write_final_csv else "no",
        "final_csv": str(final_csv) if write_final_csv else "",
        "log": str(log_path),
        "command": str(command_path),
        "manifest": str(manifest_path),
    }


def aggregate_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_case_rank: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        by_case_rank.setdefault((row["case_id"], row["mpi_ranks"]), []).append(row)

    rank1_medians: dict[str, float] = {}
    aggregates: list[dict[str, str]] = []
    for (case_id, ranks), group in sorted(by_case_rank.items()):
        passed = [row for row in group if row["passed"] == "yes"]
        driver_times = [
            value
            for row in passed
            if (value := finite_float(row["driver_wall_time_sec"])) is not None
        ]
        app_times = [
            value
            for row in passed
            if (value := finite_float(row["app_wall_time_sec"])) is not None
        ]
        driver_median = statistics.median(driver_times) if driver_times else math.nan
        if ranks == "1" and math.isfinite(driver_median):
            rank1_medians[case_id] = driver_median
        speedup = ""
        if math.isfinite(driver_median) and rank1_medians.get(case_id):
            speedup = f"{rank1_medians[case_id] / driver_median:.12g}"

        aggregates.append(
            {
                "case_id": case_id,
                "mpi_ranks": ranks,
                "repeats": str(len(group)),
                "passed_repeats": str(len(passed)),
                "driver_wall_time_median_sec": (
                    f"{driver_median:.12g}" if math.isfinite(driver_median) else ""
                ),
                "driver_wall_time_best_sec": best_text(driver_times),
                "driver_wall_time_mean_sec": mean_text(driver_times),
                "app_wall_time_median_sec": median_text(app_times),
                "app_wall_time_best_sec": best_text(app_times),
                "app_wall_time_mean_sec": mean_text(app_times),
                "speedup_vs_mpi1_driver_median": speedup,
            }
        )
    return aggregates


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError(f"no rows to write: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--case-set", choices=["smoke", "candidate"], default="candidate")
    parser.add_argument("--cases", default="", help="Comma-separated case row_ids to run.")
    parser.add_argument("--ranks", default="1,2,4")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--row-timeout-sec", type=float, default=600.0)
    parser.add_argument(
        "--output-class",
        choices=["exploratory", "candidate", "frozen"],
        default="candidate",
    )
    parser.add_argument(
        "--write-final-csv",
        action="store_true",
        help="Write final CSVs during timing rows. Disabled by default to avoid measuring output cost.",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")

    root = args.root.resolve()
    output = (
        args.output.resolve()
        if args.output
        else root.parent / f"{root.name}_mpi_explicit_rank_scan"
    )
    if output.exists():
        if not args.force:
            raise SystemExit(f"Output directory already exists: {output}. Use --force to replace it.")
        shutil.rmtree(output)
    output.mkdir(parents=True)

    mpi_exe = root / "amrex/apps/euler_compare/amrex_euler_compare2d.gnu.MPI.ex"
    if not mpi_exe.exists():
        raise SystemExit(f"MPI executable not found: {mpi_exe}")
    mpirun = shutil.which(args.mpirun)
    if not mpirun:
        raise SystemExit(f"MPI launcher not found on PATH: {args.mpirun}")

    cases = selected_cases(args.case_set, args.cases)
    ranks_list = split_csv_ints(args.ranks)
    rows: list[dict[str, str]] = []
    for case in cases:
        for ranks in ranks_list:
            for repeat in range(1, args.repeats + 1):
                row = run_rank_row(
                    root=root,
                    output=output,
                    mpi_exe=mpi_exe,
                    mpirun=mpirun,
                    case=case,
                    case_set=args.case_set,
                    ranks=ranks,
                    repeat=repeat,
                    output_class=args.output_class,
                    row_timeout_sec=args.row_timeout_sec,
                    write_final_csv=args.write_final_csv,
                )
                rows.append(row)
                print(
                    row["case_id"],
                    f"np={ranks}",
                    f"rep={repeat}",
                    row["returncode"],
                    row["status"],
                    row["passed"],
                    row["driver_wall_time_sec"],
                    flush=True,
                )

    write_csv(output / "summary.csv", rows)
    aggregates = aggregate_rows(rows)
    write_csv(output / "rank_summary.csv", aggregates)
    print(f"summary {output / 'summary.csv'}")
    print(f"rank_summary {output / 'rank_summary.csv'}")

    failed = [row for row in rows if row["passed"] != "yes"]
    if failed:
        print("failed rows: " + ", ".join(row["row_id"] for row in failed))
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
