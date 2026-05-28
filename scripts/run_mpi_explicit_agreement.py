#!/usr/bin/env python3
"""Check serial and MPI agreement for retained explicit Report 2 smoke rows."""

from __future__ import annotations

import argparse
import csv
import math
import os
import signal
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from run_manifest import expand_output_paths, shell_join, utc_now, write_manifest


STATUS_KEYS = (
    "project_amrex_euler_compare_status",
    "bdltv20_paper_t1_s2_status",
)
IGNORED_COMPARE_COLUMNS = {
    "wall_time_sec",
    "snapshot_dir",
}


@dataclass(frozen=True)
class AgreementCase:
    row_id: str
    args: tuple[str, ...]
    snapshot_times: str = ""


def explicit_cases() -> list[AgreementCase]:
    base = (
        "max_step=100000",
        "amr.plot_int=-1",
        "amr.max_grid_size=16",
    )

    riemann_common = base + (
        "amr.n_cell=128 4",
        "geometry.prob_lo=0 0",
        "geometry.prob_hi=1 0.03125",
        "geometry.is_periodic=0 0",
        "euler.problem=toro1",
        "euler.method=explicit",
        "euler.spatial_order=2",
        "euler.field_boundary=exact_dirichlet",
        "stop_time=0.02",
        "euler.cfl=0.45",
    )
    gresho_common = base + (
        "amr.n_cell=32 32",
        "geometry.prob_lo=0 0",
        "geometry.prob_hi=1 1",
        "geometry.is_periodic=0 0",
        "euler.problem=gresho_vortex",
        "euler.method=explicit",
        "euler.spatial_order=2",
        "euler.field_boundary=exact_dirichlet",
        "euler.mach=0.01",
        "stop_time=0.05",
        "euler.cfl=0.45",
    )
    advection_common = base + (
        "amr.n_cell=32 32",
        "geometry.prob_lo=0 0",
        "geometry.prob_hi=1 1",
        "geometry.is_periodic=1 1",
        "euler.problem=advection_blob",
        "euler.method=explicit",
        "euler.spatial_order=2",
        "stop_time=0.05",
        "euler.cfl=0.45",
    )
    shock_common = base + (
        "amr.n_cell=64 16",
        "geometry.prob_lo=0 0",
        "geometry.prob_hi=2 0.5",
        "geometry.is_periodic=0 0",
        "euler.problem=shock_density_bubble_2d",
        "euler.method=explicit",
        "euler.spatial_order=2",
        "euler.shock_density_bubble_snapshot_times=0,0.01",
        "stop_time=0.01",
        "euler.cfl=0.45",
    )

    return [
        AgreementCase(
            "riemann_sod_explicit_hllc",
            riemann_common + ("euler.riemann=hllc",),
        ),
        AgreementCase(
            "riemann_sod_explicit_lowmach_hllcp",
            riemann_common + ("euler.riemann=xie_am_hllc_p",),
        ),
        AgreementCase(
            "gresho_m0p01_explicit_hllc_n32",
            gresho_common + ("euler.riemann=hllc",),
        ),
        AgreementCase(
            "gresho_m0p01_explicit_lowmach_hllcp_n32",
            gresho_common + ("euler.riemann=xie_am_hllc_p",),
        ),
        AgreementCase(
            "advection_blob_explicit_hllc_n32",
            advection_common + ("euler.riemann=hllc",),
        ),
        AgreementCase(
            "advection_blob_explicit_lowmach_hllcp_n32",
            advection_common + ("euler.riemann=xie_am_hllc_p",),
        ),
        AgreementCase(
            "shock_density_bubble_explicit_hllc_64x16",
            shock_common + ("euler.riemann=hllc",),
            snapshot_times="0,0.01",
        ),
        AgreementCase(
            "shock_density_bubble_explicit_lowmach_hllcp_64x16",
            shock_common + ("euler.riemann=xie_am_hllc_p",),
            snapshot_times="0,0.01",
        ),
    ]


def parse_key_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def status_from_log(values: dict[str, str]) -> str:
    for key in STATUS_KEYS:
        if values.get(key):
            return values[key]
    return ""


def failure_from_log(values: dict[str, str], status: str) -> str:
    return values.get("failure_category", "none" if status in {"ok", "passed"} else "")


def csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        return [], []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def parse_float(text: str) -> float | None:
    try:
        return float(text)
    except ValueError:
        return None


def row_sort_key(row: dict[str, str]) -> tuple[float | str, ...]:
    if "i" in row and "j" in row:
        return (parse_float(row["i"]) or 0.0, parse_float(row["j"]) or 0.0)
    if "x" in row and "y" in row:
        return (parse_float(row["x"]) or 0.0, parse_float(row["y"]) or 0.0)
    return tuple(row.get(key, "") for key in sorted(row))


def numeric_difference(left: str, right: str) -> float | None:
    left_float = parse_float(left)
    right_float = parse_float(right)
    if left_float is None or right_float is None:
        return None
    if math.isnan(left_float) and math.isnan(right_float):
        return 0.0
    if math.isnan(left_float) or math.isnan(right_float):
        return math.inf
    return abs(left_float - right_float)


def compare_csv(left_path: Path, right_path: Path) -> dict[str, str]:
    left_header, left_rows = csv_rows(left_path)
    right_header, right_rows = csv_rows(right_path)
    result = {
        "serial_rows": str(len(left_rows)),
        "mpi_rows": str(len(right_rows)),
        "max_abs_diff": "",
        "max_abs_diff_column": "",
        "text_mismatch_count": "0",
        "csv_compare_status": "ok",
    }

    if not left_path.exists() or not right_path.exists():
        result["csv_compare_status"] = "missing_csv"
        return result
    if left_header != right_header:
        result["csv_compare_status"] = "header_mismatch"
        return result
    if len(left_rows) != len(right_rows):
        result["csv_compare_status"] = "row_count_mismatch"
        return result

    max_diff = 0.0
    max_column = ""
    text_mismatches = 0
    columns = [column for column in left_header if column not in IGNORED_COMPARE_COLUMNS]
    for left, right in zip(
        sorted(left_rows, key=row_sort_key),
        sorted(right_rows, key=row_sort_key),
    ):
        for column in columns:
            diff = numeric_difference(left[column], right[column])
            if diff is None:
                if left[column] != right[column]:
                    text_mismatches += 1
                continue
            if diff > max_diff:
                max_diff = diff
                max_column = column

    result["max_abs_diff"] = f"{max_diff:.17g}"
    result["max_abs_diff_column"] = max_column
    result["text_mismatch_count"] = str(text_mismatches)
    if not math.isfinite(max_diff):
        result["csv_compare_status"] = "nonfinite_difference"
    elif text_mismatches:
        result["csv_compare_status"] = "text_mismatch"
    return result


def run_with_timeout_manifest(
    command: list[str],
    *,
    root: Path,
    row_id: str,
    output_dir: Path,
    log_path: Path,
    command_path: Path,
    manifest_path: Path,
    output_files: list[Path],
    output_globs: list[str],
    output_class: str,
    notes: str,
    cwd: Path,
    timeout_sec: float,
) -> int | str:
    command_path.parent.mkdir(parents=True, exist_ok=True)
    command_path.write_text(shell_join(command) + "\n")

    start_utc = utc_now()
    start = time.perf_counter()
    proc = subprocess.Popen(
        command,
        cwd=str(cwd),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    timed_out = False
    try:
        stdout, _ = proc.communicate(timeout=timeout_sec)
        exit_code: int | str = proc.returncode
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            stdout, _ = proc.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            stdout, _ = proc.communicate()
        exit_code = f"timeout_after_{timeout_sec:g}s"

    wall = time.perf_counter() - start
    end_utc = utc_now()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if timed_out:
        stdout = (stdout or "") + f"\ncommand_timeout_after_seconds={timeout_sec:g}\n"
    log_path.write_text(stdout or "")

    expanded_outputs = expand_output_paths(output_files + [log_path, command_path], output_globs)
    write_manifest(
        manifest_path,
        root=root,
        row_id=row_id,
        command=command,
        start_utc=start_utc,
        end_utc=end_utc,
        wall_time_s=wall,
        exit_code=exit_code,
        output_root=output_dir,
        output_class=output_class,
        output_files=expanded_outputs,
        notes=notes,
        extra={"timeout_sec": timeout_sec, "timed_out": timed_out},
    )
    return exit_code


def with_build_env(use_mpi: bool, command: list[str], **kwargs: object) -> int | str:
    replacements = {
        "DIM": "2",
        "USE_MPI": "TRUE" if use_mpi else "FALSE",
        "USE_OMP": "FALSE",
        "USE_CUDA": "FALSE",
    }
    old_values = {key: os.environ.get(key) for key in replacements}
    os.environ.update(replacements)
    try:
        return run_with_timeout_manifest(command, **kwargs)
    finally:
        for key, old_value in old_values.items():
            if old_value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = old_value


def run_case(
    *,
    case: AgreementCase,
    root: Path,
    output: Path,
    serial_exe: Path,
    mpi_exe: Path,
    mpirun: str,
    mpi_ranks: int,
    tolerance: float,
    row_timeout_sec: float,
) -> dict[str, str]:
    fields_dir = output / "fields"
    logs_dir = output / "logs"
    commands_dir = output / "commands"
    manifests_dir = output / "manifests"
    snapshots_dir = output / "snapshots"
    for directory in (fields_dir, logs_dir, commands_dir, manifests_dir, snapshots_dir):
        directory.mkdir(parents=True, exist_ok=True)

    serial_csv = fields_dir / f"{case.row_id}_serial.csv"
    mpi_csv = fields_dir / f"{case.row_id}_mpi{mpi_ranks}.csv"
    serial_log = logs_dir / f"{case.row_id}_serial.log"
    mpi_log = logs_dir / f"{case.row_id}_mpi{mpi_ranks}.log"
    serial_manifest = manifests_dir / f"{case.row_id}_serial.manifest.json"
    mpi_manifest = manifests_dir / f"{case.row_id}_mpi{mpi_ranks}.manifest.json"
    serial_snapshots = snapshots_dir / f"{case.row_id}_serial"
    mpi_snapshots = snapshots_dir / f"{case.row_id}_mpi{mpi_ranks}"

    serial_extra: list[str] = []
    mpi_extra: list[str] = []
    serial_globs: list[str] = []
    mpi_globs: list[str] = []
    if case.snapshot_times:
        serial_extra.append(f"euler.shock_density_bubble_snapshot_dir={serial_snapshots}")
        mpi_extra.append(f"euler.shock_density_bubble_snapshot_dir={mpi_snapshots}")
        serial_globs.append(str(serial_snapshots / "*.csv"))
        mpi_globs.append(str(mpi_snapshots / "*.csv"))

    serial_command = [
        str(serial_exe),
        *case.args,
        f"euler.final_csv={serial_csv}",
        *serial_extra,
    ]
    mpi_command = [
        mpirun,
        "-np",
        str(mpi_ranks),
        str(mpi_exe),
        *case.args,
        f"euler.final_csv={mpi_csv}",
        *mpi_extra,
    ]

    serial_returncode = with_build_env(
        False,
        serial_command,
        root=root,
        row_id=f"{case.row_id}_serial",
        output_dir=output,
        log_path=serial_log,
        command_path=commands_dir / f"{case.row_id}_serial.txt",
        manifest_path=serial_manifest,
        output_files=[serial_csv],
        output_globs=serial_globs,
        output_class="exploratory",
        notes="Serial side of explicit serial/MPI agreement smoke check.",
        cwd=root,
        timeout_sec=row_timeout_sec,
    )
    mpi_returncode = with_build_env(
        True,
        mpi_command,
        root=root,
        row_id=f"{case.row_id}_mpi{mpi_ranks}",
        output_dir=output,
        log_path=mpi_log,
        command_path=commands_dir / f"{case.row_id}_mpi{mpi_ranks}.txt",
        manifest_path=mpi_manifest,
        output_files=[mpi_csv],
        output_globs=mpi_globs,
        output_class="exploratory",
        notes=f"{mpi_ranks}-rank MPI side of explicit serial/MPI agreement smoke check.",
        cwd=root,
        timeout_sec=row_timeout_sec,
    )

    serial_values = parse_key_values(serial_log)
    mpi_values = parse_key_values(mpi_log)
    serial_status = status_from_log(serial_values)
    mpi_status = status_from_log(mpi_values)
    serial_failure = failure_from_log(serial_values, serial_status)
    mpi_failure = failure_from_log(mpi_values, mpi_status)
    comparison = compare_csv(serial_csv, mpi_csv)
    max_diff = parse_float(comparison["max_abs_diff"]) if comparison["max_abs_diff"] else math.inf
    within_tolerance = (
        serial_returncode == 0
        and mpi_returncode == 0
        and serial_status in {"ok", "passed"}
        and mpi_status in {"ok", "passed"}
        and serial_failure in {"ok", "none"}
        and mpi_failure in {"ok", "none"}
        and comparison["csv_compare_status"] == "ok"
        and max_diff is not None
        and max_diff <= tolerance
    )

    return {
        "row_id": case.row_id,
        "serial_returncode": str(serial_returncode),
        "mpi_returncode": str(mpi_returncode),
        "serial_status": serial_status,
        "mpi_status": mpi_status,
        "serial_failure_category": serial_failure,
        "mpi_failure_category": mpi_failure,
        **comparison,
        "tolerance": f"{tolerance:.17g}",
        "within_tolerance": "yes" if within_tolerance else "no",
        "serial_csv": str(serial_csv),
        "mpi_csv": str(mpi_csv),
        "serial_log": str(serial_log),
        "mpi_log": str(mpi_log),
        "serial_manifest": str(serial_manifest),
        "mpi_manifest": str(mpi_manifest),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output directory. Defaults beside the source tree, not inside it.",
    )
    parser.add_argument("--mpi-ranks", type=int, default=2)
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--tolerance", type=float, default=1.0e-10)
    parser.add_argument(
        "--row-timeout-sec",
        type=float,
        default=120.0,
        help="Maximum wall time for each serial or MPI side of a smoke row.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Remove an existing output directory before running.",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    output = (
        args.output.resolve()
        if args.output
        else root.parent / f"{root.name}_mpi_explicit_agreement"
    )
    if output.exists():
        if not args.force:
            raise SystemExit(f"Output directory already exists: {output}. Use --force to replace it.")
        shutil.rmtree(output)
    output.mkdir(parents=True)

    serial_exe = root / "amrex/apps/euler_compare/amrex_euler_compare2d.gnu.ex"
    mpi_exe = root / "amrex/apps/euler_compare/amrex_euler_compare2d.gnu.MPI.ex"
    missing = [path for path in (serial_exe, mpi_exe) if not path.exists()]
    if missing:
        raise SystemExit("Missing executable(s): " + ", ".join(str(path) for path in missing))
    mpirun = shutil.which(args.mpirun)
    if not mpirun:
        raise SystemExit(f"MPI launcher not found on PATH: {args.mpirun}")

    rows: list[dict[str, str]] = []
    for case in explicit_cases():
        row = run_case(
            case=case,
            root=root,
            output=output,
            serial_exe=serial_exe,
            mpi_exe=mpi_exe,
            mpirun=mpirun,
            mpi_ranks=args.mpi_ranks,
            tolerance=args.tolerance,
            row_timeout_sec=args.row_timeout_sec,
        )
        rows.append(row)
        print(
            row["row_id"],
            row["serial_returncode"],
            row["mpi_returncode"],
            row["serial_status"],
            row["mpi_status"],
            row["max_abs_diff"],
            row["within_tolerance"],
            flush=True,
        )

    summary = output / "summary.csv"
    fields = list(rows[0].keys()) if rows else []
    with summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    failed = [row for row in rows if row["within_tolerance"] != "yes"]
    print(f"summary {summary}")
    if failed:
        print("failed rows: " + ", ".join(row["row_id"] for row in failed), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
