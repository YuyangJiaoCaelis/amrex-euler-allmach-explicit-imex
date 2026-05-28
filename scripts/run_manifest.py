#!/usr/bin/env python3
"""Run-provenance helpers for Report 2 evidence rows.

The module is importable from runner scripts and can also wrap a single command:

    python3 scripts/run_manifest.py run --row-id row --output-dir out -- command ...
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import platform
import shlex
import socket
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def shell_join(command: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in command)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, root: Path | None = None) -> dict[str, Any]:
    resolved = path.resolve() if path.exists() else path
    try:
        display = str(resolved.relative_to(root.resolve())) if root else str(resolved)
    except ValueError:
        display = str(resolved)
    record: dict[str, Any] = {
        "path": display,
        "exists": path.exists(),
    }
    if path.exists() and path.is_file():
        stat = path.stat()
        record.update(
            {
                "size_bytes": stat.st_size,
                "sha256": sha256_file(path),
            }
        )
    return record


def git_info(root: Path) -> dict[str, Any]:
    def git(*args: str) -> str:
        proc = subprocess.run(
            ["git", *args],
            cwd=root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        return proc.stdout.strip() if proc.returncode == 0 else ""

    status = git("status", "--porcelain")
    return {
        "branch": git("branch", "--show-current"),
        "commit": git("rev-parse", "HEAD"),
        "dirty": bool(status),
        "status_porcelain": status,
        "tag_nearest": git("describe", "--tags", "--always", "--dirty"),
    }


def host_info() -> dict[str, Any]:
    return {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": platform.python_version(),
        "cpu_count": os.cpu_count(),
    }


def environment_build_flags(extra: dict[str, str] | None = None) -> dict[str, str]:
    keys = [
        "DIM",
        "COMP",
        "USE_MPI",
        "USE_OMP",
        "USE_CUDA",
        "DEBUG",
        "PROFILE",
        "AMREX_HOME",
        "OMP_NUM_THREADS",
        "CUDA_VISIBLE_DEVICES",
    ]
    flags = {key: os.environ.get(key, "not_recorded") for key in keys}
    if extra:
        flags.update(extra)
    return flags


def expand_output_paths(paths: list[Path], globs: list[str]) -> list[Path]:
    out: list[Path] = list(paths)
    for pattern in globs:
        out.extend(Path(match) for match in sorted(glob.glob(pattern)))
    # Keep order stable and remove duplicates.
    seen: set[Path] = set()
    unique: list[Path] = []
    for path in out:
        key = path.resolve() if path.exists() else path
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def write_manifest(
    manifest_path: Path,
    *,
    root: Path,
    row_id: str,
    command: list[str],
    start_utc: str,
    end_utc: str,
    wall_time_s: float,
    exit_code: int | str,
    output_root: Path,
    output_class: str = "exploratory",
    input_files: list[Path] | None = None,
    output_files: list[Path] | None = None,
    build_flags: dict[str, str] | None = None,
    notes: str = "",
    extra: dict[str, Any] | None = None,
) -> None:
    input_files = input_files or []
    output_files = output_files or []
    manifest = {
        "schema": "amrex_euler_report_run_manifest_v1",
        "row_id": row_id,
        "output_class": output_class,
        "notes": notes,
        "git": git_info(root),
        "build_flags": build_flags or environment_build_flags(),
        "host": host_info(),
        "command_line": shell_join(command),
        "command_argv": command,
        "input_files": [file_record(path, root) for path in input_files],
        "start_utc": start_utc,
        "end_utc": end_utc,
        "wall_time_s": wall_time_s,
        "exit_code": exit_code,
        "output_root": str(output_root),
        "output_files": [file_record(path, root) for path in output_files],
    }
    if extra:
        manifest["extra"] = extra
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def run_with_manifest(
    command: list[str],
    *,
    root: Path,
    row_id: str,
    output_dir: Path,
    log_path: Path,
    command_path: Path | None = None,
    manifest_path: Path | None = None,
    input_files: list[Path] | None = None,
    output_files: list[Path] | None = None,
    output_globs: list[str] | None = None,
    output_class: str = "exploratory",
    notes: str = "",
    cwd: Path | None = None,
) -> int:
    command_path = command_path or output_dir / "commands" / f"{row_id}.txt"
    manifest_path = manifest_path or output_dir / "manifests" / f"{row_id}.manifest.json"
    input_files = input_files or []
    output_files = output_files or []
    output_globs = output_globs or []

    command_path.parent.mkdir(parents=True, exist_ok=True)
    command_path.write_text(shell_join(command) + "\n")

    start_utc = utc_now()
    start = time.perf_counter()
    proc = subprocess.run(
        command,
        cwd=str(cwd or root),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    wall = time.perf_counter() - start
    end_utc = utc_now()

    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(proc.stdout)

    expanded_outputs = expand_output_paths(output_files + [log_path, command_path], output_globs)
    write_manifest(
        manifest_path,
        root=root,
        row_id=row_id,
        command=command,
        start_utc=start_utc,
        end_utc=end_utc,
        wall_time_s=wall,
        exit_code=proc.returncode,
        output_root=output_dir,
        output_class=output_class,
        input_files=input_files,
        output_files=expanded_outputs,
        notes=notes,
    )
    return proc.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--root", type=Path, default=Path.cwd())
    run_parser.add_argument("--row-id", required=True)
    run_parser.add_argument("--output-dir", type=Path, required=True)
    run_parser.add_argument("--log", type=Path, required=True)
    run_parser.add_argument("--command-file", type=Path, default=None)
    run_parser.add_argument("--manifest", type=Path, default=None)
    run_parser.add_argument("--input-file", type=Path, action="append", default=[])
    run_parser.add_argument("--output-file", type=Path, action="append", default=[])
    run_parser.add_argument("--output-glob", action="append", default=[])
    run_parser.add_argument("--output-class", default=os.environ.get("RUN_OUTPUT_CLASS", "exploratory"))
    run_parser.add_argument("--notes", default="")
    run_parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    if args.mode == "run":
        command = list(args.command)
        if command and command[0] == "--":
            command = command[1:]
        if not command:
            raise SystemExit("run mode requires a command after --")
        return run_with_manifest(
            command,
            root=args.root.resolve(),
            row_id=args.row_id,
            output_dir=args.output_dir,
            log_path=args.log,
            command_path=args.command_file,
            manifest_path=args.manifest,
            input_files=args.input_file,
            output_files=args.output_file,
            output_globs=args.output_glob,
            output_class=args.output_class,
            notes=args.notes,
        )
    raise SystemExit(f"unknown mode: {args.mode}")


if __name__ == "__main__":
    raise SystemExit(main())
