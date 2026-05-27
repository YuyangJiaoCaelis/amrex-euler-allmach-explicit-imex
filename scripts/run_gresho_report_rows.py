#!/usr/bin/env python3
"""Run the uniform Gresho rows used by the report-facing figures."""

from __future__ import annotations

import argparse
import csv
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts/run_gresho_reference_matrix.py"
OUT_ROOT = ROOT / "results/amrex/project_gresho_reproduction_report_rows"


@dataclass(frozen=True)
class Group:
    name: str
    scheme: str
    mach: str
    grids: tuple[int, ...]
    max_grid_size: int


GROUPS = (
    Group("explicit_hllc_m0p001", "explicit_o2_hllc", "0.001", (32, 40, 48, 56, 64), 64),
    Group("explicit_hllc_m0p01", "explicit_o2_hllc", "0.01", (32, 48, 64, 96, 128), 128),
    Group("explicit_hllc_m0p1", "explicit_o2_hllc", "0.1", (32, 48, 64, 96, 128, 192, 256), 256),
    Group("lowmach_hllcp_m0p001", "explicit_o2_xie_am_hllc_p", "0.001", (24, 28, 32, 36, 40, 48), 48),
    Group("lowmach_hllcp_m0p01", "explicit_o2_xie_am_hllc_p", "0.01", (32, 40, 48, 56), 56),
    Group("lowmach_hllcp_m0p1", "explicit_o2_xie_am_hllc_p", "0.1", (32, 48, 64, 96), 96),
    Group("imex_t1s2_m0p001", "imex_t1s2_bdltv20", "0.001", (32, 48, 64, 96, 128, 192, 256), 256),
    Group("imex_t1s2_m0p01", "imex_t1s2_bdltv20", "0.01", (32, 48, 64, 96, 128), 128),
    Group("imex_t1s2_m0p1", "imex_t1s2_bdltv20", "0.1", (32, 48, 64, 96, 128), 128),
)


def summary_complete(path: Path, expected_count: int) -> bool:
    if not path.exists():
        return False
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    ok_rows = [row for row in rows if row.get("status") in {"ok", "failed_quality_gate"}]
    return len(ok_rows) == expected_count


def run_group(group: Group, timeout: float, force: bool) -> None:
    out_dir = OUT_ROOT / group.name
    summary = out_dir / "summary.csv"
    if not force and summary_complete(summary, len(group.grids)):
        print(f"[skip] {group.name}: existing complete summary")
        return
    cmd = [
        "python3",
        str(RUNNER),
        "--output-dir",
        str(out_dir),
        "--schemes",
        group.scheme,
        "--machs",
        group.mach,
        "--grids",
        ",".join(str(grid) for grid in group.grids),
        "--field-boundary",
        "exact_dirichlet",
        "--target-time",
        "1.2566370614359172",
        "--max-grid-size",
        str(group.max_grid_size),
        "--row-timeout-sec",
        str(timeout),
        "--explicit-max-step",
        "1000000",
        "--imex-max-step",
        "200000",
        "--clean-imex-cfl",
        "0.4",
    ]
    print("[run]", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--row-timeout-sec", type=float, default=3600.0)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--only", default="", help="Comma-separated group names to run.")
    args = parser.parse_args()

    only = {item.strip() for item in args.only.split(",") if item.strip()}
    selected = [group for group in GROUPS if not only or group.name in only]
    unknown = only - {group.name for group in GROUPS}
    if unknown:
        raise SystemExit(f"unknown group(s): {', '.join(sorted(unknown))}")
    for group in selected:
        run_group(group, args.row_timeout_sec, args.force)
    print(OUT_ROOT)


if __name__ == "__main__":
    main()
