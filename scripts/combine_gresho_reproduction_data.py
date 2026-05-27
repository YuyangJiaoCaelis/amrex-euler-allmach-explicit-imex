#!/usr/bin/env python3
"""Combine the MPhil project Gresho runner summaries for plotting."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

REPORT_ROW_ROOT = ROOT / "results/amrex/project_gresho_reproduction_report_rows"

OUTPUT = (
    ROOT
    / "project_outputs/test_packages/gresho_vortex/source_data/"
    / "amrex_gresho_efficiency_combined_data.csv"
)


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing Gresho summary CSV: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def main() -> None:
    rows: list[dict[str, str]] = []
    fieldnames: list[str] = []
    seen: set[str] = set()

    inputs = sorted(REPORT_ROW_ROOT.glob("*/summary.csv"))
    if not inputs:
        raise FileNotFoundError(
            "No Gresho summaries found. Run scripts/run_gresho_report_rows.py first."
        )

    for path in inputs:
        for row in read_rows(path):
            rows.append(row)
            for key in row.keys():
                if key not in seen:
                    seen.add(key)
                    fieldnames.append(key)

    if not rows:
        raise RuntimeError("No Gresho rows were found in the summary CSV files.")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Wrote {len(rows)} Gresho rows to {OUTPUT}")


if __name__ == "__main__":
    main()
