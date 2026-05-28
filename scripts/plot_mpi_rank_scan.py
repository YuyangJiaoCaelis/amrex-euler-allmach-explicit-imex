#!/usr/bin/env python3
"""Plot and summarise MPI rank-scan timing rows."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CASE_META = {
    "candidate_riemann_sod_n400_explicit_hllc": {
        "family": "Riemann Sod",
        "scheme": "HLLC",
        "label": "Riemann Sod HLLC",
    },
    "candidate_riemann_sod_n400_explicit_lowmach_hllcp": {
        "family": "Riemann Sod",
        "scheme": "Low-Mach HLLC-P",
        "label": "Riemann Sod low-Mach HLLC-P",
    },
    "candidate_gresho_m0p01_n32_explicit_hllc": {
        "family": "Gresho",
        "scheme": "HLLC",
        "label": "Gresho HLLC",
    },
    "candidate_gresho_m0p01_n32_explicit_lowmach_hllcp": {
        "family": "Gresho",
        "scheme": "Low-Mach HLLC-P",
        "label": "Gresho low-Mach HLLC-P",
    },
    "candidate_advection_blob_n128_explicit_hllc": {
        "family": "Advection Blob",
        "scheme": "HLLC",
        "label": "Advection blob HLLC",
    },
    "candidate_advection_blob_n128_explicit_lowmach_hllcp": {
        "family": "Advection Blob",
        "scheme": "Low-Mach HLLC-P",
        "label": "Advection blob low-Mach HLLC-P",
    },
    "candidate_shock_density_bubble_160x40_explicit_hllc": {
        "family": "Shock-Density-Bubble",
        "scheme": "HLLC",
        "label": "Shock-density-bubble HLLC",
    },
    "candidate_shock_density_bubble_160x40_explicit_lowmach_hllcp": {
        "family": "Shock-Density-Bubble",
        "scheme": "Low-Mach HLLC-P",
        "label": "Shock-density-bubble low-Mach HLLC-P",
    },
}

FAMILY_ORDER = [
    "Riemann Sod",
    "Gresho",
    "Advection Blob",
    "Shock-Density-Bubble",
]
SCHEME_ORDER = ["HLLC", "Low-Mach HLLC-P"]
SCHEME_STYLE = {
    "HLLC": {"color": "#2f6f9f", "marker": "o"},
    "Low-Mach HLLC-P": {"color": "#b85c38", "marker": "s"},
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    try:
        out = float(value)
    except ValueError:
        return math.nan
    return out if math.isfinite(out) else math.nan


def enrich_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    enriched: list[dict[str, object]] = []
    for row in rows:
        meta = CASE_META.get(row["case_id"])
        if meta is None:
            continue
        enriched.append(
            {
                **row,
                **meta,
                "mpi_ranks_int": int(row["mpi_ranks"]),
                "driver_median": as_float(row, "driver_wall_time_median_sec"),
                "driver_best": as_float(row, "driver_wall_time_best_sec"),
                "speedup": as_float(row, "speedup_vs_mpi1_driver_median"),
            }
        )
    return enriched


def grouped_by_case(rows: list[dict[str, object]]) -> dict[str, list[dict[str, object]]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["case_id"]), []).append(row)
    for group in grouped.values():
        group.sort(key=lambda row: int(row["mpi_ranks_int"]))
    return grouped


def plot_speedup(rows: list[dict[str, object]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)
    axes_by_family = {family: ax for family, ax in zip(FAMILY_ORDER, axes.ravel())}
    grouped = grouped_by_case(rows)

    for case_id, group in grouped.items():
        family = str(group[0]["family"])
        scheme = str(group[0]["scheme"])
        style = SCHEME_STYLE[scheme]
        ax = axes_by_family[family]
        ranks = [int(row["mpi_ranks_int"]) for row in group]
        speedups = [float(row["speedup"]) for row in group]
        ax.plot(
            ranks,
            speedups,
            color=style["color"],
            marker=style["marker"],
            linewidth=2.0,
            markersize=6,
            label=scheme,
        )

    for family, ax in axes_by_family.items():
        ax.plot([1, 4], [1, 4], color="#777777", linestyle="--", linewidth=1.0, label="ideal")
        ax.axhline(1.0, color="#bbbbbb", linewidth=0.8)
        ax.set_title(family)
        ax.set_xlabel("MPI ranks")
        ax.set_ylabel("Speedup vs 1 rank")
        ax.set_xticks([1, 2, 4])
        ax.set_ylim(0.75, 4.15)
        ax.grid(True, color="#dddddd", linewidth=0.8)
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("Candidate MPI Rank-Scan Speedup", y=1.04, fontsize=15)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_wall_time(rows: list[dict[str, object]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)
    axes_by_family = {family: ax for family, ax in zip(FAMILY_ORDER, axes.ravel())}
    grouped = grouped_by_case(rows)

    for case_id, group in grouped.items():
        family = str(group[0]["family"])
        scheme = str(group[0]["scheme"])
        style = SCHEME_STYLE[scheme]
        ax = axes_by_family[family]
        ranks = [int(row["mpi_ranks_int"]) for row in group]
        times = [float(row["driver_median"]) for row in group]
        ax.plot(
            ranks,
            times,
            color=style["color"],
            marker=style["marker"],
            linewidth=2.0,
            markersize=6,
            label=scheme,
        )

    for family, ax in axes_by_family.items():
        ax.set_title(family)
        ax.set_xlabel("MPI ranks")
        ax.set_ylabel("Median driver wall time (s)")
        ax.set_xticks([1, 2, 4])
        ax.set_yscale("log")
        ax.grid(True, color="#dddddd", linewidth=0.8, which="both")
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
    fig.suptitle("Candidate MPI Rank-Scan Wall Time", y=1.04, fontsize=15)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def write_report(rows: list[dict[str, object]], output: Path, input_dir: Path) -> None:
    grouped = grouped_by_case(rows)

    lines = [
        "# Candidate MPI Rank-Scan Interpretation",
        "",
        f"Input directory: `{input_dir}`",
        "",
        "This report is generated from candidate rank-scan CSVs. It is useful",
        "for deciding frozen Report 2 timing design, but it is not final",
        "cross-machine performance evidence.",
        "",
        "## Summary Table",
        "",
        "| Case | 1-rank median s | 2-rank speedup | 4-rank speedup | Interpretation |",
        "|---|---:|---:|---:|---|",
    ]

    for case_id in sorted(grouped, key=lambda item: CASE_META[item]["label"]):
        group = {int(row["mpi_ranks_int"]): row for row in grouped[case_id]}
        one = group.get(1)
        two = group.get(2)
        four = group.get(4)
        label = CASE_META[case_id]["label"]
        one_time = float(one["driver_median"]) if one else math.nan
        speed2 = float(two["speedup"]) if two else math.nan
        speed4 = float(four["speedup"]) if four else math.nan
        family = CASE_META[case_id]["family"]
        interpretation = interpretation_for(family, speed2, speed4)
        lines.append(
            "| "
            + " | ".join(
                [
                    label,
                    f"{one_time:.4g}" if math.isfinite(one_time) else "",
                    f"{speed2:.3g}" if math.isfinite(speed2) else "",
                    f"{speed4:.3g}" if math.isfinite(speed4) else "",
                    interpretation,
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Research Reading",
            "",
            "- Riemann and advection rows give a usable first MPI scaling signal.",
            "- Low-Mach Gresho is the strongest local scaling row in this candidate set.",
            "- Shock-density-bubble rows are too small at `160 x 40` for a meaningful scaling claim.",
            "- Gresho HLLC should be treated carefully: the 4-rank median is slower than the 2-rank median.",
            "",
            "## Suggested Frozen Timing Design",
            "",
            "For the next frozen evidence run, prefer at least 5 repeats and retain",
            "`1,2,4` ranks. Consider either increasing the shock grid or labelling",
            "shock timing as launch/overhead dominated. Keep IMEX timing separate",
            "as serial CPU unless the pressure-solve route changes.",
        ]
    )
    output.write_text("\n".join(lines) + "\n")


def interpretation_for(family: str, speed2: float, speed4: float) -> str:
    if family == "Shock-Density-Bubble":
        return "Too small for scaling claim"
    if math.isfinite(speed2) and math.isfinite(speed4) and speed4 < speed2:
        return "4-rank regression; inspect overhead/domain split"
    if math.isfinite(speed4) and speed4 >= 2.0:
        return "Useful candidate scaling"
    if math.isfinite(speed4) and speed4 >= 1.5:
        return "Moderate candidate scaling"
    return "Weak candidate scaling"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Rank-scan output directory.")
    parser.add_argument("--output", type=Path, default=None, help="Figure/report output directory.")
    args = parser.parse_args()

    input_dir = args.input.resolve()
    rank_summary = input_dir / "rank_summary.csv"
    detail_summary = input_dir / "summary.csv"
    if not rank_summary.exists() or not detail_summary.exists():
        raise SystemExit(f"Missing summary CSVs in {input_dir}")

    output_dir = args.output.resolve() if args.output else input_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = enrich_rows(read_csv(rank_summary))
    if not rows:
        raise SystemExit(f"No recognised candidate rows in {rank_summary}")
    detail_rows = read_csv(detail_summary)
    if any(row.get("passed") != "yes" for row in detail_rows):
        raise SystemExit("Refusing to plot: at least one detailed rank-scan row did not pass.")

    plot_speedup(rows, output_dir / "mpi_rank_scan_speedup.png")
    plot_wall_time(rows, output_dir / "mpi_rank_scan_wall_time.png")
    write_report(rows, output_dir / "mpi_rank_scan_interpretation.md", input_dir)

    print(output_dir / "mpi_rank_scan_speedup.png")
    print(output_dir / "mpi_rank_scan_wall_time.png")
    print(output_dir / "mpi_rank_scan_interpretation.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
