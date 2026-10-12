#!/usr/bin/env python3
"""Summarize IsosurfacesSweep CSV; no expected outcomes are hard-coded.

Usage:
    python tools/summarize_sweep.py outputs/parameter_sweep.csv
    python tools/summarize_sweep.py outputs/parameter_sweep.csv --markdown results.md
"""
from __future__ import annotations
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


METHODS = (
    "dual_contouring",
    "local_unfolder_secondary",
    "local_voxel_patch",
    "local_voxel_stitched",
)


def number(row: dict[str, str], key: str) -> float | None:
    raw = row.get(key, "")
    if not raw:
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def key_for(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row.get(name, "") for name in
                 ("surface", "c", "resolution", "extent", "radius", "level"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--markdown", type=Path, help="also write a Markdown summary table")
    args = parser.parse_args()
    if not args.csv_file.is_file():
        parser.error(f"CSV not found: {args.csv_file}")
    with args.csv_file.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        parser.error("CSV has no data rows")

    baseline_by_key: dict[tuple[str, ...], dict[str, str]] = {}
    by_surface: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        if row.get("method") == "classic_mc":
            baseline_by_key[key_for(row)] = row
        by_surface[row.get("surface", "")][row.get("method", "")].append(row)

    table = []
    for surface in sorted({row.get("surface", "") for row in rows}):
        methods = by_surface[surface]
        mc_rows = methods.get("classic_mc", [])
        dc_rows = methods.get("dual_contouring", [])
        stitched_rows = methods.get("local_voxel_stitched", [])
        success_rows = [row for row in stitched_rows
                        if number(row, "stitch_success") is not None]
        successes = sum(number(row, "stitch_success") == 1.0 for row in success_rows)
        success_rate = 100.0 * successes / len(success_rows) if success_rows else None
        method_reductions: dict[str, float | None] = {}
        for method in METHODS:
            reductions = []
            for row in methods.get(method, []):
                baseline = baseline_by_key.get(key_for(row))
                baseline_rms = number(baseline or {}, "local_triangle_error_rms")
                method_rms = number(row, "local_triangle_error_rms")
                if baseline_rms is not None and method_rms is not None and baseline_rms > 1e-15:
                    reductions.append(100.0 * (1.0 - method_rms / baseline_rms))
            method_reductions[method] = sum(reductions) / len(reductions) if reductions else None

        stitched_nonmanifold = [number(row, "nonmanifold_edges") for row in stitched_rows]
        stitched_nonmanifold = [value for value in stitched_nonmanifold if value is not None]
        dc_nonmanifold = [number(row, "nonmanifold_edges") for row in dc_rows]
        dc_nonmanifold = [value for value in dc_nonmanifold if value is not None]
        table.append({
            "surface": surface,
            "mc_records": len(mc_rows),
            "dc_records": len(dc_rows),
            "stitch_success_rate": success_rate,
            "stitch_attempts": len(success_rows),
            "stitch_failures": sum(number(row, "stitch_success") == 0.0 for row in success_rows),
            "dc_reduction": method_reductions["dual_contouring"],
            "unfolder_reduction": method_reductions["local_unfolder_secondary"],
            "patch_reduction": method_reductions["local_voxel_patch"],
            "stitched_reduction": method_reductions["local_voxel_stitched"],
            "max_stitched_nonmanifold": max(stitched_nonmanifold) if stitched_nonmanifold else None,
            "max_dc_nonmanifold": max(dc_nonmanifold) if dc_nonmanifold else None,
        })

    headers = [
        "surface", "MC records", "DC records", "stitch success", "stitch attempts", "stitch failures",
        "DC local-RMS reduction vs MC", "LocalUnfolder local-RMS reduction vs MC",
        "voxel-patch local-RMS reduction vs MC", "stitched local-RMS reduction vs MC",
        "max stitched non-manifold edges", "max DC non-manifold edges",
    ]
    output = ["\t".join(headers)]
    def fmt_pct(value: float | None) -> str:
        return "n/a" if value is None else f"{value:.1f}%"
    def fmt_num(value: float | None) -> str:
        return "n/a" if value is None else f"{value:.5g}"
    for row in table:
        output.append("\t".join([
            row["surface"], str(row["mc_records"]), str(row["dc_records"]),
            fmt_pct(row["stitch_success_rate"]), str(row["stitch_attempts"]), str(row["stitch_failures"]),
            fmt_pct(row["dc_reduction"]), fmt_pct(row["unfolder_reduction"]),
            fmt_pct(row["patch_reduction"]), fmt_pct(row["stitched_reduction"]),
            fmt_num(row["max_stitched_nonmanifold"]), fmt_num(row["max_dc_nonmanifold"]),
        ]))
    text = "\n".join(output)
    print(text)
    if args.markdown:
        md = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
        md.extend("| " + " | ".join(line.split("\t")) + " |" for line in output[1:])
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text("\n".join(md) + "\n", encoding="utf-8")
        print(f"Wrote Markdown summary: {args.markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
