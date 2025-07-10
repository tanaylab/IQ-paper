#!/usr/bin/env python
# (content identical to previous collect_r2.py)

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable, List

R2_LINE_RE = re.compile(r"^(?P<traj>.+?)\s+r\^2:\s+(?P<val>[-+]?\d*\.\d+|\d+)")
MEAN_RE = re.compile(r"^Mean r\^2:\s+(?P<val>[-+]?\d*\.\d+|\d+)")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Collect test R² files and output as CSV")
    p.add_argument(
        "--runs",
        nargs="+",
        required=True,
        help="One or more output directories that contain a test_correlations.txt. You can also pass glob patterns (shell-expanded)",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("r2_summary.csv"),
        help="Where to write the aggregated CSV (default: r2_summary.csv)",
    )
    return p.parse_args()


def find_correlation_file(run_dir: Path) -> Path | None:
    candidate = run_dir / "test_correlations.txt"
    return candidate if candidate.exists() else None


def parse_file(model_name: str, txt_path: Path, rows: List[dict[str, str | float]]):
    with txt_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            m = R2_LINE_RE.match(line)
            if m:
                rows.append({"model": model_name, "trajectory": m.group("traj"), "r2": float(m.group("val"))})
                continue
            m = MEAN_RE.match(line)
            if m:
                rows.append({"model": model_name, "trajectory": "MEAN", "r2": float(m.group("val"))})


def expand_globs(patterns: Iterable[str]) -> List[Path]:
    import glob as _glob

    paths: List[Path] = []
    for pat in patterns:
        paths.extend(Path(p) for p in _glob.glob(pat))
    return paths


def main() -> None:
    args = parse_args()
    run_paths = expand_globs(args.runs)
    if not run_paths:
        raise SystemExit("No run directories matched the given patterns.")

    rows: List[dict[str, str | float]] = []
    for run in run_paths:
        corr_file = find_correlation_file(run)
        if corr_file is None:
            print(f"[WARN] {run} has no test_correlations.txt – skipped")
            continue
        parse_file(model_name=run.name, txt_path=corr_file, rows=rows)

    if not rows:
        raise SystemExit("No R² values found – nothing to write.")

    with args.output.open("w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["model", "trajectory", "r2"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main() 