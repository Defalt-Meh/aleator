#!/usr/bin/env python3
"""Compares two Google Benchmark JSON result files and fails (exit 1) if any
matched benchmark regressed real_time by more than --threshold (default 5%,
matching the CI gate in CLAUDE.md #5). Benchmarks present in only one file
are reported but do not cause failure — that happens when a benchmark is
added or removed, not a regression of an existing one.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def load_benchmarks(path: Path) -> dict[str, dict]:
    data = json.loads(path.read_text())
    return {entry["name"]: entry for entry in data.get("benchmarks", [])}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("current", type=Path)
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.05,
        help="Maximum allowed fractional regression in real_time (default: 0.05 = 5%%).",
    )
    args = parser.parse_args()

    baseline = load_benchmarks(args.baseline)
    current = load_benchmarks(args.current)

    regressed = False
    print(f"{'benchmark':<40} {'baseline':>12} {'current':>12} {'delta':>10}")
    for name in sorted(set(baseline) & set(current)):
        base_time = baseline[name]["real_time"]
        cur_time = current[name]["real_time"]
        if base_time <= 0:
            continue
        delta = (cur_time - base_time) / base_time
        flag = ""
        if delta > args.threshold:
            regressed = True
            flag = "  REGRESSION"
        print(f"{name:<40} {base_time:>12.2f} {cur_time:>12.2f} {delta:>+9.1%}{flag}")

    only_baseline = sorted(set(baseline) - set(current))
    only_current = sorted(set(current) - set(baseline))
    if only_baseline:
        print(f"\nonly in baseline (removed): {', '.join(only_baseline)}")
    if only_current:
        print(f"only in current (added): {', '.join(only_current)}")

    if regressed:
        print(f"\nFAIL: one or more benchmarks regressed by more than {args.threshold:.0%}.")
        return 1

    print("\nOK: no benchmark regressed beyond threshold.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
