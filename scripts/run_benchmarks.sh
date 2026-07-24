#!/usr/bin/env bash
# Runs the Google Benchmark suite and records results as timestamped JSON
# under benchmarks/results/, plus a "latest.json" symlink, so successive
# runs can be diffed for the >5% regression gate (CLAUDE.md #5).
#
# Usage: scripts/run_benchmarks.sh [build-dir]   (default: build/bench)
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${1:-${repo_root}/build/bench}"
bench_bin="${build_dir}/benchmarks/aleator_bench"

if [ ! -x "${bench_bin}" ]; then
    echo "error: ${bench_bin} not found or not executable." >&2
    echo "Build it first: cmake --preset bench && cmake --build --preset bench" >&2
    exit 1
fi

results_dir="${repo_root}/benchmarks/results"
mkdir -p "${results_dir}"

timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
out_file="${results_dir}/${timestamp}.json"

"${bench_bin}" --benchmark_out="${out_file}" --benchmark_out_format=json

ln -sf "$(basename "${out_file}")" "${results_dir}/latest.json"
echo "wrote ${out_file} (and updated ${results_dir}/latest.json)"
