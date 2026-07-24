#!/usr/bin/env bash
# Clones and bootstraps vcpkg into .vcpkg-cibw/ at the repo root, pinned to
# the baseline commit in vcpkg-configuration.json, so that wheel builds
# (cibuildwheel) and CI use the exact same dependency versions as local
# development. Safe to re-run; it no-ops if already bootstrapped.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
vcpkg_dir="${repo_root}/.vcpkg-cibw"
baseline="$(grep -o '"baseline": *"[0-9a-f]*"' "${repo_root}/vcpkg-configuration.json" | grep -o '[0-9a-f]\{40\}')"

if [ ! -d "${vcpkg_dir}/.git" ]; then
    git clone https://github.com/microsoft/vcpkg.git "${vcpkg_dir}"
fi

git -C "${vcpkg_dir}" fetch --depth 1 origin "${baseline}"
git -C "${vcpkg_dir}" checkout "${baseline}"

if [ ! -x "${vcpkg_dir}/vcpkg" ] && [ ! -x "${vcpkg_dir}/vcpkg.exe" ]; then
    if [ -f "${vcpkg_dir}/bootstrap-vcpkg.bat" ] && [ "${OS:-}" = "Windows_NT" ]; then
        "${vcpkg_dir}/bootstrap-vcpkg.bat" -disableMetrics
    else
        "${vcpkg_dir}/bootstrap-vcpkg.sh" -disableMetrics
    fi
fi

echo "vcpkg ready at ${vcpkg_dir}"
