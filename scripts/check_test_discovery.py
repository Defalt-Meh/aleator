#!/usr/bin/env python3
"""Fails the build if CTest is not discovering every TEST_CASE.

CLAUDE.md invariant 8: "Every test is individually visible to the runner."
This is the actual guard, not the TEST_CASE renaming that used to be needed
to satisfy it -- if catch_discover_tests (or a future Catch2/CMake upgrade)
ever again collapses tests into fewer CTest entries than exist in the
source, this script is what turns that into a failed CI build instead of a
silently shrinking test suite.

Usage: check_test_discovery.py --build-dir build/dev --tests-dir tests
"""

import argparse
import pathlib
import re
import subprocess
import sys

TEST_CASE_PATTERN = re.compile(r"\bTEST_CASE\(")
TOTAL_TESTS_PATTERN = re.compile(r"^Total Tests:\s*(\d+)\s*$", re.MULTILINE)


def count_test_case_macros(tests_dir: pathlib.Path) -> int:
    total = 0
    for path in sorted(tests_dir.rglob("*.cc")):
        text = path.read_text(encoding="utf-8")
        total += len(TEST_CASE_PATTERN.findall(text))
    return total


def count_ctest_discovered(build_dir: pathlib.Path) -> int:
    result = subprocess.run(
        ["ctest", "--test-dir", str(build_dir), "-N"],
        capture_output=True,
        text=True,
        check=True,
    )
    match = TOTAL_TESTS_PATTERN.search(result.stdout)
    if match is None:
        print("check_test_discovery: could not find 'Total Tests: N' in `ctest -N` output:",
              file=sys.stderr)
        print(result.stdout, file=sys.stderr)
        sys.exit(1)
    return int(match.group(1))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", required=True, type=pathlib.Path)
    parser.add_argument("--tests-dir", required=True, type=pathlib.Path)
    args = parser.parse_args()

    source_count = count_test_case_macros(args.tests_dir)
    discovered_count = count_ctest_discovered(args.build_dir)

    print(f"TEST_CASE macros under {args.tests_dir}: {source_count}")
    print(f"Tests discovered by `ctest -N` in {args.build_dir}: {discovered_count}")

    if source_count != discovered_count:
        print(
            "check_test_discovery: MISMATCH. CTest is not seeing every TEST_CASE "
            "-- some tests are invisible to `ctest` and are not being run by CI, "
            "even though they exist in the source and may pass when the test "
            "binary is run directly. This is exactly the failure mode CLAUDE.md "
            "invariant 8 exists to catch. Do not relax this check; find out why "
            "discovery is undercounting.",
            file=sys.stderr,
        )
        return 1

    print("check_test_discovery: OK -- every TEST_CASE is individually visible to CTest.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
