// Black-box integration tests for the actual `aleator` CLI binary
// (CLAUDE.md milestone: "every example in examples/ must run end to end;
// every malformed config in tests/ must produce a helpful error"). These
// invoke the real, built executable as a subprocess and check its exit
// code and output -- the only way to genuinely test argument parsing,
// --dry-run/--json behavior, and process exit codes, none of which are
// observable by linking against the library code directly.

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#ifdef _WIN32
#include <cstdlib>
#else
#include <sys/wait.h>
#endif

namespace {

struct CommandResult {
    int exitCode;
    std::string output;
};

/// Runs `command`, capturing combined stdout+stderr, and returns its exit
/// code. Uses popen/_popen (POSIX/MSVC both provide a variant) rather than
/// pipes+fork directly -- fine for a test that just needs "did it succeed,
/// and what did it print," not a general subprocess library.
CommandResult runCommand(const std::string& command) {
    const std::string redirected = command + " 2>&1";
#ifdef _WIN32
    FILE* pipe = _popen(redirected.c_str(), "r");
#else
    FILE* pipe = popen(redirected.c_str(), "r");
#endif
    if (pipe == nullptr) {
        return {-1, "<failed to launch process>"};
    }
    std::ostringstream output;
    std::array<char, 4096> buffer{};
    while (std::fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) != nullptr) {
        output << buffer.data();
    }
#ifdef _WIN32
    const int rawStatus = _pclose(pipe);
    const int exitCode = rawStatus;
#else
    // Non-const: macOS's WIFEXITED/WEXITSTATUS macros take the address of
    // their argument internally, which trips -Wcast-qual under -Werror if
    // the argument itself is const.
    int rawStatus = pclose(pipe);
    const int exitCode = WIFEXITED(rawStatus) ? WEXITSTATUS(rawStatus) : -1;
#endif
    return {exitCode, output.str()};
}

std::filesystem::path cliBinary() { return std::filesystem::path(ALEATOR_CLI_BINARY_PATH); }
std::filesystem::path examplesDir() { return std::filesystem::path(ALEATOR_EXAMPLES_DIR); }

std::filesystem::path writeTemp(const std::string& name, const std::string& contents) {
    const auto path = std::filesystem::temp_directory_path() / name;
    std::ofstream out(path);
    out << contents;
    return path;
}

} // namespace

TEST_CASE("aleator --help exits 0 and documents only subcommands backed by working engines",
          "[integration][cli]") {
    // CLAUDE.md CLI milestone: "expose ONLY subcommands backed by working
    // engines... adding a subcommand that prints NotImplemented is worse
    // than omitting it -- it advertises capability the tool doesn't have."
    // pore/md config *schemas* still exist and `validate` still recognizes
    // them (see below), but `pore analyze`/`md run` must not appear as if
    // they were real, runnable subcommands.
    const auto result = runCommand("\"" + cliBinary().string() + "\" --help");
    REQUIRE(result.exitCode == 0);
    REQUIRE(result.output.find("gcmc run") != std::string::npos);
    REQUIRE(result.output.find("validate") != std::string::npos);
    REQUIRE(result.output.find("bench") != std::string::npos);
    // Not just absent from a command list -- not present as a runnable
    // usage line at all (the help text does mention the bare words "pore
    // analyze"/"md run" in prose, explaining why they don't exist; that's
    // fine and intended. What must never appear is the usage-line form,
    // implying they're invokable the way `aleator gcmc run <config.toml>`
    // is.)
    REQUIRE(result.output.find("aleator pore analyze") == std::string::npos);
    REQUIRE(result.output.find("aleator md run") == std::string::npos);
}

TEST_CASE("aleator with no arguments exits 1", "[integration][cli]") {
    const auto result = runCommand("\"" + cliBinary().string() + "\"");
    REQUIRE(result.exitCode == 1);
}

TEST_CASE("aleator with an unrecognized command exits 1 with a clear message",
          "[integration][cli]") {
    const auto result = runCommand("\"" + cliBinary().string() + "\" frobnicate");
    REQUIRE(result.exitCode == 1);
    REQUIRE(result.output.find("unrecognized command") != std::string::npos);
}

TEST_CASE("aleator pore analyze / aleator md run are unrecognized commands, not NotImplemented "
          "subcommands",
          "[integration][cli]") {
    // CLAUDE.md CLI milestone: these must not exist as subcommands at all
    // (see the --help test above) -- confirm the dispatcher genuinely
    // doesn't recognize them, rather than silently routing "pore"/"md" to
    // something.
    const auto poreResult = runCommand("\"" + cliBinary().string() + "\" pore analyze x.toml");
    REQUIRE(poreResult.exitCode == 1);
    REQUIRE(poreResult.output.find("unrecognized command \"pore\"") != std::string::npos);

    const auto mdResult = runCommand("\"" + cliBinary().string() + "\" md run x.toml");
    REQUIRE(mdResult.exitCode == 1);
    REQUIRE(mdResult.output.find("unrecognized command \"md\"") != std::string::npos);
}

TEST_CASE("examples/gcmc_ch4_irmof1.toml validates and dry-runs end to end, printing every "
          "default (including ones this milestone added)",
          "[integration][cli]") {
    const auto configPath = examplesDir() / "gcmc_ch4_irmof1.toml";
    REQUIRE(std::filesystem::exists(configPath));

    const auto validateResult =
        runCommand("\"" + cliBinary().string() + "\" validate \"" + configPath.string() + "\"");
    REQUIRE(validateResult.exitCode == 0);
    REQUIRE(validateResult.output.find("valid GCMC config") != std::string::npos);
    REQUIRE(validateResult.output.find("424 framework atoms") != std::string::npos);

    const auto dryRunResult = runCommand("\"" + cliBinary().string() + "\" gcmc run \"" +
                                          configPath.string() + "\" --dry-run --json");
    REQUIRE(dryRunResult.exitCode == 0);
    REQUIRE(dryRunResult.output.find("\"framework_cif\"") != std::string::npos);
    REQUIRE(dryRunResult.output.find("\"uses_peng_robinson\":true") != std::string::npos);
    // Every default filled in, not just the keys the config set explicitly:
    // output_directory and energy_grid_spacing_angstrom are both left at
    // their defaults in this example config.
    REQUIRE(dryRunResult.output.find("\"output_directory\":\"") != std::string::npos);
    REQUIRE(dryRunResult.output.find("\"energy_grid_spacing_angstrom\":null") != std::string::npos);
    // --dry-run must not write anything.
    REQUIRE_FALSE(std::filesystem::exists(examplesDir() / "out" / "resolved_config.json"));
}

TEST_CASE("examples/gcmc_ch4_irmof1.toml runs the real, validated 0.1 bar GCMC point end to "
          "end, matching the validation test's numbers, and writes reproducible artifacts",
          "[integration][cli][slow]") {
    const auto configPath = examplesDir() / "gcmc_ch4_irmof1.toml";
    const auto outputDir = examplesDir() / "out";
    std::filesystem::remove_all(outputDir);

    const auto result =
        runCommand("\"" + cliBinary().string() + "\" gcmc run \"" + configPath.string() + "\"");
    REQUIRE(result.exitCode == 0);
    REQUIRE(result.output.find("Final result") != std::string::npos);

    // This example's rng_seed (2026), pressure_bar (0.1), and step counts
    // (30000 equilibration + 90000 production) are IDENTICAL to
    // tests/validation/test_gcmc_ch4_irmof1_isotherm.cc's real 0.1 bar
    // point, and runGcmcCommand() samples <N> after every production step
    // (not once per progress update), the same estimator that test uses --
    // so this is not an approximate match, it's the same deterministic
    // computation (CLAUDE.md invariant #5) reproduced through the CLI.
    // Cross-checked directly against a fresh run of that validation test in
    // this same session/build: both report <N> = 0.2315 (4 s.f.).
    REQUIRE(result.output.find("0.231456") != std::string::npos);
    REQUIRE(result.output.find("0.0375804 mmol/g") != std::string::npos);

    // Reproducibility artifacts (CLAUDE.md CLI milestone: "a published
    // result can be reproduced from its artifacts alone").
    const auto resolvedConfigPath = outputDir / "resolved_config.json";
    REQUIRE(std::filesystem::exists(resolvedConfigPath));
    std::ifstream artifact(resolvedConfigPath);
    std::ostringstream artifactContent;
    artifactContent << artifact.rdbuf();
    REQUIRE(artifactContent.str().find("\"rng_seed\":2026") != std::string::npos);
    REQUIRE(artifactContent.str().find("\"pressure_bar\":0.10000000000000001") !=
            std::string::npos);

    std::filesystem::remove_all(outputDir);
}

TEST_CASE("aleator validate on examples/pore_irmof1.toml is schema-valid but honestly reports "
          "no engine exists for it yet",
          "[integration][cli]") {
    const auto configPath = examplesDir() / "pore_irmof1.toml";
    REQUIRE(std::filesystem::exists(configPath));

    const auto result =
        runCommand("\"" + cliBinary().string() + "\" validate \"" + configPath.string() + "\"");
    REQUIRE(result.exitCode == 2);
    REQUIRE(result.output.find("not implemented") != std::string::npos);
    REQUIRE(result.output.find("schema-valid") != std::string::npos);
}

TEST_CASE("aleator validate on examples/md_irmof1.toml is schema-valid but honestly reports no "
          "engine exists for it yet",
          "[integration][cli]") {
    const auto configPath = examplesDir() / "md_irmof1.toml";
    REQUIRE(std::filesystem::exists(configPath));

    const auto result =
        runCommand("\"" + cliBinary().string() + "\" validate \"" + configPath.string() + "\"");
    REQUIRE(result.exitCode == 2);
    REQUIRE(result.output.find("not implemented") != std::string::npos);
    REQUIRE(result.output.find("schema-valid") != std::string::npos);
}

TEST_CASE("aleator bench runs and reports timing", "[integration][cli]") {
    const auto result = runCommand("\"" + cliBinary().string() + "\" bench --json");
    REQUIRE(result.exitCode == 0);
    REQUIRE(result.output.find("\"particle_count\"") != std::string::npos);
}

TEST_CASE("a malformed config produces a line-numbered error and exit code 1 from the real CLI",
          "[integration][cli]") {
    const auto path = writeTemp("aleator_cli_bad_config.toml",
                                 "[run]\nname = \"test\"\n\n[gcmc]\nframework_cif = \"x.cif\"\n"
                                 "temperature_kelvin = \"hot\"\npressure_bar = 1.0\n\n"
                                 "[gcmc.adsorbate]\nname = \"CH4\"\nepsilon_kelvin = 158.5\n"
                                 "sigma_angstrom = 3.72\nmass_amu = 16.04246\n\n"
                                 "[[gcmc.framework_lj]]\nelement = \"Zn\"\nepsilon_kelvin = 62.4\n"
                                 "sigma_angstrom = 2.46\n");

    const auto result =
        runCommand("\"" + cliBinary().string() + "\" validate \"" + path.string() + "\"");
    REQUIRE(result.exitCode == 1);
    REQUIRE(result.output.find("gcmc.temperature_kelvin") != std::string::npos);
    REQUIRE(result.output.find("expected a number, found string") != std::string::npos);
    // A "file:line:column:" prefix -- at least two colons before the key.
    REQUIRE(result.output.find(":6:") != std::string::npos);

    std::filesystem::remove(path);
}

TEST_CASE("a framework_lj that doesn't cover every CIF element is rejected with the element name",
          "[integration][cli]") {
    // IRMOF-1 contains Zn, O, C, H -- omit H.
    const auto configPath = examplesDir() / "IRMOF-1.cif";
    const auto path = writeTemp(
        "aleator_cli_incomplete_framework_lj.toml",
        "[run]\nname = \"test\"\n\n[gcmc]\nframework_cif = \"" + configPath.string() +
            "\"\ntemperature_kelvin = 298.0\npressure_bar = 1.0\n\n"
            "[gcmc.adsorbate]\nname = \"CH4\"\nepsilon_kelvin = 158.5\nsigma_angstrom = 3.72\n"
            "mass_amu = 16.04246\n\n"
            "[[gcmc.framework_lj]]\nelement = \"Zn\"\nepsilon_kelvin = 62.3992\n"
            "sigma_angstrom = 2.46155\n\n"
            "[[gcmc.framework_lj]]\nelement = \"O\"\nepsilon_kelvin = 48.1581\n"
            "sigma_angstrom = 3.03315\n\n"
            "[[gcmc.framework_lj]]\nelement = \"C\"\nepsilon_kelvin = 47.8562\n"
            "sigma_angstrom = 3.47299\n");

    const auto result =
        runCommand("\"" + cliBinary().string() + "\" validate \"" + path.string() + "\"");
    REQUIRE(result.exitCode == 1);
    REQUIRE(result.output.find("\"H\"") != std::string::npos);

    std::filesystem::remove(path);
}
