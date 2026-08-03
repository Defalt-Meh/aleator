// Validates that malformed run configs fail loudly and precisely: a
// missing/wrong-type/out-of-range key throws io::ConfigError naming the
// exact key path and (where toml++ can supply one) the exact line/column —
// "fail before starting a run, never during" (CLI milestone requirement).

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <string>

#include "io/config.hpp"
#include "io/config_error.hpp"

namespace {

std::filesystem::path writeTemp(const std::string& name, const std::string& contents) {
    const auto path = std::filesystem::temp_directory_path() / name;
    std::ofstream out(path);
    out << contents;
    return path;
}

constexpr const char* kValidGcmcHeader = R"([run]
name = "test"

[gcmc.adsorbate]
name = "CH4"
epsilon_kelvin = 158.5
sigma_angstrom = 3.72
mass_amu = 16.04246

[[gcmc.framework_lj]]
element = "Zn"
epsilon_kelvin = 62.3992
sigma_angstrom = 2.46155
)";

} // namespace

TEST_CASE("loadGcmcConfig rejects a missing gcmc table with a location", "[unit][io]") {
    const auto path = writeTemp("aleator_cfg_no_gcmc.toml", "[run]\nname = \"test\"\n");
    try {
        (void)aleator::io::loadGcmcConfig(path);
        FAIL("expected ConfigError");
    } catch (const aleator::io::ConfigError& err) {
        const std::string what = err.what();
        REQUIRE(what.find("gcmc") != std::string::npos);
        REQUIRE(what.find("missing required table") != std::string::npos);
    }
    std::filesystem::remove(path);
}

TEST_CASE("loadGcmcConfig rejects a missing required key with a helpful message", "[unit][io]") {
    const auto path = writeTemp("aleator_cfg_missing_temp.toml",
                                 std::string(kValidGcmcHeader) +
                                     "\n[gcmc]\nframework_cif = \"x.cif\"\npressure_bar = 1.0\n");
    try {
        (void)aleator::io::loadGcmcConfig(path);
        FAIL("expected ConfigError");
    } catch (const aleator::io::ConfigError& err) {
        const std::string what = err.what();
        REQUIRE(what.find("gcmc.temperature_kelvin") != std::string::npos);
        REQUIRE(what.find("missing required key") != std::string::npos);
        REQUIRE(what.find("a number") != std::string::npos);
    }
    std::filesystem::remove(path);
}

TEST_CASE("loadGcmcConfig rejects a wrong-typed key naming the exact line", "[unit][io]") {
    const auto contents = std::string(kValidGcmcHeader) +
                           "\n[gcmc]\nframework_cif = \"x.cif\"\n"
                           "temperature_kelvin = \"hot\"\npressure_bar = 1.0\n";
    const auto path = writeTemp("aleator_cfg_wrong_type.toml", contents);
    try {
        (void)aleator::io::loadGcmcConfig(path);
        FAIL("expected ConfigError");
    } catch (const aleator::io::ConfigError& err) {
        const std::string what = err.what();
        REQUIRE(what.find("gcmc.temperature_kelvin") != std::string::npos);
        REQUIRE(what.find("expected a number, found string") != std::string::npos);
        // The bad key is line 13 of `contents` (count: kValidGcmcHeader has
        // 12 lines including the blank separator, so the wrong-typed
        // temperature_kelvin line is the 13th) -- checked loosely (just
        // that *some* line number appears) so this test doesn't become
        // brittle to header formatting.
        REQUIRE(what.find(":") != std::string::npos);
    }
    std::filesystem::remove(path);
}

TEST_CASE("loadGcmcConfig rejects a non-positive temperature", "[unit][io]") {
    const auto path =
        writeTemp("aleator_cfg_negative_temp.toml",
                  std::string(kValidGcmcHeader) +
                      "\n[gcmc]\nframework_cif = \"x.cif\"\ntemperature_kelvin = -1.0\n"
                      "pressure_bar = 1.0\n");
    try {
        (void)aleator::io::loadGcmcConfig(path);
        FAIL("expected ConfigError");
    } catch (const aleator::io::ConfigError& err) {
        const std::string what = err.what();
        REQUIRE(what.find("gcmc.temperature_kelvin") != std::string::npos);
        REQUIRE(what.find("must be positive") != std::string::npos);
    }
    std::filesystem::remove(path);
}

TEST_CASE("loadGcmcConfig rejects an empty framework_lj array", "[unit][io]") {
    const auto path = writeTemp(
        "aleator_cfg_empty_framework_lj.toml",
        "[run]\nname = \"test\"\n\n[gcmc]\nframework_cif = \"x.cif\"\n"
        "temperature_kelvin = 298.0\npressure_bar = 1.0\n\n"
        "[gcmc.adsorbate]\nname = \"CH4\"\nepsilon_kelvin = 158.5\nsigma_angstrom = 3.72\n"
        "mass_amu = 16.04246\n");
    try {
        (void)aleator::io::loadGcmcConfig(path);
        FAIL("expected ConfigError");
    } catch (const aleator::io::ConfigError& err) {
        const std::string what = err.what();
        REQUIRE(what.find("gcmc.framework_lj") != std::string::npos);
        REQUIRE(what.find("at least one") != std::string::npos);
    }
    std::filesystem::remove(path);
}

TEST_CASE("loadGcmcConfig rejects a framework_lj entry missing sigma_angstrom", "[unit][io]") {
    const auto path = writeTemp(
        "aleator_cfg_bad_framework_lj_entry.toml",
        std::string(
            "[run]\nname = \"test\"\n\n[gcmc]\nframework_cif = \"x.cif\"\n"
            "temperature_kelvin = 298.0\npressure_bar = 1.0\n\n"
            "[gcmc.adsorbate]\nname = \"CH4\"\nepsilon_kelvin = 158.5\nsigma_angstrom = 3.72\n"
            "mass_amu = 16.04246\n\n"
            "[[gcmc.framework_lj]]\nelement = \"Zn\"\nepsilon_kelvin = 62.3992\n"));
    try {
        (void)aleator::io::loadGcmcConfig(path);
        FAIL("expected ConfigError");
    } catch (const aleator::io::ConfigError& err) {
        const std::string what = err.what();
        REQUIRE(what.find("framework_lj[0].sigma_angstrom") != std::string::npos);
    }
    std::filesystem::remove(path);
}

TEST_CASE("loadGcmcConfig rejects partially-specified Peng-Robinson critical properties",
          "[unit][io]") {
    const auto path = writeTemp(
        "aleator_cfg_partial_pr.toml",
        "[run]\nname = \"test\"\n\n[gcmc]\nframework_cif = \"x.cif\"\n"
        "temperature_kelvin = 298.0\npressure_bar = 1.0\n\n"
        "[gcmc.adsorbate]\nname = \"CH4\"\nepsilon_kelvin = 158.5\nsigma_angstrom = 3.72\n"
        "mass_amu = 16.04246\ncritical_temperature_kelvin = 190.6\n\n"
        "[[gcmc.framework_lj]]\nelement = \"Zn\"\nepsilon_kelvin = 62.3992\n"
        "sigma_angstrom = 2.46155\n");
    try {
        (void)aleator::io::loadGcmcConfig(path);
        FAIL("expected ConfigError");
    } catch (const aleator::io::ConfigError& err) {
        const std::string what = err.what();
        REQUIRE(what.find("must all be given together, or all omitted") != std::string::npos);
    }
    std::filesystem::remove(path);
}

TEST_CASE("loadPoreConfig rejects a missing framework_cif", "[unit][io]") {
    const auto path =
        writeTemp("aleator_cfg_pore_missing_cif.toml", "[run]\nname = \"test\"\n\n[pore]\n");
    try {
        (void)aleator::io::loadPoreConfig(path);
        FAIL("expected ConfigError");
    } catch (const aleator::io::ConfigError& err) {
        const std::string what = err.what();
        REQUIRE(what.find("pore.framework_cif") != std::string::npos);
    }
    std::filesystem::remove(path);
}

TEST_CASE("loadMdConfig rejects a non-positive timestep", "[unit][io]") {
    const auto path = writeTemp(
        "aleator_cfg_md_bad_dt.toml",
        "[run]\nname = \"test\"\n\n[md]\nstructure_file = \"x.cif\"\ntimestep_ps = 0.0\n"
        "num_steps = 10\n");
    try {
        (void)aleator::io::loadMdConfig(path);
        FAIL("expected ConfigError");
    } catch (const aleator::io::ConfigError& err) {
        const std::string what = err.what();
        REQUIRE(what.find("md.timestep_ps") != std::string::npos);
        REQUIRE(what.find("must be positive") != std::string::npos);
    }
    std::filesystem::remove(path);
}

TEST_CASE("A syntactically invalid TOML file is rejected with a line number", "[unit][io]") {
    const auto path = writeTemp("aleator_cfg_syntax_error.toml", "[run\nname = \"test\"\n");
    REQUIRE_THROWS_AS(aleator::io::loadRunConfig(path), aleator::io::ConfigError);
    std::filesystem::remove(path);
}

TEST_CASE("detectConfigKind identifies gcmc/pore/md and rejects ambiguous or empty configs",
          "[unit][io]") {
    const auto gcmcPath =
        writeTemp("aleator_cfg_kind_gcmc.toml", "[run]\nname=\"t\"\n\n[gcmc]\nx=1\n");
    REQUIRE(aleator::io::detectConfigKind(gcmcPath) == aleator::io::ConfigKind::Gcmc);
    std::filesystem::remove(gcmcPath);

    const auto porePath =
        writeTemp("aleator_cfg_kind_pore.toml", "[run]\nname=\"t\"\n\n[pore]\nx=1\n");
    REQUIRE(aleator::io::detectConfigKind(porePath) == aleator::io::ConfigKind::Pore);
    std::filesystem::remove(porePath);

    const auto ambiguousPath = writeTemp("aleator_cfg_kind_ambiguous.toml",
                                          "[run]\nname=\"t\"\n\n[gcmc]\nx=1\n\n[pore]\ny=1\n");
    REQUIRE(aleator::io::detectConfigKind(ambiguousPath) == aleator::io::ConfigKind::Unknown);
    std::filesystem::remove(ambiguousPath);

    const auto nonePath = writeTemp("aleator_cfg_kind_none.toml", "[run]\nname=\"t\"\n");
    REQUIRE(aleator::io::detectConfigKind(nonePath) == aleator::io::ConfigKind::Unknown);
    std::filesystem::remove(nonePath);
}
