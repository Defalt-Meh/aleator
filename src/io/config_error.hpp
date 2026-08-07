#pragma once

#include <stdexcept>
#include <string>

namespace aleator::io {

/// Thrown for any problem with a run-configuration TOML file: a syntax
/// error, a missing required key, or a key present with the wrong type or
/// an out-of-range value. Always names the offending key path and (when
/// toml++ can supply one — a syntax error, or a key that exists but has the
/// wrong type) the exact line and column, so a malformed config is caught
/// before a run starts rather than partway through it.
class ConfigError : public std::runtime_error {
public:
    explicit ConfigError(const std::string& message) : std::runtime_error(message) {}

    /// `line`/`column` point at the exact key when toml++ can supply one
    /// (the key exists but has the wrong type or value); for a genuinely
    /// missing key, pass the *enclosing table's* position instead — the
    /// closest honest approximation toml++ can give, and still far more
    /// useful than no location at all — and word `problem` accordingly
    /// (callers do; see loadGcmcConfig/loadPoreConfig).
    ConfigError(const std::string& filePath, const std::string& keyPath,
                const std::string& problem, std::uint32_t line, std::uint32_t column)
        : std::runtime_error(filePath + ":" + std::to_string(line) + ":" + std::to_string(column) +
                              ": [" + keyPath + "] " + problem) {}
};

} // namespace aleator::io
