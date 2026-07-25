#pragma once

#include <stdexcept>
#include <string>

namespace aleator::io {

/// Thrown for any structurally malformed CIF input: unterminated quotes or
/// semicolon text fields, loop_ row/column mismatches, missing required
/// tags, or unparseable numeric fields. Always a clear, specific message
/// (naming the offending line where the parser has one) — never a silent
/// wrong parse or a crash.
class CifParseError : public std::runtime_error {
public:
    explicit CifParseError(const std::string& message) : std::runtime_error(message) {}

    CifParseError(const std::string& message, int lineNumber)
        : std::runtime_error(message + " (line " + std::to_string(lineNumber) + ")") {}
};

} // namespace aleator::io
