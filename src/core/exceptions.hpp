#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

namespace aleator {

/// Thrown by any declared-but-unimplemented entry point (CLAUDE.md invariant
/// #7). Never a silent no-op and never a hardcoded plausible return value —
/// if a code path can be reached before its physics is written, it must
/// throw this instead of guessing.
class NotImplemented : public std::logic_error {
public:
    explicit NotImplemented(std::string_view what)
        : std::logic_error("not implemented: " + std::string(what)) {}
};

} // namespace aleator
