#pragma once

#include <optional>
#include <string>

namespace aleator::io {

/// Standard atomic weight (amu) for element `symbol` (e.g. "H", "Si",
/// "Cu"), case-sensitive per IUPAC convention (first letter uppercase,
/// second lowercase). Returns nullopt for a symbol this table doesn't
/// recognize rather than guessing a value — CLAUDE.md #2.4: unit/identity
/// mistakes at an I/O boundary are silent and catastrophic.
[[nodiscard]] std::optional<double> atomicMass(const std::string& symbol);

/// Best-effort element symbol extracted from a CIF atom-site label such as
/// "Si1", "O12A", or "Cu2+": strips trailing digits, '+'/'-' charge marks,
/// and disorder-suffix letters, then matches the remaining leading letters
/// against known element symbols (two-letter symbols are tried before
/// one-letter ones, e.g. "Cl1" -> "Cl", not "C"). Returns nullopt if no
/// known element symbol is a prefix of the label.
[[nodiscard]] std::optional<std::string> elementFromLabel(const std::string& label);

} // namespace aleator::io
