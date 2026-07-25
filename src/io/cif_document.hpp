#pragma once

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace aleator::io {

/// One `loop_` construct: a set of columns (tag names, including the
/// leading underscore, lowercased) and the rows of raw string values under
/// them. Values are exactly as tokenized from the file (quotes stripped,
/// whitespace-delimited) — numeric parsing/esd-stripping happens later, at
/// the point each field is actually consumed.
struct CifLoop {
    std::vector<std::string> columns;
    std::vector<std::vector<std::string>> rows;

    /// Index of `tag` within `columns`, or -1 if not present in this loop.
    [[nodiscard]] int columnIndex(const std::string& tag) const;

    /// True if every tag in `requiredColumns` appears in this loop.
    [[nodiscard]] bool hasColumns(const std::vector<std::string>& requiredColumns) const;
};

/// A parsed CIF data block: the single-value tags (`_tag value`) and the
/// `loop_` tables it contains. Only the first `data_` block in a file is
/// represented — CoRE MOF / IZA structure files are single-structure CIFs,
/// and supporting multiple concatenated data blocks is out of scope here.
struct CifDocument {
    std::string blockName;
    std::unordered_map<std::string, std::string> tags;
    std::vector<CifLoop> loops;

    /// Raw string value for a single (non-loop) tag, or nullopt if absent.
    [[nodiscard]] std::optional<std::string> findTag(const std::string& tag) const;

    /// The first loop containing every tag in `requiredColumns`, or nullptr
    /// if none does.
    [[nodiscard]] const CifLoop* findLoopWithColumns(
        const std::vector<std::string>& requiredColumns) const;
};

/// Parses raw CIF text into a CifDocument. Throws CifParseError, naming the
/// offending line number where the parser has one, on: unterminated quoted
/// strings or `;`-delimited multi-line text fields, `loop_` rows whose
/// token count doesn't divide evenly into the declared column count, or
/// unrecognized content at the top level of the file.
[[nodiscard]] CifDocument parseCifDocument(const std::string& text);

/// Strips a trailing standard-uncertainty suffix in parentheses (e.g.
/// "11.9190(0)" -> "11.9190") and parses the remainder as a double. Throws
/// CifParseError if the remainder isn't a valid number.
[[nodiscard]] double parseCifNumber(const std::string& token);

/// As parseCifNumber, but returns nullopt (rather than throwing) for CIF's
/// placeholder values "?" (unknown) and "." (inapplicable), which are legal
/// wherever a field is optional.
[[nodiscard]] std::optional<double> tryParseCifNumber(const std::string& token);

} // namespace aleator::io
