#include "io/cif_document.hpp"

#include <algorithm>
#include <cctype>

#include "io/cif_error.hpp"

namespace aleator::io {

namespace {

char toLowerChar(char c) { return static_cast<char>(std::tolower(static_cast<unsigned char>(c))); }

std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), toLowerChar);
    return s;
}

bool isSpaceChar(char c) { return c == ' ' || c == '\t' || c == '\r'; }

bool isBlankLine(const std::string& line) {
    return std::all_of(line.begin(), line.end(), [](char c) {
        return std::isspace(static_cast<unsigned char>(c)) != 0;
    });
}

bool startsWithCaseInsensitive(const std::string& s, const std::string& prefix) {
    if (s.size() < prefix.size()) {
        return false;
    }
    for (std::size_t k = 0; k < prefix.size(); ++k) {
        if (toLowerChar(s[k]) != toLowerChar(prefix[k])) {
            return false;
        }
    }
    return true;
}

/// True if `line`, once leading whitespace is skipped, begins a new
/// top-level structural construct (tag, loop_, or data_) — used to decide
/// whether a loop_'s row-accumulation phase should stop.
bool looksLikeNewConstruct(const std::string& line) {
    std::size_t i = 0;
    while (i < line.size() && isSpaceChar(line[i])) {
        ++i;
    }
    if (i >= line.size() || line[i] == '#') {
        return false;
    }
    const std::string rest = line.substr(i);
    return rest[0] == '_' || startsWithCaseInsensitive(rest, "loop_") ||
           startsWithCaseInsensitive(rest, "data_");
}

std::vector<std::string> splitLines(const std::string& text) {
    std::vector<std::string> lines;
    std::string current;
    for (char c : text) {
        if (c == '\n') {
            lines.push_back(current);
            current.clear();
        } else {
            current += c;
        }
    }
    lines.push_back(current);
    for (auto& line : lines) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
    }
    return lines;
}

/// Tokenizes one line into whitespace-separated fields, honoring CIF quote
/// rules: a quoted field runs from an opening `'`/`"` to the next matching
/// quote character that is itself followed by whitespace or end-of-line
/// (so an apostrophe inside a single-quoted field, e.g. a space group
/// symbol, doesn't end the field prematurely). Everything from an
/// unquoted `#` to the end of the line is a comment.
std::vector<std::string> tokenizeLine(const std::string& line, int lineNumber) {
    std::vector<std::string> tokens;
    const std::size_t n = line.size();
    std::size_t i = 0;
    while (i < n) {
        while (i < n && isSpaceChar(line[i])) {
            ++i;
        }
        if (i >= n || line[i] == '#') {
            break;
        }
        if (line[i] == '\'' || line[i] == '"') {
            const char quote = line[i];
            ++i;
            const std::size_t start = i;
            while (i < n && !(line[i] == quote && (i + 1 >= n || isSpaceChar(line[i + 1])))) {
                ++i;
            }
            if (i >= n) {
                throw CifParseError("unterminated quoted string", lineNumber);
            }
            tokens.push_back(line.substr(start, i - start));
            ++i; // skip closing quote
        } else {
            const std::size_t start = i;
            while (i < n && !isSpaceChar(line[i]) && line[i] != '#') {
                ++i;
            }
            tokens.push_back(line.substr(start, i - start));
        }
    }
    return tokens;
}

} // namespace

int CifLoop::columnIndex(const std::string& tag) const {
    const auto it = std::find(columns.begin(), columns.end(), tag);
    return it == columns.end() ? -1 : static_cast<int>(it - columns.begin());
}

bool CifLoop::hasColumns(const std::vector<std::string>& requiredColumns) const {
    return std::all_of(requiredColumns.begin(), requiredColumns.end(),
                        [this](const std::string& tag) { return columnIndex(tag) >= 0; });
}

std::optional<std::string> CifDocument::findTag(const std::string& tag) const {
    const auto it = tags.find(toLower(tag));
    if (it == tags.end()) {
        return std::nullopt;
    }
    return it->second;
}

const CifLoop* CifDocument::findLoopWithColumns(
    const std::vector<std::string>& requiredColumns) const {
    for (const auto& loop : loops) {
        if (loop.hasColumns(requiredColumns)) {
            return &loop;
        }
    }
    return nullptr;
}

double parseCifNumber(const std::string& token) {
    const auto value = tryParseCifNumber(token);
    if (!value.has_value()) {
        throw CifParseError("expected a number but found \"" + token + "\"");
    }
    return *value;
}

std::optional<double> tryParseCifNumber(const std::string& token) {
    if (token == "?" || token == ".") {
        return std::nullopt;
    }
    // Strip a trailing standard-uncertainty suffix, e.g. "11.9190(0)" -> "11.9190".
    std::string cleaned = token;
    if (const auto paren = cleaned.find('('); paren != std::string::npos) {
        cleaned = cleaned.substr(0, paren);
    }
    if (cleaned.empty()) {
        return std::nullopt;
    }
    try {
        std::size_t consumed = 0;
        const double value = std::stod(cleaned, &consumed);
        if (consumed != cleaned.size()) {
            return std::nullopt;
        }
        return value;
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

CifDocument parseCifDocument(const std::string& text) {
    const std::vector<std::string> lines = splitLines(text);
    const std::size_t n = lines.size();

    auto stripLeadingSpace = [](const std::string& s) -> std::string {
        std::size_t k = 0;
        while (k < s.size() && isSpaceChar(s[k])) {
            ++k;
        }
        return s.substr(k);
    };

    CifDocument doc;
    bool sawDataBlock = false;
    std::size_t i = 0;

    while (i < n) {
        const std::string& line = lines[i];
        const int lineNumber = static_cast<int>(i + 1);

        if (isBlankLine(line)) {
            ++i;
            continue;
        }
        const std::string content = stripLeadingSpace(line);
        if (content.front() == '#') {
            ++i;
            continue;
        }

        if (startsWithCaseInsensitive(content, "data_")) {
            if (sawDataBlock) {
                // Only the first data_ block is supported; stop here.
                break;
            }
            const auto tokens = tokenizeLine(line, lineNumber);
            doc.blockName = tokens.empty() ? std::string{} : tokens.front().substr(5);
            sawDataBlock = true;
            ++i;
            continue;
        }

        if (startsWithCaseInsensitive(content, "loop_")) {
            ++i;
            CifLoop loop;
            // Column declarations: consecutive lines starting with '_'.
            while (i < n) {
                if (isBlankLine(lines[i])) {
                    ++i;
                    continue;
                }
                std::size_t k = 0;
                while (k < lines[i].size() && isSpaceChar(lines[i][k])) {
                    ++k;
                }
                if (k < lines[i].size() && lines[i][k] == '#') {
                    ++i;
                    continue;
                }
                if (k >= lines[i].size() || lines[i][k] != '_') {
                    break;
                }
                const auto tokens = tokenizeLine(lines[i], static_cast<int>(i + 1));
                if (tokens.empty()) {
                    throw CifParseError("empty loop_ column declaration", static_cast<int>(i + 1));
                }
                loop.columns.push_back(toLower(tokens.front()));
                ++i;
            }
            if (loop.columns.empty()) {
                throw CifParseError("loop_ with no column declarations", lineNumber);
            }

            // Data rows: accumulate tokens across lines until each group of
            // columns.size() tokens forms a complete row.
            std::vector<std::string> pending;
            std::size_t lastRowLine = static_cast<std::size_t>(lineNumber);
            while (i < n) {
                if (isBlankLine(lines[i])) {
                    ++i;
                    continue;
                }
                if (pending.empty() && looksLikeNewConstruct(lines[i])) {
                    break;
                }
                const auto tokens = tokenizeLine(lines[i], static_cast<int>(i + 1));
                lastRowLine = i + 1;
                for (const auto& token : tokens) {
                    pending.push_back(token);
                    if (pending.size() == loop.columns.size()) {
                        loop.rows.push_back(pending);
                        pending.clear();
                    }
                }
                ++i;
            }
            if (!pending.empty()) {
                throw CifParseError(
                    "loop_ row has " + std::to_string(pending.size()) + " value(s) but " +
                        std::to_string(loop.columns.size()) + " columns (" +
                        std::to_string(loop.columns.size() - pending.size()) + " short)",
                    static_cast<int>(lastRowLine));
            }
            doc.loops.push_back(std::move(loop));
            continue;
        }

        if (content.front() == '_') {
            const auto tokens = tokenizeLine(line, lineNumber);
            const std::string tag = toLower(tokens.front());
            std::string value;
            if (tokens.size() >= 2) {
                value = tokens[1];
                ++i;
            } else {
                ++i;
                if (i < n && !lines[i].empty() && lines[i][0] == ';') {
                    const std::size_t openLine = i;
                    std::string firstLine = lines[i].substr(1);
                    std::vector<std::string> textLines;
                    if (!firstLine.empty()) {
                        textLines.push_back(firstLine);
                    }
                    ++i;
                    bool closed = false;
                    while (i < n) {
                        if (!lines[i].empty() && lines[i][0] == ';') {
                            closed = true;
                            ++i;
                            break;
                        }
                        textLines.push_back(lines[i]);
                        ++i;
                    }
                    if (!closed) {
                        throw CifParseError("unterminated ';' text field for tag " + tag,
                                             static_cast<int>(openLine + 1));
                    }
                    std::string joined;
                    for (std::size_t k = 0; k < textLines.size(); ++k) {
                        if (k > 0) {
                            joined += '\n';
                        }
                        joined += textLines[k];
                    }
                    value = joined;
                } else if (i < n) {
                    const auto valueTokens = tokenizeLine(lines[i], static_cast<int>(i + 1));
                    if (valueTokens.empty()) {
                        throw CifParseError("missing value for tag " + tag, lineNumber);
                    }
                    value = valueTokens.front();
                    ++i;
                } else {
                    throw CifParseError("missing value for tag " + tag + " at end of file",
                                         lineNumber);
                }
            }
            doc.tags[tag] = value;
            continue;
        }

        throw CifParseError("unrecognized content outside of a tag, loop_, or data_ block",
                             lineNumber);
    }

    if (!sawDataBlock) {
        throw CifParseError("no data_ block found in file");
    }

    return doc;
}

} // namespace aleator::io
