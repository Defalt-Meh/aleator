#include "io/symmetry_operation.hpp"

#include <cctype>
#include <cstdlib>

#include "io/cif_error.hpp"

namespace aleator::io {

namespace {

struct AffineTerm {
    double coeffX = 0.0;
    double coeffY = 0.0;
    double coeffZ = 0.0;
    double translation = 0.0;
};

std::string stripWhitespace(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (!std::isspace(static_cast<unsigned char>(c))) {
            out += c;
        }
    }
    return out;
}

/// Parses one component (e.g. "1/2-x+y" or "z+1/3") of a symmetry
/// operation string into coefficients on x, y, z plus a constant
/// translation. `original` is the whole operation string, used only to
/// give a useful error message.
AffineTerm parseComponent(const std::string& component, const std::string& original) {
    AffineTerm term;
    std::size_t i = 0;
    const std::size_t n = component.size();
    if (n == 0) {
        throw CifParseError("empty component in symmetry operation \"" + original + "\"");
    }

    while (i < n) {
        double sign = 1.0;
        if (component[i] == '+') {
            ++i;
        } else if (component[i] == '-') {
            sign = -1.0;
            ++i;
        }
        if (i >= n) {
            throw CifParseError("trailing sign with no term in symmetry operation \"" + original +
                                 "\"");
        }

        double numeric = 1.0;
        bool hasNumeric = false;
        if (std::isdigit(static_cast<unsigned char>(component[i])) || component[i] == '.') {
            const std::size_t start = i;
            while (i < n &&
                   (std::isdigit(static_cast<unsigned char>(component[i])) || component[i] == '.')) {
                ++i;
            }
            const double numerator = std::strtod(component.substr(start, i - start).c_str(), nullptr);
            double denominator = 1.0;
            if (i < n && component[i] == '/') {
                ++i;
                const std::size_t denomStart = i;
                while (i < n && std::isdigit(static_cast<unsigned char>(component[i]))) {
                    ++i;
                }
                if (i == denomStart) {
                    throw CifParseError("malformed fraction in symmetry operation \"" + original +
                                         "\"");
                }
                denominator = std::strtod(component.substr(denomStart, i - denomStart).c_str(), nullptr);
            }
            numeric = numerator / denominator;
            hasNumeric = true;
        }

        if (i < n && (component[i] == 'x' || component[i] == 'X' || component[i] == 'y' ||
                      component[i] == 'Y' || component[i] == 'z' || component[i] == 'Z')) {
            const char variable = static_cast<char>(std::tolower(static_cast<unsigned char>(component[i])));
            const double coeff = sign * (hasNumeric ? numeric : 1.0);
            if (variable == 'x') {
                term.coeffX += coeff;
            } else if (variable == 'y') {
                term.coeffY += coeff;
            } else {
                term.coeffZ += coeff;
            }
            ++i;
        } else if (hasNumeric) {
            term.translation += sign * numeric;
        } else {
            throw CifParseError("expected x, y, z, or a number in symmetry operation \"" +
                                 original + "\"");
        }
    }
    return term;
}

} // namespace

std::array<double, 3> SymmetryOperation::apply(const std::array<double, 3>& fractional) const {
    return {
        matrix[0][0] * fractional[0] + matrix[0][1] * fractional[1] + matrix[0][2] * fractional[2] +
            translation[0],
        matrix[1][0] * fractional[0] + matrix[1][1] * fractional[1] + matrix[1][2] * fractional[2] +
            translation[1],
        matrix[2][0] * fractional[0] + matrix[2][1] * fractional[1] + matrix[2][2] * fractional[2] +
            translation[2],
    };
}

SymmetryOperation parseSymmetryOperation(const std::string& text) {
    const std::string compact = stripWhitespace(text);

    std::array<std::string, 3> components;
    std::size_t componentIndex = 0;
    std::string current;
    for (char c : compact) {
        if (c == ',') {
            if (componentIndex >= 3) {
                throw CifParseError("symmetry operation \"" + text +
                                     "\" has more than 3 comma-separated components");
            }
            components[componentIndex++] = current;
            current.clear();
        } else {
            current += c;
        }
    }
    if (componentIndex >= 3) {
        throw CifParseError("symmetry operation \"" + text +
                             "\" has more than 3 comma-separated components");
    }
    components[componentIndex++] = current;
    if (componentIndex != 3) {
        throw CifParseError("symmetry operation \"" + text + "\" must have exactly 3 components, found " +
                             std::to_string(componentIndex));
    }

    SymmetryOperation op;
    for (int row = 0; row < 3; ++row) {
        const AffineTerm term = parseComponent(components[static_cast<std::size_t>(row)], text);
        op.matrix[static_cast<std::size_t>(row)] = {term.coeffX, term.coeffY, term.coeffZ};
        op.translation[static_cast<std::size_t>(row)] = term.translation;
    }
    return op;
}

} // namespace aleator::io
