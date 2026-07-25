#include "io/element_table.hpp"

#include <array>
#include <cctype>
#include <string_view>

namespace aleator::io {

namespace {

struct ElementEntry {
    std::string_view symbol; // canonical case, e.g. "Si", "Cu"
    double mass;              // standard atomic weight, amu (IUPAC conventional values)
};

// Covers the elements realistically encountered in MOF/zeolite structure
// files (framework atoms, extra-framework cations, common metal nodes).
// Not exhaustive of the periodic table; unrecognized symbols return
// nullopt rather than a guess (CLAUDE.md #2.4).
constexpr std::array<ElementEntry, 83> kElements{{
    {"H", 1.008},      {"He", 4.003},     {"Li", 6.94},      {"Be", 9.012},
    {"B", 10.81},      {"C", 12.011},     {"N", 14.007},     {"O", 15.999},
    {"F", 18.998},     {"Ne", 20.180},    {"Na", 22.990},    {"Mg", 24.305},
    {"Al", 26.982},    {"Si", 28.085},    {"P", 30.974},     {"S", 32.06},
    {"Cl", 35.45},     {"Ar", 39.948},    {"K", 39.098},     {"Ca", 40.078},
    {"Sc", 44.956},    {"Ti", 47.867},    {"V", 50.942},     {"Cr", 51.996},
    {"Mn", 54.938},    {"Fe", 55.845},    {"Co", 58.933},    {"Ni", 58.693},
    {"Cu", 63.546},    {"Zn", 65.38},     {"Ga", 69.723},    {"Ge", 72.630},
    {"As", 74.922},    {"Se", 78.971},    {"Br", 79.904},    {"Kr", 83.798},
    {"Rb", 85.468},    {"Sr", 87.62},     {"Y", 88.906},     {"Zr", 91.224},
    {"Nb", 92.906},    {"Mo", 95.95},     {"Tc", 98.0},      {"Ru", 101.07},
    {"Rh", 102.906},   {"Pd", 106.42},    {"Ag", 107.868},   {"Cd", 112.414},
    {"In", 114.818},   {"Sn", 118.710},   {"Sb", 121.760},   {"Te", 127.60},
    {"I", 126.904},    {"Xe", 131.293},   {"Cs", 132.905},   {"Ba", 137.327},
    {"La", 138.905},   {"Ce", 140.116},   {"Pr", 140.908},   {"Nd", 144.242},
    {"Pm", 145.0},     {"Sm", 150.36},    {"Eu", 151.964},   {"Gd", 157.25},
    {"Tb", 158.925},   {"Dy", 162.500},   {"Ho", 164.930},   {"Er", 167.259},
    {"Tm", 168.934},   {"Yb", 173.045},   {"Lu", 174.967},   {"Hf", 178.49},
    {"Ta", 180.948},   {"W", 183.84},     {"Re", 186.207},   {"Os", 190.23},
    {"Ir", 192.217},   {"Pt", 195.084},   {"Au", 196.967},   {"Hg", 200.592},
    {"Tl", 204.38},    {"Pb", 207.2},     {"Bi", 208.980},
}};

bool equalsIgnoreCase(std::string_view a, std::string_view b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(a[i])) !=
            std::tolower(static_cast<unsigned char>(b[i]))) {
            return false;
        }
    }
    return true;
}

const ElementEntry* findByNameIgnoreCase(std::string_view candidate) {
    for (const auto& entry : kElements) {
        if (equalsIgnoreCase(entry.symbol, candidate)) {
            return &entry;
        }
    }
    return nullptr;
}

} // namespace

std::optional<double> atomicMass(const std::string& symbol) {
    // Case-sensitive per the public contract, but element symbols are
    // conventionally exact case, so an exact (not case-insensitive) match
    // catches genuine mistakes like "cu" or "CU" rather than silently
    // accepting them.
    for (const auto& entry : kElements) {
        if (entry.symbol == symbol) {
            return entry.mass;
        }
    }
    return std::nullopt;
}

std::optional<std::string> elementFromLabel(const std::string& label) {
    std::size_t letterCount = 0;
    while (letterCount < label.size() &&
           std::isalpha(static_cast<unsigned char>(label[letterCount]))) {
        ++letterCount;
    }
    if (letterCount == 0) {
        return std::nullopt;
    }

    if (letterCount >= 2) {
        if (const auto* twoLetter = findByNameIgnoreCase(std::string_view(label).substr(0, 2))) {
            return std::string(twoLetter->symbol);
        }
    }
    if (const auto* oneLetter = findByNameIgnoreCase(std::string_view(label).substr(0, 1))) {
        return std::string(oneLetter->symbol);
    }
    return std::nullopt;
}

} // namespace aleator::io
