// nanobind extension module. Exposes exactly one trivial function
// (vector_sum, wrapping the Highway dynamic-dispatch reference kernel from
// core/math/simd_kernels.hpp) plus the package version — nothing physics-
// related is bound yet.

#include <cstdint>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include <aleator/version.hpp>

#include "core/math/simd_kernels.hpp"

namespace nb = nanobind;

namespace {

std::int64_t vectorSum(const std::vector<std::int64_t>& values) {
    return aleator::core::vectorSum(values.data(), values.size());
}

} // namespace

NB_MODULE(_aleator, m) {
    m.doc() = "Aleator: Monte Carlo, molecular dynamics, and porous-material geometry "
              "analysis on a shared periodic-system core. This build exposes scaffolding "
              "only — no physics is implemented yet.";
    m.attr("__version__") = aleator::kVersion;
    m.def("vector_sum", &vectorSum, nb::arg("values"),
          "Sum a list of integers via the Highway runtime-SIMD-dispatch reference kernel.");
}
