#pragma once
#include <algorithm>
#include <cmath>
#include <Eigen/Core>

namespace Euclid {

//==============================================================================
// Tolerance
//==============================================================================
// Provides scale-aware tolerances for geometric computations.
// Used by primitives, intersection/trimming logic, and topology.
struct Tolerance {
    double absTol;     // Absolute tolerance
    double relTol;     // Relative tolerance (scaled by object magnitude)
    double paramTol;   // Parametric tolerance for curves/surfaces
    double evalFactor; // Scaling factor applied to all evaluations

    // Constructor with default values
    Tolerance(double absT = 1e-6, double relT = 1e-9,
              double paramT = 1e-9, double evalF = 1.0)
        : absTol(absT), relTol(relT), paramTol(paramT), evalFactor(evalF) {}

    // Returns the scale-aware epsilon for numeric comparisons, with a minimum epsilon floor
    double evaluateEpsilon(double scale) const {
        constexpr double minEpsilon = 1e-12; // Minimum allowed epsilon to prevent underflow
        double scaledTol = evalFactor * std::max(absTol, relTol * scale);
        return std::max(scaledTol, minEpsilon);
    }
};

//==============================================================================
// Helper Predicates
//==============================================================================

// Checks if two scalar values are equal within the tolerance
inline bool equalWithinTolerance(double a, double b, const Tolerance& tol, double scale) {
    return std::abs(a - b) <= tol.evaluateEpsilon(scale);
}

// Overload of equalWithinTolerance for vector/matrix comparison
template <typename Derived>
inline bool equalWithinTolerance(
    const Eigen::MatrixBase<Derived>& a,
    const Eigen::MatrixBase<Derived>& b,
    const Tolerance& tol = Tolerance())
{
    if (a.size() != b.size()) {
        return false;
    }
    double scale = std::max(a.norm(), b.norm());
    for (int i = 0; i < a.size(); ++i) {
        if (!equalWithinTolerance(a[i], b[i], tol, scale)) {
            return false;
        }
    }
    return true;
}

// Checks if the first value is definitely less than the second, considering tolerance
inline bool lessThanTolerance(double a, double b, const Tolerance& tol, double scale) {
    return (a - b) < -tol.evaluateEpsilon(scale);
}

// Checks if the first value is definitely greater than the second, considering tolerance
inline bool greaterThanTolerance(double a, double b, const Tolerance& tol, double scale) {
    return (a - b) > tol.evaluateEpsilon(scale);
}

} // namespace Euclid
