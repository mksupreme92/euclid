#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <Eigen/Dense>

namespace euclid::tests {

// Print a test result with emoji
inline void printTest(const std::string& name, bool result) {
    std::cout << (result ? "✅ " : "❌ ") << name << "\n";
}

// Compare a measured time to a benchmark time, show both times and % difference, emoji shows pass/fail
inline void benchmarkTest(const std::string& name, double measuredTime, double benchmarkTime, double maxDivergencePercent) {
    double divergence = (measuredTime / benchmarkTime - 1.0) * 100.0;
    bool pass = std::abs(divergence) <= maxDivergencePercent; // abs so negatives don't fail

    auto round3 = [](double x) { return std::round(x * 1000.0) / 1000.0; };
    
    std::cout << name << "\n"
              << "measured = " << round3(measuredTime) << "s \n"
              << "benchmark = " << round3(benchmarkTime) << "s \n"
              << "diff = " << round3(divergence) << "% "
              << (pass ? "✅" : "❌") << "\n\n";
}

// Scalar approximate equality check
template<typename T>
inline bool approxEqual(T a, T b, T tol = static_cast<T>(1e-9)) {
    return std::abs(a - b) <= tol;
}

// Check if two Eigen vectors are approximately equal within a tolerance
template<typename T, int N>
inline bool approxEqualVector(const Eigen::Matrix<T, N, 1>& a, const Eigen::Matrix<T, N, 1>& b, T tol = static_cast<T>(1e-9)) {
    for (int i = 0; i < N; ++i) {
        if (!approxEqual(a(i), b(i), tol)) return false;
    }
    return true;
}

// Compare two Eigen matrices approximately (element-wise)
template<typename T, int Rows, int Cols>
inline bool approxEqualMatrix(const Eigen::Matrix<T, Rows, Cols>& A,
                              const Eigen::Matrix<T, Rows, Cols>& B,
                              T tol = static_cast<T>(1e-9)) {
    for (int i = 0; i < Rows; ++i)
        for (int j = 0; j < Cols; ++j)
            if (!approxEqual(A(i,j), B(i,j), tol)) return false;
    return true;
}

} // namespace euclid::tests
