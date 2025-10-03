#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <Eigen/Dense>
//#include "geometry/tolerance.hpp"

namespace Euclid::Tests {

// -----------------------
// Test printing helpers
// -----------------------
inline void printTest(const std::string& name, bool result) {
    std::cout << (result ? "✅ " : "❌ ") << name << "\n";
}

inline void benchmarkTest(const std::string& name, double measuredTime, double benchmarkTime, double maxDivergencePercent) {
    double divergence = (measuredTime / benchmarkTime - 1.0) * 100.0;
    bool pass = std::abs(divergence) <= maxDivergencePercent;
    auto round3 = [](double x) { return std::round(x * 1000.0) / 1000.0; };

    std::cout << name << "\n"
              << "measured = " << round3(measuredTime) << "s \n"
              << "benchmark = " << round3(benchmarkTime) << "s \n"
              << "diff = " << round3(divergence) << "% "
              << (pass ? "✅" : "❌") << "\n\n";
}


} // namespace Euclid::Tests
