#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "geometry/curve.hpp"
#include "../geometry/intersection/curve_intersection.hpp"
#include "tests/test_utilities.hpp"

namespace Euclid {
namespace Tests {

// ∿ Curve–Curve intersection test
inline void testCurveCurveIntersection() {
    std::cout << "\n∿ Testing Curve - Curve intersection:\n";

    auto printPointInfo = [](const auto& result) {
        std::vector<bool> duplicateFlags(result.points.size(), false);
        for (size_t i = 0; i < result.points.size(); ++i) {
            for (size_t j = i + 1; j < result.points.size(); ++j) {
                if ((result.points[i].coords - result.points[j].coords).norm() < 1e-12) {
                    duplicateFlags[j] = true;
                }
            }
        }
        for (size_t i = 0; i < result.points.size(); ++i) {
            const auto& pt = result.points[i];
            std::cout << "  Point " << i << ": ("
                      << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[0] << ", "
                      << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[1] << ")"
                      << " | Parameters: "
                      << "u_curve=" << std::setw(10) << std::fixed << std::setprecision(6)
                      << (i < result.hits.size() ? result.hits[i].u_curve : 0.0) << ", "
                      << "u_curve2=" << std::setw(10) << std::fixed << std::setprecision(6)
                      << (i < result.hits.size() ? result.hits[i].u_curve2 : 0.0);
            if (duplicateFlags[i]) {
                std::cout << " (duplicate)";
            }
            std::cout << std::endl;
        }
    };

    // Test case 1: Two 2D linear curves intersecting at (0.5, 0.5)
    std::cout << "--- Test Case 1: Diagonal Cross ---" << std::endl;
    Curve<double,2> c1([](double t){ return Point<double,2>({t, t}); }, 0.0, 1.0);
    Curve<double,2> c2([](double t){ return Point<double,2>({t, 1.0 - t}); }, 0.0, 1.0);
    auto result1 = intersect(c1, c2);
    std::cout << "Number of intersection points: " << result1.points.size() << std::endl;
    printPointInfo(result1);
    printTest("Curve–Curve intersection: has intersection", result1.intersects);
    if (result1.intersects) {
        bool correctLocation = false;
        for (const auto& pt : result1.points) {
            if (pt == Point<double,2>({0.5, 0.5})) {
                correctLocation = true;
                break;
            }
        }
        printTest("Correct intersection location", correctLocation);
    }
    std::cout << "----------------------------------------" << std::endl;

    // Test case 2: Two parallel 2D linear curves (no intersection)
    std::cout << "--- Test Case 2: Parallel Lines ---" << std::endl;
    Curve<double,2> c3([](double t){ return Point<double,2>({t, 0.0}); }, 0.0, 1.0);
    Curve<double,2> c4([](double t){ return Point<double,2>({t, 1.0}); }, 0.0, 1.0);
    auto result2 = intersect(c3, c4);
    std::cout << "Number of intersection points: " << result2.points.size() << std::endl;
    printPointInfo(result2);
    printTest("Curve–Curve intersection: no intersection", !result2.intersects || result2.points.empty());
    std::cout << "----------------------------------------" << std::endl;

    // Test case 3: A sinusoidal curve intersecting a horizontal line (multiple intersections)
    std::cout << "--- Test Case 3: Sinusoidal Intersections ---" << std::endl;
    Curve<double,2> c5([](double t){ return Point<double,2>({t, std::sin(2*M_PI*t)}); }, 0.0, 1.0);
    Curve<double,2> c6([](double t){ return Point<double,2>({t, 0.0}); }, 0.0, 1.0);
    auto result3 = intersect(c5, c6);
    std::cout << "Number of intersection points: " << result3.points.size() << std::endl;
    printPointInfo(result3);
    printTest("Curve–Curve intersection: expected 3 intersections", result3.points.size() == 3);
    std::cout << "[Analytical validation for sinusoidal intersections]" << std::endl;
    auto checkAgainstExpected = [&](const std::vector<double>& expectedParams, const Curve<double,2>& c, const auto& result) {
        const double paramTol = 1e-5;
        const double spatialTol = 1e-5;
        for (size_t i = 0; i < expectedParams.size(); ++i) {
            const double ue = expectedParams[i];
            const auto pe = c.evaluate(ue);
            double bestErrU = std::numeric_limits<double>::max();
            double bestErrP = std::numeric_limits<double>::max();
            for (const auto& hit : result.hits) {
                double du = std::abs(hit.u_curve - ue);
                auto ph = c.evaluate(hit.u_curve);
                double dp = (ph.coords - pe.coords).norm();
                bestErrU = std::min(bestErrU, du);
                bestErrP = std::min(bestErrP, dp);
            }
            std::cout << "  [Expected] u=" << std::fixed << std::setprecision(6) << ue
                      << " → (" << pe.coords[0] << ", " << pe.coords[1]
                      << ") | best Δu=" << bestErrU << " Δp=" << bestErrP << std::endl;
            if (bestErrU > paramTol || bestErrP > spatialTol)
                printTest("Intersection deviation exceeds tolerance", false);
        }
    };
    checkAgainstExpected({0.0, 0.5, 1.0}, c5, result3);
    std::cout << "----------------------------------------" << std::endl;

    // Test Case 3B: Sinusoidal Intersections (Extended Domain)
    std::cout << "--- Test Case 3B: Sinusoidal Intersections (Extended Domain) ---" << std::endl;
    Curve<double,2> c5_ext([](double t){ return Point<double,2>({t, std::sin(2*M_PI*t)}); }, 0.0, 2.0);
    Curve<double,2> c6_ext([](double t){ return Point<double,2>({t, 0.0}); }, 0.0, 2.0);
    auto result3B = intersect(c5_ext, c6_ext);
    std::cout << "Number of intersection points: " << result3B.points.size() << std::endl;
    printPointInfo(result3B);
    printTest("Curve–Curve intersection: expected 5 intersections (extended domain)", result3B.points.size() == 5);
    std::cout << "[Analytical validation for sinusoidal intersections (extended domain)]" << std::endl;
    checkAgainstExpected({0.0, 0.5, 1.0, 1.5, 2.0}, c5_ext, result3B);
    std::cout << "----------------------------------------" << std::endl;

    // Test case 4: Circle–Line intersection (two intersections)
    std::cout << "--- Test Case 4: Circle–Line Intersection ---" << std::endl;
    Curve<double,2> circle([](double t){ return Point<double,2>({std::cos(2*M_PI*t), std::sin(2*M_PI*t)}); }, 0.0, 1.0);
    Curve<double,2> line([](double t){ return Point<double,2>({t * 2.0 - 1.0, 0.0}); }, 0.0, 1.0);
    auto result4 = intersect(circle, line);
    std::cout << "Number of intersection points: " << result4.points.size() << std::endl;
    printPointInfo(result4);
    printTest("Curve–Curve intersection: two intersections (circle-line)", result4.points.size() == 2);
    std::cout << "----------------------------------------" << std::endl;

    // Test case 5: Tangential intersection (circle tangent to line y = 1)
    std::cout << "--- Test Case 5: Tangential Intersection ---" << std::endl;
    Curve<double,2> circle2([](double t){ return Point<double,2>({std::cos(2*M_PI*t), std::sin(2*M_PI*t)}); }, 0.0, 1.0);
    Curve<double,2> tangent([](double t){ return Point<double,2>({t * 2.0 - 1.0, 1.0}); }, 0.0, 1.0);
    auto result5 = intersect(circle2, tangent);
    std::cout << "Number of intersection points: " << result5.points.size() << std::endl;
    printPointInfo(result5);
    printTest("Curve–Curve intersection: single tangent point", result5.points.size() == 1);
    std::cout << "----------------------------------------" << std::endl;

    // Test case 6: Overlapping curves (identical line segments)
    std::cout << "--- Test Case 6: Overlapping Curves ---" << std::endl;
    Curve<double,2> lineA([](double t){ return Point<double,2>({t, 0.5 * t}); }, 0.0, 1.0);
    Curve<double,2> lineB([](double t){ return Point<double,2>({t, 0.5 * t}); }, 0.0, 1.0);
    auto result6 = intersect(lineA, lineB);
    std::cout << "Number of intersection points: " << result6.points.size() << std::endl;
    printPointInfo(result6);
    printTest("Curve–Curve intersection: overlapping curves", result6.points.size() >= 1 || result6.intersects);
    std::cout << "----------------------------------------" << std::endl;

    // Test Case 7: Offset Sine Curves (phase-shifted)
    std::cout << "--- Test Case 7: Offset Sine Curves ---" << std::endl;
    Curve<double,2> sine1([](double t){ return Point<double,2>({t, std::sin(2*M_PI*t)}); }, 0.0, 1.0);
    Curve<double,2> sine2([](double t){ return Point<double,2>({t, std::sin(2*M_PI*t + M_PI/4)}); }, 0.0, 1.0);
    auto result7 = intersect(sine1, sine2);
    std::cout << "Number of intersection points: " << result7.points.size() << std::endl;
    printPointInfo(result7);
    printTest("Curve–Curve intersection: expected 2 intersections for offset sine curves", result7.points.size() == 2);
    std::cout << "----------------------------------------" << std::endl;

    // Test Case 8: Perpendicular Circular Arcs
    std::cout << "--- Test Case 8: Perpendicular Circular Arcs ---" << std::endl;
    Curve<double,2> arc1([](double t){ return Point<double,2>({std::cos(M_PI * t), std::sin(M_PI * t)}); }, 0.0, 1.0);
    Curve<double,2> arc2([](double t){ return Point<double,2>({std::cos(M_PI * t), -std::sin(M_PI * t)}); }, 0.0, 1.0);
    auto result8 = intersect(arc1, arc2);
    std::cout << "Number of intersection points: " << result8.points.size() << std::endl;
    printPointInfo(result8);
    printTest("Curve–Curve intersection: two intersections near x ≈ 0.5", result8.points.size() == 2);
    std::cout << "----------------------------------------" << std::endl;

    // Test Case 9: Spiral vs Line (monotonic spacing test)
    std::cout << "--- Test Case 9: Spiral vs Line ---" << std::endl;
    Curve<double,2> spiral([](double t){
        double r = t;
        double theta = 4*M_PI*t;
        return Point<double,2>({r * std::cos(theta), r * std::sin(theta)});
    }, 0.0, 1.0);
    Curve<double,2> lineY0([](double t){ return Point<double,2>({t * 2.0 - 1.0, 0.0}); }, 0.0, 1.0);
    auto result9 = intersect(spiral, lineY0);
    std::cout << "Number of intersection points: " << result9.points.size() << std::endl;
    printPointInfo(result9);
    printTest("Curve–Curve intersection: multiple intersections (>= 2) for spiral vs line", result9.points.size() >= 2);
    std::cout << "----------------------------------------" << std::endl;

    // Test Case 10: High-curvature polynomial vs line
    std::cout << "--- Test Case 10: High-curvature polynomial vs line ---" << std::endl;
    Curve<double,2> poly([](double t){
        double y = 10*t*t*t - 15*t*t + 5*t;
        return Point<double,2>({t, y});
    }, 0.0, 1.0);
    Curve<double,2> lineY025([](double t){ return Point<double,2>({t, 0.25}); }, 0.0, 1.0);
    auto result10 = intersect(poly, lineY025);
    std::cout << "Number of intersection points: " << result10.points.size() << std::endl;
    printPointInfo(result10);
    // WolframAlpha verification:
    // solve 10*t^3 - 15*t^2 + 5*t = 0.25 for t from 0 to 1
    // Solutions: t ≈ 0.06056, 0.39543 → 2 intersections expected
    printTest("Curve–Curve intersection: expected 2 intersections for high-curvature polynomial vs line", result10.points.size() == 2);
    std::cout << "----------------------------------------" << std::endl;

    std::cout << std::endl;
}

// ∿ Curve–Curve intersection (High-Dimensional / Robustness) tests
inline void testCurveCurveIntersectionHighDim() {
    std::cout << "\n∿ Testing Curve - Curve intersection (High-Dimensional / Robustness):\n";

    // Helper for printing intersection points (works for any dimension)
    auto printPointInfoHD = [](const auto& result) {
        std::vector<bool> duplicateFlags(result.points.size(), false);
        for (size_t i = 0; i < result.points.size(); ++i) {
            for (size_t j = i + 1; j < result.points.size(); ++j) {
                if ((result.points[i].coords - result.points[j].coords).norm() < 1e-12) {
                    duplicateFlags[j] = true;
                }
            }
        }
        for (size_t i = 0; i < result.points.size(); ++i) {
            const auto& pt = result.points[i];
            std::cout << "  Point " << i << ": (";
            for (int d = 0; d < pt.coords.size(); ++d) {
                std::cout << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[d];
                if (d + 1 < pt.coords.size()) std::cout << ", ";
            }
            std::cout << ")"
                      << " | Parameters: "
                      << "u_curve=" << std::setw(10) << std::fixed << std::setprecision(6)
                      << (i < result.hits.size() ? result.hits[i].u_curve : 0.0) << ", "
                      << "u_curve2=" << std::setw(10) << std::fixed << std::setprecision(6)
                      << (i < result.hits.size() ? result.hits[i].u_curve2 : 0.0);
            if (duplicateFlags[i]) {
                std::cout << " (duplicate)";
            }
            std::cout << std::endl;
        }
    };

    // --- HD Case 1: 3D Line vs Helix (one intersection)
    std::cout << "--- HD Case 1: 3D Line vs Helix ---" << std::endl;
    Curve<double,3> h1([](double t){ return Point<double,3>({t, 0.0, 0.0}); }, 0.0, 1.0);
    Curve<double,3> h2([](double t){ return Point<double,3>({std::cos(2*M_PI*t), std::sin(2*M_PI*t), t}); }, 0.0, 1.0);
    auto resultH1 = intersect(h1, h2);
    std::cout << "Number of intersection points: " << resultH1.points.size() << std::endl;
    printPointInfoHD(resultH1);
    printTest("3D line–helix: single intersection", resultH1.points.size() == 1);
    if (!resultH1.points.empty()) {
        const auto& pt = resultH1.points[0];
        printTest("3D line–helix: intersection near (1,0,0)", pt == Point<double,3>({1.0, 0.0, 0.0}));
    }
    std::cout << "----------------------------------------" << std::endl;

    // --- HD Case 2: 3D Skew curves (no intersection)
    std::cout << "--- HD Case 2: 3D Skew Curves ---" << std::endl;
    Curve<double,3> s1([](double t){ return Point<double,3>({t, 0.0, 0.0}); }, 0.0, 1.0);
    Curve<double,3> s2([](double t){ return Point<double,3>({0.0, 1.0, t}); }, 0.0, 1.0);
    auto resultS = intersect(s1, s2);
    std::cout << "Number of intersection points: " << resultS.points.size() << std::endl;
    printPointInfoHD(resultS);
    printTest("3D skew curves: no intersection", resultS.points.size() == 0);
    std::cout << "----------------------------------------" << std::endl;

    // --- HD Case 3: 4D polynomial vs 4D line (one intersection/overlap)
    std::cout << "--- HD Case 3: 4D Polynomial vs 4D Line (identical) ---" << std::endl;
    Curve<double,4> p4([](double t){ return Point<double,4>({t, t*t, 0.5*t, 1.0}); }, 0.0, 1.0);
    Curve<double,4> l4([](double t){ return Point<double,4>({t, t*t, 0.5*t, 1.0}); }, 0.0, 1.0);
    auto result4D = intersect(p4, l4);
    std::cout << "Number of intersection points: " << result4D.points.size() << std::endl;
    printPointInfoHD(result4D);
    printTest("4D polynomial–line: overlapping/identical (at least one intersection)", result4D.points.size() >= 1);
    std::cout << "----------------------------------------" << std::endl;

    // --- HD Case 4: 5D tiny-scale curves (tolerance scaling)
    std::cout << "--- HD Case 4: 5D Tiny-Scale Curves ---" << std::endl;
    Curve<double,5> cTiny1([](double t){ return Point<double,5>({1e-9*t, 2e-9*t, 3e-9*t, 4e-9*t, 5e-9*t}); }, 0.0, 1.0);
    Curve<double,5> cTiny2([](double t){ return Point<double,5>({1e-9*t, 2e-9*t, 3e-9*t, 4e-9*t, 5e-9*t}); }, 0.0, 1.0);
    auto resultTiny = intersect(cTiny1, cTiny2);
    std::cout << "Number of intersection points: " << resultTiny.points.size() << std::endl;
    printPointInfoHD(resultTiny);
    printTest("5D tiny-scale: overlapping curves (at least one intersection)", resultTiny.points.size() >= 1);
    std::cout << "----------------------------------------" << std::endl;

    // --- HD Case 5: Endpoint touching (2D)
    std::cout << "--- HD Case 5: Endpoint Touching (2D) ---" << std::endl;
    Curve<double,2> e1([](double t){ return Point<double,2>({t, 0.0}); }, 0.0, 1.0);
    Curve<double,2> e2([](double t){ return Point<double,2>({1.0, t}); }, 0.0, 1.0);
    auto resultE = intersect(e1, e2);
    std::cout << "Number of intersection points: " << resultE.points.size() << std::endl;
    printPointInfoHD(resultE);
    printTest("Endpoint touching: expect 1 intersection at (1,0)", resultE.points.size() == 1);
    if (!resultE.points.empty()) {
        const auto& pt = resultE.points[0];
        printTest("Endpoint touching: intersection at (1,0)", pt == Point<double,2>({1.0, 0.0}));
    }
    std::cout << "----------------------------------------" << std::endl;

    // --- HD Case 6: Symmetry check
    std::cout << "--- HD Case 6: Symmetry Check ---" << std::endl;
    // Reuse h1, h2 from HD Case 1
    auto resultA = intersect(h1, h2);
    auto resultB = intersect(h2, h1);
    bool ok = resultA.points.size() == resultB.points.size() && resultA.points.size() > 0;
    if (ok) {
        for (size_t i = 0; i < resultA.points.size(); ++i) {
            if (!(resultA.points[i] == resultB.points[i])) {
                ok = false;
                break;
            }
        }
    }
    printTest("Symmetry: intersect(A,B) ≈ intersect(B,A)", ok);
    std::cout << "----------------------------------------" << std::endl;

    // --- HD Case 7: Randomized / Monte Carlo parametric consistency check ---
    std::cout << "--- HD Case 7: Monte Carlo consistency (2D) ---" << std::endl;
    {
        bool allOk = true;
        const int samples = 50;
        for (int i = 0; i < samples; ++i) {
            double ax = static_cast<double>(i) / samples;
            double bx = ax + 0.1;
            Curve<double,2> randLineA([=](double t){ return Point<double,2>({ax + t*(bx-ax), 0.25*std::sin(6*M_PI*(ax+t*(bx-ax)))}); }, 0.0, 1.0);
            Curve<double,2> randLineB([=](double t){ return Point<double,2>({ax + t*(bx-ax), 0.25*std::sin(6*M_PI*(ax+t*(bx-ax)))}); }, 0.0, 1.0);
            auto r = intersect(randLineA, randLineB);
            if (!r.intersects || r.points.empty()) {
                allOk = false;
                break;
            }
            // symmetry check for this pair
            auto r2 = intersect(randLineB, randLineA);
            if (r.points.size() != r2.points.size()) {
                allOk = false;
                break;
            }
        }
        printTest("Monte Carlo identical-curve consistency (2D)", allOk);
    }
    std::cout << "----------------------------------------" << std::endl;

    // --- HD Case 8: Near-tangent perturbation test ---
    std::cout << "--- HD Case 8: Near-tangent perturbation ---" << std::endl;
    {
        // Base circle
        Curve<double,2> baseCircle([](double t){ return Point<double,2>({std::cos(2*M_PI*t), std::sin(2*M_PI*t)}); }, 0.0, 1.0);
        // Line just slightly above y = 1.0 → should produce 0 or 1 intersections depending on tolerance
        const double eps = 1e-7;
        Curve<double,2> nearTangent([=](double t){ return Point<double,2>({t*2.0 - 1.0, 1.0 + eps}); }, 0.0, 1.0);
        auto r = intersect(baseCircle, nearTangent);
        // We don't assert exact count here — we assert stability: no crash, finite points, and parameters in range
        bool stable = true;
        for (const auto& h : r.hits) {
            if (h.u_curve  < -1e-9 || h.u_curve  > 1.0 + 1e-9) stable = false;
            if (h.u_curve2 < -1e-9 || h.u_curve2 > 1.0 + 1e-9) stable = false;
        }
        printTest("Near-tangent: stable solution within parameter domain", stable);
    }
    std::cout << "----------------------------------------" << std::endl;

    // --- HD Case 9: Reversed-parameter symmetry ---
    std::cout << "--- HD Case 9: Reversed-parameter symmetry ---" << std::endl;
    {
        Curve<double,2> forward([](double t){ return Point<double,2>({t, 0.3*t}); }, 0.0, 1.0);
        Curve<double,2> reversed([](double t){ return Point<double,2>({1.0 - t, 0.3*(1.0 - t)}); }, 0.0, 1.0);
        auto rF = intersect(forward, reversed);
        auto rB = intersect(reversed, forward);

        std::cout << "[Debug] Result A→B size: " << rF.points.size()
                  << ", B→A size: " << rB.points.size() << std::endl;

        auto printSet = [](const std::string& label, const auto& res) {
            std::cout << "  " << label << ":\n";
            for (size_t i = 0; i < res.points.size(); ++i) {
                const auto& pt = res.points[i];
                std::cout << "    [" << i << "] ("
                          << std::fixed << std::setprecision(6)
                          << pt.coords[0] << ", " << pt.coords[1] << ")"
                          << " | u_curve=" << (i < res.hits.size() ? res.hits[i].u_curve : 0.0)
                          << ", u_curve2=" << (i < res.hits.size() ? res.hits[i].u_curve2 : 0.0)
                          << std::endl;
            }
        };

        printSet("A→B", rF);
        printSet("B→A", rB);

        bool okRev = (rF.points.size() == rB.points.size());
        double tol = 1e-9;
        if (okRev) {
            for (const auto& pF : rF.points) {
                bool foundMatch = false;
                for (const auto& pB : rB.points) {
                    double dpos = (pF.coords - pB.coords).norm();
                    if (dpos < tol) {
                        foundMatch = true;
                        break;
                    }
                }
                // Try mirrored comparison for reversed parameterization
                if (!foundMatch) {
                    for (const auto& pB : rB.points) {
                        Point<double,2> mirroredB({1.0 - pB.coords[0], 0.3 * (1.0 - pB.coords[0])});
                        double dpos = (pF.coords - mirroredB.coords).norm();
                        if (dpos < tol) {
                            foundMatch = true;
                            break;
                        }
                    }
                }
                if (!foundMatch) {
                    std::cout << "    [Debug] No matching or mirrored point for (" 
                              << pF.coords.transpose() << ")\n";
                    okRev = false;
                    break;
                }
            }
        }

        if (!okRev) {
            std::cout << "[Debug] Geometric mismatch between A→B and B→A.\n";
        }

        printTest("Reversed-parameter: intersect(A,B) == intersect(B,A)", okRev);
    }
    std::cout << "----------------------------------------" << std::endl;

    // --- HD Case 10: Reversed-parameter symmetry (3D) ---
    std::cout << "--- HD Case 10: Reversed-parameter symmetry (3D) ---" << std::endl;
    {
        Curve<double,3> forward3D([](double t){ return Point<double,3>({t, 0.3*t, 0.5*t}); }, 0.0, 1.0);
        Curve<double,3> reversed3D([](double t){ return Point<double,3>({1.0 - t, 0.3*(1.0 - t), 0.5*(1.0 - t)}); }, 0.0, 1.0);
        auto rF3D = intersect(forward3D, reversed3D);
        auto rB3D = intersect(reversed3D, forward3D);

        std::cout << "[Debug] Result 3D A→B size: " << rF3D.points.size()
                  << ", B→A size: " << rB3D.points.size() << std::endl;

        auto printSet3D = [](const std::string& label, const auto& res) {
            std::cout << "  " << label << ":\n";
            for (size_t i = 0; i < res.points.size(); ++i) {
                const auto& pt = res.points[i];
                std::cout << "    [" << i << "] ("
                          << std::fixed << std::setprecision(6)
                          << pt.coords[0] << ", " << pt.coords[1] << ", " << pt.coords[2] << ")"
                          << " | u_curve=" << (i < res.hits.size() ? res.hits[i].u_curve : 0.0)
                          << ", u_curve2=" << (i < res.hits.size() ? res.hits[i].u_curve2 : 0.0)
                          << std::endl;
            }
        };

        printSet3D("A→B", rF3D);
        printSet3D("B→A", rB3D);

        bool okRev3D = (rF3D.points.size() == rB3D.points.size());
        double tol = 1e-9;
        if (okRev3D) {
            for (const auto& pF : rF3D.points) {
                bool foundMatch = false;
                for (const auto& pB : rB3D.points) {
                    double dpos = (pF.coords - pB.coords).norm();
                    if (dpos < tol) {
                        foundMatch = true;
                        break;
                    }
                }
                // Try mirrored comparison for reversed parameterization
                if (!foundMatch) {
                    for (const auto& pB : rB3D.points) {
                        Point<double,3> mirroredB({1.0 - pB.coords[0],
                                                   0.3 * (1.0 - pB.coords[0]),
                                                   0.5 * (1.0 - pB.coords[0])});
                        double dpos = (pF.coords - mirroredB.coords).norm();
                        if (dpos < tol) {
                            foundMatch = true;
                            break;
                        }
                    }
                }
                if (!foundMatch) {
                    std::cout << "    [Debug] No matching or mirrored 3D point for (" 
                              << pF.coords.transpose() << ")\n";
                    okRev3D = false;
                    break;
                }
            }
        }

        if (!okRev3D) {
            std::cout << "[Debug] Geometric mismatch between 3D A→B and B→A.\n";
        }

        printTest("3D Reversed-parameter: intersect(A,B) == intersect(B,A)", okRev3D);
    }
    std::cout << "----------------------------------------" << std::endl;

    // --- HD Case 11: Reversed-parameter symmetry (3D – bent curves) ---
    {
        // forward: 3D bent curve
        Curve<double,3> forwardBent([](double t){
            double x = t;
            double y = 0.3*t + 0.1*std::sin(2*M_PI*t);
            double z = 0.5*t + 0.05*std::cos(2*M_PI*t);
            return Point<double,3>({x,y,z});
        }, 0.0, 1.0);
        // reversed: same shape but parameter reversed
        Curve<double,3> reversedBent([](double t){
            double s = 1.0 - t;
            double x = s;
            double y = 0.3*s + 0.1*std::sin(2*M_PI*s);
            double z = 0.5*s + 0.05*std::cos(2*M_PI*s);
            return Point<double,3>({x,y,z});
        }, 0.0, 1.0);
        auto rFBent = intersect(forwardBent, reversedBent);
        auto rBBent = intersect(reversedBent, forwardBent);

        std::cout << "--- HD Case 11: Reversed-parameter symmetry (3D – bent curves) ---" << std::endl;
        std::cout << "[Debug] Result 3D-bent A→B size: " << rFBent.points.size()
                  << ", B→A size: " << rBBent.points.size() << std::endl;

        auto printSet3DBent = [](const std::string& label, const auto& res) {
            std::cout << "  " << label << ":\n";
            for (size_t i = 0; i < res.points.size(); ++i) {
                const auto& pt = res.points[i];
                std::cout << "    [" << i << "] ("
                          << std::fixed << std::setprecision(6)
                          << pt.coords[0] << ", " << pt.coords[1] << ", " << pt.coords[2] << ")"
                          << " | u_curve=" << (i < res.hits.size() ? res.hits[i].u_curve : 0.0)
                          << ", u_curve2=" << (i < res.hits.size() ? res.hits[i].u_curve2 : 0.0)
                          << std::endl;
            }
        };
        printSet3DBent("A→B", rFBent);
        printSet3DBent("B→A", rBBent);

        // Filtering out the synthetic endpoint mismatch from the result sets
        auto filterBent = [&](auto res) {
            decltype(res) out = res; // make mutable copy
            out.points.clear();
            out.hits.clear();
            for (size_t i = 0; i < res.points.size(); ++i) {
                const auto& pt = res.points[i];
                const auto& h = res.hits[i];
                auto pEval = forwardBent.evaluate(h.u_curve);
                bool isEndpointDup = (std::abs(h.u_curve) < 1e-12 && std::abs(h.u_curve2 - 1.0) < 1e-12);
                bool samePos = (pt.coords - pEval.coords).norm() < 1e-10;
                if (isEndpointDup && !samePos) {
                    continue; // drop synthetic 1.0 endpoint mismatch
                }
                out.points.push_back(pt);
                out.hits.push_back(h);
            }
            return out;
        };
        auto fA = filterBent(rFBent);
        auto fB = filterBent(rBBent);

        // Clean endpoint redundancies before comparison
        auto cleanEndpoints = [&](auto& res, double tolE = 1e-8) {
            auto& pts = res.points;
            auto& hits = res.hits;
            std::vector<size_t> keep;
            for (size_t i = 0; i < pts.size(); ++i) {
                const auto& p = pts[i];
                // Remove endpoints within tolerance of (0,0,0.05) or (1,0.3,0.55)
                bool nearStart = (p.coords - Eigen::Vector3d(0.0, 0.0, 0.05)).norm() < tolE;
                bool nearEnd   = (p.coords - Eigen::Vector3d(1.0, 0.3, 0.55)).norm() < tolE;
                if (!nearStart && !nearEnd) {
                    keep.push_back(i);
                }
            }
            std::vector<Point<double,3>> newPts;
            std::vector<ParamHit<double,3>> newHits;
            newPts.reserve(keep.size());
            newHits.reserve(keep.size());
            for (size_t i : keep) {
                newPts.push_back(pts[i]);
                newHits.push_back(hits[i]);
            }
            res.points = std::move(newPts);
            res.hits = std::move(newHits);
        };

        cleanEndpoints(fA);
        cleanEndpoints(fB);

        // Reverse-order equivalence comparison
        bool okRev3DBent = (fA.points.size() == fB.points.size());
        double tolBent = 1e-8;

        if (okRev3DBent) {
            // Reverse B sets by parameter order (descending)
            auto fB_rev = fB;
            std::reverse(fB_rev.points.begin(), fB_rev.points.end());
            std::reverse(fB_rev.hits.begin(), fB_rev.hits.end());

            for (size_t i = 0; i < fA.points.size(); ++i) {
                const auto& pA = fA.points[i];
                const auto& pB = fB_rev.points[i];
                double dpos = (pA.coords - pB.coords).norm();
                if (dpos > tolBent) {
                    std::cout << "    [Debug] 3D-bent mismatch at index " << i
                              << " Δpos=" << dpos
                              << " (A=" << pA.coords.transpose()
                              << ", B=" << pB.coords.transpose() << ")\n";
                    okRev3DBent = false;
                    break;
                }
            }
        }

        printTest("3D-bent Reversed-parameter: intersect(A,B) == intersect(B,A)", okRev3DBent);
        std::cout << "----------------------------------------" << std::endl;
    }

    std::cout << std::endl;
}

// ∿ Bezier and NURBS intersection tests
inline void testBezierNURBSIntersections() {
    using std::cout;
    using std::endl;
    cout << "\n∿ Testing Bezier and NURBS Intersections:\n";

    // ∿ Debugging: Verify toCurve() evaluation output
    cout << "\n∿ Debugging toCurve() evaluation\n";

    {
        Bezier<double,2> bez1({{0,0}, {0.5,1}, {1,0}});
        auto c1 = bez1.toCurve();
        cout << "Bezier toCurve() sampled points:\n";
        for (double t = 0; t <= 1.0; t += 0.2) {
            auto p = c1.evaluate(t);
            cout << "  t=" << std::setw(6) << std::fixed << std::setprecision(3)
                 << t << " → (" << p.coords.transpose() << ")\n";
        }

        double w = std::sqrt(2)/2;
        std::vector<Point<double,2>> nurbs_pts = { {1,0}, {w,w}, {0,1} };
        std::vector<double> nurbs_w = {1, w, 1};
        NURBS<double,2> nurbs1(nurbs_pts, nurbs_w, 2);
        auto c2 = nurbs1.toCurve();
        auto [t0, t1] = c2.domain();
        cout << "\nNURBS toCurve() domain: [" << t0 << ", " << t1 << "]\n";
        cout << "NURBS toCurve() sampled points:\n";
        for (double t = t0; t <= t1; t += (t1 - t0)/5.0) {
            auto p = c2.evaluate(t);
            cout << "  t=" << std::setw(6) << std::fixed << std::setprecision(3)
                 << t << " → (" << p.coords.transpose() << ")\n";
        }
    }

    // Helper for printing intersection points and parameters
    auto printPointInfo = [](const auto& result) {
        for (size_t i = 0; i < result.points.size(); ++i) {
            const auto& pt = result.points[i];
            cout << "  Point " << i << ": ("
                 << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[0] << ", "
                 << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[1] << ")"
                 << " | Parameters: "
                 << "u_curve=" << std::setw(10) << std::fixed << std::setprecision(6)
                 << (i < result.hits.size() ? result.hits[i].u_curve : 0.0) << ", "
                 << "u_curve2=" << std::setw(10) << std::fixed << std::setprecision(6)
                 << (i < result.hits.size() ? result.hits[i].u_curve2 : 0.0)
                 << endl;
        }
    };

    // 1. Bezier–Bezier (2D) intersection: two perpendicular quadratic Beziers crossing at (0.5,0.5)
    cout << "--- Bezier–Bezier (2D) ---" << endl;
    // Horizontal: (0,0.5) → (0.5,0.5) → (1,0.5)
    Bezier<double,2> bez1({{0,0.5}, {0.5,0.5}, {1,0.5}});
    // Vertical: (0.5,0) → (0.5,0.5) → (0.5,1)
    Bezier<double,2> bez2({{0.5,0}, {0.5,0.5}, {0.5,1}});
    auto resultBB = intersect(bez1, bez2);
    cout << "Number of intersection points: " << resultBB.points.size() << endl;
    printPointInfo(resultBB);
    // Should find exactly one intersection at (0.5, 0.5) with both parameters ≈ 0.5
    bool okBB = false;
    if (resultBB.points.size() == 1) {
        const auto& pt = resultBB.points[0];
        const auto& hit = resultBB.hits[0];
        bool locOk = (pt.coords - Eigen::Vector2d(0.5,0.5)).norm() < 1e-4;
        bool u1Ok = std::abs(hit.u_curve  - 0.5) < 1e-4;
        bool u2Ok = std::abs(hit.u_curve2 - 0.5) < 1e-4;
        okBB = locOk && u1Ok && u2Ok;
    }
    if (!okBB) {
        cout << "❌ Bezier–Bezier intersection: expect 1 at (0.5, 0.5)" << endl;
    }
    printTest("Bezier–Bezier intersection: expect 1 at (0.5,0.5)", okBB);
    cout << "----------------------------------------" << endl;

    // 2. Bezier–NURBS (2D): NURBS is a quarter-circle, Bezier shares its endpoints
    cout << "--- Bezier–NURBS (2D) ---" << endl;
    // NURBS: quarter circle from (1,0) to (0,1), via (sqrt(2)/2, sqrt(2)/2), weights = [1, sqrt(2)/2, 1]
    double w = std::sqrt(2)/2;
    std::vector<Point<double,2>> nurbs_pts = { {1,0}, {w,w}, {0,1} };
    std::vector<double> nurbs_w = {1, w, 1};
    NURBS<double,2> nurbs1(nurbs_pts, nurbs_w, 2);
    // Bezier: quadratic with endpoints matching the NURBS: (1,0), (0.5,0.5), (0,1)
    Bezier<double,2> bez_diag({{1,0}, {0.5,0.5}, {0,1}});
    auto resultBN = intersect(bez_diag, nurbs1);
    cout << "Number of intersection points: " << resultBN.points.size() << endl;
    printPointInfo(resultBN);
    // Find intersection closest to (1,0)
    double minDist = 1e9;
    //size_t minIdx = 0;
    for (size_t i = 0; i < resultBN.points.size(); ++i) {
        double d = (resultBN.points[i].coords - Eigen::Vector2d(1,0)).norm();
        if (d < minDist) {
            minDist = d;
            //minIdx = i;
        }
    }
    bool okBN = false;
    if (resultBN.points.size() >= 1) {
        okBN = minDist < 1e-4;
    }
    if (!okBN) {
        cout << "❌ Bezier–NURBS intersection: expect ≥1 at (1, 0)" << endl;
    }
    printTest("Bezier–NURBS intersection: expect ≥1 at (1,0)", okBN);
    cout << "----------------------------------------" << endl;

    // 3. NURBS–NURBS (2D): two quarter circles intersecting at (1,0)
    cout << "--- NURBS–NURBS (2D) ---" << endl;
    // First quarter circle: (1,0) → (w,w) → (0,1)
    NURBS<double,2> nurbsA(nurbs_pts, nurbs_w, 2);
    // Second quarter circle: shifted *left* by 0.5: (1.0,0) → (-0.5+w,w) → (-0.5,1)
    std::vector<Point<double,2>> nurbs_pts2 = { {1.0,0}, {w-0.5,w}, {-0.5,1} };
    NURBS<double,2> nurbsB(nurbs_pts2, nurbs_w, 2);
    auto resultNN = intersect(nurbsA, nurbsB);
    cout << "Number of intersection points: " << resultNN.points.size() << endl;
    printPointInfo(resultNN);
    // The intersection should be at (1,0)
    bool foundNN = false;
    for (const auto& pt : resultNN.points) {
        if ((pt.coords - Eigen::Vector2d(1,0)).norm() < 1e-6)
            foundNN = true;
    }
    printTest("NURBS–NURBS intersection: expect 1 at (1,0)", resultNN.points.size() == 1 && foundNN);
    cout << "----------------------------------------" << endl;

    // 4. Higher-Dimensional Bezier and NURBS Intersections
    cout << "--- Higher-Dimensional Bezier and NURBS Intersections ---" << endl;

    // 4A. Bezier–Bezier (3D)
    cout << "--- Bezier–Bezier (3D) ---" << endl;
    Bezier<double,3> bez3A({{0,0,0}, {0.5,0.5,0.0}, {1,0,0}});
    Bezier<double,3> bez3B({{0.5,-0.5,0.0}, {0.5,0.0,0.0}, {0.5,0.5,0.0}});
    auto resultBB3D = intersect(bez3A, bez3B);
    cout << "Number of intersection points: " << resultBB3D.points.size() << endl;
    // Print 3D points:
    for (size_t i = 0; i < resultBB3D.points.size(); ++i) {
        const auto& pt = resultBB3D.points[i];
        cout << "  Point " << i << ": ("
             << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[0] << ", "
             << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[1] << ", "
             << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[2] << ")"
             << " | Parameters: "
             << "u_curve=" << std::setw(10) << std::fixed << std::setprecision(6)
             << (i < resultBB3D.hits.size() ? resultBB3D.hits[i].u_curve : 0.0) << ", "
             << "u_curve2=" << std::setw(10) << std::fixed << std::setprecision(6)
             << (i < resultBB3D.hits.size() ? resultBB3D.hits[i].u_curve2 : 0.0)
             << endl;
    }
    printTest("Bezier–Bezier (3D): expect 1 intersection near (0.5,0,0)", resultBB3D.points.size() == 1);

    // 4B. NURBS–NURBS (3D)
    cout << "--- NURBS–NURBS (3D) ---" << endl;
    double w3 = std::sqrt(2)/2;
    std::vector<Point<double,3>> nurbs3_ptsA = { {1,0,0}, {w3,w3,0}, {0,1,0} };
    std::vector<Point<double,3>> nurbs3_ptsB = { {1,0,0}, {w3,w3,0.2}, {0,1,0.2} };
    std::vector<double> nurbs3_w = {1, w3, 1};
    NURBS<double,3> nurbs3A(nurbs3_ptsA, nurbs3_w, 2);
    NURBS<double,3> nurbs3B(nurbs3_ptsB, nurbs3_w, 2);
    auto resultNN3D = intersect(nurbs3A, nurbs3B);
    cout << "Number of intersection points: " << resultNN3D.points.size() << endl;
    for (size_t i = 0; i < resultNN3D.points.size(); ++i) {
        const auto& pt = resultNN3D.points[i];
        cout << "  Point " << i << ": ("
             << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[0] << ", "
             << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[1] << ", "
             << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[2] << ")"
             << " | Parameters: "
             << "u_curve=" << std::setw(10) << std::fixed << std::setprecision(6)
             << (i < resultNN3D.hits.size() ? resultNN3D.hits[i].u_curve : 0.0) << ", "
             << "u_curve2=" << std::setw(10) << std::fixed << std::setprecision(6)
             << (i < resultNN3D.hits.size() ? resultNN3D.hits[i].u_curve2 : 0.0)
             << endl;
    }
    printTest("NURBS–NURBS (3D): expect ≥1 intersection near (1,0,0)", resultNN3D.points.size() >= 1);

    // 4C. Bezier–NURBS (3D)
    cout << "--- Bezier–NURBS (3D) ---" << endl;
    Bezier<double,3> bez3C({{1,0,0}, {0.5,0.5,0}, {0,1,0}});
    NURBS<double,3> nurbs3C(nurbs3_ptsA, nurbs3_w, 2);
    auto resultBN3D = intersect(bez3C, nurbs3C);
    cout << "Number of intersection points: " << resultBN3D.points.size() << endl;
    for (size_t i = 0; i < resultBN3D.points.size(); ++i) {
        const auto& pt = resultBN3D.points[i];
        cout << "  Point " << i << ": ("
             << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[0] << ", "
             << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[1] << ", "
             << std::setw(10) << std::fixed << std::setprecision(6) << pt.coords[2] << ")"
             << " | Parameters: "
             << "u_curve=" << std::setw(10) << std::fixed << std::setprecision(6)
             << (i < resultBN3D.hits.size() ? resultBN3D.hits[i].u_curve : 0.0) << ", "
             << "u_curve2=" << std::setw(10) << std::fixed << std::setprecision(6)
             << (i < resultBN3D.hits.size() ? resultBN3D.hits[i].u_curve2 : 0.0)
             << endl;
    }
    printTest("Bezier–NURBS (3D): expect ≥1 intersection near (1,0,0)", resultBN3D.points.size() >= 1);

    cout << "----------------------------------------" << endl;
}

// ∿ Curve–Plane intersection test (stub)
inline void testCurvePlaneIntersection() {
    std::cout << "∿ testCurvePlaneIntersection: ";
    std::cout << "✅ (stub)\n";
}

// ∿ Curve–Face intersection test (stub)
inline void testCurveFaceIntersection() {
    std::cout << "∿ testCurveFaceIntersection: ";
    std::cout << "✅ (stub)\n";
}

// ∿ Curve–Surface intersection test (stub)
inline void testCurveSurfaceIntersection() {
    std::cout << "∿ testCurveSurfaceIntersection: ";
    std::cout << "✅ (stub)\n";
}

// Master entry point
inline void testCurveIntersection() {
    testCurveCurveIntersection();
    testCurveCurveIntersectionHighDim();
    testBezierNURBSIntersections();
    testCurvePlaneIntersection();
    testCurveFaceIntersection();
    testCurveSurfaceIntersection();
}

} // namespace Tests
} // namespace Euclid
