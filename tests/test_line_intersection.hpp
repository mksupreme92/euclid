#pragma once


#include "geometry/geometry.hpp"
#include "geometry/intersection/line_intersection.hpp"
#include "test_utilities.hpp"

namespace Euclid {
namespace Tests {

inline void testLineLineIntersection() {
    // 1. Line-Line intersection at a point
    {
        Line l1(Point<double,3>({0, 0, 0}), Eigen::Vector3d({1, 0, 0}));
        Line l2(Point<double,3>({0, 0, 0}), Eigen::Vector3d({0, 1, 0}));
        auto result = intersect(l1, l2);
        printTest("Line-Line intersecting at a point: has intersection", result.intersects);
        if (result.intersects) {
            printTest("Line-Line intersection result is a Point", !result.points.empty());
            if (!result.points.empty()) {
                Point<double,3> p = result.points[0];
                printTest("Line-Line intersection at origin", p == Point<double,3>({0,0,0}));
            }
        }
    }

    // 2. Line-Line skew (no intersection)
    {
        Line l1(Point<double,3>({0, 0, 0}), Eigen::Vector3d({1, 0, 0}));
        Line l2(Point<double,3>({0, 1, 1}), Eigen::Vector3d({0, 1, 0}));
        auto result = intersect(l1, l2);
        printTest("Line-Line skew: no intersection", !result.intersects);
    }

    // 3. Line-Line coincident
    {
        Line l1(Point<double,3>({1, 2, 3}), Eigen::Vector3d({2, 2, 2}));
        Line l2(Point<double,3>({3, 4, 5}), Eigen::Vector3d({1, 1, 1}));
        auto result = intersect(l1, l2);
        printTest("Line-Line coincident: has intersection", result.intersects);
        if (result.intersects) {
            printTest("Line-Line coincident: result is a Line", !result.lines.empty());
        }
    }

    // Additional tests from test_line.hpp migrated here:

    // 6. 2D Line-Line intersection
    {
        Line l1(Point<double,2>({0,0}), Eigen::Vector2d({1,1}));
        Line l2(Point<double,2>({0,1}), Eigen::Vector2d({1,-1}));
        auto result = intersect(l1, l2);
        printTest("2D Line-Line intersection: has intersection", result.intersects);
        if (result.intersects) {
            printTest("2D Line-Line intersection result is a Point", !result.points.empty());
            if (!result.points.empty()) {
                Point<double,2> p = result.points[0];
                printTest("2D Line-Line intersection at (0.5,0.5)", p == Point<double,2>({0.5,0.5}));
            }
        }
    }

    // 7. 3D Line-Line intersection (skew lines)
    {
        Line l1(Point<double,3>({0,0,0}), Eigen::Vector3d({1,0,0}));
        Line l2(Point<double,3>({0,1,1}), Eigen::Vector3d({0,1,0}));
        auto result = intersect(l1, l2);
        printTest("3D Line-Line skew: no intersection", !result.intersects);
    }

    // 8. 4D Line-Line intersection (coincident lines)
    {
        Line l1(Point<double,4>({1,2,3,4}), Eigen::Vector4d({1,1,1,1}));
        Line l2(Point<double,4>({2,3,4,5}), Eigen::Vector4d({2,2,2,2}));
        auto result = intersect(l1, l2);
        printTest("4D Line-Line coincident: has intersection", result.intersects);
        if (result.intersects) {
            printTest("4D Line-Line coincident: result is a Line", !result.lines.empty());
        }
    }

    // 9. 5D Line-Line intersection (parallel lines, no intersection)
    {
        Eigen::Matrix<double,5,1> dir1; dir1 << 1,0,0,0,0;
        Eigen::Matrix<double,5,1> dir2; dir2 << 1,0,0,0,0;
        Line l1(Point<double,5>({0,0,0,0,0}), dir1);
        Line l2(Point<double,5>({0,1,0,0,0}), dir2);
        auto result = intersect(l1, l2);
        printTest("5D Line-Line parallel: no intersection", !result.intersects);
    }
}

inline void testLinePlaneIntersection() {
    // 4. Line-Plane intersection
    {
        Line l(Point<double,3>({0,0,0}), Eigen::Vector3d({0,0,1}));
        Plane pl(Point<double,3>({0,0,5}), Eigen::Vector3d({0,0,1}));
        auto result = intersect(l, pl);
        printTest("Line-Plane: has intersection", result.intersects);
        if (result.intersects) {
            printTest("Line-Plane: result is a Point", !result.points.empty());
            if (!result.points.empty()) {
                Point<double,3> p = result.points[0];
                printTest("Line-Plane intersection at (0,0,5)", p == Point<double,3>({0,0,5}));
            }
        }
    }
}

// Trivial tests for ND (4D, 5D) line-plane intersection
inline void testLinePlaneIntersectionND() {
    using namespace Euclid;
    using std::cout;
    cout << "ND Line-Plane Intersection Tests\n";

    // 4D: line intersects plane
    {
        // Plane: passes through (1,2,3,4), normal (0,0,0,1)
        Point<double,4> plane_point({1,2,3,4});
        Eigen::Matrix<double,4,1> plane_normal; plane_normal << 0,0,0,1;
        Plane pl(plane_point, plane_normal);
        // Line: passes through (1,2,3,0), direction (0,0,0,1)
        Point<double,4> line_point({1,2,3,0});
        Eigen::Matrix<double,4,1> line_dir; line_dir << 0,0,0,1;
        Line l(line_point, line_dir);
        auto result = intersect(l, pl);
        printTest("4D Line-Plane: has intersection", result.intersects);
        if (result.intersects && !result.points.empty()) {
            const auto& p = result.points[0];
            std::cout << "  4D intersection point: (";
            for (int i = 0; i < 4; ++i) {
                std::cout << p[i];
                if (i < 3) std::cout << ", ";
            }
            std::cout << ")\n";
        }
    }

    // 5D: intersection case
    {
        Point<double,5> plane_point({0,0,0,0,0});
        Eigen::Matrix<double,5,1> plane_normal; plane_normal << 1,0,0,0,0;
        Plane pl(plane_point, plane_normal);

        Point<double,5> line_point({0,1,2,3,4});
        Eigen::Matrix<double,5,1> line_dir; line_dir << 1,0,0,0,0; // direction along plane normal

        Line l(line_point, line_dir);
        auto result = intersect(l, pl);

        std::cout << "5D Line-Plane Intersection Case Debug:\n";
        std::cout << "Line point: "; for (int i=0;i<5;i++){std::cout<<line_point[i]<<(i<4?",":"");} std::cout<<"\n";
        std::cout << "Line direction: "; for (int i=0;i<5;i++){std::cout<<line_dir[i]<<(i<4?",":"");} std::cout<<"\n";
        std::cout << "Plane point: "; for (int i=0;i<5;i++){std::cout<<plane_point[i]<<(i<4?",":"");} std::cout<<"\n";
        std::cout << "Plane normal: "; for (int i=0;i<5;i++){std::cout<<plane_normal[i]<<(i<4?",":"");} std::cout<<"\n";
        double dot_product = plane_normal.dot(line_dir);
        std::cout << "Dot product (normal • direction) = " << dot_product << "\n";

        printTest("5D Line-Plane intersection: has intersection", result.intersects);
        if(result.intersects && !result.points.empty()) {
            const auto& p = result.points[0];
            std::cout << "  Intersection point: ("; for(int i=0;i<5;i++){std::cout<<p[i]<<(i<4?", ":"");} std::cout<<")\n";
        }
    }

    // 5D: no intersection case
    {
        Point<double,5> plane_point({0,0,0,0,0});
        Eigen::Matrix<double,5,1> plane_normal; plane_normal << 1,0,0,0,0;
        Plane pl(plane_point, plane_normal);

        Point<double,5> line_point({1,1,2,3,4}); // shifted along normal
        Eigen::Matrix<double,5,1> line_dir; line_dir << 0,1,0,0,0; // direction parallel to plane

        Line l(line_point, line_dir);
        auto result = intersect(l, pl);

        std::cout << "5D Line-Plane No-Intersection Case Debug:\n";
        std::cout << "Line point: "; for (int i=0;i<5;i++){std::cout<<line_point[i]<<(i<4?",":"");} std::cout<<"\n";
        std::cout << "Line direction: "; for (int i=0;i<5;i++){std::cout<<line_dir[i]<<(i<4?",":"");} std::cout<<"\n";
        std::cout << "Plane point: "; for (int i=0;i<5;i++){std::cout<<plane_point[i]<<(i<4?",":"");} std::cout<<"\n";
        std::cout << "Plane normal: "; for (int i=0;i<5;i++){std::cout<<plane_normal[i]<<(i<4?",":"");} std::cout<<"\n";
        double dot_product = plane_normal.dot(line_dir);
        std::cout << "Dot product (normal • direction) = " << dot_product << "\n";

        printTest("5D Line-Plane: no intersection", !result.intersects);
        if(result.intersects && !result.points.empty()) {
            const auto& p = result.points[0];
            std::cout << "  Intersection point: ("; for(int i=0;i<5;i++){std::cout<<p[i]<<(i<4?", ":"");} std::cout<<")\n";
        }
    }
}


inline void testLineCurveIntersection() {
    // 1. Line intersecting a 3D curve (parabola)
    {
        Curve<double, 3> parabola(
            [](double t) { return Point<double,3>({t, t*t, 0}); },
            -1.0, 2.0
        );
        Line<double, 3> l(Point<double,3>({0,0,0}), Eigen::Vector3d({1,1,0}));
        auto result = intersect(l, parabola);
        printTest("Line-Curve (parabola) intersection: has intersection", result.intersects);
        if (result.intersects && !result.points.empty()) {
            const auto& p = result.points[0];
            printTest("Line-Curve intersection point", true);
            std::cout << "Intersection at: (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
        }
    }

    // 2. Line not intersecting the curve
    {
        Curve<double, 3> parabola(
            [](double t) { return Point<double,3>({t, t*t, 0}); },
            -1.0, 2.0
        );
        Line<double, 3> l(Point<double,3>({0,0,1}), Eigen::Vector3d({1,1,0})); // Line is offset in z=1 plane
        auto result = intersect(l, parabola);
        printTest("Line-Curve (parabola) no intersection", !result.intersects);
    }

    // 2b. Analytical resolution limit test
    {
        std::cout << "---- Analytical Frequency Resolution Limit Test ----" << std::endl;

        Line<double, 2> axis(Point<double, 2>({0, 0}), Eigen::Vector2d({1, 0}));
        double freq = 1.0;

        while (true) {
            // Define a sine curve with `freq` oscillations over [0, 1]
            Curve<double, 2> wave(
                [freq](double t) {
                    return Point<double, 2>({t, std::sin(2.0 * M_PI * freq * t)});
                },
                0.0, 1.0
            );

            // Analytical number of zero-crossings:
            // Each full period contributes two crossings (up/down),
            // plus one at t=0 and one at t=1.
            int expected = static_cast<int>(2 * freq + 1);

            auto result = intersect(axis, wave);
            int found = static_cast<int>(result.points.size());

            bool ok = std::abs(found - expected) <= 1;
            std::cout << "freq=" << freq
                      << "  expected=" << expected
                      << "  found=" << found
                      << (ok ? " ✅" : " ❌") << std::endl;

            if (!ok) {
                std::cout << "❌ Analytical mismatch: intersection resolution limit reached at frequency = "
                          << freq << std::endl;
                break;
            }

            // Double frequency until failure
            freq *= 2.0;

            // Safety cap to avoid infinite loop in case of perfect numerical behavior
            if (freq > 1e9) {
                std::cout << "✅ No analytical mismatch detected up to freq = " << freq << std::endl;
                break;
            }
        }
    }

    // 2c. Frequency resolution test with line just below sine peaks
    {
        std::cout << "---- Frequency Test: Line below sine peaks ----" << std::endl;

        double freq = 1.0;
        while (true) {
            // Define a sine curve with `freq` oscillations over [0, 1]
            Curve<double, 2> wave(
                [freq](double t) {
                    return Point<double, 2>({t, std::sin(2.0 * M_PI * freq * t)});
                },
                0.0, 1.0
            );

            // Define a horizontal line slightly below the sine peaks (y = 0.9)
            Line<double, 2> lineBelowPeak(Point<double, 2>({0.0, 0.9}), Eigen::Vector2d({1.0, 0.0}));

            // Expected crossings: two per period where sin(x)=0.9 (both up and down)
            // Solve sin(2πf t)=0.9 ⇒ crossings occur at ±arcsin(0.9)/(2πf) in each half-period.
            // Analytical number of crossings ≈ 2 * freq (same as zero-crossings)
            int expected = static_cast<int>(2 * freq);

            auto result = intersect(lineBelowPeak, wave);
            int found = static_cast<int>(result.points.size());

            bool ok = std::abs(found - expected) <= 2;
            std::cout << "freq=" << freq
                      << "  expected≈" << expected
                      << "  found=" << found
                      << (ok ? " ✅" : " ❌") << std::endl;

            if (!ok) {
                std::cout << "❌ Analytical mismatch (below peaks): resolution limit reached at frequency = "
                          << freq << std::endl;
                break;
            }

            freq *= 2.0;
            if (freq > 1e9) {
                std::cout << "✅ No mismatch detected up to freq = " << freq << std::endl;
                break;
            }
        }
    }

    /*
    // 2d. Empirical plot: missed intersections vs tolerance (zero line and below-peak line)
    {
        std::cout << "---- Empirical Plot: Missed vs Tolerance ----" << std::endl;

        const double freq = 4096.0; // near the resolution limit
        // Sine curve over [0,1]
        Curve<double, 2> wave(
            [freq](double t) { return Point<double, 2>({t, std::sin(2.0 * M_PI * freq * t)}); },
            0.0, 1.0
        );

        // Lines: y=0 (zero-crossings) and y=0.9 (near peaks)
        Line<double, 2> axis(Point<double, 2>({0.0, 0.0}), Eigen::Vector2d({1.0, 0.0}));
        Line<double, 2> belowPeak(Point<double, 2>({0.0, 0.9}), Eigen::Vector2d({1.0, 0.0}));

        const int expectedZero = static_cast<int>(2 * freq + 1);
        const int expectedPeak = static_cast<int>(2 * freq);

        std::vector<double> tolerances = {1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12};

        // CSV for zero line
        {
            std::ofstream csv("missed_vs_tol_zero.csv");
            csv << "tolerance,found,expected,missed\n";
            for (double tval : tolerances) {
                Euclid::Tolerance tol(1e-6, 1e-9, tval, 1.0);//( absTol=1e-6, 1e-9, paramTol=tval, evalFactor=1.0)
                auto res = intersect(axis, wave, tol);
                int found = static_cast<int>(res.points.size());
                int missed = std::max(0, expectedZero - found);
                std::cout << "ZERO tol=" << tval << "  found=" << found << "  missed=" << missed << "\n";
                csv << tval << "," << found << "," << expectedZero << "," << missed << "\n";
            }
        }

        // CSV for line below peaks
        {
            std::ofstream csv("missed_vs_tol_peak.csv");
            csv << "tolerance,found,expected,missed\n";
            for (double tval : tolerances) {
                Euclid::Tolerance tol(1e-6, 1e-9, tval, 1.0); //( absTol=1e-6, 1e-9, paramTol=tval, evalFactor=1.0)
                auto res = intersect(belowPeak, wave, tol);
                int found = static_cast<int>(res.points.size());
                int missed = std::max(0, expectedPeak - found);
                std::cout << "PEAK tol=" << tval << "  found=" << found << "  missed=" << missed << "\n";
                csv << tval << "," << found << "," << expectedPeak << "," << missed << "\n";
            }
        }

        std::cout << "✅ Wrote missed_vs_tol_zero.csv and missed_vs_tol_peak.csv" << std::endl;
    }
   */

    // 3. Line intersecting sine curve (multiple crossings)
    {
        Curve<double, 2> sineCurve(
            [](double t) { return Point<double,2>({t, std::sin(M_PI * t)}); },
            -1.0, 3.0
        );
        Line<double, 2> l(Point<double,2>({-1.0, 0.0}), Eigen::Vector2d({1.0, 0.0}));
        auto result = intersect(l, sineCurve);
        printTest("Line-Curve (sine wave) intersection: has intersection", result.intersects);
        std::cout << "DEBUG: sine curve intersection count = " << result.points.size() << "\n";
        for (size_t i = 0; i < result.points.size(); ++i) {
            const auto& p = result.points[i];
            std::cout << "  sine point[" << i << "] = (" << p[0] << ", " << p[1] << ")\n";
        }
    }

    // 4. Line intersecting circle (two intersections)
    {
        Curve<double, 2> circle(
            [](double theta) { return Point<double,2>({std::cos(theta), std::sin(theta)}); },
            0.0, 2*M_PI
        );
        Line<double, 2> l(Point<double,2>({0.5, -2.0}), Eigen::Vector2d({0.0, 1.0}));
        auto result = intersect(l, circle);
        printTest("Line-Curve (circle) intersection: has intersection", result.intersects);
        std::cout << "DEBUG: circle intersection count = " << result.points.size() << "\n";
        for (size_t i = 0; i < result.points.size(); ++i) {
            const auto& p = result.points[i];
            std::cout << "  circle point[" << i << "] = (" << p[0] << ", " << p[1] << ")\n";
        }
    }
}

inline void testLineSurfaceIntersection() {
    using namespace Euclid;
    using std::cout;
    // 1. Line intersecting planar surface z=0 (as parametric surface)
    {
        auto surface = Surface<double,3>(
            [](double u, double v) { return Point<double,3>({u, v, 0.0}); },
            std::make_pair(std::make_pair(-10.0, 10.0), std::make_pair(-10.0, 10.0))
        );
        Line l(Point<double,3>({0,0,-1}), Eigen::Vector3d({0,0,1}));
        std::cout << "DEBUG: Line start: (0,0,-1), direction: (0,0,1)\n";
        std::cout << "DEBUG: Surface point: (0,0,0), normal: (0,0,1)\n";
        auto result = intersect(l, surface);
        printTest("Line-Surface (plane z=0) intersection: has intersection", result.intersects);
        if (result.intersects) {
            printTest("Line-Surface intersection result is a Point", !result.points.empty());
            if (!result.points.empty()) {
                Point<double,3> p = result.points[0];
                printTest("Line-Surface intersection at (0,0,0)", p == Point<double,3>({0,0,0}));
            }
        }
    }

    // 2. Line missing the surface (parallel, above plane)
    {
        auto surface = Surface<double,3>(
            [](double u, double v) { return Point<double,3>({u, v, 0.0}); },
            std::make_pair(std::make_pair(-10.0, 10.0), std::make_pair(-10.0, 10.0))
        );
        Line l(Point<double,3>({0,0,1}), Eigen::Vector3d({1,0,0})); // Parallel and above plane
        std::cout << "DEBUG: Line start: (0,0,1), direction: (1,0,0)\n";
        std::cout << "DEBUG: Surface point: (0,0,0), normal: (0,0,1)\n";
        auto result = intersect(l, surface);
        printTest("Line-Surface no intersection", !result.intersects);
    }

    // 2a. Additional test: Line through (0,0,0) (segment start point in segment-surface tests)
    {
        auto surface = Surface<double,3>(
            [](double u, double v) { return Point<double,3>({u, v, 0.0}); },
            std::make_pair(std::make_pair(-10.0, 10.0), std::make_pair(-10.0, 10.0))
        );
        Line l(Point<double,3>({0,0,0}), Eigen::Vector3d({1,1,1}));
        std::cout << "DEBUG: Line start: (0,0,0), direction: (1,1,1)\n";
        std::cout << "DEBUG: Surface point: (0,0,0), normal: (0,0,1)\n";
        auto result = intersect(l, surface);
        printTest("Line-Surface (line through segment start): has intersection", result.intersects);
        if (result.intersects && !result.points.empty()) {
            Point<double,3> p = result.points[0];
            std::cout << "DEBUG: Intersection point: (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
            printTest("Line-Surface intersection at (0,0,0)", p == Point<double,3>({0,0,0}));
        }
    }

    // 2b. Additional test: Line through (1,1,1) (segment end point in segment-surface tests)
    {
        auto surface = Surface<double,3>(
            [](double u, double v) { return Point<double,3>({u, v, 0.0}); },
            std::make_pair(std::make_pair(-10.0, 10.0), std::make_pair(-10.0, 10.0))
        );
        Line l(Point<double,3>({1,1,1}), Eigen::Vector3d({1,1,1}));
        std::cout << "DEBUG: Line start: (1,1,1), direction: (1,1,1)\n";
        std::cout << "DEBUG: Surface point: (0,0,0), normal: (0,0,1)\n";
        auto result = intersect(l, surface);
        printTest("Line-Surface (line through segment end): has intersection", result.intersects);
        if (result.intersects && !result.points.empty()) {
            Point<double,3> p = result.points[0];
            std::cout << "DEBUG: Intersection point: (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
            printTest("Line-Surface intersection at (0,0,0)", p == Point<double,3>({0,0,0}));
        }
    }

    // 3. 4D: Line-Surface (w=0) intersection (as parametric surface)
    {
        // 4D surface w = 0, domain (-1,1)x(-1,1)
        auto surface4 = Surface<double,4>(
            [](double u, double v) { return Point<double,4>({u, v, 0.0, 0.0}); },
            std::make_pair(std::make_pair(-1.0, 1.0), std::make_pair(-1.0, 1.0))
        );
        // Line: passes through (0.5,0.5,-1,0), direction (0,0,1,0)
        Point<double,4> line_point({0.5,0.5,-1,0});
        Eigen::Matrix<double,4,1> line_dir; line_dir << 0,0,1,0;
        Line l(line_point, line_dir);
        auto result = intersect(l, surface4);
        printTest("4D Line-Surface intersection: has intersection", result.intersects);
        if (result.intersects && !result.points.empty()) {
            const auto& p = result.points[0];
            std::cout << "  4D Line-Surface intersection point: (";
            for (int i = 0; i < 4; ++i) {
                std::cout << p[i];
                if (i < 3) std::cout << ", ";
            }
            std::cout << ")\n";
        }
    }

    // 4. 4D: Line missing surface (parallel, offset)
    {
        // 4D surface w = 0, domain (-1,1)x(-1,1)
        auto surface4 = Surface<double,4>(
            [](double u, double v) { return Point<double,4>({u, v, 0.0, 0.0}); },
            std::make_pair(std::make_pair(-1.0, 1.0), std::make_pair(-1.0, 1.0))
        );
        // Line: passes through (0.5,0.5,-1,1), direction (0,0,1,0) (parallel, but offset in w)
        Point<double,4> line_point({0.5,0.5,-1,1});
        Eigen::Matrix<double,4,1> line_dir; line_dir << 0,0,1,0;
        Line l(line_point, line_dir);
        auto result = intersect(l, surface4);
        printTest("4D Line-Surface no intersection", !result.intersects);
    }

    // 5. 5D: Line-Surface intersection (z=0 plane as surface)
    {
        // 5D surface z = 0.0, domain (-1,1)x(-1,1)
        auto surface5 = Surface<double,5>(
            [](double u, double v) { return Point<double,5>({u, v, 0.0, 0.0, 0.0}); },
            std::make_pair(std::make_pair(-1.0, 1.0), std::make_pair(-1.0, 1.0))
        );
        // Line: passes through (0.5,0.5,-1,0,0), direction (0,0,1,0,0)
        Point<double,5> line_point({0.5,0.5,-1,0,0});
        Eigen::Matrix<double,5,1> line_dir; line_dir << 0,0,1,0,0;
        Line l(line_point, line_dir);
        auto result = intersect(l, surface5);
        printTest("5D Line-Surface intersection: has intersection", result.intersects);
        if (result.intersects && !result.points.empty()) {
            const auto& p = result.points[0];
            std::cout << "  5D Line-Surface intersection point: (";
            for (int i = 0; i < 5; ++i) {
                std::cout << p[i];
                if (i < 4) std::cout << ", ";
            }
            std::cout << ")\n";
        }
    }

    // 6. 5D: Line missing surface (parallel, offset)
    {
        // 5D surface z = 0.0, domain (-1,1)x(-1,1)
        auto surface5 = Surface<double,5>(
            [](double u, double v) { return Point<double,5>({u, v, 0.0, 0.0, 0.0}); },
            std::make_pair(std::make_pair(-1.0, 1.0), std::make_pair(-1.0, 1.0))
        );
        // Line: passes through (0.5,0.5,2,0,0), direction (0,1,0,0,0) (parallel to z=0 plane, offset in z)
        Point<double,5> line_point({0.5,0.5,2,0,0});
        Eigen::Matrix<double,5,1> line_dir; line_dir << 0,1,0,0,0;
        Line l(line_point, line_dir);
        auto result = intersect(l, surface5);
        printTest("5D Line-Surface no intersection", !result.intersects);
    }

    // --- Additional tests: Lines from start/end of previously failing segment-surface intersections ---
    // 3D: start (0.5,0.5,-1), end (0.5,0.5,1), surface z=0 (parametric surface)
    {
        // 3D surface z = 0
        auto surface = Surface<double,3>(
            [](double u, double v) { return Point<double,3>({u, v, 0.0}); },
            std::make_pair(std::make_pair(-1.0, 1.0), std::make_pair(-1.0, 1.0))
        );

        Point<double,3> seg_start({0.5,0.5,-1});
        Point<double,3> seg_end({0.5,0.5,1});
        Eigen::Vector3d dir = seg_end.coords - seg_start.coords;
        Line l(seg_start, dir);
        std::cout << "DEBUG: 3D extra test - Line start: (" << seg_start[0] << "," << seg_start[1] << "," << seg_start[2] << "), ";
        std::cout << "direction: (" << dir[0] << "," << dir[1] << "," << dir[2] << ")\n";
        std::cout << "DEBUG: Parametric surface z=0\n";
        auto result = intersect(l, surface);
        std::cout << "DEBUG: Intersection result: intersects=" << result.intersects << ", #points=" << result.points.size() << "\n";
        if (result.intersects && !result.points.empty()) {
            const auto& p = result.points[0];
            std::cout << "DEBUG: Intersection point: (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
        }
        printTest("3D Line-Surface (segment start/end test): has intersection", result.intersects);
    }
    // 4D: start (0.5,0.5,-1,0), end (0.5,0.5,1,0), surface w=0 (parametric surface)
    {
        // 4D surface w = 0
        auto surface4 = Surface<double,4>(
            [](double u, double v) { return Point<double,4>({u, v, 0.0, 0.0}); },
            std::make_pair(std::make_pair(-1.0, 1.0), std::make_pair(-1.0, 1.0))
        );
        Point<double,4> seg_start4({0.5,0.5,-1,0});
        Point<double,4> seg_end4({0.5,0.5,1,0});
        Eigen::Matrix<double,4,1> dir4 = seg_end4.coords - seg_start4.coords;
        Line l4(seg_start4, dir4);
        std::cout << "DEBUG: 4D extra test - Line start: (" << seg_start4[0] << "," << seg_start4[1] << "," << seg_start4[2] << "," << seg_start4[3] << "), ";
        std::cout << "direction: (" << dir4[0] << "," << dir4[1] << "," << dir4[2] << "," << dir4[3] << ")\n";
        std::cout << "DEBUG: Parametric surface w=0\n";
        auto result4 = intersect(l4, surface4);
        std::cout << "DEBUG: Intersection result: intersects=" << result4.intersects << ", #points=" << result4.points.size() << "\n";
        if (result4.intersects && !result4.points.empty()) {
            const auto& p = result4.points[0];
            std::cout << "DEBUG: Intersection point: (" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << ")\n";
        }
        printTest("4D Line-Surface (segment start/end test): has intersection", result4.intersects);
    }
    // 5D: start (0.5,0.5,-1,0,0), end (0.5,0.5,1,0,0), surface z=0 (parametric surface)
    {
        // 5D surface z = 0
        auto surface5 = Surface<double,5>(
            [](double u, double v) { return Point<double,5>({u, v, 0.0, 0.0, 0.0}); },
            std::make_pair(std::make_pair(-1.0, 1.0), std::make_pair(-1.0, 1.0))
        );
        Point<double,5> seg_start5({0.5,0.5,-1,0,0});
        Point<double,5> seg_end5({0.5,0.5,1,0,0});
        Eigen::Matrix<double,5,1> dir5 = seg_end5.coords - seg_start5.coords;
        Line l5(seg_start5, dir5);
        std::cout << "DEBUG: 5D extra test - Line start: (" << seg_start5[0] << "," << seg_start5[1] << "," << seg_start5[2] << "," << seg_start5[3] << "," << seg_start5[4] << "), ";
        std::cout << "direction: (" << dir5[0] << "," << dir5[1] << "," << dir5[2] << "," << dir5[3] << "," << dir5[4] << ")\n";
        std::cout << "DEBUG: Parametric surface z=0\n";
        auto result5 = intersect(l5, surface5);
        std::cout << "DEBUG: Intersection result: intersects=" << result5.intersects << ", #points=" << result5.points.size() << "\n";
        if (result5.intersects && !result5.points.empty()) {
            const auto& p = result5.points[0];
            std::cout << "DEBUG: Intersection point: (" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << ", " << p[4] << ")\n";
        }
        printTest("5D Line-Surface (segment start/end test): has intersection", result5.intersects);
    }
}


inline void testLineSegmentIntersection() {
    // 1. Line intersects segment
    {
        Line l(Point<double,3>({0.5, 0.5, 0.5}), Eigen::Vector3d({1,0,0}));
        Point<double,3> a({0,0,0});
        Point<double,3> b({1,1,1});
        Segment s(a, b);
        auto result = intersect(l, s);
        printTest("Line-Segment: intersection exists", result.intersects);
        if (result.intersects) {
            printTest("Line-Segment: result is a Point", !result.points.empty());
            if (!result.points.empty()) {
                Point<double,3> p = result.points[0];
                printTest("Line-Segment: intersection at (0.5,0.5,0.5)", p == Point<double,3>({0.5,0.5,0.5}));
            }
        }
    }
    // 2. Line misses segment
    {
        Line l(Point<double,3>({2,2,2}), Eigen::Vector3d({1,0,0}));
        Point<double,3> a({0,0,0});
        Point<double,3> b({1,1,1});
        Segment s(a, b);
        auto result = intersect(l, s);
        printTest("Line-Segment: no intersection", !result.intersects);
    }
}

inline void testLineFaceIntersection() {
    // 1. Line intersects face (square in XY plane)
    {
        std::vector<Point<double,3>> vertices = {
            Point<double,3>({0,0,0}),
            Point<double,3>({1,0,0}),
            Point<double,3>({1,1,0}),
            Point<double,3>({0,1,0})
        };
        Eigen::Vector3d normal = (vertices[1].coords - vertices[0].coords).cross(vertices[2].coords - vertices[0].coords).normalized();
        Face<double,3> face(vertices[0], normal, vertices);
        Line l(Point<double,3>({0.5,0.5,1}), Eigen::Vector3d({0,0,-1}));
        auto result = intersect(l, face);
        printTest("Line-Face: intersection exists", result.intersects);
        if (result.intersects) {
            printTest("Line-Face: result is a Point", !result.points.empty());
            if (!result.points.empty()) {
                Point<double,3> p = result.points[0];
                printTest("Line-Face: intersection at (0.5,0.5,0)", p == Point<double,3>({0.5,0.5,0}));
            }
        }
    }
    // 2. Line misses face (parallel above)
    {
        std::vector<Point<double,3>> vertices = {
            Point<double,3>({0,0,0}),
            Point<double,3>({1,0,0}),
            Point<double,3>({1,1,0}),
            Point<double,3>({0,1,0})
        };
        Eigen::Vector3d normal = (vertices[1].coords - vertices[0].coords).cross(vertices[2].coords - vertices[0].coords).normalized();
        Face<double,3> face(vertices[0], normal, vertices);
        Line l(Point<double,3>({0.5,0.5,1}), Eigen::Vector3d({1,0,0}));
        auto result = intersect(l, face);
        printTest("Line-Face: no intersection", !result.intersects);
    }
}

inline void testLineMultiIntersectionSurface() {
    using namespace Euclid;
    using std::cout;

    // 1. Sphere (2 intersections)
    {
        auto sphere = Surface<double,3>(
                                        [](double u, double v){
                                            double r = 1.0;
                                            return Point<double,3>({
                                                r*std::sin(u)*std::cos(v),
                                                r*std::sin(u)*std::sin(v),
                                                r*std::cos(u)
                                            });
                                        },
                                        {{0.0, M_PI}, {0.0, 2*M_PI}}
                                        );
        Line l(Point<double,3>({0,0,-2}), Eigen::Vector3d({0,0,1}));
        auto result = intersect(l, sphere);
        
        std::cout << "DEBUG: Sphere intersection count: " << result.points.size() << "\n";
        
        for (size_t i=0;i<result.points.size();++i){
            const auto& p = result.points[i];
            std::cout << "DEBUG: Sphere point " << i << ": (" << p[0] << "," << p[1] << "," << p[2] << ")\n";
        }
        
        printTest("Sphere: line through center intersects twice", result.intersects && result.points.size() == 2);
        
    }

    // 2. Large Torus (multiple intersections through tube)
    {
        double R = 20.0, r = 5.0;
        auto torus = Surface<double,3>(
            [R,r](double u,double v){
                return Point<double,3>({
                    (R + r*std::cos(v))*std::cos(u),
                    (R + r*std::cos(v))*std::sin(u),
                    r*std::sin(v)
                });
            },
            {{0.0, 2*M_PI}, {0.0, 2*M_PI}}
        );
        // Line through torus tube along x-axis at y=0, z=0
        Line l(Point<double,3>({1.0,0,0}), Eigen::Vector3d({1,0,0}));
        auto result = intersect(l, torus);
        std::cout << "DEBUG: Torus intersection count: " << result.points.size() << "\n";
        for (size_t i=0;i<result.points.size();++i){
            const auto& p = result.points[i];
            std::cout << "DEBUG: Torus point " << i << ": (" << p[0] << "," << p[1] << "," << p[2] << ")\n";
        }
        printTest("Large Torus: line through tube intersects multiple times", result.intersects && result.points.size() >= 2);
    }
    
    // 2. Small Torus (multiple intersections through tube)
    {
        double R = 4.0, r = 0.5;
        auto torus = Surface<double,3>(
            [R,r](double u,double v){
                return Point<double,3>({
                    (R + r*std::cos(v))*std::cos(u),
                    (R + r*std::cos(v))*std::sin(u),
                    r*std::sin(v)
                });
            },
            {{0.0, 2*M_PI}, {0.0, 2*M_PI}}
        );
        // Line through torus tube along x-axis at y=0, z=0
        Line l(Point<double,3>({1.0,0,0}), Eigen::Vector3d({1,0,0}));
        auto result = intersect(l, torus);
        std::cout << "DEBUG: Torus intersection count: " << result.points.size() << "\n";
        for (size_t i=0;i<result.points.size();++i){
            const auto& p = result.points[i];
            std::cout << "DEBUG: Torus point " << i << ": (" << p[0] << "," << p[1] << "," << p[2] << ")\n";
        }
        printTest("Torus: line through tube intersects multiple times", result.intersects && result.points.size() >= 2);
    }
    
    // 2. Very small Torus (multiple intersections through tube)
    {
        double R = 0.1, r = 0.092;
        auto torus = Surface<double,3>(
            [R,r](double u,double v){
                return Point<double,3>({
                    (R + r*std::cos(v))*std::cos(u),
                    (R + r*std::cos(v))*std::sin(u),
                    r*std::sin(v)
                });
            },
            {{0.0, 2*M_PI}, {0.0, 2*M_PI}}
        );
        // Line through torus tube along x-axis at y=0, z=0
        Line l(Point<double,3>({1.0,0,0}), Eigen::Vector3d({1,0,0}));
        auto result = intersect(l, torus);
        std::cout << "DEBUG: Torus intersection count: " << result.points.size() << "\n";
        for (size_t i=0;i<result.points.size();++i){
            const auto& p = result.points[i];
            std::cout << "DEBUG: Torus point " << i << ": (" << p[0] << "," << p[1] << "," << p[2] << ")\n";
        }
        printTest("Very Small Torus: line through tube intersects multiple times", result.intersects && result.points.size() >= 2);
    }

    // 3. 3D sinusoidal surface (ribbon with multiple intersections)
    {
        using namespace Euclid;
        // Line along x-axis at y=0.0 (center of ribbon)
        Line l(Point<double,3>({0,0.0,0}), Eigen::Vector3d({1.0,0.0,0.0}));

        // Proper 2D ribbon surface: width 0.5, v spreads across y axis
        auto surface3D = Surface<double,3>(
            [](double u,double v){
                double amplitude = 1.0;
                double frequency = 4.0;
                double width = 1;
                Point<double,3> p{u, (v - 0.5) * 2.0 * width, amplitude * std::sin(frequency*M_PI*u)};
                //std::cout << "DEBUG: surface3D(u=" << u << ", v=" << v << ") = ("
                //          << p[0] << "," << p[1] << "," << p[2] << ")\n";
                return p;
            },
            {{0.0, 4.0}, {0.0, 1.0}}
        );
        // Line-surface intersection
        auto result = intersect(l, surface3D);
        std::cout << "DEBUG: 3D sinusoidal ribbon surface intersection count: " << result.points.size() << "\n";
        for (size_t i=0;i<result.points.size();++i){
            const auto& p = result.points[i];
            std::cout << "DEBUG: 3D surface point " << i << ": (" << p[0] << "," << p[1] << "," << p[2] << ")\n";
        }
        printTest("3D sinusoidal ribbon surface: line intersects multiple peaks", result.intersects && result.points.size() >= 2);

        // Export mesh for inspection (high resolution, 2D surface)
        std::string filename = "3d_sinusoidal_ribbon.obj";
        int M = 400; // samples along u (higher resolution)
        int N = 40;  // ribbon width resolution along v (higher resolution)
        auto mesh = generateSurfaceMesh(surface3D, M, N);
        if (auto objData = mesh.exportOBJ()) {
            std::ofstream file(filename);
            file << *objData;
            file.close();
            std::cout << "Exported 3d_sinusoidal_ribbon.obj successfully.\n";
        }
        // Also output a simple .obj for intersection points as a point cloud for inspection
        if (!result.points.empty()) {
            std::ofstream pts("3d_sinusoidal_intersections.obj");
            for (const auto& p : result.points) {
                pts << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
            }
            pts.close();
            std::cout << "Exported 3d_sinusoidal_intersections.obj (intersection points).\n";
        }
    }

    // 3a. 3D high-frequency sinusoidal ribbon surface (multiple intersections)
    {
        using namespace Euclid;
        // Line along x-axis at y=0.0 (center of ribbon)
        Line l(Point<double,3>({0,0.0,0}), Eigen::Vector3d({1.0,0.0,0.0}));

        // 2D ribbon surface, high frequency = 8, width = 1, v spreads across y axis
        auto surface3D_highfreq = Surface<double,3>(
            [](double u,double v){
                double amplitude = 1.0;
                double frequency = 8.0;
                double width = 1;
                Point<double,3> p{u, (v - 0.5) * 2.0 * width, amplitude * std::sin(frequency*M_PI*u)};
                return p;
            },
            {{0.0, 4.0}, {0.0, 1.0}}
        );
        // Line-surface intersection
        auto result = intersect(l, surface3D_highfreq);
        std::cout << "DEBUG: 3D high-frequency sinusoidal ribbon surface intersection count: " << result.points.size() << "\n";
        for (size_t i=0;i<result.points.size();++i){
            const auto& p = result.points[i];
            std::cout << "DEBUG: 3D highfreq surface point " << i << ": (" << p[0] << "," << p[1] << "," << p[2] << ")\n";
        }
        printTest("3D high-frequency sinusoidal ribbon surface: line intersects multiple peaks", result.intersects && result.points.size() >= 2);

        // Export mesh for inspection (high resolution, 2D surface)
        std::string filename = "3d_sinusoidal_ribbon_highfreq.obj";
        int M = 400; // samples along u (higher resolution)
        int N = 40;  // ribbon width resolution along v (higher resolution)
        auto mesh = generateSurfaceMesh(surface3D_highfreq, M, N);
        if (auto objData = mesh.exportOBJ()) {
            std::ofstream file(filename);
            file << *objData;
            file.close();
            std::cout << "Exported 3d_sinusoidal_ribbon_highfreq.obj successfully.\n";
        }
        // Output intersection points as a point cloud for inspection (optional)
        if (!result.points.empty()) {
            std::ofstream pts("3d_sinusoidal_highfreq_intersections.obj");
            for (const auto& p : result.points) {
                pts << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
            }
            pts.close();
            std::cout << "Exported 3d_sinusoidal_highfreq_intersections.obj (intersection points).\n";
        }
    }

    // 3b. 3D even higher-frequency sinusoidal ribbon surface (frequency = 16, 32, 64, 128)
    for (double freq : {16.0, 32.0, 64.0, 128.0}) {
        using namespace Euclid;
        Line l(Point<double,3>({0,0.0,0}), Eigen::Vector3d({1.0,0.0,0.0}));
        auto surface3D = Surface<double,3>(
            [freq](double u, double v) {
                double amplitude = 1.0;
                double width = 1.0;
                return Point<double,3>({u, (v - 0.5) * 2.0 * width, amplitude * std::sin(freq * M_PI * u)});
            },
            {{0.0, 4.0}, {0.0, 1.0}}
        );
        auto result = intersect(l, surface3D);
        std::cout << "DEBUG: 3D sinusoidal ribbon freq=" << freq
                  << " intersection count: " << result.points.size() << "\n";
        size_t n_show = std::min<size_t>(result.points.size(), 5);
        for (size_t i=0;i<n_show;++i){
            const auto& p = result.points[i];
            std::cout << "  freq=" << freq << " pt" << i << ": ("
                      << p[0] << "," << p[1] << "," << p[2] << ")\n";
        }
        if (result.points.size() > n_show)
            std::cout << "  ... (" << result.points.size() - n_show << " more points)\n";
        printTest(("3D sinusoidal ribbon freq=" + std::to_string(int(freq)) + ": intersects").c_str(), result.intersects && result.points.size() >= 2);
        // No OBJ export for these high-frequency tests
    }

    // 4. 4D sinusoidal surface (ribbon with multiple intersections)
    {
        using namespace Euclid;

        // Line along x-axis at y=0.25, w=0
        Line l4(Point<double,4>({0,0.25,0,0}), Eigen::Matrix<double,4,1>({1.0,0.0,0.0,0.0}));

        // Ribbon surface: small thickness along v
        auto surface4D = Surface<double,4>(
            [](double u,double v){
                double amplitude = 1.0;
                double frequency = 4.0;      // number of peaks
                double width = 0.05;         // ribbon half-width
                return Point<double,4>({
                    u,
                    0.25 + width * (v - 0.5),    // ribbon along y
                    amplitude * std::sin(frequency*M_PI*u),
                    0.0                          // keep w=0 for simplicity
                });
            },
            {{0.0, 4.0}, {0.0, 1.0}}        // u and v domains
        );

        // Line-surface intersection
        auto result4 = intersect(l4, surface4D);
        std::cout << "DEBUG: 4D sinusoidal ribbon surface intersection count: " << result4.points.size() << "\n";
        for (size_t i=0;i<result4.points.size();++i){
            const auto& p = result4.points[i];
            std::cout << "DEBUG: 4D surface point " << i << ": ("
                      << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << ")\n";
        }
        printTest("4D sinusoidal ribbon surface: line intersects multiple peaks", result4.intersects && result4.points.size() >= 2);

        // Export mesh for inspection (project 4D -> 3D by dropping w)
        std::string filename = "4d_sinusoidal_ribbon.obj";
        int M = 100; // samples along u
        int N = 10;  // ribbon thickness resolution along v

        // Simple projection lambda: drop w
        auto projectedSurface = Surface<double,3>(
            [&surface4D](double u,double v){
                Point<double,4> p4 = surface4D.evaluate(u,v);   // <-- use evaluate() instead of operator()
                return Point<double,3>({p4[0], p4[1], p4[2]});
            },
            {{0.0, 4.0}, {0.0, 1.0}}   // use same domain
        );

        auto mesh = generateSurfaceMesh(projectedSurface, M, N);
        if (auto objData = mesh.exportOBJ()) {
            std::ofstream file(filename);
            file << *objData;
            file.close();
            std::cout << "Exported 4d_sinusoidal_ribbon.obj successfully.\n";
        }
    }

 
}

inline void testLineIntersection() {
    using namespace Euclid;
    using std::cout;
    cout << "\nLine Intersection Tests\n" << std::endl;
    
    //testLineLineIntersection();
    //testLineSegmentIntersection();
    testLineCurveIntersection();
    //testLinePlaneIntersection();
    //testLinePlaneIntersectionND();
    //testLineFaceIntersection();
    //testLineSurfaceIntersection();
    //testLineMultiIntersectionSurface();
}

} // namespace Tests
} // namespace Euclid
