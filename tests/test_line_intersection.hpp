#pragma once


#include "geometry/geometry.hpp"
#include "geometry/intersection/line_intersection.hpp"
#include "test_utilities.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>

inline std::string currentTimestamp() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

namespace Euclid {
namespace Tests {

// Frequency resolution test for line–surface intersection (sinusoidal surface)
inline void testLineRibbonSurfaceFrequencyResolution() {
    using namespace Euclid;
    using std::cout;

    cout << "---- Line–Surface Frequency Resolution Test ----" << std::endl;

    // Log test start time and date
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    cout << "🕒 Test started at: " << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S") << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    // Use Euclid's tolerance model
    Tolerance tol(1e-6, 1e-9, 1e-6, 1.0);
    double expectedLimit = 1.0 / (2.0 * M_PI * tol.paramTol);
    cout << "Expected frequency limit ≈ " << expectedLimit
         << " (based on paramTol = " << tol.paramTol << ")\n";

    Line<double,3> line(Point<double,3>({0, 0, 0}), Eigen::Vector3d({1, 0, 0}));
    double freq = 1.0;
    [[maybe_unused]] double failFreq = -1.0;
    [[maybe_unused]] int iteration_count = 0;

    // For profiling and scaling analysis
    std::vector<double> freqs, times, avg_times;
    std::vector<int> expected_vals, found_vals;
    std::vector<double> alphas;

    cout << std::setw(10) << "freq"
         << std::setw(12) << "expected"
         << std::setw(10) << "found"
         << std::setw(10) << "status"
         << std::setw(15) << "runtime (ms)"
         << std::setw(15) << "avg (ms/pt)"
         << std::setw(10) << "α"
         << std::endl;

    double prev_time = 0.0;
    double prev_freq = 0.0;
    bool first = true;
    while (true) {
        auto surface = Surface<double,3>(
            [freq](double u, double v) {
                return Point<double,3>({
                    u,
                    v - 0.5, // small strip around y=0
                    std::sin(2.0 * M_PI * freq * u)
                });
            },
            {{0.0, 1.0}, {0.0, 1.0}}
        );

        int expected = static_cast<int>(2 * freq + 1);
        auto t0 = std::chrono::high_resolution_clock::now();
        auto result = intersect(line, surface, tol);
        auto t1 = std::chrono::high_resolution_clock::now();
        int found = static_cast<int>(result.points.size());
        bool ok = (found == expected);
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        double avg_ms = found > 0 ? ms / found : 0.0;

        freqs.push_back(freq);
        times.push_back(ms);
        avg_times.push_back(avg_ms);
        expected_vals.push_back(expected);
        found_vals.push_back(found);

        double alpha = 0.0;
        if (!first && prev_time > 0.0 && prev_freq > 0.0 && ms > 0.0) {
            alpha = std::log(ms / prev_time) / std::log(freq / prev_freq);
            alphas.push_back(alpha);
        }
        cout << std::setw(10) << freq
             << std::setw(12) << expected
             << std::setw(10) << found
             << std::setw(8)
             << (found == expected ? "[OK]" : (std::abs(found - expected) <= 1 ? "[NEAR]" : "[!!]"))
             << std::setw(15) << std::fixed << std::setprecision(3) << ms
             << std::setw(15) << std::fixed << std::setprecision(6) << avg_ms;
        if (!first && prev_time > 0.0 && prev_freq > 0.0 && ms > 0.0)
            cout << std::setw(10) << std::fixed << std::setprecision(3) << alpha;
        else
            cout << std::setw(10) << "-";
        cout << std::endl;

        ++iteration_count;

        if (!ok) {
            failFreq = freq;
            cout << "[FAIL] Numerical resolution limit reached at freq = " << freq
                 << "\nExpected theoretical limit ≈ " << expectedLimit
                 << "\nRatio actual/theoretical = " << freq / expectedLimit << std::endl;
            break;
        }

        prev_time = ms;
        prev_freq = freq;
        first = false;

        freq *= 2.0; // double each iteration

        if (freq > 1e5) {
            cout << "[PASS] No mismatch detected up to freq = " << freq << std::endl;
            break;
        }
    }
    // Compute average scaling exponent α
    double alpha_sum = 0.0;
    int alpha_count = 0;
    for (double a : alphas) {
        if (std::isfinite(a)) {
            alpha_sum += a;
            ++alpha_count;
        }
    }
    if (alpha_count > 0) {
        double avg_alpha = alpha_sum / alpha_count;
        cout << "Average scaling exponent α ≈ " << std::fixed << std::setprecision(4) << avg_alpha << std::endl;
    } else {
        cout << "Average scaling exponent α: not available (insufficient data)" << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    auto end_time = std::chrono::system_clock::now();
    std::time_t end_c = std::chrono::system_clock::to_time_t(end_time);
    cout << "✅ Test completed at: " << std::put_time(std::localtime(&end_c), "%Y-%m-%d %H:%M:%S")
         << " | Runtime: " << duration.count() << " seconds\n";
}

// Scaling test: empirical vs theoretical frequency limit across paramTol values
inline void testLineSurfaceFrequencyScaling() {
    using namespace Euclid;
    using std::cout;

    cout << "---- Line–Surface Frequency Scaling Test ----" << std::endl;

    // Log test start time and date
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    cout << "🕒 Test started at: " << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S") << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<double> paramTols = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
    cout << std::setw(10) << "paramTol" << std::setw(15) << "f_theory"
         << std::setw(15) << "f_empirical" << std::setw(15) << "ratio" << std::endl;

    int iteration_count = 0;
    for (double paramTol : paramTols) {
        Tolerance tol(1e-6, 1e-9, paramTol, 1.0);
        double expectedLimit = 1.0 / (2.0 * M_PI * tol.paramTol);
        double freq = 1.0;
        double failFreq = -1.0;

        Line<double,3> line(Point<double,3>({0, 0, 0}), Eigen::Vector3d({1, 0, 0}));

        while (freq < expectedLimit * 10) {
            const double thisFreq = freq; // isolate capture
            Tolerance localTol(1e-6, 1e-9, paramTol, 1.0);
            Line<double,3> line(Point<double,3>({0, 0, 0}), Eigen::Vector3d({1, 0, 0}));

            Surface<double,3> surface(
                [thisFreq](double u, double v) {
                    return Point<double,3>({
                        u,
                        v - 0.5,
                        std::sin(2.0 * M_PI * thisFreq * u)
                    });
                },
                {{0.0, 1.0}, {0.0, 1.0}}
            );

            int expected = static_cast<int>(2 * thisFreq + 1);
            auto result = intersect(line, surface, localTol);
            int found = static_cast<int>(result.points.size());
            bool ok = std::abs(found - expected) <= 1;

            ++iteration_count;

            if (!ok) {
                failFreq = thisFreq;
                break;
            }

            freq *= 2.0;
        }

        double empiricalLimit = (failFreq > 0 ? failFreq : freq);
        double ratio = empiricalLimit / expectedLimit;

        cout << std::setw(10) << std::scientific << paramTol
             << std::setw(15) << std::fixed << std::setprecision(0) << expectedLimit
             << std::setw(15) << std::fixed << std::setprecision(0) << empiricalLimit
             << std::setw(15) << std::setprecision(6) << ratio
             << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    auto end_time = std::chrono::system_clock::now();
    std::time_t end_c = std::chrono::system_clock::to_time_t(end_time);
    cout << "✅ Test completed at: " << std::put_time(std::localtime(&end_c), "%Y-%m-%d %H:%M:%S")
         << " | Runtime: " << duration.count() << " seconds\n";
    cout << "⚙️  Approx complexity ~ O(" << iteration_count << " intersections)\n";
    cout << "✅ Completed frequency scaling sweep across paramTol values.\n";
}

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

    // 1b. Oblique intersection (non-perpendicular line)
    {
        Line l(Point<double,3>({0,0,0}), Eigen::Vector3d({1,1,1}));
        Plane pl(Point<double,3>({0,0,5}), Eigen::Vector3d({0,0,1}));
        auto result = intersect(l, pl);
        printTest("Line-Plane (oblique): has intersection", result.intersects);
        if (result.intersects && !result.points.empty()) {
            const auto& p = result.points[0];
            printTest("Line-Plane (oblique): z≈5", std::abs(p[2] - 5.0) < 1e-9);
        }
    }

    // 1c. Parallel non-coincident (above plane)
    {
        Line l(Point<double,3>({0,0,1}), Eigen::Vector3d({1,0,0}));
        Plane pl(Point<double,3>({0,0,0}), Eigen::Vector3d({0,0,1}));
        auto result = intersect(l, pl);
        printTest("Line-Plane (parallel above): no intersection", !result.intersects);
    }

    // 1d. Coincident case
    {
        Line l(Point<double,3>({0,0,5}), Eigen::Vector3d({1,0,0}));
        Plane pl(Point<double,3>({0,0,5}), Eigen::Vector3d({0,0,1}));
        auto result = intersect(l, pl);
        std::cout << "DEBUG: Coincident Test -- dot(line_dir, plane_normal)="
                  << l.direction().dot(pl.normal) << ", offset="
                  << (pl.base.coords - l.point1().coords).dot(pl.normal)
                  << std::endl;
        printTest("Line-Plane (coincident): has intersection", result.intersects);
        printTest("Line-Plane (coincident): result is a Line", !result.lines.empty());
    }


    // 1f. Scaled geometry (tolerance consistency)
    {
        double scale = 1e6;
        Line l(Point<double,3>({0,0,0}), Eigen::Vector3d({0,0,scale}));
        Plane pl(Point<double,3>({0,0,5*scale}), Eigen::Vector3d({0,0,1}));
        auto result = intersect(l, pl);
        printTest("Line-Plane (scaled units): has intersection", result.intersects);
        if (result.intersects && !result.points.empty()) {
            printTest("Line-Plane (scaled units): z≈5e6", std::abs(result.points[0][2] - 5*scale) < 1e-3);
        }
    }
}

inline void testLinePlaneParallelThreshold() {
    using namespace Euclid;
    using std::cout;
    cout << "---- Line–Plane Parallelism Logarithmic Sweep (micro angles in radians, forward) ----" << std::endl;

    Plane<double,3> plane(Point<double,3>({0,0,0}), Eigen::Vector3d({0,0,1}));
    double z0 = 1.0;
    // Logarithmic sweep: start at 1e-10 rad (nearly parallel), up to 1e-1 rad (steep)
    const int N = 40;
    double minAngle = 1e-10;
    double maxAngle = 1e-1;
    double logMin = std::log10(minAngle);
    double logMax = std::log10(maxAngle);
    double lastAngle = -1.0;
    bool found_transition = false;
    double first_success_angle = -1.0;
    for (int i = 0; i <= N; ++i) {
        double logAngle = logMin + (logMax - logMin) * (double(i) / N);
        double theta = std::pow(10.0, logAngle);
        // Avoid duplicate theta for i=0 and i=N due to floating point
        if (theta == lastAngle) continue;
        lastAngle = theta;
        // Direction: (cos(theta), 0, sin(theta)), so theta ≈ 0 is parallel to plane
        Eigen::Vector3d dir(std::cos(theta), 0, std::sin(theta));
        Line<double,3> line(Point<double,3>({0,0,z0}), dir);

        auto result = intersect(line, plane);
        bool hit = result.intersects;

        double dotVal = std::abs(dir.normalized().dot(plane.normal));
        std::cout << "angle (rad)=" << std::scientific << std::setprecision(2) << theta
                  << ", dot=" << std::setprecision(10) << dotVal
                  << "  -> intersects=" << (hit ? "✅" : "❌") << std::endl;

        // Look for first transition from ❌ to ✅ (first hit after a miss)
        if (!found_transition && hit) {
            first_success_angle = theta;
            found_transition = true;
            // Optionally, break here if you want only the first transition
            break;
        }
    }

    if (found_transition)
        std::cout << "✅ Intersection first succeeds at angle ≈ " << std::scientific << std::setprecision(2) << first_success_angle << " rad from parallel\n";
    else
        std::cout << "❌ No intersection transition detected (never succeeded) in tested angle range\n";
}

// Logarithmic sweep of micro angles for line-segment intersection parallel threshold
inline void testLineSegmentParallelThreshold() {
    using namespace Euclid;
    using std::cout;
    cout << "---- Line–Segment Parallelism Logarithmic Sweep (micro angles in radians, backward from intersecting to parallel) ----" << std::endl;

    // Base segment along X-axis, z=0
    Point<double,3> A({0, 0, 0});
    Point<double,3> B({1, 0, 0});
    Segment<double,3> seg(A, B);

    const int N = 80;
    double minAngle = 1e-10;
    double maxAngle = 1e-1;
    double logMin = std::log10(minAngle);
    double logMax = std::log10(maxAngle);
    double lastAngle = -1.0;
    bool found_hit = false;
    bool lost_transition = false;
    double first_failure_angle = -1.0;

    Eigen::Vector3d axisY(0, 1, 0); // rotate around Y axis to tilt downward

    // Sweep from steep (logMax) to nearly parallel (logMin)
    for (int i = 0; i <= N; ++i) {
        double logAngle = logMax - (logMax - logMin) * (double(i) / N);
        double theta = std::pow(10.0, logAngle);
        if (theta == lastAngle) continue;
        lastAngle = theta;

        // Compute dynamic offset proportional to actual vertical displacement
        double segment_length = (B.coords - A.coords).norm();
        // Start with negative z0 so the initial configuration is below the segment (ensures intersection at steep angles)
        double z0 = -0.9 * segment_length * std::sin(theta);  // proportional to actual vertical displacement

        // construct line using same segment endpoints but tilted around Y axis
        Eigen::AngleAxisd rot(-theta, axisY);
        Eigen::Vector3d A_off = A.coords + Eigen::Vector3d(0, 0, z0);
        Eigen::Vector3d B_off = B.coords + Eigen::Vector3d(0, 0, z0);
        Eigen::Vector3d A_rot = rot * A_off;
        Eigen::Vector3d B_rot = rot * B_off;

        // build line directly through tilted segment endpoints
        Line<double,3> line{Point<double,3>(A_rot), Point<double,3>(B_rot)};

        auto result = intersect(line, seg);
        bool hit = result.intersects && !result.points.empty();
        if (result.intersects && result.points.empty() && !result.lines.empty()) {
            std::cout << "  (coincident case detected)\n";
            if (found_hit && !lost_transition) {
                first_failure_angle = theta;
                lost_transition = true;
            }
            break; // stop after marking transition
        }

        // Compute minimal distance between line and segment
        Eigen::Vector3d lineDir = line.direction().normalized();
        Eigen::Vector3d segDir = seg.end.coords - seg.start.coords;
        Eigen::Vector3d w0 = line.point1().coords - seg.start.coords;
        double a = lineDir.dot(lineDir);
        double b = lineDir.dot(segDir);
        double c = segDir.dot(segDir);
        double d = lineDir.dot(w0);
        double e = segDir.dot(w0);
        double denom = a*c - b*b;
        double sc, tc;
        if (std::abs(denom) < 1e-15) {
            sc = 0.0;
            tc = (b > c ? d/b : e/c);
        } else {
            sc = (b*e - c*d) / denom;
            tc = (a*e - b*d) / denom;
        }
        Eigen::Vector3d dP = w0 + sc*lineDir - tc*segDir;
        double dist = dP.norm();

        double dotVal = std::abs(lineDir.dot(segDir.normalized()));
        std::cout << "angle (rad)=" << std::scientific << std::setprecision(2) << theta
                  << ", dot=" << std::setprecision(10) << dotVal
                  << ", dist=" << std::setprecision(4) << dist
                  << "  -> intersects=" << (hit ? "✅" : "❌") << std::endl;

        // Track transition: after seeing any hit, record first non-hit as the loss point
        if (!found_hit && hit) {
            found_hit = true;
        }
        if (found_hit && !hit && !lost_transition) {
            first_failure_angle = theta;
            lost_transition = true;
            // Do not break; keep printing all values for the sweep
        }
    }

    if (lost_transition)
        std::cout << "✅ Line–Segment intersection first lost (from intersecting to parallel) at angle ≈ "
                  << std::scientific << std::setprecision(2) << first_failure_angle
                  << " rad from parallel\n";
    else
        std::cout << "❌ No loss-of-intersection transition detected (always intersected)\n";
}

// Logarithmic sweep of micro angles for line-face intersection parallel threshold
inline void testLineFaceParallelThreshold() {
    using namespace Euclid;
    using std::cout;
    cout << "---- Line–Face Parallelism Logarithmic Sweep (micro angles in radians, backward from intersecting to parallel) ----" << std::endl;

    // Define square face in XY plane, z=0, normal (0,0,1)
    std::vector<Point<double,3>> vertices = {
        Point<double,3>({0,0,0}),
        Point<double,3>({1,0,0}),
        Point<double,3>({1,1,0}),
        Point<double,3>({0,1,0})
    };
    Eigen::Vector3d normal = (vertices[1].coords - vertices[0].coords).cross(vertices[2].coords - vertices[0].coords).normalized();
    Face<double,3> face(vertices[0], normal, vertices);

    const int N = 80;
    double minAngle = 1e-10;
    double maxAngle = 1e-1;
    double logMin = std::log10(minAngle);
    double logMax = std::log10(maxAngle);
    double lastAngle = -1.0;
    bool found_hit = false;
    bool lost_transition = false;
    double first_failure_angle = -1.0;

    Eigen::Vector3d axisY(0, 1, 0); // rotate about Y to tilt toward parallel

    for (int i = 0; i <= N; ++i) {
        double logAngle = logMax - (logMax - logMin) * (double(i) / N);
        double theta = std::pow(10.0, logAngle);
        if (theta == lastAngle) continue;
        lastAngle = theta;

        double face_size = 1.0;
        double z0 = -0.9 * face_size * std::sin(theta);

        // Construct line tilted toward the face
        Eigen::AngleAxisd rot(-theta, axisY);
        Eigen::Vector3d origin(0.5, 0.5, z0);
        Eigen::Vector3d dir = rot * Eigen::Vector3d(0, 0, 1);
        Line<double,3> line(Point<double,3>(origin), dir);

        auto result = intersect(line, face);
        bool hit = result.intersects && !result.points.empty();

        if (result.intersects && result.points.empty() && !result.lines.empty()) {
            std::cout << "  (coincident case detected)\n";
            if (found_hit && !lost_transition) {
                first_failure_angle = theta;
                lost_transition = true;
            }
            break;
        }

        // Compute minimal distance from line origin to face plane
        double dist = std::abs(origin.dot(normal));
        double dotVal = std::abs(dir.normalized().dot(normal));

        std::cout << "angle (rad)=" << std::scientific << std::setprecision(2) << theta
                  << ", dot=" << std::setprecision(10) << dotVal
                  << ", dist=" << std::setprecision(4) << dist
                  << "  -> intersects=" << (hit ? "✅" : "❌") << std::endl;

        if (!found_hit && hit)
            found_hit = true;
        if (found_hit && !hit && !lost_transition) {
            first_failure_angle = theta;
            lost_transition = true;
            break; // stop sweeping after first loss
        }
    }

    if (lost_transition)
        std::cout << "✅ Line–Face intersection first lost (from intersecting to parallel) at angle ≈ "
                  << std::scientific << std::setprecision(2) << first_failure_angle
                  << " rad from parallel\n";
    else
        std::cout << "❌ No loss-of-intersection transition detected (always intersected)\n";
}

// Logarithmic sweep of micro angles for line-line intersection parallel threshold
inline void testLineLineParallelThreshold() {
    using namespace Euclid;
    using std::cout;
    cout << "---- Line–Line Parallelism Logarithmic Sweep (micro angles in radians, forward) ----" << std::endl;

    const int N = 40;
    double z0 = 1.0;
    double minAngle = 1e-10;
    double maxAngle = 1e-1;
    double logMin = std::log10(minAngle);
    double logMax = std::log10(maxAngle);
    double lastAngle = -1.0;
    bool found_transition = false;
    double first_success_angle = -1.0;

    for (int i = 0; i <= N; ++i) {
        double logAngle = logMin + (logMax - logMin) * (double(i) / N);
        double theta = std::pow(10.0, logAngle);
        if (theta == lastAngle) continue;
        lastAngle = theta;

        // Base line along x-axis
        Line<double,3> base(Point<double,3>({0,0,0}), Eigen::Vector3d({1,0,0}));
        // Second line rotated slightly upward by theta, starting at z=z0
        Eigen::Vector3d dir(std::cos(theta), 0, std::sin(theta));
        Line<double,3> l(Point<double,3>({0,0,z0}), dir);

        auto result = intersect(base, l);
        bool hit = result.intersects;

        double dotVal = std::abs(dir.normalized().dot(base.direction().normalized()));
        std::cout << "angle (rad)=" << std::scientific << std::setprecision(2) << theta
                  << ", dot=" << std::setprecision(10) << dotVal
                  << "  -> intersects=" << (hit ? "✅" : "❌") << std::endl;

        if (!found_transition && hit) {
            first_success_angle = theta;
            found_transition = true;
            break;
        }
    }

    if (found_transition)
        std::cout << "✅ Line-Line intersection first succeeds at angle ≈ " << std::scientific << std::setprecision(2) << first_success_angle << " rad from parallel\n";
    else
        std::cout << "❌ No intersection transition detected (never succeeded) in tested angle range\n";
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
    
    // 4D sinusoidal surface (ribbon with multiple intersections)
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
    
    // Sphere (2 intersections)
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

    // Large Torus (multiple intersections through tube)
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


inline void testLineTorusScaleResolution() {
    using namespace Euclid::Geometry;
    using T = double;
    constexpr int N = 3;

    std::cout << "---- Line–Torus Scale Resolution Test ----\n";
    std::cout << "🕒 Test started at: " << currentTimestamp() << "\n";

    // Fixed line through origin
    Line<T, N> line(Point<T, N>({0.0, 0.0, 0.0}), Eigen::Matrix<T, N, 1>(1.0, 0.0, 0.0));

    // Torus definition parameters
    const T R0 = 10.0;  // major radius
    const T r0 = 2.5;   // minor radius
    const T shrinkFactor = 0.8; // geometric shrink step

    int expected = 4;
    int found = 0;
    int step = 0;
    while (true) {
        T R = R0 * std::pow(shrinkFactor, step);
        T r = r0 * std::pow(shrinkFactor, step);

        Surface<T, N> torus(
            [R, r](double u, double v) {
                return Point<T, N>({
                    (R + r * std::cos(v)) * std::cos(u),
                    (R + r * std::cos(v)) * std::sin(u),
                    r * std::sin(v)
                });
            },
            {{0.0, 2 * M_PI}, {0.0, 2 * M_PI}}
        );
        auto t0 = std::chrono::high_resolution_clock::now();
        auto res = intersect(line, torus);
        auto t1 = std::chrono::high_resolution_clock::now();

        found = static_cast<int>(res.points.size());
        std::cout << std::scientific << std::setprecision(8) << R << " " << r
                  << "  expected=" << std::setw(3) << expected
                  << "  found=" << std::setw(3) << found
                  << "  runtime(ms)=";
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << std::fixed << std::setprecision(3) << std::setw(10) << ms << "\n";

        // Failure if found < expected
        if (found < expected) {
            std::cout << "[FAIL] Resolution limit reached at R="
                      << std::scientific << std::setprecision(8) << R
                      << ", r=" << std::scientific << std::setprecision(8) << r << "\n";
            break;
        }
        // Failure if found > expected * 2
        if (found > expected * 2) {
            std::cout << "[FAIL] Excess intersections detected at R="
                      << std::scientific << std::setprecision(8) << R
                      << ", r=" << std::scientific << std::setprecision(8) << r
                      << " (found=" << found << ", expected=" << expected << ")\n";
            break;
        }
        ++step;
    }
    std::cout << "✅ Test completed at: " << currentTimestamp() << "\n";
}

inline void testLineSphereScaleResolution() {
    using namespace Euclid::Geometry;
    using T = double;
    constexpr int N = 3;

    std::cout << "---- Line–Sphere Scale Resolution Test ----\n";
    std::cout << "🕒 Test started at: " << currentTimestamp() << "\n";

    Line<T, N> line(Point<T, N>({0.0, 0.0, -2.0}), Eigen::Matrix<T, N, 1>(0.0, 0.0, 1.0));

    const T R0 = 1.0;
    const T shrinkFactor = 0.8;
    int expected = 2;
    int found = 0;
    int step = 0;

    while (true) {
        T R = R0 * std::pow(shrinkFactor, step);

        Surface<T, N> sphere(
            [R](double u, double v) {
                return Point<T, N>({
                    R * std::sin(u) * std::cos(v),
                    R * std::sin(u) * std::sin(v),
                    R * std::cos(u)
                });
            },
            {{0.0, M_PI}, {0.0, 2 * M_PI}}
        );

        auto t0 = std::chrono::high_resolution_clock::now();
        auto res = intersect(line, sphere);
        auto t1 = std::chrono::high_resolution_clock::now();

        found = static_cast<int>(res.points.size());
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << std::scientific << std::setprecision(8)
                  << R << "  expected=" << std::setw(3) << expected
                  << "  found=" << std::setw(3) << found
                  << "  runtime(ms)=" << std::fixed << std::setprecision(3) << std::setw(10) << ms << "\n";

        if (found < expected) {
            std::cout << "[FAIL] Resolution limit reached at R=" << std::scientific << R << "\n";
            break;
        }
        if (found > expected * 2) {
            std::cout << "[FAIL] Excess intersections detected at R=" << std::scientific << R
                      << " (found=" << found << ", expected=" << expected << ")\n";
            break;
        }

        ++step;
    }

    std::cout << "✅ Test completed at: " << currentTimestamp() << "\n";
}

inline void testLineSaddleScaleResolution() {
    using namespace Euclid::Geometry;
    using T = double;
    constexpr int N = 3;

    std::cout << "---- Line–Saddle Scale Resolution Test ----\n";
    std::cout << "🕒 Test started at: " << currentTimestamp() << "\n";

    // Line passing through both sides of the saddle, yielding two intersections within the domain
    Line<T, N> line(Point<T, N>({-1.0, 0.0, -1.0}), Point<T, N>({1.0, 0.0, 1.0}));

    const T scale0 = 1.0;
    const T shrinkFactor = 0.8;
    int expected = 2; // usually 2 crossings
    int found = 0;
    int step = 0;

    while (true) {
        T s = scale0 * std::pow(shrinkFactor, step);

        Surface<T, N> saddle(
            [s](double u, double v) {
                return Point<T, N>({
                    s * u,
                    s * v,
                    s * (u*u - v*v)
                });
            },
            {{-1.0, 1.0}, {-1.0, 1.0}}
        );

        auto t0 = std::chrono::high_resolution_clock::now();
        auto res = intersect(line, saddle);
        auto t1 = std::chrono::high_resolution_clock::now();

        found = static_cast<int>(res.points.size());
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << std::scientific << std::setprecision(8)
                  << s << "  expected=" << std::setw(3) << expected
                  << "  found=" << std::setw(3) << found
                  << "  runtime(ms)=" << std::fixed << std::setprecision(3) << std::setw(10) << ms << "\n";

        if (found < expected) {
            std::cout << "[FAIL] Resolution limit reached at scale=" << std::scientific << s << "\n";
            break;
        }
        if (found > expected * 2) {
            std::cout << "[FAIL] Excess intersections detected at scale=" << std::scientific << s
                      << " (found=" << found << ", expected=" << expected << ")\n";
            break;
        }

        ++step;
    }

    std::cout << "✅ Test completed at: " << currentTimestamp() << "\n";
}

inline void testLineIntersection() {
    using namespace Euclid;
    using std::cout;
    cout << "\nLine Intersection Tests\n" << std::endl;
    
    testLineLineIntersection();
    testLineLineParallelThreshold();
    testLineSegmentIntersection();
    testLineSegmentParallelThreshold();
    testLineCurveIntersection();
    testLinePlaneIntersection();
    testLinePlaneIntersectionND();
    testLinePlaneParallelThreshold();
    testLineFaceIntersection();
    testLineFaceParallelThreshold();
    testLineSurfaceIntersection();
    testLineRibbonSurfaceFrequencyResolution();
    testLineTorusScaleResolution();
    testLineSphereScaleResolution();
    testLineSaddleScaleResolution();
    testLineSurfaceFrequencyScaling();
}

} // namespace Tests
} // namespace Euclid
