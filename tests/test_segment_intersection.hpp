// Modular submethod for segment-segment intersection tests
#include "../geometry/geometry.hpp"
#include "../geometry/intersection/segment_intersection.hpp"
#include "../geometry/intersection/line_intersection.hpp"

namespace Euclid::Tests {

inline void testSegmentSegmentIntersection() {
    // 1. Intersecting segments in 2D
    {
        Segment<double, 2> s1(Point<double, 2>({0,0}), Point<double, 2>({1,1}));
        Segment<double, 2> s2(Point<double, 2>({0,1}), Point<double, 2>({1,0}));
        auto r = intersect(s1, s2);
        printTest("2D Segment-Segment intersecting\n", r.intersects && r.points.size() == 1);
        if (r.intersects) {
            std::cout << "Intersection at: ("
                      << r.points[0][0] << ", " << r.points[0][1] << ")\n";
        }
    }

    // 2. Parallel but non-overlapping segments
    {
        Segment<double, 2> s1(Point<double, 2>({0,0}), Point<double, 2>({1,0}));
        Segment<double, 2> s2(Point<double, 2>({0,1}), Point<double, 2>({1,1}));
        auto r = intersect(s1, s2);
        printTest("2D Segment-Segment parallel no intersection\n", !r.intersects);
    }

    // 3. Colinear overlapping segments
    {
        Segment<double, 2> s1(Point<double, 2>({0,0}), Point<double, 2>({2,0}));
        Segment<double, 2> s2(Point<double, 2>({1,0}), Point<double, 2>({3,0}));
        auto r = intersect(s1, s2);
        printTest("2D Segment-Segment colinear overlap\n", r.intersects && r.points.size() >= 2);
    }

    // 4. Intersecting segments in 3D
    {
        Segment<double, 3> s1(Point<double, 3>({0,0,0}), Point<double, 3>({1,1,0}));
        Segment<double, 3> s2(Point<double, 3>({0,1,0}), Point<double, 3>({1,0,0}));
        auto r = intersect(s1, s2);
        printTest("3D Segment-Segment intersecting\n", r.intersects);
        if (r.intersects) {
            std::cout << "Intersection at: ("
                      << r.points[0][0] << ", " << r.points[0][1] << ", " << r.points[0][2] << ")\n";
        }
    }

    // 5. Non-intersecting skew segments in 3D
    {
        Segment<double, 3> s1(Point<double, 3>({0,0,0}), Point<double, 3>({1,0,0}));
        Segment<double, 3> s2(Point<double, 3>({0,1,1}), Point<double, 3>({1,1,1}));
        auto r = intersect(s1, s2);
        printTest("3D Segment-Segment skew no intersection\n", !r.intersects);
    }

    // 6. Intersecting segments in 4D
    {
        Segment<double, 4> s1(Point<double, 4>({0,0,0,0}), Point<double, 4>({1,1,0,0}));
        Segment<double, 4> s2(Point<double, 4>({0,1,0,0}), Point<double, 4>({1,0,0,0}));
        auto r = intersect(s1, s2);
        printTest("4D Segment-Segment intersecting\n", r.intersects);
        if (r.intersects) {
            std::cout << "Intersection at: ("
                      << r.points[0][0] << ", " << r.points[0][1] << ", " << r.points[0][2] << ", " << r.points[0][3] << ")\n";
        }
    }

    // 7. Non-intersecting segments in 4D
    {
        Segment<double, 4> s1(Point<double, 4>({0,0,0,0}), Point<double, 4>({1,0,0,0}));
        Segment<double, 4> s2(Point<double, 4>({0,1,1,1}), Point<double, 4>({1,1,1,1}));
        auto r = intersect(s1, s2);
        printTest("4D Segment-Segment no intersection\n", !r.intersects);
    }

    // 8. Intersecting segments in 5D
    {
        Segment<double, 5> s1(Point<double, 5>({0,0,0,0,0}), Point<double, 5>({1,1,0,0,0}));
        Segment<double, 5> s2(Point<double, 5>({0,1,0,0,0}), Point<double, 5>({1,0,0,0,0}));
        auto r = intersect(s1, s2);
        printTest("5D Segment-Segment intersecting\n", r.intersects);
        if (r.intersects) {
            std::cout << "Intersection at: ("
                      << r.points[0][0] << ", " << r.points[0][1] << ", " << r.points[0][2] << ", "
                      << r.points[0][3] << ", " << r.points[0][4] << ")\n";
        }
    }

    // 9. Non-intersecting segments in 5D
    {
        Segment<double, 5> s1(Point<double, 5>({0,0,0,0,0}), Point<double, 5>({1,0,0,0,0}));
        Segment<double, 5> s2(Point<double, 5>({0,1,1,1,1}), Point<double, 5>({1,1,1,1,1}));
        auto r = intersect(s1, s2);
        printTest("5D Segment-Segment no intersection\n", !r.intersects);
    }
}

// Modular submethod for segment-curve intersection tests
inline void testSegmentCurveIntersection() {
    using namespace Euclid::Geometry;
    using std::cout;
    using std::endl;
    using std::setw;
    using std::fixed;
    using std::setprecision;
    std::cout << "---- Segment–Curve Intersection Tests ----\n\n";

    // 1. Parabola intersection/no-intersection tests
    std::cout << "[Parabola Intersection Tests]\n";
    {
        // Parabola: y = x^2, x in [-1,1]
        auto parabola = [](double t) {
            return Point<double,2>({t, t*t});
        };
        Curve<double,2> curve(parabola, -1.0, 1.0);

        // Segment crossing parabola (through y=0): from (-1,0) to (1,0)
        Segment<double,2> seg1({-1.0, 0.0}, {1.0, 0.0});
        auto r1 = intersect(seg1, curve);
        printTest("Segment crosses parabola (expect 2)", r1.intersects && r1.points.size() == 2);
        if (r1.intersects) {
            for (const auto& p : r1.points) {
                std::cout << "Intersection at: (" << p[0] << ", " << p[1] << ")\n";
            }
        }
        // Segment above parabola (no intersection): from (-1,2) to (1,2)
        Segment<double,2> seg2({-1.0, 2.0}, {1.0, 2.0});
        auto r2 = intersect(seg2, curve);
        printTest("Segment above parabola (no intersection)", !r2.intersects);
    }

    // 2. Analytical frequency resolution tests (y=0 segment, expect 2f+1 crossings)
    std::cout << "[Sine Curve Frequency Resolution]\n";
    std::cout << std::right
              << setw(10) << "freq" << "  "
              << setw(10) << "expected" << "  "
              << setw(7)  << "found" << "\n";
    std::cout << "-------------------------------\n";
    {
        double freq = 1.0;
        while (true) {
            // Sine curve: y = sin(2π f x), x in [0,1]
            auto sine = [freq](double t) {
                return Point<double,2>({t, std::sin(2.0*M_PI*freq*t)});
            };
            Curve<double,2> curve(sine, 0.0, 1.0);
            // Segment along y=0, from (0,0) to (1,0)
            Segment<double,2> seg({0.0, 0.0}, {1.0, 0.0});
            int expected = static_cast<int>(2*freq + 1);
            auto r = intersect(seg, curve);
            size_t found = r.points.size();
            std::cout << fixed << setprecision(3)
                      << setw(10) << freq << "  "
                      << setw(10) << expected << "  "
                      << setw(7) << found << "\n";
            if (found < static_cast<size_t>(expected)) {
                std::cout << "[FAIL] Resolution limit reached at freq=" << freq
                          << " (found=" << found << ", expected=" << expected << ")\n";
                break;
            }
            if (found > static_cast<size_t>(expected * 2)) {
                std::cout << "[FAIL] Excess intersections detected at freq=" << freq
                          << " (found=" << found << ", expected=" << expected << ")\n";
                break;
            }
            freq *= 2.0;
            if (freq > 128) break; // avoid endless loop
        }
    }

    // 3. Frequency tests below sine peaks (expect ≈2f crossings)
    std::cout << "[Sine Curve Below Peaks]\n";
    std::cout << std::right
              << setw(10) << "freq" << "  "
              << setw(10) << "expected" << "  "
              << setw(7)  << "found" << "\n";
    std::cout << "-------------------------------\n";
    {
        double freq = 1.0;
        while (true) {
            // Sine curve: y = sin(2π f x), x in [0,1]
            auto sine = [freq](double t) {
                return Point<double,2>({t, std::sin(2.0*M_PI*freq*t)});
            };
            Curve<double,2> curve(sine, 0.0, 1.0);
            // Segment along y=0.5, from (0,0.5) to (1,0.5), which cuts sine at about 2f times
            Segment<double,2> seg({0.0, 0.5}, {1.0, 0.5});
            int expected = static_cast<int>(2*freq); // Not always exact, but close
            auto r = intersect(seg, curve);
            size_t found = r.points.size();
            std::cout << fixed << setprecision(3)
                      << setw(10) << freq << "  "
                      << setw(10) << expected << "  "
                      << setw(7) << found << "\n";
            if (found < static_cast<size_t>(expected)) {
                std::cout << "[FAIL] Resolution limit reached at freq=" << freq
                          << " (found=" << found << ", expected≈" << expected << ")\n";
                break;
            }
            if (found > static_cast<size_t>(expected * 2 + 1)) {
                std::cout << "[FAIL] Excess intersections detected at freq=" << freq
                          << " (found=" << found << ", expected≈" << expected << ")\n";
                break;
            }
            freq *= 2.0;
            if (freq > 128) break;
        }
    }

    // 4. Sine wave extended-domain sanity test
    std::cout << "[Sine Curve Extended Domain]\n";
    {
        // Sine curve: y = sin(2π x), x in [-2,2]
        auto sine = [](double t) {
            return Point<double,2>({t, std::sin(2.0*M_PI*t)});
        };
        Curve<double,2> curve(sine, -2.0, 2.0);
        // Segment along y=0 from (-2,0) to (2,0): should cross at x = n*0.5, n=-4...4, so 9 crossings
        Segment<double,2> seg({-2.0, 0.0}, {2.0, 0.0});
        auto r = intersect(seg, curve);
        printTest("Segment–sine extended domain (expect 9)", r.intersects && r.points.size() == 9);
        if (r.intersects) {
            for (const auto& p : r.points) {
                std::cout << "Intersection at: (" << p[0] << ", " << p[1] << ")\n";
            }
        }
    }

    // 5. Circle intersection test (two hits expected)
    std::cout << "[Circle Intersection Test]\n";
    {
        // Unit circle: x = cos(t), y = sin(t), t in [0,2π]
        auto circle = [](double t) {
            return Point<double,2>({std::cos(t), std::sin(t)});
        };
        Curve<double,2> curve(circle, 0.0, 2.0*M_PI);
        // Segment along x-axis from (-2,0) to (2,0): should hit circle at x=-1 and x=1
        Segment<double,2> seg({-2.0, 0.0}, {2.0, 0.0});
        auto r = intersect(seg, curve);
        printTest("Segment–circle intersection (expect 2)", r.intersects && r.points.size() == 2);
        if (r.intersects) {
            for (const auto& p : r.points) {
                std::cout << "Intersection at: (" << p[0] << ", " << p[1] << ")\n";
            }
        }
    }
    std::cout << "---- End Segment–Curve Intersection Tests ----\n\n";
}

// Modular submethod for segment-plane intersection tests
inline void testSegmentPlaneIntersection() {
    // 1. 2D: Segment intersects plane (line)
    {
        Segment<double, 2> seg(Point<double, 2>({0, 0}), Point<double, 2>({2, 2}));
        Plane<double, 2> plane(Point<double, 2>({1, 1}), Plane<double, 2>::VectorType({1, -1}));
        std::cout << "2D Segment start: (" << seg.start[0] << ", " << seg.start[1] << ")\n";
        std::cout << "2D Segment end: (" << seg.end[0] << ", " << seg.end[1] << ")\n";
        std::cout << "2D Plane base: (" << plane.base[0] << ", " << plane.base[1] << ")\n";
        std::cout << "2D Plane normal: (" << plane.normal[0] << ", " << plane.normal[1] << ")\n";
        auto r = intersect(seg, plane);
        printTest("2D Segment-Plane intersection", r.intersects && !r.points.empty());
        if (r.intersects) {
            std::cout << "Intersection at: (" << r.points[0][0] << ", " << r.points[0][1] << ")\n";
        }
    }
    // 2. 2D: Segment does not intersect plane (line)
    {
        Segment<double, 2> seg(Point<double, 2>({0, 2}), Point<double, 2>({2, 4}));
        Plane<double, 2> plane(Point<double, 2>({1, 1}), Plane<double, 2>::VectorType({1, -1}));
        std::cout << "2D Segment start: (" << seg.start[0] << ", " << seg.start[1] << ")\n";
        std::cout << "2D Segment end: (" << seg.end[0] << ", " << seg.end[1] << ")\n";
        std::cout << "2D Plane base: (" << plane.base[0] << ", " << plane.base[1] << ")\n";
        std::cout << "2D Plane normal: (" << plane.normal[0] << ", " << plane.normal[1] << ")\n";
        auto r = intersect(seg, plane);
        printTest("2D Segment-Plane no intersection", !r.intersects);
    }
    // 3. 3D: Segment intersects plane
    {
        Segment<double, 3> seg(Point<double, 3>({0, 0, 0}), Point<double, 3>({0, 0, 2}));
        Plane<double, 3> plane(Point<double, 3>({0, 0, 1}), Plane<double, 3>::VectorType({0, 0, 1}));
        std::cout << "3D Segment start: (" << seg.start[0] << ", " << seg.start[1] << ", " << seg.start[2] << ")\n";
        std::cout << "3D Segment end: (" << seg.end[0] << ", " << seg.end[1] << ", " << seg.end[2] << ")\n";
        std::cout << "3D Plane base: (" << plane.base[0] << ", " << plane.base[1] << ", " << plane.base[2] << ")\n";
        std::cout << "3D Plane normal: (" << plane.normal[0] << ", " << plane.normal[1] << ", " << plane.normal[2] << ")\n";
        auto r = intersect(seg, plane);
        printTest("3D Segment-Plane intersection", r.intersects && !r.points.empty());
        if (r.intersects) {
            std::cout << "Intersection at: (" << r.points[0][0] << ", " << r.points[0][1] << ", " << r.points[0][2] << ")\n";
        }
    }
    // 4. 3D: Segment does not intersect plane
    {
        Segment<double, 3> seg(Point<double, 3>({0, 0, 2}), Point<double, 3>({0, 0, 4}));
        Plane<double, 3> plane(Point<double, 3>({0, 0, 1}), Plane<double, 3>::VectorType({0, 0, 1}));
        std::cout << "3D Segment start: (" << seg.start[0] << ", " << seg.start[1] << ", " << seg.start[2] << ")\n";
        std::cout << "3D Segment end: (" << seg.end[0] << ", " << seg.end[1] << ", " << seg.end[2] << ")\n";
        std::cout << "3D Plane base: (" << plane.base[0] << ", " << plane.base[1] << ", " << plane.base[2] << ")\n";
        std::cout << "3D Plane normal: (" << plane.normal[0] << ", " << plane.normal[1] << ", " << plane.normal[2] << ")\n";
        auto r = intersect(seg, plane);
        printTest("3D Segment-Plane no intersection", !r.intersects);
    }
    // 5. 4D: Segment intersects hyperplane
    {
        Segment<double, 4> seg(Point<double, 4>({0, 0, 0, 0}), Point<double, 4>({0, 0, 0, 2}));
        Plane<double, 4> plane(Point<double, 4>({0, 0, 0, 1}), Plane<double, 4>::VectorType({0, 0, 0, 1}));
        std::cout << "4D Segment start: (" << seg.start[0] << ", " << seg.start[1] << ", " << seg.start[2] << ", " << seg.start[3] << ")\n";
        std::cout << "4D Segment end: (" << seg.end[0] << ", " << seg.end[1] << ", " << seg.end[2] << ", " << seg.end[3] << ")\n";
        std::cout << "4D Plane base: (" << plane.base[0] << ", " << plane.base[1] << ", " << plane.base[2] << ", " << plane.base[3] << ")\n";
        std::cout << "4D Plane normal: (" << plane.normal[0] << ", " << plane.normal[1] << ", " << plane.normal[2] << ", " << plane.normal[3] << ")\n";
        auto r = intersect(seg, plane);
        printTest("4D Segment-Plane intersection", r.intersects && !r.points.empty());
        if (r.intersects) {
            std::cout << "Intersection at: (" << r.points[0][0] << ", " << r.points[0][1] << ", " << r.points[0][2] << ", " << r.points[0][3] << ")\n";
        }
    }
    // 6. 4D: Segment does not intersect hyperplane
    {
        Segment<double, 4> seg(Point<double, 4>({0, 0, 0, 2}), Point<double, 4>({0, 0, 0, 3}));
        Plane<double, 4> plane(Point<double, 4>({0, 0, 0, 1}), Plane<double, 4>::VectorType({0, 0, 0, 1}));
        std::cout << "4D Segment start: (" << seg.start[0] << ", " << seg.start[1] << ", " << seg.start[2] << ", " << seg.start[3] << ")\n";
        std::cout << "4D Segment end: (" << seg.end[0] << ", " << seg.end[1] << ", " << seg.end[2] << ", " << seg.end[3] << ")\n";
        std::cout << "4D Plane base: (" << plane.base[0] << ", " << plane.base[1] << ", " << plane.base[2] << ", " << plane.base[3] << ")\n";
        std::cout << "4D Plane normal: (" << plane.normal[0] << ", " << plane.normal[1] << ", " << plane.normal[2] << ", " << plane.normal[3] << ")\n";
        auto r = intersect(seg, plane);
        printTest("4D Segment-Plane no intersection", !r.intersects);
    }
}


inline void testSegmentFaceIntersection() {
    std::cout << "Testing Segment-Face intersections...\n";

    // 2D face (line segment as face)
    {
        // Face defined by two points (line segment): {0,0} to {2,0}, normal {0,1}
        std::vector<Point<double, 2>> verts2 = {{0,0}, {2,0}};
        Point<double, 2> base2 = verts2[0];
        Face<double, 2>::VectorType normal2({0,1});
        Face<double, 2> face2(base2, normal2, verts2);
        Segment<double, 2> seg2({1, -1}, {1, 1});
        auto r2 = intersect(seg2, face2);
        printTest("2D Segment intersects face (line)", r2.intersects);
        if (r2.intersects && !r2.points.empty()) {
            std::cout << "Intersection at: (" << r2.points[0][0] << ", " << r2.points[0][1] << ")\n";
        }
        // Segment outside face (parallel above)
        Segment<double, 2> seg2_out({0, 2}, {2, 2});
        auto r2_out = intersect(seg2_out, face2);
        printTest("2D Segment outside face (line)", !r2_out.intersects);
    }

    // Simple 3D triangle
    {
        std::vector<Point<double, 3>> verts3 = {{0,0,0}, {1,0,0}, {0,1,0}};
        Point<double, 3> base3 = verts3[0];
        Face<double, 3>::VectorType normal3({0,0,1});
        Face<double, 3> face3(base3, normal3, verts3);
        Segment<double, 3> seg3({0.2,0.2,-1}, {0.2,0.2,1});
        auto r3 = intersect(seg3, face3);
        printTest("3D Segment intersects face\n", r3.intersects);
        if (r3.intersects && !r3.points.empty()) {
            std::cout << "Intersection at: (" << r3.points[0][0] << ", " << r3.points[0][1] << ", " << r3.points[0][2] << ")\n";
        }
        // Segment lying completely in face plane
        Segment<double, 3> seg3_plane({0.1,0.1,0}, {0.2,0.2,0});
        auto r3_plane = intersect(seg3_plane, face3);
        printTest("3D Segment in face plane", r3_plane.intersects && r3_plane.points.size() == 2);
        if (r3_plane.intersects && r3_plane.points.size() == 2) {
            std::cout << "Intersection at: (" << r3_plane.points[0][0] << ", " << r3_plane.points[0][1] << ", " << r3_plane.points[0][2] << ")\n";
            std::cout << "Intersection at: (" << r3_plane.points[1][0] << ", " << r3_plane.points[1][1] << ", " << r3_plane.points[1][2] << ")\n";
        }
        // Segment outside face
        Segment<double, 3> seg3_out({-1,-1,-1}, {-0.5,-0.5,-0.5});
        auto r3_out = intersect(seg3_out, face3);
        printTest("3D Segment outside face\n", !r3_out.intersects);
    }

    // 4D triangle face
    {
        std::vector<Point<double, 4>> verts4 = {{0,0,0,0}, {1,0,0,0}, {0,1,0,0}};
        Point<double, 4> base4 = verts4[0];
        Face<double, 4>::VectorType normal4({0,0,1,0});
        Face<double, 4> face4(base4, normal4, verts4);
        Segment<double, 4> seg4({0.1,0.1,-1,0}, {0.1,0.1,1,0});
        auto r4 = intersect(seg4, face4);
        printTest("4D Segment intersects face", r4.intersects);
        if (r4.intersects && !r4.points.empty()) {
            std::cout << "Intersection at: (" << r4.points[0][0] << ", " << r4.points[0][1] << ", " << r4.points[0][2] << ", " << r4.points[0][3] << ")\n";
        }
        // Segment outside face
        Segment<double, 4> seg4_out({-1,-1,-1,-1}, {-0.5,-0.5,-0.5,-0.5});
        auto r4_out = intersect(seg4_out, face4);
        printTest("4D Segment outside face\n", !r4_out.intersects);
    }

    // 5D triangle face
    {
        std::vector<Point<double, 5>> verts5 = {{0,0,0,0,0}, {1,0,0,0,0}, {0,1,0,0,0}};
        Point<double, 5> base5 = verts5[0];
        Face<double, 5>::VectorType normal5({0,0,1,0,0});
        Face<double, 5> face5(base5, normal5, verts5);
        Segment<double, 5> seg5({0.1,0.1,-1,0,0}, {0.1,0.1,1,0,0});
        auto r5 = intersect(seg5, face5);
        printTest("5D Segment intersects face\n", r5.intersects);
        if (r5.intersects && !r5.points.empty()) {
            std::cout << "Intersection at: (" << r5.points[0][0] << ", " << r5.points[0][1] << ", " << r5.points[0][2] << ", " << r5.points[0][3] << ", " << r5.points[0][4] << ")\n";
        }
        // Segment outside face
        Segment<double, 5> seg5_out({-1,-1,-1,-1,-1}, {-0.5,-0.5,-0.5,-0.5,-0.5});
        auto r5_out = intersect(seg5_out, face5);
        printTest("5D Segment outside face\n", !r5_out.intersects);
    }

}

// Modular submethod for segment-surface intersection tests (3D, 4D, 5D)
inline void testSegmentSurfaceIntersection() {
    std::cout << "---- Segment–Surface Intersection Tests ----\n\n";

    // 1. Plane surface (z = 0)
    std::cout << "[Plane Surface]\n";
    {
        auto planeFn = [](double u, double v) {
            return Point<double,3>({u, v, 0.0});
        };
        Surface<double,3> surf(planeFn, {{-1.0,1.0},{-1.0,1.0}});
        // Segment through the plane
        Segment<double,3> seg(Point<double,3>({0.0, 0.0, -1.0}), Point<double,3>({0.0, 0.0, 1.0}));
        auto r = intersect(seg, surf);
        printTest("Segment intersects plane surface", r.intersects && !r.points.empty());
        if (r.intersects && !r.points.empty()) {
            std::cout << "Intersection at: (" << r.points[0][0] << ", " << r.points[0][1] << ", " << r.points[0][2] << ")\n";
        }
        // Segment parallel above the plane
        Segment<double,3> seg_above(Point<double,3>({0.0, 0.0, 2.0}), Point<double,3>({1.0, 1.0, 2.0}));
        auto r2 = intersect(seg_above, surf);
        printTest("Segment above plane surface (no intersection)", !r2.intersects);
    }

    // 2. Sphere surface (unit sphere)
    std::cout << "[Sphere Surface]\n";
    {
        auto sphereFn = [](double u, double v) {
            double theta = 2*M_PI*u, phi = M_PI*v;
            double x = std::sin(phi)*std::cos(theta);
            double y = std::sin(phi)*std::sin(theta);
            double z = std::cos(phi);
            return Point<double,3>({x, y, z});
        };
        Surface<double,3> surf(sphereFn, {{-1.0,1.0},{-1.0,1.0}});
        // Segment through the sphere center
        Segment<double,3> seg(Point<double,3>({0.0, 0.0, -2.0}), Point<double,3>({0.0, 0.0, 2.0}));
        auto r = intersect(seg, surf);
        printTest("Segment intersects sphere surface (2 points expected)", r.intersects && r.points.size() == 2);
        if (r.intersects && !r.points.empty()) {
            for (size_t i = 0; i < r.points.size(); ++i) {
                std::cout << "Intersection at: (" << r.points[i][0] << ", " << r.points[i][1] << ", " << r.points[i][2] << ")\n";
            }
        }
        // Segment outside the sphere
        Segment<double,3> seg_out(Point<double,3>({2.0, 2.0, 2.0}), Point<double,3>({3.0, 3.0, 3.0}));
        auto r2 = intersect(seg_out, surf);
        printTest("Segment outside sphere surface (no intersection)", !r2.intersects);
    }

    // 3. Saddle surface (z = x^2 - y^2)
    std::cout << "[Saddle Surface]\n";
    {
        auto saddleFn = [](double u, double v) {
            double x = u, y = v;
            double z = x*x - y*y;
            return Point<double,3>({x, y, z});
        };
        Surface<double,3> surf(saddleFn, {{-1.0,1.0},{-1.0,1.0}});
        // Segment along the x axis at y=0: matches the line test for saddle
        Segment<double,3> seg(Point<double,3>({-1.0, 0.0, -1.0}), Point<double,3>({1.0, 0.0, 1.0}));
        auto r = intersect(seg, surf);
        printTest("Segment intersects saddle surface (2 points expected)", r.intersects && r.points.size() == 2);
        if (r.intersects && !r.points.empty()) {
            for (size_t i = 0; i < r.points.size(); ++i) {
                std::cout << "Intersection at: (" << r.points[i][0] << ", " << r.points[i][1] << ", " << r.points[i][2] << ")\n";
            }
        }
        // Segment away from the saddle
        Segment<double,3> seg_off(Point<double,3>({2.0, 2.0, -2.0}), Point<double,3>({2.0, 2.0, 2.0}));
        auto r2 = intersect(seg_off, surf);
        printTest("Segment outside saddle surface (no intersection)", !r2.intersects);
    }

    std::cout << "---- End Segment–Surface Intersection Tests ----\n \n";
}


// Reuse-based frequency and scaling tests for Segment–Surface intersections
inline void testSegmentRibbonSurfaceFrequencyResolution() {
    std::cout << "---- Segment–Surface Frequency Resolution Test ----\n\n";
    using namespace std::chrono;

    // Pretty, fixed-width header identical to line–surface ribbon frequency test
    std::cout << std::right
              << std::setw(10) << "freq" << "  "
              << std::setw(10) << "expected" << "  "
              << std::setw(7)  << "found" << "  "
              << std::setw(12) << "runtime(ms)" << "\n";
    std::cout << "-----------------------------------------------\n";

    // The ribbon surface domain is u in [0,1], v in [0,1]
    double freq = 1.0;
    while (true) {
        // Ribbon surface: z = sin(2π*freq*u), v in [0,1], centered at y=0
        auto ribbonFn = [freq](double u, double v) {
            return Point<double,3>({
                u,
                v - 0.5, // centered ribbon
                std::sin(2.0 * M_PI * freq * u)
            });
        };
        Surface<double,3> surf(ribbonFn, {{0.0, 1.0}, {0.0, 1.0}});

        // Segment from (0,0,0) to (1,0,0) (matches the line ribbon test geometry)
        Segment<double,3> seg(Point<double,3>({0.0, 0.0, 0.0}),
                              Point<double,3>({1.0, 0.0, 0.0}));

        // For each frequency, expect 2*freq+1 crossings (analytic formula, matches line–surface test).
        int expected = static_cast<int>(2 * freq + 1);

        auto t0 = high_resolution_clock::now();
        auto r_seg = intersect(seg, surf);
        auto t1 = high_resolution_clock::now();

        size_t n = r_seg.points.size();
        double ms = duration<double, std::milli>(t1-t0).count();

        std::cout << std::fixed << std::setprecision(3)
                  << std::setw(10) << freq << "  "
                  << std::setw(10) << expected << "  "
                  << std::setw(7)  << n << "  "
                  << std::setw(12) << ms << "\n";

        if (n < static_cast<size_t>(expected)) {
            std::cout << "[FAIL] Resolution limit reached at freq=" << freq
                      << " (found=" << n << ", expected=" << expected << ")\n";
            break;
        }
        if (n > static_cast<size_t>(expected * 2)) {
            std::cout << "[FAIL] Excess intersections detected at freq=" << freq
                      << " (found=" << n << ", expected=" << expected << ")\n";
            break;
        }
        freq *= 2.0;
    }
    std::cout << "Test completed\n";
    std::cout << "---- End Segment–Surface Frequency Resolution Test ----\n\n";
}



inline void testSegmentTorusScaleResolution() {
    std::cout << "---- Segment–Torus Scale Resolution Test ----\n";
    using namespace std::chrono;
    using namespace Euclid::Geometry;

    const double R0 = 10.0;      // major radius
    const double r0 = 2.5;       // minor radius
    const double shrink = 0.5;   // scaling factor
    const int expected = 4;      // four crossings for a full-through segment

    // Pretty, fixed-width header identical to line-torus test (tabular spacing)
    std::cout << std::right
              << std::setw(14) << "R" << "  "
              << std::setw(14) << "r" << "  "
              << std::setw(10) << "expected" << "  "
              << std::setw(7)  << "found" << "  "
              << std::setw(12) << "runtime(ms)" << "\n";
    std::cout << "---------------------------------------------------------\n";

    int step = 0;
    while (true) {
        double R = R0 * std::pow(shrink, step);
        double r = r0 * std::pow(shrink, step);

        // Torus identical to line test
        Surface<double,3> torus(
            [R, r](double u, double v) {
                return Point<double,3>({
                    (R + r * std::cos(v)) * std::cos(u),
                    (R + r * std::cos(v)) * std::sin(u),
                    r * std::sin(v)
                });
            },
            {{0.0, 2*M_PI}, {0.0, 2*M_PI}}
        );

        // Segment through torus tube center along x-axis with ±1.2 margin
        double span = 1.2 * (R + r);
        Segment<double,3> seg(Point<double,3>({-span, 0.0, 0.0}),
                              Point<double,3>({ span, 0.0, 0.0}));

        auto t0 = high_resolution_clock::now();
        auto result = intersect(seg, torus);
        auto t1 = high_resolution_clock::now();

        int found = static_cast<int>(result.points.size());
        double ms = duration<double, std::milli>(t1 - t0).count();

        // Row output: perfectly aligned scientific and fixed columns
        std::cout << std::scientific << std::setprecision(8)
                  << std::setw(14) << std::right << R << "  "
                  << std::setw(14) << std::right << r << "  "
                  << std::setw(10) << std::right << expected << "  "
                  << std::setw(7)  << std::right << found << "  "
                  << std::fixed << std::setprecision(3)
                  << std::setw(12) << std::right << ms << "\n";

        // Stop conditions with a single clearly formatted FAIL line, identical style to line test
        if (found < expected) {
            std::cout << "[FAIL] Resolution limit reached at R="
                      << std::scientific << std::setprecision(3) << R
                      << ", r=" << r
                      << " (found=" << found << ", expected=" << expected << ")\n";
            break;
        }
        if (found > expected * 2) {
            std::cout << "[FAIL] Excess intersections detected at R="
                      << std::scientific << std::setprecision(3) << R
                      << ", r=" << r
                      << " (found=" << found << ", expected=" << expected << ")\n";
            break;
        }

        ++step;
        if (R < std::numeric_limits<double>::epsilon() || r < std::numeric_limits<double>::epsilon()) {
            std::cout << "[FAIL] Reached numeric underflow before detecting resolution limit\n";
            break;
        }
    }

    std::cout << "Test completed\n";
    std::cout << "---- End Segment–Torus Scale Resolution Test ----\n\n";
}

inline void testSegmentSphereScaleResolution() {
    std::cout << "---- Segment–Sphere Scale Resolution Test ----\n\n";
    using namespace std::chrono;

    double radius = 1.0;
    double scale = 1.0;

    // Pretty, fixed-width header, identical to torus/saddle test
    std::cout << std::right
              << std::setw(14) << "scale" << "  "
              << std::setw(10) << "expected" << "  "
              << std::setw(7)  << "found" << "  "
              << std::setw(12) << "runtime(ms)" << "\n";
    std::cout << "-----------------------------------------------\n";

    // Segment through sphere center
    Segment<double,3> seg(Point<double,3>({0, 0, -2.0}), Point<double,3>({0, 0, 2.0}));

    while (scale > 1e-5) {
        double rs = radius*scale;
        auto sphereF = [rs](double u, double v) {
            double theta = 2*M_PI*u, phi = M_PI*v;
            double x = rs*std::sin(phi)*std::cos(theta);
            double y = rs*std::sin(phi)*std::sin(theta);
            double z = rs*std::cos(phi);
            return Point<double,3>({x, y, z});
        };
        Surface<double,3> sphere(sphereF, {{-1.0,1.0},{-1.0,1.0}});

        // Ground truth from Line–Surface with identical geometry
        Line<double,3> line(seg.start, seg.end - seg.start);
        auto r_line = intersect(line, sphere);
        size_t expected = r_line.points.size();

        auto t0 = high_resolution_clock::now();
        auto r_seg = intersect(seg, sphere);
        auto t1 = high_resolution_clock::now();

        size_t n = r_seg.points.size();
        double ms = duration<double, std::milli>(t1-t0).count();

        std::cout << std::scientific << std::setprecision(8)
                  << std::setw(14) << scale << "  "
                  << std::setw(10) << expected << "  "
                  << std::setw(7)  << n << "  "
                  << std::fixed << std::setprecision(3)
                  << std::setw(12) << ms << "\n";

        if (n < expected) {
            std::cout << "[FAIL] Resolution limit reached at scale="
                      << std::scientific << std::setprecision(3) << scale
                      << " (found=" << n << ", expected=" << expected << ")\n";
            break;
        }
        if (n > expected * 2) {
            std::cout << "[FAIL] Excess intersections detected at scale="
                      << std::scientific << std::setprecision(3) << scale
                      << " (found=" << n << ", expected=" << expected << ")\n";
            break;
        }
        scale /= 2.0;
        if (scale < std::numeric_limits<double>::epsilon()) {
            std::cout << "[FAIL] Reached numeric underflow before detecting resolution limit\n";
            break;
        }
    }
    std::cout << "Test completed\n";
    std::cout << "---- End Segment–Sphere Scale Resolution Test ----\n\n";
}

inline void testSegmentSaddleScaleResolution() {
    std::cout << "---- Segment–Saddle Scale Resolution Test ----\n";
    using namespace std::chrono;
    using namespace Euclid::Geometry;

    // Pretty, fixed-width header, identical to torus test but for saddle
    std::cout << std::right
              << std::setw(14) << "scale" << "  "
              << std::setw(10) << "expected" << "  "
              << std::setw(7)  << "found" << "  "
              << std::setw(12) << "runtime(ms)" << "\n";
    std::cout << "-----------------------------------------------\n";

    double scale = 1.0;
    const int expected = 2; // Always 2 for the symmetric diagonal segment through the saddle
    while (true) {
        // Saddle surface: z = x^2 - y^2, x = scale*u, y = scale*v
        auto saddleF = [scale](double u, double v) {
            double x = scale * u;
            double y = scale * v;
            double z = scale * (u*u - v*v);
            return Point<double,3>({x, y, z});
        };
        Surface<double,3> saddle(saddleF, {{-1.0,1.0},{-1.0,1.0}});

        // Symmetric diagonal segment through saddle origin: x = z, y = 0
        Segment<double,3> seg(Point<double,3>({-1.0, 0.0, -1.0}), Point<double,3>({1.0, 0.0, 1.0}));

        auto t0 = high_resolution_clock::now();
        auto r_seg = intersect(seg, saddle);
        auto t1 = high_resolution_clock::now();

        int found = static_cast<int>(r_seg.points.size());
        double ms = duration<double, std::milli>(t1 - t0).count();

        std::cout << std::scientific << std::setprecision(8)
                  << std::setw(14) << std::right << scale << "  "
                  << std::setw(10) << std::right << expected << "  "
                  << std::setw(7)  << std::right << found << "  "
                  << std::fixed << std::setprecision(3)
                  << std::setw(12) << std::right << ms << "\n";

        if (found < expected) {
            std::cout << "[FAIL] Resolution limit reached at scale="
                      << std::scientific << std::setprecision(3) << scale
                      << " (found=" << found << ", expected=" << expected << ")\n";
            break;
        }
        if (found > expected * 2) {
            std::cout << "[FAIL] Excess intersections detected at scale="
                      << std::scientific << std::setprecision(3) << scale
                      << " (found=" << found << ", expected=" << expected << ")\n";
            break;
        }
        scale /= 2.0;
        if (scale < std::numeric_limits<double>::epsilon()) {
            std::cout << "[FAIL] Reached numeric underflow before detecting resolution limit\n";
            break;
        }
    }
    std::cout << "Test completed\n";
    std::cout << "---- End Segment–Saddle Scale Resolution Test ----\n\n";
}

inline void testSegmentIntersection() {
    testSegmentSegmentIntersection();
    testSegmentCurveIntersection();
    testSegmentPlaneIntersection();
    testSegmentFaceIntersection();
    testSegmentSurfaceIntersection();
    
    
    testSegmentRibbonSurfaceFrequencyResolution();
    testSegmentSphereScaleResolution();
    testSegmentTorusScaleResolution();
    testSegmentSaddleScaleResolution();
}

} // namespace Euclid::Tests
