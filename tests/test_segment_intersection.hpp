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
        printTest("2D Segment-Segment intersecting", r.intersects && r.points.size() == 1);
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
        printTest("2D Segment-Segment parallel no intersection", !r.intersects);
    }

    // 3. Colinear overlapping segments
    {
        Segment<double, 2> s1(Point<double, 2>({0,0}), Point<double, 2>({2,0}));
        Segment<double, 2> s2(Point<double, 2>({1,0}), Point<double, 2>({3,0}));
        auto r = intersect(s1, s2);
        printTest("2D Segment-Segment colinear overlap", r.intersects && r.points.size() >= 2);
    }

    // 4. Intersecting segments in 3D
    {
        Segment<double, 3> s1(Point<double, 3>({0,0,0}), Point<double, 3>({1,1,0}));
        Segment<double, 3> s2(Point<double, 3>({0,1,0}), Point<double, 3>({1,0,0}));
        auto r = intersect(s1, s2);
        printTest("3D Segment-Segment intersecting", r.intersects);
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
        printTest("3D Segment-Segment skew no intersection", !r.intersects);
    }

    // 6. Intersecting segments in 4D
    {
        Segment<double, 4> s1(Point<double, 4>({0,0,0,0}), Point<double, 4>({1,1,0,0}));
        Segment<double, 4> s2(Point<double, 4>({0,1,0,0}), Point<double, 4>({1,0,0,0}));
        auto r = intersect(s1, s2);
        printTest("4D Segment-Segment intersecting", r.intersects);
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
        printTest("4D Segment-Segment no intersection", !r.intersects);
    }

    // 8. Intersecting segments in 5D
    {
        Segment<double, 5> s1(Point<double, 5>({0,0,0,0,0}), Point<double, 5>({1,1,0,0,0}));
        Segment<double, 5> s2(Point<double, 5>({0,1,0,0,0}), Point<double, 5>({1,0,0,0,0}));
        auto r = intersect(s1, s2);
        printTest("5D Segment-Segment intersecting", r.intersects);
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
        printTest("5D Segment-Segment no intersection", !r.intersects);
    }
}

// Modular submethod for segment-curve intersection tests
inline void testSegmentCurveIntersection() {

    // 4c. Trivial 4D curve intersection (curve varies in all coordinates, segment intersects)
    {
        auto curve = Curve<double, 4>(
            [](double t) -> Point<double, 4> {
                // Line from (0,0,0,0) to (1,1,1,1)
                return Point<double, 4>({t, t, t, t});
            },
            0.0, 1.0
        );
        // Segment from (0,1,1,1) to (1,0,0,0), which passes through (0.5,0.5,0.5,0.5) at t=0.5
        Segment<double, 4> seg(Point<double, 4>({0, 1, 1, 1}), Point<double, 4>({1, 0, 0, 0}));
        std::cout << "Segment start: (" << seg.start[0] << ", " << seg.start[1] << ", " << seg.start[2] << ", " << seg.start[3] << ")\n";
        std::cout << "Segment end: (" << seg.end[0] << ", " << seg.end[1] << ", " << seg.end[2] << ", " << seg.end[3] << ")\n";
        auto [t_min, t_max] = curve.domain();
        for (double t = t_min; t <= t_max; t += 0.1) {
            auto pt = curve.evaluate(t);
            std::cout << "Curve point at t=" << t << ": (" << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << ")\n";
        }
        auto r = intersect(seg, curve);
        printTest("4D Trivial Segment-ParametricCurve intersection", r.intersects && !r.points.empty());
        if (r.intersects) {
            for (size_t i = 0; i < r.points.size(); ++i) {
                std::cout << "Intersection at: ("
                          << r.points[i][0] << ", " << r.points[i][1] << ", " << r.points[i][2] << ", " << r.points[i][3] << ")\n";
            }
        }
    }

    // 4d. Trivial 5D curve intersection (curve varies in all coordinates, segment intersects)
    {
        auto curve = Curve<double, 5>(
            [](double t) -> Point<double, 5> {
                // Line from (0,0,0,0,0) to (1,1,1,1,1)
                return Point<double, 5>({t, t, t, t, t});
            },
            0.0, 1.0
        );
        // Segment from (0,1,1,1,1) to (1,0,0,0,0), which passes through (0.5,0.5,0.5,0.5,0.5) at t=0.5
        Segment<double, 5> seg(Point<double, 5>({0, 1, 1, 1, 1}), Point<double, 5>({1, 0, 0, 0, 0}));
        std::cout << "Segment start: (" << seg.start[0] << ", " << seg.start[1] << ", " << seg.start[2] << ", " << seg.start[3] << ", " << seg.start[4] << ")\n";
        std::cout << "Segment end: (" << seg.end[0] << ", " << seg.end[1] << ", " << seg.end[2] << ", " << seg.end[3] << ", " << seg.end[4] << ")\n";
        auto [t_min, t_max] = curve.domain();
        for (double t = t_min; t <= t_max; t += 0.1) {
            auto pt = curve.evaluate(t);
            std::cout << "Curve point at t=" << t << ": (" << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << ", " << pt[4] << ")\n";
        }
        auto r = intersect(seg, curve);
        printTest("5D Trivial Segment-ParametricCurve intersection", r.intersects && !r.points.empty());
        if (r.intersects) {
            for (size_t i = 0; i < r.points.size(); ++i) {
                std::cout << "Intersection at: ("
                          << r.points[i][0] << ", " << r.points[i][1] << ", " << r.points[i][2] << ", "
                          << r.points[i][3] << ", " << r.points[i][4] << ")\n";
            }
        }
    }
    // 1. Segment intersects a quadratic Bezier curve in 2D
    {
        auto curve = Curve<double, 2>(
            [](double t) -> Point<double, 2> {
                double u = 1.0 - t;
                return Point<double, 2>({u*u*0 + 2*u*t*1 + t*t*2, u*u*0 + 2*u*t*1 + t*t*0});
            },
            0.0, 1.0
        );
        // Segment crossing the curve
        Segment<double, 2> seg(Point<double, 2>({1, -1}), Point<double, 2>({1, 2}));
        auto r = intersect(seg, curve);
        printTest("2D Segment-QuadraticCurve intersection", r.intersects && !r.points.empty());
        if (r.intersects) {
            for (size_t i = 0; i < r.points.size(); ++i) {
                std::cout << "Intersection at: ("
                          << r.points[i][0] << ", " << r.points[i][1] << ")\n";
            }
        }
    }

    // 2. Segment not intersecting the curve
    {
        auto curve = Curve<double, 2>(
            [](double t) -> Point<double, 2> {
                double u = 1.0 - t;
                return Point<double, 2>({u*u*0 + 2*u*t*1 + t*t*2, u*u*0 + 2*u*t*1 + t*t*0});
            },
            0.0, 1.0
        );
        Segment<double, 2> seg(Point<double, 2>({-1, 2}), Point<double, 2>({-1, 3}));
        auto r = intersect(seg, curve);
        printTest("2D Segment-QuadraticCurve no intersection", !r.intersects);
    }

    // 3. Segment tangent to the curve
    {
        auto curve = Curve<double, 2>(
            [](double t) -> Point<double, 2> {
                double u = 1.0 - t;
                return Point<double, 2>({u*u*0 + 2*u*t*1 + t*t*2, u*u*0 + 2*u*t*1 + t*t*0});
            },
            0.0, 1.0
        );
        // Tangent segment at (1,1) with direction matching curve derivative
        Segment<double, 2> seg(Point<double, 2>({1, 1}), Point<double, 2>({2, 0}));
        auto r = intersect(seg, curve);
        printTest("2D Segment-QuadraticCurve tangent", r.intersects && !r.points.empty());
        if (r.intersects) {
            for (size_t i = 0; i < r.points.size(); ++i) {
                std::cout << "Intersection at: ("
                          << r.points[i][0] << ", " << r.points[i][1] << ")\n";
            }
        }
    }

    // 4. Segment intersects a simple parametric curve in 3D
    {
        auto curve = Curve<double, 3>(
            [](double t) -> Point<double, 3> {
                double u = 1.0 - t;
                return Point<double, 3>({
                    u*u*0 + 2*u*t*1 + t*t*2,
                    u*u*0 + 2*u*t*1 + t*t*0,
                    0.0
                });
            },
            0.0, 1.0
        );
        Segment<double, 3> seg(Point<double, 3>({1, -1, 0}), Point<double, 3>({1, 2, 0}));
        // Debug prints of segment endpoints and curve points
        std::cout << "Segment start: (" << seg.start[0] << ", " << seg.start[1] << ", " << seg.start[2] << ")\n";
        std::cout << "Segment end: (" << seg.end[0] << ", " << seg.end[1] << ", " << seg.end[2] << ")\n";
        auto [t_min, t_max] = curve.domain();
        for (double t = t_min; t <= t_max; t += 0.1) {
            auto pt = curve.evaluate(t);
            std::cout << "Curve point at t=" << t << ": (" << pt[0] << ", " << pt[1] << ", " << pt[2] << ")\n";
        }
        auto r = intersect(seg, curve);
        printTest("3D Segment-QuadraticCurve intersection", r.intersects && !r.points.empty());
        if (r.intersects) {
            for (size_t i = 0; i < r.points.size(); ++i) {
                std::cout << "Intersection at: ("
                          << r.points[i][0] << ", " << r.points[i][1] << ", " << r.points[i][2] << ")\n";
                // Debug print of curve.evaluate(t) at intersection parameter could be added here if parameter known
            }
        }
    }

    // 4b. Trivial 3D curve intersection (curve varies in all coordinates, segment intersects)
    {
        auto curve = Curve<double, 3>(
            [](double t) -> Point<double, 3> {
                // Line from (0,0,0) to (1,1,1)
                return Point<double, 3>({t, t, t});
            },
            0.0, 1.0
        );
        // Segment from (0,1,1) to (1,0,0), which passes through (0.5,0.5,0.5) at t=0.5
        Segment<double, 3> seg(Point<double, 3>({0, 1, 1}), Point<double, 3>({1, 0, 0}));
        std::cout << "Segment start: (" << seg.start[0] << ", " << seg.start[1] << ", " << seg.start[2] << ")\n";
        std::cout << "Segment end: (" << seg.end[0] << ", " << seg.end[1] << ", " << seg.end[2] << ")\n";
        auto [t_min, t_max] = curve.domain();
        for (double t = t_min; t <= t_max; t += 0.1) {
            auto pt = curve.evaluate(t);
            std::cout << "Curve point at t=" << t << ": (" << pt[0] << ", " << pt[1] << ", " << pt[2] << ")\n";
        }
        auto r = intersect(seg, curve);
        printTest("3D Trivial Segment-ParametricCurve intersection", r.intersects && !r.points.empty());
        if (r.intersects) {
            for (size_t i = 0; i < r.points.size(); ++i) {
                std::cout << "Intersection at: ("
                          << r.points[i][0] << ", " << r.points[i][1] << ", " << r.points[i][2] << ")\n";
            }
        }
    }

    // 5. Segment intersects a simple parametric curve in 4D
    {
        auto curve = Curve<double, 4>(
            [](double t) -> Point<double, 4> {
                double u = 1.0 - t;
                return Point<double, 4>({
                    u*u*0 + 2*u*t*1 + t*t*2,
                    u*u*0 + 2*u*t*1 + t*t*0,
                    0.0,
                    0.0
                });
            },
            0.0, 1.0
        );
        Segment<double, 4> seg(Point<double, 4>({1, -1, 0, 0}), Point<double, 4>({1, 2, 0, 0}));
        // Debug prints of segment endpoints and curve points
        std::cout << "Segment start: (" << seg.start[0] << ", " << seg.start[1] << ", " << seg.start[2] << ", " << seg.start[3] << ")\n";
        std::cout << "Segment end: (" << seg.end[0] << ", " << seg.end[1] << ", " << seg.end[2] << ", " << seg.end[3] << ")\n";
        auto [t_min, t_max] = curve.domain();
        for (double t = t_min; t <= t_max; t += 0.1) {
            auto pt = curve.evaluate(t);
            std::cout << "Curve point at t=" << t << ": (" << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << ")\n";
        }
        auto r = intersect(seg, curve);
        printTest("4D Segment-QuadraticCurve intersection", r.intersects && !r.points.empty());
        if (r.intersects) {
            for (size_t i = 0; i < r.points.size(); ++i) {
                std::cout << "Intersection at: ("
                          << r.points[i][0] << ", " << r.points[i][1] << ", " << r.points[i][2] << ", " << r.points[i][3] << ")\n";
                // Debug print of curve.evaluate(t) at intersection parameter could be added here if parameter known
            }
        }
    }

    // 6. Segment intersects a simple parametric curve in 5D
    {
        auto curve = Curve<double, 5>(
            [](double t) -> Point<double, 5> {
                double u = 1.0 - t;
                return Point<double, 5>({
                    u*u*0 + 2*u*t*1 + t*t*2,
                    u*u*0 + 2*u*t*1 + t*t*0,
                    0.0,
                    0.0,
                    0.0
                });
            },
            0.0, 1.0
        );
        Segment<double, 5> seg(Point<double, 5>({1, -1, 0, 0, 0}), Point<double, 5>({1, 2, 0, 0, 0}));
        // Debug prints of segment endpoints and curve points
        std::cout << "Segment start: (" << seg.start[0] << ", " << seg.start[1] << ", " << seg.start[2] << ", " << seg.start[3] << ", " << seg.start[4] << ")\n";
        std::cout << "Segment end: (" << seg.end[0] << ", " << seg.end[1] << ", " << seg.end[2] << ", " << seg.end[3] << ", " << seg.end[4] << ")\n";
        auto [t_min, t_max] = curve.domain();
        for (double t = t_min; t <= t_max; t += 0.1) {
            auto pt = curve.evaluate(t);
            std::cout << "Curve point at t=" << t << ": (" << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << ", " << pt[4] << ")\n";
        }
        auto r = intersect(seg, curve);
        printTest("5D Segment-QuadraticCurve intersection", r.intersects && !r.points.empty());
        if (r.intersects) {
            for (size_t i = 0; i < r.points.size(); ++i) {
                std::cout << "Intersection at: ("
                          << r.points[i][0] << ", " << r.points[i][1] << ", " << r.points[i][2] << ", "
                          << r.points[i][3] << ", " << r.points[i][4] << ")\n";
                // Debug print of curve.evaluate(t) at intersection parameter could be added here if parameter known
            }
        }
    }
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
        printTest("3D Segment intersects face", r3.intersects);
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
        printTest("3D Segment outside face", !r3_out.intersects);
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
        printTest("4D Segment outside face", !r4_out.intersects);
    }

    // 5D triangle face
    {
        std::vector<Point<double, 5>> verts5 = {{0,0,0,0,0}, {1,0,0,0,0}, {0,1,0,0,0}};
        Point<double, 5> base5 = verts5[0];
        Face<double, 5>::VectorType normal5({0,0,1,0,0});
        Face<double, 5> face5(base5, normal5, verts5);
        Segment<double, 5> seg5({0.1,0.1,-1,0,0}, {0.1,0.1,1,0,0});
        auto r5 = intersect(seg5, face5);
        printTest("5D Segment intersects face", r5.intersects);
        if (r5.intersects && !r5.points.empty()) {
            std::cout << "Intersection at: (" << r5.points[0][0] << ", " << r5.points[0][1] << ", " << r5.points[0][2] << ", " << r5.points[0][3] << ", " << r5.points[0][4] << ")\n";
        }
        // Segment outside face
        Segment<double, 5> seg5_out({-1,-1,-1,-1,-1}, {-0.5,-0.5,-0.5,-0.5,-0.5});
        auto r5_out = intersect(seg5_out, face5);
        printTest("5D Segment outside face", !r5_out.intersects);
    }

}

// Modular submethod for segment-surface intersection tests (3D, 4D, 5D)
inline void testSegmentSurfaceIntersection() {
    std::cout << "Testing Segment-Surface intersections...\n";

    // 3D ND-compatible parametric surface: the unit square in XY at z=0
    {
        auto surf3 = Surface<double, 3>(
            [](double u, double v) -> Point<double, 3> {
                // (u, v, 0), u,v in [0,1]
                return Point<double, 3>({u, v, 0});
            },
            std::pair<std::pair<double,double>, std::pair<double,double>>{{0.0, 0.0}, {1.0, 1.0}}
        );
        // Generate mesh and export to OBJ for inspection
        auto mesh3 = Euclid::Geometry::generateSurfaceMesh(surf3, 20, 20);
        if (auto objData = mesh3.exportOBJ()) {
            std::ofstream("./surf3.obj") << *objData;
        }
        std::cout << "3D Surface mesh written to ./surf3.obj\n";
        // Segment piercing the center of the surface at (0.5,0.5,0)
        Segment<double, 3> seg3({0.5,0.5,-1}, {0.5,0.5,1});
        std::cout << "3D Segment start: (" << seg3.start[0] << ", " << seg3.start[1] << ", " << seg3.start[2] << ")\n";
        std::cout << "3D Segment end:   (" << seg3.end[0] << ", " << seg3.end[1] << ", " << seg3.end[2] << ")\n";
        std::cout << "3D Surface evaluation grid (u,v):\n";
        for (double u = 0.0; u <= 1.001; u += 0.25) {
            for (double v = 0.0; v <= 1.001; v += 0.25) {
                auto pt = surf3.evaluate(u, v);
                std::cout << "  surf(" << u << "," << v << ") = ("
                          << pt[0] << ", " << pt[1] << ", " << pt[2] << ")\n";
            }
        }
        // Line-surface intersection test
        auto line_test = Line<double, 3>(seg3.start, seg3.end);
        auto lr = intersect(line_test, surf3);
        printTest("3D Line intersects parametric surface", lr.intersects);
        if (lr.intersects) {
            for (const auto& pt : lr.points) {
                std::cout << "Line intersection at: (" << pt[0] << ", " << pt[1] << ", " << pt[2] << ")\n";
            }
        }
        auto r3 = intersect(seg3, surf3);
        printTest("3D Segment intersects parametric surface", r3.intersects);
        if (r3.intersects) {
            for (size_t i = 0; i < r3.points.size(); ++i) {
                std::cout << "Intersection at: (" << r3.points[i][0] << ", " << r3.points[i][1] << ", " << r3.points[i][2] << ")\n";
            }
        }
        // Segment outside surface
        Segment<double, 3> seg3_out({2,2,-1}, {2,2,1});
        std::cout << "3D OUT Segment start: (" << seg3_out.start[0] << ", " << seg3_out.start[1] << ", " << seg3_out.start[2] << ")\n";
        std::cout << "3D OUT Segment end:   (" << seg3_out.end[0] << ", " << seg3_out.end[1] << ", " << seg3_out.end[2] << ")\n";
        auto r3_out = intersect(seg3_out, surf3);
        printTest("3D Segment outside parametric surface", !r3_out.intersects);
    }

    // 4D ND-compatible parametric surface: unit square in (x,y) at z=0, w=0
    {
        auto surf4 = Surface<double, 4>(
            [](double u, double v) -> Point<double, 4> {
                // (u, v, 0, 0), u,v in [0,1]
                return Point<double, 4>({u, v, 0, 0});
            },
            std::pair<std::pair<double,double>, std::pair<double,double>>{{0.0, 0.0}, {1.0, 1.0}}
        );
        // Segment piercing the center of the surface at (0.5,0.5,0,0)
        Segment<double, 4> seg4({0.5,0.5,-1,0}, {0.5,0.5,1,0});
        std::cout << "4D Segment start: (" << seg4.start[0] << ", " << seg4.start[1] << ", " << seg4.start[2] << ", " << seg4.start[3] << ")\n";
        std::cout << "4D Segment end:   (" << seg4.end[0] << ", " << seg4.end[1] << ", " << seg4.end[2] << ", " << seg4.end[3] << ")\n";
        std::cout << "4D Surface evaluation grid (u,v):\n";
        for (double u = 0.0; u <= 1.001; u += 0.25) {
            for (double v = 0.0; v <= 1.001; v += 0.25) {
                auto pt = surf4.evaluate(u, v);
                std::cout << "  surf(" << u << "," << v << ") = ("
                          << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << ")\n";
            }
        }
        // Line-surface intersection test
        auto line_test = Line<double, 4>(seg4.start, seg4.end);
        auto lr = intersect(line_test, surf4);
        printTest("4D Line intersects parametric surface", lr.intersects);
        if (lr.intersects) {
            for (const auto& pt : lr.points) {
                std::cout << "Line intersection at: (" << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << ")\n";
            }
        }
        auto r4 = intersect(seg4, surf4);
        printTest("4D Segment intersects parametric surface", r4.intersects);
        if (r4.intersects) {
            for (size_t i = 0; i < r4.points.size(); ++i) {
                std::cout << "Intersection at: (" << r4.points[i][0] << ", " << r4.points[i][1] << ", " << r4.points[i][2] << ", " << r4.points[i][3] << ")\n";
            }
        }
        // Segment outside surface
        Segment<double, 4> seg4_out({2,2,-1,0}, {2,2,1,0});
        std::cout << "4D OUT Segment start: (" << seg4_out.start[0] << ", " << seg4_out.start[1] << ", " << seg4_out.start[2] << ", " << seg4_out.start[3] << ")\n";
        std::cout << "4D OUT Segment end:   (" << seg4_out.end[0] << ", " << seg4_out.end[1] << ", " << seg4_out.end[2] << ", " << seg4_out.end[3] << ")\n";
        auto r4_out = intersect(seg4_out, surf4);
        printTest("4D Segment outside parametric surface", !r4_out.intersects);
    }

    // 5D ND-compatible parametric surface: unit square in (x,y) at z=0, w=0, v=0
    {
        auto surf5 = Surface<double, 5>(
            [](double u, double v) -> Point<double, 5> {
                // (u, v, 0, 0, 0), u,v in [0,1]
                return Point<double, 5>({u, v, 0, 0, 0});
            },
            std::pair<std::pair<double,double>, std::pair<double,double>>{{0.0, 0.0}, {1.0, 1.0}}
        );
        // Segment piercing the center of the surface at (0.5,0.5,0,0,0)
        Segment<double, 5> seg5({0.5,0.5,-1,0,0}, {0.5,0.5,1,0,0});
        std::cout << "5D Segment start: (" << seg5.start[0] << ", " << seg5.start[1] << ", " << seg5.start[2] << ", " << seg5.start[3] << ", " << seg5.start[4] << ")\n";
        std::cout << "5D Segment end:   (" << seg5.end[0] << ", " << seg5.end[1] << ", " << seg5.end[2] << ", " << seg5.end[3] << ", " << seg5.end[4] << ")\n";
        std::cout << "5D Surface evaluation grid (u,v):\n";
        for (double u = 0.0; u <= 1.001; u += 0.25) {
            for (double v = 0.0; v <= 1.001; v += 0.25) {
                auto pt = surf5.evaluate(u, v);
                std::cout << "  surf(" << u << "," << v << ") = ("
                          << pt[0] << ", " << pt[1] << ", " << pt[2] << ", " << pt[3] << ", " << pt[4] << ")\n";
            }
        }
        // Line-surface intersection test
        auto line_test = Line<double, 5>(seg5.start, seg5.end);
        auto lr = intersect(line_test, surf5);
        printTest("5D Line intersects parametric surface", lr.intersects);
        if (lr.intersects) {
            for (const auto& pt : lr.points) {
                std::cout << "Line intersection at: (" << pt[0] << ", " << pt[1] << ", " << pt[2]
                          << ", " << pt[3] << ", " << pt[4] << ")\n";
            }
        }
        auto r5 = intersect(seg5, surf5);
        printTest("5D Segment intersects parametric surface", r5.intersects);
        if (r5.intersects) {
            for (size_t i = 0; i < r5.points.size(); ++i) {
                std::cout << "Intersection at: (" << r5.points[i][0] << ", " << r5.points[i][1] << ", " << r5.points[i][2]
                          << ", " << r5.points[i][3] << ", " << r5.points[i][4] << ")\n";
            }
        }
        // Segment outside surface
        Segment<double, 5> seg5_out({2,2,-1,0,0}, {2,2,1,0,0});
        std::cout << "5D OUT Segment start: (" << seg5_out.start[0] << ", " << seg5_out.start[1] << ", " << seg5_out.start[2] << ", " << seg5_out.start[3] << ", " << seg5_out.start[4] << ")\n";
        std::cout << "5D OUT Segment end:   (" << seg5_out.end[0] << ", " << seg5_out.end[1] << ", " << seg5_out.end[2] << ", " << seg5_out.end[3] << ", " << seg5_out.end[4] << ")\n";
        auto r5_out = intersect(seg5_out, surf5);
        printTest("5D Segment outside parametric surface", !r5_out.intersects);
    }
}


inline void testSegmentIntersection() {
    testSegmentSegmentIntersection();
    testSegmentCurveIntersection();
    testSegmentPlaneIntersection();
    testSegmentFaceIntersection();
    testSegmentSurfaceIntersection();
}

} // namespace Euclid::Tests
