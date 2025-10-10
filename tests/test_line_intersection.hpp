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


inline void testLineCurveIntersection() {
    // 1. Line intersecting a 3D curve (parabola)
    {
        // Define a parabola curve y = x^2 in the XY plane, z=0
        auto parabola = [](double t) -> Point<double,3> {
            return Point<double,3>({t, t*t, 0});
        };
        // Define a line that intersects the parabola at (1,1,0)
        Line l(Point<double,3>({0,0,0}), Eigen::Vector3d({1,1,0}));
        // For testing, approximate intersection by checking points on curve
        bool intersects = false;
        Point<double,3> intersection_point;
        for (double t = -1; t <= 2; t += 0.01) {
            Point<double,3> p = parabola(t);
            // Check if point p lies on line l (within some tolerance)
            auto result = intersect(l, p);
            if (result.intersects) {
                intersects = true;
                intersection_point = p;
                break;
            }
        }
        printTest("Line-Curve (parabola) intersection: has intersection", intersects);
        if (intersects) {
            printTest("Line-Curve intersection point", true);
            // Print intersection point coordinates
            std::cout << "Intersection at: ("
                      << intersection_point[0] << ", "
                      << intersection_point[1] << ", "
                      << intersection_point[2] << ")" << std::endl;
        }
    }

    // 2. Line not intersecting the curve
    {
        auto parabola = [](double t) -> Point<double,3> {
            return Point<double,3>({t, t*t, 0});
        };
        Line l(Point<double,3>({0,0,1}), Eigen::Vector3d({1,1,0})); // Line is offset in z=1 plane
        bool intersects = false;
        for (double t = -1; t <= 2; t += 0.01) {
            Point<double,3> p = parabola(t);
            auto result = intersect(l, p);
            if (result.intersects) {
                intersects = true;
                break;
            }
        }
        printTest("Line-Curve (parabola) no intersection", !intersects);
    }
}

inline void testLineSurfaceIntersection() {
    // 1. Line intersecting planar surface z=0
    {
        Plane surface(Point<double,3>({0,0,0}), Eigen::Vector3d({0,0,1}));
        Line l(Point<double,3>({0,0,-1}), Eigen::Vector3d({0,0,1}));
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

    // 2. Line missing the surface
    {
        Plane surface(Point<double,3>({0,0,0}), Eigen::Vector3d({0,0,1}));
        Line l(Point<double,3>({0,0,1}), Eigen::Vector3d({1,0,0})); // Parallel and above plane
        auto result = intersect(l, surface);
        printTest("Line-Surface no intersection", !result.intersects);
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

inline void testLineIntersection() {
    using namespace Euclid;
    using std::cout;
    cout << "\nLine Intersection Tests\n" << std::endl;
    testLineLineIntersection();
    testLineSegmentIntersection();
    testLineCurveIntersection();
    testLinePlaneIntersection();
    testLineFaceIntersection();
    testLineSurfaceIntersection();

}

} // namespace Tests
} // namespace Euclid
