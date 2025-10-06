#pragma once


#include "geometry/geometry.hpp"
#include "geometry/intersection/line_intersection.hpp"
#include "test_utilities.hpp"

namespace Euclid {
namespace Tests {

inline void testLineIntersection() {
    using namespace Euclid;
    using std::cout;
    cout << "Running Line Intersection Tests..." << std::endl;

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

    // 5. Line-Point intersection (point on line and point off line)
    {
        Line l(Point<double,3>({1,1,1}), Eigen::Vector3d({1,1,1}));
        Point<double,3> p_on({2,2,2});
        Point<double,3> p_off({1,2,3});
        auto result_on = intersect(l, p_on);
        printTest("Line-Point (on line): has intersection", result_on.intersects);
        if (result_on.intersects) {
            printTest("Line-Point (on line): result is Point", !result_on.points.empty());
            if (!result_on.points.empty()) {
                Point<double,3> p = result_on.points[0];
                printTest("Line-Point (on line): intersection is the point", p == p_on);
            }
        }
        auto result_off = intersect(l, p_off);
        printTest("Line-Point (off line): no intersection", !result_off.intersects);
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

} // namespace Tests
} // namespace Euclid
