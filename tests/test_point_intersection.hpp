#pragma once

#include "../geometry/intersection/intersection.hpp"
#include "test_utilities.hpp"

namespace Euclid::Tests {

void testPointIntersection() {
    using namespace Euclid::Geometry;
    
    std::cout << "\nTesting Point Intersections\n\n";

    // -----------------------
    // Point–Point Tests
    // -----------------------
    {
        Point<double, 3> p1({0.0, 0.0, 0.0});
        Point<double, 3> p2({0.0, 0.0, 0.0});
        Point<double, 3> p3({1.0, 0.0, 0.0});

        auto r1 = intersect(p1, p2);
        printTest("Point–Point coincident", r1.intersects && r1.points.size() == 1 && r1.points[0].isEqual(p1));
        std::cout << "Result: intersects = " << r1.intersects << ", points = ";
        for (const auto& pt : r1.points) {
            std::cout << "(";
            for (int i = 0; i < 3; ++i) {
                std::cout << pt[i];
                if (i < 2) std::cout << ", ";
            }
            std::cout << ") ";
        }
        std::cout << std::endl;

        auto r2 = intersect(p1, p3);
        printTest("Point–Point distinct", !r2.intersects && r2.points.empty());
        std::cout << "Result: intersects = " << r2.intersects << ", points = ";
        for (const auto& pt : r2.points) {
            std::cout << "(";
            for (int i = 0; i < 3; ++i) {
                std::cout << pt[i];
                if (i < 2) std::cout << ", ";
            }
            std::cout << ") ";
        }
        std::cout << std::endl;
    }

    // -----------------------
    // Point–Line Tests
    // -----------------------
    {
        Point<double, 3> p({1.0, 1.0, 0.0});
        Line<double, 3> l(Point<double, 3>({0.0, 0.0, 0.0}),
                           Point<double, 3>({2.0, 2.0, 0.0}));

        auto r = intersect(p, l);
        printTest("Point–Line closest point",
                  r.points.size() == 1 &&
                  r.points[0].isEqual(Point<double, 3>({1.0, 1.0, 0.0})));
        std::cout << "Result: intersects = " << r.intersects << ", points = ";
        for (const auto& pt : r.points) {
            std::cout << "(";
            for (int i = 0; i < 3; ++i) {
                std::cout << pt[i];
                if (i < 2) std::cout << ", ";
            }
            std::cout << ") ";
        }
        std::cout << std::endl;
    }

    // -----------------------
    // Point–Plane Tests
    // -----------------------
    {
        Point<double, 3> p({1.0, 1.0, 5.0});
        Plane<double, 3> plane(Point<double, 3>({0.0, 0.0, 0.0}),
                                Eigen::Matrix<double, 3, 1>({0.0, 0.0, 1.0})); // z = 0 plane

        auto r = intersect(p, plane);
        printTest("Point–Plane projection",
                  r.points.size() == 1 &&
                  r.points[0].isEqual(Point<double, 3>({1.0, 1.0, 0.0})));
        std::cout << "Result: intersects = " << r.intersects << ", points = ";
        for (const auto& pt : r.points) {
            std::cout << "(";
            for (int i = 0; i < 3; ++i) {
                std::cout << pt[i];
                if (i < 2) std::cout << ", ";
            }
            std::cout << ") ";
        }
        std::cout << std::endl;
    }

    // -----------------------
    // Point–Plane True Intersection Test
    // -----------------------
    {
        Point<double, 3> p({2.0, 3.0, 0.0});
        Plane<double, 3> plane(Point<double, 3>({0.0, 0.0, 0.0}),
                                Eigen::Matrix<double, 3, 1>({0.0, 0.0, 1.0})); // z = 0 plane

        auto r = intersect(p, plane);
        printTest("Point–Plane true intersection",
                  r.intersects &&
                  r.points.size() == 1 &&
                  r.points[0].isEqual(p));
        std::cout << "Result: intersects = " << r.intersects << ", points = ";
        for (const auto& pt : r.points) {
            std::cout << "(";
            for (int i = 0; i < 3; ++i) {
                std::cout << pt[i];
                if (i < 2) std::cout << ", ";
            }
            std::cout << ") ";
        }
        std::cout << std::endl;
    }

    // -----------------------
    // Point–Line Non-Intersecting Tests
    // -----------------------
    {
        Point<double, 3> p({3.0, 3.0, 3.0});
        Line<double, 3> l(Point<double, 3>({0.0, 0.0, 0.0}),
                          Point<double, 3>({1.0, 0.0, 0.0}));

        auto r = intersect(p, l);
        printTest("Point–Line non-intersecting",
                  !r.intersects);
        std::cout << "Result: intersects = " << r.intersects << ", points = ";
        for (const auto& pt : r.points) {
            std::cout << "(";
            for (int i = 0; i < 3; ++i) {
                std::cout << pt[i];
                if (i < 2) std::cout << ", ";
            }
            std::cout << ") ";
        }
        std::cout << std::endl;
    }

    // -----------------------
    // Point–Plane Non-Intersecting Tests
    // -----------------------
    {
        Point<double, 3> p({0.0, 0.0, 5.0});
        Plane<double, 3> plane(Point<double, 3>({0.0, 0.0, 0.0}),
                                Eigen::Matrix<double, 3, 1>({0.0, 0.0, 1.0})); // z = 0 plane

        auto r = intersect(p, plane);
        printTest("Point–Plane non-intersecting",
                  !r.intersects);
        std::cout << "Result: intersects = " << r.intersects << ", points = ";
        for (const auto& pt : r.points) {
            std::cout << "(";
            for (int i = 0; i < 3; ++i) {
                std::cout << pt[i];
                if (i < 2) std::cout << ", ";
            }
            std::cout << ") ";
        }
        std::cout << std::endl;
    }

    // -----------------------
    // Point–Curve Tests (placeholder)
    // -----------------------
    {
        // Placeholder test: no intersection expected yet
        printTest("Point–Curve placeholder test (expected fail)", false);
        std::cout << "Result: intersects = false, points = " << std::endl;
    }

    // -----------------------
    // Point–Surface Tests (placeholder)
    // -----------------------
    {
        // Placeholder test: no intersection expected yet
        printTest("Point–Surface placeholder test (expected fail)", false);
        std::cout << "Result: intersects = false, points = " << std::endl;
    }
}

} // namespace Euclid::Tests
