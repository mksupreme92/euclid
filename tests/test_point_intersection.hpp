#pragma once

#include "../geometry/intersection/intersection.hpp"
#include "test_utilities.hpp"

namespace Euclid::Tests {

// Sub-functions for each test category
namespace {
    using namespace Euclid::Geometry;

    void testPointPoint() {
        // -----------------------
        // Point–Point Tests
        // -----------------------
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

    void testPointLine() {
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
    }

    void testPointSegment() {
        // -----------------------
        // Point–Segment Tests
        // -----------------------
        {
            // Point exactly on the segment
            Point<double, 3> p({1.0, 1.0, 0.0});
            Segment<double, 3> seg(Point<double, 3>({0.0, 0.0, 0.0}),
                                   Point<double, 3>({2.0, 2.0, 0.0}));

            auto r1 = intersect(p, seg);
            printTest("Point–Segment on segment (p, seg)",
                      r1.intersects && r1.points.size() == 1 && r1.points[0].isEqual(p));
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

            auto r1b = intersect(seg, p);
            printTest("Point–Segment on segment (seg, p)",
                      r1b.intersects && r1b.points.size() == 1 && r1b.points[0].isEqual(p));
            std::cout << "Result: intersects = " << r1b.intersects << ", points = ";
            for (const auto& pt : r1b.points) {
                std::cout << "(";
                for (int i = 0; i < 3; ++i) {
                    std::cout << pt[i];
                    if (i < 2) std::cout << ", ";
                }
                std::cout << ") ";
            }
            std::cout << std::endl;
        }
        {
            // Point outside segment bounds (closest point should be endpoint)
            Point<double, 3> p({4.0, 4.0, 0.0});
            Segment<double, 3> seg(Point<double, 3>({0.0, 0.0, 0.0}),
                                   Point<double, 3>({2.0, 2.0, 0.0}));

            auto r2 = intersect(p, seg);
            printTest("Point–Segment outside bounds (p, seg)",
                      !r2.intersects && r2.points.size() == 1 &&
                      r2.points[0].isEqual(Point<double, 3>({2.0, 2.0, 0.0})));
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

            auto r2b = intersect(seg, p);
            printTest("Point–Segment outside bounds (seg, p)",
                      !r2b.intersects && r2b.points.size() == 1 &&
                      r2b.points[0].isEqual(Point<double, 3>({2.0, 2.0, 0.0})));
            std::cout << "Result: intersects = " << r2b.intersects << ", points = ";
            for (const auto& pt : r2b.points) {
                std::cout << "(";
                for (int i = 0; i < 3; ++i) {
                    std::cout << pt[i];
                    if (i < 2) std::cout << ", ";
                }
                std::cout << ") ";
            }
            std::cout << std::endl;
        }
        {
            // Degenerate segment (start = end)
            Point<double, 3> p({1.0, 1.0, 1.0});
            Segment<double, 3> seg(Point<double, 3>({1.0, 1.0, 1.0}),
                                   Point<double, 3>({1.0, 1.0, 1.0}));

            auto r3 = intersect(p, seg);
            printTest("Point–Degenerate Segment coincident (p, seg)",
                      r3.intersects && r3.points.size() == 1 && r3.points[0].isEqual(p));
            std::cout << "Result: intersects = " << r3.intersects << ", points = ";
            for (const auto& pt : r3.points) {
                std::cout << "(";
                for (int i = 0; i < 3; ++i) {
                    std::cout << pt[i];
                    if (i < 2) std::cout << ", ";
                }
                std::cout << ") ";
            }
            std::cout << std::endl;

            Point<double, 3> p2({2.0, 2.0, 2.0});
            auto r3b = intersect(seg, p2);
            printTest("Point–Degenerate Segment distinct (seg, p2)",
                      !r3b.intersects && r3b.points.size() == 1 &&
                      r3b.points[0].isEqual(Point<double, 3>({1.0, 1.0, 1.0})));
            std::cout << "Result: intersects = " << r3b.intersects << ", points = ";
            for (const auto& pt : r3b.points) {
                std::cout << "(";
                for (int i = 0; i < 3; ++i) {
                    std::cout << pt[i];
                    if (i < 2) std::cout << ", ";
                }
                std::cout << ") ";
            }
            std::cout << std::endl;
        }
    }

    void testPointCurve() {
        // -----------------------
        // Point–Curve Tests
        // -----------------------
        {
            using Curve3d = Curve<double, 3>;
            // Create a linear curve from (0,0,0) to (2,2,0)
            std::vector<Point<double, 3>> ctrlpts = {
                Point<double, 3>({0.0, 0.0, 0.0}),
                Point<double, 3>({2.0, 2.0, 0.0})
            };
            auto lineCurve = Curve3d::linearCurve(ctrlpts[0], ctrlpts[1]);

            // 1. Point exactly on the curve
            Point<double, 3> p_on({1.0, 1.0, 0.0});
            auto r_on = intersect(p_on, lineCurve);
            printTest("Point–Curve on curve",
                      r_on.intersects && r_on.points.size() == 1 && r_on.points[0].isEqual(p_on));
            std::cout << "Result: intersects = " << r_on.intersects << ", closest point = ";
            for (const auto& pt : r_on.points) {
                std::cout << "(";
                for (int i = 0; i < 3; ++i) {
                    std::cout << pt[i];
                    if (i < 2) std::cout << ", ";
                }
                std::cout << ") ";
            }
            std::cout << std::endl;

            // 2. Point off the curve but near
            Point<double, 3> p_near({1.0, 1.0, 1.0});
            auto r_near = intersect(p_near, lineCurve);
            printTest("Point–Curve near but not on",
                      !r_near.intersects && r_near.points.size() == 1 &&
                      r_near.points[0].isEqual(Point<double, 3>({1.0, 1.0, 0.0})));
            std::cout << "Result: intersects = " << r_near.intersects << ", closest point = ";
            for (const auto& pt : r_near.points) {
                std::cout << "(";
                for (int i = 0; i < 3; ++i) {
                    std::cout << pt[i];
                    if (i < 2) std::cout << ", ";
                }
                std::cout << ") ";
            }
            std::cout << std::endl;

            // 3. Point far from the curve
            Point<double, 3> p_far({10.0, -5.0, 3.0});
            auto r_far = intersect(p_far, lineCurve);
            printTest("Point–Curve far from curve",
                      !r_far.intersects && r_far.points.size() == 1 &&
                      r_far.points[0].isEqual(Point<double, 3>({2.0, 2.0, 0.0})));
            std::cout << "Result: intersects = " << r_far.intersects << ", closest point = ";
            for (const auto& pt : r_far.points) {
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
        // Point–Complex Curve Tests
        // -----------------------
        {
            using Curve3d = Curve<double, 3>;
            auto parabolaFunc = [](double t) -> Point<double, 3> {
                return Point<double, 3>({t, t*t, 0.0});
            };
            Curve3d parabola(parabolaFunc, 0.0, 2.0);

            // Point exactly on the curve
            Point<double, 3> p_on({1.0, 1.0, 0.0});
            auto r_on = intersect(p_on, parabola);
            printTest("Point–Parabola on curve",
                      r_on.intersects &&
                      r_on.points.size() == 1 &&
                      r_on.points[0].isEqual(p_on));
            std::cout << "Result: intersects = " << r_on.intersects << ", closest point = ";
            for (const auto& pt : r_on.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;

            // Point off the curve
            Point<double, 3> p_off({1.0, 1.5, 0.0});
            auto r_off = intersect(p_off, parabola);
            printTest("Point–Parabola near",
                      !r_off.intersects &&
                      r_off.points.size() == 1);
            std::cout << "Result: intersects = " << r_off.intersects << ", closest point = ";
            for (const auto& pt : r_off.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;
        }
    }

    void testPointPlane() {
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
    }

    void testPointFace() {
        // -----------------------
        // Point–Face Tests
        // -----------------------
        // Define a face (triangle in the z=0 plane)
        std::vector<Point<double, 3>> faceVerts = {
            Point<double, 3>({0.0, 0.0, 0.0}),
            Point<double, 3>({2.0, 0.0, 0.0}),
            Point<double, 3>({1.0, 2.0, 0.0})
        };
        Face<double, 3> face(
            faceVerts[0],                            // base point
            Eigen::Matrix<double, 3, 1>({0, 0, 1}), // normal
            faceVerts                                // vertices
        );

        // 1. Point inside the face
        Point<double, 3> p_inside({1.0, 1.0, 0.0});
        auto r_inside = intersect(p_inside, face);
        printTest("Point–Face inside",
                  r_inside.intersects &&
                  r_inside.points.size() == 1 &&
                  r_inside.points[0].isEqual(p_inside));
        
        /*
        std::cout << "Result: intersects = " << r_inside.intersects << ", points = ";
        for (const auto& pt : r_inside.points) {
            std::cout << "(";
            for (int i = 0; i < 3; ++i) {
                std::cout << pt[i];
                if (i < 2) std::cout << ", ";
            }
            std::cout << ") ";
        }
        std::cout << std::endl;
        */

        // 2. Point outside the face
        Point<double, 3> p_out({2.0, 2.0, 0.0});
        auto r_out = intersect(p_out, face);
        printTest("Point–Face outside",
                  !r_out.intersects);
        /*
        std::cout << "Result: intersects = " << r_out.intersects << ", points = ";
        for (const auto& pt : r_out.points) {
            std::cout << "(";
            for (int i = 0; i < 3; ++i) {
                std::cout << pt[i];
                if (i < 2) std::cout << ", ";
            }
            std::cout << ") ";
        
        }
        std::cout << std::endl;
        std::cout << "[DEBUG] Point–Face outside: projected = ("
                  << r_out.points[0][0] << ", " << r_out.points[0][1] << ", " << r_out.points[0][2]
                  << "), intersects = " << r_out.intersects << std::endl;
         */

        // 3. Point above the face (should project to inside, but not intersect)
        Point<double, 3> p_above({1.0, 1.0, 3.0});
        auto r_above = intersect(p_above, face);
        printTest("Point–Face above (z offset)",
                  !r_above.intersects);
        
        /*
        std::cout << "Result: intersects = " << r_above.intersects << ", points = ";
        for (const auto& pt : r_above.points) {
            std::cout << "(";
            for (int i = 0; i < 3; ++i) {
                std::cout << pt[i];
                if (i < 2) std::cout << ", ";
            }
            std::cout << ") ";
        }
        std::cout << std::endl;
        std::cout << "[DEBUG] Point–Face above: projected = ("
                  << r_above.points[0][0] << ", " << r_above.points[0][1] << ", " << r_above.points[0][2]
                  << "), intersects = " << r_above.intersects << std::endl;
         */
    }

    void testPointSurface() {
        // -----------------------
        // Point–Surface Tests
        // -----------------------
        {
            // Define a simple parametric surface: z = 0 plane
            using Surface3d = Surface<double, 3>;
            auto planeFunc = [](double u, double v) -> Point<double, 3> {
                return Point<double, 3>({u, v, 0.0});
            };
            auto domain = std::make_pair(std::make_pair(0.0, 2.0), std::make_pair(0.0, 2.0));
            Surface3d surf(planeFunc, domain); // u,v in [0,2]

            // 1. Point exactly on the surface
            Point<double, 3> p_on({1.0, 1.0, 0.0});
            auto r_on = intersect(p_on, surf);
            printTest("Point–Surface on surface",
                      r_on.intersects &&
                      r_on.points.size() == 1 &&
                      r_on.points[0].isEqual(p_on));
            std::cout << "Result: intersects = " << r_on.intersects << ", projected point = ";
            for (const auto& pt : r_on.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;

            // 2. Point above the surface (should not intersect; get projected point)
            Point<double, 3> p_above({1.0, 1.0, 2.0});
            auto r_above = intersect(p_above, surf);
            Point<double, 3> expected_proj({1.0, 1.0, 0.0});
            double practical_tol = 0.015; // slightly larger to account for residuals
            printTest("Point–Surface above surface",
                      !r_above.intersects &&
                      r_above.points.size() == 1 &&
                      (r_above.points[0] - expected_proj).norm() <= practical_tol);
            std::cout << "Result: intersects = " << r_above.intersects << ", projected point = ";
            for (const auto& pt : r_above.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;

            // 3. Point far from the surface (should not intersect; get closest point)
            Point<double, 3> p_far({10.0, -5.0, 3.0});
            auto r_far = intersect(p_far, surf);
            // Closest point should be at (2,0,0) (since surface is [0,2]x[0,2])
            Point<double, 3> closest({2.0, 0.0, 0.0});
            printTest("Point–Surface far from surface",
                      !r_far.intersects &&
                      r_far.points.size() == 1 &&
                      r_far.points[0].isEqual(closest));
            std::cout << "Result: intersects = " << r_far.intersects << ", closest point = ";
            for (const auto& pt : r_far.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;
        }

        // Additional Point–Surface tests: Paraboloid
        {
            using Surface3d = Surface<double, 3>;
            auto paraboloidFunc = [](double u, double v) -> Point<double, 3> {
                return Point<double, 3>({u, v, u*u + v*v});
            };
            auto paraboloidDomain = std::make_pair(std::make_pair(-1.0, 1.0), std::make_pair(-1.0, 1.0));
            Surface3d paraboloid(paraboloidFunc, paraboloidDomain);
            double u0 = paraboloidDomain.first.first;
            double u1 = paraboloidDomain.first.second;
            double v0 = paraboloidDomain.second.first;
            double v1 = paraboloidDomain.second.second;
            double adaptive_tol = std::max(Tolerance().evaluateEpsilon(3), 0.01 * std::max(u1-u0, v1-v0));

            Point<double,3> p_on_parab({0.5, 0.5, 0.5*0.5 + 0.5*0.5});
            auto r_on_parab = intersect(p_on_parab, paraboloid);
            printTest("Point–Paraboloid on surface",
                      r_on_parab.intersects &&
                      r_on_parab.points.size() == 1 &&
                      (r_on_parab.points[0] - p_on_parab).norm() <= adaptive_tol);
            std::cout << "Result: intersects = " << r_on_parab.intersects << ", projected/closest point = ";
            for (const auto& pt : r_on_parab.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;

            Point<double,3> p_above_parab({0.5, 0.5, 1.0});
            auto r_above_parab = intersect(p_above_parab, paraboloid);
            // Analytical computation for symmetric point above paraboloid z = u^2 + v^2
            // Set u = v = t, P = (t0, t0, z0)
            // Solve cubic: 4 t^3 + (-2*z0 + 1)*t - t0 = 0
            {
                double t0 = 0.5, z0 = 1.0;
                // Coefficients: 4 t^3 + (-2*z0 + 1)*t - t0 = 0
                // Use Cardano's formula for cubic: 4 t^3 + a1 t + a0 = 0
                double a = 4.0;
                double c = (-2.0 * z0 + 1.0);
                double d = -t0;
                // t^3 + (c/a) t + (d/a) = 0
                double p = c / a;
                double q = d / a;
                // Depressed cubic: t^3 + p t + q = 0
                double discrim = (q*q)/4.0 + (p*p*p)/27.0;
                double t = 0.0;
                if (discrim >= 0.0) {
                    double sqrt_disc = std::sqrt(discrim);
                    double A = std::cbrt(-q/2.0 + sqrt_disc);
                    double B = std::cbrt(-q/2.0 - sqrt_disc);
                    t = A + B;
                } else {
                    double r = std::sqrt(-(p*p*p)/27.0);
                    double phi = std::acos(-q/(2.0*r));
                    t = 2*std::cbrt(r) * std::cos(phi/3.0);
                }
                Point<double,3> expected_proj({t, t, 2*t*t});
                printTest("Point–Paraboloid above surface",
                          !r_above_parab.intersects &&
                          r_above_parab.points.size() == 1 &&
                          (r_above_parab.points[0] - expected_proj).norm() <= adaptive_tol);
            }
            std::cout << "Result: intersects = " << r_above_parab.intersects << ", projected/closest point = ";
            for (const auto& pt : r_above_parab.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;

            Point<double,3> p_far_parab({2.0, 2.0, 5.0});
            auto r_far_parab = intersect(p_far_parab, paraboloid);
            printTest("Point–Paraboloid far from surface",
                      !r_far_parab.intersects &&
                      r_far_parab.points.size() == 1);
            std::cout << "Result: intersects = " << r_far_parab.intersects << ", projected/closest point = ";
            for (const auto& pt : r_far_parab.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;
        }

        // Additional Point–Surface tests: Saddle (Hyperbolic paraboloid)
        {
            using Surface3d = Surface<double, 3>;
            auto saddleFunc = [](double u, double v) -> Point<double, 3> {
                return Point<double, 3>({u, v, u*v});
            };
            auto saddleDomain = std::make_pair(std::make_pair(-1.0,1.0), std::make_pair(-1.0,1.0));
            Surface3d saddle(saddleFunc, saddleDomain);
            double u0 = saddleDomain.first.first;
            double u1 = saddleDomain.first.second;
            double v0 = saddleDomain.second.first;
            double v1 = saddleDomain.second.second;
            double adaptive_tol = std::max(Tolerance().evaluateEpsilon(3), 0.01 * std::max(u1-u0, v1-v0));

            Point<double,3> p_on_saddle({0.5,0.5,0.25});
            auto r_on_saddle = intersect(p_on_saddle, saddle);
            printTest("Point–Saddle on surface",
                      r_on_saddle.intersects &&
                      r_on_saddle.points.size() == 1 &&
                      (r_on_saddle.points[0] - p_on_saddle).norm() <= adaptive_tol);
            std::cout << "Result: intersects = " << r_on_saddle.intersects << ", projected/closest point = ";
            for (const auto& pt : r_on_saddle.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;

            Point<double,3> p_above_saddle({0.5,0.5,0.5});
            auto r_above_saddle = intersect(p_above_saddle, saddle);
            // Tolerance-based validation: check closest point is on the surface within adaptive_tol
            bool valid_proj = false;
            if (r_above_saddle.points.size() == 1) {
                // Find (u,v) such that saddleFunc(u,v) == r_above_saddle.points[0]
                // For the test, we can invert: u = r[0], v = r[1], z = r[2], and check z == u*v
                double u_proj = r_above_saddle.points[0][0];
                double v_proj = r_above_saddle.points[0][1];
                double z_proj = r_above_saddle.points[0][2];
                double surf_val = u_proj * v_proj;
                valid_proj = std::abs(z_proj - surf_val) <= adaptive_tol;
            }
            printTest("Point–Saddle above surface",
                      !r_above_saddle.intersects &&
                      r_above_saddle.points.size() == 1 &&
                      valid_proj);
            std::cout << "Result: intersects = " << r_above_saddle.intersects << ", projected/closest point = ";
            for (const auto& pt : r_above_saddle.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;

            Point<double,3> p_far_saddle({2.0,2.0,0.0});
            auto r_far_saddle = intersect(p_far_saddle, saddle);
            printTest("Point–Saddle far from surface",
                      !r_far_saddle.intersects &&
                      r_far_saddle.points.size() == 1);
            std::cout << "Result: intersects = " << r_far_saddle.intersects << ", projected/closest point = ";
            for (const auto& pt : r_far_saddle.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;
        }

        // Additional Point–Surface tests: Sinusoidal surface
        {
            using Surface3d = Surface<double, 3>;
            auto sinFunc = [](double u, double v) -> Point<double,3> {
                return Point<double,3>({u, v, std::sin(u)*std::cos(v)});
            };
            auto sinDomain = std::make_pair(std::make_pair(0.0, 2*M_PI), std::make_pair(0.0, 2*M_PI));
            Surface3d sinSurface(sinFunc, sinDomain);
            double u0 = sinDomain.first.first;
            double u1 = sinDomain.first.second;
            double v0 = sinDomain.second.first;
            double v1 = sinDomain.second.second;
            double adaptive_tol = std::max(Tolerance().evaluateEpsilon(3), 0.01 * std::max(u1-u0, v1-v0));

            Point<double,3> p_on_sin({M_PI/2, M_PI/2, std::sin(M_PI/2)*std::cos(M_PI/2)});
            auto r_on_sin = intersect(p_on_sin, sinSurface);
            printTest("Point–SinSurface on surface",
                      r_on_sin.intersects &&
                      r_on_sin.points.size() == 1 &&
                      (r_on_sin.points[0] - p_on_sin).norm() <= adaptive_tol);
            std::cout << "Result: intersects = " << r_on_sin.intersects << ", projected/closest point = ";
            for (const auto& pt : r_on_sin.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;

            Point<double,3> p_above_sin({M_PI/2, M_PI/2, 2.0});
            auto r_above_sin = intersect(p_above_sin, sinSurface);
            printTest("Point–SinSurface above surface",
                      !r_above_sin.intersects &&
                      r_above_sin.points.size() == 1);
            std::cout << "Result: intersects = " << r_above_sin.intersects << ", projected/closest point = ";
            for (const auto& pt : r_above_sin.points) {
                std::cout << "(" << pt[0] << ", " << pt[1] << ", " << pt[2] << ") ";
            }
            std::cout << std::endl;
        }
    }
} // end anonymous namespace

void testPointIntersection() {
    std::cout << "\nTesting Point Intersections\n\n";
    testPointPoint();
    testPointLine();
    testPointSegment();
    testPointCurve();
    testPointPlane();
    testPointFace();
    testPointSurface();
}

} // namespace Euclid::Tests
