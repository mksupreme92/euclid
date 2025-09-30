#pragma once

#include <iostream>
#include <array>
#include "geometry/point.hpp"
#include "geometry/plane.hpp"
#include "test_utilities.hpp"

using namespace euclid::geometry;
using namespace euclid::tests;

inline void testPlane() {
    std::cout << "\n⧉ Testing Plane Primitive\n";


    // --- 2D Plane (line) ---
    using Point2 = Point<float,2>;
    using Plane2 = Plane<float,2>;

    Point2 p0{0.0f,0.0f};
    Point2 p1{1.0f,0.0f};
    std::array<Point2,2> points2 = {p0,p1};
    auto maybePlane2 = Plane2::fromPoints(points2);
    if (maybePlane2) {
        Plane2 plane2 = *maybePlane2;
        printTest("2D Plane contains p0", plane2.contains(p0));
        printTest("2D Plane contains p1", plane2.contains(p1));
        printTest("2D Plane does not contain off-line point", !plane2.contains(Point2{0.5f,1.0f}));

        // Signed distance tests for 2D
        float dist_p0 = plane2.signedDistance(p0);
        float dist_p1 = plane2.signedDistance(p1);
        float dist_off = plane2.signedDistance(Point2{0.5f,1.0f});
        printTest("2D Plane signed distance is zero on p0", std::abs(dist_p0) < 1e-6f);
        printTest("2D Plane signed distance is zero on p1", std::abs(dist_p1) < 1e-6f);
        printTest("2D Plane signed distance is non-zero off-line point", std::abs(dist_off) > 1e-6f);
    } else {
        std::cerr << "Failed to construct 2D plane from points\n";
    }

    // --- 3D Plane ---
    using Point3 = Point<float,3>;
    using Plane3 = Plane<float,3>;

    Point3 pa{0.0f,0.0f,0.0f};
    Point3 pb{1.0f,0.0f,0.0f};
    Point3 pc{0.0f,1.0f,0.0f};

    std::array<Point3,3> points3 = {pa,pb,pc};
    auto maybePlane3 = Plane3::fromPoints(points3);
    if (maybePlane3) {
        Plane3 plane3 = *maybePlane3;
        printTest("3D Plane contains pa", plane3.contains(pa));
        printTest("3D Plane contains pb", plane3.contains(pb));
        printTest("3D Plane contains pc", plane3.contains(pc));
        printTest("3D Plane does not contain off-plane point", !plane3.contains(Point3{0.0f,0.0f,1.0f}));

        // Signed distance tests for 3D
        float dist_pa = plane3.signedDistance(pa);
        float dist_pb = plane3.signedDistance(pb);
        float dist_pc = plane3.signedDistance(pc);
        float dist_off = plane3.signedDistance(Point3{0.0f,0.0f,1.0f});
        printTest("3D Plane signed distance is zero on pa", std::abs(dist_pa) < 1e-6f);
        printTest("3D Plane signed distance is zero on pb", std::abs(dist_pb) < 1e-6f);
        printTest("3D Plane signed distance is zero on pc", std::abs(dist_pc) < 1e-6f);
        printTest("3D Plane signed distance is non-zero off-plane point", std::abs(dist_off) > 1e-6f);
    } else {
        std::cerr << "Failed to construct 3D plane from points\n";
    }

    // --- 4D Hyperplane ---
    using Point4 = Point<float,4>;
    using Plane4 = Plane<float,4>;

    Point4 p4a{0.0f,0.0f,0.0f,0.0f};
    Point4 p4b{1.0f,0.0f,0.0f,0.0f};
    Point4 p4c{0.0f,1.0f,0.0f,0.0f};
    Point4 p4d{0.0f,0.0f,1.0f,0.0f};

    std::array<Point4,4> points4 = {p4a,p4b,p4c,p4d};
    auto maybePlane4 = Plane4::fromPoints(points4);
    if (maybePlane4) {
        Plane4 plane4 = *maybePlane4;
        printTest("4D Plane contains p4a", plane4.contains(p4a));
        printTest("4D Plane contains p4b", plane4.contains(p4b));
        printTest("4D Plane contains p4c", plane4.contains(p4c));
        printTest("4D Plane contains p4d", plane4.contains(p4d));
        printTest("4D Plane does not contain off-plane point", !plane4.contains(Point4{0.0f,0.0f,0.0f,1.0f}));

        // Signed distance tests for 4D
        float dist_p4a = plane4.signedDistance(p4a);
        float dist_p4b = plane4.signedDistance(p4b);
        float dist_p4c = plane4.signedDistance(p4c);
        float dist_p4d = plane4.signedDistance(p4d);
        float dist_off = plane4.signedDistance(Point4{0.0f,0.0f,0.0f,1.0f});
        printTest("4D Plane signed distance is zero on p4a", std::abs(dist_p4a) < 1e-6f);
        printTest("4D Plane signed distance is zero on p4b", std::abs(dist_p4b) < 1e-6f);
        printTest("4D Plane signed distance is zero on p4c", std::abs(dist_p4c) < 1e-6f);
        printTest("4D Plane signed distance is zero on p4d", std::abs(dist_p4d) < 1e-6f);
        printTest("4D Plane signed distance is non-zero off-plane point", std::abs(dist_off) > 1e-6f);
    } else {
        std::cerr << "Failed to construct 4D plane from points\n";
    }

    // --- 5D Hyperplane ---
    using Point5 = Point<float,5>;
    using Plane5 = Plane<float,5>;

    Point5 p5a{0,0,0,0,0};
    Point5 p5b{1,0,0,0,0};
    Point5 p5c{0,1,0,0,0};
    Point5 p5d{0,0,1,0,0};
    Point5 p5e{0,0,0,1,0};

    std::array<Point5,5> points5 = {p5a,p5b,p5c,p5d,p5e};
    auto maybePlane5 = Plane5::fromPoints(points5);
    if (maybePlane5) {
        Plane5 plane5 = *maybePlane5;
        printTest("5D Plane contains p5a", plane5.contains(p5a));
        printTest("5D Plane contains p5b", plane5.contains(p5b));
        printTest("5D Plane contains p5c", plane5.contains(p5c));
        printTest("5D Plane contains p5d", plane5.contains(p5d));
        printTest("5D Plane contains p5e", plane5.contains(p5e));
        printTest("5D Plane does not contain off-plane point", !plane5.contains(Point5{0,0,0,0,1}));

        // Signed distance tests for 5D
        float dist_p5a = plane5.signedDistance(p5a);
        float dist_p5b = plane5.signedDistance(p5b);
        float dist_p5c = plane5.signedDistance(p5c);
        float dist_p5d = plane5.signedDistance(p5d);
        float dist_p5e = plane5.signedDistance(p5e);
        float dist_off = plane5.signedDistance(Point5{0,0,0,0,1});
        printTest("5D Plane signed distance is zero on p5a", std::abs(dist_p5a) < 1e-6f);
        printTest("5D Plane signed distance is zero on p5b", std::abs(dist_p5b) < 1e-6f);
        printTest("5D Plane signed distance is zero on p5c", std::abs(dist_p5c) < 1e-6f);
        printTest("5D Plane signed distance is zero on p5d", std::abs(dist_p5d) < 1e-6f);
        printTest("5D Plane signed distance is zero on p5e", std::abs(dist_p5e) < 1e-6f);
        printTest("5D Plane signed distance is non-zero off-plane point", std::abs(dist_off) > 1e-6f);
    } else {
        std::cerr << "Failed to construct 5D plane from points\n";
    }


    // --- Plane from point + direction (fixed 3D) ---
    using Point3 = Point<float,3>;
    using Plane3 = Plane<float,3>;
    Eigen::Vector3f dir3;
    dir3 << 0.0f, 0.0f, 1.0f;
    Point3 base3{0.0f, 0.0f, 0.0f};
    auto maybePlaneFromPointNormal = Plane3::fromPointAndNormal(base3, dir3);
    if (maybePlaneFromPointNormal) {
        Plane3 plane3 = *maybePlaneFromPointNormal;
        printTest("Plane3 (from point+normal) contains base3", plane3.contains(base3));
        printTest("Plane3 (from point+normal) does not contain off-plane point",
                  !plane3.contains(Point3{1.0f, 0.0f, 1.0f}));
    } else {
        std::cerr << "Failed to construct Plane3 from point + direction\n";
    }

    // --- Plane from point + direction (2D) ---
    using Point2 = Point<float,2>;
    using Plane2 = Plane<float,2>;
    Eigen::Vector2f dir2;
    dir2 << 0.0f, 1.0f;
    Point2 base2{0.0f, 0.0f};
    auto maybePlaneFromPointNormal2 = Plane2::fromPointAndNormal(base2, dir2);
    if (maybePlaneFromPointNormal2) {
        Plane2 plane2 = *maybePlaneFromPointNormal2;
        printTest("Plane2 (from point+normal) contains base2", plane2.contains(base2));
        printTest("Plane2 (from point+normal) does not contain off-plane point",
                  !plane2.contains(Point2{1.0f, 1.0f}));
    } else {
        std::cerr << "Failed to construct Plane2 from point + direction\n";
    }

    // --- Plane from point + direction (4D) ---
    using Point4 = Point<float,4>;
    using Plane4 = Plane<float,4>;
    Eigen::Vector4f dir4;
    dir4 << 0.0f, 0.0f, 1.0f, 0.0f;
    Point4 base4{0.0f, 0.0f, 0.0f, 0.0f};
    auto maybePlaneFromPointNormal4 = Plane4::fromPointAndNormal(base4, dir4);
    if (maybePlaneFromPointNormal4) {
        Plane4 plane4 = *maybePlaneFromPointNormal4;
        printTest("Plane4 (from point+normal) contains base4", plane4.contains(base4));
        printTest("Plane4 (from point+normal) does not contain off-plane point",
                  !plane4.contains(Point4{0.0f, 0.0f, 1.0f, 1.0f}));
    } else {
        std::cerr << "Failed to construct Plane4 from point + direction\n";
    }

    // --- Plane from point + direction (5D) ---
    using Point5 = Point<float,5>;
    using Plane5 = Plane<float,5>;
    Eigen::Matrix<float,5,1> dir5;
    dir5 << 0.0f, 0.0f, 0.0f, 1.0f, 0.0f;
    Point5 base5{0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    auto maybePlaneFromPointNormal5 = Plane5::fromPointAndNormal(base5, dir5);
    if (maybePlaneFromPointNormal5) {
        Plane5 plane5 = *maybePlaneFromPointNormal5;
        printTest("Plane5 (from point+normal) contains base5", plane5.contains(base5));
        printTest("Plane5 (from point+normal) does not contain off-plane point",
                  !plane5.contains(Point5{0.0f, 0.0f, 0.0f, 1.0f, 1.0f}));
    } else {
        std::cerr << "Failed to construct Plane5 from point + direction\n";
    }

    // --- Plane Transform Tests ---
    std::cout << "\nTesting Plane Transformations\n";
    using Point3 = Point<float,3>;
    using Plane3 = Plane<float,3>;

    std::array<Point3,3> pts3 = { Point3{0,0,0}, Point3{1,0,0}, Point3{0,1,0} };
    auto maybePlane = Plane3::fromPoints(pts3);
    if (maybePlane) {
        Plane3 plane = *maybePlane;

        // Define affine transform: rotation + translation
        Eigen::AngleAxisf rotation(M_PI/2.0f, Eigen::Vector3f::UnitZ());
        Eigen::Vector3f translation(1.0f, 2.0f, 3.0f);
        euclid::algebra::Affine<float,3> affineTransform(rotation.toRotationMatrix(), translation);

        Plane3 transformedPlane = plane.applyTransform(affineTransform);

        for (const auto& p : pts3) {
            Point3 pTrans = affineTransform.apply(p);
            printTest("Transformed plane contains transformed point", transformedPlane.contains(pTrans));
            float sd = transformedPlane.signedDistance(pTrans);
            printTest("Signed distance to transformed plane is near zero", std::abs(sd) < 1e-5f);
        }
    } else {
        std::cerr << "Failed to construct Plane3 for transform test\n";
    }

}
