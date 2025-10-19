#pragma once

#include <iostream>
#include "geometry/face.hpp"
#include "geometry/point.hpp"

using namespace Euclid::Geometry;
using namespace Euclid::Algebra;
using namespace Euclid::Tests;

inline void testFaceTransforms() {
    std::cout << "\nTesting Face Transformations\n";


    // --- 3D Face ---
    using Point3 = Point<double,3>;
    using Face3  = Face<double,3>;

    std::vector<Point3> faceVerts3 = { {0,0,0}, {1,0,0}, {0,1,0} };
    Face3 face3 = Face3::fromPoints(faceVerts3);

    Eigen::Matrix3d rot3;
    rot3 << 0.0, -1.0, 0.0,
            1.0,  0.0, 0.0,
            0.0,  0.0, 1.0;
    Eigen::Vector3d trans3(1.0, 2.0, 3.0);
    Affine<double,3> transform3(rot3, trans3);

    Face3 face3Trans = face3.applyTransform(transform3);

    bool verticesOk3 = true;
    for (size_t i = 0; i < faceVerts3.size(); ++i) {
        Point3 expected(faceVerts3[i].coords);
        expected.coords = rot3 * expected.coords + trans3;
        verticesOk3 &= face3Trans.vertices[i].isEqual(expected, Euclid::Tolerance());
    }

    Eigen::Vector3d expectedNormal3 = rot3 * face3.normal;
    expectedNormal3.normalize();
    bool normalOk3 = Euclid::equalWithinTolerance(face3Trans.normal, expectedNormal3, Euclid::Tolerance());

    printTest("Face3 transform: vertices", verticesOk3);
    printTest("Face3 transform: normal", normalOk3);

    // --- 5D Face ---
    using Point5 = Point<double,5>;
    using Face5  = Face<double,5>;

    std::vector<Point5> faceVerts5 = { {0,0,0,0,0}, {1,0,0,0,0}, {0,1,0,0,0}, {0,0,1,0,0}, {0,0,0,1,0} };
    Face5 face5 = Face5::fromPoints(faceVerts5);

    Eigen::Matrix<double,5,5> linear5 = Eigen::Matrix<double,5,5>::Identity();
    linear5(0,0) = 2.0; linear5(1,1) = 0.5; linear5(2,2) = -1.0; linear5(3,3) = 1.5; linear5(4,4) = 0.25;
    Eigen::Matrix<double,5,1> trans5;
    trans5 << 1.0, -1.0, 0.5, 2.0, -0.5;

    Affine<double,5> transform5(linear5, trans5);

    Face5 face5Trans = face5.applyTransform(transform5);

    bool verticesOk5 = true;
    for (size_t i = 0; i < faceVerts5.size(); ++i) {
        Point5 expected(faceVerts5[i].coords);
        expected.coords = linear5 * expected.coords + trans5;
        verticesOk5 &= face5Trans.vertices[i].isEqual(expected, Euclid::Tolerance());
    }

    bool allInside = true;
    for (const auto& v : face5Trans.vertices)
        allInside &= intersect(v, face5Trans).intersects;

    printTest("Face5 transform: vertices", verticesOk5);
    printTest("Face5 transform: plane consistency", allInside);

}

template<typename T>
void testFace2D() {

    using Point2 = Point<T, 2>;
    using Face2  = Face<T, 2>;

    std::vector<Point2> pts2 = { Point2{0,0}, Point2{1,0}, Point2{1,1}, Point2{0,1} };
    Face2 face2 = Face2::fromPoints(pts2);

    printTest("2D Face: vertex inside", intersect(pts2[0], face2).intersects);
    printTest("2D Face: edge vertex inside", intersect(pts2[1], face2).intersects);
    printTest("2D Face: outside point", !intersect(Point2{0.5, 1.5}, face2).intersects);
}

template<typename T>
void testFace3D() {

    using Point3 = Point<T, 3>;
    using Face3  = Face<T, 3>;

    std::vector<Point3> pts3 = { Point3{0,0,0}, Point3{1,0,0}, Point3{0,1,0} };
    Face3 face3 = Face3::fromPoints(pts3);

    printTest("3D Face: vertex inside", intersect(pts3[0], face3).intersects);
    printTest("3D Face: vertex inside", intersect(pts3[1], face3).intersects);
    printTest("3D Face: vertex inside", intersect(pts3[2], face3).intersects);
    printTest("3D Face: outside point", !intersect(Point3{0,0,1}, face3).intersects);
}

template<typename T>
void testFace5D() {

    using Point5 = Point<T, 5>;
    using Face5  = Face<T, 5>;

    std::vector<Point5> pts5 = { Point5{0,0,0,0,0}, Point5{1,0,0,0,0}, Point5{0,1,0,0,0}, Point5{0,0,1,0,0}, Point5{0,0,0,1,0} };
    Face5 face5 = Face5::fromPoints(pts5);

    printTest("5D Face: vertex inside", intersect(pts5[0], face5).intersects);
    printTest("5D Face: vertex inside", intersect(pts5[1], face5).intersects);
    printTest("5D Face: vertex inside", intersect(pts5[2], face5).intersects);
    printTest("5D Face: vertex inside", intersect(pts5[3], face5).intersects);
    printTest("5D Face: vertex inside", intersect(pts5[4], face5).intersects);
    printTest("5D Face: outside point", !intersect(Point5{0,0,0,0,1}, face5).intersects);
}

inline void testFace() {
    std::cout << "\n⛛ Testing Face Primitive\n";
    testFace2D<double>();
    testFace3D<double>();
    testFace5D<double>();

    // --- Face Transformation Tests ---
    testFaceTransforms();
}
