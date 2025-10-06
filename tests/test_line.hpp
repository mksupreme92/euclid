#include <iostream>
#include "geometry/line.hpp"
#include "algebra/transform.hpp"
#include "test_line_intersection.hpp"
#include "test_utilities.hpp"

using namespace Euclid::Geometry;
using namespace Euclid::Algebra;
using namespace Euclid::Tests;


inline void testLineMeasureAngle() {
    std::cout << "\nTesting Line MeasureAngle Logic\n";


    using Point2 = Point<float,2>;
    using Line2 = Line<float,2>;

    Point2 p0{0.0f, 0.0f};
    Point2 p1{1.0f, 0.0f};
    Point2 p2{0.0f, 1.0f};
    Point2 p3{2.0f, 0.0f};

    Line2 l_horiz(p0, p1);   // horizontal
    Line2 l_vert(p0, p2);    // vertical
    Line2 l_sameDir(p0, p3); // same direction as l_horiz
    Line2 l_oppDir(p1, p0);  // opposite direction

    float angleOrth = l_horiz.measureAngle(l_vert);
    float angleZero = l_horiz.measureAngle(l_sameDir);
    float anglePi   = l_horiz.measureAngle(l_oppDir);

    float tol = Euclid::Tolerance().evaluateEpsilon(M_PI/2.0f);
    printTest("Angle between orthogonal lines", std::abs(angleOrth - M_PI/2.0f) <= tol);

    tol = Euclid::Tolerance().evaluateEpsilon(0.0f);
    printTest("Angle between same direction lines", std::abs(angleZero - 0.0f) <= tol);

    tol = Euclid::Tolerance().evaluateEpsilon(M_PI);
    printTest("Angle between opposite direction lines", std::abs(anglePi - M_PI) <= tol);

}

inline void testLine() {
    std::cout << "\n⟷ Testing Line Primitive\n\n";
    
    
    // --- 2D Tests ---
    using Point2 = Point<float, 2>;
    using Vec2 = Eigen::Matrix<float, 2, 1>;
    using Line2 = Line<float, 2>;
    
    Point2 p2a{0.0f, 0.0f};
    Point2 p2b{1.0f, 1.0f};
    Vec2 dir2{1.0f, 0.0f};
    
    Line2 l2_pts(p2a, p2b);
    printTest("2D Line from two points", l2_pts.point1() == p2a && l2_pts.point2() == p2b);
    
    Line2 l2_dir(p2a, dir2);
    printTest("2D Line from point + direction", l2_dir.point1() == p2a && Euclid::equalWithinTolerance(l2_dir.direction(), dir2.normalized()));

    // --- 2D Intersection Tests ---
    {
        // Line2 l2_1(Point2{0.0f, 0.0f}, Vec2{1.0f, 1.0f});
        // Line2 l2_2(Point2{0.0f, 1.0f}, Vec2{1.0f, -1.0f});
        // Line2::LineIntersection result = l2_1.intersect(l2_2);
        // printTest("2D Lines intersect", result.intersects);
        // if (result.intersects) {
        //     printTest("2D Intersection point on line 1", result.pointOnThis == result.pointOnThis);
        //     printTest("2D Intersection point on line 2", result.pointOnOther == result.pointOnOther);
        // }

        // Parallel lines
        // Line2 l2_3(Point2{0.0f, 0.0f}, Vec2{1.0f, 0.0f});
        // Line2 l2_4(Point2{0.0f, 1.0f}, Vec2{1.0f, 0.0f});
        // result = l2_3.intersect(l2_4);
        // printTest("2D Parallel lines do not intersect", !result.intersects);
    }
    
    // --- 3D Tests ---
    using Point3 = Point<double, 3>;
    using Vec3 = Eigen::Matrix<double, 3, 1>;
    using Line3 = Line<double, 3>;
    
    Point3 p3a{0.0, 0.0, 0.0};
    Point3 p3b{1.0, 1.0, 1.0};
    Vec3 dir3{0.0, 1.0, 0.0};
    
    Line3 l3_pts(p3a, p3b);
    printTest("3D Line from two points", l3_pts.point1() == p3a && l3_pts.point2() == p3b);
    
    Line3 l3_dir(p3a, dir3);
    printTest("3D Line from point + direction", l3_dir.point1() == p3a && Euclid::equalWithinTolerance(l3_dir.direction(), dir3.normalized()));

    // --- 3D Intersection Tests ---
    {
        // Line3 l3_1(Point3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
        // Line3 l3_2(Point3{0.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});
        // Line3::LineIntersection result = l3_1.intersect(l3_2);
        // printTest("3D Lines intersect (at origin)", result.intersects);
        // if (result.intersects) {
        //     printTest("3D Intersection point on line 1", result.pointOnThis == result.pointOnThis);
        //     printTest("3D Intersection point on line 2", result.pointOnOther == result.pointOnOther);
        //     printTest("3D Intersection points equal", result.pointOnThis == result.pointOnOther);
        // }

        // Skew lines (no intersection)
        // Line3 l3_3(Point3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
        // Line3 l3_4(Point3{0.0, 1.0, 1.0}, Vec3{0.0, 1.0, 0.0});
        // result = l3_3.intersect(l3_4);
        // printTest("3D Skew lines do not intersect", !result.intersects);
    }
    
    // --- 4D Tests ---
    using Point4 = Point<float, 4>;
    using Vec4 = Eigen::Matrix<float, 4, 1>;
    using Line4 = Line<float, 4>;
    
    Point4 p4a{1.0f, 2.0f, 3.0f, 4.0f};
    Point4 p4b{5.0f, 6.0f, 7.0f, 8.0f};
    Vec4 dir4{0.0f, 1.0f, 0.0f, 0.0f};
    
    Line4 l4_pts(p4a, p4b);
    printTest("4D Line from two points", l4_pts.point1() == p4a && l4_pts.point2() == p4b);
    
    Line4 l4_dir(p4a, dir4);
    printTest("4D Line from point + direction", Euclid::equalWithinTolerance(l4_dir.direction(), dir4.normalized()) && l4_dir.point1() == p4a);

    // --- 4D Intersection Tests ---
    {
        // Line4 l4_1(Point4{0.0f, 0.0f, 0.0f, 0.0f}, Vec4{1.0f, 0.0f, 0.0f, 0.0f});
        // Line4 l4_2(Point4{0.0f, 0.0f, 0.0f, 0.0f}, Vec4{0.0f, 1.0f, 0.0f, 0.0f});
        // Line4::LineIntersection result = l4_1.intersect(l4_2);
        // printTest("4D Lines intersect (at origin)", result.intersects);
        // if (result.intersects) {
        //     printTest("4D Intersection point on line 1", result.pointOnThis == result.pointOnThis);
        //     printTest("4D Intersection point on line 2", result.pointOnOther == result.pointOnOther);
        //     printTest("4D Intersection points equal", result.pointOnThis == result.pointOnOther);
        // }

        // Parallel lines in 4D
        // Line4 l4_3(Point4{0.0f, 0.0f, 0.0f, 0.0f}, Vec4{1.0f, 0.0f, 0.0f, 0.0f});
        // Line4 l4_4(Point4{0.0f, 1.0f, 0.0f, 0.0f}, Vec4{1.0f, 0.0f, 0.0f, 0.0f});
        // result = l4_3.intersect(l4_4);
        // printTest("4D Parallel lines do not intersect", !result.intersects);
    }
    
    // --- 5D Tests ---
    using Point5 = Point<double, 5>;
    using Vec5 = Eigen::Matrix<double, 5, 1>;
    using Line5 = Line<double, 5>;
    
    Point5 p5a{1.0, 2.0, 3.0, 4.0, 5.0};
    Point5 p5b{6.0, 7.0, 8.0, 9.0, 10.0};
    Vec5 dir5{1.0, 0.0, 0.0, 0.0, 0.0};
    
    Line5 l5_pts(p5a, p5b);
    printTest("5D Line from two points", l5_pts.point1() == p5a && l5_pts.point2() == p5b);
    
    Line5 l5_dir(p5a, dir5);
    printTest("5D Line from point + direction", Euclid::equalWithinTolerance(l5_dir.direction(), dir5.normalized()) && l5_dir.point1() == p5a);

    // --- 5D Intersection Tests ---
    {
        // Line5 l5_1(Point5{0.0, 0.0, 0.0, 0.0, 0.0}, Vec5{1.0, 0.0, 0.0, 0.0, 0.0});
        // Line5 l5_2(Point5{0.0, 0.0, 0.0, 0.0, 0.0}, Vec5{0.0, 1.0, 0.0, 0.0, 0.0});
        // Line5::LineIntersection result = l5_1.intersect(l5_2);
        // printTest("5D Lines intersect (at origin)", result.intersects);
        // if (result.intersects) {
        //     printTest("5D Intersection point on line 1", result.pointOnThis == result.pointOnThis);
        //     printTest("5D Intersection point on line 2", result.pointOnOther == result.pointOnOther);
        //     printTest("5D Intersection points equal", result.pointOnThis == result.pointOnOther);
        // }

        // Parallel lines in 5D
        // Line5 l5_3(Point5{0.0, 0.0, 0.0, 0.0, 0.0}, Vec5{1.0, 0.0, 0.0, 0.0, 0.0});
        // Line5 l5_4(Point5{0.0, 1.0, 0.0, 0.0, 0.0}, Vec5{1.0, 0.0, 0.0, 0.0, 0.0});
        // result = l5_3.intersect(l5_4);
        // printTest("5D Parallel lines do not intersect", !result.intersects);
    }
    
    std::cout << "\nTesting Line Transform Logic\n";
    
    using Point2 = Point<float,2>;
    using Vec2 = Eigen::Matrix<float,2,1>;
    using Line2 = Line<float,2>;
    
    Point2 p0{1.0f, 2.0f};
    Point2 p1{3.0f, 5.0f};
    Vec2 dir = p1 - p0;
    Line2 l(p0, dir);
    
    // --- 1. Scale + translate ---
    Eigen::Matrix2f linear;
    linear << 2.0f, 0.0f,
    0.0f, 3.0f;
    Eigen::Vector2f translation;
    translation << 1.0f, -1.0f;
    Affine<float,2> transform(linear, translation);
    
    Line2 l_transformed = l.applyTransform(transform);
    Point2 expectedPoint{3.0f, 5.0f};
    Vec2 expectedDir = (linear * dir).normalized();
    printTest("Line transform (scale + translate) point", l_transformed.point1() == expectedPoint);
    printTest("Line transform (scale + translate) direction", Euclid::equalWithinTolerance(l_transformed.direction(), expectedDir));
    
    // --- 2. Pure translation ---
    Eigen::Matrix2f identity = Eigen::Matrix2f::Identity();
    Eigen::Vector2f translation2;
    translation2 << -2.0f, 3.0f;
    Affine<float,2> transformTranslate(identity, translation2);
    
    Line2 l_translated = l.applyTransform(transformTranslate);
    Point2 expectedTransPoint{-1.0f, 5.0f};
    Vec2 expectedTransDir = dir.normalized();
    printTest("Line transform (pure translation) point", l_translated.point1() == expectedTransPoint);
    printTest("Line transform (pure translation) direction", Euclid::equalWithinTolerance(l_translated.direction(), expectedTransDir));
    
    // --- 3. Rotation about pivot ---
    Point2 pivot{2.0f, 1.0f};
    float angle = M_PI / 2.0f;
    Eigen::Matrix2f rotMat;
    rotMat << std::cos(angle), -std::sin(angle),
    std::sin(angle),  std::cos(angle);
    SpecialOrthogonal<float,2> rotation(rotMat);
    Eigen::Vector2f translationPivot = pivot.coords - rotation.A * pivot.coords;
    Affine<float,2> rotatePivot(rotation.A, translationPivot);
    
    Line2 l_rotated = l.applyTransform(rotatePivot);
    
    // Expected points computed manually (API-based)
    Vec2 relP0 = p0.coords - pivot.coords;
    Vec2 relP1 = p1.coords - pivot.coords;
    Eigen::Vector2f expectedDirRot = (rotMat * dir).normalized();
    Point2 expectedPointRot = rotatePivot.apply(p0);
    
    printTest("Line transform (rotation about pivot) point", l_rotated.point1() == expectedPointRot);
    printTest("Line transform (rotation about pivot) direction", Euclid::equalWithinTolerance(l_rotated.direction(), expectedDirRot));
    
    testLineMeasureAngle();
    testLineIntersection();
    
}
    
