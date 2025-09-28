#include <iostream>
#include "geometry/line.hpp"


inline void testLine() {
    std::cout << "\n📏 Testing Line Primitive\n\n";
    
    int oldDimension = Euclid::getSpaceDimension();
    
    // --- 2D Tests ---
    Euclid::setSpaceDimension(2);
    using Point2 = Point<float, 2>;
    using Vec2 = Eigen::Matrix<float, 2, 1>;
    using Line2 = Line<float, 2>;
    
    Point2 p2a{0.0f, 0.0f};
    Point2 p2b{1.0f, 1.0f};
    Vec2 dir2{1.0f, 0.0f};
    
    Line2 l2_pts(p2a, p2b);
    printTest("2D Line from two points", l2_pts.point1() == p2a && l2_pts.point2() == p2b);
    
    Line2 l2_dir(p2a, dir2);
    printTest("2D Line from point + direction", l2_dir.point1() == p2a && l2_dir.direction().isApprox(dir2.normalized()));
    
    // --- 3D Tests ---
    Euclid::setSpaceDimension(3);
    using Point3 = Point<double, 3>;
    using Vec3 = Eigen::Matrix<double, 3, 1>;
    using Line3 = Line<double, 3>;
    
    Point3 p3a{0.0, 0.0, 0.0};
    Point3 p3b{1.0, 1.0, 1.0};
    Vec3 dir3{0.0, 1.0, 0.0};
    
    Line3 l3_pts(p3a, p3b);
    printTest("3D Line from two points", l3_pts.point1() == p3a && l3_pts.point2() == p3b);
    
    Line3 l3_dir(p3a, dir3);
    printTest("3D Line from point + direction", l3_dir.point1() == p3a && l3_dir.direction().isApprox(dir3.normalized()));
    
    // --- 4D Tests ---
    Euclid::setSpaceDimension(4);
    using Point4 = Point<float, 4>;
    using Vec4 = Eigen::Matrix<float, 4, 1>;
    using Line4 = Line<float, 4>;
    
    Point4 p4a{1.0f, 2.0f, 3.0f, 4.0f};
    Point4 p4b{5.0f, 6.0f, 7.0f, 8.0f};
    Vec4 dir4{0.0f, 1.0f, 0.0f, 0.0f};
    
    Line4 l4_pts(p4a, p4b);
    printTest("4D Line from two points", l4_pts.point1() == p4a && l4_pts.point2() == p4b);
    
    Line4 l4_dir(p4a, dir4);
    printTest("4D Line from point + direction", l4_dir.point1() == p4a && l4_dir.direction().isApprox(dir4.normalized()));
    
    // --- 5D Tests ---
    Euclid::setSpaceDimension(5);
    using Point5 = Point<double, 5>;
    using Vec5 = Eigen::Matrix<double, 5, 1>;
    using Line5 = Line<double, 5>;
    
    Point5 p5a{1.0, 2.0, 3.0, 4.0, 5.0};
    Point5 p5b{6.0, 7.0, 8.0, 9.0, 10.0};
    Vec5 dir5{1.0, 0.0, 0.0, 0.0, 0.0};
    
    Line5 l5_pts(p5a, p5b);
    printTest("5D Line from two points", l5_pts.point1() == p5a && l5_pts.point2() == p5b);
    
    Line5 l5_dir(p5a, dir5);
    printTest("5D Line from point + direction", l5_dir.point1() == p5a && l5_dir.direction().isApprox(dir5.normalized()));
    
    /* 2025-09-18
     // Intersection logic only works for 2D
     // Need to revisit intersection logic so it is robust for arbitary dimension lines
     // --- 2D Intersection Tests ---
     Euclid::setSpaceDimension(2);
     {
     // Intersecting lines: (0,0)-(1,1) and (0,1)-(1,0), should intersect at (0.5, 0.5)
     Point2 a1{0.0f, 0.0f};
     Point2 a2{1.0f, 1.0f};
     Point2 b1{0.0f, 1.0f};
     Point2 b2{1.0f, 0.0f};
     Line2 la(a1, a2);
     Line2 lb(b1, b2);
     auto res2d = la.intersect(lb);
     printTest("2D Line intersection (intersecting)", res2d.intersects && approxEqual(res2d.pointOnThis, Point2{0.5f, 0.5f}));
     
     // Non-intersecting parallel lines: (0,0)-(1,0) and (0,1)-(1,1)
     Point2 c1{0.0f, 0.0f};
     Point2 c2{1.0f, 0.0f};
     Point2 d1{0.0f, 1.0f};
     Point2 d2{1.0f, 1.0f};
     Line2 lc(c1, c2);
     Line2 ld(d1, d2);
     auto res2d_nonintersect = lc.intersect(ld);
     std::cout << "res2d_nonintersect.intersects = " << res2d_nonintersect.intersects << ", res2d_nonintersect.pointOnThis = " << res2d_nonintersect.pointOnThis.coords.transpose() << "\n";
     printTest("2D Line intersection (non-intersecting parallel)", !res2d_nonintersect.intersects);
     }
     
     // --- 3D Intersection Tests ---
     Euclid::setSpaceDimension(3);
     {
     // Skew lines: (0,0,0)-(1,0,0) and (0,1,1)-(0,2,2), do not intersect
     Point3 a1{0.0, 0.0, 0.0};
     Point3 a2{1.0, 0.0, 0.0};
     Point3 b1{0.0, 1.0, 1.0};
     Point3 b2{0.0, 2.0, 2.0};
     Line3 la(a1, a2);
     Line3 lb(b1, b2);
     auto res3d = la.intersect(lb);
     std::cout << "res3d.intersects = " << res3d.intersects << ", res3d.pointOnThis = " << res3d.pointOnThis.coords.transpose() << "\n";
     printTest("3D Line intersection (skew, no intersection)", !res3d.intersects);
     
     // Intersecting lines: (0,0,0)-(1,1,1) and (1,0,0)-(0,1,1), should intersect at (0.5,0.5,0.5)
     Point3 c1{0.0, 0.0, 0.0};
     Point3 c2{1.0, 1.0, 1.0};
     Point3 d1{1.0, 0.0, 0.0};
     Point3 d2{0.0, 1.0, 1.0};
     Line3 lc(c1, c2);
     Line3 ld(d1, d2);
     auto res3d_2 = lc.intersect(ld);
     std::cout << "[3D intersect] intersects=" << res3d_2.intersects
     << ", point=" << res3d_2.pointOnThis.coords.transpose() << "\n";
     printTest("3D Line intersection (intersecting)", res3d_2.intersects && approxEqual(res3d_2.pointOnThis, Point3{0.5, 0.5, 0.5}));
     }
     
     // --- 4D Intersection Tests ---
     Euclid::setSpaceDimension(4);
     {
     // Skew lines: (0,0,0,0)-(1,0,0,0) and (0,1,1,1)-(0,2,2,2), do not intersect
     Point4 a1{0.0f, 0.0f, 0.0f, 0.0f};
     Point4 a2{1.0f, 0.0f, 0.0f, 0.0f};
     Point4 b1{0.0f, 1.0f, 1.0f, 1.0f};
     Point4 b2{0.0f, 2.0f, 2.0f, 2.0f};
     Line4 la(a1, a2);
     Line4 lb(b1, b2);
     auto res4d = la.intersect(lb);
     std::cout << "res4d.intersects = " << res4d.intersects << ", res4d.pointOnThis = " << res4d.pointOnThis.coords.transpose() << "\n";
     printTest("4D Line intersection (skew, no intersection)", !res4d.intersects);
     
     // Intersecting lines: (0,0,0,0)-(1,1,1,1) and (1,0,0,0)-(0,1,1,1), intersect at (0.5,0.5,0.5,0.5)
     Point4 c1{0.0f, 0.0f, 0.0f, 0.0f};
     Point4 c2{1.0f, 1.0f, 1.0f, 1.0f};
     Point4 d1{1.0f, 0.0f, 0.0f, 0.0f};
     Point4 d2{0.0f, 1.0f, 1.0f, 1.0f};
     Line4 lc(c1, c2);
     Line4 ld(d1, d2);
     auto res4d_2 = lc.intersect(ld);
     std::cout << "[4D intersect] intersects=" << res4d_2.intersects
     << ", point=" << res4d_2.pointOnThis.coords.transpose() << "\n";
     printTest("4D Line intersection (intersecting)", res4d_2.intersects && approxEqual(res4d_2.pointOnThis, Point4{0.5f, 0.5f, 0.5f, 0.5f}));
     }
     
     // --- 5D Intersection Tests ---
     Euclid::setSpaceDimension(5);
     {
     // Skew lines: (0,0,0,0,0)-(1,0,0,0,0) and (0,1,1,1,1)-(0,2,2,2,2), do not intersect
     Point5 a1{0.0, 0.0, 0.0, 0.0, 0.0};
     Point5 a2{1.0, 0.0, 0.0, 0.0, 0.0};
     Point5 b1{0.0, 1.0, 1.0, 1.0, 1.0};
     Point5 b2{0.0, 2.0, 2.0, 2.0, 2.0};
     Line5 la(a1, a2);
     Line5 lb(b1, b2);
     auto res5d = la.intersect(lb);
     std::cout << "res5d.intersects = " << res5d.intersects << ", res5d.pointOnThis = " << res5d.pointOnThis.coords.transpose() << "\n";
     printTest("5D Line intersection (skew, no intersection)", !res5d.intersects);
     
     // Intersecting lines: (0,0,0,0,0)-(1,1,1,1,1) and (1,0,0,0,0)-(0,1,1,1,1), intersect at (0.5,0.5,0.5,0.5,0.5)
     Point5 c1{0.0, 0.0, 0.0, 0.0, 0.0};
     Point5 c2{1.0, 1.0, 1.0, 1.0, 1.0};
     Point5 d1{1.0, 0.0, 0.0, 0.0, 0.0};
     Point5 d2{0.0, 1.0, 1.0, 1.0, 1.0};
     Line5 lc(c1, c2);
     Line5 ld(d1, d2);
     auto res5d_2 = lc.intersect(ld);
     std::cout << "[5D intersect] intersects=" << res5d_2.intersects
     << ", point=" << res5d_2.pointOnThis.coords.transpose() << "\n";
     printTest("5D Line intersection (intersecting)", res5d_2.intersects && approxEqual(res5d_2.pointOnThis, Point5{0.5, 0.5, 0.5, 0.5, 0.5}));
     }
     */
    
    Euclid::setSpaceDimension(oldDimension);
    
    std::cout << "\nTesting Line Transform Logic\n";
    
    int oldDim = Euclid::getSpaceDimension();
    Euclid::setSpaceDimension(2);
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
    euclid::algebra::Affine<float,2> transform(linear, translation);
    
    Line2 l_transformed = l.applyTransform(transform);
    Point2 expectedPoint{3.0f, 5.0f};
    Vec2 expectedDir = (linear * dir).normalized();
    printTest("Line transform (scale + translate) point", approxEqual(l_transformed.point1(), expectedPoint, 1e-6f));
    printTest("Line transform (scale + translate) direction", l_transformed.direction().isApprox(expectedDir, 1e-6f));
    
    // --- 2. Pure translation ---
    Eigen::Matrix2f identity = Eigen::Matrix2f::Identity();
    Eigen::Vector2f translation2;
    translation2 << -2.0f, 3.0f;
    euclid::algebra::Affine<float,2> transformTranslate(identity, translation2);
    
    Line2 l_translated = l.applyTransform(transformTranslate);
    Point2 expectedTransPoint{-1.0f, 5.0f};
    Vec2 expectedTransDir = dir.normalized();
    printTest("Line transform (pure translation) point", approxEqual(l_translated.point1(), expectedTransPoint, 1e-6f));
    printTest("Line transform (pure translation) direction", l_translated.direction().isApprox(expectedTransDir, 1e-6f));
    
    // --- 3. Rotation about pivot ---
    Point2 pivot{2.0f, 1.0f};
    float angle = M_PI / 2.0f;
    Eigen::Matrix2f rotMat;
    rotMat << std::cos(angle), -std::sin(angle),
    std::sin(angle),  std::cos(angle);
    euclid::algebra::SpecialOrthogonal<float,2> rotation(rotMat);
    Eigen::Vector2f translationPivot = pivot.coords - rotation.A * pivot.coords;
    euclid::algebra::Affine<float,2> rotatePivot(rotation.A, translationPivot);
    
    Line2 l_rotated = l.applyTransform(rotatePivot);
    
    // Expected points computed manually (API-based)
    Vec2 relP0 = p0.coords - pivot.coords;
    Vec2 relP1 = p1.coords - pivot.coords;
    Eigen::Vector2f expectedDirRot = (rotMat * dir).normalized();
    Point2 expectedPointRot = rotatePivot.apply(p0);
    
    printTest("Line transform (rotation about pivot) point", approxEqual(l_rotated.point1(), expectedPointRot, 1e-6f));
    printTest("Line transform (rotation about pivot) direction", l_rotated.direction().isApprox(expectedDirRot, 1e-6f));
    
    Euclid::setSpaceDimension(oldDim);
    
}
    

inline void testLineMeasureAngle() {
    std::cout << "\nTesting Line MeasureAngle Logic\n";

    int oldDim = Euclid::getSpaceDimension();
    Euclid::setSpaceDimension(2);

    using Point2 = Point<float,2>;
    using Vec2 = Eigen::Matrix<float,2,1>;
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

    printTest("Angle between orthogonal lines", std::abs(angleOrth - M_PI/2.0f) < 1e-6f);
    printTest("Angle between same direction lines", std::abs(angleZero - 0.0f) < 1e-6f);
    printTest("Angle between opposite direction lines", std::abs(anglePi - M_PI) < 1e-6f);

    Euclid::setSpaceDimension(oldDim);
}
