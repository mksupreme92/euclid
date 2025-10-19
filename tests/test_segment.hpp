#include <iostream>
#include "geometry/segment.hpp"
#include "algebra/transform.hpp"
#include "test_utilities.hpp"
#include "test_segment_intersection.hpp"

using namespace Euclid::Geometry;
using namespace Euclid::Algebra;
using namespace Euclid::Tests;

inline void testSegment() {
    std::cout << "\n⎯ Testing Segment Primitive\n\n";

    // --- 2D Tests ---
    using Point2 = Point<float,2>;
    using Segment2 = Segment<float,2>;

    Point2 p0{0.0f, 0.0f};
    Point2 p1{1.0f, 0.0f};
    Point2 p2{0.5f, 0.0f};
    Point2 p3{2.0f, 0.0f};
    Point2 pOff{0.5f, 0.5f};

    Segment2 seg2(p0, p1);
    printTest("2D Segment length", std::abs(seg2.length() - 1.0f) <= Euclid::Tolerance().evaluateEpsilon(1));
    printTest("2D Segment midpoint", seg2.midpoint() == Point2{0.5f,0.0f});
    printTest("2D Segment vector", seg2.vector() == p1.coords - p0.coords);

    printTest("Point on segment", seg2.contains(Point2{p2}));
    printTest("Point before segment start", !seg2.contains(Point2{p0 - Point2{1.0f,0.0f}}));
    printTest("Point after segment end", !seg2.contains(Point2{p3}));
    printTest("Point off segment line", !seg2.contains(Point2{pOff}));

    // --- 3D Tests ---
    using Point3 = Point<double,3>;
    using Segment3 = Segment<double,3>;

    Point3 q0{0.0,0.0,0.0};
    Point3 q1{0.0,1.0,0.0};
    Point3 qMid{0.0,0.5,0.0};
    Segment3 seg3(q0,q1);

    printTest("3D Segment length", std::abs(seg3.length() - 1.0) <= Euclid::Tolerance().evaluateEpsilon(1));
    printTest("3D Segment midpoint", seg3.midpoint() == qMid);
    printTest("3D Point on segment", seg3.contains(Point3{qMid}));
    printTest("3D Point off segment", !seg3.contains(Point3{0.0,1.5,0.0}));

    // --- 4D Tests ---
    using Point4 = Point<float,4>;
    using Segment4 = Segment<float,4>;
    Point4 r0{1.0f,2.0f,3.0f,4.0f};
    Point4 r1{2.0f,3.0f,4.0f,5.0f};
    Point4 rMid{1.5f,2.5f,3.5f,4.5f};
    Segment4 seg4(r0, r1);
    printTest("4D Segment length", std::abs(seg4.length() - std::sqrt(4.0f)) <= Euclid::Tolerance().evaluateEpsilon(1));
    printTest("4D Segment midpoint", seg4.midpoint() == rMid);
    printTest("4D Point on segment", seg4.contains(Point4{rMid}));
    printTest("4D Point off segment", !seg4.contains(Point4{2.5f,3.5f,4.5f,5.5f}));

    // --- 5D Tests ---
    using Point5 = Point<double,5>;
    using Segment5 = Segment<double,5>;
    Point5 s0{0.0,0.0,0.0,0.0,0.0};
    Point5 s1{1.0,2.0,3.0,4.0,5.0};
    Point5 sMid{0.5,1.0,1.5,2.0,2.5};
    Segment5 seg5(s0, sMid);
    printTest("5D Segment length", std::abs(seg5.length() - std::sqrt(0.25+1+2.25+4+6.25)) <= Euclid::Tolerance().evaluateEpsilon(1));
    printTest("5D Segment midpoint", seg5.midpoint() == Point5{0.25,0.5,0.75,1.0,1.25});
    printTest("5D Point on segment", seg5.contains(Point5{sMid}));
    printTest("5D Point off segment", !seg5.contains(Point5{1.0,2.0,3.0,4.0,6.0}));

    // --- Degenerate Segment Tests (start == end) ---
    Segment2 deg2(p0, p0);
    printTest("Degenerate 2D segment length zero", deg2.length() == 0.0f);
    printTest("Degenerate 2D segment contains start", deg2.contains(Point2{p0}));
    printTest("Degenerate 2D segment does not contain other", !deg2.contains(Point2{1.0f,0.0f}));

    Segment3 deg3(q0, q0);
    printTest("Degenerate 3D segment length zero", deg3.length() == 0.0);
    printTest("Degenerate 3D segment contains start", deg3.contains(Point3{q0}));
    printTest("Degenerate 3D segment does not contain other", !deg3.contains(Point3{0.0,1.0,0.0}));

    // --- Transform Tests (similar to testLine) ---
    // 2D
    Eigen::Matrix2f scale;
    scale << 2.0f, 0.0f,
             0.0f, 2.0f;
    Segment2 seg2Scaled = seg2;
    seg2Scaled.applyTransform(scale);
    printTest("Segment transform start", seg2Scaled.start == Point2{0.0f,0.0f});
    printTest("Segment transform end", seg2Scaled.end == Point2{2.0f,0.0f});

    // 3D
    Eigen::Matrix3d rot = Eigen::Matrix3d::Identity();
    rot(0,0) = 0; rot(0,1) = -1; rot(1,0) = 1; rot(1,1) = 0; // 90 deg rotation in XY
    Segment3 seg3Rot = seg3;
    seg3Rot.applyTransform(rot);
    printTest("3D Segment transform start", seg3Rot.start == Point3{0.0,0.0,0.0});
    printTest("3D Segment transform end", seg3Rot.end == Point3{-1.0,0.0,0.0});
    
    
    testSegmentIntersection();
}
