#pragma once

#include <iostream>
#include "geometry/point.hpp"
#include "geometry/line.hpp"
#include "geometry/curve.hpp"
#include "test_utilities.hpp"


using namespace euclid::geometry;
using namespace euclid::tests;

inline void testCurve() {
    std::cout << "\n∿ Testing Curve Primitive\n";


    // --- 2D Linear Curve ---
    using Point2 = Point<float,2>;
    using Curve2 = Curve<float,2>;

    Point2 p0{0.0f,0.0f};
    Point2 p1{1.0f,1.0f};
    Curve2 linear = Curve2::linearCurve(p0,p1);

    Point2 eval0 = linear.evaluate(0.0f);
    Point2 eval1 = linear.evaluate(1.0f);
    Point2 evalHalf = linear.evaluate(0.5f);

    printTest("Linear 2D curve t=0", approxEqual(eval0,p0,1e-6f));
    printTest("Linear 2D curve t=1", approxEqual(eval1,p1,1e-6f));
    printTest("Linear 2D curve t=0.5", approxEqual(evalHalf,Point2{0.5f,0.5f},1e-6f));

    // --- 3D Linear Curve ---
    using Point3 = Point<float,3>;
    using Curve3 = Curve<float,3>;
    Point3 p3_0{0.0f,0.0f,0.0f};
    Point3 p3_1{1.0f,2.0f,3.0f};
    Curve3 linear3 = Curve3::linearCurve(p3_0,p3_1);

    Point3 eval3Half = linear3.evaluate(0.5f);
    printTest("Linear 3D curve t=0.5", approxEqual(eval3Half,Point3{0.5f,1.0f,1.5f},1e-6f));

    // --- Nonlinear Curve ---
    auto nonlinearFunc = [](float t) -> Point2 { return Point2{t, t*t}; };
    Curve2 nonlinear(nonlinearFunc,0.0f,1.0f);

    Point2 evalNonlin = nonlinear.evaluate(0.5f);
    printTest("Nonlinear 2D curve t=0.5", approxEqual(evalNonlin,Point2{0.5f,0.25f},1e-6f));

    // --- 3D Nonlinear Curve ---
    auto nonlinear3D = [](float t) -> Point3 { return Point3{t, t*t, t*t*t}; };
    Curve3 nonlinearCurve3D(nonlinear3D, 0.0f, 1.0f);
    Point3 eval3D_0 = nonlinearCurve3D.evaluate(0.0f);
    Point3 eval3D_half = nonlinearCurve3D.evaluate(0.5f);
    Point3 eval3D_1 = nonlinearCurve3D.evaluate(1.0f);
    printTest("Nonlinear 3D curve t=0.0", approxEqual(eval3D_0, Point3{0.0f, 0.0f, 0.0f}, 1e-6f));
    printTest("Nonlinear 3D curve t=0.5", approxEqual(eval3D_half, Point3{0.5f, 0.25f, 0.125f}, 1e-6f));
    printTest("Nonlinear 3D curve t=1.0", approxEqual(eval3D_1, Point3{1.0f, 1.0f, 1.0f}, 1e-6f));

    // --- 2D Sinusoidal Curve: y = sin(2*pi*x), x = t in [0,1] ---
    auto sinusoidalFunc = [](float t) -> Point2 {
        return Point2{t, std::sin(2.0f * float(M_PI) * t)};
    };
    Curve2 sinusoidal(sinusoidalFunc, 0.0f, 1.0f);
    printTest("Sinusoidal 2D curve t=0", approxEqual(sinusoidal.evaluate(0.0f), Point2{0.0f, 0.0f}, 1e-6f));
    printTest("Sinusoidal 2D curve t=0.25", approxEqual(sinusoidal.evaluate(0.25f), Point2{0.25f, 1.0f}, 1e-6f));
    printTest("Sinusoidal 2D curve t=0.5", approxEqual(sinusoidal.evaluate(0.5f), Point2{0.5f, 0.0f}, 1e-6f));
    printTest("Sinusoidal 2D curve t=0.75", approxEqual(sinusoidal.evaluate(0.75f), Point2{0.75f, -1.0f}, 1e-6f));
    printTest("Sinusoidal 2D curve t=1.0", approxEqual(sinusoidal.evaluate(1.0f), Point2{1.0f, 0.0f}, 1e-6f));

    // --- Transform on Curve ---
    Eigen::Matrix2f scaleMat; scaleMat << 2.0f,0.0f,0.0f,3.0f;
    Eigen::Vector2f translation; translation << 1.0f,-1.0f;
    euclid::algebra::Affine<float,2> transform(scaleMat, translation);

    Curve2 transformedCurve = nonlinear.applyTransform(transform);
    Point2 evalTransformed = transformedCurve.evaluate(0.5f);
    Point2 expectedTransformed{2.0f, -0.25f + -1.0f}; // x scaled by 2, y by 3, then translation applied
    printTest("Curve transform evaluation", approxEqual(evalTransformed,Point2{2.0f, -0.25f},1e-5f));

    // --- Multi-dimensional Nonlinear Curve Tests (2D–5D) ---

    // 2D Nonlinear Curve
    auto func2D = [](float t) -> Point<float,2> { return Point<float,2>{-t, t * float(M_PI)}; };
    Curve<float,2> curve2(func2D, 0.0f, 1.0f);
    printTest("Nonlinear 2D t=0", approxEqual(curve2.evaluate(0.0f), Point<float,2>{0.0f,0.0f},1e-6f));
    printTest("Nonlinear 2D t=0.5", approxEqual(curve2.evaluate(0.5f), Point<float,2>{-0.5f, 0.5f*float(M_PI)},1e-6f));
    printTest("Nonlinear 2D t=1.0", approxEqual(curve2.evaluate(1.0f), Point<float,2>{-1.0f, float(M_PI)},1e-6f));
    printTest("Nonlinear 2D t=-0.5", approxEqual(curve2.evaluate(-0.5f), Point<float,2>{0.5f, -0.5f*float(M_PI)},1e-6f));
    printTest("Nonlinear 2D t=1.5", approxEqual(curve2.evaluate(1.5f), Point<float,2>{-1.5f, 1.5f*float(M_PI)},1e-6f));

    // 3D Nonlinear Curve
    auto func3D = [](float t) -> Point<float,3> { return Point<float,3>{-t, t * float(M_PI), t*t}; };
    Curve<float,3> curve3(func3D, 0.0f, 1.0f);
    printTest("Nonlinear 3D t=0.5", approxEqual(curve3.evaluate(0.5f), Point<float,3>{-0.5f, 0.5f*float(M_PI), 0.25f},1e-6f));
    printTest("Nonlinear 3D t=1.5", approxEqual(curve3.evaluate(1.5f), Point<float,3>{-1.5f, 1.5f*float(M_PI), 2.25f},1e-6f));

    // 4D Nonlinear Curve
    auto func4D = [](float t) -> Point<float,4> { return Point<float,4>{-t, t * float(M_PI), t*t, t*t*t}; };
    Curve<float,4> curve4(func4D, 0.0f, 1.0f);
    printTest("Nonlinear 4D t=0.5", approxEqual(curve4.evaluate(0.5f), Point<float,4>{-0.5f, 0.5f*float(M_PI), 0.25f, 0.125f},1e-6f));
    printTest("Nonlinear 4D t=-0.5", approxEqual(curve4.evaluate(-0.5f), Point<float,4>{0.5f, -0.5f*float(M_PI), 0.25f, -0.125f},1e-6f));

    // 5D Nonlinear Curve
    auto func5D = [](float t) -> Point<float,5> { return Point<float,5>{-t, t * float(M_PI), t*t, t*t*t, std::sqrt(t+1.0f)}; };
    Curve<float,5> curve5(func5D, 0.0f, 1.0f);
    printTest("Nonlinear 5D t=0.5", approxEqual(curve5.evaluate(0.5f), Point<float,5>{-0.5f, 0.5f*float(M_PI), 0.25f, 0.125f, std::sqrt(1.5f)},1e-6f));
    printTest("Nonlinear 5D t=1.5", approxEqual(curve5.evaluate(1.5f), Point<float,5>{-1.5f, 1.5f*float(M_PI), 2.25f, 3.375f, std::sqrt(2.5f)},1e-6f));

}
