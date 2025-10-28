#pragma once

#include <iostream>
#include "geometry/point.hpp"
#include "geometry/line.hpp"
#include "geometry/curve.hpp"
#include "test_utilities.hpp"


using namespace Euclid::Geometry;
using namespace Euclid::Algebra;
using namespace Euclid::Tests;

inline void testCurveDerivative() {
    std::cout << "\n∿ Testing Curve Derivative (Adaptive Tolerance)\n";

    Euclid::Tolerance tol;
    // Unified adaptive tolerance lambda
    auto adaptiveTolForCurve = [&](const auto& expectedVec, const auto& posVec) {
        float geomScale = posVec.norm();
        float sumAbs = posVec.cwiseAbs().sum();
        float slopeEstimate = geomScale + sumAbs;
        float baseTol = tol.evaluateEpsilon(std::max(geomScale, expectedVec.norm()));
        return baseTol * (10.0f + slopeEstimate);
    };

    using Point2 = Point<float,2>;
    using Curve2 = Curve<float,2>;
    Point2 p0{0.0f, 0.0f};
    Point2 p1{1.0f, 1.0f};
    Curve2 linear = Curve2::linearCurve(p0, p1);

    float tMid = 0.5f;
    auto d = linear.evaluateDerivative(tMid);
    Eigen::Vector2f expected(1.0f, 1.0f);
    Eigen::Vector2f posExpected(0.5f, 0.5f);
    float derivTol = adaptiveTolForCurve(expected, posExpected);
    bool ok = (d - expected).norm() <= 3.0f * derivTol;
    std::cout << "[DEBUG] expected: " << expected.transpose()
              << ", got: " << d.transpose()
              << ", tol: " << derivTol << std::endl;
    printTest("Linear curve derivative (adaptive)", ok);

    // --- Quadratic Curve Derivative Test ---
    auto quadraticFunc = [](float t) -> Point2 { return Point2{t, t * t}; };
    Curve2 quadratic(quadraticFunc, 0.0f, 1.0f);
    auto dq = quadratic.evaluateDerivative(tMid);
    Eigen::Vector2f expectedQ(1.0f, 2.0f * tMid);
    Eigen::Vector2f posExpectedQ(0.5f, 0.25f);
    float derivTolQ = adaptiveTolForCurve(expectedQ, posExpectedQ);
    bool okQ = (dq - expectedQ).norm() <= 3.0f * derivTolQ;
    std::cout << "[DEBUG] expected: " << expectedQ.transpose()
              << ", got: " << dq.transpose()
              << ", tol: " << derivTolQ << std::endl;
    printTest("Quadratic curve derivative (adaptive)", okQ);

    // --- Exponential Curve Derivative Test ---
    auto exponentialFunc = [](float t) -> Point2 { return Point2{t, std::exp(t)}; };
    Curve2 exponential(exponentialFunc, 0.0f, 1.0f);
    auto de = exponential.evaluateDerivative(tMid);
    Eigen::Vector2f expectedE(1.0f, std::exp(tMid));
    Eigen::Vector2f posExpectedE(0.5f, std::exp(tMid));
    float derivTolE = adaptiveTolForCurve(expectedE, posExpectedE);
    bool okE = (de - expectedE).norm() <= 3.0f * derivTolE;
    std::cout << "[DEBUG] expected: " << expectedE.transpose()
              << ", got: " << de.transpose()
              << ", tol: " << derivTolE << std::endl;
    printTest("Exponential curve derivative (adaptive)", okE);

    // --- Sinusoidal Curve Derivative Test ---
    auto sinusoidalFunc = [](float t) -> Point2 { return Point2{t, std::sin(2.0f * float(M_PI) * t)}; };
    Curve2 sinusoidal(sinusoidalFunc, 0.0f, 1.0f);
    auto ds = sinusoidal.evaluateDerivative(tMid);
    Eigen::Vector2f expectedS(1.0f, 2.0f * float(M_PI) * std::cos(2.0f * float(M_PI) * tMid));
    Eigen::Vector2f posExpectedS(tMid, std::sin(2.0f * float(M_PI) * tMid));
    float derivTolS = adaptiveTolForCurve(expectedS, posExpectedS);
    float absErr = (ds - expectedS).norm();
    float relErr = absErr / std::max(expectedS.norm(), 1e-8f);
    bool okS = absErr <= 5.0f * derivTolS;
    std::cout << "[DEBUG] expected: " << expectedS.transpose()
              << ", got: " << ds.transpose()
              << ", absErr: " << absErr
              << ", relErr: " << relErr
              << ", tol: " << derivTolS << std::endl;
    printTest("Sinusoidal curve derivative (adaptive)", okS);

    // --- Circular Curve Derivative Test ---
    tMid = 0.25f;
    auto circularFunc = [](float t) -> Point2 {
        return Point2{std::cos(2.0f * float(M_PI) * t), std::sin(2.0f * float(M_PI) * t)};
    };
    Curve2 circular(circularFunc, 0.0f, 1.0f);
    auto dc = circular.evaluateDerivative(tMid);
    Eigen::Vector2f expectedC(-2.0f * float(M_PI) * std::sin(2.0f * float(M_PI) * tMid),
                              2.0f * float(M_PI) * std::cos(2.0f * float(M_PI) * tMid));
    Eigen::Vector2f posExpectedC(std::cos(2.0f * float(M_PI) * tMid), std::sin(2.0f * float(M_PI) * tMid));
    float derivTolC = adaptiveTolForCurve(expectedC, posExpectedC);
    absErr = (dc - expectedC).norm();
    relErr = absErr / std::max(expectedC.norm(), 1e-8f);
    bool okC = absErr <= 5.0f * derivTolC;
    std::cout << "[DEBUG] expected: " << expectedC.transpose()
              << ", got: " << dc.transpose()
              << ", absErr: " << absErr
              << ", relErr: " << relErr
              << ", tol: " << derivTolC << std::endl;
    printTest("Circular curve derivative (adaptive)", okC);

    // --- Cubic Curve Derivative Test ---
    auto cubicFunc = [](float t) -> Point2 { return Point2{t, t*t*t - 3*t*t + 2*t}; };
    Curve2 cubic(cubicFunc, 0.0f, 1.0f);
    auto dcub = cubic.evaluateDerivative(0.5f);
    float tCub = 0.5f;
    Eigen::Vector2f expectedCub(1.0f, 3*tCub*tCub - 6*tCub + 2);
    Eigen::Vector2f posExpectedCub(tCub, tCub*tCub*tCub - 3*tCub*tCub + 2*tCub);
    float derivTolCub = adaptiveTolForCurve(expectedCub, posExpectedCub);
    bool okCub = (dcub - expectedCub).norm() <= 3.0f * derivTolCub;
    std::cout << "[DEBUG] expected: " << expectedCub.transpose()
              << ", got: " << dcub.transpose()
              << ", tol: " << derivTolCub << std::endl;
    printTest("Cubic curve derivative (adaptive)", okCub);

    // --- Rational Quadratic Derivative Test ---
    auto rationalFunc = [](float t) -> Point2 {
        float denom = 1.0f + t*t;
        return Point2{t / denom, t*t / denom};
    };
    Curve2 rational(rationalFunc, -1.0f, 1.0f);
    float tR = 0.5f;
    auto dr = rational.evaluateDerivative(tR);
    float denom = 1.0f + tR*tR;
    Eigen::Vector2f expectedR(
        (1.0f*denom - tR*(2*tR)) / (denom*denom),
        (2*tR*denom - tR*tR*(2*tR)) / (denom*denom)
    );
    Eigen::Vector2f posExpectedR(tR / denom, tR*tR / denom);
    float derivTolR = adaptiveTolForCurve(expectedR, posExpectedR);
    bool okR = (dr - expectedR).norm() <= 3.0f * derivTolR;
    std::cout << "[DEBUG] expected: " << expectedR.transpose()
              << ", got: " << dr.transpose()
              << ", tol: " << derivTolR << std::endl;
    printTest("Rational quadratic curve derivative (adaptive)", okR);

    // --- 3D Helix Derivative Test ---
    using Point3 = Point<float,3>;
    using Curve3 = Curve<float,3>;
    auto helixFunc = [](float t) -> Point3 {
        return Point3{static_cast<float>(std::cos(2*M_PI*t)), static_cast<float>(std::sin(2*M_PI*t)), t};
    };
    Curve3 helix(helixFunc, 0.0f, 1.0f);
    float tHelix = 0.25f;
    auto dh = helix.evaluateDerivative(tHelix);
    Eigen::Vector3f expectedH(-2*M_PI*std::sin(2*M_PI*tHelix), 2*M_PI*std::cos(2*M_PI*tHelix), 1.0f);
    Eigen::Vector3f posExpectedH(std::cos(2*M_PI*tHelix), std::sin(2*M_PI*tHelix), tHelix);
    float derivTolH = adaptiveTolForCurve(expectedH, posExpectedH);
    bool okH = (dh - expectedH).norm() <= 3.0f * derivTolH;
    std::cout << "[DEBUG] expected: " << expectedH.transpose()
              << ", got: " << dh.transpose()
              << ", tol: " << derivTolH << std::endl;
    printTest("Helix curve derivative (adaptive)", okH);

    // --- Sigmoid Curve Derivative Test ---
    auto sigmoidFunc = [](float t) -> Point2 {
        return Point2{t, 1.0f / (1.0f + std::exp(-10.0f * (t - 0.5f)))};
    };
    Curve2 sigmoid(sigmoidFunc, 0.0f, 1.0f);
    float tSig = 0.5f;
    auto dsig = sigmoid.evaluateDerivative(tSig);
    float expTerm = std::exp(-10.0f * (tSig - 0.5f));
    float denom2 = std::pow(1.0f + expTerm, 2);
    Eigen::Vector2f expectedSig(1.0f, (10.0f * expTerm) / denom2);
    Eigen::Vector2f posExpectedSig(tSig, 1.0f / (1.0f + std::exp(-10.0f * (tSig - 0.5f))));
    float derivTolSig = adaptiveTolForCurve(expectedSig, posExpectedSig);
    bool okSig = (dsig - expectedSig).norm() <= 3.0f * derivTolSig;
    std::cout << "[DEBUG] expected: " << expectedSig.transpose()
              << ", got: " << dsig.transpose()
              << ", tol: " << derivTolSig << std::endl;
    printTest("Sigmoid curve derivative (adaptive)", okSig);

    // --- Torus Knot Derivative Test ---
    auto torusKnotFunc = [](float t) -> Point3 {
        const float R = 2.0f, r = 0.5f;
        const int p = 2, q = 3;
        float angleP = p * 2 * M_PI * t;
        float angleQ = q * 2 * M_PI * t;
        return Point3{(R + r * std::cos(angleQ)) * std::cos(angleP),
                      (R + r * std::cos(angleQ)) * std::sin(angleP),
                      r * std::sin(angleQ)};
    };
    Curve3 torusKnot(torusKnotFunc, 0.0f, 1.0f);
    float tTK = 0.25f, R = 2.0f, r = 0.5f, p = 2, q = 3;
    float angleP = p * 2*M_PI*tTK, angleQ = q * 2*M_PI*tTK;
    auto dtk = torusKnot.evaluateDerivative(tTK);
    Eigen::Vector3f expectedTK(
        -p*2*M_PI*(R + r*std::cos(angleQ))*std::sin(angleP) - r*q*2*M_PI*std::sin(angleQ)*std::cos(angleP),
         p*2*M_PI*(R + r*std::cos(angleQ))*std::cos(angleP) - r*q*2*M_PI*std::sin(angleQ)*std::sin(angleP),
         r*q*2*M_PI*std::cos(angleQ)
    );
    Eigen::Vector3f posExpectedTK((R + r*std::cos(angleQ))*std::cos(angleP),
                                  (R + r*std::cos(angleQ))*std::sin(angleP),
                                  r*std::sin(angleQ));
    absErr = (dtk - expectedTK).norm();
    relErr = absErr / std::max(expectedTK.norm(), 1e-8f);
    float geomScaleTK = posExpectedTK.norm();
    float slopeEstimateTK = expectedTK.norm() + geomScaleTK;

    // Compute base tolerance from the unified tolerance model
    float baseTolTK = tol.evaluateEpsilon(std::max(geomScaleTK, expectedTK.norm()));

    // Scale tolerance dynamically by curvature and slope (robust for high-curvature)
    float scaledTolTK = baseTolTK * (200.0f + 5.0f * slopeEstimateTK);
    float relTolTK    = baseTolTK * (250.0f + 6.0f * slopeEstimateTK);

    bool okTK = (absErr <= scaledTolTK) || (relErr <= relTolTK);

    std::cout << "[DEBUG] expected: " << expectedTK.transpose()
              << ", got: " << dtk.transpose()
              << ", absErr: " << absErr
              << ", relErr: " << relErr
              << ", scaledTol: " << scaledTolTK
              << ", relTol: " << relTolTK
              << std::endl;

    printTest("Torus knot curve derivative (adaptive)", okTK);
}

// Test evaluateIntegral() for representative curves
inline void testCurveIntegral() {
    std::cout << "\n∿ Testing Curve Integral (Adaptive Tolerance)\n";

    Euclid::Tolerance tol;
    using Point2 = Point<float,2>;
    using Curve2 = Curve<float,2>;
    using Point3 = Point<float,3>;
    using Curve3 = Curve<float,3>;

    // --- Linear Curve: from (0,0) to (1,1) ---
    Point2 p0{0.0f, 0.0f};
    Point2 p1{1.0f, 1.0f};
    Curve2 linear = Curve2::linearCurve(p0, p1);
    float expectedLinear = std::sqrt(2.0f); // length of diagonal in unit square
    float computedLinear = linear.evaluateIntegral();
    float diffLinear = std::abs(computedLinear - expectedLinear);
    float tolLinear = tol.evaluateEpsilon(expectedLinear) * 10.0f;
    std::cout << "[DEBUG] Linear curve: expected=" << expectedLinear
              << ", computed=" << computedLinear
              << ", diff=" << diffLinear
              << ", tol=" << tolLinear << std::endl;
    printTest("Linear curve integral (arc length)", diffLinear <= 3.0f * tolLinear);

    // --- Quadratic Curve: (t, t^2), t in [0,1] ---
    auto quadraticFunc = [](float t) -> Point2 { return Point2{t, t*t}; };
    Curve2 quadratic(quadraticFunc, 0.0f, 1.0f);
    // Analytical: ∫₀¹ sqrt(1 + 4t^2) dt
    // This can be computed: (sinh⁻¹(2t)/4) + (t/2) * sqrt(1 + 4t^2) |₀¹
    auto F = [](float t) {
        return 0.25f * std::asinh(2.0f * t) + 0.5f * t * std::sqrt(1.0f + 4.0f * t * t);
    };
    float expectedQuad = F(1.0f) - F(0.0f);
    float computedQuad = quadratic.evaluateIntegral();
    float diffQuad = std::abs(computedQuad - expectedQuad);
    float tolQuad = tol.evaluateEpsilon(expectedQuad) * 15.0f;
    std::cout << "[DEBUG] Quadratic curve: expected=" << expectedQuad
              << ", computed=" << computedQuad
              << ", diff=" << diffQuad
              << ", tol=" << tolQuad << std::endl;
    printTest("Quadratic curve integral (arc length)", diffQuad <= 5.0f * tolQuad);

    // --- Circular Curve: unit circle x=cos(2πt), y=sin(2πt), t in [0,1] ---
    auto circularFunc = [](float t) -> Point2 {
        return Point2{std::cos(2.0f * float(M_PI) * t), std::sin(2.0f * float(M_PI) * t)};
    };
    Curve2 circular(circularFunc, 0.0f, 1.0f);
    float expectedCircle = 2.0f * float(M_PI); // circumference of unit circle
    float computedCircle = circular.evaluateIntegral();
    float diffCircle = std::abs(computedCircle - expectedCircle);
    float tolCircle = tol.evaluateEpsilon(expectedCircle) * 15.0f;
    std::cout << "[DEBUG] Circular curve: expected=" << expectedCircle
              << ", computed=" << computedCircle
              << ", diff=" << diffCircle
              << ", tol=" << tolCircle << std::endl;
    printTest("Circular curve integral (arc length)", diffCircle <= 5.0f * tolCircle);

    // --- Sinusoidal Curve: (t, sin(2πt)), t in [0,1] ---
    auto sinusoidalFunc = [](float t) -> Point2 { return Point2{t, std::sin(2.0f * float(M_PI) * t)}; };
    Curve2 sinusoidal(sinusoidalFunc, 0.0f, 1.0f);
    // Corrected analytical reference:
    // ∫₀¹ sqrt(1 + (2π cos(2πt))^2) dt ≈ 4.18828
    // Verified numerically with WolframAlpha and matches Euclid output
    float expectedSine = 4.18828f;
    float computedSine = sinusoidal.evaluateIntegral();
    float diffSine = std::abs(computedSine - expectedSine);
    float tolSine = tol.evaluateEpsilon(expectedSine) * 30.0f;
    std::cout << "[DEBUG] Sinusoidal curve: expected=" << expectedSine
              << ", computed=" << computedSine
              << ", diff=" << diffSine
              << ", tol=" << tolSine << std::endl;
    printTest("Sinusoidal curve integral (arc length)", diffSine <= 5.0f * tolSine);

    // --- 3D Helix: x=cos(2πt), y=sin(2πt), z=t, t in [0,1] ---
    auto helixFunc = [](float t) -> Point3 {
        return Point3{static_cast<float>(std::cos(2*M_PI*t)), static_cast<float>(std::sin(2*M_PI*t)), t};
    };
    Curve3 helix(helixFunc, 0.0f, 1.0f);
    // Arc length for one turn: ∫₀¹ sqrt((−2π sin(2πt))^2 + (2π cos(2πt))^2 + 1) dt
    // The sum of squares of the first two terms is (2π)^2, so integrand is sqrt((2π)^2 + 1)
    float expectedHelix = std::sqrt(4.0f * float(M_PI) * float(M_PI) + 1.0f); // constant integrand
    float computedHelix = helix.evaluateIntegral();
    float diffHelix = std::abs(computedHelix - expectedHelix);
    float tolHelix = tol.evaluateEpsilon(expectedHelix) * 20.0f;
    std::cout << "[DEBUG] Helix curve: expected=" << expectedHelix
              << ", computed=" << computedHelix
              << ", diff=" << diffHelix
              << ", tol=" << tolHelix << std::endl;
    printTest("Helix curve integral (arc length)", diffHelix <= 5.0f * tolHelix);

    // --- Exponential Curve: (t, e^t), t in [0,1] ---
    auto exponentialFunc = [](float t) -> Point2 { return Point2{t, std::exp(t)}; };
    Curve2 exponential(exponentialFunc, 0.0f, 1.0f);
    // Arc length: ∫₀¹ sqrt(1 + (e^t)^2) dt ≈ 2.00345
    float expectedExp = 2.00345f;
    float computedExp = exponential.evaluateIntegral();
    float diffExp = std::abs(computedExp - expectedExp);
    float tolExp = tol.evaluateEpsilon(expectedExp) * 20.0f;
    std::cout << "[DEBUG] Exponential curve: expected=" << expectedExp
              << ", computed=" << computedExp
              << ", diff=" << diffExp
              << ", tol=" << tolExp << std::endl;
    printTest("Exponential curve integral (arc length)", diffExp <= 5.0f * tolExp);

    // --- Cubic Curve: (t, t^3 - 3t^2 + 2t), t in [0,1] ---
    auto cubicFunc = [](float t) -> Point2 { return Point2{t, t*t*t - 3*t*t + 2*t}; };
    Curve2 cubic(cubicFunc, 0.0f, 1.0f);
    // Approximate reference from numerical integration ≈ 1.31136
    float expectedCubic = 1.31136f;
    float computedCubic = cubic.evaluateIntegral();
    float diffCubic = std::abs(computedCubic - expectedCubic);
    float tolCubic = tol.evaluateEpsilon(expectedCubic) * 25.0f;
    std::cout << "[DEBUG] Cubic curve: expected=" << expectedCubic
              << ", computed=" << computedCubic
              << ", diff=" << diffCubic
              << ", tol=" << tolCubic << std::endl;
    printTest("Cubic curve integral (arc length)", diffCubic <= 5.0f * tolCubic);

    // --- Rational Quadratic Curve: (t/(1+t^2), t^2/(1+t^2)), t in [-1,1] ---
    auto rationalFunc = [](float t) -> Point2 {
        float denom = 1.0f + t*t;
        return Point2{t / denom, t*t / denom};
    };
    Curve2 rational(rationalFunc, -1.0f, 1.0f);
    // Approximate analytical length ≈ pi/2
    float expectedRational = M_PI / 2.0f;
    float computedRational = rational.evaluateIntegral();
    float diffRational = std::abs(computedRational - expectedRational);
    float tolRational = tol.evaluateEpsilon(expectedRational) * 25.0f;
    std::cout << "[DEBUG] Rational quadratic curve: expected=" << expectedRational
              << ", computed=" << computedRational
              << ", diff=" << diffRational
              << ", tol=" << tolRational << std::endl;
    printTest("Rational quadratic curve integral (arc length)", diffRational <= 5.0f * tolRational);

    // --- Sigmoid Curve: (t, 1/(1+e^{-10(t-0.5)})), t in [0,1] ---
    auto sigmoidFunc = [](float t) -> Point2 {
        return Point2{t, 1.0f / (1.0f + std::exp(-10.0f * (t - 0.5f)))};
    };
    Curve2 sigmoid(sigmoidFunc, 0.0f, 1.0f);
    // Numerical reference ≈ 1.52326
    float expectedSigmoid = 1.52326f;
    float computedSigmoid = sigmoid.evaluateIntegral();
    float diffSigmoid = std::abs(computedSigmoid - expectedSigmoid);
    float tolSigmoid = tol.evaluateEpsilon(expectedSigmoid) * 25.0f;
    std::cout << "[DEBUG] Sigmoid curve: expected=" << expectedSigmoid
              << ", computed=" << computedSigmoid
              << ", diff=" << diffSigmoid
              << ", tol=" << tolSigmoid << std::endl;
    printTest("Sigmoid curve integral (arc length)", diffSigmoid <= 5.0f * tolSigmoid);

    // --- 3D Torus Knot (p=2, q=3): R=2, r=0.5, t in [0,1] ---
    auto torusKnotFunc = [](float t) -> Point3 {
        const float R = 2.0f, r = 0.5f;
        const int p = 2, q = 3;
        float angleP = p * 2 * M_PI * t;
        float angleQ = q * 2 * M_PI * t;
        return Point3{(R + r * std::cos(angleQ)) * std::cos(angleP),
                      (R + r * std::cos(angleQ)) * std::sin(angleP),
                      r * std::sin(angleQ)};
    };
    Curve3 torusKnot(torusKnotFunc, 0.0f, 1.0f);
    // Reference from numerical integration ≈ 26.8887
    float expectedTorus = 26.8887f;
    float computedTorus = torusKnot.evaluateIntegral();
    float diffTorus = std::abs(computedTorus - expectedTorus);
    float tolTorus = tol.evaluateEpsilon(expectedTorus) * 50.0f;
    std::cout << "[DEBUG] Torus knot curve: expected=" << expectedTorus
              << ", computed=" << computedTorus
              << ", diff=" << diffTorus
              << ", tol=" << tolTorus << std::endl;
    printTest("Torus knot curve integral (arc length)", diffTorus <= 5.0f * tolTorus);
}

// Test evaluateCurvature() for representative curves
inline void testCurveCurvature() {
    std::cout << "\n∿ Testing Curve Local Curvature (Adaptive Tolerance)\n";

    Euclid::Tolerance tol;
    using Point2 = Euclid::Geometry::Point<float, 2>;
    using Curve2 = Euclid::Geometry::Curve<float, 2>;
    using Point3 = Euclid::Geometry::Point<float, 3>;
    using Curve3 = Euclid::Geometry::Curve<float, 3>;

    // --- Linear Curve: (t, t), t in [0,1] ---
    Curve2 linear = Curve2::linearCurve(Point2{0.0f, 0.0f}, Point2{1.0f, 1.0f});
    float tLin = 0.5f;
    float expectedLin = 0.0f;
    float computedLin = linear.evaluateCurvature(tLin);
    float tolLin = tol.evaluateEpsilon(1.0f) * 10.0f;
    float scaledTolLin = tolLin * 100.0f;
    float diffLin = std::abs(computedLin - expectedLin);
    std::cout << "[DEBUG] Linear: expected=" << expectedLin
              << ", computed=" << computedLin
              << ", diff=" << diffLin
              << ", tol=" << std::scientific << std::setprecision(6) << tolLin
              << ", scaledTol=" << std::scientific << std::setprecision(6) << scaledTolLin
              << std::defaultfloat << std::endl;
    printTest("Linear curve curvature", diffLin <= scaledTolLin);

    // --- Circular Curve: unit circle x=cos(2πt), y=sin(2πt), t in [0,1] ---
    float tCirc = 0.25f;
    auto circularFunc = [](float t) -> Point2 {
        return Point2{std::cos(2.0f * float(M_PI) * t), std::sin(2.0f * float(M_PI) * t)};
    };
    Curve2 circular(circularFunc, 0.0f, 1.0f);
    // For a unit circle parametrized by t in [0,1], the curvature is always 1
    float expectedCirc = 1.0f;
    float computedCirc = circular.evaluateCurvature(tCirc);
    float tolCirc = tol.evaluateEpsilon(expectedCirc) * 10.0f;
    float scaledTolCirc = tolCirc * 100.0f;
    float diffCirc = std::abs(computedCirc - expectedCirc);
    std::cout << "[DEBUG] Circular: expected=" << expectedCirc
              << ", computed=" << computedCirc
              << ", diff=" << diffCirc
              << ", tol=" << std::scientific << std::setprecision(6) << tolCirc
              << ", scaledTol=" << std::scientific << std::setprecision(6) << scaledTolCirc
              << std::defaultfloat << std::endl;
    printTest("Circular curve curvature", diffCirc <= scaledTolCirc);

    // --- Quadratic Curve: (t, t^2), t in [0,1] ---
    float tQuad = 0.5f;
    auto quadraticFunc = [](float t) -> Point2 { return Point2{t, t * t}; };
    Curve2 quadratic(quadraticFunc, 0.0f, 1.0f);
    // Analytical curvature: kappa = |x'y'' - y'x''| / (x'^2 + y'^2)^(3/2)
    // x = t, y = t^2
    // x' = 1, x'' = 0, y' = 2t, y'' = 2
    // Numerator: |1*2 - 2t*0| = 2
    // Denominator: (1^2 + (2t)^2)^(3/2) = (1 + 4t^2)^(3/2)
    float denomQuad = std::pow(1.0f + 4.0f * tQuad * tQuad, 1.5f);
    float expectedQuad = 2.0f / denomQuad;
    float computedQuad = quadratic.evaluateCurvature(tQuad);
    float tolQuad = tol.evaluateEpsilon(std::abs(expectedQuad)) * 10.0f;
    float scaledTolQuad = tolQuad * 100.0f;
    float diffQuad = std::abs(computedQuad - expectedQuad);
    std::cout << "[DEBUG] Quadratic: expected=" << expectedQuad
              << ", computed=" << computedQuad
              << ", diff=" << diffQuad
              << ", tol=" << std::scientific << std::setprecision(6) << tolQuad
              << ", scaledTol=" << std::scientific << std::setprecision(6) << scaledTolQuad
              << std::defaultfloat << std::endl;
    printTest("Quadratic curve curvature", diffQuad <= scaledTolQuad);

    // --- Helix Curve: x=cos(2πt), y=sin(2πt), z=t, t in [0,1] ---
    float tHelix = 0.25f;
    auto helixFunc = [](float t) -> Point3 {
        return Point3{static_cast<float>(std::cos(2*M_PI*t)), static_cast<float>(std::sin(2*M_PI*t)), t};
    };
    Curve3 helix(helixFunc, 0.0f, 1.0f);
    float expectedHelix = 0.975295f;
    float computedHelix = helix.evaluateCurvature(tHelix);
    float tolHelix = tol.evaluateEpsilon(std::abs(expectedHelix)) * 10.0f;
    float scaledTolHelix = tolHelix * 200.0f;
    float diffHelix = std::abs(computedHelix - expectedHelix);
    std::cout << "[DEBUG] Helix: expected=" << expectedHelix
              << ", computed=" << computedHelix
              << ", diff=" << diffHelix
              << ", tol=" << std::scientific << std::setprecision(6) << tolHelix
              << ", scaledTol=" << std::scientific << std::setprecision(6) << scaledTolHelix
              << std::defaultfloat << std::endl;
    printTest("Helix curve curvature", diffHelix <= scaledTolHelix);

    // --- Sinusoidal Curve: (t, sin(2πt)), t in [0,1] ---
    float tSine = 0.25f;
    auto sinusoidalFunc = [](float t) -> Point2 { return Point2{t, std::sin(2.0f * float(M_PI) * t)}; };
    Curve2 sinusoidal(sinusoidalFunc, 0.0f, 1.0f);
    // x = t, y = sin(2πt)
    // x' = 1, x'' = 0
    // y' = 2π cos(2πt), y'' = -4π^2 sin(2πt)
    // kappa = |x'y'' - y'x''| / (x'^2 + y'^2)^(3/2)
    float x1 = 1.0f;
    float x2 = 0.0f;
    float y1 = 2.0f * float(M_PI) * std::cos(2.0f * float(M_PI) * tSine);
    float y2 = -4.0f * float(M_PI) * float(M_PI) * std::sin(2.0f * float(M_PI) * tSine);
    float numSine = std::abs(x1 * y2 - y1 * x2);
    float denomSine = std::pow(x1 * x1 + y1 * y1, 1.5f);
    float expectedSine = denomSine > 1e-12f ? numSine / denomSine : 0.0f;
    float computedSine = sinusoidal.evaluateCurvature(tSine);
    float tolSine = tol.evaluateEpsilon(std::abs(expectedSine)) * 10.0f;
    float scaledTolSine = tolSine * 500.0f;
    float diffSine = std::abs(computedSine - expectedSine);
    std::cout << "[DEBUG] Sinusoidal: expected=" << expectedSine
              << ", computed=" << computedSine
              << ", diff=" << diffSine
              << ", tol=" << std::scientific << std::setprecision(6) << tolSine
              << ", scaledTol=" << std::scientific << std::setprecision(6) << scaledTolSine
              << std::defaultfloat << std::endl;
    printTest("Sinusoidal curve curvature", diffSine <= scaledTolSine);
}


// Resolution sweep for curve integral (adaptive tolerance)
inline void testCurveIntegralResolutionSweep() {
    std::cout << "\n∿ Curve Integral Resolution Sweep (Adaptive Tolerance)\n";

    Euclid::Tolerance tol;
    using Point2 = Euclid::Geometry::Point<float, 2>;
    using Curve2 = Euclid::Geometry::Curve<float, 2>;
    //using Point3 = Euclid::Geometry::Point<float, 3>;
    //using Curve3 = Euclid::Geometry::Curve<float, 3>;

    // Representative test curves with verified arc lengths
    std::vector<std::tuple<std::string, std::function<Point2(float)>, float>> tests = {
        {"Linear",      [](float t){ return Point2{t, t}; }, std::sqrt(2.0f)},
        {"Quadratic",   [](float t){ return Point2{t, t*t}; }, 1.47894f},
        {"Exponential", [](float t){ return Point2{t, std::exp(t)}; }, 2.00349f},
        {"Sinusoidal",  [](float t){ return Point2{t, static_cast<float>(std::sin(2.0 * M_PI * t))}; }, 4.18828f},
        {"Cubic",       [](float t){ return Point2{t, t*t*t - 3*t*t + 2*t}; }, 1.31135f},
        {"Sigmoid",     [](float t){ return Point2{t, 1.0f/(1.0f + std::exp(-10*(t - 0.5f)))}; }, 1.52326f}
    };

    std::vector<float> tolScales = {0.1f, 0.3f, 1.0f, 3.0f, 10.0f};

    for (auto& [name, func, expected] : tests) {
        Curve2 c(func, 0.0f, 1.0f);
        std::cout << "→ " << name << std::endl;

        for (float scale : tolScales) {
            float scaledTol = tol.evaluateEpsilon(expected) * scale * 10.0f;
            float computed = c.evaluateIntegral();
            float absErr = std::abs(computed - expected);
            float relErr = absErr / std::max(expected, 1e-8f);

            std::cout << "   scale=" << scale
                      << "  absErr=" << absErr
                      << "  relErr=" << relErr
                      << "  tol=" << scaledTol
                      << std::endl;
        }
    }
}


inline void testCurveDerivativeResolution() {
    std::cout << "\n∿ Testing Curve Derivative Resolution (Confidence-Aware)\n";

    using Point2 = Euclid::Geometry::Point<float, 2>;
    using Curve2 = Euclid::Geometry::Curve<float, 2>;

    Euclid::Tolerance tol;
    float t = 0.5f;

    std::vector<std::pair<std::string, std::function<Point2(float)>>> funcs = {
        {"Linear",    [](float t){ return Point2{t, t}; }},
        {"Quadratic", [](float t){ return Point2{t, t*t}; }},
        {"Cubic",     [](float t){ return Point2{t, t*t*t}; }},
        {"Exponential",[](float t){ return Point2{t, std::exp(3.0f * t)}; }},
        {"Sinusoidal",[](float t){ return Point2{t, std::sin(10.0f * float(M_PI) * t)}; }},
        {"Torus-like", [](float t){
            return Point2{
                static_cast<float>(std::cos(8.0f*M_PI*t) * (2.0f + 0.5f*std::cos(12.0f*M_PI*t))),
                static_cast<float>(std::sin(8.0f*M_PI*t) * (2.0f + 0.5f*std::cos(12.0f*M_PI*t)))
            };
        }},
    };

    for (auto& [name, func] : funcs) {
        Curve2 c(func, 0.0f, 1.0f);
        auto [d, conf] = c.evaluateDerivativeWithConfidence(t);

        // Analytical derivative for comparison
        Eigen::Vector2f expected;
        if (name == "Linear") expected = {1, 1};
        else if (name == "Quadratic") expected = {1, 2*t};
        else if (name == "Cubic") expected = {1, 3*t*t};
        else if (name == "Exponential") expected = {1, 3.0f*std::exp(3.0f*t)};
        else if (name == "Sinusoidal") expected = {1, 10.0f*M_PI*std::cos(10.0f*M_PI*t)};
        else if (name == "Torus-like") {
            float θ = 8.0f*M_PI*t, φ = 12.0f*M_PI*t;
            expected = {
                -std::sin(θ)*(2.0f + 0.5f*std::cos(φ))*8.0f*M_PI
                + std::cos(θ)*(-0.5f*std::sin(φ))*12.0f*M_PI,
                 std::cos(θ)*(2.0f + 0.5f*std::cos(φ))*8.0f*M_PI
                + std::sin(θ)*(-0.5f*std::sin(φ))*12.0f*M_PI
            };
        }

        float absErr = (d - expected).norm();
        float relErr = absErr / std::max(expected.norm(), 1e-8f);

        // Compute adaptive tol from geometry
        Point2 pos = c.evaluate(t);
        float geomScale = pos.coords.norm();
        float baseTol = tol.evaluateEpsilon(std::max(geomScale, expected.norm()));
        float adaptiveTol = baseTol * (10.0f + geomScale + expected.norm());

        std::cout << "→ " << name
                  << " | absErr=" << absErr
                  << "  relErr=" << relErr
                  << "  conf=" << conf
                  << "  tol=" << adaptiveTol
                  << std::endl;
    }
}


// Sweep tolerance scaling for confidence-aware derivative estimation
inline void testCurveDerivativeResolutionSweep() {
    std::cout << "\n∿ Curve Derivative Tolerance Scaling Sweep (Confidence-Aware)\n";

    using Point2 = Euclid::Geometry::Point<float, 2>;
    using Curve2 = Euclid::Geometry::Curve<float, 2>;
    Euclid::Tolerance tol;
    float t = 0.5f;

    // Representative curves
    std::vector<std::pair<std::string, std::function<Point2(float)>>> funcs = {
        {"Linear",    [](float t){ return Point2{t, t}; }},
        {"Quadratic", [](float t){ return Point2{t, t*t}; }},
        {"Sinusoidal",[](float t){ return Point2{t, std::sin(10.0f * float(M_PI) * t)}; }},
        {"Torus-like", [](float t){
            return Point2{
                static_cast<float>(std::cos(8.0f*M_PI*t) * (2.0f + 0.5f*std::cos(12.0f*M_PI*t))),
                static_cast<float>(std::sin(8.0f*M_PI*t) * (2.0f + 0.5f*std::cos(12.0f*M_PI*t)))
            };
        }},
    };

    std::vector<float> scaleFactors = {0.1f, 0.3f, 1.0f, 3.0f, 10.0f};

    for (auto& [name, func] : funcs) {
        Curve2 c(func, 0.0f, 1.0f);

        // Analytical derivative for comparison
        Eigen::Vector2f expected;
        if (name == "Linear") expected = {1, 1};
        else if (name == "Quadratic") expected = {1, 2*t};
        else if (name == "Sinusoidal") expected = {1, 10.0f*M_PI*std::cos(10.0f*M_PI*t)};
        else if (name == "Torus-like") {
            float θ = 8.0f*M_PI*t, φ = 12.0f*M_PI*t;
            expected = {
                -std::sin(θ)*(2.0f + 0.5f*std::cos(φ))*8.0f*M_PI
                + std::cos(θ)*(-0.5f*std::sin(φ))*12.0f*M_PI,
                 std::cos(θ)*(2.0f + 0.5f*std::cos(φ))*8.0f*M_PI
                + std::sin(θ)*(-0.5f*std::sin(φ))*12.0f*M_PI
            };
        }

        Point2 pos = c.evaluate(t);
        float geomScale = pos.coords.norm();
        float baseTol = tol.evaluateEpsilon(std::max(geomScale, expected.norm()));

        for (float scale : scaleFactors) {
            float scaledTol = baseTol * (10.0f + geomScale + expected.norm()) * scale;
            // Here, we assume evaluateDerivativeWithConfidence uses adaptive tolerance internally,
            // but for this sweep, we only change the tolerance used for test/inspection.
            auto [d, conf] = c.evaluateDerivativeWithConfidence(t);
            float absErr = (d - expected).norm();
            float relErr = absErr / std::max(expected.norm(), 1e-8f);
            std::cout << "→ " << name
                      << "  scale=" << scale
                      << "  absErr=" << absErr
                      << "  relErr=" << relErr
                      << "  conf=" << conf
                      << "  scaledTol=" << scaledTol
                      << std::endl;
        }
    }
}

// Correlation between curvature and confidence for representative curves
inline void testCurveDerivativeCurvatureCorrelation() {
    std::cout << "\n∿ Curve Derivative Curvature–Confidence Correlation\n";

    using Point2 = Euclid::Geometry::Point<float, 2>;
    using Curve2 = Euclid::Geometry::Curve<float, 2>;
    float eps = 1e-4f;
    constexpr float kappaFloor = 1e-4f; // clamp to avoid division by ~0 in efficiency metrics
    std::vector<std::pair<std::string, std::function<Point2(float)>>> funcs = {
        {"Linear",    [](float t){ return Point2{t, t}; }},
        {"Quadratic", [](float t){ return Point2{t, t*t}; }},
        {"Sinusoidal",[](float t){ return Point2{t, std::sin(2.0f * float(M_PI) * t)}; }},
        {"Torus-like", [](float t){
            return Point2{
                static_cast<float>(std::cos(8.0f*M_PI*t) * (2.0f + 0.5f*std::cos(12.0f*M_PI*t))),
                static_cast<float>(std::sin(8.0f*M_PI*t) * (2.0f + 0.5f*std::cos(12.0f*M_PI*t)))
            };
        }},
    };
    std::vector<float> tSamples = {0.1f, 0.3f, 0.5f, 0.7f, 0.9f};
    // Vector to store (name, meanEfficiency)
    std::vector<std::pair<std::string, float>> meanEfficiencies;
    std::vector<std::pair<std::string, float>> meanConfVec;
    for (auto& [name, func] : funcs) {
        Curve2 c(func, 0.0f, 1.0f);
        float effSum = 0.0f;
        int effCount = 0;
        for (float t : tSamples) {
            // First derivative (x', y')
            auto [d, conf] = c.evaluateDerivativeWithConfidence(t);
            // Second derivative (x'', y'') via finite difference
            Eigen::Vector2f d_left = c.evaluateDerivative(std::max(0.0f, t - eps));
            Eigen::Vector2f d_right = c.evaluateDerivative(std::min(1.0f, t + eps));
            Eigen::Vector2f d2 = (d_right - d_left) / (2.0f * eps);
            float xp = d[0], yp = d[1];
            float xpp = d2[0], ypp = d2[1];
            float denom = std::pow(xp * xp + yp * yp, 1.5f);
            float curvature = denom > 1e-12f ? std::abs(xp * ypp - yp * xpp) / denom : 0.0f;
            float stability = (curvature > 1e-6f) ? conf / curvature : 0.0f;
            std::cout << "→ " << name
                      << "  t=" << t
                      << "  curvature=" << curvature
                      << "  conf=" << conf
                      << "  stability=" << stability
                      << std::endl;
            // Efficiency: conf / max(curvature, kappaFloor)
            effSum  += (conf / std::max(curvature, kappaFloor));
            effCount += 1;
            // Track confidence for normalization
            // (we'll divide mean efficiency by mean confidence to get the normalized index)
        }
        float meanEfficiency = (effCount > 0) ? (effSum / effCount) : 0.0f;
        // Compute mean confidence across samples for this curve
        // We need the per-sample confidences; recompute lightweight mean here
        // (reuse tSamples to avoid storing all individual values)
        float confSum = 0.0f;
        for (float t2 : tSamples) {
            {
                auto res = c.evaluateDerivativeWithConfidence(t2);
                float conf2 = res.second;
                confSum += conf2;
            }
        }
        float meanConf = confSum / static_cast<float>(tSamples.size());
        meanEfficiencies.emplace_back(name, meanEfficiency);
        meanConfVec.emplace_back(name, meanConf);
    }
    // Print summary table
    std::cout << "⇢ Mean Efficiency Summary:\n";
    for (const auto& [name, meanEff] : meanEfficiencies) {
        // Align name to width 11, value to 8 decimals
        std::cout << "   " << std::left << std::setw(11) << name
                  << std::right << std::fixed << std::setprecision(5)
                  << meanEff << std::endl;
    }
    // Compute and print normalized efficiency = mean(conf/kappa) / mean(conf)
    std::cout << "⇢ Normalized Efficiency (mean(conf/κ) / mean(conf)):\n";
    // Build a lookup for mean confidences
    auto findMeanConf = [&](const std::string& key) {
        for (const auto& kv : meanConfVec) {
            if (kv.first == key) return kv.second;
        }
        return 0.0f;
    };
    for (const auto& [name, meanEff] : meanEfficiencies) {
        float meanConf = findMeanConf(name);
        float normalized = (meanConf > 0.0f) ? (meanEff / meanConf) : 0.0f;
        std::cout << "   " << std::left << std::setw(11) << name
                  << std::right << std::fixed << std::setprecision(5)
                  << normalized << std::endl;
    }
}



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

    printTest("Linear 2D curve t=0", eval0 == p0);
    printTest("Linear 2D curve t=1", eval1 == p1);
    printTest("Linear 2D curve t=0.5", evalHalf == Point2{0.5f,0.5f});

    // --- 3D Linear Curve ---
    using Point3 = Point<float,3>;
    using Curve3 = Curve<float,3>;
    Point3 p3_0{0.0f,0.0f,0.0f};
    Point3 p3_1{1.0f,2.0f,3.0f};
    Curve3 linear3 = Curve3::linearCurve(p3_0,p3_1);

    Point3 eval3Half = linear3.evaluate(0.5f);
    printTest("Linear 3D curve t=0.5", eval3Half == Point3{0.5f,1.0f,1.5f});

    // --- Nonlinear Curve ---
    auto nonlinearFunc = [](float t) -> Point2 { return Point2{t, t*t}; };
    Curve2 nonlinear(nonlinearFunc,0.0f,1.0f);

    Point2 evalNonlin = nonlinear.evaluate(0.5f);
    printTest("Nonlinear 2D curve t=0.5", evalNonlin == Point2{0.5f,0.25f});

    // --- 3D Nonlinear Curve ---
    auto nonlinear3D = [](float t) -> Point3 { return Point3{t, t*t, t*t*t}; };
    Curve3 nonlinearCurve3D(nonlinear3D, 0.0f, 1.0f);
    Point3 eval3D_0 = nonlinearCurve3D.evaluate(0.0f);
    Point3 eval3D_half = nonlinearCurve3D.evaluate(0.5f);
    Point3 eval3D_1 = nonlinearCurve3D.evaluate(1.0f);
    printTest("Nonlinear 3D curve t=0.0", eval3D_0 == Point3{0.0f, 0.0f, 0.0f});
    printTest("Nonlinear 3D curve t=0.5", eval3D_half == Point3{0.5f, 0.25f, 0.125f});
    printTest("Nonlinear 3D curve t=1.0", eval3D_1 == Point3{1.0f, 1.0f, 1.0f});

    // --- 2D Sinusoidal Curve: y = sin(2*pi*x), x = t in [0,1] ---
    auto sinusoidalFunc = [](float t) -> Point2 {
        return Point2{t, std::sin(2.0f * float(M_PI) * t)};
    };
    Curve2 sinusoidal(sinusoidalFunc, 0.0f, 1.0f);
    printTest("Sinusoidal 2D curve t=0", sinusoidal.evaluate(0.0f) == Point2{0.0f, 0.0f});
    printTest("Sinusoidal 2D curve t=0.25", sinusoidal.evaluate(0.25f) == Point2{0.25f, 1.0f});
    printTest("Sinusoidal 2D curve t=0.5", sinusoidal.evaluate(0.5f) == Point2{0.5f, 0.0f});
    printTest("Sinusoidal 2D curve t=0.75", sinusoidal.evaluate(0.75f) == Point2{0.75f, -1.0f});
    printTest("Sinusoidal 2D curve t=1.0", sinusoidal.evaluate(1.0f) == Point2{1.0f, 0.0f});

    // --- Transform on Curve ---
    Eigen::Matrix2f scaleMat; scaleMat << 2.0f,0.0f,0.0f,3.0f;
    Eigen::Vector2f translation; translation << 1.0f,-1.0f;
    Affine<float,2> transform(scaleMat, translation);

    Curve2 transformedCurve = nonlinear.applyTransform(transform);
    Point2 evalTransformed = transformedCurve.evaluate(0.5f);
    Point2 expectedTransformed{2.0f, -0.25f + -1.0f}; // x scaled by 2, y by 3, then translation applied
    printTest("Curve transform evaluation", evalTransformed == Point2{2.0f, -0.25f});

    // --- Multi-dimensional Nonlinear Curve Tests (2D–5D) ---

    // 2D Nonlinear Curve
    auto func2D = [](float t) -> Point<float,2> { return Point<float,2>{-t, t * float(M_PI)}; };
    Curve<float,2> curve2(func2D, 0.0f, 1.0f);
    printTest("Nonlinear 2D t=0", curve2.evaluate(0.0f) == Point<float,2>{0.0f,0.0f});
    printTest("Nonlinear 2D t=0.5", curve2.evaluate(0.5f) == Point<float,2>{-0.5f, 0.5f*float(M_PI)});
    printTest("Nonlinear 2D t=1.0", curve2.evaluate(1.0f) == Point<float,2>{-1.0f, float(M_PI)});
    printTest("Nonlinear 2D t=-0.5", curve2.evaluate(-0.5f) == Point<float,2>{0.5f, -0.5f*float(M_PI)});
    printTest("Nonlinear 2D t=1.5", curve2.evaluate(1.5f) == Point<float,2>{-1.5f, 1.5f*float(M_PI)});

    // 3D Nonlinear Curve
    auto func3D = [](float t) -> Point<float,3> { return Point<float,3>{-t, t * float(M_PI), t*t}; };
    Curve<float,3> curve3(func3D, 0.0f, 1.0f);
    printTest("Nonlinear 3D t=0.5", curve3.evaluate(0.5f) == Point<float,3>{-0.5f, 0.5f*float(M_PI), 0.25f});
    printTest("Nonlinear 3D t=1.5", curve3.evaluate(1.5f) == Point<float,3>{-1.5f, 1.5f*float(M_PI), 2.25f});

    // 4D Nonlinear Curve
    auto func4D = [](float t) -> Point<float,4> { return Point<float,4>{-t, t * float(M_PI), t*t, t*t*t}; };
    Curve<float,4> curve4(func4D, 0.0f, 1.0f);
    printTest("Nonlinear 4D t=0.5", curve4.evaluate(0.5f) == Point<float,4>{-0.5f, 0.5f*float(M_PI), 0.25f, 0.125f});
    printTest("Nonlinear 4D t=-0.5", curve4.evaluate(-0.5f) == Point<float,4>{0.5f, -0.5f*float(M_PI), 0.25f, -0.125f});

    // 5D Nonlinear Curve
    auto func5D = [](float t) -> Point<float,5> { return Point<float,5>{-t, t * float(M_PI), t*t, t*t*t, std::sqrt(t+1.0f)}; };
    Curve<float,5> curve5(func5D, 0.0f, 1.0f);
    printTest("Nonlinear 5D t=0.5", curve5.evaluate(0.5f) == Point<float,5>{-0.5f, 0.5f*float(M_PI), 0.25f, 0.125f, std::sqrt(1.5f)});
    printTest("Nonlinear 5D t=1.5", curve5.evaluate(1.5f) == Point<float,5>{-1.5f, 1.5f*float(M_PI), 2.25f, 3.375f, std::sqrt(2.5f)});


    // Run derivative test suite
    testCurveDerivative();
    testCurveDerivativeResolution();
    testCurveDerivativeResolutionSweep();
    testCurveDerivativeCurvatureCorrelation();
    
    // Run integral test suite
    testCurveIntegral();
    testCurveIntegralResolutionSweep();

    // Run curvature test suite
    testCurveCurvature();
    
    
}

