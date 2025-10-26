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
}

