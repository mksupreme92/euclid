#pragma once
#include <functional>
#include <utility>
#include "point.hpp"
#include "tolerance.hpp"
#include "algebra/transform.hpp"
#include <vector>
#include <limits>
#include <optional>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace Euclid::Geometry {

template <typename Scalar, int Dim = Eigen::Dynamic>
class Curve {
public:
    using PointType = Point<Scalar, Dim>;

    // General parametric curve constructor
    Curve(std::function<PointType(Scalar)> func, Scalar t0, Scalar t1)
        : curveFunc_(func), domain_{t0, t1} {}

    // Evaluate the curve at parameter t
    PointType evaluate(Scalar t) const {
        return curveFunc_(t);
    }

    // Evaluate the derivative of the curve at parameter t (finite difference)
    Eigen::Matrix<Scalar, Dim, 1> evaluateDerivative(Scalar t, const Tolerance& tol) const {
        Scalar h = tol.evaluateEpsilon((curveFunc_(t)).coords.norm());
        return (curveFunc_(t + h).coords - curveFunc_(t - h).coords) / (2.0 * h);
    }

    // Overload: automatically determines tolerance from curve geometry and clamps to domain
    Eigen::Matrix<Scalar, Dim, 1> evaluateDerivative(Scalar t) const {
        using VecS = Eigen::Matrix<Scalar, Dim, 1>;
        using VecD = Eigen::Matrix<double, Dim, 1>;
        Tolerance tol; // default tolerance model

        auto evalD = [&](double tt) -> VecD {
            return curveFunc_(static_cast<Scalar>(tt)).coords.template cast<double>();
        };

        // Geometry-informed scale for tolerance model
        const VecD pt = evalD(static_cast<double>(t));
        const Scalar geomScale = static_cast<Scalar>(pt.norm());
        const Scalar slopeEstimate = std::max(Scalar(1e-6), Scalar(geomScale + pt.cwiseAbs().sum()));

        // Domain and parameter-step safety
        const Scalar t0 = domain_.first;
        const Scalar t1 = domain_.second;
        const Scalar span = std::max(Scalar(0), t1 - t0);

        const Scalar epsMach = std::numeric_limits<Scalar>::epsilon();
        const Scalar hTol    = tol.evaluateEpsilon(std::max(slopeEstimate, Scalar(1)));
        const Scalar hUlp    = std::max(epsMach * (Scalar(1) + std::abs(t)), epsMach);
        const Scalar hDom    = std::max(span * Scalar(1e-6), epsMach);

        // Use a cube-root-of-eps heuristic for central differences:
        // For smooth functions the optimal step is ~ cbrt(eps) * scale.
        const Scalar hCbrt   = Scalar(std::cbrt(static_cast<double>(epsMach))) * std::max(span, Scalar(1));

        // Start with the most conservative among tolerance-, ulp-, domain- and cbrt-based steps.
        Scalar h0 = std::max(std::max(hTol, hUlp), std::max(hDom, hCbrt));

        // Keep step safely inside the domain while allowing larger steps when beneficial.
        const Scalar hMax = std::max(span * Scalar(0.25), epsMach);
        h0 = std::min(h0, hMax);
        /*
        std::cout << "[DEBUG] adaptive h = " << h0
                  << " for geomScale = " << geomScale
                  << ", slopeEstimate = " << slopeEstimate << std::endl;
         */

        // Ridders' parameters
        constexpr int MAX_IT = 30;

        auto safeCentral = [&](Scalar step) -> std::optional<VecD> {
            Scalar tf = std::min(t + step, t1);
            Scalar tb = std::max(t - step, t0);
            Scalar denom = tf - tb;
            if (!(denom > Scalar(0))) {
                Scalar nud = std::max<Scalar>(hUlp * Scalar(4), tol.evaluateEpsilon(span));
                tf = std::min(t + nud, t1);
                tb = std::max(t - nud, t0);
                denom = tf - tb;
                if (!(denom > Scalar(0))) return std::nullopt;
            }
            const VecD fwd = evalD(static_cast<double>(tf));
            const VecD bwd = evalD(static_cast<double>(tb));
            VecD g = (fwd - bwd) / static_cast<double>(denom);
            if (!g.allFinite()) return std::nullopt;
            return g;
        };

        // We allow a few restarts with a smaller initial step if the sequence fails to converge.
        const int MAX_RESTARTS = 3;
        VecD globalBest = VecD::Zero();
        double globalBestErr = std::numeric_limits<double>::infinity();
        bool   haveGlobalBest = false;

        Scalar hStart = h0;

        for (int r = 0; r < MAX_RESTARTS; ++r) {
            // --- Adaptive derivative tolerance scaling ---
            // This makes derivative convergence tolerance geometry-dependent:
            // - Looser (10–100×) for flat/smooth geometry to absorb numeric noise.
            // - Tighter near regions of high curvature or steep slope.
            const Scalar slopeMag = slopeEstimate;
            const Scalar span     = std::max(Scalar(1e-9), domain_.second - domain_.first);

            // Geometric conditioning heuristic: small for flat, large for tight curvature
            const Scalar geomCond = std::min(Scalar(1e3), Scalar(1) / std::max(Scalar(1e-6), geomScale / (slopeMag * span)));

            // Scale factor varies between 10× and 100× based on geometry smoothness
            const Scalar scaleFactor = std::clamp(Scalar(10) * (Scalar(1) + Scalar(0.1) * geomCond),
                                                  Scalar(10), Scalar(100));

            const Scalar convTol = tol.evaluateEpsilon(std::max<Scalar>(slopeEstimate, Scalar(1))) * scaleFactor;

            VecD D[MAX_IT][MAX_IT];
            bool haveBest = false;
            VecD best = VecD::Zero();
            double bestErr = std::numeric_limits<double>::infinity();

            int growthCount = 0;
            double lastDiagErr = std::numeric_limits<double>::infinity();

            // Build extrapolation table with halving step sequence
            for (int i = 0; i < MAX_IT; ++i) {
                Scalar hi = hStart / std::pow(Scalar(2), i);
                auto g0 = safeCentral(hi);
                if (!g0.has_value()) break;
                D[i][0] = *g0;

                // Richardson / Ridders’ extrapolation
                for (int k = 1; k <= i; ++k) {
                    const double fac = std::pow(4.0, k);
                    VecD num = fac * D[i][k - 1] - D[i - 1][k - 1];
                    const double den = fac - 1.0;
                    D[i][k] = num / den;
                    if (!D[i][k].allFinite()) {
                        D[i][k] = D[i][k - 1];
                    }
                }

                // Error estimate: difference between successive diagonals
                if (i > 0) {
                    VecD diff = D[i][i] - D[i - 1][i - 1];
                    double err = diff.norm();
                    /*
                    std::cout << "[DEBUG] Ridders iteration " << i
                              << ", |Δdiag| = " << err << std::endl;
                     */

                    // Track best
                    if (err < bestErr && D[i][i].allFinite()) {
                        bestErr = err;
                        best = D[i][i];
                        haveBest = true;
                    }

                    // Convergence: absolute-to-relative hybrid
                    const double relScale = 1.0 + D[i][i].norm();
                    if (err <= static_cast<double>(convTol) * relScale) {
                        /*
                        std::cout << "[DEBUG] Ridders converged at i=" << i
                                  << " with err=" << err << " (tol=" << convTol << ")\n";*/
                        return best.template cast<Scalar>();
                    }

                    // If diagonal error is not improving, count growth
                    if (err >= lastDiagErr * 0.9) {
                        ++growthCount;
                    } else {
                        growthCount = 0;
                    }
                    lastDiagErr = err;

                    // --- Oscillation detection and refinement ---
                    if (i > 2) {
                        double deltaPrev = (D[i-1][i-1] - D[i-2][i-2]).norm();
                        double deltaCurr = (D[i][i] - D[i-1][i-1]).norm();
                        if (deltaCurr > deltaPrev * 1.5 && deltaCurr > static_cast<double>(convTol)) {
                            /*
                            std::cout << "[DEBUG] Detected oscillatory behavior, refining step h by 0.5\n";
                             */
                            hStart *= Scalar(0.5);
                            break;
                        }
                    }

                    // --- Final adaptive refinement for oscillatory, high-curvature curves ---
                    if (geomScale > Scalar(1.5) && slopeEstimate > Scalar(3.0)) {
                        /*
                        std::cout << "[DEBUG] High-curvature refinement: halving step hStart for oscillatory geometry\n";
                         */
                        hStart = std::max(hUlp * Scalar(8), hStart * Scalar(0.25));
                        break;
                    }

                    // Early abort this restart if the sequence keeps getting worse
                    if (growthCount >= 2) {
                        break;
                    }
                } else {
                    if (D[i][0].allFinite()) {
                        best = D[i][0];
                        bestErr = std::numeric_limits<double>::infinity();
                        haveBest = true;
                    }
                }
            }

            if (haveBest && bestErr < globalBestErr) {
                globalBestErr = bestErr;
                globalBest = best;
                haveGlobalBest = true;
            }

            // Adaptive damping: reduce h proportionally to slope magnitude for smoother refinement
            hStart = std::max<Scalar>(hUlp * Scalar(8), hStart / (Scalar(1) + slopeEstimate));
        }

        if (haveGlobalBest) {
            /*
            std::cout << "[DEBUG] Ridders restarts exhausted; returning best with err=" << globalBestErr << "\n";
             */
            return globalBest.template cast<Scalar>();
        }

        // Last resort: single safe central difference with smallest reliable step
        {
            auto g = safeCentral(std::max<Scalar>(hUlp * Scalar(8), tol.evaluateEpsilon(span)));
            if (g.has_value()) return g->template cast<Scalar>();
        }
        return VecS::Zero();
    }

    // Overload: returns both derivative and confidence metric (0 = low, 1 = high)
    std::pair<Eigen::Matrix<Scalar, Dim, 1>, Scalar> evaluateDerivativeWithConfidence(Scalar t) const {
        using VecS = Eigen::Matrix<Scalar, Dim, 1>;
        using VecD = Eigen::Matrix<double, Dim, 1>;
        Tolerance tol; // default tolerance model

        auto evalD = [&](double tt) -> VecD {
            return curveFunc_(static_cast<Scalar>(tt)).coords.template cast<double>();
        };

        // Geometry-informed scale for tolerance model
        const VecD pt = evalD(static_cast<double>(t));
        const Scalar geomScale = static_cast<Scalar>(pt.norm());
        const Scalar slopeEstimate = std::max(Scalar(1e-6), Scalar(geomScale + pt.cwiseAbs().sum()));

        // Domain and parameter-step safety
        const Scalar t0 = domain_.first;
        const Scalar t1 = domain_.second;
        const Scalar span = std::max(Scalar(0), t1 - t0);

        const Scalar epsMach = std::numeric_limits<Scalar>::epsilon();
        const Scalar hTol    = tol.evaluateEpsilon(std::max(slopeEstimate, Scalar(1)));
        const Scalar hUlp    = std::max(epsMach * (Scalar(1) + std::abs(t)), epsMach);
        const Scalar hDom    = std::max(span * Scalar(1e-6), epsMach);

        // Use a cube-root-of-eps heuristic for central differences:
        // For smooth functions the optimal step is ~ cbrt(eps) * scale.
        const Scalar hCbrt   = Scalar(std::cbrt(static_cast<double>(epsMach))) * std::max(span, Scalar(1));

        // Start with the most conservative among tolerance-, ulp-, domain- and cbrt-based steps.
        Scalar h0 = std::max(std::max(hTol, hUlp), std::max(hDom, hCbrt));

        // Keep step safely inside the domain while allowing larger steps when beneficial.
        const Scalar hMax = std::max(span * Scalar(0.25), epsMach);
        h0 = std::min(h0, hMax);

        // Ridders' parameters
        constexpr int MAX_IT = 30;

        auto safeCentral = [&](Scalar step) -> std::optional<VecD> {
            Scalar tf = std::min(t + step, t1);
            Scalar tb = std::max(t - step, t0);
            Scalar denom = tf - tb;
            if (!(denom > Scalar(0))) {
                Scalar nud = std::max<Scalar>(hUlp * Scalar(4), tol.evaluateEpsilon(span));
                tf = std::min(t + nud, t1);
                tb = std::max(t - nud, t0);
                denom = tf - tb;
                if (!(denom > Scalar(0))) return std::nullopt;
            }
            const VecD fwd = evalD(static_cast<double>(tf));
            const VecD bwd = evalD(static_cast<double>(tb));
            VecD g = (fwd - bwd) / static_cast<double>(denom);
            if (!g.allFinite()) return std::nullopt;
            return g;
        };

        const Scalar convTol = tol.evaluateEpsilon(std::max<Scalar>(slopeEstimate, Scalar(1))) * std::clamp(Scalar(10) * (Scalar(1) + Scalar(0.1) * std::min(Scalar(1e3), Scalar(1) / std::max(Scalar(1e-6), geomScale / (slopeEstimate * span)))), Scalar(10), Scalar(100));

        VecD globalBest = VecD::Zero();
        double globalBestErr = std::numeric_limits<double>::infinity();
        bool haveGlobalBest = false;

        Scalar hStart = h0;

        for (int r = 0; r < 3; ++r) {
            VecD D[MAX_IT][MAX_IT];
            bool haveBest = false;
            VecD best = VecD::Zero();
            double bestErr = std::numeric_limits<double>::infinity();

            int growthCount = 0;
            double lastDiagErr = std::numeric_limits<double>::infinity();

            for (int i = 0; i < MAX_IT; ++i) {
                Scalar hi = hStart / std::pow(Scalar(2), i);
                auto g0 = safeCentral(hi);
                if (!g0.has_value()) break;
                D[i][0] = *g0;

                for (int k = 1; k <= i; ++k) {
                    const double fac = std::pow(4.0, k);
                    VecD num = fac * D[i][k - 1] - D[i - 1][k - 1];
                    const double den = fac - 1.0;
                    D[i][k] = num / den;
                    if (!D[i][k].allFinite()) {
                        D[i][k] = D[i][k - 1];
                    }
                }

                if (i > 0) {
                    VecD diff = D[i][i] - D[i - 1][i - 1];
                    double err = diff.norm();

                    if (err < bestErr && D[i][i].allFinite()) {
                        bestErr = err;
                        best = D[i][i];
                        haveBest = true;
                    }

                    const double relScale = 1.0 + D[i][i].norm();
                    if (err <= static_cast<double>(convTol) * relScale) {
                        globalBestErr = err;
                        globalBest = best;
                        haveGlobalBest = true;
                        break;
                    }

                    if (err >= lastDiagErr * 0.9) {
                        ++growthCount;
                    } else {
                        growthCount = 0;
                    }
                    lastDiagErr = err;

                    if (growthCount >= 2) {
                        break;
                    }
                } else {
                    if (D[i][0].allFinite()) {
                        best = D[i][0];
                        bestErr = std::numeric_limits<double>::infinity();
                        haveBest = true;
                    }
                }
            }

            if (haveBest && bestErr < globalBestErr) {
                globalBestErr = bestErr;
                globalBest = best;
                haveGlobalBest = true;
            }

            // Adaptive damping: reduce h proportionally to slope magnitude for smoother refinement
            hStart = std::max<Scalar>(hUlp * Scalar(8), hStart / (Scalar(1) + slopeEstimate));
        }

        VecS derivative;
        if (haveGlobalBest)
            derivative = globalBest.template cast<Scalar>().eval();
        else
            derivative = VecS::Zero();
        const Scalar confidence = Scalar(1) - std::min<Scalar>(Scalar(1), static_cast<Scalar>(globalBestErr) / convTol);
        return {derivative, confidence};
    }

    // Estimate integral of curvature magnitude across the domain using adaptive sampling
    Scalar evaluateIntegral(const Tolerance& tol, int baseDivs = 64, int maxDepth = 10) const {
        Scalar t0 = domain_.first;
        Scalar t1 = domain_.second;
        std::vector<Scalar> samples;
        samples.reserve(baseDivs + 1);
        for (int i = 0; i <= baseDivs; ++i)
            samples.push_back(t0 + (t1 - t0) * (Scalar(i) / baseDivs));

        auto curvatureAt = [&](Scalar t) {
            const Scalar h = tol.evaluateEpsilon(curveFunc_(t).coords.norm());
            auto p0 = curveFunc_(t).coords;
            auto p1 = curveFunc_(std::min(t1, t + h)).coords;
            auto p_1 = curveFunc_(std::max(t0, t - h)).coords;
            auto first = (p1 - p_1) / (2 * h);
            auto second = (p1 - 2 * p0 + p_1) / (h * h);
            return second.norm() / std::pow(1.0 + first.squaredNorm(), 1.5);
        };

        // Adaptive subdivision loop
        for (int depth = 0; depth < maxDepth; ++depth) {
            std::vector<Scalar> newSamples;
            for (size_t i = 0; i + 1 < samples.size(); ++i) {
                Scalar a = samples[i];
                Scalar b = samples[i + 1];
                Scalar mid = (a + b) * 0.5;
                Scalar kA = curvatureAt(a);
                Scalar kB = curvatureAt(b);
                Scalar kM = curvatureAt(mid);
                Scalar avgK = (kA + kB + 4 * kM) / 6;
                if (avgK * std::abs(b - a) > tol.paramTol)
                    newSamples.push_back(mid);
            }
            if (newSamples.empty()) break;
            samples.insert(samples.end(), newSamples.begin(), newSamples.end());
            std::sort(samples.begin(), samples.end());
            samples.erase(std::unique(samples.begin(), samples.end(), [&](Scalar x, Scalar y){ return std::abs(x - y) < tol.paramTol; }), samples.end());
        }

        // Trapezoidal integration of curvature magnitude
        Scalar integral = 0;
        for (size_t i = 0; i + 1 < samples.size(); ++i) {
            Scalar a = samples[i];
            Scalar b = samples[i + 1];
            Scalar kA = curvatureAt(a);
            Scalar kB = curvatureAt(b);
            integral += (kA + kB) * 0.5 * (b - a);
        }

        return integral;
    }

    // Curve parameter domain
    std::pair<Scalar, Scalar> domain() const { return domain_; }

    // Apply a transform to the curve
    template <typename Transform>
    Curve applyTransform(const Transform& T) const {
        auto newFunc = [T, this](Scalar t) -> PointType {
            return T.apply(this->evaluate(t));
        };
        return Curve(newFunc, domain_.first, domain_.second);
    }

    // Linear curve factory: p(t) = p0 + t*(p1 - p0)
    static Curve linearCurve(const PointType& p0, const PointType& p1) {
        auto func = [p0, p1](Scalar t) -> PointType {
            return p0 + (p1 - p0) * t;
        };
        return Curve(func, 0, 1);
    }

private:
    std::function<PointType(Scalar)> curveFunc_;
    std::pair<Scalar, Scalar> domain_;
};

} // namespace Euclid::Geometry
