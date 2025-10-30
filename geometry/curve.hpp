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
#include <numeric>

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
        // Reduce the maximum finite difference step for oscillatory or high-curvature curves
        const Scalar hMax = std::max(span * Scalar(1e-2), epsMach);
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

    // Cached variant: avoids recomputing expensive adaptive derivative
    Eigen::Matrix<Scalar, Dim, 1> evaluateDerivativeCached(Scalar t) const {
        // If cache hits and param matches exactly, return cached value
        if (derivativeCache_.has_value() && derivativeCache_->first == t) {
            return derivativeCache_->second;
        }
        auto d = evaluateDerivative(t);
        derivativeCache_ = std::make_pair(t, d);
        return d;
    }

    // Utility: adaptive 1D sampling helper for intersection-style solvers
    template <typename Func>
    Scalar adaptiveParametricSolve(Func f, Scalar t0, Scalar t1, Scalar targetTol) const {
        Tolerance tol;
        Scalar a = std::max(domain_.first, t0);
        Scalar b = std::min(domain_.second, t1);
        if (!(a < b)) return a;

        Scalar fa = f(a);
        Scalar fb = f(b);

        // If either endpoint is already good, return it
        if (std::abs(fa) <= targetTol) return a;
        if (std::abs(fb) <= targetTol) return b;

        // Fallback to bisection-like refinement (no assumptions about monotonicity)
        const int maxIter = 32;
        for (int i = 0; i < maxIter; ++i) {
            Scalar m = (a + b) * Scalar(0.5);
            Scalar fm = f(m);
            if (std::abs(fm) <= targetTol) {
                return m;
            }

            // Choose side with smaller magnitude to continue refinement
            if (std::abs(fa) < std::abs(fb)) {
                b = m;
                fb = fm;
            } else {
                a = m;
                fa = fm;
            }

            if (std::abs(b - a) <= targetTol * Scalar(2)) {
                return (a + b) * Scalar(0.5);
            }
        }
        return (a + b) * Scalar(0.5);
    }

    // Evaluate the second derivative of the curve at parameter t (consistent with evaluateDerivative)
    Eigen::Matrix<Scalar, Dim, 1> evaluateSecondDerivative(Scalar t) const {
        using VecS = Eigen::Matrix<Scalar, Dim, 1>;
        using VecD = Eigen::Matrix<double, Dim, 1>;
        Tolerance tol;

        auto evalD = [&](double tt) -> VecD {
            return curveFunc_(static_cast<Scalar>(tt)).coords.template cast<double>();
        };

        const Scalar t0 = domain_.first;
        const Scalar t1 = domain_.second;
        const Scalar span = std::max(Scalar(1e-9), t1 - t0);
        const Scalar epsMach = std::numeric_limits<Scalar>::epsilon();

        // Compute adaptive step similar to evaluateDerivative
        const VecD pt = evalD(static_cast<double>(t));
        const Scalar geomScale = static_cast<Scalar>(pt.norm());
        const Scalar slopeEstimate = std::max(Scalar(1e-6), Scalar(geomScale + pt.cwiseAbs().sum()));
        const Scalar hTol = tol.evaluateEpsilon(std::max(slopeEstimate, Scalar(1)));
        const Scalar hCbrt = Scalar(std::cbrt(static_cast<double>(epsMach))) * std::max(span, Scalar(1));
        Scalar hStart = std::max(hTol, hCbrt);
        hStart = std::clamp(hStart, span * Scalar(1e-6), span * Scalar(1e-2));

        auto safeEval = [&](Scalar tVal) -> VecD {
            Scalar tc = std::clamp(tVal, t0, t1);
            return evalD(static_cast<double>(tc));
        };

        constexpr int MAX_IT = 20;
        VecD best = VecD::Zero();
        double bestErr = std::numeric_limits<double>::infinity();

        for (int i = 0; i < MAX_IT; ++i) {
            Scalar hi = hStart / std::pow(2.0, i);
            Scalar tp = std::min(t + hi, t1);
            Scalar tm = std::max(t - hi, t0);

            VecD fp = safeEval(tp);
            VecD fm = safeEval(tm);
            VecD f0 = safeEval(t);

            VecD g2 = (fp - 2.0 * f0 + fm) / (hi * hi);

            if (i > 0) {
                VecD diff = g2 - best;
                double err = diff.norm();
                if (err < bestErr && g2.allFinite()) {
                    bestErr = err;
                    best = g2;
                }
                if (err < tol.evaluateEpsilon(span)) break;
            } else {
                best = g2;
            }
        }

        if (!best.allFinite()) return VecS::Zero();
        return best.template cast<Scalar>();
    }

    // Evaluate the normalized tangent vector at parameter t
    Eigen::Matrix<Scalar, Dim, 1> evaluateTangent(Scalar t) const {
        using VecS = Eigen::Matrix<Scalar, Dim, 1>;
        Tolerance tol; // use the default adaptive tolerance model

        VecS v = evaluateDerivative(t);
        Scalar n = v.norm();

        // Avoid division by zero or subnormal magnitudes
        if (n < tol.paramTol) {
            return VecS::Zero();
        }

        return v / n;
    }

    // Overload: returns both derivative and confidence metric (0 = low, 1 = high)
    std::pair<Eigen::Matrix<Scalar, Dim, 1>, Scalar> evaluateDerivativeWithConfidence(Scalar t) const {
        using VecS = Eigen::Matrix<Scalar, Dim, 1>;
        Tolerance tol;

        // 1) Get our best derivative once (expensive Ridders path)
        VecS d0 = evaluateDerivative(t);
        const Scalar d0norm = d0.norm();

        // 2) Do two *cheap* directional checks using a small param perturbation, but
        //    compute derivatives via a *single-step central difference* rather than
        //    calling the full adaptive evaluator again. This preserves the intent of the
        //    confidence metric while avoiding 2 extra Ridders runs.
        const Scalar t0 = domain_.first;
        const Scalar t1 = domain_.second;
        const Scalar span = std::max(Scalar(1e-9), t1 - t0);
        const Scalar epsMach = std::numeric_limits<Scalar>::epsilon();
        // base perturbation tied to tolerance and span
        const Scalar baseDelta = std::max(Scalar(tol.paramTol), std::cbrt(epsMach) * span);
        const Scalar delta = std::min(baseDelta, span * Scalar(1e-3));

        auto cheapDerivative = [&](Scalar tc) -> VecS {
            // clamp inside domain
            const Scalar tp = std::min(tc + delta, t1);
            const Scalar tm = std::max(tc - delta, t0);
            const Scalar denom = tp - tm;
            if (denom <= Scalar(0)) {
                return d0; // fall back to main derivative
            }
            // evaluate curve once per side (cheap)
            const VecS fp = evaluate(tp).coords;
            const VecS fm = evaluate(tm).coords;
            return (fp - fm) / denom;
        };

        VecS dPlus  = cheapDerivative(std::min(t + delta, t1));
        VecS dMinus = cheapDerivative(std::max(t - delta, t0));

        // 3) Local fluctuation between the two cheap estimates
        const VecS fluct = (dPlus - dMinus) * Scalar(0.5);
        const Scalar localErr = fluct.norm();

        // 4) Map to confidence similarly to previous version, but cheaper
        const Scalar baseTol = tol.evaluateEpsilon(std::max(Scalar(1), d0norm));
        const Scalar floorTol = std::max(baseTol * Scalar(50), Scalar(1e-6));
        const Scalar scale = std::max(floorTol, localErr + floorTol);
        const Scalar ratio = localErr / scale;

        Scalar conf;
        if (ratio < Scalar(0.5)) {
            conf = Scalar(1) - ratio * Scalar(0.1);
        } else if (ratio < Scalar(2)) {
            conf = Scalar(0.9) - (ratio - Scalar(0.5)) * Scalar(0.3);
        } else {
            conf = Scalar(0.3) * std::exp(-ratio * Scalar(0.5));
        }
        conf = std::clamp(conf, Scalar(0.05), Scalar(1.0));

        return {d0, conf};
    }

    // Compute local curvature at parameter t:
    // κ(t) = ||a_perp|| / ||v||^2, where a_perp = a - ((a·v)/||v||^2) v
    Scalar evaluateCurvature(Scalar t) const {
        using VecS = Eigen::Matrix<Scalar, Dim, 1>;
        Tolerance tol;

        const Scalar t0 = domain_.first;
        const Scalar t1 = domain_.second;
        const Scalar tc = std::clamp(t, t0, t1);

        const Scalar eps = std::max(static_cast<Scalar>(tol.paramTol),
                                    static_cast<Scalar>(std::cbrt(std::numeric_limits<Scalar>::epsilon())));
        const Scalar h = std::clamp((domain_.second - domain_.first) * Scalar(1e-4),
                                    eps, (domain_.second - domain_.first) * Scalar(1e-2));

        // Evaluate points
        VecS f0  = evaluate(tc).coords;
        VecS fp  = evaluate(std::min(tc + h, domain_.second)).coords;
        VecS fm  = evaluate(std::max(tc - h, domain_.first)).coords;

        // First and second derivatives (standard central difference)
        VecS v = (fp - fm) / (Scalar(2) * h);
        VecS a = (fp - Scalar(2) * f0 + fm) / (h * h);

        Scalar vnorm = v.norm();
        if (vnorm < tol.paramTol) return Scalar(0);

        // Project out parallel component
        Scalar vdotv = vnorm * vnorm;
        VecS a_perp = a - (a.dot(v) / vdotv) * v;

        Scalar kappa = a_perp.norm() / vdotv;
        return kappa;
    }


    // Overload: Estimate arc length integral adaptively using internal tolerance model (panel-doubling Simpson-like)
    Scalar evaluateIntegral() const {
        Tolerance tol;
        const Scalar t0 = domain_.first;
        const Scalar t1 = domain_.second;
        if (!(t0 < t1)) return Scalar(0);

        // speed(t) = |p'(t)|, same as before
        auto speed = [&](Scalar t) -> Scalar {
            Scalar v = evaluateDerivative(t).norm();
            if (!std::isfinite(v)) v = Scalar(0);
            return v;
        };

        const Scalar span = t1 - t0;

        // Estimate oscillation by sampling speed at 8 intervals
        Scalar osc = Scalar(0);
        Scalar prev = speed(t0);
        const int oscSamples = 8;
        for (int j = 1; j <= oscSamples; ++j) {
            Scalar tj = t0 + span * Scalar(j) / Scalar(oscSamples);
            Scalar s = speed(tj);
            osc += std::abs(s - prev);
            prev = s;
        }
        // Oscillation threshold heuristic: tighten tolerance for rougher curves
        const Scalar oscThreshold = Scalar(1.0);
        Scalar oscFactor = Scalar(1.0);
        if (osc > oscThreshold) {
            oscFactor = Scalar(0.1) / std::min(Scalar(10.0), osc);
        }

        // target accuracy for whole interval (geometry-aware and oscillation-aware)
        const Scalar baseTol = tol.evaluateEpsilon(std::max(span, Scalar(1)));
        const Scalar target = std::max(baseTol * Scalar(10) * oscFactor, Scalar(1e-6));

        // basic adaptive panel-doubling (composite Simpson-like)
        int panels = 8;               // start cheap
        const int maxPanels = 1 << 12; // 4096 max
        Scalar prevInt = std::numeric_limits<Scalar>::quiet_NaN();

        while (panels <= maxPanels) {
            const Scalar h = span / Scalar(panels);
            Scalar integral = Scalar(0);

            // Composite Simpson-style accumulation
            Scalar f0 = speed(t0);
            Scalar f1 = speed(t1);
            integral += f0 + f1;

            for (int i = 1; i < panels; ++i) {
                const Scalar ti = t0 + h * Scalar(i);
                const Scalar fi = speed(ti);
                // Simpson weights 4,2,4,2,...
                integral += (i % 2 ? Scalar(4) : Scalar(2)) * fi;
            }
            integral *= (h / Scalar(3));

            if (std::isfinite(prevInt)) {
                Scalar diff = std::abs(integral - prevInt);
                if (diff <= target) {
                    return integral;
                }
            }

            prevInt = integral;
            panels *= 2;
        }

        // if we exit the loop, return the last/best estimate
        return prevInt;
    }
    
    // Compute an axis-aligned bounding box of the curve adaptively using tolerance model
    std::pair<PointType, PointType> boundingBox() const {
        Tolerance tol;
        const Scalar t0 = domain_.first;
        const Scalar t1 = domain_.second;
        const Scalar span = std::max(Scalar(1e-9), t1 - t0);

        // Initial bounding box using endpoints
        PointType p0 = evaluate(t0);
        Eigen::Matrix<Scalar, Dim, 1> minC = p0.coords;
        Eigen::Matrix<Scalar, Dim, 1> maxC = p0.coords;

        // Estimate curvature-driven refinement
        Scalar k0 = evaluateCurvature(t0 + 0.1 * span);
        Scalar k1 = evaluateCurvature(t0 + 0.5 * span);
        Scalar k2 = evaluateCurvature(t0 + 0.9 * span);
        Scalar avgK = (k0 + k1 + k2) / 3.0;

        // Compute adaptive sampling density based on tolerance and curvature
        Scalar res = tol.evaluateEpsilon(span);
        int samples = std::clamp(
            int(16 + std::log10(1.0 + avgK * span / res) * 32),
            16, 512
        );

        // Sample adaptively along curve
        for (int i = 1; i <= samples; ++i) {
            Scalar t = t0 + (t1 - t0) * (Scalar(i) / Scalar(samples));
            PointType p = evaluate(t);
            minC = minC.cwiseMin(p.coords);
            maxC = maxC.cwiseMax(p.coords);
        }

        return { PointType(minC), PointType(maxC) };
    }

    // Subdivide curve into two segments at parameter tSplit
    std::pair<Curve, Curve> subDivide(Scalar tSplit) const {
        const Scalar tMin = domain_.first;
        const Scalar tMax = domain_.second;
        const Scalar t0 = std::max(tMin, std::min(tSplit, tMax));

        // Reparameterize sub-curves into [0, 1] domain
        auto func1 = [=, this](Scalar t) -> PointType {
            Scalar u = tMin + (t0 - tMin) * t;
            return this->curveFunc_(u);
        };
        auto func2 = [=, this](Scalar t) -> PointType {
            Scalar u = t0 + (tMax - t0) * t;
            return this->curveFunc_(u);
        };

        Curve c1(func1, 0.0, 1.0);
        Curve c2(func2, 0.0, 1.0);

        return {c1, c2};
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
    mutable std::optional<std::pair<Scalar, Eigen::Matrix<Scalar, Dim, 1>>> derivativeCache_;
};

} // namespace Euclid::Geometry
