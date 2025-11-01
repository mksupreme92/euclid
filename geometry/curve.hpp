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

// ---- Optional performance timing instrumentation ----
#ifdef EUCLID_DEBUG_TIMING_CURVE
#  include <chrono>
#  define CURVE_TIME_START() auto __curve_time_start = std::chrono::high_resolution_clock::now();
#  define CURVE_TIME_END(label) \
    do { \
        auto __curve_time_end = std::chrono::high_resolution_clock::now(); \
        auto __curve_time_dur = std::chrono::duration_cast<std::chrono::microseconds>(__curve_time_end - __curve_time_start).count(); \
        std::cerr << "[CURVE_TIMING] " << (label) << " took " << __curve_time_dur << " us\n"; \
    } while(0)
#else
#  define CURVE_TIME_START()
#  define CURVE_TIME_END(label)
#endif

namespace Euclid::Geometry {

// Curve now takes its function as a template parameter, for zero runtime overhead.
// The default (type-erased) variant is still supported for backwards compatibility.
template <typename Scalar, int Dim = Eigen::Dynamic, typename FuncT = void>
class Curve;

// Primary template for inline function type (FuncT != void)
template <typename Scalar, int Dim, typename FuncT>
class Curve {
public:
    using PointType = Point<Scalar, Dim>;
    // General parametric curve constructor: stores callable directly (no std::function)
    Curve(FuncT func, Scalar t0, Scalar t1)
        : curveFunc_(std::move(func)), domain_{t0, t1} {}

    // Evaluate the curve at parameter t
    PointType evaluate(Scalar t) const {
        CURVE_TIME_START();
        auto result = curveFunc_(t);
        CURVE_TIME_END("evaluate");
        return result;
    }


    // Overload: automatically determines tolerance from curve geometry and clamps to domain
    Eigen::Matrix<Scalar, Dim, 1> evaluateDerivative(Scalar t) const {
        CURVE_TIME_START();
        using VecS = Eigen::Matrix<Scalar, Dim, 1>;
        using VecD = Eigen::Matrix<double, Dim, 1>;
        Tolerance tol; // default tolerance model

        // Inline evaluation function (direct call, no std::function overhead)
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

        // Ridders' parameters
        constexpr int MAX_IT = 30;

        // --- Local memoization cache for f(t+h), f(t-h) to avoid redundant evaluations ---
        struct FCache {
            double tp = std::numeric_limits<double>::quiet_NaN();
            double tm = std::numeric_limits<double>::quiet_NaN();
            VecD fp;
            VecD fm;
        };
        FCache cache;
        // Returns (fp, fm) for given step h, using cache if possible.
        auto getSymmetric = [&](Scalar step) -> std::pair<VecD, VecD> {
            double tp = std::min(static_cast<double>(t + step), static_cast<double>(t1));
            double tm = std::max(static_cast<double>(t - step), static_cast<double>(t0));
            VecD fp, fm;
            // Only recompute if not cached
            if (tp == cache.tp) {
                fp = cache.fp;
            } else {
                fp = evalD(tp);
                cache.tp = tp;
                cache.fp = fp;
            }
            if (tm == cache.tm) {
                fm = cache.fm;
            } else {
                fm = evalD(tm);
                cache.tm = tm;
                cache.fm = fm;
            }
            return {fp, fm};
        };

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
            // Use cache for symmetric evaluations
            auto [fwd, bwd] = getSymmetric(step);
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
                    // Track best
                    if (err < bestErr && D[i][i].allFinite()) {
                        bestErr = err;
                        best = D[i][i];
                        haveBest = true;
                    }

                    // Convergence: absolute-to-relative hybrid
                    const double relScale = 1.0 + D[i][i].norm();
                    if (err <= static_cast<double>(convTol) * relScale) {
                        // Return best estimate
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
                            // Detected oscillatory behavior, refining step h by 0.5
                            hStart *= Scalar(0.5);
                            break;
                        }
                    }

                    // --- Final adaptive refinement for oscillatory, high-curvature curves ---
                    if (geomScale > Scalar(1.5) && slopeEstimate > Scalar(3.0)) {
                        // High-curvature refinement: halving step hStart for oscillatory geometry
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
            CURVE_TIME_END("evaluateDerivative");
            return globalBest.template cast<Scalar>();
        }

        // Last resort: single safe central difference with smallest reliable step
        {
            auto g = safeCentral(std::max<Scalar>(hUlp * Scalar(8), tol.evaluateEpsilon(span)));
            if (g.has_value()) {
                CURVE_TIME_END("evaluateDerivative");
                return g->template cast<Scalar>();
            }
        }
        CURVE_TIME_END("evaluateDerivative");
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

    // Evaluate the second derivative of the curve at parameter t (efficient adaptive implementation)
    Eigen::Matrix<Scalar, Dim, 1> evaluateSecondDerivative(Scalar t) const {
        CURVE_TIME_START();
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

        VecD best = VecD::Zero();
        double bestErr = std::numeric_limits<double>::infinity();

        Scalar hi = hStart;
        VecD f0 = safeEval(t);
        VecD fp = VecD::Zero(), fm = VecD::Zero();
        for (int i = 0; i < 5; ++i) {  // increased iterations for stability
            Scalar tp = std::min(t + hi, t1);
            Scalar tm = std::max(t - hi, t0);
            // Cache symmetric evaluations if possible
            if (i == 0) {
                fp = safeEval(tp);
                fm = safeEval(tm);
            } else {
                fp = safeEval(tp);
                fm = safeEval(tm);
            }
            VecD g2 = (fp - 2.0 * f0 + fm) / (hi * hi);
            if (i > 0) {
                VecD diff = g2 - best;
                double err = diff.norm();
                if (err < bestErr && g2.allFinite()) {
                    bestErr = err;
                    best = g2;
                }
                if (err < tol.evaluateEpsilon(span)) break;
                if (i >= 3 && err > bestErr * 0.5) {
                    break; // stop if not improving
                }
            } else {
                best = g2;
            }
            hi *= Scalar(0.5);
        }
        if (!best.allFinite()) {
            CURVE_TIME_END("evaluateSecondDerivative");
            return VecS::Zero();
        }
        CURVE_TIME_END("evaluateSecondDerivative");
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
        CURVE_TIME_START();
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
        CURVE_TIME_END("evaluateCurvature");
        return kappa;
    }


    // Overload: Estimate arc length integral adaptively using internal tolerance model (panel-doubling Simpson-like)
    Scalar evaluateIntegral() const {
        CURVE_TIME_START();
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
                    CURVE_TIME_END("evaluateIntegral");
                    return integral;
                }
            }

            prevInt = integral;
            panels *= 2;
        }

        // if we exit the loop, return the last/best estimate
        CURVE_TIME_END("evaluateIntegral");
        return prevInt;
    }
    
    // Compute an axis-aligned bounding box of the curve adaptively using tolerance model
    std::pair<PointType, PointType> boundingBox() const {
        // If cache is valid and we have sample points, reuse both
        if (bboxValid_ && !bboxSamples_.empty()) {
            return cachedBox_;
        }
        CURVE_TIME_START();
        Tolerance tol;
        const Scalar t0 = domain_.first;
        const Scalar t1 = domain_.second;
        const Scalar span = std::max(Scalar(1e-9), t1 - t0);

        // --- Timing: curvature ---
#ifdef EUCLID_DEBUG_TIMING_CURVE
        auto __curve_curvature_start = std::chrono::high_resolution_clock::now();
#endif
        // Estimate curvature-driven refinement
        Scalar k0 = evaluateCurvature(t0 + 0.1 * span);
        Scalar k1 = evaluateCurvature(t0 + 0.5 * span);
        Scalar k2 = evaluateCurvature(t0 + 0.9 * span);
        Scalar avgK = (k0 + k1 + k2) / 3.0;
#ifdef EUCLID_DEBUG_TIMING_CURVE
        auto __curve_curvature_end = std::chrono::high_resolution_clock::now();
        auto __curve_curvature_dur = std::chrono::duration_cast<std::chrono::microseconds>(__curve_curvature_end - __curve_curvature_start).count();
        std::cerr << "[CURVE_TIMING] boundingBox:curvature took " << __curve_curvature_dur << " us\n";
#endif

        // --- Timing: density ---
#ifdef EUCLID_DEBUG_TIMING_CURVE
        auto __curve_density_start = std::chrono::high_resolution_clock::now();
#endif
        // Compute adaptive sampling density based on tolerance and curvature
        Scalar res = tol.evaluateEpsilon(span);
        int samples = std::clamp(
            int(8 + 24 * std::log2(1.0 + avgK * span / res)),
            8, 256
        );
#ifdef EUCLID_DEBUG_TIMING_CURVE
        auto __curve_density_end = std::chrono::high_resolution_clock::now();
        auto __curve_density_dur = std::chrono::duration_cast<std::chrono::microseconds>(__curve_density_end - __curve_density_start).count();
        std::cerr << "[CURVE_TIMING] boundingBox:density took " << __curve_density_dur << " us\n";
#endif

        // --- Two-pass adaptive sampling ---
#ifdef EUCLID_DEBUG_TIMING_CURVE
        auto __curve_sampling_start = std::chrono::high_resolution_clock::now();
#endif
        bboxSamples_.clear();
        // First pass: coarse sampling
        constexpr int coarseSamples = 32;
        std::vector<Scalar> coarseT(coarseSamples + 1);
        std::vector<Eigen::Matrix<Scalar, Dim, 1>> coarsePts(coarseSamples + 1);
        // Always include the left endpoint (t0)
        for (int i = 0; i <= coarseSamples; ++i) {
            Scalar t = t0 + (t1 - t0) * (Scalar(i) / Scalar(coarseSamples));
            coarseT[i] = t;
            coarsePts[i] = evaluate(t).coords;
        }
        // Compute initial min/max from coarse
        Eigen::Matrix<Scalar, Dim, 1> minC = coarsePts[0];
        Eigen::Matrix<Scalar, Dim, 1> maxC = coarsePts[0];
        for (int i = 1; i <= coarseSamples; ++i) {
            minC = minC.cwiseMin(coarsePts[i]);
            maxC = maxC.cwiseMax(coarsePts[i]);
        }
        // Compute box diagonal for thresholding
        Eigen::Matrix<Scalar, Dim, 1> boxDiag = maxC - minC;
        Scalar diagLen = boxDiag.norm();
        Scalar deltaThresh = diagLen * Scalar(0.05); // 5% of diagonal
        // Second pass: refine intervals with large deltas
        // Prepare the final sample points: start with coarse samples
        std::vector<std::pair<Scalar, Eigen::Matrix<Scalar, Dim, 1>>> allSamples;
        allSamples.reserve((coarseSamples + 1) * 2); // estimate upper bound
        for (int i = 0; i <= coarseSamples; ++i) {
            allSamples.emplace_back(coarseT[i], coarsePts[i]);
        }
        // Mark intervals to refine
        std::vector<std::pair<Scalar, Scalar>> refineIntervals;
        for (int i = 0; i < coarseSamples; ++i) {
            Scalar delta = (coarsePts[i + 1] - coarsePts[i]).norm();
            if (delta > deltaThresh) {
                refineIntervals.emplace_back(coarseT[i], coarseT[i + 1]);
            }
        }
        // For each interval, sample at higher density
        int fineSamples = std::max(samples, coarseSamples * 2); // use at least as many as main
        for (const auto& interval : refineIntervals) {
            Scalar ta = interval.first;
            Scalar tb = interval.second;
            // Avoid duplicating endpoints (they are already in allSamples)
            int nFine = std::max(4, fineSamples / coarseSamples);
            for (int j = 1; j < nFine; ++j) { // skip endpoints
                Scalar t = ta + (tb - ta) * (Scalar(j) / Scalar(nFine));
                Eigen::Matrix<Scalar, Dim, 1> pt = evaluate(t).coords;
                allSamples.emplace_back(t, pt);
            }
        }
        // Sort allSamples by t to preserve order
        std::sort(allSamples.begin(), allSamples.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        // Remove duplicate t values (from overlapping intervals)
        auto last = std::unique(allSamples.begin(), allSamples.end(),
                                [](const auto& a, const auto& b) {
                                    return std::abs(a.first - b.first) < Scalar(1e-12);
                                });
        allSamples.erase(last, allSamples.end());
        // Compute final min/max from allSamples
        minC = allSamples[0].second;
        maxC = allSamples[0].second;
        bboxSamples_.clear();
        bboxSamples_.reserve(allSamples.size());
        for (const auto& samp : allSamples) {
            minC = minC.cwiseMin(samp.second);
            maxC = maxC.cwiseMax(samp.second);
            bboxSamples_.push_back(samp.second);
        }
#ifdef EUCLID_DEBUG_TIMING_CURVE
        auto __curve_sampling_end = std::chrono::high_resolution_clock::now();
        auto __curve_sampling_dur = std::chrono::duration_cast<std::chrono::microseconds>(__curve_sampling_end - __curve_sampling_start).count();
        std::cerr << "[CURVE_TIMING] boundingBox:sampling took " << __curve_sampling_dur << " us\n";
#endif

        auto result = std::make_pair(PointType(minC), PointType(maxC));
        cachedBox_ = result;
        bboxValid_ = true;
        CURVE_TIME_END("boundingBox");
        return result;
    }

    // Invalidate the cached bounding box (call when control points or transform change)
    inline void invalidateBoundingBoxCache() const {
        bboxValid_ = false;
        bboxSamples_.clear();
    }

    // Subdivide curve into two segments at parameter tSplit
    std::pair<Curve, Curve> subDivide(Scalar tSplit) const {
        const Scalar tMin = domain_.first;
        const Scalar tMax = domain_.second;
        const Scalar t0 = std::max(tMin, std::min(tSplit, tMax));
        // Reparameterize sub-curves into [0, 1] domain
        auto func1 = [=, *this](Scalar t) -> PointType {
            Scalar u = tMin + (t0 - tMin) * t;
            return this->curveFunc_(u);
        };
        auto func2 = [=, *this](Scalar t) -> PointType {
            Scalar u = t0 + (tMax - t0) * t;
            return this->curveFunc_(u);
        };
        Curve<Scalar, Dim, decltype(func1)> c1(func1, 0.0, 1.0);
        Curve<Scalar, Dim, decltype(func2)> c2(func2, 0.0, 1.0);
        return {c1, c2};
    }

    // Curve parameter domain
    std::pair<Scalar, Scalar> domain() const { return domain_; }

    // Apply a transform to the curve
    template <typename Transform>
    auto applyTransform(const Transform& T) const {
        auto newFunc = [T, *this](Scalar t) -> PointType {
            return T.apply(this->evaluate(t));
        };
        return Curve<Scalar, Dim, decltype(newFunc)>(newFunc, domain_.first, domain_.second);
    }

    // Linear curve factory: p(t) = p0 + t*(p1 - p0)
    static Curve linearCurve(const PointType& p0, const PointType& p1) {
        auto func = [p0, p1](Scalar t) -> PointType {
            return p0 + (p1 - p0) * t;
        };
        return Curve<Scalar, Dim, decltype(func)>(func, 0, 1);
    }

private:
    FuncT curveFunc_;
    std::pair<Scalar, Scalar> domain_;
    mutable std::optional<std::pair<Scalar, Eigen::Matrix<Scalar, Dim, 1>>> derivativeCache_;
    // Bounding box cache mechanism
    mutable bool bboxValid_ = false;
    mutable std::pair<PointType, PointType> cachedBox_;
    // Store sampled points for bounding box computation (for reuse)
    mutable std::vector<Eigen::Matrix<Scalar, Dim, 1>> bboxSamples_;
};

// Specialization for legacy type-erased std::function (FuncT = void)
template <typename Scalar, int Dim>
class Curve<Scalar, Dim, void> {
public:
    using PointType = Point<Scalar, Dim>;
    Curve(std::function<PointType(Scalar)> func, Scalar t0, Scalar t1)
        : curveFunc_(std::move(func)), domain_{t0, t1} {}
    PointType evaluate(Scalar t) const {
        CURVE_TIME_START();
        auto result = curveFunc_(t);
        CURVE_TIME_END("evaluate");
        return result;
    }
    Eigen::Matrix<Scalar, Dim, 1> evaluateDerivative(Scalar t) const {
        // Use the new inline version with std::function as the callable.
        // This will be a bit slower, but preserves API for legacy code.
        Curve<Scalar, Dim, std::function<PointType(Scalar)>> tmp(curveFunc_, domain_.first, domain_.second);
        return tmp.evaluateDerivative(t);
    }
    Eigen::Matrix<Scalar, Dim, 1> evaluateDerivativeCached(Scalar t) const {
        if (derivativeCache_.has_value() && derivativeCache_->first == t) {
            return derivativeCache_->second;
        }
        auto d = evaluateDerivative(t);
        derivativeCache_ = std::make_pair(t, d);
        return d;
    }
    // All other methods forward to the inline version for std::function
#define CURVE_FORWARD(NAME, ...) \
    Curve<Scalar, Dim, std::function<PointType(Scalar)>> tmp(curveFunc_, domain_.first, domain_.second); \
    return tmp.NAME(__VA_ARGS__);
    template <typename Func>
    Scalar adaptiveParametricSolve(Func f, Scalar t0, Scalar t1, Scalar targetTol) const {
        CURVE_FORWARD(adaptiveParametricSolve, f, t0, t1, targetTol);
    }
    Eigen::Matrix<Scalar, Dim, 1> evaluateSecondDerivative(Scalar t) const {
        CURVE_FORWARD(evaluateSecondDerivative, t);
    }
    Eigen::Matrix<Scalar, Dim, 1> evaluateTangent(Scalar t) const {
        CURVE_FORWARD(evaluateTangent, t);
    }
    std::pair<Eigen::Matrix<Scalar, Dim, 1>, Scalar> evaluateDerivativeWithConfidence(Scalar t) const {
        CURVE_FORWARD(evaluateDerivativeWithConfidence, t);
    }
    Scalar evaluateCurvature(Scalar t) const {
        CURVE_FORWARD(evaluateCurvature, t);
    }
    Scalar evaluateIntegral() const {
        CURVE_FORWARD(evaluateIntegral);
    }
    std::pair<PointType, PointType> boundingBox() const {
        CURVE_FORWARD(boundingBox);
    }
    inline void invalidateBoundingBoxCache() const {
        bboxValid_ = false;
        bboxSamples_.clear();
    }
    std::pair<Curve, Curve> subDivide(Scalar tSplit) const {
        // Use legacy constructor
        const Scalar tMin = domain_.first;
        const Scalar tMax = domain_.second;
        const Scalar t0 = std::max(tMin, std::min(tSplit, tMax));
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
    std::pair<Scalar, Scalar> domain() const { return domain_; }
    template <typename Transform>
    Curve applyTransform(const Transform& T) const {
        auto newFunc = [T, this](Scalar t) -> PointType {
            return T.apply(this->evaluate(t));
        };
        return Curve(newFunc, domain_.first, domain_.second);
    }
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
    mutable bool bboxValid_ = false;
    mutable std::pair<PointType, PointType> cachedBox_;
    mutable std::vector<Eigen::Matrix<Scalar, Dim, 1>> bboxSamples_;
};


// ---- Bezier curve class ----
template <typename Scalar, int Dim = Eigen::Dynamic>
class Bezier {
public:
    using PointType = Point<Scalar, Dim>;
    using VecType = Eigen::Matrix<Scalar, Dim, 1>;
    using ScalarType = Scalar;
    static constexpr int Dimension = Dim;

    Bezier() = default;
    Bezier(const std::vector<PointType>& controlPoints)
        : controlPoints_(controlPoints) {}

    // Evaluate at parameter t in [0,1] using de Casteljau's algorithm
    PointType evaluate(Scalar t) const {
        if (controlPoints_.empty()) return PointType(VecType::Zero());
        std::vector<VecType> pts;
        pts.reserve(controlPoints_.size());
        for (const auto& p : controlPoints_) pts.push_back(p.coords);
        size_t n = pts.size();
        for (size_t k = 1; k < n; ++k) {
            for (size_t i = 0; i < n - k; ++i) {
                pts[i] = (Scalar(1) - t) * pts[i] + t * pts[i + 1];
            }
        }
        return PointType(pts[0]);
    }

    // Convert to a Curve object over [0,1]
    Curve<Scalar, Dim> toCurve() const {
        auto func = [*this](Scalar t) { return this->evaluate(t); };
        return Curve<Scalar, Dim>(func, Scalar(0), Scalar(1));
    }

    // Evaluate derivative at parameter t in [0,1]
    VecType evaluateDerivative(Scalar t) const {
        const std::size_t n = controlPoints_.size();
        if (n <= 1) {
            return VecType::Zero();
        }

        // Build derivative control points in-place (stack/local)
        std::vector<VecType> dctrl;
        dctrl.reserve(n - 1);
        const Scalar deg = static_cast<Scalar>(n - 1);
        for (std::size_t i = 0; i + 1 < n; ++i) {
            dctrl.push_back(deg * (controlPoints_[i + 1].coords - controlPoints_[i].coords));
        }

        // de Casteljau on derivative control points
        std::size_t m = dctrl.size();
        for (std::size_t k = 1; k < m; ++k) {
            for (std::size_t i = 0; i < m - k; ++i) {
                dctrl[i] = (Scalar(1) - t) * dctrl[i] + t * dctrl[i + 1];
            }
        }
        return dctrl[0];
    }

    // Second derivative: degree (n-2) Bezier built from first-derivative control points
    VecType evaluateSecondDerivative(Scalar t) const {
        const std::size_t n = controlPoints_.size();
        if (n <= 2) {
            return VecType::Zero();
        }
        const Scalar deg  = static_cast<Scalar>(n - 1);
        const Scalar deg2 = static_cast<Scalar>(n - 2);

        // First-derivative control points
        std::vector<VecType> d1;
        d1.reserve(n - 1);
        for (std::size_t i = 0; i + 1 < n; ++i) {
            d1.push_back(deg * (controlPoints_[i + 1].coords - controlPoints_[i].coords));
        }

        // Second-derivative control points from d1
        std::vector<VecType> d2;
        d2.reserve(n - 2);
        for (std::size_t i = 0; i + 1 < d1.size(); ++i) {
            d2.push_back(deg2 * (d1[i + 1] - d1[i]));
        }

        // de Casteljau on d2
        std::size_t m = d2.size();
        for (std::size_t k = 1; k < m; ++k) {
            for (std::size_t i = 0; i < m - k; ++i) {
                d2[i] = (Scalar(1) - t) * d2[i] + t * d2[i + 1];
            }
        }
        return d2[0];
    }

    // Access control points
    const std::vector<PointType>& controlPoints() const { return controlPoints_; }
    std::vector<PointType>& controlPoints() { return controlPoints_; }

private:
    std::vector<PointType> controlPoints_;
};

// ---- NURBS curve class ----
template <typename Scalar, int Dim = Eigen::Dynamic>
class NURBS {
public:
    using PointType = Point<Scalar, Dim>;
    using VecType = Eigen::Matrix<Scalar, Dim, 1>;
    using ScalarType = Scalar;
    static constexpr int Dimension = Dim;

    NURBS() = default;
    NURBS(const std::vector<PointType>& controlPoints,
          const std::vector<Scalar>& weights,
          const std::vector<Scalar>& knotVector,
          int degree)
        : controlPoints_(controlPoints),
          weights_(weights),
          knotVector_(knotVector),
          degree_(degree) {}

    // Overload: generate clamped uniform knot vector automatically
    NURBS(const std::vector<PointType>& controlPoints,
          const std::vector<Scalar>& weights,
          int degree)
        : controlPoints_(controlPoints),
          weights_(weights),
          degree_(degree)
    {
        const int n = static_cast<int>(controlPoints_.size());
        const int m = n + degree + 1;
        knotVector_.resize(m);
        for (int i = 0; i < m; ++i) {
            if (i <= degree) knotVector_[i] = Scalar(0);
            else if (i >= n) knotVector_[i] = Scalar(1);
            else knotVector_[i] = Scalar(i - degree) / Scalar(n - degree);
        }
    }

    // Evaluate at parameter t using rational Cox–de Boor
    PointType evaluate(Scalar t) const {
        int n = static_cast<int>(controlPoints_.size());
        if (n == 0 || weights_.size() != size_t(n) || knotVector_.size() < size_t(n + degree_ + 1))
            return PointType(VecType::Zero());

        // Find span
        int span = findKnotSpan(t);
        // Compute basis
        std::vector<Scalar> N = basisFunctions(span, t);
        VecType num = VecType::Zero();
        Scalar denom = Scalar(0);
        for (int i = 0; i <= degree_; ++i) {
            Scalar wN = weights_[span - degree_ + i] * N[i];
            num += wN * controlPoints_[span - degree_ + i].coords;
            denom += wN;
        }
        if (std::abs(denom) > Scalar(0)) {
            return PointType(num / denom);
        }
        return PointType(VecType::Zero());
    }

    // Convert to a Curve object over [knotVector_[degree_], knotVector_[n]]
    Curve<Scalar, Dim> toCurve() const {
        Scalar t0 = knotVector_[degree_];
        Scalar t1 = knotVector_[controlPoints_.size()];
        auto func = [*this](Scalar t) { return this->evaluate(t); };
        return Curve<Scalar, Dim>(func, t0, t1);
    }

    // Accessors
    const std::vector<PointType>& controlPoints() const { return controlPoints_; }
    const std::vector<Scalar>& weights() const { return weights_; }
    const std::vector<Scalar>& knotVector() const { return knotVector_; }
    int degree() const { return degree_; }

private:
    // Find span such that knotVector_[span] <= t < knotVector_[span+1]
    int findKnotSpan(Scalar t) const {
        int n = static_cast<int>(controlPoints_.size());
        if (t >= knotVector_[n]) return n - 1;
        if (t <= knotVector_[degree_]) return degree_;
        int low = degree_, high = n, mid = (low + high) / 2;
        while (t < knotVector_[mid] || t >= knotVector_[mid + 1]) {
            if (t < knotVector_[mid])
                high = mid;
            else
                low = mid;
            mid = (low + high) / 2;
        }
        return mid;
    }

    // Compute nonrational basis functions N[0..degree_] at t
    std::vector<Scalar> basisFunctions(int span, Scalar t) const {
        std::vector<Scalar> N(degree_ + 1, Scalar(0));
        std::vector<Scalar> left(degree_ + 1), right(degree_ + 1);
        N[0] = Scalar(1);
        for (int j = 1; j <= degree_; ++j) {
            left[j] = t - knotVector_[span + 1 - j];
            right[j] = knotVector_[span + j] - t;
            Scalar saved = Scalar(0);
            for (int r = 0; r < j; ++r) {
                Scalar denom = right[r + 1] + left[j - r];
                Scalar temp = N[r] / (denom != Scalar(0) ? denom : Scalar(1));
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }
        return N;
    }

    // Compute basis functions and their first derivatives at parameter t
    std::pair<std::vector<Scalar>, std::vector<Scalar>>
    basisFunctionsAndDerivatives(int span, Scalar t) const {
        const int p = degree_;
        std::vector<Scalar> N(p + 1, Scalar(0));
        std::vector<Scalar> dN(p + 1, Scalar(0));
        std::vector<Scalar> left(p + 1), right(p + 1);

        N[0] = Scalar(1);
        for (int j = 1; j <= p; ++j) {
            left[j] = t - knotVector_[span + 1 - j];
            right[j] = knotVector_[span + j] - t;
            Scalar saved = Scalar(0);
            for (int r = 0; r < j; ++r) {
                Scalar denom = right[r + 1] + left[j - r];
                Scalar temp  = N[r] / (denom != Scalar(0) ? denom : Scalar(1));
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }

        // Derivatives (first derivative only)
        for (int r = 0; r <= p; ++r) {
            Scalar dr = Scalar(0);

            if (r > 0) {
                Scalar denom = knotVector_[span + r] - knotVector_[span + r - p];
                if (std::abs(denom) > Scalar(0)) {
                    dr += (Scalar)p * N[r - 1] / denom;
                }
            }
            if (r < p) {
                Scalar denom = knotVector_[span + r + 1] - knotVector_[span + r + 1 - p];
                if (std::abs(denom) > Scalar(0)) {
                    dr -= (Scalar)p * N[r] / denom;
                }
            }
            dN[r] = dr;
        }

        return {N, dN};
    }

    // Evaluate first derivative C'(t)
    VecType evaluateDerivative(Scalar t) const {
        int n = static_cast<int>(controlPoints_.size());
        if (n == 0 ||
            weights_.size() != static_cast<std::size_t>(n) ||
            knotVector_.size() < static_cast<std::size_t>(n + degree_ + 1)) {
            return VecType::Zero();
        }

        int span = findKnotSpan(t);
        auto [N, dN] = basisFunctionsAndDerivatives(span, t);

        VecType Cw = VecType::Zero();   // weighted point sum
        VecType dCw = VecType::Zero();  // derivative of weighted point sum
        Scalar  wSum = Scalar(0);
        Scalar  dwSum = Scalar(0);

        for (int i = 0; i <= degree_; ++i) {
            int idx = span - degree_ + i;
            const Scalar w   = weights_[idx];
            const Scalar Ni  = N[i];
            const Scalar dNi = dN[i];

            Cw  += (w * Ni) * controlPoints_[idx].coords;
            dCw += (w * dNi) * controlPoints_[idx].coords;
            wSum  += w * Ni;
            dwSum += w * dNi;
        }

        if (std::abs(wSum) < Scalar(1e-12)) {
            return VecType::Zero();
        }

        VecType num = dCw * wSum - Cw * dwSum;
        Scalar  den = wSum * wSum;
        return num / den;
    }

    std::vector<PointType> controlPoints_;
    std::vector<Scalar> weights_;
    std::vector<Scalar> knotVector_;
    int degree_ = 0;
};

} // namespace Euclid::Geometry
