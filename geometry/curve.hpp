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

    // Compute local curvature at parameter t:
    // κ(t) = ||a_perp|| / ||v||^2, where a_perp = a - ((a·v)/||v||^2) v
    Scalar evaluateCurvature(Scalar t) const {
        using VecD = Eigen::Matrix<double, Dim, 1>;
        Tolerance tol;
        const Scalar t0 = domain_.first;
        const Scalar t1 = domain_.second;
        // Clamp t to domain
        Scalar tc = t;
        if (tc < t0) tc = t0;
        if (tc > t1) tc = t1;

        // --- Use same adaptive step logic as in evaluateDerivative() ---
        const Scalar epsMach = std::numeric_limits<Scalar>::epsilon();
        const Scalar span    = std::max(Scalar(1e-12), t1 - t0);
        auto safeClampToDomain = [&](Scalar tval) -> Scalar {
            if (tval <= t0) return t0;
            if (tval >= t1) return t1;
            return tval;
        };
        auto choose_h = [&](Scalar tval) -> Scalar {
            VecD pd = curveFunc_(tval).coords.template cast<double>();
            const double geomScale = pd.norm();
            const double base = std::max(1.0, geomScale + pd.cwiseAbs().sum());
            const Scalar hTol = tol.evaluateEpsilon(static_cast<Scalar>(base));
            const Scalar hUlp = std::max(epsMach * (Scalar(1) + std::abs(tval)), epsMach);
            const Scalar hDom = std::max(span * Scalar(1e-6), epsMach);
            const Scalar hCbrt = Scalar(std::cbrt(static_cast<double>(epsMach))) * std::max(span, Scalar(1));
            Scalar h0 = std::max(std::max(hTol, hUlp), std::max(hDom, hCbrt));
            const Scalar hMax = std::max(span * Scalar(1e-2), epsMach);
            h0 = std::min(h0, hMax);
            h0 = std::max(h0, hUlp * Scalar(8));
            const Scalar toLeft  = tval - t0;
            const Scalar toRight = t1 - tval;
            const Scalar hBound  = std::max(std::min(toLeft, toRight), epsMach);
            return std::min(h0, hBound);
        };
        // --- Compute both velocity and acceleration using a single finite-difference block ---
        Scalar h = choose_h(tc);
        Scalar tp = safeClampToDomain(tc + h);
        Scalar tm = safeClampToDomain(tc - h);
        double denom = std::max(static_cast<double>(tp - tm), static_cast<double>(epsMach));
        double h2 = std::max(static_cast<double>(h) * static_cast<double>(h), static_cast<double>(epsMach));
        VecD fp = curveFunc_(tp).coords.template cast<double>();
        VecD fm = curveFunc_(tm).coords.template cast<double>();
        VecD f0 = curveFunc_(tc).coords.template cast<double>();
        VecD v = (fp - fm) / denom;
        if (!v.allFinite()) return Scalar(0);
        double vnorm = v.norm();
        if (!(vnorm > 0.0)) return Scalar(0);
        VecD a = (fp - 2.0 * f0 + fm) / h2;
        if (!a.allFinite()) return Scalar(0);
        double v2 = std::max(vnorm * vnorm, static_cast<double>(epsMach));
        double proj = (a.dot(v)) / v2;
        VecD a_perp = a - proj * v;
        double numer = a_perp.norm();
        if (!std::isfinite(numer)) return Scalar(0);
        double kappa = numer / v2;
        return static_cast<Scalar>(kappa);
    }


    // Overload: Estimate arc length integral adaptively using internal tolerance model (Gauss–Kronrod)
    Scalar evaluateIntegral() const {
        // Use a default Tolerance model as in evaluateDerivative()
        Tolerance tol;
        Scalar t0 = domain_.first;
        Scalar t1 = domain_.second;
        constexpr int maxDepth = 10;
        if (!(t0 < t1)) return Scalar(0);

        auto speed = [&](Scalar t) -> Scalar {
            Scalar v = evaluateDerivative(t).norm();
            if (!std::isfinite(v)) v = Scalar(0);
            return v;
        };

        // --- Step 1: Periodicity / oscillation detection on speed(t) (correlation-based) ---
        const int M0 = 256;
        std::vector<Scalar> s(M0);
        for (int i = 0; i < M0; ++i) {
            Scalar t = t0 + (t1 - t0) * Scalar(i) / Scalar(M0 - 1);
            s[i] = speed(t);
        }

        // Basic stats
        const Scalar meanS = static_cast<Scalar>(std::accumulate(s.begin(), s.end(), static_cast<double>(0)) / static_cast<double>(M0));

        // Demean for correlation
        std::vector<Scalar> sd(M0);
        for (int i = 0; i < M0; ++i) sd[i] = s[i] - meanS;

        // Autocorrelation (naive O(N^2), small N)
        auto autocorr_at = [&](int lag) -> Scalar {
            // Corr(lag) = sum_i sd[i]*sd[i+lag]
            const int L = M0 - lag;
            if (L <= 1) return Scalar(0);
            Scalar num = Scalar(0);
            Scalar den = Scalar(0);
            for (int i = 0; i < L; ++i) {
                num += sd[i] * sd[i + lag];
                den += sd[i] * sd[i];
            }
            if (den == Scalar(0)) return Scalar(0);
            return num / den; // normalized (not strictly Pearson but enough as indicator)
        };

        // Search for a strong nonzero-lag peak up to 1/3 of the window
        Scalar bestR = Scalar(0);
        const int maxLag = M0 / 3;
        for (int lag = 1; lag <= maxLag; ++lag) {
            Scalar r = autocorr_at(lag);
            if (r > bestR) {
                bestR = r;
            }
        }

        // Zero-cross count (secondary cue) -- removed as unused

        // The periodicStrict and oscillatory flags are detected here for diagnostics,
        // but the integration path is always Gauss–Kronrod below.

        // Arc length integrand: velocity magnitude
        auto velocityAt = speed;
        // Gauss–Kronrod 15-point nodes (abscissae) and weights (on [-1,1])
        static constexpr double GK15_x[8] = {
            0.0000000000000000,
            0.2077849550078985,
            0.4058451513773972,
            0.5860872354676911,
            0.7415311855993945,
            0.8648644233597691,
            0.9491079123427585,
            0.9914553711208126
        };
        static constexpr double GK15_w[8] = {
            0.2094821410847278,
            0.2044329400752989,
            0.1903505780647854,
            0.1690047266392679,
            0.1406532597155259,
            0.1047900103222502,
            0.0630920926299786,
            0.0229353220105292
        };
        static constexpr double G7_x[4] = {
            0.0000000000000000,
            0.4058451513773972,
            0.7415311855993945,
            0.9491079123427585
        };
        static constexpr double G7_w[4] = {
            0.4179591836734694,
            0.3818300505051189,
            0.2797053914892766,
            0.1294849661688697
        };
        // Adaptive recursive Gauss–Kronrod quadrature on [a,b]
        std::function<Scalar(Scalar, Scalar, int)> integrateGK;
        integrateGK = [&](Scalar a, Scalar b, int depth) -> Scalar {
            const Scalar mid  = (a + b) * Scalar(0.5);
            const Scalar half = (b - a) * Scalar(0.5);
            const Scalar width = std::abs(b - a);

            // Only fm is needed for quadrature; fa/fb are not used except in osc (diagnostic, but can be omitted)
            const Scalar fm = velocityAt(mid);

            // --- Gauss–Kronrod 15-point integral for arc length ---
            Scalar I15 = 0;
            {
                I15 += Scalar(GK15_w[0]) * fm; // center once
                for (int i = 1; i < 8; ++i) {
                    const Scalar xi = Scalar(GK15_x[i]);
                    const Scalar fp = velocityAt(mid + half * xi);
                    const Scalar fn = velocityAt(mid - half * xi);
                    I15 += Scalar(GK15_w[i]) * (fp + fn);
                }
                I15 *= half;
            }

            // --- Embedded 7-point Gauss subset ---
            Scalar I7 = 0;
            {
                I7 += Scalar(G7_w[0]) * fm; // center once
                for (int i = 1; i < 4; ++i) {
                    const Scalar xi = Scalar(G7_x[i]);
                    const Scalar fp = velocityAt(mid + half * xi);
                    const Scalar fn = velocityAt(mid - half * xi);
                    I7 += Scalar(G7_w[i]) * (fp + fn);
                }
                I7 *= half;
            }

            // Error and tolerance for this interval
            const Scalar err     = std::abs(I15 - I7);
            // tolHere is only used once, so can be omitted

            // Decide whether to split: only if error large OR significant normalized oscillation
            const bool needSplit = (err >= tol.paramTol * width);

            // Hard stops: recursion depth and minimal width relative to total span
            if (!needSplit || depth >= maxDepth || width < (t1 - t0) * Scalar(1e-6)) {
                return I15;
            }

            const Scalar left  = integrateGK(a,  mid, depth + 1);
            const Scalar right = integrateGK(mid, b,  depth + 1);
            return left + right;
        };
        return integrateGK(t0, t1, 0);
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
