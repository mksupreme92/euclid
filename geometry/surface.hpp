#pragma once
#include <array>
#include <functional>
#include <optional>
#include <cmath>
#include <numeric>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <Eigen/Dense>
#include <limits>
#include "point.hpp"

namespace Euclid::Geometry {

template<typename T>
struct NumericalLimits {
    T eps_param;
    T eps_det;
    T eps_cross;

    static NumericalLimits<T> derive(
        const Euclid::Tolerance& tol, T J = T(1))
    {
        T floor = std::max<T>(tol.paramTol, std::numeric_limits<T>::epsilon());
        T eps_param = std::max<T>(tol.paramTol, tol.evaluateEpsilon(std::max<T>(J, floor)));
        T eps_det   = std::max<T>(tol.evaluateEpsilon(std::max<T>(J, floor)), T(1e-30));
        T eps_cross = std::max<T>(tol.evaluateEpsilon(std::max<T>(J, floor)), T(1e-20));
        return {eps_param, eps_det, eps_cross};
    }
};

template<typename T, int N>
class Surface {
public:
    using PointType = Point<T, N>;

    // Parametric function: S(u,v) -> ℝⁿ
    std::function<PointType(T, T)> surfaceFunc;

    // Optional analytic SDF: SDF(P) -> ℝ
    std::function<T(const Eigen::Matrix<T, N, 1>&)> analyticSDF = nullptr;

    // Domain: ((umin, vmin), (umax, vmax))
    std::pair<std::pair<T,T>, std::pair<T,T>> surfaceDomain;

    struct SurfacePatchAABB {
        T umin, umax, vmin, vmax;
        PointType minPoint;
        PointType maxPoint;
    };

    // Holds the result of projecting a point onto the surface
    struct PointProjection {
        bool ok{false};
        T u{T(0)}, v{T(0)};
        T signed_distance{std::numeric_limits<T>::infinity()};
        PointType surface_point{};
        int iterations{0};
    };
    // Gather a small set of seed (u,v) guesses from the BVH near a point P
    std::vector<std::pair<T,T>> gatherSeedsFromBVHNearPoint(const Eigen::Matrix<T, N, 1>& P,
                                                           const Euclid::Tolerance& tol) const
    {
        std::vector<std::pair<T,T>> seeds;
        if (!bvh_valid) generateBVH(40, tol);
        const auto& patches = bvh_patches;
        // Compute adaptive max_seeds based on tolerance
        T surface_domain_scale_est = std::max<T>(
            std::abs(surfaceDomain.first.second - surfaceDomain.first.first),
            std::abs(surfaceDomain.second.second - surfaceDomain.second.first)
        );
        auto limits = NumericalLimits<T>::derive(tol, surface_domain_scale_est);
        int max_seeds = std::clamp<int>(
            static_cast<int>(6 / std::max<T>(tol.paramTol, limits.eps_param)),
            4, 32
        );
        // Compute distances from P to each patch AABB and keep a few best
        struct Item { T d2; T umin, umax, vmin, vmax; };
        std::vector<Item> items;
        items.reserve(patches.size());
        for (const auto& patch : patches) {
            Eigen::Matrix<T,N,1> q;
            for (int k=0;k<N;++k) {
                T a = patch.minPoint.coords(k);
                T b = patch.maxPoint.coords(k);
                T x = P(k);
                if (x < a) q(k) = a;
                else if (x > b) q(k) = b;
                else q(k) = x;
            }
            T d2 = (q - P).squaredNorm();
            items.push_back({d2, patch.umin, patch.umax, patch.vmin, patch.vmax});
        }
        std::nth_element(items.begin(), items.begin() + std::min<size_t>(items.size(), (size_t)max_seeds), items.end(),
                         [](const Item& A, const Item& B){ return A.d2 < B.d2; });
        size_t take = std::min<size_t>(items.size(), (size_t)max_seeds);
        seeds.reserve(take);
        for (size_t i=0;i<take;++i) {
            const auto& it = items[i];
            seeds.emplace_back((it.umin + it.umax)/T(2), (it.vmin + it.vmax)/T(2));
        }
        // Always include domain center as a fallback
        seeds.emplace_back((surfaceDomain.first.first + surfaceDomain.first.second)/T(2),
                           (surfaceDomain.second.first + surfaceDomain.second.second)/T(2));
        // Add tolerance-aware fallback seeds along domain boundaries if extent is large
        if ((surfaceDomain.first.second - surfaceDomain.first.first) > tol.paramTol * 10)
            seeds.emplace_back(surfaceDomain.first.first, (surfaceDomain.second.first + surfaceDomain.second.second)/2);
        if ((surfaceDomain.second.second - surfaceDomain.second.first) > tol.paramTol * 10)
            seeds.emplace_back((surfaceDomain.first.first + surfaceDomain.first.second)/2, surfaceDomain.second.first);
        return seeds;
    }

private:
    mutable std::vector<SurfacePatchAABB> bvh_patches;
    mutable bool bvh_valid = false;
    // Optional direction hint for directional BVH refinement
    mutable std::optional<Eigen::Matrix<T, N, 1>> bvh_dir_hint;
    // Cache of last successful projection seed (u,v)
    mutable std::optional<std::pair<T,T>> last_seed_uv;
    // Cache for curvature gradient magnitudes at (u,v)
    mutable std::map<std::pair<T,T>, T> curvature_cache;

    // --- Curvature proxy and gradient helpers (tolerance-aware) ---
    // Dimension-agnostic "curvature" proxy: sine of angle between Su and Sv (area of parallelogram normalized)
    T curvatureProxy(T u, T v, const Euclid::Tolerance& tol) const {
        auto [Su, Sv] = evaluatePartials(u, v);
        T su2 = Su.squaredNorm();
        T sv2 = Sv.squaredNorm();
        T dot = Su.dot(Sv);
        T area2 = std::max<T>(0, su2 * sv2 - dot * dot); // = ||Su x Sv||^2 in 3D; Gram determinant in ND
        T denom = std::sqrt(std::max<T>(su2, tol.evaluateEpsilon(su2))) *
                  std::sqrt(std::max<T>(sv2, tol.evaluateEpsilon(sv2)));
        if (denom <= T(0)) return T(0);
        return std::sqrt(area2) / denom; // in [0,1]
    }

    // Approximate |∇kappa| via central differences in param-space, with conservative monotone-max caching
    T curvatureGradientMagnitude(T u, T v, T du, T dv, const Euclid::Tolerance& tol) const {
        // Conservative caching: store the MAX observed gradient at quantized (u,v).
        // Quantize keys to reduce cache misses from near-identical parameter pairs.
        auto limits_q = NumericalLimits<T>::derive(tol, std::max<T>(tol.paramTol, std::max(du, dv)));
        T qstep = std::max<T>(tol.paramTol, limits_q.eps_param);
        T uq = std::round(u / qstep) * qstep;
        T vq = std::round(v / qstep) * qstep;
        auto key = std::make_pair(uq, vq);
        Euclid::Tolerance tlocal = tol;
        const T epsu = std::max<T>(std::abs(du), tlocal.paramTol);
        const T epsv = std::max<T>(std::abs(dv), tlocal.paramTol);

        T ku_p = curvatureProxy(u + epsu, v, tol);
        T ku_m = curvatureProxy(u - epsu, v, tol);
        T kv_p = curvatureProxy(u, v + epsv, tol);
        T kv_m = curvatureProxy(u, v - epsv, tol);

        T gu = std::abs(ku_p - ku_m) / (T(2) * epsu);
        T gv = std::abs(kv_p - kv_m) / (T(2) * epsv);
        T computed = gu + gv; // L1 blend; robust and cheap
        // Optionally, smooth conservative value to ensure not rounded to zero prematurely
        T scale_ref = std::max<T>(computed, tol.paramTol);
        computed = std::max<T>(computed, tol.evaluateEpsilon(scale_ref));

        auto it = curvature_cache.find(key);
        if (it != curvature_cache.end()) {
            if (computed < it->second) {
                // Return the larger cached value to remain conservative
                return it->second;
            } else {
                // Tighten the cache to the higher (more conservative) value
                it->second = computed;
                return computed;
            }
        } else {
            curvature_cache[key] = computed;
            return computed;
        }
    }

public:
    Surface(const std::function<PointType(T,T)>& f,
            const std::pair<std::pair<T,T>, std::pair<T,T>>& domain,
            std::function<T(const Eigen::Matrix<T, N, 1>&)> sdf = nullptr)
        : surfaceFunc(f), analyticSDF(sdf), surfaceDomain(domain) {}

    // Optional: provide a direction hint to bias BVH refinement along that direction
    void setBVHDirectionHint(const Eigen::Matrix<T, N, 1>& dir_hint) const {
        if (dir_hint.norm() > T(0)) {
            bvh_dir_hint = dir_hint.normalized();
            bvh_valid = false; // force rebuild on next getBVH()
        }
    }

    // Clear any previously set BVH direction hint
    void clearBVHDirectionHint() const {
        bvh_dir_hint.reset();
        bvh_valid = false;
    }

    // Evaluate at (u,v)
    PointType evaluate(T u, T v) const {
        return surfaceFunc(u, v);
    }

    // Evaluate (Su, Sv) partial derivatives at (u,v) with default tolerance
    std::pair<Eigen::Matrix<T, N, 1>, Eigen::Matrix<T, N, 1>> evaluatePartials(T u, T v) const {
        Euclid::Tolerance tol;
        return evaluatePartials(u, v, tol);
    }

    // Evaluate (Su, Sv) partial derivatives at (u,v) with provided tolerance
    std::pair<Eigen::Matrix<T, N, 1>, Eigen::Matrix<T, N, 1>> evaluatePartials(T u, T v, const Euclid::Tolerance& tol) const {
        if constexpr (std::is_member_function_pointer_v<decltype(&Surface<T, N>::analyticalPartials)>) {
            // If analyticalPartials is overridden, use it
            return analyticalPartials(u, v);
        } else {
            // Fallback to numerical
            return numericalPartials(u, v, tol);
        }
    }

    // Analytical partial derivatives (override in derived classes for analytical surfaces)
    virtual std::pair<Eigen::Matrix<T, N, 1>, Eigen::Matrix<T, N, 1>> analyticalPartials(T u, T v) const {
        // Default: not implemented, fallback to numerical
        Euclid::Tolerance tol;
        return numericalPartials(u, v, tol);
    }

    // Numerical partial derivatives using central finite differences
    std::pair<Eigen::Matrix<T, N, 1>, Eigen::Matrix<T, N, 1>> numericalPartials(T u, T v, const Euclid::Tolerance& tol) const {
        // Adaptive step size based on tolerance and domain size
        T umin = surfaceDomain.first.first;
        T umax = surfaceDomain.first.second;
        T vmin = surfaceDomain.second.first;
        T vmax = surfaceDomain.second.second;
        T u_extent = std::abs(umax - umin);
        T v_extent = std::abs(vmax - vmin);
        // Step size is tol.paramTol * param extent, but not smaller than tol.paramTol itself
        T du = std::max(tol.paramTol * u_extent, tol.paramTol);
        T dv = std::max(tol.paramTol * v_extent, tol.paramTol);
        // Clamp du/dv to not exceed a tolerance-driven fraction of the domain, for robustness
        auto limits = NumericalLimits<T>::derive(tol, std::max(u_extent, v_extent));
        // derive a maximum fraction dynamically from tolerance and numeric precision
        T max_frac = std::max<T>(
            tol.paramTol,
            limits.eps_param * T(10)
        );
        du = std::min(du, u_extent * max_frac);
        dv = std::min(dv, v_extent * max_frac);
        // Central finite difference for each partial
        auto Su = (evaluate(u + du, v).coords - evaluate(u - du, v).coords) / (T(2) * du);
        auto Sv = (evaluate(u, v + dv).coords - evaluate(u, v - dv).coords) / (T(2) * dv);
        return {Su, Sv};
    }

    // Evaluate surface normal at (u,v) (only for N==3), with provided tolerance
    PointType evaluateNormal(T u, T v, const Euclid::Tolerance& tol) const {
        static_assert(N == 3, "evaluateNormal only makes sense for 3D surfaces");
        auto [Su, Sv] = evaluatePartials(u, v, tol);
        auto nvec = Su.cross(Sv);
        // Use tolerance to check for degeneracy
        if (nvec.norm() < tol.evaluateEpsilon(nvec.norm())) {
            // Degenerate, return PointType(Eigen::Matrix<T, N, 1>::Zero());
            return PointType(Eigen::Matrix<T, N, 1>::Zero());
        }
        return PointType(nvec.normalized());
    }

    // Generalized principal normal in R^N using Gram–Schmidt and optional reference direction
    Eigen::Matrix<T, N, 1> principalNormal(T u, T v,
                                           const Euclid::Tolerance& tol,
                                           const std::optional<Eigen::Matrix<T, N, 1>>& ref_dir_opt = std::nullopt) const
    {
        // Evaluate first-order partials with tolerance to keep step sizes consistent
        auto [Su_raw, Sv_raw] = evaluatePartials(u, v, tol);

        // Local derivative scales
        T su = Su_raw.norm();
        T sv = Sv_raw.norm();
        T floor = tol.evaluateEpsilon(std::max({su, sv, tol.paramTol}));
        T scale = std::max({su, sv, floor});

        // Tolerance-scaled epsilons driven by local conditioning
        T eps_dir   = tol.evaluateEpsilon(scale);                    // for vector magnitudes
        T eps_cross = tol.evaluateEpsilon(std::max<T>(su * sv, tol.paramTol));  // for area / cross terms

        // Guard against degeneracy
        if (su <= eps_dir && sv <= eps_dir) {
            // Fallback: pick a canonical axis or ref_dir if provided
            Eigen::Matrix<T,N,1> n = ref_dir_opt ? *ref_dir_opt : Eigen::Matrix<T,N,1>::UnitX();
            if (n.norm() == T(0)) n = Eigen::Matrix<T,N,1>::UnitX();
            return n.normalized();
        }

        // Gram–Schmidt to build an orthonormal basis of the tangent plane span{Su,Sv}
        Eigen::Matrix<T,N,1> e1, e2;
        e1 = (su > eps_dir) ? (Su_raw / su) : Su_raw;
        // remove component of Sv along e1
        Eigen::Matrix<T,N,1> v2 = Sv_raw - (Sv_raw.dot(e1)) * e1;
        T v2n = v2.norm();
        if (v2n <= eps_dir) {
            // Su and Sv nearly colinear; try a numerical perturbation by swapping roles
            // or create an arbitrary vector not parallel to e1
            Eigen::Matrix<T,N,1> arbit = Eigen::Matrix<T,N,1>::Zero();
            // Choose the axis with smallest alignment with e1
            int best_k = 0;
            T best_val = std::abs(e1(0));
            for (int k=1; k<N; ++k) {
                T val = std::abs(e1(k));
                if (val < best_val) { best_val = val; best_k = k; }
            }
            arbit(best_k) = tol.paramTol;
            v2 = arbit - (arbit.dot(e1)) * e1;
            v2n = v2.norm();
            if (v2n <= eps_dir) v2 = Eigen::Matrix<T,N,1>::UnitY(); // final fallback where available
            v2n = v2.norm();
        }
        e2 = (v2n > eps_dir) ? (v2 / v2n) : v2;

        if constexpr (N == 3) {
            auto n3 = e1.cross(e2);
            T nn = n3.norm();
            if (nn <= eps_cross) {
                // fallback to ref_dir or canonical axis
                Eigen::Matrix<T,3,1> n = ref_dir_opt ? *ref_dir_opt : Eigen::Matrix<T,3,1>::UnitZ();
                return n.normalized();
            }
            n3 /= nn;
            if (ref_dir_opt) {
                if (n3.dot(*ref_dir_opt) < T(0)) n3 = -n3;
            }
            return n3;
        } else {
            // Project a preferred direction into the orthogonal complement of span{e1,e2}
            Eigen::Matrix<T,N,1> pref = ref_dir_opt ? *ref_dir_opt : Eigen::Matrix<T,N,1>::UnitX();
            if (pref.norm() <= T(0)) pref = Eigen::Matrix<T,N,1>::UnitX();

            // Remove tangent components
            Eigen::Matrix<T,N,1> n = pref - (pref.dot(e1))*e1 - (pref.dot(e2))*e2;
            T nn = n.norm();
            if (nn <= eps_dir) {
                // Try canonical axes to find a stable orthogonal direction
                T best_norm = T(0);
                Eigen::Matrix<T,N,1> best = Eigen::Matrix<T,N,1>::Zero();
                for (int k=0; k<N; ++k) {
                    Eigen::Matrix<T,N,1> axis = Eigen::Matrix<T,N,1>::Zero();
                    axis(k) = tol.paramTol;
                    Eigen::Matrix<T,N,1> cand = axis - (axis.dot(e1))*e1 - (axis.dot(e2))*e2;
                    T cn = cand.norm();
                    if (cn > best_norm) { best_norm = cn; best = cand; }
                }
                if (best_norm > eps_dir) n = best; else n = pref; // ultimate fallback
            }
            n.normalize();
            // Align with preferred direction if provided
            if (ref_dir_opt && n.dot(*ref_dir_opt) < T(0)) n = -n;
            return n;
        }
    }

    // Project a point P onto the surface starting from a seed (u0, v0) using Gauss–Newton
    PointProjection projectPointFromSeed(const Eigen::Matrix<T, N, 1>& P,
                                         T u0, T v0,
                                         const Euclid::Tolerance& tol,
                                         int max_iters = 20) const
    {
        PointProjection r;
        r.u = u0; r.v = v0;
        T prev_norm = std::numeric_limits<T>::infinity();
        for (int it=0; it<max_iters; ++it) {
            auto S = evaluate(r.u, r.v).coords;
            auto R = S - P; // residual
            T rn = R.norm();
            T ref_scale = std::max<T>(rn, tol.paramTol);
            if (rn <= tol.evaluateEpsilon(ref_scale)) {
                r.ok = true; r.signed_distance = rn;
                r.surface_point = PointType(S);
                r.iterations = it;
                // Signed distance via principal normal aligned with any BVH hint
                Eigen::Matrix<T,N,1> n = principalNormal(r.u, r.v, tol, bvh_dir_hint);
                if (n.norm() > tol.evaluateEpsilon(std::max<T>(n.norm(), tol.paramTol))) {
                    r.signed_distance = (R.dot(n));
                } else {
                    r.signed_distance = rn;
                }
                return r;
            }
            // Build J (N x 2)
            auto [Su, Sv] = evaluatePartials(r.u, r.v);
            Eigen::Matrix<T,2,2> JTJ;
            JTJ(0,0) = Su.dot(Su);
            JTJ(0,1) = Su.dot(Sv);
            JTJ(1,0) = JTJ(0,1);
            JTJ(1,1) = Sv.dot(Sv);
            Eigen::Matrix<T,2,1> JTr;
            JTr(0) = Su.dot(R);
            JTr(1) = Sv.dot(R);
            // Solve JTJ * delta = -JTr (with tolerance-scaled Tikhonov for robustness)
            const T scale_J = std::max<T>(JTJ(0,0) + JTJ(1,1), tol.paramTol);
            T lam = std::max<T>(tol.evaluateEpsilon(scale_J), T(1e-14));
            JTJ(0,0) += lam;
            JTJ(1,1) += lam;

            // Determinant guard using tolerance-scaled epsilon
            const T det_eps = std::max<T>(tol.evaluateEpsilon(scale_J), T(1e-18));
            T det = JTJ(0,0)*JTJ(1,1) - JTJ(0,1)*JTJ(1,0);
            if (std::abs(det) < std::max<T>(det_eps, tol.evaluateEpsilon(scale_J))) break;
            Eigen::Matrix<T,2,1> delta;
            delta(0) = (-JTr(0)*JTJ(1,1) + JTr(1)*JTJ(0,1)) / det;
            delta(1) = (-JTr(1)*JTJ(0,0) + JTr(0)*JTJ(1,0)) / det;
            // Damped step with tolerance-driven adaptive scaling
            T step_scale = T(1);
            T domain_extent = std::max<T>(
                std::abs(surfaceDomain.first.second - surfaceDomain.first.first),
                std::abs(surfaceDomain.second.second - surfaceDomain.second.first)
            );
            T step_min = std::max<T>(
                tol.paramTol * T(0.25),
                tol.evaluateEpsilon(domain_extent)
            );
            int max_backtrack = std::clamp<int>(int(std::ceil(-std::log2(double(std::max<T>(tol.paramTol * T(10), T(1e-16)))))), 3, 10);

            for (int bs = 0; bs < max_backtrack; ++bs) {
                T nu = r.u - step_scale * delta(0);
                T nv = r.v - step_scale * delta(1);
                auto Stry = evaluate(nu, nv).coords;
                T rn_try = (Stry - P).norm();

                // Accept if residual decreases or the change is within tolerance noise (domain-scaled)
                T domain_extent = std::max<T>(
                    std::abs(surfaceDomain.first.second - surfaceDomain.first.first),
                    std::abs(surfaceDomain.second.second - surfaceDomain.second.first)
                );
                if (rn_try < rn || std::abs(rn_try - rn) < tol.evaluateEpsilon(std::max<T>(rn, domain_extent))) {
                    r.u = nu;
                    r.v = nv;
                    rn = rn_try;
                    break;
                }
                step_scale *= T(0.5);
                if (step_scale < step_min) {
                    // Prevent stalling below parametric tolerance
                    break;
                }
            }
            // Convergence checks
            T ref_scale2 = std::max<T>(rn, tol.paramTol);
            if (std::abs(prev_norm - rn) <= tol.evaluateEpsilon(ref_scale2)) {
                prev_norm = rn; r.iterations = it+1;
                r.ok = true; r.signed_distance = rn;
                r.surface_point = PointType(evaluate(r.u, r.v).coords);
                Eigen::Matrix<T,N,1> n = principalNormal(r.u, r.v, tol, bvh_dir_hint);
                if (n.norm() > tol.evaluateEpsilon(std::max<T>(n.norm(), tol.paramTol))) {
                    r.signed_distance = ((evaluate(r.u, r.v).coords - P).dot(n));
                } else {
                    r.signed_distance = rn;
                }
                return r;
            }
            prev_norm = rn;
        }
        // Finalize (may be not-ok)
        r.surface_point = PointType(evaluate(r.u, r.v).coords);
        {
            Eigen::Matrix<T,N,1> n = principalNormal(r.u, r.v, tol, bvh_dir_hint);
            if (n.norm() > T(0)) {
                r.signed_distance = ((r.surface_point.coords - P).dot(n));
            } else {
                r.signed_distance = (r.surface_point.coords - P).norm();
            }
        }
        return r;
    }

    // Try multiple BVH-based seeds and return the best projection
    PointProjection projectPoint(const Eigen::Matrix<T, N, 1>& P,
                                 const Euclid::Tolerance& tol,
                                 int max_iters = 20) const
    {
        // Ensure BVH exists for seeding
        if (!bvh_valid) {
            generateBVH(40, tol);
        }

        PointProjection best; // best so far

        // 1) Warm-start from last successful (u,v) if available
        if (last_seed_uv.has_value()) {
            auto [u0, v0] = *last_seed_uv;
            auto r = projectPointFromSeed(P, u0, v0, tol, max_iters);
            if (r.ok) {
                best = r;
            }
        }

        // 2) Gather BVH seeds and try them (skip if they are very close to the warm-start to avoid duplicates)
        auto seeds = gatherSeedsFromBVHNearPoint(P, tol);
        for (auto [u0, v0] : seeds) {
            // Skip near-duplicate seeds relative to the warm-start using domain-scaled tolerance
            T domain_extent = std::max<T>(
                std::abs(surfaceDomain.first.second - surfaceDomain.first.first),
                std::abs(surfaceDomain.second.second - surfaceDomain.second.first)
            );
            if (best.ok && std::abs(u0 - best.u) < tol.evaluateEpsilon(domain_extent) &&
                           std::abs(v0 - best.v) < tol.evaluateEpsilon(domain_extent)) {
                continue;
            }
            auto r = projectPointFromSeed(P, u0, v0, tol, max_iters);
            if (!best.ok || std::abs(r.signed_distance) < std::abs(best.signed_distance)) {
                best = r;
            }
        }

        // 3) Update warm-start cache
        if (best.ok) {
            last_seed_uv = std::make_pair(best.u, best.v);
        }

        return best;
    }

    // Signed distance convenience wrapper (returns signed in 3D, unsigned otherwise)
    T signedDistanceToSurface(const Eigen::Matrix<T, N, 1>& P,
                              const Euclid::Tolerance& tol,
                              T* out_u = nullptr,
                              T* out_v = nullptr,
                              int max_iters = 20) const
    {
        if (analyticSDF) {
            // Use analytic SDF if provided
            return analyticSDF(P);
        }
        auto r = projectPoint(P, tol, max_iters);
        if (out_u) *out_u = r.u;
        if (out_v) *out_v = r.v;

        // Use principal normal with optional global sign hint for consistent signed distances
        Eigen::Matrix<T,N,1> n = principalNormal(r.u, r.v, tol, bvh_dir_hint);
        // Enforce normal alignment with BVH direction hint if present
        if (bvh_dir_hint && n.dot(*bvh_dir_hint) < T(0)) {
            n = -n;
        }
        if (n.norm() > tol.evaluateEpsilon(std::max<T>(n.norm(), tol.paramTol))) {
            return (P - r.surface_point.coords).dot(n);
        } else {
            // Fallback to unsigned distance if normal is degenerate
            return (P - r.surface_point.coords).norm();
        }
    }

    // Apply a transformation to the surface
    template<typename Transform>
    Surface<T,N> applyTransform(const Transform& transform) const {
        auto newFunc = [transform, f = surfaceFunc](T uVal, T vVal) -> PointType {
            return transform.apply(f(uVal, vVal));
        };
        return Surface<T,N>(newFunc, surfaceDomain);
    }



    // Helper: Adaptive recursive patch subdivision for BVH (tolerance-aware, with optional direction hint)
    void subdividePatch(T umin, T umax, T vmin, T vmax,
                        const Euclid::Tolerance& tol,
                        T min_param_extent, int max_depth, int depth,
                        int budget,
                        bool have_dir, const Eigen::Matrix<T, N, 1>& dir_hint) const
    {
        T uc = (umin + umax) / T(2);
        T vc = (vmin + vmax) / T(2);
        auto [Su, Sv] = evaluatePartials(uc, vc);
        T J = std::max(Su.norm(), Sv.norm());
        auto limits = NumericalLimits<T>::derive(tol, J);
        T du = std::abs(umax - umin);
        T dv = std::abs(vmax - vmin);
        T surface_extent = std::max<T>(du, dv);
        T local_scale = std::max<T>(surface_extent * J, tol.paramTol);
        T eps_local = tol.evaluateEpsilon(local_scale);

        // Per-branch work budget: stop subdividing if exhausted
        if (budget <= 0) {
            // Adaptive sampling density based on curvature and tolerance
            int base_samples = 5;
            T G_center = curvatureGradientMagnitude(uc, vc, std::max<T>(du*0.5, tol.paramTol), std::max<T>(dv*0.5, tol.paramTol), tol);
            T curvature_scale = std::clamp<T>(1 + G_center * (1 / std::sqrt(std::max<T>(tol.paramTol, limits.eps_param))), 1, 4);
            T tol_scale = std::clamp<T>(1 / std::sqrt(std::max<T>(tol.paramTol, limits.eps_param)), 1, 4);
            int adaptive_samples = std::clamp<int>(int(base_samples * 0.5 * (curvature_scale + tol_scale)), 5, 25);

            std::vector<PointType> samples;
            samples.reserve(adaptive_samples);
            int grid_n = std::max<int>(2, std::round(std::sqrt(adaptive_samples)));
            for (int i = 0; i < grid_n; ++i) {
                for (int j = 0; j < grid_n; ++j) {
                    T u_sample = umin + (umax - umin) * (T(i) / T(grid_n - 1));
                    T v_sample = vmin + (vmax - vmin) * (T(j) / T(grid_n - 1));
                    samples.push_back(evaluate(u_sample, v_sample));
                }
            }
            PointType minPt = samples[0];
            PointType maxPt = samples[0];
            for (const auto& p : samples) {
                for (int k = 0; k < N; ++k) {
                    if (p.coords(k) < minPt.coords(k)) minPt.coords(k) = p.coords(k);
                    if (p.coords(k) > maxPt.coords(k)) maxPt.coords(k) = p.coords(k);
                }
            }
            bvh_patches.push_back(SurfacePatchAABB{umin, umax, vmin, vmax, minPt, maxPt});
            return;
        }

        // Global cap to prevent exponential explosion in highly oscillatory regions
        size_t MAX_PATCHES = static_cast<size_t>(
            std::clamp<T>(
                T(4096) / std::max<T>(tol.paramTol, limits.eps_param),
                T(1024), T(32768)
            ));
        if (bvh_patches.size() >= MAX_PATCHES) {
            int base_samples = 5;
            T G_center = curvatureGradientMagnitude(uc, vc, std::max<T>(du*0.5, tol.paramTol), std::max<T>(dv*0.5, tol.paramTol), tol);
            T curvature_scale = std::clamp<T>(1 + G_center * (1 / std::sqrt(std::max<T>(tol.paramTol, limits.eps_param))), 1, 4);
            T tol_scale = std::clamp<T>(1 / std::sqrt(std::max<T>(tol.paramTol, limits.eps_param)), 1, 4);
            int adaptive_samples = std::clamp<int>(int(base_samples * 0.5 * (curvature_scale + tol_scale)), 5, 25);

            std::vector<PointType> samples;
            samples.reserve(adaptive_samples);
            int grid_n = std::max<int>(2, std::round(std::sqrt(adaptive_samples)));
            for (int i = 0; i < grid_n; ++i) {
                for (int j = 0; j < grid_n; ++j) {
                    T u_sample = umin + (umax - umin) * (T(i) / T(grid_n - 1));
                    T v_sample = vmin + (vmax - vmin) * (T(j) / T(grid_n - 1));
                    samples.push_back(evaluate(u_sample, v_sample));
                }
            }
            PointType minPt = samples[0];
            PointType maxPt = samples[0];
            for (const auto& p : samples) {
                for (int k = 0; k < N; ++k) {
                    if (p.coords(k) < minPt.coords(k)) minPt.coords(k) = p.coords(k);
                    if (p.coords(k) > maxPt.coords(k)) maxPt.coords(k) = p.coords(k);
                }
            }
            bvh_patches.push_back(SurfacePatchAABB{umin, umax, vmin, vmax, minPt, maxPt});
            return;
        }

        int base_samples = 5;
        T G_center = curvatureGradientMagnitude(uc, vc, std::max<T>(du*0.5, tol.paramTol), std::max<T>(dv*0.5, tol.paramTol), tol);
        T curvature_scale = std::clamp<T>(1 + G_center * (1 / std::sqrt(std::max<T>(tol.paramTol, limits.eps_param))), 1, 4);
        T tol_scale = std::clamp<T>(1 / std::sqrt(std::max<T>(tol.paramTol, limits.eps_param)), 1, 4);
        int adaptive_samples = std::clamp<int>(int(base_samples * 0.5 * (curvature_scale + tol_scale)), 5, 25);
        std::vector<PointType> samples;
        samples.reserve(adaptive_samples);
        int grid_n = std::max<int>(2, std::round(std::sqrt(adaptive_samples)));
        for (int i = 0; i < grid_n; ++i) {
            for (int j = 0; j < grid_n; ++j) {
                T u_sample = umin + (umax - umin) * (T(i) / T(grid_n - 1));
                T v_sample = vmin + (vmax - vmin) * (T(j) / T(grid_n - 1));
                samples.push_back(evaluate(u_sample, v_sample));
            }
        }
        PointType minPt = samples[0];
        PointType maxPt = samples[0];
        for (const auto& p : samples) {
            for (int k = 0; k < N; ++k) {
                if (p.coords(k) < minPt.coords(k)) minPt.coords(k) = p.coords(k);
                if (p.coords(k) > maxPt.coords(k)) maxPt.coords(k) = p.coords(k);
            }
        }
        T bbox_diag = (maxPt.coords - minPt.coords).norm();

        // Param-step for gradient probing: half patch, but not smaller than tolerance-induced param eps
        T param_eps = std::max(min_param_extent * T(0.5), tol.paramTol);
        T step_u = std::max<T>(du * T(0.5), param_eps);
        T step_v = std::max<T>(dv * T(0.5), param_eps);

        // Curvature gradient magnitude at center
        T G = curvatureGradientMagnitude(uc, vc, step_u, step_v, tol);

        // --- Tolerance-normalized flatness tests ---
        // Parametric span (dimensionless)
        T s_param = std::sqrt(du*du + dv*dv);

        // Curvature change across patch (dimensionless): G ~ |∇kappa| per param
        // Threshold blends a small absolute tolerance with param scaled epsilon
        T curv_thresh = tol.evaluateEpsilon(T(1)) + (eps_local / J);
        bool curvatureFlat = (G * s_param) <= curv_thresh;

        bool geomFlat = bbox_diag <= eps_local * (1 + J * (du + dv));

        bool depthStop = (depth >= max_depth);
        bool paramStop = (du <= min_param_extent && dv <= min_param_extent);

        if (depthStop || paramStop || (geomFlat && curvatureFlat)) {
            bvh_patches.push_back(SurfacePatchAABB{umin, umax, vmin, vmax, minPt, maxPt});
            return;
        }

        bool split_u = true;
        bool split_v = true;

        if (have_dir) {
            // Project direction hint onto param space using J = [Su Sv] at the patch center
            // (Su, Sv already computed above)
            Eigen::Matrix<T,2,2> JTJ;
            JTJ(0,0) = Su.dot(Su);
            JTJ(0,1) = Su.dot(Sv);
            JTJ(1,0) = JTJ(0,1);
            JTJ(1,1) = Sv.dot(Sv);
            Eigen::Matrix<T,2,1> JTd;
            JTd(0) = Su.dot(dir_hint);
            JTd(1) = Sv.dot(dir_hint);

            T scale_J = std::max<T>(JTJ(0,0) + JTJ(1,1), tol.paramTol);
            T det_eps = std::max<T>(tol.evaluateEpsilon(scale_J), limits.eps_det);
            T det = JTJ(0,0)*JTJ(1,1) - JTJ(0,1)*JTJ(1,0);
            if (std::abs(det) > det_eps) {
                Eigen::Matrix<T,2,1> ab;
                ab(0) = ( JTJ(1,1)*JTd(0) - JTJ(0,1)*JTd(1)) / det;
                ab(1) = (-JTJ(1,0)*JTd(0) + JTJ(0,0)*JTd(1)) / det;
                T au = std::abs(ab(0));
                T bv = std::abs(ab(1));
                // Adaptive anisotropy threshold based on tolerance scale and local_scale
                T anisotropyThresh = T(2) + tol.evaluateEpsilon(local_scale) * T(1e3);
                T ratio = (std::max(au, bv)) / (std::max(limits.eps_cross, std::min(au, bv)));
                if (ratio > anisotropyThresh) {
                    if (au > bv) { split_u = true; split_v = false; }
                    else         { split_u = false; split_v = true; }
                }
            }
        }

        T uc_mid = (umin + umax) / 2;
        T vc_mid = (vmin + vmax) / 2;

        if (split_u && split_v) {
            subdividePatch(umin, uc_mid, vmin, vc_mid, tol, min_param_extent, max_depth, depth + 1, budget - 1, have_dir, dir_hint);
            subdividePatch(uc_mid, umax, vmin, vc_mid, tol, min_param_extent, max_depth, depth + 1, budget - 1, have_dir, dir_hint);
            subdividePatch(umin, uc_mid, vc_mid, vmax, tol, min_param_extent, max_depth, depth + 1, budget - 1, have_dir, dir_hint);
            subdividePatch(uc_mid, umax, vc_mid, vmax, tol, min_param_extent, max_depth, depth + 1, budget - 1, have_dir, dir_hint);
        } else if (split_u) {
            subdividePatch(umin, uc_mid, vmin, vmax, tol, min_param_extent, max_depth, depth + 1, budget - 1, have_dir, dir_hint);
            subdividePatch(uc_mid, umax, vmin, vmax, tol, min_param_extent, max_depth, depth + 1, budget - 1, have_dir, dir_hint);
        } else {
            subdividePatch(umin, umax, vmin, vc_mid, tol, min_param_extent, max_depth, depth + 1, budget - 1, have_dir, dir_hint);
            subdividePatch(umin, umax, vc_mid, vmax, tol, min_param_extent, max_depth, depth + 1, budget - 1, have_dir, dir_hint);
        }
    }

    // Adaptive recursive BVH generation (tolerance-driven)
    void generateBVH(int grid, const Euclid::Tolerance& tol, int max_depth = -1) const {
        // Will set J after partials below
        bvh_patches.clear();
        curvature_cache.clear();
        // Pre-reserve BVH patch vector capacity to avoid excessive reallocations during recursion
        auto limits_prealloc = NumericalLimits<T>::derive(tol, tol.paramTol);
        size_t reserve_cap = static_cast<size_t>(
            std::clamp<T>(4096 / std::max<T>(tol.paramTol, limits_prealloc.eps_param), T(1024), T(32768))
        );
        bvh_patches.reserve(reserve_cap);
        // Domain
        T umin = surfaceDomain.first.first;
        T umax = surfaceDomain.first.second;
        T vmin = surfaceDomain.second.first;
        T vmax = surfaceDomain.second.second;

        // Estimate scale if not provided: bbox diagonal of corners + center
        auto p00 = evaluate(umin, vmin).coords;
        auto p10 = evaluate(umax, vmin).coords;
        auto p01 = evaluate(umin, vmax).coords;
        auto p11 = evaluate(umax, vmax).coords;
        auto pc  = evaluate((umin+umax)/T(2), (vmin+vmax)/T(2)).coords;
        Eigen::Matrix<T,N,1> bbmin = p00, bbmax = p00;
        auto upd = [&](const Eigen::Matrix<T,N,1>& q){
            for (int k=0;k<N;++k){ bbmin(k)=std::min(bbmin(k),q(k)); bbmax(k)=std::max(bbmax(k),q(k)); }
        };
        upd(p10); upd(p01); upd(p11); upd(pc);
        T est_scale = (bbmax - bbmin).norm();
        if (est_scale == T(0)) est_scale = T(1);

        // Parametric step lower bound derived adaptively from surface extent and Jacobian magnitude at center
        T uc = (umin + umax) / T(2);
        T vc = (vmin + vmax) / T(2);
        auto [Su0, Sv0] = evaluatePartials(uc, vc);
        T J = std::max(Su0.norm(), Sv0.norm());
        auto limits = NumericalLimits<T>::derive(tol, J);

        // Use local surface extent as geometric scale
        T surface_extent = std::max<T>(std::abs(umax - umin), std::abs(vmax - vmin));
        T local_scale = std::max<T>(surface_extent * J, tol.paramTol);

        // Compute adaptive epsilon using the tolerance model
        T adaptive_eps = tol.evaluateEpsilon(local_scale);

        // Map adaptive epsilon (based on local surface extent) to param-space resolution
        T param_eps_local = adaptive_eps / std::max<T>(J, tol.paramTol);

        // Grid-based cap in param space
        T grid_eps_u = std::abs(umax-umin) / T(std::max(1, grid));
        T grid_eps_v = std::abs(vmax-vmin) / T(std::max(1, grid));
        T min_param_extent = std::max(param_eps_local, T(0.25) * std::min(grid_eps_u, grid_eps_v));

        int effective_depth = max_depth;
        if (effective_depth <= 0) {
            // Estimate curvature-driven refinement need at domain center
            T uc_mid = (umin + umax) / T(2);
            T vc_mid = (vmin + vmax) / T(2);
            T G_center = curvatureGradientMagnitude(uc_mid, vc_mid,
                                                    std::max<T>(std::abs(umax - umin) * T(0.25), tol.paramTol),
                                                    std::max<T>(std::abs(vmax - vmin) * T(0.25), tol.paramTol),
                                                    tol);

            // Compute domain ratio and tolerance scaling
            T span = std::max(std::abs(umax - umin), std::abs(vmax - vmin));
            T ratio = std::max<T>(T(1), span / std::max<T>(min_param_extent, limits.eps_param));
            T tol_scale = std::clamp<T>(T(1.0) / std::sqrt(std::max<T>(tol.paramTol, limits.eps_param)), T(1), T(64));

            // Adaptive depth scaling: curvature and tolerance drive deeper recursion
            T curvature_scale = std::clamp<T>(1 + G_center * 20, 1, 8);
            int depth_est = int(std::ceil(std::log2(double(ratio * curvature_scale * tol_scale))));
            effective_depth = std::clamp<int>(depth_est, 3, 18);
        }

        // Pass optional direction hint (if any)
        bool have_dir = bvh_dir_hint.has_value();
        Eigen::Matrix<T,N,1> dir_hint = have_dir ? *bvh_dir_hint : Eigen::Matrix<T,N,1>::Zero();

        // Per-branch work budget derived from adaptive tolerance scaling and resolution ratio
        T tol_clamped = std::max<T>(tol.paramTol, limits.eps_param);
        T domain_extent = std::max<T>(std::abs(umax - umin), std::abs(vmax - vmin));
        T ratio = domain_extent / std::max<T>(min_param_extent, tol_clamped);
        T adaptive_scale = std::max<T>(1, std::log2(std::max<T>(T(2), ratio)));
        int base_budget = 6 + int(std::ceil(std::log10(double(1.0 / tol_clamped))));
        int per_branch_budget = std::clamp<int>(int(base_budget * adaptive_scale), 8, 24);

        subdividePatch(umin, umax, vmin, vmax, tol, min_param_extent, effective_depth, 0, per_branch_budget, have_dir, dir_hint);
        bvh_valid = true;
    }

    // Ray vs AABB slab test in N-D; returns (hit, tmin, tmax)
    std::tuple<bool, T, T> rayAABBInterval(const Eigen::Matrix<T, N, 1>& rayOrigin,
                                           const Eigen::Matrix<T, N, 1>& rayDir,
                                           const Eigen::Matrix<T, N, 1>& bmin,
                                           const Eigen::Matrix<T, N, 1>& bmax,
                                           const Euclid::Tolerance& tol,
                                           T local_scale = T(1)) const
    {
        T tmin = -std::numeric_limits<T>::infinity();
        T tmax =  std::numeric_limits<T>::infinity();

        for (int k = 0; k < N; ++k) {
            T o = rayOrigin(k);
            T d = rayDir(k);
            // Use a tolerance-scaled threshold based on the slab width in this dimension
            T slab_extent = std::max<T>(std::abs(bmax(k) - bmin(k)), tol.evaluateEpsilon(local_scale));
            T eps_d = tol.evaluateEpsilon(slab_extent);
            T invD;
            if (std::abs(d) <= eps_d) {
                // Ray effectively parallel to slab: origin must lie within (padded externally)
                if (o < bmin(k) || o > bmax(k)) return {false, T(0), T(0)};
                continue;
            } else {
                invD = T(1) / d;
            }
            T t0 = (bmin(k) - o) * invD;
            T t1 = (bmax(k) - o) * invD;
            if (t0 > t1) std::swap(t0, t1);
            tmin = std::max(tmin, t0);
            tmax = std::min(tmax, t1);
            if (tmax < tmin) return {false, T(0), T(0)};
        }
        return {true, tmin, tmax};
    }

    // Lightweight predicate using the slab test (with epsilon padding)
    bool rayIntersectsAABBND(const Eigen::Matrix<T, N, 1>& rayOrigin,
                             const Eigen::Matrix<T, N, 1>& rayDir,
                             const Eigen::Matrix<T, N, 1>& bmin,
                             const Eigen::Matrix<T, N, 1>& bmax,
                             const Euclid::Tolerance& tol,
                             T local_scale = T(1)) const
    {
        Eigen::Matrix<T, N, 1> pad = Eigen::Matrix<T, N, 1>::Constant(tol.evaluateEpsilon(local_scale));
        auto [hit, t0, t1] = rayAABBInterval(rayOrigin, rayDir, bmin - pad, bmax + pad, tol, local_scale);
        (void)t0; (void)t1;
        return hit;
    }

    // Simplified, normal-oriented signed distance: uses surface normal for sign and removes extra near-surface sampling.
    T rayPointSignedDistanceToSurface(const Eigen::Matrix<T, N, 1>& P,
                                      const Euclid::Tolerance& tol) const
    {
        // Prefer analytic SDF when available
        if (analyticSDF) {
            return analyticSDF(P);
        }

        // Project to the surface to get nearest point and parameters
        auto proj = projectPoint(P, tol);
        const Eigen::Matrix<T,N,1> Q = proj.surface_point.coords;
        const Eigen::Matrix<T,N,1> diff = P - Q;

        // Compute an oriented normal from partials; fall back to principalNormal in ND
        Eigen::Matrix<T,N,1> n;
        auto [Su, Sv] = evaluatePartials(proj.u, proj.v);
        if constexpr (N == 3) {
            n = Su.cross(Sv);
            T nn = n.norm();
            if (nn > T(0)) n /= nn; else n = principalNormal(proj.u, proj.v, tol, bvh_dir_hint);
        } else {
            n = principalNormal(proj.u, proj.v, tol, bvh_dir_hint);
        }

        // If we still have a degenerate normal, use the direction to the surface
        if (n.norm() <= T(0)) {
            if (diff.norm() > T(0)) n = diff.normalized();
            else n = Eigen::Matrix<T,N,1>::UnitX();
        }


        // Signed distance is dot with oriented normal
        return diff.dot(n);
    }

    // Compute a conservative bbox for a param patch by sampling its corners and center.
    SurfacePatchAABB makePatchAABB(T umin, T umax, T vmin, T vmax) const
    {
        T uc = (umin + umax) / T(2);
        T vc = (vmin + vmax) / T(2);
        std::array<PointType,5> samples = {
            evaluate(umin, vmin),
            evaluate(umax, vmin),
            evaluate(umin, vmax),
            evaluate(umax, vmax),
            evaluate(uc,   vc  )
        };
        PointType minPt = samples[0];
        PointType maxPt = samples[0];
        for (const auto& p : samples) {
            for (int k=0;k<N;++k) {
                if (p.coords(k) < minPt.coords(k)) minPt.coords(k) = p.coords(k);
                if (p.coords(k) > maxPt.coords(k)) maxPt.coords(k) = p.coords(k);
            }
        }
        return SurfacePatchAABB{umin, umax, vmin, vmax, minPt, maxPt};
    }

    // Recursive refinement of a single patch guided by the ray and tolerance.
    void refinePatchNearRayRecursive(const Eigen::Matrix<T, N, 1>& rayOrigin,
                                     const Eigen::Matrix<T, N, 1>& rayDir,
                                     SurfacePatchAABB& patch,
                                     const Euclid::Tolerance& tol,
                                     int depth,
                                     int max_depth,
                                     int budget) const
    {
        // Derive local geometry and tolerance-driven scales
        T du = std::abs(patch.umax - patch.umin);
        T dv = std::abs(patch.vmax - patch.vmin);
        T surface_extent = std::max<T>(du, dv);
        T uc = (patch.umin + patch.umax) / T(2);
        T vc = (patch.vmin + patch.vmax) / T(2);
        auto [Su, Sv] = evaluatePartials(uc, vc);
        T J = std::max(Su.norm(), Sv.norm());
        T local_scale = std::max<T>(surface_extent * J, tol.paramTol);
        T eps_local = tol.evaluateEpsilon(local_scale);
        T param_eps_local = eps_local / std::max<T>(J, tol.paramTol);
        auto limits_J = NumericalLimits<T>::derive(tol, std::max<T>(J, tol.paramTol));

        // Hard cap to avoid explosions in pathological cases, dynamically scaled by tolerance
        size_t MAX_EXTRA = std::clamp<size_t>(
            static_cast<size_t>(8192 / std::max<T>(tol.paramTol, limits_J.eps_param)),
            size_t(2048),
            size_t(32768)
        );
        if (bvh_patches.size() >= MAX_EXTRA) return;
        if (budget <= 0) return;

        // Quick reject: ray vs AABB
        Eigen::Matrix<T,N,1> bmin = patch.minPoint.coords;
        Eigen::Matrix<T,N,1> bmax = patch.maxPoint.coords;
        if (!rayIntersectsAABBND(rayOrigin, rayDir, bmin, bmax, tol, local_scale)) return;

        // If the bbox is already within epsilon, stop.
        T bbox_diag = (bmax - bmin).norm();
        bool small_geom = (bbox_diag <= eps_local);
        bool small_param = (std::max(du, dv) <= param_eps_local);

        // Stop if depth exhausted, or if either space has converged sufficiently and curvature is mild
        if (depth >= max_depth) return;
        // Lightweight curvature proxy to avoid needless splitting when already tiny
        T G_pre = curvatureGradientMagnitude(uc, vc, std::max(du * T(0.5), param_eps_local),
                                             std::max(dv * T(0.5), param_eps_local), tol);
        bool mild_curv = (G_pre * (du + dv)) <= (tol.evaluateEpsilon(std::max<T>(tol.paramTol, limits_J.eps_param)) + eps_local / std::max<T>(J, tol.paramTol));
        if ((small_geom || small_param) && mild_curv) return;

        // Find t-interval of ray within the AABB
        auto [hit, t0, t1] = rayAABBInterval(rayOrigin, rayDir, bmin, bmax, tol, local_scale);
        if (!hit) return;

        // Local plane proxy at patch center to avoid iterative projection during sampling
        const Eigen::Matrix<T,N,1> S0 = evaluate(uc, vc).coords;
        auto [Su_loc, Sv_loc] = evaluatePartials(uc, vc);
        Eigen::Matrix<T,N,1> n_loc;
        if constexpr (N == 3) {
            n_loc = Su_loc.cross(Sv_loc);
            T nn = n_loc.norm();
            if (nn > T(0)) n_loc /= nn; else n_loc = principalNormal(uc, vc, tol, bvh_dir_hint);
        } else {
            n_loc = principalNormal(uc, vc, tol, bvh_dir_hint);
        }
        if (n_loc.norm() <= T(0)) {
            // fallback: use ray direction as a stable proxy
            n_loc = (rayDir.norm() > T(0)) ? rayDir.normalized() : Eigen::Matrix<T,N,1>::UnitX();
        }
        auto proxySigned = [&](const Eigen::Matrix<T,N,1>& P) -> T {
            return (P - S0).dot(n_loc);
        };

        // Adaptive sampling along the ray segment inside the AABB, using the local proxy
        auto sample_segment_adaptive = [&](int samples, bool& sign_change_out, bool& near_zero_out){
            T dist_prev = std::numeric_limits<T>::quiet_NaN();
            sign_change_out = false; near_zero_out = false;
            for (int i = 0; i < samples; ++i) {
                T s = (samples == 1) ? T(0.5) : T(i) / T(samples - 1);
                T t = t0 + (t1 - t0) * s;
                Eigen::Matrix<T,N,1> P = rayOrigin + t * rayDir;
                T d = proxySigned(P);
                if (std::abs(d) <= eps_local) near_zero_out = true;
                if (!std::isnan(dist_prev) && d * dist_prev < T(0)) { sign_change_out = true; break; }
                dist_prev = d;
            }
        };

        // Choose number of samples based on how large the bbox segment is relative to epsilon
        // Larger boxes (relative to eps_local) get more samples, up to a modest cap.
        T rel_floor = std::max<T>(tol.paramTol, limits_J.eps_param);
        T rel = std::max<T>(rel_floor, bbox_diag / std::max<T>(eps_local, rel_floor));
        int base_samples = 3 + int(std::ceil(std::log2(double(rel))));
        int samples_coarse = std::clamp(base_samples, 3, 7);   // 3,5,7 typical
        // Adaptive refinement logic (purely tolerance-adaptive, no hardcoded constants)
        T scale_factor = T(1) / std::sqrt(std::max<T>(tol.paramTol, limits_J.eps_param));
        T refine_thresh_geom = eps_local * scale_factor;
        T refine_thresh_curv = std::max<T>(limits_J.eps_param, tol.paramTol);

        bool refine_needed =
            (bbox_diag > refine_thresh_geom) &&
            (G_pre > refine_thresh_curv);

        int samples_refine = refine_needed
            ? std::clamp(base_samples + 2, 5, 9)
            : samples_coarse;

        bool sign_change = false, near_zero = false;
        sample_segment_adaptive(samples_coarse, sign_change, near_zero);

        // If uncertain and bbox is not tiny, do a denser pass
        if (!sign_change && !near_zero && bbox_diag > eps_local * T(2) && samples_refine > samples_coarse) {
            sample_segment_adaptive(samples_refine, sign_change, near_zero);
        }

        // If proxy was inconclusive on a sizable box, do a single exact mid-point check to avoid false negatives
        if (!sign_change && !near_zero && bbox_diag > eps_local * T(4)) {
            T t_mid = (t0 + t1) * T(0.5);
            Eigen::Matrix<T,N,1> Pmid = rayOrigin + t_mid * rayDir;
            T d_exact = rayPointSignedDistanceToSurface(Pmid, tol);
            if (std::abs(d_exact) <= eps_local) near_zero = true;
        }

        // Curvature gradient at center (parametric)
        T G = curvatureGradientMagnitude(uc, vc, std::max(du * T(0.5), param_eps_local),
                                         std::max(dv * T(0.5), param_eps_local), tol);

        // Adaptive curvature threshold and geometric flatness
        T curv_thresh = std::max<T>(limits_J.eps_param, limits_J.eps_param / std::max<T>(J, tol.paramTol));
        bool high_curvature = (G * (du + dv)) > curv_thresh;

        bool geomFlat = bbox_diag <= eps_local * (1 + J * (du + dv));

        // Early exit: if there is no sign change or near-zero crossing and the patch is both geom-flat and not high curvature, stop.
        if (!sign_change && !near_zero && geomFlat && !high_curvature) return;

        // Otherwise, subdivide the patch and recurse.
        T um = (patch.umin + patch.umax) / T(2);
        T vm = (patch.vmin + patch.vmax) / T(2);

        SurfacePatchAABB q00 = makePatchAABB(patch.umin, um, patch.vmin, vm);
        SurfacePatchAABB q10 = makePatchAABB(um, patch.umax, patch.vmin, vm);
        SurfacePatchAABB q01 = makePatchAABB(patch.umin, um, vm, patch.vmax);
        SurfacePatchAABB q11 = makePatchAABB(um, patch.umax, vm, patch.vmax);

        refinePatchNearRayRecursive(rayOrigin, rayDir, q00, tol, depth+1, max_depth, budget - 1);
        refinePatchNearRayRecursive(rayOrigin, rayDir, q10, tol, depth+1, max_depth, budget - 1);
        refinePatchNearRayRecursive(rayOrigin, rayDir, q01, tol, depth+1, max_depth, budget - 1);
        refinePatchNearRayRecursive(rayOrigin, rayDir, q11, tol, depth+1, max_depth, budget - 1);

        // Replace current patch with one child and append the others
        patch = q00;
        bvh_patches.push_back(q10);
        bvh_patches.push_back(q01);
        bvh_patches.push_back(q11);
    }

public:
    // Public entry: refine existing BVH only where the ray traverses, up to a limited depth
    const std::vector<SurfacePatchAABB>& refineBVHNearRay(
        const Eigen::Matrix<T, N, 1>& rayOrigin,
        const Eigen::Matrix<T, N, 1>& rayDir,
        const Euclid::Tolerance& tol,
        int max_depth = 9) const
    {
        curvature_cache.clear();
        // Ensure BVH exists; reuse scale estimated during build
        if (!bvh_valid) {
            generateBVH(40, tol);
        }

        // Estimate a local geometric extent from current BVH bounds (quick pass)
        Eigen::Matrix<T,N,1> bbmin = bvh_patches.front().minPoint.coords;
        Eigen::Matrix<T,N,1> bbmax = bvh_patches.front().maxPoint.coords;
        for (const auto& p : bvh_patches) {
            for (int k=0;k<N;++k) {
                bbmin(k) = std::min(bbmin(k), p.minPoint.coords(k));
                bbmax(k) = std::max(bbmax(k), p.maxPoint.coords(k));
            }
        }
        T geom_extent = (bbmax - bbmin).norm();
        geom_extent = std::max<T>(geom_extent, tol.paramTol);

        // Estimate a representative Jacobian magnitude at the domain center
        T u_mid = (surfaceDomain.first.first + surfaceDomain.first.second)/T(2);
        T v_mid = (surfaceDomain.second.first + surfaceDomain.second.second)/T(2);
        auto [Su_c, Sv_c] = evaluatePartials(u_mid, v_mid, tol);
        T J_bvh = std::max(Su_c.norm(), Sv_c.norm());
        J_bvh = std::max<T>(J_bvh, tol.paramTol);

        // Pre-reserve to avoid reallocation that would invalidate references during refinement
        {
            size_t start_sz = bvh_patches.size();
            // A soft upper bound for new patches created during refinement
            // (matches recursive guard caps and gives slack).
            size_t soft_budget = start_sz + size_t(65536);
            if (bvh_patches.capacity() < soft_budget) {
                bvh_patches.reserve(soft_budget);
            }
        }

        // Important: iterate over a snapshot of current size because vector grows during refinement.
        // Work on a local copy per element to avoid reference invalidation if pushes occur.
        auto limits_rb = NumericalLimits<T>::derive(tol, J_bvh);
        T tol_clamped = std::max<T>(tol.paramTol, limits_rb.eps_param);
        int refine_budget = std::clamp(6 + int(std::ceil(std::log10(double(1.0 / tol_clamped)))) + max_depth/2, 8, 24);
        size_t initial = bvh_patches.size();
        for (size_t i = 0; i < initial; ++i) {
            SurfacePatchAABB local = bvh_patches[i];                   // copy
            refinePatchNearRayRecursive(rayOrigin, rayDir, local, tol, /*depth*/0, max_depth, refine_budget);
            bvh_patches[i] = local;                                    // write back
        }
        // BVH remains valid; patches were refined in-place.
        return bvh_patches;
    }

    const std::vector<SurfacePatchAABB>& getBVH(
        int grid,
        const Euclid::Tolerance& tol,
        int max_depth = 9) const
    {
        if (!bvh_valid) {
            generateBVH(grid, tol, max_depth);
        }
        return bvh_patches;
    }

    // ======================================================
    // Candidate Sampling using Cached BVH
    // ======================================================
    std::vector<std::pair<T, T>> generateCandidateSamples(
        const Euclid::Tolerance& tol,
        int max_grid = 40,
        std::function<bool(const Eigen::Matrix<T, N, 1>& minpt,
                           const Eigen::Matrix<T, N, 1>& maxpt)> filter = nullptr) const
    {
        auto domain = surfaceDomain;
        T u0 = domain.first.first;
        T u1 = domain.first.second;
        T v0 = domain.second.first;
        T v1 = domain.second.second;
        if (u1 <= u0 || v1 <= v0) return {};

        T extent = std::max(std::abs(u1 - u0), std::abs(v1 - v0));
        const auto& patches = getBVH(max_grid > 0 ? max_grid : 40, tol);

        std::vector<std::pair<T, T>> samples;
        samples.reserve(patches.size() * 36 + 16);

        // Helper to add a unique sample (u,v) to samples, using tolerance on extent
        auto add_unique = [&](T uu, T vv) {
            for (const auto& s : samples)
                if (std::abs(s.first - uu) < tol.evaluateEpsilon(extent) &&
                    std::abs(s.second - vv) < tol.evaluateEpsilon(extent))
                    return;
            samples.emplace_back(uu, vv);
        };

        for (const auto& patch : patches) {
            if (filter && !filter(patch.minPoint.coords, patch.maxPoint.coords))
                continue;

            // Local param extents
            T du = std::abs(patch.umax - patch.umin);
            T dv = std::abs(patch.vmax - patch.vmin);

            // Always insert center and corners (cheap and robust)
            T uc = (patch.umin + patch.umax) / 2;
            T vc = (patch.vmin + patch.vmax) / 2;
            add_unique(uc, vc);
            add_unique(patch.umin, patch.vmin);
            add_unique(patch.umin, patch.vmax);
            add_unique(patch.umax, patch.vmin);
            add_unique(patch.umax, patch.vmax);

            // Check if patch is large in param-space relative to tolerance
            bool large_extent = (du > 2 * tol.evaluateEpsilon(std::max(std::abs(u1 - u0), std::abs(v1 - v0)))) ||
                                (dv > 2 * tol.evaluateEpsilon(std::max(std::abs(u1 - u0), std::abs(v1 - v0))));

            // Curvature/conditioning indicator at patch center (tolerance-adaptive)
            auto [Su, Sv] = evaluatePartials(uc, vc);
            T su = Su.norm();
            T sv = Sv.norm();
            T mixed = std::abs(Su.dot(Sv));
            // Use adaptive indicator: normalize by tolerance and numerical limits
            auto limits = NumericalLimits<T>::derive(tol, std::max<T>(du, dv));
            T indicator = (su * sv + mixed) / std::max<T>(limits.eps_param, tol.paramTol);

            // Adaptive refinement thresholds based on tolerance and limits
            T thresh_lo  = std::max<T>(tol.paramTol * 5,   limits.eps_param * 10);
            T thresh_mid = std::max<T>(tol.paramTol * 25,  limits.eps_param * 50);
            T thresh_hi  = std::max<T>(tol.paramTol * 100, limits.eps_param * 200);

            // Choose refinement level based on adaptive indicator
            // 1 => only center/corners, 3 => 3x3, 5 => 5x5, 7 => 7x7
            int refine = 1;
            if (large_extent || indicator > thresh_lo)  refine = std::max(refine, 3);
            if (indicator > thresh_mid)                 refine = std::max(refine, 5);
            if (indicator > thresh_hi)                  refine = std::max(refine, 7);

            if (refine > 1) {
                // Place a regular (refine+1)x(refine+1) grid including boundaries to
                // improve seeding density in oscillatory regions.
                for (int i = 0; i <= refine; ++i) {
                    T fu = T(i) / T(refine);
                    T uu = patch.umin + du * fu;
                    for (int j = 0; j <= refine; ++j) {
                        T fv = T(j) / T(refine);
                        T vv = patch.vmin + dv * fv;
                        add_unique(uu, vv);
                    }
                }
            }
        }

        // Always add the domain center and corners
        add_unique((u0 + u1) / 2, (v0 + v1) / 2);
        add_unique(u0, v0);
        add_unique(u0, v1);
        add_unique(u1, v0);
        add_unique(u1, v1);

        return samples;
    }
};

// Mesh representation
template<typename T, int N>
class SurfaceMesh {
public:
    using PointType = Point<T,N>;

    std::vector<PointType> vertices;
    std::vector<std::vector<size_t>> faces;

    SurfaceMesh() = default;

    bool validate() const {
        size_t numVertices = vertices.size();
        for (const auto& f : faces) {
            if (f.size() < 3) return false;
            for (size_t idx : f) {
                if (idx >= numVertices) return false;
            }
        }
        return true;
    }

    // Apply a transformation to the mesh
    template<typename Transform>
    SurfaceMesh<T,N> applyTransform(const Transform& transform) const {
        SurfaceMesh<T,N> meshT;
        meshT.vertices.reserve(vertices.size());
        for (auto& v : vertices)
            meshT.vertices.push_back(transform.apply(v));
        meshT.faces = faces;
        return meshT;
    }

    // OBJ export only supports 3D
    std::optional<std::string> exportOBJ() const {
        if constexpr (N != 3) {
            return std::nullopt;
        } else {
            std::string s;
            for (auto& v : vertices) {
                s += "v ";
                for (int i=0;i<3;++i) s += std::to_string(v.coords(i)) + " ";
                s += "\n";
            }
            for (auto& f : faces) {
                s += "f ";
                for (auto idx : f) s += std::to_string(idx+1) + " ";
                s += "\n";
            }
            return s;
        }
    }

    // Approximate area for 3D surfaces
    T area() const {
        if constexpr (N != 3) return T(0);
        T total = 0;
        for (auto& f : faces) {
            if (f.size() != 3) continue;
            auto& a = vertices[f[0]].coords;
            auto& b = vertices[f[1]].coords;
            auto& c = vertices[f[2]].coords;
            auto ab = b - a;
            auto ac = c - a;
            total += (ab.cross(ac)).norm() / 2;
        }
        return total;
    }

    struct TopologyCheckResult {
        bool isClosed;
        std::vector<std::vector<size_t>> boundarySubfaces;
    };

    TopologyCheckResult hasClosedTopology(const Euclid::Tolerance& tol) const {
        TopologyCheckResult result{true, {}};
        if constexpr (N < 2) {
            // Not defined for dimension less than 2 (no faces)
            result.isClosed = false;
            return result;
        }

        // Determine subface size: (N-2)-faces have (N-1) vertices
        const int subfaceSize = N - 1;

        // Map from canonical subface (sorted indices) to count
        std::map<std::vector<size_t>, int> subfaceCount;

        // Helper lambda to check if face is degenerate (repeated vertices or near-zero measure)
        auto isDegenerateFace = [&](const std::vector<size_t>& face) -> bool {
            std::set<size_t> uniqueVerts(face.begin(), face.end());
            if (uniqueVerts.size() < face.size()) return true;

            if constexpr (N == 3) {
                // For 3D triangles, also check near-zero area using tolerance-driven epsilon
                if (face.size() == 3) {
                    const auto& a = vertices[face[0]].coords;
                    const auto& b = vertices[face[1]].coords;
                    const auto& c = vertices[face[2]].coords;
                    auto ab = b - a;
                    auto ac = c - a;
                    auto cross = ab.cross(ac);
                    // Scale area tolerance by local edge lengths to remain scale-aware
                    T len_ab = ab.norm();
                    T len_ac = ac.norm();
                    T scale = std::max<T>(tol.paramTol, len_ab * len_ac);
                    if (cross.norm() < tol.evaluateEpsilon(scale)) return true;
                }
            }
            return false;
        };

        for (const auto& face : faces) {
            if (face.size() < subfaceSize) {
                // Skip degenerate or too small faces
                continue;
            }
            if (isDegenerateFace(face)) continue;

            if constexpr (N == 3) {
                // For 3D triangular mesh, (N-2)-faces are edges (pairs of vertices)
                // Each face has 3 edges
                for (size_t i = 0; i < 3; ++i) {
                    std::vector<size_t> edge = {face[i], face[(i+1)%3]};
                    std::sort(edge.begin(), edge.end());
                    subfaceCount[edge]++;
                }
            } else {
                // For general ND, generate all (N-2)-faces by removing one vertex from the face
                // Each (N-2)-face has subfaceSize = N-1 vertices
                for (size_t i = 0; i < face.size(); ++i) {
                    std::vector<size_t> subface;
                    subface.reserve(subfaceSize);
                    for (size_t j = 0; j < face.size(); ++j) {
                        if (j != i) {
                            subface.push_back(face[j]);
                        }
                    }
                    std::sort(subface.begin(), subface.end());
                    subfaceCount[subface]++;
                }
            }
        }

        // Identify boundary subfaces (those with count != 2)
        for (const auto& [subface, count] : subfaceCount) {
            if (count != 2) {
                result.isClosed = false;
                result.boundarySubfaces.push_back(subface);
            }
        }

        return result;
    }

    // Convenience overload that uses a default tolerance
    TopologyCheckResult hasClosedTopology() const {
        Euclid::Tolerance tol;
        return hasClosedTopology(tol);
    }
};

// Generate a grid-based mesh
template<typename T, int N>
SurfaceMesh<T,N> generateSurfaceMesh(const Surface<T,N>& surf, int uSteps, int vSteps) {
    SurfaceMesh<T,N> mesh;

    // Correctly interpret domain
    T umin = surf.surfaceDomain.first.first;
    T umax = surf.surfaceDomain.first.second;
    T vmin = surf.surfaceDomain.second.first;
    T vmax = surf.surfaceDomain.second.second;

    // Steps
    T stepU = (umax - umin)/T(uSteps-1);
    T stepV = (vmax - vmin)/T(vSteps-1);

    // Sample points
    std::vector<T> us(uSteps);
    std::vector<T> vs(vSteps);
    for(int i=0;i<uSteps;++i) us[i] = umin + stepU*i;
    for(int j=0;j<vSteps;++j) vs[j] = vmin + stepV*j;

    // Vertices (row-major: v first, then u)
    for(int j=0;j<vSteps;++j){
        for(int i=0;i<uSteps;++i){
            mesh.vertices.push_back(surf.evaluate(us[i], vs[j]));
        }
    }

    // Faces (triangles from quads)
    for(int j=0;j<vSteps-1;++j){
        for(int i=0;i<uSteps-1;++i){
            size_t idx = i + j*uSteps;
            size_t idxRight = idx+1;
            size_t idxDown = idx+uSteps;
            size_t idxDiag = idx+uSteps+1;
            mesh.faces.push_back({idx, idxDiag, idxRight});
            mesh.faces.push_back({idx, idxDown, idxDiag});
        }
    }

    return mesh;
}

template<typename T>
SurfaceMesh<T,3> generatePeriodicWrappedMesh(
    const Surface<T,3>& surf,
    int uSteps,
    int vSteps,
    bool uPeriodic = false,
    bool vPeriodic = false
) {
    SurfaceMesh<T,3> mesh;
    auto& domain = surf.surfaceDomain;
    T umin = domain.first.first;
    T vmin = domain.second.first;
    T umax = domain.first.second;
    T vmax = domain.second.second;

    T stepU = (umax - umin) / T(uSteps - 1);
    T stepV = (vmax - vmin) / T(vSteps - 1);

    // Determine effective number of vertices in u and v directions
    int vertCountU = uPeriodic ? uSteps - 1 : uSteps;
    int vertCountV = vPeriodic ? vSteps - 1 : vSteps;

    // Generate vertices with seam vertices shared if periodic
    for (int j = 0; j < vertCountV; ++j) {
        T v = vmin + stepV * j;
        for (int i = 0; i < vertCountU; ++i) {
            T u = umin + stepU * i;
            mesh.vertices.push_back(surf.evaluate(u, v));
        }
    }

    // Helper to get vertex index with wrapping
    auto vertexIndex = [&](int i, int j) -> size_t {
        if (uPeriodic) {
            i = i % vertCountU;
            if (i < 0) i += vertCountU;
        }
        if (vPeriodic) {
            j = j % vertCountV;
            if (j < 0) j += vertCountV;
        }
        return i + j * vertCountU;
    };

    // Generate faces
    int faceCountU = uPeriodic ? vertCountU : vertCountU - 1;
    int faceCountV = vPeriodic ? vertCountV : vertCountV - 1;

    for (int j = 0; j < faceCountV; ++j) {
        for (int i = 0; i < faceCountU; ++i) {
            size_t idx00 = vertexIndex(i, j);
            size_t idx10 = vertexIndex(i + 1, j);
            size_t idx01 = vertexIndex(i, j + 1);
            size_t idx11 = vertexIndex(i + 1, j + 1);

            mesh.faces.push_back({idx00, idx10, idx11});
            mesh.faces.push_back({idx00, idx11, idx01});
        }
    }

    return mesh;
}

} // namespace Euclid::Geometry
