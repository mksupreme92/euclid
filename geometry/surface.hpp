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
                                                           const Euclid::Tolerance& tol,
                                                           int max_seeds = 6) const
    {
        std::vector<std::pair<T,T>> seeds;
        if (!bvh_valid) generateBVH(40, tol);
        const auto& patches = bvh_patches;
        // Compute distances from P to each patch AABB and keep a few best
        struct Item { T d2; T umin, umax, vmin, vmax; };
        std::vector<Item> items;
        items.reserve(patches.size());
        //auto sqr = [](T x){ return x*x; };
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
        return seeds;
    }

private:
    mutable std::vector<SurfacePatchAABB> bvh_patches;
    mutable bool bvh_valid = false;
    // Optional world-space direction hint for directional BVH refinement
    mutable std::optional<Eigen::Matrix<T, N, 1>> bvh_dir_hint_world;
    // Cache of last successful projection seed (u,v)
    mutable std::optional<std::pair<T,T>> last_seed_uv;

    // --- Curvature proxy and gradient helpers (tolerance-aware) ---
    // Dimension-agnostic "curvature" proxy: sine of angle between Su and Sv (area of parallelogram normalized)
    T curvatureProxy(T u, T v) const {
        auto [Su, Sv] = evaluatePartials(u, v);
        T su2 = Su.squaredNorm();
        T sv2 = Sv.squaredNorm();
        T dot = Su.dot(Sv);
        T area2 = std::max<T>(0, su2 * sv2 - dot * dot); // = ||Su x Sv||^2 in 3D; Gram determinant in ND
        T denom = std::sqrt(std::max<T>(su2, T(1e-30))) * std::sqrt(std::max<T>(sv2, T(1e-30)));
        if (denom <= T(0)) return T(0);
        return std::sqrt(area2) / denom; // in [0,1]
    }

    // Approximate |∇kappa| via central differences in param-space
    T curvatureGradientMagnitude(T u, T v, T du, T dv) const {
        const T epsu = std::max<T>(std::abs(du), T(1e-8));
        const T epsv = std::max<T>(std::abs(dv), T(1e-8));
        T ku_p = curvatureProxy(u + epsu, v);
        T ku_m = curvatureProxy(u - epsu, v);
        T kv_p = curvatureProxy(u, v + epsv);
        T kv_m = curvatureProxy(u, v - epsv);
        T gu = std::abs(ku_p - ku_m) / (T(2) * epsu);
        T gv = std::abs(kv_p - kv_m) / (T(2) * epsv);
        return gu + gv; // L1 blend; robust and cheap
    }

public:
    Surface(const std::function<PointType(T,T)>& f,
            const std::pair<std::pair<T,T>, std::pair<T,T>>& domain,
            std::function<T(const Eigen::Matrix<T, N, 1>&)> sdf = nullptr)
        : surfaceFunc(f), analyticSDF(sdf), surfaceDomain(domain) {}

    // Optional: provide a world-space direction hint to bias BVH refinement along that direction
    void setBVHDirectionHint(const Eigen::Matrix<T, N, 1>& dir_world) const {
        if (dir_world.norm() > T(0)) {
            bvh_dir_hint_world = dir_world.normalized();
            bvh_valid = false; // force rebuild on next getBVH()
        }
    }

    // Clear any previously set BVH direction hint
    void clearBVHDirectionHint() const {
        bvh_dir_hint_world.reset();
        bvh_valid = false;
    }

    // Evaluate at (u,v)
    PointType evaluate(T u, T v) const {
        return surfaceFunc(u, v);
    }

    // Evaluate (Su, Sv) partial derivatives at (u,v)
    std::pair<Eigen::Matrix<T, N, 1>, Eigen::Matrix<T, N, 1>> evaluatePartials(T u, T v) const {
        if constexpr (std::is_member_function_pointer_v<decltype(&Surface<T, N>::analyticalPartials)>) {
            // If analyticalPartials is overridden, use it
            return analyticalPartials(u, v);
        } else {
            // Fallback to numerical
            return numericalPartials(u, v);
        }
    }

    // Analytical partial derivatives (override in derived classes for analytical surfaces)
    virtual std::pair<Eigen::Matrix<T, N, 1>, Eigen::Matrix<T, N, 1>> analyticalPartials(T u, T v) const {
        // Default: not implemented, fallback to numerical
        return numericalPartials(u, v);
    }

    // Numerical partial derivatives using central finite differences
    std::pair<Eigen::Matrix<T, N, 1>, Eigen::Matrix<T, N, 1>> numericalPartials(T u, T v, T h = T(1e-6)) const {
        // Central finite difference for each partial
        auto Su = (evaluate(u + h, v).coords - evaluate(u - h, v).coords) / (T(2) * h);
        auto Sv = (evaluate(u, v + h).coords - evaluate(u, v - h).coords) / (T(2) * h);
        return {Su, Sv};
    }

    // Evaluate surface normal at (u,v) (only for N==3)
    PointType evaluateNormal(T u, T v) const {
        static_assert(N == 3, "evaluateNormal only makes sense for 3D surfaces");
        auto [Su, Sv] = evaluatePartials(u, v);
        auto nvec = Su.cross(Sv);
        if (nvec.norm() == T(0)) {
            // Degenerate, return PointType(Eigen::Matrix<T, N, 1>::Zero());
            return PointType(Eigen::Matrix<T, N, 1>::Zero());
        }
        return PointType(nvec.normalized());
    }

    // Generalized principal normal in R^N using Gram–Schmidt and optional reference direction
    Eigen::Matrix<T, N, 1> principalNormal(T u, T v,
                                           const std::optional<Eigen::Matrix<T, N, 1>>& ref_dir_opt = std::nullopt) const
    {
        // Evaluate first-order partials
        auto [Su_raw, Sv_raw] = evaluatePartials(u, v);

        // Guard against degeneracy
        T su = Su_raw.norm();
        T sv = Sv_raw.norm();
        if (su < T(1e-16) && sv < T(1e-16)) {
            // Fallback: pick a canonical axis or ref_dir if provided
            Eigen::Matrix<T,N,1> n = ref_dir_opt ? *ref_dir_opt : Eigen::Matrix<T,N,1>::UnitX();
            if (n.norm() == T(0)) n = Eigen::Matrix<T,N,1>::UnitX();
            return n.normalized();
        }

        // Gram–Schmidt to build an orthonormal basis of the tangent plane span{Su,Sv}
        Eigen::Matrix<T,N,1> e1, e2;
        e1 = (su > T(1e-16)) ? (Su_raw / su) : Su_raw;
        // remove component of Sv along e1
        Eigen::Matrix<T,N,1> v2 = Sv_raw - (Sv_raw.dot(e1)) * e1;
        T v2n = v2.norm();
        if (v2n < T(1e-16)) {
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
            arbit(best_k) = T(1);
            v2 = arbit - (arbit.dot(e1)) * e1;
            v2n = v2.norm();
            if (v2n < T(1e-16)) v2 = Eigen::Matrix<T,N,1>::UnitY(); // final fallback where available
            v2n = v2.norm();
        }
        e2 = (v2n > T(1e-16)) ? (v2 / v2n) : v2;

        if constexpr (N == 3) {
            auto n3 = e1.cross(e2);
            T nn = n3.norm();
            if (nn <= T(0)) {
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
            if (nn < T(1e-16)) {
                // Try canonical axes to find a stable orthogonal direction
                T best_norm = T(0);
                Eigen::Matrix<T,N,1> best = Eigen::Matrix<T,N,1>::Zero();
                for (int k=0; k<N; ++k) {
                    Eigen::Matrix<T,N,1> axis = Eigen::Matrix<T,N,1>::Zero();
                    axis(k) = T(1);
                    Eigen::Matrix<T,N,1> cand = axis - (axis.dot(e1))*e1 - (axis.dot(e2))*e2;
                    T cn = cand.norm();
                    if (cn > best_norm) { best_norm = cn; best = cand; }
                }
                if (best_norm > T(0)) n = best; else n = pref; // ultimate fallback
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
            if (rn <= tol.evaluateEpsilon(std::max<T>(rn, T(1)))) {
                r.ok = true; r.signed_distance = rn;
                r.surface_point = PointType(S);
                r.iterations = it;
                // Signed distance via principal normal aligned with any BVH hint
                Eigen::Matrix<T,N,1> n = principalNormal(r.u, r.v, bvh_dir_hint_world);
                if (n.norm() > T(0)) {
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
            // Solve JTJ * delta = -JTr  (with tiny Tikhonov for robustness)
            T lam = std::max<T>(tol.evaluateEpsilon(T(1)), T(1e-12));
            JTJ(0,0) += lam; JTJ(1,1) += lam;
            T det = JTJ(0,0)*JTJ(1,1) - JTJ(0,1)*JTJ(1,0);
            if (std::abs(det) < T(1e-30)) break;
            Eigen::Matrix<T,2,1> delta;
            delta(0) = (-JTr(0)*JTJ(1,1) + JTr(1)*JTJ(0,1)) / det;
            delta(1) = (-JTr(1)*JTJ(0,0) + JTr(0)*JTJ(1,0)) / det;
            // Damped step
            T step_scale = T(1);
            // Backtracking if residual grows
            for (int bs=0; bs<5; ++bs) {
                T nu = r.u - step_scale*delta(0);
                T nv = r.v - step_scale*delta(1);
                auto Stry = evaluate(nu, nv).coords;
                T rn_try = (Stry - P).norm();
                if (rn_try < rn) { r.u = nu; r.v = nv; rn = rn_try; break; }
                step_scale *= T(0.5);
            }
            // Convergence checks
            if (std::abs(prev_norm - rn) <= tol.evaluateEpsilon(std::max<T>(rn, T(1)))) {
                prev_norm = rn; r.iterations = it+1;
                r.ok = true; r.signed_distance = rn;
                r.surface_point = PointType(evaluate(r.u, r.v).coords);
                Eigen::Matrix<T,N,1> n = principalNormal(r.u, r.v, bvh_dir_hint_world);
                if (n.norm() > T(0)) {
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
            Eigen::Matrix<T,N,1> n = principalNormal(r.u, r.v, bvh_dir_hint_world);
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
                                 int max_iters = 20,
                                 int max_seeds = 6) const
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
        auto seeds = gatherSeedsFromBVHNearPoint(P, tol, max_seeds);
        for (auto [u0, v0] : seeds) {
            // Skip near-duplicate seeds relative to the warm-start
            if (best.ok && std::abs(u0 - best.u) < tol.evaluateEpsilon(T(1)) &&
                           std::abs(v0 - best.v) < tol.evaluateEpsilon(T(1))) {
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
        Eigen::Matrix<T,N,1> n = principalNormal(r.u, r.v, bvh_dir_hint_world);
        // Enforce normal alignment with BVH direction hint if present
        if (bvh_dir_hint_world && n.dot(*bvh_dir_hint_world) < T(0)) {
            n = -n;
        }
        if (n.norm() > T(0)) {
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

    // Sweep a cross-section along a spine curve to create a surface
    template<typename Curve, typename CrossSectionFunc>
    [[deprecated("sweepSurface is deprecated. Use B-rep extrude/sweep instead.")]]
    static Surface<T, 3> sweepSurface(
        const Curve& spineCurve,
        const CrossSectionFunc& crossSectionFunc,
        std::pair<std::pair<T, T>, std::pair<T, T>> domain = {})
    {
        T umin = spineCurve.domain().first;
        T umax = spineCurve.domain().second;

        if (domain.first.first == T() && domain.first.second == T() &&
            domain.second.first == T() && domain.second.second == T()) {
            domain = {{umin, T(0)}, {umax, T(1)}};
        }

        // --- Rotation-Minimizing Frame (RMF) ---
        int samples = 500;
        std::vector<T> uSamples(samples);
        std::vector<Eigen::Matrix<T,3,1>> spinePoints(samples);
        for(int i=0;i<samples;++i){
            T uVal = umin + (umax-umin)*T(i)/(samples-1);
            uSamples[i] = uVal;
            spinePoints[i] = spineCurve.evaluate(uVal).coords;
        }

        std::vector<Eigen::Matrix<T,3,1>> tangents(samples);
        std::vector<Eigen::Matrix<T,3,1>> normals(samples);
        std::vector<Eigen::Matrix<T,3,1>> binormals(samples);

        // Compute tangents
        for(int i=0;i<samples;++i){
            if(i<samples-1)
                tangents[i] = (spinePoints[i+1]-spinePoints[i]).normalized();
            else
                tangents[i] = (spinePoints[i]-spinePoints[i-1]).normalized();
        }

        // Initialize first normal
        Eigen::Matrix<T,3,1> arbitrary = Eigen::Matrix<T,3,1>::UnitY();
        if(std::abs(tangents[0].dot(arbitrary)) > T(0.9)) arbitrary = Eigen::Matrix<T,3,1>::UnitX();
        normals[0] = (arbitrary - tangents[0]*(tangents[0].dot(arbitrary))).normalized();
        binormals[0] = tangents[0].cross(normals[0]).normalized();

        // RMF propagation
        for(int i=1;i<samples;++i){
            Eigen::Matrix<T,3,1> vPrev = tangents[i-1];
            Eigen::Matrix<T,3,1> vCurr = tangents[i];
            Eigen::Matrix<T,3,1> axis = vPrev.cross(vCurr);
            T axisNorm = axis.norm();
            Euclid::Tolerance tol;
            if(Euclid::equalWithinTolerance(axisNorm, T(0), tol, axisNorm)){
                normals[i] = normals[i-1];
                binormals[i] = binormals[i-1];
                continue;
            }
            axis /= axisNorm;
            T angle = std::acos(std::clamp(vPrev.dot(vCurr), T(-1), T(1)));
            auto rotate = [&](const Eigen::Matrix<T,3,1>& vec){
                Eigen::Matrix<T,3,1> rotated = vec*std::cos(angle) + axis.cross(vec)*std::sin(angle) + axis*(axis.dot(vec))*(1-std::cos(angle));
                return rotated.normalized();
            };
            normals[i] = rotate(normals[i-1]);
            binormals[i] = rotate(binormals[i-1]);
        }

        // Frame interpolation function
        auto frameAt = [=](T uVal){
            T tNorm = (uVal-umin)/(umax-umin);
            int idx = std::min(int(tNorm*(samples-1)), samples-1);
            return std::make_tuple(tangents[idx], normals[idx], binormals[idx]);
        };

        auto surfFunc = [spineCurve, crossSectionFunc, frameAt](T uVal, T vVal){
            auto [tangent, normal, binormal] = frameAt(uVal);
            Eigen::Matrix<T,3,1> center = spineCurve.evaluate(uVal).coords;

            Point<T,2> cs = crossSectionFunc(vVal);

            Eigen::Matrix<T,3,1> offset = normal*cs.coords(0) + binormal*cs.coords(1);
            return Point<T,3>(center + offset);
        };

        return Surface<T,3>(surfFunc, domain);
    }

    // Helper: Adaptive recursive patch subdivision for BVH (tolerance-aware, with optional direction hint)
    void subdividePatch(T umin, T umax, T vmin, T vmax,
                        const Euclid::Tolerance& tol, T world_scale,
                        T min_param_extent, int max_depth, int depth,
                        bool have_dir, const Eigen::Matrix<T, N, 1>& dir_world) const
    {
        // Compute patch center and extents
        T uc = (umin + umax) / T(2);
        T vc = (vmin + vmax) / T(2);
        T du = std::abs(umax - umin);
        T dv = std::abs(vmax - vmin);

        // Global cap to prevent exponential explosion in highly oscillatory regions
        const size_t MAX_PATCHES = 4096;
        if (bvh_patches.size() >= MAX_PATCHES) {
            std::array<PointType, 5> samples = {
                evaluate(umin, vmin),
                evaluate(umax, vmin),
                evaluate(umin, vmax),
                evaluate(umax, vmax),
                evaluate(uc,   vc  )
            };
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

        // World-space flatness test using a 5-point bbox
        std::array<PointType, 5> samples = {
            evaluate(umin, vmin),
            evaluate(umax, vmin),
            evaluate(umin, vmax),
            evaluate(umax, vmax),
            evaluate(uc,   vc  )
        };
        PointType minPt = samples[0];
        PointType maxPt = samples[0];
        for (const auto& p : samples) {
            for (int k = 0; k < N; ++k) {
                if (p.coords(k) < minPt.coords(k)) minPt.coords(k) = p.coords(k);
                if (p.coords(k) > maxPt.coords(k)) maxPt.coords(k) = p.coords(k);
            }
        }
        T bbox_diag = (maxPt.coords - minPt.coords).norm();
        T world_eps = tol.evaluateEpsilon(world_scale > T(0) ? world_scale : T(1));

        // Param-step for gradient probing: half patch, but not smaller than tolerance-induced param eps
        T param_eps = std::max(min_param_extent * T(0.5), T(1e-8));
        T step_u = std::max<T>(du * T(0.5), param_eps);
        T step_v = std::max<T>(dv * T(0.5), param_eps);

        // Curvature gradient magnitude at center
        T G = curvatureGradientMagnitude(uc, vc, step_u, step_v);

        // Stopping criteria
        bool depthStop = (depth >= max_depth);
        bool paramStop = (du <= min_param_extent && dv <= min_param_extent);

        // If world variation is tiny AND curvature does not change meaningfully across the patch, stop.
        //T param_eps_flat = std::max(min_param_extent * T(0.5), T(1e-8));
        bool curvatureFlat = (G * (du + dv)) <= tol.evaluateEpsilon(T(1));

        if (depthStop || paramStop || (bbox_diag <= world_eps && curvatureFlat)) {
            bvh_patches.push_back(SurfacePatchAABB{umin, umax, vmin, vmax, minPt, maxPt});
            return;
        }

        // --- Directional splitting logic (bias toward the projected line direction) ---
        bool split_u = true;
        bool split_v = true;

        if (have_dir) {
            // Project world direction onto param space using J = [Su Sv] at the patch center
            auto [Su, Sv] = evaluatePartials(uc, vc);
            // Solve min_{a,b} || Su*a + Sv*b - dir_world || using normal equations
            Eigen::Matrix<T,2,2> JTJ;
            JTJ(0,0) = Su.dot(Su);
            JTJ(0,1) = Su.dot(Sv);
            JTJ(1,0) = JTJ(0,1);
            JTJ(1,1) = Sv.dot(Sv);
            Eigen::Matrix<T,2,1> JTd;
            JTd(0) = Su.dot(dir_world);
            JTd(1) = Sv.dot(dir_world);

            T det = JTJ(0,0)*JTJ(1,1) - JTJ(0,1)*JTJ(1,0);
            if (std::abs(det) > T(1e-20)) {
                Eigen::Matrix<T,2,1> ab;
                // Inverse of 2x2
                ab(0) = ( JTJ(1,1)*JTd(0) - JTJ(0,1)*JTd(1)) / det; // a
                ab(1) = (-JTJ(1,0)*JTd(0) + JTJ(0,0)*JTd(1)) / det; // b
                T au = std::abs(ab(0));
                T bv = std::abs(ab(1));

                // If the parametric direction is strongly anisotropic, split only along dominant axis
                T ratio = (std::max(au, bv)) / (std::max(T(1e-20), std::min(au, bv)));
                if (ratio > T(2)) {
                    if (au > bv) { split_u = true; split_v = false; }
                    else         { split_u = false; split_v = true; }
                }
            }
        }

        // Otherwise, subdivide accordingly
        T uc_mid = (umin + umax) / 2;
        T vc_mid = (vmin + vmax) / 2;

        if (split_u && split_v) {
            // 4-way split (default)
            subdividePatch(umin, uc_mid, vmin, vc_mid, tol, world_scale, min_param_extent, max_depth, depth + 1, have_dir, dir_world);
            subdividePatch(uc_mid, umax, vmin, vc_mid, tol, world_scale, min_param_extent, max_depth, depth + 1, have_dir, dir_world);
            subdividePatch(umin, uc_mid, vc_mid, vmax, tol, world_scale, min_param_extent, max_depth, depth + 1, have_dir, dir_world);
            subdividePatch(uc_mid, umax, vc_mid, vmax, tol, world_scale, min_param_extent, max_depth, depth + 1, have_dir, dir_world);
        } else if (split_u) {
            // 2-way split along u
            subdividePatch(umin, uc_mid, vmin, vmax, tol, world_scale, min_param_extent, max_depth, depth + 1, have_dir, dir_world);
            subdividePatch(uc_mid, umax, vmin, vmax, tol, world_scale, min_param_extent, max_depth, depth + 1, have_dir, dir_world);
        } else { // split_v only
            subdividePatch(umin, umax, vmin, vc_mid, tol, world_scale, min_param_extent, max_depth, depth + 1, have_dir, dir_world);
            subdividePatch(umin, umax, vc_mid, vmax, tol, world_scale, min_param_extent, max_depth, depth + 1, have_dir, dir_world);
        }
    }

    // Adaptive recursive BVH generation (tolerance-driven)
    void generateBVH(int grid, const Euclid::Tolerance& tol, T world_scale = T(0), int max_depth = -1) const {
        bvh_patches.clear();
        // Domain
        T umin = surfaceDomain.first.first;
        T umax = surfaceDomain.first.second;
        T vmin = surfaceDomain.second.first;
        T vmax = surfaceDomain.second.second;

        // Estimate world scale if not provided: bbox diagonal of corners + center
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
        if (world_scale > T(0)) est_scale = world_scale; // caller override

        // Parametric step lower bound derived from world epsilon and Jacobian magnitude at center
        T uc = (umin+umax)/T(2), vc = (vmin+vmax)/T(2);
        auto [Su0,Sv0] = evaluatePartials(uc, vc);
        T J = std::max(Su0.norm(), Sv0.norm());
        if (J <= T(0)) J = T(1);
        T world_eps = tol.evaluateEpsilon(est_scale);
        T param_eps_from_world = world_eps / J;

        // Grid-based cap in param space
        T grid_eps_u = std::abs(umax-umin) / T(std::max(1, grid));
        T grid_eps_v = std::abs(vmax-vmin) / T(std::max(1, grid));
        T min_param_extent = std::max(param_eps_from_world, T(0.25) * std::min(grid_eps_u, grid_eps_v));

        int effective_depth = max_depth;
        if (effective_depth <= 0) {
            // Derive depth from ratio of domain span to minimal param extent; clamp to [3, 18]
            T span = std::max(std::abs(umax - umin), std::abs(vmax - vmin));
            T ratio = std::max<T>(T(1), span / std::max<T>(min_param_extent, T(1e-8)));
            effective_depth = std::clamp<int>(int(std::ceil(std::log2(double(ratio)))), 3, 18);
        }

        // Pass optional direction hint (if any)
        bool have_dir = bvh_dir_hint_world.has_value();
        Eigen::Matrix<T,N,1> dir_world = have_dir ? *bvh_dir_hint_world : Eigen::Matrix<T,N,1>::Zero();

        subdividePatch(umin, umax, vmin, vmax, tol, est_scale, min_param_extent, effective_depth, 0, have_dir, dir_world);
        bvh_valid = true;
    }

    // Ray vs AABB slab test in N-D; returns (hit, tmin, tmax)
    std::tuple<bool, T, T> rayAABBInterval(const Eigen::Matrix<T, N, 1>& rayOrigin,
                                           const Eigen::Matrix<T, N, 1>& rayDir,
                                           const Eigen::Matrix<T, N, 1>& bmin,
                                           const Eigen::Matrix<T, N, 1>& bmax) const
    {
        T tmin = -std::numeric_limits<T>::infinity();
        T tmax =  std::numeric_limits<T>::infinity();

        for (int k = 0; k < N; ++k) {
            T o = rayOrigin(k);
            T d = rayDir(k);
            T invD;
            if (std::abs(d) < T(1e-30)) {
                // Ray parallel to slab: must be within bounds
                if (o < bmin(k) || o > bmax(k)) return {false, T(0), T(0)};
                else continue;
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
                             T world_scale_hint = T(1)) const
    {
        Eigen::Matrix<T, N, 1> pad = Eigen::Matrix<T, N, 1>::Constant(tol.evaluateEpsilon(world_scale_hint));
        auto [hit, t0, t1] = rayAABBInterval(rayOrigin, rayDir, bmin - pad, bmax + pad);
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
            if (nn > T(0)) n /= nn; else n = principalNormal(proj.u, proj.v, bvh_dir_hint_world);
        } else {
            n = principalNormal(proj.u, proj.v, bvh_dir_hint_world);
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
                                     T world_scale,
                                     int depth,
                                     int max_depth)
    const
    {
        // Hard cap to avoid explosions in pathological cases
        const size_t MAX_EXTRA = 8192;
        if (bvh_patches.size() >= MAX_EXTRA) return;

        // Quick reject: ray vs AABB
        Eigen::Matrix<T,N,1> bmin = patch.minPoint.coords;
        Eigen::Matrix<T,N,1> bmax = patch.maxPoint.coords;
        if (!rayIntersectsAABBND(rayOrigin, rayDir, bmin, bmax, tol, world_scale)) return;

        // Intersect parametric box size
        T du = std::abs(patch.umax - patch.umin);
        T dv = std::abs(patch.vmax - patch.vmin);

        // If the world bbox is already within epsilon, stop.
        T bbox_diag = (bmax - bmin).norm();
        T world_eps = tol.evaluateEpsilon(world_scale > T(0) ? world_scale : T(1));
        bool small_world = (bbox_diag <= world_eps);

        // Param-based lower bound derived from generateBVH heuristic
        T uc = (patch.umin + patch.umax) / T(2);
        T vc = (patch.vmin + patch.vmax) / T(2);
        auto [Su, Sv] = evaluatePartials(uc, vc);
        T J = std::max(Su.norm(), Sv.norm());
        if (J <= T(0)) J = T(1);
        T param_eps_from_world = world_eps / J;
        bool small_param = (std::max(du, dv) <= param_eps_from_world);

        if (depth >= max_depth || (small_world && small_param)) return;

        // Find t-interval of ray within the AABB
        auto [hit, t0, t1] = rayAABBInterval(rayOrigin, rayDir, bmin, bmax);
        if (!hit) return;

        // Sample signed distances to the surface at a few points on that segment (adaptive)
        auto sample_segment = [&](const std::vector<T>& sgrid, bool& sign_change_out, bool& near_zero_out){
            T dist_prev = std::numeric_limits<T>::quiet_NaN();
            sign_change_out = false; near_zero_out = false;
            for (T s : sgrid) {
                T t = t0 + (t1 - t0) * s;
                Eigen::Matrix<T,N,1> P = rayOrigin + t * rayDir;
                T d = rayPointSignedDistanceToSurface(P, tol);
                if (std::abs(d) <= world_eps) near_zero_out = true;
                if (!std::isnan(dist_prev) && d * dist_prev < T(0)) sign_change_out = true;
                dist_prev = d;
            }
        };

        bool sign_change = false, near_zero = false;
        // Coarse 3-sample pass
        sample_segment({ T(0.0), T(0.5), T(1.0) }, sign_change, near_zero);

        // If uncertain and bbox is not tiny, refine to 5 samples
        if (!sign_change && !near_zero && bbox_diag > world_eps * T(2)) {
            sample_segment({ T(0.0), T(0.25), T(0.5), T(0.75), T(1.0) }, sign_change, near_zero);
        }

        // Check local curvature magnitude to guide refinement
        T G = curvatureGradientMagnitude(uc, vc, du * T(0.5), dv * T(0.5));
        bool high_curvature = (G * (du + dv)) > tol.evaluateEpsilon(world_scale > T(0) ? world_scale : T(1));

        // If neither sign change, proximity, nor high curvature is observed, we can stop for this patch.
        if (!sign_change && !near_zero && !high_curvature) return;

        // Otherwise, subdivide the patch and recurse.
        T um = (patch.umin + patch.umax) / T(2);
        T vm = (patch.vmin + patch.vmax) / T(2);

        SurfacePatchAABB q00 = makePatchAABB(patch.umin, um, patch.vmin, vm);
        SurfacePatchAABB q10 = makePatchAABB(um, patch.umax, patch.vmin, vm);
        SurfacePatchAABB q01 = makePatchAABB(patch.umin, um, vm, patch.vmax);
        SurfacePatchAABB q11 = makePatchAABB(um, patch.umax, vm, patch.vmax);

        refinePatchNearRayRecursive(rayOrigin, rayDir, q00, tol, world_scale, depth+1, max_depth);
        refinePatchNearRayRecursive(rayOrigin, rayDir, q10, tol, world_scale, depth+1, max_depth);
        refinePatchNearRayRecursive(rayOrigin, rayDir, q01, tol, world_scale, depth+1, max_depth);
        refinePatchNearRayRecursive(rayOrigin, rayDir, q11, tol, world_scale, depth+1, max_depth);

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
        // Ensure BVH exists; reuse world-scale estimated during build
        if (!bvh_valid) {
            generateBVH(40, tol);
        }

        // Estimate global scale from current BVH bounds (quick pass)
        Eigen::Matrix<T,N,1> bbmin = bvh_patches.front().minPoint.coords;
        Eigen::Matrix<T,N,1> bbmax = bvh_patches.front().maxPoint.coords;
        for (const auto& p : bvh_patches) {
            for (int k=0;k<N;++k) {
                bbmin(k) = std::min(bbmin(k), p.minPoint.coords(k));
                bbmax(k) = std::max(bbmax(k), p.maxPoint.coords(k));
            }
        }
        T world_scale = (bbmax - bbmin).norm();
        if (world_scale <= T(0)) world_scale = T(1e-6);

        // Important: iterate over a snapshot of current size because vector grows during refinement
        size_t initial = bvh_patches.size();
        for (size_t i = 0; i < initial; ++i) {
            refinePatchNearRayRecursive(rayOrigin, rayDir, bvh_patches[i], tol, world_scale, /*depth*/0, max_depth);
        }
        // BVH remains valid; patches were refined in-place.
        return bvh_patches;
    }

    const std::vector<SurfacePatchAABB>& getBVH(int grid, const Euclid::Tolerance& tol, T world_scale = T(0), int max_depth = 9) const {
        if (!bvh_valid) {
            generateBVH(grid, tol, world_scale, max_depth);
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

            // Curvature/conditioning indicator at patch center
            auto [Su, Sv] = evaluatePartials(uc, vc);
            T su = Su.norm();
            T sv = Sv.norm();
            T mixed = std::abs(Su.dot(Sv));
            // A slightly stronger indicator than just su*sv
            T indicator = su * sv + mixed;

            // Choose refinement level based on indicator
            // 1 => only center/corners, 3 => 3x3, 5 => 5x5, 7 => 7x7
            int refine = 1;
            if (large_extent || indicator > T(10)) refine = std::max(refine, 3);
            if (indicator > T(40))               refine = std::max(refine, 5);
            if (indicator > T(160))              refine = std::max(refine, 7);

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

    TopologyCheckResult hasClosedTopology() const {
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

        // Helper lambda to check if face is degenerate (repeated vertices)
        auto isDegenerateFace = [&](const std::vector<size_t>& face) -> bool {
            std::set<size_t> uniqueVerts(face.begin(), face.end());
            if (uniqueVerts.size() < face.size()) return true;
            if constexpr (N == 3) {
                // For 3D triangles, also check zero area
                if (face.size() == 3) {
                    const auto& a = vertices[face[0]].coords;
                    const auto& b = vertices[face[1]].coords;
                    const auto& c = vertices[face[2]].coords;
                    auto ab = b - a;
                    auto ac = c - a;
                    auto cross = ab.cross(ac);
                    if (cross.norm() < 1e-12) return true;
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
