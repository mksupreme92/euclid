#pragma once

#include "../geometry.hpp"
#include "../tolerance.hpp"
#include "intersection_result.hpp"
#include <optional>
#include <vector>
#include <string>
#include <algorithm>

#include <tuple>
#include <iostream>

namespace Euclid::Geometry {

// ======================================================
// Ray–AABB intersection helper (used for BVH patch filtering)
// ======================================================
template <typename T, int N>
inline bool rayIntersectsAABB(const Eigen::Matrix<T, N, 1>& rayOrigin,
                              const Eigen::Matrix<T, N, 1>& rayDir,
                              const Eigen::Matrix<T, N, 1>& aabbMin,
                              const Eigen::Matrix<T, N, 1>& aabbMax) {
    T tmin = -std::numeric_limits<T>::infinity();
    T tmax =  std::numeric_limits<T>::infinity();
    for (int i = 0; i < N; ++i) {
        if (std::abs(rayDir[i]) < Tolerance().evaluateEpsilon(std::max(aabbMax[i] - aabbMin[i], (T)1))) {
            if (rayOrigin[i] < aabbMin[i] || rayOrigin[i] > aabbMax[i])
                return false; // Parallel and outside slab
        } else {
            T invD = T(1) / rayDir[i];
            T t0 = (aabbMin[i] - rayOrigin[i]) * invD;
            T t1 = (aabbMax[i] - rayOrigin[i]) * invD;
            if (t0 > t1) std::swap(t0, t1);
            tmin = std::max(tmin, t0);
            tmax = std::min(tmax, t1);
            if (tmax < tmin) return false;
        }
    }
    return tmax >= tmin;
}

// Ray–AABB interval helper (returns whether hit, and [t_enter, t_exit])
template <typename T, int N>
inline std::tuple<bool, T, T> rayAABBInterval(const Eigen::Matrix<T, N, 1>& rayOrigin,
                                              const Eigen::Matrix<T, N, 1>& rayDir,
                                              const Eigen::Matrix<T, N, 1>& aabbMin,
                                              const Eigen::Matrix<T, N, 1>& aabbMax,
                                              const Tolerance& tol,
                                              T world_scale)
{
    T tmin = -std::numeric_limits<T>::infinity();
    T tmax =  std::numeric_limits<T>::infinity();
    const T eps = tol.evaluateEpsilon(std::max(world_scale, T(1)));
    for (int i = 0; i < N; ++i) {
        const T d = rayDir[i];
        const T o = rayOrigin[i];
        const T mn = aabbMin[i];
        const T mx = aabbMax[i];
        if (std::abs(d) < eps) {
            // Ray parallel to slab; reject if origin not within
            if (o < mn - eps || o > mx + eps) return {false, T(0), T(0)};
            continue;
        }
        const T invD = T(1) / d;
        T t0 = (mn - o) * invD;
        T t1 = (mx - o) * invD;
        if (t0 > t1) std::swap(t0, t1);
        tmin = std::max(tmin, t0);
        tmax = std::min(tmax, t1);
        if (tmax < tmin) return {false, T(0), T(0)};
    }
    // Always return the interval, even if it is behind the origin (full line support)
    return {true, tmin, tmax};
}

// ======================================================
// Line Intersection Result
// ======================================================
template <typename T, int N>
struct LineIntersectionResult : public IntersectionResult<T, N> {
    std::vector<Point<T, N>> points; // intersection points (if any)
    std::vector<Line<T, N>> lines;   // intersection lines (if any)

    LineIntersectionResult() = default;

    LineIntersectionResult(bool intersects, const std::vector<Point<T, N>>& pts = {},
                           const std::vector<Line<T, N>>& lns = {}, std::string desc = "")
        : IntersectionResult<T, N>{intersects, desc}, points(pts), lines(lns) {}
};

// ======================================================
// Line–Line Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& l1, const Line<T, N>& l2,
                                      const Tolerance& tol = Tolerance()) {
    LineIntersectionResult<T, N> result;

    // Adaptive scale for tolerance model (mirrors line–plane scale logic)
    T scale = l1.direction().norm() + l2.direction().norm() +
              (l1.point1().coords - l2.point1().coords).norm();

    if constexpr (N == 2) {
        // 2D intersection using angular tolerance model (sinTheta)
        auto d1n = l1.direction().normalized();
        auto d2n = l2.direction().normalized();
        T angularTol = tol.evaluateEpsilon(std::max(scale, T(1)));
        T cosTheta = std::abs(d1n.dot(d2n));
        T sinTheta = std::sqrt(std::max(T(0), T(1) - cosTheta * cosTheta));
        if (sinTheta < angularTol) {
            // Lines are nearly parallel — check if coincident
            auto diff = l2.point1().coords - l1.point1().coords;
            T eps = tol.absTol + tol.evaluateEpsilon(scale);
            bool coincident = false;
            Eigen::Matrix<T, N, 1> proj = diff - diff.dot(d1n) * d1n;
            if (proj.norm() < eps) coincident = true;
            if (coincident) {
                result.intersects = true;
                result.lines = {l1};
                result.description = "Lines are coincident";
            } else {
                result.description = "Lines are parallel";
            }
            return result;
        }
        // Not parallel: solve for intersection using determinants
        T dx1 = l1.direction()[0], dy1 = l1.direction()[1];
        T dx2 = l2.direction()[0], dy2 = l2.direction()[1];
        auto diff = l2.point1().coords - l1.point1().coords;
        T det = dx1 * dy2 - dy1 * dx2;
        T t1 = (diff[0] * dy2 - diff[1] * dx2) / det;
        T t2 = (diff[0] * dy1 - diff[1] * dx1) / det; // param on l2
        Point<T, N> pt(l1.point1().coords + t1 * l1.direction());
        result.intersects = true;
        result.points = {pt};
        result.description = "Lines intersect at a point";
        ParamHit<T, N> hit;
        hit.t_line = t1;
        hit.u_curve = t2;
        hit.v_surface = T(0);
        hit.p = pt;
        hit.tangential = false;
        result.addHit(hit);
        result.hits.push_back(hit);
        return result;
    }

    if constexpr (N == 3) {
        // 3D intersection using angular tolerance model (sinTheta)
        auto d1n = l1.direction().normalized();
        auto d2n = l2.direction().normalized();
        T angularTol = tol.evaluateEpsilon(std::max(scale, T(1)));
        T cosTheta = std::abs(d1n.dot(d2n));
        T sinTheta = std::sqrt(std::max(T(0), T(1) - cosTheta * cosTheta));
        if (sinTheta < angularTol) {
            // Lines are nearly parallel — check if coincident
            auto diff = l2.point1().coords - l1.point1().coords;
            T eps = tol.absTol + tol.evaluateEpsilon(scale);
            bool coincident = false;
            Eigen::Matrix<T, N, 1> proj = diff - diff.dot(d1n) * d1n;
            if (proj.norm() < eps) coincident = true;
            if (coincident) {
                result.intersects = true;
                result.lines = {l1};
                result.description = "Lines are coincident";
            } else {
                result.description = "Lines are parallel";
            }
            return result;
        }
        // Not parallel: solve for intersection using closest point method
        auto p1 = l1.point1().coords;
        auto d1 = l1.direction();
        auto p2 = l2.point1().coords;
        auto d2 = l2.direction();
        auto r = p1 - p2;
        T a = d1.dot(d1), b = d1.dot(d2), c = d2.dot(d2);
        T d = d1.dot(r), e = d2.dot(r);
        T denom = a * c - b * b;
        T t1 = (b * e - c * d) / denom;
        T t2 = (a * e - b * d) / denom;
        Point<T, N> pt1(p1 + t1 * d1);
        Point<T, N> pt2(p2 + t2 * d2);
        T dist = (pt1 - pt2).norm();
        if (dist > tol.absTol + tol.evaluateEpsilon(scale)) {
            // Nearly skew lines: stabilize using midpoint
            Point<T, N> mid((pt1.coords + pt2.coords) * T(0.5));
            if (dist < tol.absTol * T(10) + tol.evaluateEpsilon(scale)) {
                result.intersects = true;
                result.points = {mid};
                result.description = "Lines nearly intersect (stabilized midpoint)";
                ParamHit<T, N> hit;
                hit.t_line = t1;
                hit.u_curve = t2;
                hit.v_surface = T(0);
                hit.p = mid;
                hit.tangential = false;
                result.addHit(hit);
                result.hits.push_back(hit);
                return result;
            }
            result.description = "Lines are skew (do not intersect)";
            return result;
        }
        result.intersects = true;
        Point<T, N> mid((pt1.coords + pt2.coords) * T(0.5));
        result.points = {mid};
        result.description = "Lines intersect at a point";
        ParamHit<T, N> hit;
        hit.t_line = t1;
        hit.u_curve = t2;
        hit.v_surface = T(0);
        hit.p = mid;
        hit.tangential = false;
        result.addHit(hit);
        result.hits.push_back(hit);
        return result;
    }

    // Generic ND implementation: for N != 2 and N != 3
    if constexpr (N != 2 && N != 3) {
        // Check if lines are parallel or coincident
        auto p1 = l1.point1().coords;
        auto d1 = l1.direction().normalized();
        auto p2 = l2.point1().coords;
        auto d2 = l2.direction().normalized();
        auto diff = p2 - p1;

        T dot_dir = d1.dot(d2);
        T eps = tol.absTol + tol.evaluateEpsilon(scale);

        if (std::abs(std::abs(dot_dir) - 1) < eps) {
            // Lines are parallel, check if coincident
            T dist = diff.norm();
            if (dist < eps) {
                // Lines are coincident
                result.intersects = true;
                result.lines = {l1};
                result.description = "Lines are coincident";
            } else {
                // Check if diff is along the direction line (collinear)
                T proj = diff.dot(d1);
                Point<T, N> projected_point = Point<T, N>(p1 + proj * d1);
                if ((projected_point.coords - p2).norm() < eps) {
                    result.intersects = true;
                    result.lines = {l1};
                    result.description = "Lines are coincident";
                } else {
                    result.description = "Lines are parallel";
                }
            }
            return result;
        }

        // For ND, solve for intersection parameters if possible
        // Using least squares to solve for parameters t1 and t2 where:
        // p1 + t1*d1 = p2 + t2*d2

        // Construct system: t1*d1 - t2*d2 = p2 - p1
        // Rewrite as A * [t1; t2] = b, where A = [d1, -d2], b = diff

        Eigen::Matrix<T, N, 2> A;
        for (int i = 0; i < N; ++i) {
            A(i, 0) = d1[i];
            A(i, 1) = -d2[i];
        }
        Eigen::Matrix<T, N, 1> b;
        for (int i = 0; i < N; ++i) {
            b(i, 0) = diff[i];
        }

        // Solve least squares
        Eigen::Matrix<T, 2, 1> t = A.colPivHouseholderQr().solve(b);

        Point<T, N> pt1 = Point<T, N>(p1 + t(0, 0) * d1);
        Point<T, N> pt2 = Point<T, N>(p2 + t(1, 0) * d2);
        T dist = (pt1 - pt2).norm();

        if (dist > tol.absTol + tol.evaluateEpsilon(scale)) {
            result.description = "Lines do not intersect";
            return result;
        }

        result.intersects = true;
        result.points = {pt1};
        result.description = "Lines intersect at a point";
        return result;
    }
}

// ======================================================
// Line–Curve Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& line, const Curve<T, N>& c,
                                       const Tolerance& tol = Tolerance()) {
    LineIntersectionResult<T, N> result;

    // Domain of curve
    auto domain = c.domain();
    T u0 = domain.first;
    T u1 = domain.second;
    if (u1 <= u0) { result.description = "Empty curve domain"; return result; }

    const auto& P0  = line.point1().coords;
    const auto  dir = line.direction();
    const T dir_norm2 = dir.squaredNorm();
    if (dir_norm2 < tol.evaluateEpsilon(std::max(dir.norm(), T(1)))) {
        result.description = "Line direction too small (degenerate)";
        return result;
    }

    // ----------------------------------------------
    // Build a unit vector n orthogonal to dir for signed distance s(u) = n·(C(u)-P0)
    // ----------------------------------------------
    Eigen::Matrix<T, N, 1> e = dir / std::sqrt(dir_norm2);
    Eigen::Matrix<T, N, 1> a = Eigen::Matrix<T, N, 1>::Zero();
    // pick an axis least aligned with e
    int bestAxis = 0; T minAbs = std::numeric_limits<T>::infinity();
    for (int i = 0; i < N; ++i) {
        T val = std::abs(e[i]);
        if (val < minAbs) { minAbs = val; bestAxis = i; }
    }
    a[bestAxis] = T(1);
    Eigen::Matrix<T, N, 1> n = a - (a.dot(e)) * e; // Gram-Schmidt
    T n_norm = n.norm();
    if (n_norm < std::numeric_limits<T>::epsilon()) {
        // fallback: choose next axis
        int alt = (bestAxis + 1) % N; a.setZero(); a[alt] = T(1);
        n = a - (a.dot(e)) * e; n_norm = n.norm();
    }
    if (n_norm == T(0)) { n = e; n_norm = T(1); } // degenerate fallback
    n /= n_norm;

    // ----------------------------------------------
    // Hybrid baseline + progressive curvature refinement
    // ----------------------------------------------
    T extent = std::abs(u1 - u0);

    // Efficient curvature-aware adaptive cap (bounded)
    T curvatureIntegral = std::clamp(c.evaluateIntegral(tol, 64), T(1e-6), T(1e3));

    // Hybrid curvature + tolerance scaling (Parasolid-style heuristic)
    int base_divs = std::max(64, int(std::min(
        extent / tol.paramTol / 4,
        std::clamp(512.0 * (1.0 + std::sqrt(curvatureIntegral)), 64.0, 8192.0)
    )));

    std::vector<T> samplePoints;
    samplePoints.reserve(base_divs + 1);
    for (int i = 0; i <= base_divs; ++i) {
        samplePoints.push_back(u0 + (u1 - u0) * (T(i) / base_divs));
    }

    // Curvature estimator for adaptive refinement
    auto curvatureAt = [&](T u) {
        const T h = std::max(tol.paramTol * extent, T(1e-6));
        auto C0 = c.evaluate(u).coords;
        auto Cp = c.evaluate(std::min(u1, u + h)).coords;
        auto Cm = c.evaluate(std::max(u0, u - h)).coords;
        Eigen::Matrix<T, N, 1> second = (Cp - 2 * C0 + Cm) / (h * h);
        return second.norm();
    };

    auto shouldRefine = [&](T a, T b) {
        T ka = curvatureAt(a);
        T kb = curvatureAt(b);
        T kmax = std::max(ka, kb);
        T segLen = std::abs(b - a);
        return (kmax * segLen) > (tol.paramTol * extent * tol.evalFactor);
    };

    // Progressive curvature-based refinement (3 passes)
    for (int pass = 0; pass < 3; ++pass) {
        std::vector<T> newPoints;
        newPoints.reserve(samplePoints.size() / 2);
        for (size_t i = 0; i + 1 < samplePoints.size(); ++i) {
            T a = samplePoints[i];
            T b = samplePoints[i + 1];
            if (shouldRefine(a, b)) {
                T mid = (a + b) * T(0.5);
                newPoints.push_back(mid);
            }
        }
        if (newPoints.empty()) break;
        samplePoints.insert(samplePoints.end(), newPoints.begin(), newPoints.end());
        std::sort(samplePoints.begin(), samplePoints.end());
        samplePoints.erase(std::unique(samplePoints.begin(), samplePoints.end(), [&](T x, T y){ return std::abs(x - y) < tol.paramTol; }), samplePoints.end());
    }

    // Compute signed distances s(u)
    std::vector<T> svals;
    svals.reserve(samplePoints.size());
    for (T uu : samplePoints) {
        svals.push_back(n.dot(c.evaluate(uu).coords - P0));
    }

    // Bracket intervals where sign changes or near-zero values occur (robust version)
    struct Bracket { T a, b; };
    std::vector<Bracket> brackets;
    T sTol = tol.evaluateEpsilon(std::max({extent, P0.norm(), dir.norm()}));
    for (size_t i = 0; i + 1 < samplePoints.size(); ++i) {
        T uA = samplePoints[i], uB = samplePoints[i + 1];
        T sA = svals[i], sB = svals[i + 1];
        bool crosses = (sA * sB < T(0));
        bool nearZero = (std::abs(sA) < sTol || std::abs(sB) < sTol);
        if (crosses || nearZero) {
            brackets.push_back({uA, uB});
        }
    }

    // --- NEW: collapse overlapping/nearby brackets by midpoint in u ---
    if (!brackets.empty()) {
        // Compute an estimate of the sampling step for clustering
        const T sample_step_u = (u1 - u0) / std::max<T>(T(1), static_cast<T>(samplePoints.size()));
        const T u_cluster_from_brackets = std::max(T(8) * tol.paramTol, T(2.5) * sample_step_u);

        // Sort by midpoint and keep only the first in each cluster
        struct MidB { T mid; T a; T b; };
        std::vector<MidB> mids; mids.reserve(brackets.size());
        for (const auto& br : brackets) {
            T mid = (br.a + br.b) * T(0.5);
            mids.push_back({mid, br.a, br.b});
        }
        std::sort(mids.begin(), mids.end(), [](const MidB& x, const MidB& y){ return x.mid < y.mid; });

        std::vector<Bracket> collapsed; collapsed.reserve(mids.size());
        for (const auto& mb : mids) {
            if (collapsed.empty()) {
                collapsed.push_back({mb.a, mb.b});
            } else {
                T prev_mid = (collapsed.back().a + collapsed.back().b) * T(0.5);
                if (std::abs(mb.mid - prev_mid) > u_cluster_from_brackets) {
                    collapsed.push_back({mb.a, mb.b});
                } else {
                    // Merge by enlarging the latest bracket for better Newton seeding stability
                    collapsed.back().a = std::min(collapsed.back().a, mb.a);
                    collapsed.back().b = std::max(collapsed.back().b, mb.b);
                }
            }
        }
        brackets.swap(collapsed);
    }

    if (brackets.empty()) {
        // fallback to midpoint if nothing found
        brackets.push_back({(u0 + u1) * 0.5, (u0 + u1) * 0.5});
    }

    // ----------------------------------------------
    // Refine each bracket via Gauss–Newton on (u,t)
    // ----------------------------------------------
    auto geometric_scale = [&](T u_seed){
        return std::max<T>({T(1), extent, c.evaluate(u_seed).coords.norm(), line.direction().norm()});
    };

    auto tangency_metric = [&](const Eigen::Matrix<T,N,1>& Cup) {
        // Gram determinant of [Cup, dir]
        T a2 = Cup.squaredNorm();
        T b2 = dir_norm2;
        T dot = Cup.dot(dir);
        T area2 = a2 * b2 - dot * dot; // squared area of parallelogram
        return std::sqrt(std::max<T>(area2, T(0)));
    };

    std::vector<ParamHit<T,N>> hits;
    std::vector<Point<T,N>>    pts; // for dedup and legacy points

    int max_iters = 25;
    for (const auto& br : brackets) {
        // Seed u at midpoint (or endpoint if degenerate)
        T u = (br.a == br.b) ? br.a : (br.a + br.b) * T(0.5);
        Eigen::Matrix<T,N,1> C0 = c.evaluate(u).coords;
        T t = dir.dot(C0 - P0) / dir_norm2;

        T prevDist = (C0 - (P0 + t * dir)).norm();
        T scale    = geometric_scale(u);

        for (int iter = 0; iter < max_iters; ++iter) {
            // Numerical derivative C'(u)
            T h = std::max(tol.paramTol, tol.evaluateEpsilon(std::abs(u) + T(1)));
            T up = std::min(u1, u + h);
            T um = std::max(u0, u - h);
            Eigen::Matrix<T,N,1> Cup = (c.evaluate(up).coords - c.evaluate(um).coords) / (up - um);

            Eigen::Matrix<T,N,1> Cu  = c.evaluate(u).coords;
            Eigen::Matrix<T,N,1> L   = P0 + t * dir;
            Eigen::Matrix<T,N,1> R   = Cu - L; // residual

            // J = [C'(u), -dir]
            Eigen::Matrix<T,N,2> J;
            for (int k = 0; k < N; ++k) { J(k,0) = Cup[k]; J(k,1) = -dir[k]; }

            Eigen::Matrix<T,2,2> JTJ = J.transpose() * J;
            Eigen::Matrix<T,2,1> JTR = J.transpose() * R;

            Eigen::Matrix<T,2,1> delta;
            Eigen::FullPivLU<Eigen::Matrix<T,2,2>> lu(JTJ);
            if (lu.isInvertible()) delta = lu.solve(-JTR);
            else {
                Eigen::Matrix<T,2,2> JTJ_d = JTJ;
                T damping = tol.evaluateEpsilon(std::max(scale, T(1)));
                JTJ_d(0,0) += damping; JTJ_d(1,1) += damping;
                delta = JTJ_d.colPivHouseholderQr().solve(-JTR);
            }

            T du = delta(0,0);
            T dt = delta(1,0);

            // Simple adaptive damping
            T alpha = (iter < 5) ? T(0.7) : T(0.9);
            u += alpha * du;
            t += alpha * dt;

            // Clamp u to domain
            if (u < u0) u = u0; if (u > u1) u = u1;

            Eigen::Matrix<T,N,1> Cn = c.evaluate(u).coords;
            Eigen::Matrix<T,N,1> Ln = P0 + t * dir;
            T dist = (Cn - Ln).norm();
            if (std::abs(dist - prevDist) < tol.evaluateEpsilon(std::max(dist, T(1)))) break;
            prevDist = dist;
        }

        // Acceptance test
        Eigen::Matrix<T,N,1> Cfin = c.evaluate(u).coords;
        Eigen::Matrix<T,N,1> Lfin = P0 + t * dir;
        T finalDist = (Cfin - Lfin).norm();
        T scale_fin = geometric_scale(u);
        T acceptTol = tol.absTol + tol.evaluateEpsilon(scale_fin) * T(5);
        if (finalDist <= acceptTol) {
            // Dedup by param and by world distance
            Point<T,N> Pfin((Cfin + Lfin) * T(0.5));
            bool duplicate = false;
            // Use tolerance-model–based scaling for deduplication
            T dedup_param = tol.paramTol * std::max<T>(T(2), std::cbrt(scale_fin));
            T dedup_world = tol.evaluateEpsilon(scale_fin) * std::max<T>(T(2), std::cbrt(scale_fin));
            for (size_t i = 0; i < hits.size(); ++i) {
                if (std::abs(hits[i].u_curve - u) <= dedup_param) { duplicate = true; break; }
                if ((pts[i].coords - Pfin.coords).norm() <= dedup_world) { duplicate = true; break; }
            }
            if (!duplicate) {
                // Tangency classification
                T h = std::max(tol.paramTol, tol.evaluateEpsilon(std::abs(u) + T(1)));
                T up = std::min(u1, u + h);
                T um = std::max(u0, u - h);
                Eigen::Matrix<T,N,1> Cup = (c.evaluate(up).coords - c.evaluate(um).coords) / (up - um);
                bool isTangential = tangency_metric(Cup) <= tol.evaluateEpsilon(scale_fin);

                ParamHit<T,N> hit; hit.t_line = t; hit.u_curve = u; hit.v_surface = T(0); hit.p = Pfin; hit.tangential = isTangential;
                // Fill both the new hits vector and legacy points for compatibility
                result.addHit(hit);           // fills base points
                hits.push_back(hit);
                pts.push_back(Pfin);
            }
        }
    }

    // ======================================================
    // Unified deduplication on hits (param-space + world-space)
    //   Rationale:
    //   - Many near-tangent configurations (e.g., line y=0 and parabola y=x^2)
    //     produce clusters of Newton solutions around the true root (u≈0).
    //   - Dedup must operate on "hits" (which carry u parameters), not only on
    //     result.points, otherwise later syncing would reintroduce duplicates.
    // ======================================================
    if (!hits.empty()) {
        // Estimate average sampling step to scale param clustering
        T sample_step_u = (u1 - u0) / std::max<T>(T(1), static_cast<T>(samplePoints.size()));

        // Param-space and world-space clustering thresholds (more conservative)
        const T u_cluster_tol = std::max(T(8) * tol.paramTol, sample_step_u * T(3.0));
        const T world_cluster_tol = std::max(
            tol.absTol * T(100),
            tol.evaluateEpsilon(std::max<T>({T(1), dir.norm(), P0.norm()})) * T(75)
        );

        std::vector<ParamHit<T,N>> filtered_hits;
        filtered_hits.reserve(hits.size());

        for (const auto& h : hits) {
            bool dup = false;
            for (auto& f : filtered_hits) {
                const T du = std::abs(h.u_curve - f.u_curve);
                const T dw = (h.p.coords - f.p.coords).norm();
                if (du <= u_cluster_tol || dw <= world_cluster_tol) {
                    // Prefer the representative closer to the line (smaller residual along Newton)
                    if ((h.p.coords - (P0 + h.t_line * dir)).norm() <=
                        (f.p.coords - (P0 + f.t_line * dir)).norm()) {
                        f = h; // replace with better representative
                    }
                    dup = true;
                    break;
                }
            }
            if (!dup) filtered_hits.push_back(h);
        }

        hits.swap(filtered_hits);
    }

    // Sync deduped hits back to result.points
    if (!hits.empty()) {
        result.points.clear();
        for (const auto& h : hits)
            result.points.push_back(h.p);
    }

    if (!hits.empty()) {
        result.intersects = true;
        if (hits.size() > 1) result.description = "Line intersects curve (multi-hit)";
        else                 result.description = "Line intersects curve (single hit)";
    } else {
        result.description = "Line does not intersect curve";
    }
    return result;
}


// ======================================================
// Public utility to precompute reusable (u,v) samples for repeated Line–Surface queries.
// ======================================================
template <typename T, int N>
std::vector<std::pair<T, T>> prepareSurfaceLineSamples(
    const Surface<T, N>& s,
    const Tolerance& tol = Tolerance(),
    int max_grid = -1)
{
    return s.generateCandidateSamples(tol, max_grid);
}

// Overload that uses precomputed samples (avoids rebuilding/visiting BVH for each line).
template <typename T, int N>
LineIntersectionResult<T, N> intersect(
    const Line<T, N>& line,
    const Surface<T, N>& s,
    const std::vector<std::pair<T, T>>& samples,
    const Tolerance& tol = Tolerance())
{
    LineIntersectionResult<T, N> result;
    if (samples.empty()) { result.description = "No candidate samples provided"; return result; }

    auto domain = s.surfaceDomain;
    T u0 = domain.first.first;
    T u1 = domain.first.second;
    T v0 = domain.second.first;
    T v1 = domain.second.second;
    if (u1 <= u0 || v1 <= v0) { result.description = "Empty surface domain"; return result; }

    // Precompute derivative steps once (cached inside Newton loop scope)
    T u_extent = std::abs(u1 - u0);
    T v_extent = std::abs(v1 - v0);
    T extent   = std::max(u_extent, v_extent);
    T eps      = tol.evaluateEpsilon(extent);

    const auto& P0  = line.point1().coords;
    const auto  dir = line.direction();
    T dir_norm2     = dir.squaredNorm();

    // Scale for tolerance checks
    T scale = 0;
    for (const auto& samp : samples) {
        Point<T, N> spt = s.evaluate(samp.first, samp.second);
        scale = std::max(scale, spt.coords.norm());
    }
    T coarse_thresh = std::max(tol.evaluateEpsilon(scale), tol.paramTol * extent);

    struct Candidate { T u, v, t, dist; };
    std::vector<Candidate> candidates;
    candidates.reserve(samples.size());

    T best_dist = std::numeric_limits<T>::infinity();
    size_t best_idx = 0;
    for (size_t i = 0; i < samples.size(); ++i) {
        T uu = samples[i].first;
        T vv = samples[i].second;
        Point<T, N> spt = s.evaluate(uu, vv);
        T proj_t = dir.dot(spt.coords - P0) / dir_norm2;
        Point<T, N> lpt(P0 + proj_t * dir);
        T dist = (spt - lpt).norm();
        if (dist <= coarse_thresh) {
            candidates.push_back({uu, vv, proj_t, dist});
        }
        if (dist < best_dist) { best_dist = dist; best_idx = i; }
    }
    if (candidates.empty()) {
        T uu = samples[best_idx].first;
        T vv = samples[best_idx].second;
        Point<T, N> spt = s.evaluate(uu, vv);
        T proj_t = dir.dot(spt.coords - P0) / dir_norm2;
        T dist = (spt - Point<T, N>(P0 + proj_t * dir)).norm();
        candidates.push_back({uu, vv, proj_t, dist});
    }

    // Cache derivative step sizes once (per call), not per Newton iteration
    const T uh = std::max(tol.paramTol * u_extent, std::numeric_limits<T>::epsilon() * (T)10);
    const T vh = std::max(tol.paramTol * v_extent, std::numeric_limits<T>::epsilon() * (T)10);

    std::vector<Eigen::Matrix<T, N, 1>> intersection_coords;
    intersection_coords.reserve(candidates.size());

    int max_iters = 30;
    T alpha = T(1.0);

    for (const auto& cand : candidates) {
        T u = cand.u, v = cand.v, t = cand.t;
        T prevDist = cand.dist;

        for (int iter = 0; iter < max_iters; ++iter) {
            Point<T, N> Spt = s.evaluate(u, v);

            // Use cached uh, vh — only clamp the evaluation points
            T up = std::min(u1, u + uh);
            T um = std::max(u0, u - uh);
            T vp = std::min(v1, v + vh);
            T vm = std::max(v0, v - vh);

            Point<T, N> Sup = s.evaluate(up, v);
            Point<T, N> Sum = s.evaluate(um, v);
            Point<T, N> Svp = s.evaluate(u, vp);
            Point<T, N> Svm = s.evaluate(u, vm);

            Eigen::Matrix<T, N, 1> Su = (Sup.coords - Sum.coords) / (up - um);
            Eigen::Matrix<T, N, 1> Sv = (Svp.coords - Svm.coords) / (vp - vm);

            Eigen::Matrix<T, N, 1> R;
            for (int k = 0; k < N; ++k) R(k, 0) = Spt.coords[k] - (P0[k] + t * dir[k]);

            Eigen::Matrix<T, N, 3> J;
            for (int k = 0; k < N; ++k) { J(k,0)=Su[k]; J(k,1)=Sv[k]; J(k,2)=-dir[k]; }

            Eigen::Matrix<T, 3, 1> delta;
            if constexpr (N == 3) {
                Eigen::FullPivLU<Eigen::Matrix<T, 3, 3>> lu(J);
                if (lu.isInvertible()) {
                    delta = lu.solve(-R);
                } else {
                    Eigen::Matrix<T, 3, 3> JTJ = J.transpose() * J;
                    Eigen::Matrix<T, 3, 1> JTR = J.transpose() * R;
                    T damping = tol.evaluateEpsilon(std::max(scale, (T)1));
                    for (int d = 0; d < 3; ++d) JTJ(d, d) += damping;
                    delta = JTJ.colPivHouseholderQr().solve(-JTR);
                }
            } else {
                delta = J.colPivHouseholderQr().solve(-R);
            }

            T du = delta(0,0), dv = delta(1,0), dt = delta(2,0);
            u += alpha * du; v += alpha * dv; t += alpha * dt;

            // clamp to domain
            if (u < u0) u = u0; if (u > u1) u = u1;
            if (v < v0) v = v0; if (v > v1) v = v1;

            Point<T, N> Snew = s.evaluate(u, v);
            Point<T, N> Lnew(P0 + t * dir);
            T dist = (Snew - Lnew).norm();
            if (std::abs(dist - prevDist) < eps * (T)0.5) break;
            prevDist = dist;
        }

        Point<T, N> Sfinal = s.evaluate(u, v);
        Point<T, N> Lfinal(P0 + t * dir);
        T finalDist = (Sfinal - Lfinal).norm();
        if (finalDist <= tol.evaluateEpsilon(std::max(scale, (T)1))) {
            intersection_coords.push_back((Sfinal.coords + Lfinal.coords) * (T)0.5);
        }
    }

    // Deduplicate
    std::vector<Point<T, N>> intersections;
    std::vector<bool> taken(intersection_coords.size(), false);
    T dedup_tol = tol.evaluateEpsilon(std::max(scale, (T)1));
    for (size_t i = 0; i < intersection_coords.size(); ++i) {
        if (taken[i]) continue;
        intersections.emplace_back(intersection_coords[i]);
        for (size_t j = i + 1; j < intersection_coords.size(); ++j) {
            if (!taken[j] && (intersection_coords[i] - intersection_coords[j]).norm() < dedup_tol) {
                taken[j] = true;
            }
        }
    }

    if (!intersections.empty()) {
        result.intersects = true;
        result.points = std::move(intersections);
        result.description = "Line intersects surface (precomputed BVH samples)";
    } else {
        result.description = "Line does not intersect surface";
    }
    return result;
}

// ======================================================
// Line–Surface Intersection (with BVH patch filtering)
// ======================================================
// ======================================================
// Line–Surface Intersection (adaptive ray-marching + Newton refinement)
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& line, const Surface<T, N>& s,
                                      const Tolerance& tol = Tolerance(), int max_grid = -1) {
    LineIntersectionResult<T, N> result;

    /*
    // [DEBUG] Entry parameters
    std::cout << "[DEBUG][Line–Surface] ENTER "
              << "origin=" << line.point1().coords.transpose()
              << " dir=" << line.direction().transpose()
              << " |dir|=" << line.direction().norm()
              << std::endl;
     */

    // Get some scale for tolerance
    const auto& P0 = line.point1().coords;
    const auto dir = line.direction();
    T dir_norm = dir.norm();
    if (dir_norm == 0) {
        result.description = "Line direction is zero";
        return result;
    }
    auto domain = s.surfaceDomain;
    T u0 = domain.first.first, u1 = domain.first.second;
    T v0 = domain.second.first, v1 = domain.second.second;
    T u_extent = std::abs(u1 - u0);
    T v_extent = std::abs(v1 - v0);

    /*
    std::cout << "[DEBUG][Line–Surface] domain u:[" << u0 << "," << u1 << "] v:[" << v0 << "," << v1 << "]"
              << " u_extent=" << u_extent << " v_extent=" << v_extent << std::endl;
     */

    // Estimate world bounding box to derive t-range
    std::array<std::pair<T,T>,5> uvProbe = {{{u0,v0},{u1,v0},{u0,v1},{u1,v1},{(u0+u1)/2,(v0+v1)/2}}};
    Eigen::Matrix<T,N,1> gmin = s.evaluate(uvProbe[0].first, uvProbe[0].second).coords;
    Eigen::Matrix<T,N,1> gmax = gmin;
    for (size_t i=1;i<uvProbe.size();++i){
        auto p = s.evaluate(uvProbe[i].first, uvProbe[i].second).coords;
        gmin = gmin.cwiseMin(p); gmax = gmax.cwiseMax(p);
    }
    T diag = (gmax - gmin).norm();

    //std::cout << "[DEBUG][Line–Surface] bbox diag=" << diag << std::endl;

    // Set up min step sizes based on tolerance model
    T min_step_u = tol.paramTol;
    T min_step_v = tol.paramTol;
    min_step_u = std::max(min_step_u, tol.evaluateEpsilon(std::max(u_extent, T(1))));
    min_step_v = std::max(min_step_v, tol.evaluateEpsilon(std::max(v_extent, T(1))));

    // Bounded grid seeding in (u,v) with optional refinement pass
    struct Candidate {T u,v,t,dist;};
    std::vector<Candidate> candidates;

    T curvature_factor = std::clamp(diag / (tol.absTol + tol.evaluateEpsilon(diag)), T(1), T(32));
    int grid = max_grid > 0 ? max_grid : std::clamp<int>(int(16 * curvature_factor), 16, 256);

    auto seedCandidates = [&](int g){
        for(int i=0;i<=g;++i){
            for(int j=0;j<=g;++j){
                T u = u0 + (u1-u0)*i/g;
                T v = v0 + (v1-v0)*j/g;
                Point<T,N> S = s.evaluate(u,v);
                T t = dir.dot(S.coords - P0)/dir.squaredNorm();
                Point<T,N> L(P0 + t*dir);
                T dist = (S - L).norm();
                if(dist < tol.absTol + tol.evaluateEpsilon(diag)*T(10)){
                    candidates.push_back({u,v,t,dist});
                }
            }
        }
    };

    seedCandidates(grid);
    /*
    std::cout << "[DEBUG][Line–Surface] seeded candidates (coarse)=" << candidates.size()
              << " grid=" << grid << std::endl;
     */

    // If we appear under-seeded for highly oscillatory geometry, add a second, denser pass.
    if ((int)candidates.size() < grid/2 || curvature_factor > T(8)) {
        int grid_refined = std::min(grid * 2, 256);
        seedCandidates(grid_refined);
        /*
        std::cout << "[DEBUG][Line–Surface] seeded candidates (refined)=" << candidates.size()
                  << " grid_refined=" << grid_refined << std::endl;
         */
    }

    // Ensure candidates are processed in a stable order by (u,v)
    if (!candidates.empty()) {
        std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b){
            if (a.u == b.u) return a.v < b.v; return a.u < b.u;
        });
    }

    // --------------------------------------------------
    // Multi‑resolution refinement: detect oscillatory residuals and re‑seed locally
    // --------------------------------------------------
    if (!candidates.empty() && curvature_factor > T(4)) {
        std::vector<Candidate> fine_candidates;
        fine_candidates.reserve(candidates.size() * 2);
        for (size_t i = 1; i < candidates.size(); ++i) {
            const auto& c0 = candidates[i-1];
            const auto& c1 = candidates[i];
            // Detect oscillatory residual: alternating distance slope sign
            T slope0 = (c1.dist - c0.dist);
            if (i+1 < candidates.size()) {
                T slope1 = (candidates[i+1].dist - c1.dist);
                T osc_thresh = std::max<T>(tol.absTol * T(5), tol.evaluateEpsilon(std::max(diag, T(1))));
                if (slope0 * slope1 < 0 && std::abs(slope0) > osc_thresh) {
                    T umid = (c0.u + c1.u) * T(0.5);
                    T vmid = (c0.v + c1.v) * T(0.5);
                    Point<T,N> S = s.evaluate(umid, vmid);
                    T t = dir.dot(S.coords - P0) / dir.squaredNorm();
                    Point<T,N> L(P0 + t * dir);
                    T dist = (S - L).norm();
                    if (dist < tol.absTol * T(50)) {
                        fine_candidates.push_back({umid, vmid, t, dist});
                    }
                }
            }
        }
        if (!fine_candidates.empty()) {
            candidates.insert(candidates.end(), fine_candidates.begin(), fine_candidates.end());
        }
    }

    // Param/world-space dedup containers and tolerances
    struct Accepted { T u, v, t; Eigen::Matrix<T,N,1> mid; };
    std::vector<Accepted> accepted;

    // Dedup radii: use both param-space and world-space criteria
    T dedup_world = std::max<T>(
        tol.evaluateEpsilon(std::max(diag, T(1))) * T(2),
        tol.absTol * T(20));
    T dedup_u = std::max<T>(tol.paramTol * T(4), tol.evaluateEpsilon(std::max(u_extent, T(1))) * T(10));
    T dedup_v = std::max<T>(tol.paramTol * T(4), tol.evaluateEpsilon(std::max(v_extent, T(1))) * T(10));

    // Cache derivative step sizes once (per call)
    const T uh = std::max(min_step_u * T(0.0625), std::numeric_limits<T>::epsilon() * 10);
    const T vh = std::max(min_step_v * T(0.0625), std::numeric_limits<T>::epsilon() * 10);

    std::vector<Point<T,N>> intersections;

    //std::cout << "[DEBUG][Line–Surface] processing candidates count=" << candidates.size() << std::endl;

    for(const auto& cand : candidates){
        T u = cand.u, v = cand.v, t = cand.t;
        T prevDist = cand.dist;
        int max_newton = (curvature_factor > T(8)) ? 40 : 25;
        T alpha = T(1.0);
        for(int iter=0; iter<max_newton; ++iter){
            Point<T,N> Spt = s.evaluate(u, v);
            T up = std::min(u1, u+uh);
            T um = std::max(u0, u-uh);
            T vp = std::min(v1, v+vh);
            T vm = std::max(v0, v-vh);
            Point<T,N> Sup = s.evaluate(up, v);
            Point<T,N> Sum = s.evaluate(um, v);
            Point<T,N> Svp = s.evaluate(u, vp);
            Point<T,N> Svm = s.evaluate(u, vm);
            Eigen::Matrix<T,N,1> Su = (Sup.coords - Sum.coords) / (up - um);
            Eigen::Matrix<T,N,1> Sv = (Svp.coords - Svm.coords) / (vp - vm);
            Eigen::Matrix<T,N,1> R = Spt.coords - (P0 + t*dir);
            Eigen::Matrix<T,N,3> J;
            for(int k=0;k<N;++k){ J(k,0)=Su[k]; J(k,1)=Sv[k]; J(k,2)=-dir[k]; }
            Eigen::Matrix<T,3,1> delta;
            if constexpr (N == 3) {
                Eigen::FullPivLU<Eigen::Matrix<T, 3, 3>> lu(J);
                if (lu.isInvertible()) {
                    delta = lu.solve(-R);
                } else {
                    Eigen::Matrix<T, 3, 3> JTJ = J.transpose() * J;
                    Eigen::Matrix<T, 3, 1> JTR = J.transpose() * R;
                    T damping = tol.evaluateEpsilon(std::max(diag, (T)1));
                    for (int d = 0; d < 3; ++d) JTJ(d, d) += damping;
                    delta = JTJ.colPivHouseholderQr().solve(-JTR);
                }
            } else {
                // Generic ND fallback
                delta = J.colPivHouseholderQr().solve(-R);
            }
            T du = delta(0,0), dv = delta(1,0), dt = delta(2,0);
            // Adaptive alpha update
            if (iter > 3 && prevDist > 0 && (Spt.coords - (P0 + t*dir)).norm() > prevDist * T(0.995))
                alpha *= T(0.85);
            else
                alpha = (iter<3) ? T(0.9) : T(0.97);
            u += alpha * du; v += alpha * dv; t += alpha * dt;
            u = std::clamp(u, u0, u1);
            v = std::clamp(v, v0, v1);
            Point<T,N> Snew = s.evaluate(u, v);
            Point<T,N> Lnew(P0 + t*dir);
            T dist = (Snew - Lnew).norm();
            bool converged = (std::abs(dist - prevDist) < tol.evaluateEpsilon(std::max(diag, T(1))) * T(0.5));
            prevDist = dist;
            if(converged) break;
        }
        Point<T,N> Sfinal = s.evaluate(u, v);
        Point<T,N> Lfinal(P0 + t*dir);
        T finalDist = (Sfinal - Lfinal).norm();
        T curvature_boost = std::min<T>(T(8) + diag * T(0.1), T(32));
        if(finalDist <= tol.absTol + tol.evaluateEpsilon(std::max(diag, T(1))) * curvature_boost){
            Point<T,N> mid((Sfinal.coords + Lfinal.coords) * T(0.5));

            // Param-space and world-space deduplication
            bool duplicate = false;
            for (const auto& a : accepted) {
                if (std::abs(a.u - u) <= dedup_u && std::abs(a.v - v) <= dedup_v) { duplicate = true; break; }
                if ((a.mid - mid.coords).norm() <= dedup_world) { duplicate = true; break; }
            }
            if (!duplicate) {
                accepted.push_back({u, v, t, mid.coords});
                intersections.push_back(mid);

                // Record detailed hit for downstream consumers (tests may ignore)
                ParamHit<T,N> hit; hit.t_line = t; hit.u_curve = u; hit.v_surface = v; hit.p = mid; hit.tangential = false;
                result.addHit(hit);
                /*
                std::cout << "[DEBUG][Line–Surface] accepted hit u=" << u << " v=" << v << " t=" << t
                          << " |mid|=" << mid.coords.norm() << std::endl;
                 */
            }
        }
    }

    //std::cout << "[DEBUG][Line–Surface] total accepted=" << intersections.size() << std::endl;

    if (!intersections.empty()) {
        result.intersects = true;
        result.points = intersections;
        result.description = "Line intersects surface (bounded grid Newton)";
        /*
        std::cout << "[DEBUG][Line–Surface] EXIT intersects=1 points=" << result.points.size()
                  << " desc=" << result.description << std::endl;
         */
    } else {
        result.description = "Line does not intersect surface (bounded grid Newton)";
        //std::cout << "[DEBUG][Line–Surface] EXIT intersects=0 points=0 desc=" << result.description << std::endl;
    }
    return result;
}

// ======================================================
// Line–Plane Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& line, const Plane<T, N>& plane,
                                      const Tolerance& tol = Tolerance()) {
    LineIntersectionResult<T, N> result;

    const auto& n = plane.normal;
    const auto& p0 = line.point1().coords;
    const auto& d = line.direction();
    const auto& b = plane.base.coords;

    T scale = d.norm() + n.norm();

    // Degenerate normal or direction check
    if (n.norm() < tol.evaluateEpsilon(scale) || d.norm() < tol.evaluateEpsilon(scale)) {
        result.description = "Degenerate plane or line direction";
        return result;
    }

    // Use tolerance model for angular threshold
    T angularTol = tol.absTol / std::max(scale, T(1));
    T distTol = tol.absTol * scale * tol.evalFactor * T(10);

    if constexpr (N == 3) {
        // 3D fast path: geometric intersection
        T denom = n.dot(d);
        if (std::abs(n.dot(d.normalized())) < angularTol) {
            // Near-parallel or coincident handling
            T dist = n.dot(p0 - b);
            // 1. Coincident case: explicit handling
            if (std::abs(dist) < tol.absTol + tol.evaluateEpsilon(scale)) {
                result.points.clear();
                result.lines = { line };
                result.intersects = true;
                result.description = "Line lies in plane (coincident)";
                return result;
            }
            // 2. Near-parallel case: accept only if both nearly parallel and within distance tolerance
            if (std::abs(dist) < distTol) {
                // Accept very shallow intersection within tolerance and proceed
            } else {
                result.description = "Line is nearly parallel but offset";
                result.intersects = false;
                return result;
            }
        }
        T t = -n.dot(p0 - b) / denom;
        Point<T, 3> pt(p0 + t * d);
        result.intersects = true;
        result.points = {pt};
        result.description = "Line intersects plane at a point (3D fast path)";
        return result;
    } else {
        // Generic ND implementation
        T denom = n.dot(d);
        if (std::abs(n.dot(d.normalized())) < angularTol) {
            // Near-parallel or coincident handling
            T dist = n.dot(p0 - b);
            // 1. Coincident case: explicit handling
            if (std::abs(dist) < tol.absTol + tol.evaluateEpsilon(scale)) {
                result.points.clear();
                result.lines = { line };
                result.intersects = true;
                result.description = "Line lies in plane (coincident)";
                return result;
            }
            // 2. Near-parallel case: accept only if both nearly parallel and within distance tolerance
            if (std::abs(dist) < distTol) {
                // Accept very shallow intersection within tolerance and proceed
            } else {
                result.description = "Line is nearly parallel but offset";
                result.intersects = false;
                return result;
            }
        }
        T t = -n.dot(p0 - b) / denom;
        Point<T, N> pt(p0 + t * d);
        result.intersects = true;
        result.points = {pt};
        result.description = "Line intersects plane at a point (ND fallback)";
        return result;
    }
}

// Reverse overload: Plane–Line
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Plane<T, N>& plane, const Line<T, N>& line,
                                      const Tolerance& tol = Tolerance()) {
    return intersect(line, plane, tol);
}

// ======================================================
// Line–Segment Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& line, const Segment<T, N>& seg,
                                      const Tolerance& tol = Tolerance()) {
    LineIntersectionResult<T, N> result;

    const auto& P0 = line.point1().coords;
    const auto dir = line.direction();
    const auto& A = seg.start.coords;
    const auto& B = seg.end.coords;
    const auto segDir = B - A;
    T dirNorm = dir.norm();
    T segNorm = segDir.norm();
    T sepNorm = (P0 - A).norm();
    // Adaptive scale for tolerance model (mirrors line–line/plane)
    T scale = dirNorm + segNorm + sepNorm;
    T angularTol = tol.evaluateEpsilon(std::max(scale, T(1)));
    // Compute normalized directions
    Eigen::Matrix<T, N, 1> d1n = dir;
    if (dirNorm > 0) d1n /= dirNorm;

    Eigen::Matrix<T, N, 1> d2n = segDir;
    if (segNorm > 0) d2n /= segNorm;

    T cosTheta = std::abs(d1n.dot(d2n));
    T sinTheta = std::sqrt(std::max(T(0), T(1) - cosTheta * cosTheta));
    if (sinTheta < angularTol) {
        // Nearly parallel: project segment start onto line and check separation
        Eigen::Matrix<T, N, 1> diff = A - P0;
        Eigen::Matrix<T, N, 1> proj = diff - diff.dot(d1n) * d1n;
        T eps = tol.absTol + tol.evaluateEpsilon(scale);
        if (proj.norm() < eps) {
            // Coincident within tolerance
            result.intersects = true;
            result.lines = { line };
            result.description = "Line and segment are coincident within tolerance";
        } else {
            result.description = "Lines nearly parallel (no intersection)";
        }
        return result;
    }
    // Not parallel: solve for intersection
    Eigen::Matrix<T, N, 2> A_mat;
    for (int i = 0; i < N; ++i) {
        A_mat(i, 0) = dir[i];
        A_mat(i, 1) = -segDir[i];
    }
    Eigen::Matrix<T, N, 1> b = A - P0;
    Eigen::Matrix<T, 2, 1> x = A_mat.colPivHouseholderQr().solve(b);
    T t = x(0, 0);
    T u = x(1, 0);
    // Check if intersection lies within segment
    if (u < -tol.evaluateEpsilon(std::max(scale, T(1))) || u > 1.0 + tol.evaluateEpsilon(std::max(scale, T(1)))) {
        result.description = "Intersection outside segment bounds";
        return result;
    }
    Point<T, N> P_int(P0 + t * dir);
    Point<T, N> S_int(A + u * segDir);
    T dist = (P_int - S_int).norm();
    if (dist <= tol.absTol + tol.evaluateEpsilon(scale)) {
        result.intersects = true;
        result.points = {P_int};
        result.description = "Line intersects segment at a point";
    } else {
        result.description = "Line and segment are skew (no intersection)";
    }
    return result;
}

// ======================================================
// Line–Face Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& line, const Face<T, N>& face,
                                      const Tolerance& tol = Tolerance()) {
    LineIntersectionResult<T, N> result;
    // Compute scale for tolerance
    const auto& n = face.normal;
    const auto& d = line.direction();
    const auto& p0 = line.point1().coords;
    const auto& b = face.base.coords;
    T scale = d.norm() + n.norm() + (p0 - b).norm();
    T angularTol = tol.evaluateEpsilon(std::max(scale, T(1)));
    // Compute angle between line and face normal
    Eigen::Matrix<T, N, 1> d_hat = (d.norm() > 0) ? d.normalized() : d;
    Eigen::Matrix<T, N, 1> n_hat = (n.norm() > 0) ? n.normalized() : n;
    T cosTheta = std::abs(n_hat.dot(d_hat));
    T sinTheta = std::sqrt(std::max(T(0), T(1) - cosTheta * cosTheta));
    // Intersect with face plane using angular tolerance model
    Plane<T, N> facePlane(face.base, face.normal);
    // Use same logic as line-plane: check for near-parallel
    if (sinTheta < angularTol) {
        // Near-parallel: check if line lies in plane (coincident)
        T dist = n.dot(p0 - b);
        T eps = tol.absTol + tol.evaluateEpsilon(scale);
        if (std::abs(dist) < eps) {
            // Coincident: infinite intersection (line in plane)
            // Only mark as intersecting if the segment of the line overlaps the face polygon
            // We can treat as "coincident within tolerance"
            result.intersects = true;
            result.lines = { line };
            result.description = "Line is coincident with face plane (within tolerance)";
            return result;
        } else {
            result.description = "Lines nearly parallel (no intersection with face plane)";
            return result;
        }
    }
    // Not parallel: intersect line with plane
    auto planeRes = intersect(line, facePlane, tol);
    if (!planeRes.intersects || planeRes.points.empty()) {
        result.description = "Line does not intersect face plane";
        return result;
    }
    Point<T, N> candidate = planeRes.points[0];
    // If line is nearly parallel, only accept if offset is within tolerance
    if (sinTheta < angularTol * T(10)) {
        // Extra check: only continue to polygon test if offset is within tol.absTol + tol.evaluateEpsilon(scale)
        T dist = n.dot(p0 - b);
        if (std::abs(dist) > tol.absTol + tol.evaluateEpsilon(scale)) {
            result.description = "Line nearly parallel to face plane and offset exceeds tolerance";
            return result;
        }
    }
    // Check if intersection point lies within the face polygon
    if (intersect(candidate, face, tol).intersects) {
        result.intersects = true;
        result.points = {candidate};
        result.description = "Line intersects face at a point";
    } else {
        result.description = "Line intersects face plane but outside face bounds";
    }
    return result;
}


// Reverse overload: Curve–Line
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Curve<T, N>& c, const Line<T, N>& line,
                                       const Tolerance& tol = Tolerance(), int max_samples = -1) {
    return intersect(line, c, tol, max_samples);
}

// Reverse overload: Surface–Line
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Surface<T, N>& s, const Line<T, N>& line,
                                       const Tolerance& tol = Tolerance(), int max_grid = -1) {
    return intersect(line, s, tol, max_grid);
}

// Reverse overload: Segment–Line
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Segment<T, N>& seg, const Line<T, N>& line,
                                       const Tolerance& tol = Tolerance()) {
    return intersect(line, seg, tol);
}

// Reverse overload: Face–Line
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Face<T, N>& face, const Line<T, N>& line,
                                       const Tolerance& tol = Tolerance()) {
    return intersect(line, face, tol);
}

} // namespace Euclid::Geometry
