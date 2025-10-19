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
struct LineIntersectionResult : public IntersectionResult {
    std::vector<Point<T, N>> points; // intersection points (if any)
    std::vector<Line<T, N>> lines;   // intersection lines (if any)

    LineIntersectionResult() = default;

    LineIntersectionResult(bool intersects, const std::vector<Point<T, N>>& pts = {},
                           const std::vector<Line<T, N>>& lns = {}, std::string desc = "")
        : IntersectionResult{intersects, desc}, points(pts), lines(lns) {}
};

// ======================================================
// Line–Line Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& l1, const Line<T, N>& l2,
                                      const Tolerance& tol = Tolerance()) {
    LineIntersectionResult<T, N> result;

    T scale = l1.direction().norm() + l2.direction().norm();
    if constexpr (N == 2) {
        // 2D intersection using determinant
        T dx1 = l1.direction()[0], dy1 = l1.direction()[1];
        T dx2 = l2.direction()[0], dy2 = l2.direction()[1];
        T det = dx1 * dy2 - dy1 * dx2;

        if (std::abs(det) < tol.evaluateEpsilon(scale)) {
            // Lines are parallel — check if they are coincident
            auto diff = l2.point1().coords - l1.point1().coords;
            T eps = tol.evaluateEpsilon(scale);

            bool coincident = false;
            if (diff.norm() < eps) {
                coincident = true;
            } else {
                T dot = diff.normalized().dot(l1.direction().normalized());
                coincident = std::abs(std::abs(dot) - 1.0) < eps;
            }

            if (coincident) {
                result.intersects = true;
                result.lines = {l1};
                result.description = "Lines are coincident";
            } else {
                result.description = "Lines are parallel";
            }
            return result;
        }
        auto diff = l2.point1().coords - l1.point1().coords;
        T t1 = (diff[0] * dy2 - diff[1] * dx2) / det;
        Point<T, N> pt(l1.point1().coords + t1 * l1.direction());
        result.intersects = true;
        result.points = {pt};
        result.description = "Lines intersect at a point";
        return result;
    }

    if constexpr (N == 3) {
        // 3D intersection using closest point method
        auto p1 = l1.point1().coords;
        auto d1 = l1.direction();
        auto p2 = l2.point1().coords;
        auto d2 = l2.direction();
        auto r = p1 - p2;
        T a = d1.dot(d1), b = d1.dot(d2), c = d2.dot(d2);
        T d = d1.dot(r), e = d2.dot(r);
        T denom = a * c - b * b;

        if (std::abs(denom) < tol.evaluateEpsilon(scale)) {
            // Lines are parallel — check if they are coincident
            auto diff = l2.point1().coords - l1.point1().coords;
            T eps = tol.evaluateEpsilon(scale);

            bool coincident = false;
            if (diff.norm() < eps) {
                coincident = true;
            } else {
                T dot = diff.normalized().dot(l1.direction().normalized());
                coincident = std::abs(std::abs(dot) - 1.0) < eps;
            }

            if (coincident) {
                result.intersects = true;
                result.lines = {l1};
                result.description = "Lines are coincident";
            } else {
                result.description = "Lines are parallel";
            }
            return result;
        }
        T t1 = (b * e - c * d) / denom;
        T t2 = (a * e - b * d) / denom;
        Point<T, N> pt1(p1 + t1 * d1);
        Point<T, N> pt2(p2 + t2 * d2);
        T dist = (pt1 - pt2).norm();
        if (dist > tol.evaluateEpsilon(dist)) {
            result.description = "Lines are skew (do not intersect)";
            return result;
        }
        result.intersects = true;
        result.points = {pt1};
        result.description = "Lines intersect at a point";
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
        T eps = tol.evaluateEpsilon(scale);

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

        if (dist > tol.evaluateEpsilon(dist)) {
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
                                      const Tolerance& tol = Tolerance(), int max_samples = -1) {
    LineIntersectionResult<T, N> result;

    // domain of curve
    auto domain = c.domain();
    T t0 = domain.first.first;
    T t1 = domain.first.second;
    if (t1 <= t0) return result;

    // adaptive sample count based on curve extent and tolerance
    T extent = std::abs(t1 - t0);
    T eps = tol.evaluateEpsilon(extent);
    int samples = max_samples > 0 ? max_samples : std::max(20, static_cast<int>(std::min(200.0, std::max(20.0, (double)(extent / eps) * 5.0))));

    // ensure the three canonical samples exist: t0, mid, t1
    std::vector<T> sampled_ts;
    sampled_ts.reserve(samples);
    sampled_ts.push_back(t0);
    sampled_ts.push_back((t0 + t1) / 2);
    sampled_ts.push_back(t1);

    for (int i = 0; i < samples; ++i) {
        T uu = t0 + (t1 - t0) * static_cast<T>(i) / static_cast<T>(samples - 1);
        // avoid duplicating the three canonical samples
        if (std::abs(uu - t0) < eps || std::abs(uu - t1) < eps || std::abs(uu - (t0 + t1) / 2) < eps) continue;
        sampled_ts.push_back(uu);
    }

    // Closest projection onto line helper
    const auto& P0 = line.point1().coords;
    const auto dir = line.direction();

    T bestDist = std::numeric_limits<T>::infinity();
    T best_u = t0;
    T best_t = T(0);
    Point<T, N> best_pt = c.evaluate(t0);

    for (size_t i = 0; i < sampled_ts.size(); ++i) {
        T u = sampled_ts[i];
        Point<T, N> cpt = c.evaluate(u);
        // projector parameter t along line
        T proj_t = dir.dot(cpt.coords - P0) / dir.squaredNorm();
        Point<T, N> lpt(P0 + proj_t * dir);
        T dist = (cpt - lpt).norm();
        if (dist < bestDist) {
            bestDist = dist;
            best_u = u;
            best_t = proj_t;
            best_pt = cpt;
        }
    }

    // Newton-like refinement solving for (u, t) that minimize ||C(u) - L(t)||^2.
    // Use numerical derivative for C'(u).
    int max_iters = 20;
    T alpha = T(0.2);
    T prevDist = bestDist;
    T u = best_u;
    T t = best_t;

    for (int iter = 0; iter < max_iters; ++iter) {
        Point<T, N> Cu = c.evaluate(u);
        // numerical derivative C'(u)
        T h = std::max((T)1e-7, (t1 - t0) * (T)1e-6);
        T up = std::min(t1, u + h);
        T um = std::max(t0, u - h);
        Eigen::Matrix<T, N, 1> deriv = (c.evaluate(up).coords - c.evaluate(um).coords) / (up - um);

        // residual and jacobian
        Eigen::Matrix<T, N, 1> R;
        for (int i = 0; i < N; ++i) R(i, 0) = Cu.coords[i] - (P0[i] + t * dir[i]);

        // J = [C'(u), -dir]
        Eigen::Matrix<T, N, 2> J;
        for (int i = 0; i < N; ++i) {
            J(i, 0) = deriv[i];
            J(i, 1) = -dir[i];
        }

        Eigen::Matrix<T, 2, 2> JTJ = J.transpose() * J;
        Eigen::Matrix<T, 2, 1> JTR = J.transpose() * R;

        // solve JTJ * delta = -JTR with damping fallback
        Eigen::Matrix<T, 2, 1> delta;
        Eigen::FullPivLU<Eigen::Matrix<T, 2, 2>> lu(JTJ);
        if (lu.isInvertible()) {
            delta = lu.solve(-JTR);
        } else {
            // damped least squares fallback
            Eigen::Matrix<T, 2, 2> JTJ_damped = JTJ;
            JTJ_damped(0, 0) += eps;
            JTJ_damped(1, 1) += eps;
            delta = JTJ_damped.colPivHouseholderQr().solve(-JTR);
        }

        T du = delta(0, 0);
        T dt = delta(1, 0);

        u += alpha * du;
        t += alpha * dt;

        // clamp
        if (u < t0) u = t0;
        if (u > t1) u = t1;

        // recompute distance
        Point<T, N> Cnew = c.evaluate(u);
        Point<T, N> Lnew(P0 + t * dir);
        T dist = (Cnew - Lnew).norm();

        if (std::abs(dist - prevDist) < eps * 0.5) break;
        prevDist = dist;
    }

    // final check
    Point<T, N> Cfinal = c.evaluate(u);
    Point<T, N> Lfinal(P0 + t * dir);
    T finalDist = (Cfinal - Lfinal).norm();
    if (finalDist <= tol.evaluateEpsilon(std::max((T)1.0, finalDist))) {
        result.intersects = true;
        result.points = {Point<T, N>( (Cfinal.coords + Lfinal.coords) * (T)0.5 )};
        result.description = "Line intersects curve (within tolerance)";
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

    // Get some scale for tolerance
    const auto& P0 = line.point1().coords;
    const auto dir = line.direction();
    T dir_norm = dir.norm();
    if (dir_norm == 0) {
        result.description = "Line direction is zero";
        return result;
    }
    const auto dir_hat = dir / dir_norm;
    // Hint the surface BVH/normal orientation with the ray direction for stable signed distances
    s.setBVHDirectionHint(dir_hat);

    // Estimate extent using surface bounding box if available, else use domain extents
    auto domain = s.surfaceDomain;
    T u0 = domain.first.first, u1 = domain.first.second;
    T v0 = domain.second.first, v1 = domain.second.second;
    T u_extent = std::abs(u1 - u0);
    T v_extent = std::abs(v1 - v0);
    T surface_extent = std::max(u_extent, v_extent);
    // Try to get a bounding box (if implemented for Surface)
    T t_min = -10 * surface_extent, t_max = 10 * surface_extent;

    // Hint the surface BVH/normal orientation with the ray direction for stable signed distances
    s.setBVHDirectionHint(dir_hat);

    // Refine BVH w.r.t. BOTH directions of the infinite line so patches behind P0 are not culled.
    // This ensures torus-style 4-hit configurations (±(R±r)) are not lost by directional pruning.
    s.refineBVHNearRay(P0, dir_hat, tol);
    s.refineBVHNearRay(P0, -dir_hat, tol);

    // Fetch BVH patches after symmetric refinement
    auto patches = s.getBVH(max_grid > 0 ? max_grid : 40, tol, surface_extent, 9);

    // Build ray intervals from BVH patch AABBs
    std::vector<std::pair<T, T>> t_intervals;
    t_intervals.reserve(patches.size());
    for (const auto& p : patches) {
        bool hit; T t_enter, t_exit;
        std::tie(hit, t_enter, t_exit) = rayAABBInterval<T, N>(P0, dir, p.minPoint.coords, p.maxPoint.coords, tol, surface_extent);
        if (hit && t_exit >= t_min && t_enter <= t_max) {
            t_intervals.emplace_back(
                std::max(t_enter - tol.evaluateEpsilon(surface_extent), t_min),
                std::min(t_exit + tol.evaluateEpsilon(surface_extent), t_max)
            );
        }
    }
    // Sort and merge overlapping intervals
    std::sort(t_intervals.begin(), t_intervals.end(), [](const auto& a, const auto& b){ return a.first < b.first; });
    std::vector<std::pair<T, T>> merged_intervals;
    for (const auto& iv : t_intervals) {
        if (merged_intervals.empty() || iv.first > merged_intervals.back().second) {
            merged_intervals.push_back(iv);
        } else {
            merged_intervals.back().second = std::max(merged_intervals.back().second, iv.second);
        }
    }
    
    /*
    // Debug: print BVH t-intervals for the ray
    if constexpr (N == 3) { // only for 3D for clarity
        std::cout << "DEBUG: BVH t-intervals for line-surface:" << '\n';
        if (merged_intervals.empty()) {
            std::cout << "DEBUG: No intervals found." << '\n';
        } else {
            for (const auto& iv : merged_intervals) {
                std::cout << "DEBUG: [" << iv.first << ", " << iv.second << "]" << '\n';
            }
        }
    }
     */

    // ======================================================
    // (u,v,t) Newton solver based on BVH patch midpoints
    // ======================================================
    auto dedupHits = [&](std::vector<Point<T,N>>& hits, const Point<T,N>& newPt, T tol) {
        for (const auto& h : hits) {
            if ((h.coords - newPt.coords).norm() < tol) return false;
        }
        hits.emplace_back(newPt);
        return true;
    };

    T tol_val = tol.evaluateEpsilon(surface_extent);

    auto polishHit = [&](Point<T,N>& Spt, T& u, T& v, T& t) {
        for (int i = 0; i < 2; ++i) {
            Point<T,N> Lpt(P0 + t * dir);
            Eigen::Matrix<T,N,1> R = Spt.coords - Lpt.coords;
            if (R.norm() < tol_val) break;
            T uh = std::max(tol_val, std::numeric_limits<T>::epsilon() * 10);
            T vh = std::max(tol_val, std::numeric_limits<T>::epsilon() * 10);
            Point<T,N> Sup = s.evaluate(std::min(u1,u+uh), v);
            Point<T,N> Sum = s.evaluate(std::max(u0,u-uh), v);
            Point<T,N> Svp = s.evaluate(u, std::min(v1,v+vh));
            Point<T,N> Svm = s.evaluate(u, std::max(v0,v-vh));
            Eigen::Matrix<T,N,1> Su = (Sup.coords - Sum.coords) / (2*uh);
            Eigen::Matrix<T,N,1> Sv = (Svp.coords - Svm.coords) / (2*vh);
            Eigen::Matrix<T,N,3> J;
            for (int k = 0; k < N; ++k) { J(k,0)=Su[k]; J(k,1)=Sv[k]; J(k,2)=-dir[k]; }
            Eigen::Matrix<T,N,1> negR = -R;
            Eigen::Matrix<T,3,1> delta = J.colPivHouseholderQr().solve(negR);
            u += delta(0,0); v += delta(1,0); t += delta(2,0);
            Spt = s.evaluate(u,v);
        }
    };

    std::vector<Point<T,N>> intersections;
    T dedup_tol = tol_val * 1.5;
    int max_iters = 20;

    for (const auto& iv : merged_intervals) {
        // Adaptive seeding: number of seeds based on interval length and surface extent
        int num_seeds = std::max(5, (int)((iv.second - iv.first) / (surface_extent * 0.1)));
        for (int i = 0; i < num_seeds; ++i) {
            T t_seed = iv.first + (iv.second - iv.first) * (i + 0.5) / num_seeds;
            // Add random offset to seed point
            Eigen::Matrix<T,N,1> offset = Eigen::Matrix<T,N,1>::Random().normalized() * (tol_val * 10);
            Point<T,N> P_seed(P0 + t_seed * dir + offset);
            // Point<T,N> P_seed(P0 + t_seed * dir);
            auto proj = s.projectPoint(P_seed.coords, tol);
            T u = proj.u;
            T v = proj.v;
            T t = t_seed;

            for (int iter = 0; iter < max_iters; ++iter) {
                Point<T,N> Spt = s.evaluate(u, v);
                Eigen::Matrix<T,N,1> R = Spt.coords - (P0 + t * dir);
                if (R.norm() < tol_val) break;

                // Partial derivatives
                T uh = std::max(tol_val, std::numeric_limits<T>::epsilon() * 10);
                T vh = std::max(tol_val, std::numeric_limits<T>::epsilon() * 10);
                Point<T,N> Sup = s.evaluate(std::min(u1,u+uh), v);
                Point<T,N> Sum = s.evaluate(std::max(u0,u-uh), v);
                Point<T,N> Svp = s.evaluate(u, std::min(v1,v+vh));
                Point<T,N> Svm = s.evaluate(u, std::max(v0,v-vh));
                Eigen::Matrix<T,N,1> Su = (Sup.coords - Sum.coords) / (2*uh);
                Eigen::Matrix<T,N,1> Sv = (Svp.coords - Svm.coords) / (2*vh);

                Eigen::Matrix<T,N,3> J;
                for (int k = 0; k < N; ++k) {
                    J(k,0) = Su[k];
                    J(k,1) = Sv[k];
                    J(k,2) = -dir[k];
                }

                Eigen::Matrix<T,3,1> delta = J.colPivHouseholderQr().solve(-R);
                T du = delta(0,0), dv = delta(1,0), dtau = delta(2,0);

                // Damping for stability
                T alpha = 0.8;
                u += alpha * du;
                v += alpha * dv;
                t += alpha * dtau;

                // Clamp to domain
                u = std::clamp(u, u0, u1);
                v = std::clamp(v, v0, v1);
            }

            // Check final residual
            Point<T,N> Sfinal = s.evaluate(u,v);
            Point<T,N> Lfinal(P0 + t * dir);
            T dist = (Sfinal - Lfinal).norm();
            if (dist <= tol_val * 5.0) {
                polishHit(Sfinal, u, v, t);
                Point<T,N> polished((Sfinal.coords + Lfinal.coords) * 0.5);
                dedupHits(intersections, polished, dedup_tol);
            }
        }
    }

    // After all BVH operations are done, clear the BVH direction hint
    s.clearBVHDirectionHint();
    if (!intersections.empty()) {
        result.intersects = true;
        result.points = intersections;
        result.description = "Line intersects surface (u,v,t Newton solver)";
    } else {
        result.description = "Line does not intersect surface (u,v,t Newton)";
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
    T eps = tol.evaluateEpsilon(scale);

    if constexpr (N == 3) {
        // 3D fast path: geometric intersection
        T denom = n.dot(d);
        if (std::abs(denom) < eps) {
            T dist = n.dot(p0 - b);
            if (std::abs(dist) < eps) {
                result.intersects = true;
                result.lines = {line};
                result.description = "Line lies in plane (coincident)";
            } else {
                result.description = "Line is parallel to plane";
            }
            return result;
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
        if (std::abs(denom) < eps) {
            T dist = n.dot(p0 - b);
            if (std::abs(dist) < eps) {
                result.intersects = true;
                result.lines = {line};
                result.description = "Line lies in plane (coincident)";
            } else {
                result.description = "Line is parallel to plane";
            }
            return result;
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

    // Solve for t (line) and u (segment) such that P0 + t*dir = A + u*segDir
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
    if (u < -tol.evaluateEpsilon(1.0) || u > 1.0 + tol.evaluateEpsilon(1.0)) {
        result.description = "Line and segment do not intersect (outside segment bounds)";
        return result;
    }

    Point<T, N> P_int(P0 + t * dir);
    Point<T, N> S_int(A + u * segDir);

    if ((P_int - S_int).norm() <= tol.evaluateEpsilon((P_int - S_int).norm())) {
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

    // First, intersect with the face’s supporting plane
    Plane<T, N> facePlane(face.base, face.normal);
    auto planeRes = intersect(line, facePlane, tol);
    if (!planeRes.intersects || planeRes.points.empty()) {
        result.description = "Line does not intersect face plane";
        return result;
    }

    Point<T, N> candidate = planeRes.points[0];

    // Check if intersection point lies within the face polygon
    if (intersect(candidate, face, tol).intersects) {
        result.intersects = true;
        result.points = {candidate};
        result.description = "Line intersects face at a point";
    } else {
        result.description = "Line intersects face plane but outside bounds";
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
