#pragma once

#include "../geometry.hpp"
#include "../tolerance.hpp"
#include "intersection_result.hpp"

#include <vector>

namespace Euclid::Geometry {

// Specialized intersection result for point intersections
template <typename T, int N>
struct PointIntersectionResult : public IntersectionResult {
    std::vector<Point<T, N>> points;
    PointIntersectionResult(bool inter, const std::string& desc)
        : IntersectionResult(inter, desc) {}
    PointIntersectionResult(bool inter, const std::string& desc, const std::vector<Point<T, N>>& pts)
        : IntersectionResult(inter, desc), points(pts) {}
};

// ======================================================
// Point–Point Intersection
// ======================================================
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Point<T, N>& a, const Point<T, N>& b,
                             const Tolerance& tol = Tolerance()) {
    bool equal = a.isEqual(b, tol);
    if (equal) {
        return PointIntersectionResult<T, N>(true, "Points coincide", {a});
    } else {
        return PointIntersectionResult<T, N>(false, "Points are distinct");
    }
}

// ======================================================
// Point–Line Intersection (projection / closest point)
// ======================================================
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Point<T, N>& p, const Line<T, N>& l,
                             const Tolerance& tol = Tolerance()) {
    using Vec = Eigen::Matrix<T, N, 1>;
    Vec lineDir = l.direction();
    Vec diff = p.coords - l.point1().coords;

    // Project point onto line
    T t = lineDir.dot(diff);
    Vec projected = l.point1().coords + t * lineDir;
    Point<T, N> closest(projected);

    bool onLine = Euclid::equalWithinTolerance((closest - p).norm(), 0.0, tol, p.coordinateMagnitude());
    if (onLine) {
        return PointIntersectionResult<T, N>(true, "Point lies on line", {p});
    } else {
        return PointIntersectionResult<T, N>(false, "Closest point on line", {closest});
    }
}

// Added reverse overload for Line-Point intersection
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Line<T, N>& l, const Point<T, N>& p,
                             const Tolerance& tol = Tolerance()) {
    return intersect(p, l, tol);
}

// ======================================================
// Point–Segment Intersection (checks if point lies on segment)
// ======================================================
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Point<T, N>& p, const Segment<T, N>& seg,
                             const Tolerance& tol = Tolerance()) {
    using Vec = Eigen::Matrix<T, N, 1>;
    const Vec& a = seg.start.coords;
    const Vec& b = seg.end.coords;
    const Vec& x = p.coords;
    Vec ab = b - a;
    Vec ap = x - a;

    T abNorm = ab.norm();

    // Check for degenerate segment (start == end)
    if (abNorm == 0) {
        if (p.isEqual(seg.start, tol)) {
            return PointIntersectionResult<T, N>(true, "Point coincides with degenerate segment point", {p});
        } else {
            return PointIntersectionResult<T, N>(false, "Point does not coincide with degenerate segment point", {seg.start});
        }
    }

    T apNorm = ap.norm();
    T bpNorm = (x - b).norm();

    bool colinear = false;
    if constexpr (N == 3) {
        // Use cross product norm for 3D
        T crossNorm = (ab.cross(ap)).norm();
        colinear = Euclid::equalWithinTolerance(crossNorm, 0.0, tol, abNorm + apNorm);
    } else {
        // Use projection residual for ND (including 2D)
        Vec proj = ab * (ab.dot(ap) / ab.squaredNorm());
        Vec residual = ap - proj;
        T residualNorm = residual.norm();
        colinear = residualNorm <= tol.evaluateEpsilon(N) * (abNorm + apNorm);
    }

    // Check if projection is within segment bounds
    T dot = ab.dot(ap);
    bool within = dot >= -tol.evaluateEpsilon(N) * abNorm * apNorm && dot <= ab.squaredNorm() + tol.evaluateEpsilon(N) * abNorm * apNorm;

    // Also check if point is close to segment endpoints (for degenerate/short segments)
    bool closeToA = Euclid::equalWithinTolerance(apNorm, 0.0, tol, abNorm);
    bool closeToB = Euclid::equalWithinTolerance(bpNorm, 0.0, tol, abNorm);

    if ((colinear && within) || closeToA || closeToB) {
        return PointIntersectionResult<T, N>(true, "Point lies on segment", {p});
    } else {
        // Find closest point on segment to p
        T t = ab.dot(ap) / ab.squaredNorm();
        t = std::max(T(0), std::min(T(1), t));
        Vec closestCoords = a + t * ab;
        Point<T, N> closest(closestCoords);
        return PointIntersectionResult<T, N>(false, "Closest point on segment", {closest});
    }
}

// Added reverse overload for Segment-Point intersection
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Segment<T, N>& seg, const Point<T, N>& p,
                             const Tolerance& tol = Tolerance()) {
    return intersect(p, seg, tol);
}


// ======================================================
// Point–Plane Intersection (projection)
// ======================================================
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Point<T, N>& p, const Plane<T, N>& plane,
                             const Tolerance& tol = Tolerance()) {
    using Vec = Eigen::Matrix<T, N, 1>;
    Vec n = plane.normal;
    T dist = n.dot(p.coords - plane.base.coords);
    Vec projected = p.coords - dist * n;
    Point<T, N> projPt(projected);

    bool onPlane = Euclid::equalWithinTolerance(std::abs(dist), 0.0, tol, p.coordinateMagnitude());
    if (onPlane) {
        return PointIntersectionResult<T, N>(true, "Point lies on plane", {p});
    } else {
        return PointIntersectionResult<T, N>(false, "Projected point on plane", {projPt});
    }
}

// Added reverse overload for Plane-Point intersection
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Plane<T, N>& plane, const Point<T, N>& p,
                             const Tolerance& tol = Tolerance()) {
    return intersect(p, plane, tol);
}

// ======================================================
// Point–Curve Intersection (coarse sampling + local refinement)
// ======================================================
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Point<T, N>& p, const Curve<T, N>& c,
                             const Tolerance& tol = Tolerance(), int N_samples = -1) {
    using Vec = Eigen::Matrix<T, N, 1>;
    auto [t0, t1] = c.domain();
    T extent = std::abs(t1 - t0);

    // --- Adaptive coarse sampling ---
    // Try to estimate a "length" or bounding box extent as basis for adaptive sampling
    T curve_length_est = extent;
    {
        // Optionally, sample the curve at a few points to estimate length
        constexpr int quick_samples = 7;
        T length_sum = 0;
        Point<T, N> prev = c.evaluate(t0);
        for (int i = 1; i < quick_samples; ++i) {
            T t = t0 + (t1 - t0) * T(i) / T(quick_samples - 1);
            Point<T, N> curr = c.evaluate(t);
            length_sum += (curr.coords - prev.coords).norm();
            prev = curr;
        }
        curve_length_est = std::max(length_sum, extent);
    }
    T tol_val = tol.evaluateEpsilon(N);
    int samples = N_samples;
    if (samples <= 0) {
        samples = std::max(5, std::min(50, static_cast<int>(curve_length_est / tol_val)));
    }
    if (samples < 20) samples = 20;

    // --- Coarse sampling: always sample t0, (t0+t1)/2, t1 first ---
    T minDist = std::numeric_limits<T>::max();
    T best_t = t0;
    Point<T, N> best_pt = c.evaluate(t0);
    struct SampleInfo {
        int index;
        T t;
        T dist;
        Point<T, N> pt;
    };
    std::vector<SampleInfo> coarse_samples;
    // Three mandatory coarse samples
    std::vector<std::pair<int, T>> fixed_samples = {
        {0, t0},
        {1, (t0 + t1) / T(2)},
        {2, t1}
    };
    for (size_t i = 0; i < fixed_samples.size(); ++i) {
        int idx = fixed_samples[i].first;
        T t = fixed_samples[i].second;
        Point<T, N> cpt = c.evaluate(t);
        T dist = (cpt.coords - p.coords).norm();
#ifdef DEBUG
        std::cout << "[DEBUG] Coarse sample " << idx << ": t=" << t << ", dist=" << dist << std::endl;
#endif
        coarse_samples.push_back({idx, t, dist, cpt});
        if (dist < minDist) {
            minDist = dist;
            best_t = t;
            best_pt = cpt;
        }
    }
    // Now continue with adaptive coarse sampling for other samples, skipping t0, (t0+t1)/2, t1
    for (int i = 0; i < samples; ++i) {
        T t = t0 + (t1 - t0) * T(i) / T(samples - 1);
        // skip if t matches one of the three already sampled points (avoid double sample)
        bool is_fixed = false;
        for (const auto& s : fixed_samples) {
            if (std::abs(t - s.second) <= std::numeric_limits<T>::epsilon() * 10 || (samples <= 3 && t == s.second)) {
                is_fixed = true;
                break;
            }
        }
        if (is_fixed) continue;
        Point<T, N> cpt = c.evaluate(t);
        T dist = (cpt.coords - p.coords).norm();
#ifdef DEBUG
        std::cout << "[DEBUG] Coarse sample " << i << ": t=" << t << ", dist=" << dist << std::endl;
#endif
        if (dist < minDist) {
            minDist = dist;
            best_t = t;
            best_pt = cpt;
        }
        // Early exit if within tolerance
        if (dist < tol_val) {
            break;
        }
    }
    // --- Newton-Raphson refinement ---
    // F(t) = ||c(t) - p||^2, dF/dt = 2 * (c(t) - p).dot(c'(t))
    T eps = tol_val;
    int max_iters = 15;
    T alpha = 0.3;
    T t = best_t;
    T dist_old = (c.evaluate(t).coords - p.coords).norm();
    for (int iter = 0; iter < max_iters; ++iter) {
        Vec c_t = c.evaluate(t).coords;
        Vec diff = c_t - p.coords;
        // Always use numerical derivative
        T h = std::max(T(1e-7), (t1 - t0) * T(1e-6));
        T t_ph = std::min(t1, t + h);
        T t_mh = std::max(t0, t - h);
        Vec deriv = (c.evaluate(t_ph).coords - c.evaluate(t_mh).coords) / (t_ph - t_mh);
        T grad = 2 * diff.dot(deriv);
        // Take step
        t = t - alpha * grad;
        // Clamp to domain
        t = std::max(t0, std::min(t1, t));
        // Check convergence
        T dist_new = (c.evaluate(t).coords - p.coords).norm();
        if (std::abs(dist_new - dist_old) < eps) {
            break;
        }
        dist_old = dist_new;
    }
    // Set result to refined value
    best_t = t;
    best_pt = c.evaluate(best_t);
    minDist = (best_pt.coords - p.coords).norm();

    bool coincides = Euclid::equalWithinTolerance(minDist, 0.0, tol, p.coordinateMagnitude());
    if (coincides) {
        return PointIntersectionResult<T, N>(true, "Point lies on curve", {p});
    } else {
        return PointIntersectionResult<T, N>(false, "Closest point on curve", {best_pt});
    }
}

// Added reverse overload for Curve-Point intersection
template <typename T, int N>
PointIntersectionResult<T, N> intersect(
    const Curve<T, N>& c, const Point<T, N>& p,
    const Tolerance& tol = Tolerance(), int N_samples = -1)
{
    return intersect(p, c, tol, N_samples);
}

// ======================================================
// Point–Surface Intersection (projection / closest point search)
// ======================================================
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Point<T, N>& p, const Surface<T, N>& s,
                             const Tolerance& tol = Tolerance()) {
    using Vec = Eigen::Matrix<T, N, 1>;
    auto domainPair = s.surfaceDomain;
    auto [u0, u1] = domainPair.first;
    auto [v0, v1] = domainPair.second;
    T u_extent = std::abs(u1 - u0);
    T v_extent = std::abs(v1 - v0);

    // --- Adaptive coarse grid ---
    T eps = tol.evaluateEpsilon(N);
    int M = static_cast<int>((u_extent + v_extent) / eps * 20);
    if (M < 50) M = 50;
    if (M > 500) M = 500;

    T minDist = std::numeric_limits<T>::max();
    T best_u = u0;
    T best_v = v0;
    Point<T, N> best_pt;

    for (int i = 0; i < M; ++i) {
        T u = u0 + (u1 - u0) * T(i) / T(M - 1);
        for (int j = 0; j < M; ++j) {
            T v = v0 + (v1 - v0) * T(j) / T(M - 1);
            Point<T, N> spt = s.evaluate(u, v);
            T dist = (spt.coords - p.coords).norm();
            if (dist < minDist) {
                minDist = dist;
                best_u = u;
                best_v = v;
                best_pt = spt;
            }
        }
    }

    // --- Local refinement ---
    int refinement_samples = 9;
    T refinement_u_range = std::max(u_extent * 0.15, eps);
    T refinement_v_range = std::max(v_extent * 0.15, eps);
    T half_u_range = refinement_u_range / 2;
    T half_v_range = refinement_v_range / 2;

    T refine_u_start = std::max(u0, best_u - half_u_range);
    T refine_u_end = std::min(u1, best_u + half_u_range);
    T refine_v_start = std::max(v0, best_v - half_v_range);
    T refine_v_end = std::min(v1, best_v + half_v_range);

    for (int i = 0; i < refinement_samples; ++i) {
        T u = refine_u_start + (refine_u_end - refine_u_start) * T(i) / T(refinement_samples - 1);
        for (int j = 0; j < refinement_samples; ++j) {
            T v = refine_v_start + (refine_v_end - refine_v_start) * T(j) / T(refinement_samples - 1);
            Point<T, N> spt = s.evaluate(u, v);
            T dist = (spt.coords - p.coords).norm();
            if (dist < minDist) {
                minDist = dist;
                best_u = u;
                best_v = v;
                best_pt = spt;
            }
        }
    }

    // --- Newton-Raphson iterative solver for closest point projection ---
    int max_iters = 30;
    T alpha = 0.2;

    auto S = [&](T u, T v) {
        return s.evaluate(u, v).coords;
    };

    auto partialDerivativeU = [&](T u, T v) {
        T h = eps;
        T u_minus = std::max(u0, u - h);
        T u_plus = std::min(u1, u + h);
        Vec diff = (S(u_plus, v) - S(u_minus, v)) / (u_plus - u_minus);
        return diff;
    };

    auto partialDerivativeV = [&](T u, T v) {
        T h = eps;
        T v_minus = std::max(v0, v - h);
        T v_plus = std::min(v1, v + h);
        Vec diff = (S(u, v_plus) - S(u, v_minus)) / (v_plus - v_minus);
        return diff;
    };

    T u = best_u;
    T v = best_v;
    T prevDist = minDist; // initialize with distance from coarse grid

    for (int iter = 0; iter < max_iters; ++iter) {
        Vec Suv = S(u, v);
        Vec dSdu = partialDerivativeU(u, v);
        Vec dSdv = partialDerivativeV(u, v);
        Vec diff = Suv - p.coords;

        T grad_u = 2 * diff.dot(dSdu);
        T grad_v = 2 * diff.dot(dSdv);

        u = u - alpha * grad_u;
        v = v - alpha * grad_v;

        u = std::min(u1, std::max(u0, u));
        v = std::min(v1, std::max(v0, v));

        T dist = (S(u, v) - p.coords).norm();
        if (std::abs(dist - prevDist) < eps) {
            break;
        }
        prevDist = dist;
    }

    best_u = u;
    best_v = v;
    best_pt = s.evaluate(best_u, best_v);
    minDist = (best_pt.coords - p.coords).norm();

    // Final intersection test
    T thres = std::max(eps, 0.01 * std::max(u_extent, v_extent)); // scale tolerance relative to surface size
    bool intersects = minDist <= thres;
    if (intersects) {
        return PointIntersectionResult<T, N>(true, "Point lies on surface", {p});
    } else {
        return PointIntersectionResult<T, N>(false, "Closest point on surface", {best_pt});
    }
}

// ======================================================
// Point–Face Intersection (projection and containment test)
// ======================================================
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Point<T, N>& p, const Face<T, N>& f,
                             const Tolerance& tol = Tolerance()) {
    using Vec = Eigen::Matrix<T, N, 1>;
    // Project point onto face plane
    Vec n = f.normal;
    T dist = n.dot(p.coords - f.base.coords);
    Vec projected = p.coords - dist * n;
    Point<T, N> projPt(projected);

    bool onPlane = std::abs(dist) <= tol.evaluateEpsilon(N);
    bool inside = true;
    T eps = tol.evaluateEpsilon(N);
    for (size_t i = 0; i < f.vertices.size(); ++i) {
        const Vec& v_i = f.vertices[i].coords;
        const Vec& v_next = f.vertices[(i + 1) % f.vertices.size()].coords;
        Vec edge = v_next - v_i;
        Vec toProj = projected - v_i;
        T crossMag;
        if constexpr (N == 2) {
            crossMag = edge[0] * toProj[1] - edge[1] * toProj[0];
        } else if constexpr (N == 3) {
            crossMag = edge.cross(toProj).dot(n);
        } else {
            Eigen::Matrix<T, N, 1> proj = toProj - (edge.dot(toProj) / edge.squaredNorm()) * edge;
            crossMag = n.dot(proj);
        }
        if (crossMag < -eps) {
            inside = false;
            break;
        }
    }
    bool intersects = onPlane && inside;
    return PointIntersectionResult<T, N>(
        intersects,
        intersects ? "Point lies on face" : "Point does not intersect face",
        {projPt}
    );
}

// Added reverse overload for Face–Point intersection
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Face<T, N>& f, const Point<T, N>& p,
                             const Tolerance& tol = Tolerance()) {
    return intersect(p, f, tol);
}

} // namespace Euclid::Geometry

// Added reverse overload for Surface-Point intersection
template <typename T, int N>
Euclid::Geometry::PointIntersectionResult<T, N> intersect(
    const Euclid::Geometry::Surface<T, N>& s,
    const Euclid::Geometry::Point<T, N>& p,
    const Euclid::Tolerance& tol = Euclid::Tolerance())
{
    return intersect(p, s, tol);
}
