#pragma once

#include "../geometry.hpp"
#include "intersection_result.hpp"
#include "geometry/tolerance.hpp"

#include <vector>
#include <Eigen/Core>

namespace Euclid::Geometry {

// Helper function for 2D cross product (scalar)
template <typename T>
T cross2d(const Eigen::Matrix<T, 2, 1>& a, const Eigen::Matrix<T, 2, 1>& b) {
    return a.x() * b.y() - a.y() * b.x();
}

// Segment-Segment intersection
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Segment<T, N>& seg1, const Segment<T, N>& seg2,
                             const Tolerance& tol = Tolerance()) {
    using Vector = Eigen::Matrix<T, N, 1>;

    Vector p = seg1.start.coords;
    Vector r = seg1.end.coords - seg1.start.coords;
    Vector q = seg2.start.coords;
    Vector s = seg2.end.coords - seg2.start.coords;

    T rxs;
    T q_pxr;
    if constexpr (N == 2) {
        rxs = cross2d(r.eval(), s.eval());
        q_pxr = cross2d((q - p).eval(), r.eval());
    } else if constexpr (N == 3) {
        rxs = r.cross(s).squaredNorm();
        q_pxr = (q - p).cross(r).norm();
    } else {
        // General N-dimensional approach using projection onto the plane spanned by r and s
        T r_norm = r.norm();
        T s_norm = s.norm();
        T r_dot_r = r.dot(r);
        T s_dot_s = s.dot(s);
        T r_dot_s = r.dot(s);
        Vector qp = q - p;
        T qp_dot_r = qp.dot(r);
        T qp_dot_s = qp.dot(s);

        // Compute denominator for parameters t and u
        T denom = r_dot_r * s_dot_s - r_dot_s * r_dot_s;

        if (std::abs(denom) < tol.evaluateEpsilon(std::max(r_norm, s_norm))) {
            // Lines are parallel or degenerate
            T dist = (qp - r * (qp_dot_r / r_dot_r)).norm();
            if (dist < tol.evaluateEpsilon(std::max(r_norm, s_norm))) {
                // Colinear case
                T t0 = qp_dot_r / r_dot_r;
                T t1 = t0 + r_dot_s / r_dot_r;

                T tmin = std::min(t0, t1);
                T tmax = std::max(t0, t1);

                if (tmin > 1 || tmax < 0) {
                    // No overlap
                    return PointIntersectionResult<T, N>(false, "No intersection", {});
                }

                // Overlapping segment
                T overlap_start = std::max(T(0), tmin);
                T overlap_end = std::min(T(1), tmax);

                Vector overlap_seg_start = p + r * overlap_start;
                Vector overlap_seg_end = p + r * overlap_end;

                std::vector<Point<T, N>> points;
                if ((overlap_seg_start - overlap_seg_end).norm() < tol.evaluateEpsilon((overlap_seg_start - overlap_seg_end).norm())) {
                    points.push_back(Point<T, N>{overlap_seg_start});
                } else {
                    points.push_back(Point<T, N>{overlap_seg_start});
                    points.push_back(Point<T, N>{overlap_seg_end});
                }
                return PointIntersectionResult<T, N>(true, "Segments intersect", std::move(points));
            } else {
                // Parallel and non-intersecting
                return PointIntersectionResult<T, N>(false, "No intersection", {});
            }
        } else {
            // Lines are not parallel, compute parameters t and u
            T t = (qp_dot_r * s_dot_s - qp_dot_s * r_dot_s) / denom;
            T u = (qp_dot_r * r_dot_s - qp_dot_s * r_dot_r) / denom;

            if (t >= -tol.evaluateEpsilon(r_norm) && t <= 1 + tol.evaluateEpsilon(r_norm) &&
                u >= -tol.evaluateEpsilon(s_norm) && u <= 1 + tol.evaluateEpsilon(s_norm)) {
                Vector intersection_point = p + r * t;
                std::vector<Point<T, N>> points{Point<T, N>{intersection_point}};
                return PointIntersectionResult<T, N>(true, "Segments intersect", std::move(points));
            } else {
                return PointIntersectionResult<T, N>(false, "No intersection", {});
            }
        }
    }

    if (std::abs(rxs) < tol.evaluateEpsilon(r.norm())) {
        if (std::abs(q_pxr) < tol.evaluateEpsilon(r.norm())) {
            // Segments are colinear
            T t0 = (q - p).dot(r) / r.dot(r);
            T t1 = t0 + s.dot(r) / r.dot(r);

            T tmin = std::min(t0, t1);
            T tmax = std::max(t0, t1);

            if (tmin > 1 || tmax < 0) {
                // No overlap
                return PointIntersectionResult<T, N>(false, "No intersection", {});
            }

            // Overlapping segment
            T overlap_start = std::max(T(0), tmin);
            T overlap_end = std::min(T(1), tmax);

            Vector overlap_seg_start = p + r * overlap_start;
            Vector overlap_seg_end = p + r * overlap_end;

            std::vector<Point<T, N>> points;
            if ((overlap_seg_start - overlap_seg_end).norm() < tol.evaluateEpsilon((overlap_seg_start - overlap_seg_end).norm())) {
                points.push_back(Point<T, N>{overlap_seg_start});
            } else {
                points.push_back(Point<T, N>{overlap_seg_start});
                points.push_back(Point<T, N>{overlap_seg_end});
            }
            return PointIntersectionResult<T, N>(true, "Segments intersect", std::move(points));
        } else {
            // Parallel and non-intersecting
            return PointIntersectionResult<T, N>(false, "No intersection", {});
        }
    } else {
        T t;
        T u;
        if constexpr (N == 2) {
            t = cross2d((q - p).eval(), s.eval()) / rxs;
            u = cross2d((q - p).eval(), r.eval()) / rxs;
        } else if constexpr (N == 3) {
            Vector r_cross_s = r.cross(s);
            T denom = r_cross_s.squaredNorm();
            t = (q - p).cross(s).dot(r_cross_s) / denom;
            u = (q - p).cross(r).dot(r_cross_s) / denom;
        }

        if (t >= -tol.evaluateEpsilon(r.norm()) && t <= 1 + tol.evaluateEpsilon(r.norm()) && u >= -tol.evaluateEpsilon(r.norm()) && u <= 1 + tol.evaluateEpsilon(r.norm())) {
            Vector intersection_point = p + r * t;
            std::vector<Point<T, N>> points{Point<T, N>{intersection_point}};
            return PointIntersectionResult<T, N>(true, "Segments intersect", std::move(points));
        } else {
            return PointIntersectionResult<T, N>(false, "No intersection", {});
        }
    }
}

// Segment-Curve intersection using sampling and Newton refinement
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Segment<T, N>& seg, const Curve<T, N>& curve,
                             const Tolerance& tol = Tolerance()) {
    using Vector = Eigen::Matrix<T, N, 1>;
    constexpr int NUM_SAMPLES = 200; // increase sampling for robustness in higher dimensions
    std::vector<Point<T, N>> intersection_points;

    // Parametric representation of the segment: S(t) = seg.start.coords + t * (seg.end.coords - seg.start.coords), t in [0,1]
    Vector p0 = seg.start.coords;
    Vector p1 = seg.end.coords;
    Vector dir = p1 - p0;
    T seg_length = dir.norm();
    if (seg_length < tol.evaluateEpsilon(seg_length)) {
        return PointIntersectionResult<T, N>(false, "Degenerate segment", {});
    }

    // Curve domain
    auto curve_domain = curve.domain();
    T t_min = curve_domain.first;
    T t_max = curve_domain.second;
    if (t_max <= t_min) t_max = t_min + T(1);

    // Sampling step
    T curve_param_step = (t_max - t_min) / static_cast<T>(NUM_SAMPLES);

    auto is_duplicate = [&](const Vector &pt)->bool{
        for (const auto &existing : intersection_points) {
            if ((existing.coords - pt).norm() < tol.evaluateEpsilon(seg_length)) return true;
        }
        return false;
    };

    for (int i = 0; i <= NUM_SAMPLES; ++i) {
        T t_curve = t_min + static_cast<T>(i) * curve_param_step;
        Vector curve_pt = curve.evaluate(t_curve).coords;

        // Project curve_pt onto the infinite line of the segment
        Vector v = curve_pt - p0;
        T t_seg = dir.dot(v) / (dir.squaredNorm());

        // Only consider points whose projection falls (approximately) within the segment bounds
        if (t_seg >= -tol.evaluateEpsilon(seg_length) && t_seg <= T(1) + tol.evaluateEpsilon(seg_length)) {
            Vector seg_pt = p0 + t_seg * dir;
            T dist = (curve_pt - seg_pt).norm();
            if (dist < tol.evaluateEpsilon(seg_length) * (T)1.0) {
                // Attempt local refinement on curve parameter to reduce distance to the segment
                T t_c = t_curve;
                constexpr int MAX_NEWTON_ITER = 16;

                for (int iter = 0; iter < MAX_NEWTON_ITER; ++iter) {
                    Vector cpt = curve.evaluate(t_c).coords;
                    Vector cder = curve.evaluateDerivative(t_c, tol); // use tolerance-aware derivative

                    // projection onto the segment
                    Vector v2 = cpt - p0;
                    T t_s = dir.dot(v2) / dir.squaredNorm();
                    t_s = std::clamp(t_s, T(0), T(1));
                    Vector spt = p0 + t_s * dir;
                    Vector diff = cpt - spt;
                    T dist2 = diff.norm();

                    if (dist2 < tol.evaluateEpsilon(seg_length)) {
                        break;
                    }

                    // Solve 2x2 normal equations for simultaneous update of (t_c, t_s):
                    // J = [c'(t)  -dir]
                    // (J^T J) * [dt; ds] = - J^T * f  where f = c(t) - s(t_s)
                    T a00 = cder.dot(cder);
                    T a01 = -cder.dot(dir);
                    T a11 = dir.dot(dir);

                    T b0 = -cder.dot(diff);
                    T b1 = dir.dot(diff);

                    T det = a00 * a11 - a01 * a01;
                    if (std::abs(det) < tol.evaluateEpsilon(seg_length)) {
                        // ill-conditioned; bail out
                        break;
                    }

                    T dt = (b0 * a11 - a01 * b1) / det;
                    T ds = (a00 * b1 - b0 * a01) / det;

                    // safeguard step sizes
                    T max_dt = (t_max - t_min) * (T)0.5;
                    if (std::abs(dt) > max_dt) dt = std::copysign(max_dt, dt);
                    if (std::abs(ds) > (T)0.5) ds = std::copysign((T)0.5, ds);

                    t_c += dt;
                    t_s += ds;

                    if (t_c < t_min) t_c = t_min;
                    if (t_c > t_max) t_c = t_max;
                    // keep t_s within segment bounds but allow a bit of overshoot for convergence
                    if (t_s < -tol.evaluateEpsilon(seg_length)) t_s = -tol.evaluateEpsilon(seg_length);
                    if (t_s > T(1) + tol.evaluateEpsilon(seg_length)) t_s = T(1) + tol.evaluateEpsilon(seg_length);

                    // recompute spt and diff for next iteration (loop does at top)
                }

                Vector final_pt = curve.evaluate(t_c).coords;
                // Final projection and distance check
                Vector vfinal = final_pt - p0;
                T final_tseg = dir.dot(vfinal) / dir.squaredNorm();
                if (final_tseg >= -tol.evaluateEpsilon(seg_length) && final_tseg <= T(1) + tol.evaluateEpsilon(seg_length)) {
                    Vector final_spt = p0 + std::clamp(final_tseg, T(0), T(1)) * dir;
                    T final_dist = (final_pt - final_spt).norm();
                    if (final_dist < tol.evaluateEpsilon(seg_length)) {
                        if (!is_duplicate(final_pt)) {
                            intersection_points.push_back(Point<T, N>{final_pt});
                        }
                    }
                }
            }
        }
    }

    if (!intersection_points.empty()) {
        return PointIntersectionResult<T, N>(true, "Segment and curve intersect", std::move(intersection_points));
    }
    return PointIntersectionResult<T, N>(false, "No intersection", {});
}


// Segment-Plane intersection
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Segment<T, N>& seg, const Plane<T, N>& plane,
                             const Tolerance& tol = Tolerance()) {
    using Vector = typename Plane<T, N>::VectorType;
    const Vector& p0 = seg.start.coords;
    const Vector& p1 = seg.end.coords;
    const Vector& n = plane.normal;
    const Vector& p_plane = plane.base.coords;
    T d0 = n.dot(p0 - p_plane);
    T d1 = n.dot(p1 - p_plane);
    T denom = d1 - d0;
    std::vector<Point<T, N>> result_points;
    // Check if both points are on the plane
    if (std::abs(d0) < tol.evaluateEpsilon(p0.norm()) && std::abs(d1) < tol.evaluateEpsilon(p1.norm())) {
        // The segment lies in the plane
        result_points.push_back(Point<T, N>{p0});
        result_points.push_back(Point<T, N>{p1});
        return PointIntersectionResult<T, N>(true, "Segment lies in the plane", std::move(result_points));
    }
    // If denom is zero, the segment is parallel to the plane
    if (std::abs(denom) < tol.evaluateEpsilon(std::max(std::abs(d0), std::abs(d1)))) {
        // No intersection (unless both points on plane, handled above)
        return PointIntersectionResult<T, N>(false, "Segment is parallel to the plane", {});
    }
    // Compute t such that p = (1-t)p0 + t*p1 is on the plane
    T t = -d0 / denom;
    if (t >= -tol.evaluateEpsilon(T(1)) && t <= 1 + tol.evaluateEpsilon(T(1))) {
        Vector intersection = p0 + (p1 - p0) * t;
        result_points.push_back(Point<T, N>{intersection});
        return PointIntersectionResult<T, N>(true, "Segment intersects plane", std::move(result_points));
    }
    return PointIntersectionResult<T, N>(false, "No intersection", {});
}

// Segment-Face intersection
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Segment<T, N>& seg, const Face<T, N>& face,
                             const Tolerance& tol = Tolerance()) {
    std::vector<Point<T, N>> intersection_points;

    // 1. Intersect with supporting plane
    Plane<T, N> plane(face.base, face.normal);
    auto plane_result = intersect(seg, plane, tol);
    if (!plane_result.intersects) {
        // The segment does not intersect the face's plane
        return PointIntersectionResult<T, N>(false, "No intersection", {});
    }

    // 2. For each intersection point with the plane, check if inside the face polygon
    for (const auto& pt : plane_result.points) {
        auto pt_in_face = intersect(pt, face, tol);
        if (pt_in_face.intersects) {
            intersection_points.push_back(pt);
        }
    }

    // 3. Clip segment endpoints to the face (if they are inside)
    auto r0 = intersect(seg.start, face, tol);
    auto r1 = intersect(seg.end, face, tol);
    if (r0.intersects) {
        for (const auto& pt : r0.points)
            intersection_points.push_back(pt);
    }
    if (r1.intersects) {
        for (const auto& pt : r1.points)
            intersection_points.push_back(pt);
    }

    // Remove duplicates (points may be found via both plane intersection and endpoint check)
    std::vector<Point<T, N>> unique_points;
    for (const auto& pt : intersection_points) {
        bool found = false;
        for (const auto& upt : unique_points) {
            if ((upt.coords - pt.coords).norm() < tol.evaluateEpsilon(pt.coords.norm())) {
                found = true;
                break;
            }
        }
        if (!found)
            unique_points.push_back(pt);
    }

    if (!unique_points.empty())
        return PointIntersectionResult<T, N>(true, "Segment intersects face", std::move(unique_points));
    return PointIntersectionResult<T, N>(false, "No intersection", {});
}

// Segment-Surface intersection (ND-compatible, uses line-surface intersection and bounds to segment)
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Segment<T, N>& seg, const Surface<T, N>& surface,
                             const Tolerance& tol = Tolerance()) {
    using Vector = Eigen::Matrix<T, N, 1>;
    std::vector<Point<T, N>> intersection_points;
    Vector p0 = seg.start.coords;
    Vector p1 = seg.end.coords;
    Vector dir = p1 - p0;
    T seg_length = dir.norm();
    if (seg_length < tol.evaluateEpsilon(seg_length)) {
        // Degenerate segment
        return PointIntersectionResult<T, N>(false, "Degenerate segment", {});
    }

    // Fast path for 3D: use Surface's specialized line intersection if available
    if constexpr (N == 3) {
        auto line_result = intersect(Line<T, 3>(seg.start, seg.end), surface, tol);
        if (!line_result.intersects) {
            // No intersection with infinite line, so none with segment
            return PointIntersectionResult<T, N>(false, "No intersection", {});
        }
        // For each intersection point, check if it lies on the segment
        for (const auto& pt : line_result.points) {
            Vector v = pt.coords - p0;
            T t = dir.dot(v) / dir.squaredNorm();
            if (t >= -tol.evaluateEpsilon(seg_length) && t <= T(1) + tol.evaluateEpsilon(seg_length)) {
                // Clamp to segment and check actual distance
                Vector segpt = p0 + std::clamp(t, T(0), T(1)) * dir;
                if ((pt.coords - segpt).norm() < tol.evaluateEpsilon(seg_length)) {
                    // Avoid duplicates
                    bool dup = false;
                    for (const auto& ex : intersection_points) {
                        if ((ex.coords - pt.coords).norm() < tol.evaluateEpsilon(seg_length)) { dup = true; break; }
                    }
                    if (!dup) intersection_points.push_back(pt);
                }
            }
        }
    } else {
        // ND: use line-surface intersection and bound to segment
        auto line_result = intersect(Line<T, N>(seg.start, seg.end), surface, tol);
        if (!line_result.intersects) {
            return PointIntersectionResult<T, N>(false, "No intersection", {});
        }
        for (const auto& pt : line_result.points) {
            Vector v = pt.coords - p0;
            T t = dir.dot(v) / dir.squaredNorm();
            if (t >= -tol.evaluateEpsilon(seg_length) && t <= T(1) + tol.evaluateEpsilon(seg_length)) {
                Vector segpt = p0 + std::clamp(t, T(0), T(1)) * dir;
                if ((pt.coords - segpt).norm() < tol.evaluateEpsilon(seg_length)) {
                    bool dup = false;
                    for (const auto& ex : intersection_points) {
                        if ((ex.coords - pt.coords).norm() < tol.evaluateEpsilon(seg_length)) { dup = true; break; }
                    }
                    if (!dup) intersection_points.push_back(pt);
                }
            }
        }
    }
    // Optionally also check endpoints for intersection, e.g. if endpoints lie on the surface
    auto r0 = intersect(seg.start, surface, tol);
    auto r1 = intersect(seg.end, surface, tol);
    if (r0.intersects) {
        for (const auto& pt : r0.points) {
            bool dup = false;
            for (const auto& ex : intersection_points) {
                if ((ex.coords - pt.coords).norm() < tol.evaluateEpsilon(seg_length)) { dup = true; break; }
            }
            if (!dup) intersection_points.push_back(pt);
        }
    }
    if (r1.intersects) {
        for (const auto& pt : r1.points) {
            bool dup = false;
            for (const auto& ex : intersection_points) {
                if ((ex.coords - pt.coords).norm() < tol.evaluateEpsilon(seg_length)) { dup = true; break; }
            }
            if (!dup) intersection_points.push_back(pt);
        }
    }
    if (!intersection_points.empty()) {
        return PointIntersectionResult<T, N>(true, "Segment intersects surface", std::move(intersection_points));
    }
    return PointIntersectionResult<T, N>(false, "No intersection", {});
}

// Symmetric overloads for reverse calls
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Curve<T, N>& curve, const Segment<T, N>& seg,
                             const Tolerance& tol = Tolerance()) {
    return intersect(seg, curve, tol);
}

template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Plane<T, N>& plane, const Segment<T, N>& seg,
                             const Tolerance& tol = Tolerance()) {
    return intersect(seg, plane, tol);
}

template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Face<T, N>& face, const Segment<T, N>& seg,
                             const Tolerance& tol = Tolerance()) {
    return intersect(seg, face, tol);
}

template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Surface<T, N>& surface, const Segment<T, N>& seg,
                             const Tolerance& tol = Tolerance()) {
    return intersect(seg, surface, tol);
}

} // namespace Euclid::Geometry
