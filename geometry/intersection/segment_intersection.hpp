#pragma once

#include "../geometry.hpp"
#include "intersection_result.hpp"
#include "geometry/tolerance.hpp"
#include "line_intersection.hpp"

#include <vector>
#include <Eigen/Core>

namespace Euclid::Geometry {


// Utility: project a point onto a segment and decide if it lies (within tolerance) in [0,1]
template <typename T, int N>
static inline bool point_on_segment_domain(const Eigen::Matrix<T, N, 1>& pt,
                                           const Segment<T, N>& seg,
                                           const Tolerance& tol) {
    using Vector = Eigen::Matrix<T, N, 1>;
    const Vector p0 = seg.start.coords;
    const Vector dir = seg.end.coords - seg.start.coords;
    const T seg_len = dir.norm();
    if (seg_len < tol.evaluateEpsilon(seg_len)) return false; // degenerate
    const T t = dir.dot(pt - p0) / dir.squaredNorm();
    return (t >= -tol.evaluateEpsilon(seg_len)) && (t <= T(1) + tol.evaluateEpsilon(seg_len));
}

// Utility: push a point into a vector if not a duplicate (w.r.t. tolerance)
template <typename T, int N>
static inline void push_unique_point(std::vector<Point<T, N>>& out,
                                     const Point<T, N>& pt,
                                     T seg_len,
                                     const Tolerance& tol) {
    for (const auto& ex : out) {
        if ((ex.coords - pt.coords).norm() < tol.evaluateEpsilon(seg_len)) return;
    }
    out.push_back(pt);
}

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

// Segment-Curve intersection using line–curve logic and segment-domain filter
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Segment<T, N>& seg, const Curve<T, N>& curve,
                             const Tolerance& tol = Tolerance()) {
    using Vector = Eigen::Matrix<T, N, 1>;
    std::vector<Point<T, N>> intersection_points;

    const Vector p0 = seg.start.coords;
    const Vector dir = seg.end.coords - seg.start.coords;
    const T seg_length = dir.norm();
    if (seg_length < tol.evaluateEpsilon(seg_length)) {
        return PointIntersectionResult<T, N>(false, "Degenerate segment", {});
    }

    // Delegate to the line–curve intersection using unit direction
    const Vector dir_unit = (seg_length > T(0)) ? (dir / seg_length) : dir;
    auto line_result = intersect(Line<T, N>(seg.start, dir_unit), curve, tol);
    if (!line_result.intersects) {
        return PointIntersectionResult<T, N>(false, "No intersection", {});
    }

    // Filter line intersections to the segment domain
    for (const auto& pt : line_result.points) {
        if (point_on_segment_domain<T, N>(pt.coords, seg, tol)) {
            // additionally clamp to the segment and verify distance (paranoia)
            const T t = dir.dot(pt.coords - p0) / dir.squaredNorm();
            const Vector segpt = p0 + std::clamp(t, T(0), T(1)) * dir;
            if ((pt.coords - segpt).norm() < tol.evaluateEpsilon(seg_length)) {
                push_unique_point<T, N>(intersection_points, pt, seg_length, tol);
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
    // NOTE: Could be simplified by delegating to Line–Plane and filtering to [0,1], kept explicit for clarity and identical behavior to prior tests.
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

    // Reuse Line–Surface intersection and clip hits to the segment domain [0,1]
    const Vector p0 = seg.start.coords;
    const Vector dir = seg.end.coords - seg.start.coords;
    const T seg_len = dir.norm();

    if (seg_len < tol.evaluateEpsilon(seg_len)) {
        return PointIntersectionResult<T, N>(false, "Degenerate segment", {});
    }

    /*
    // [DEBUG] print direction info
    std::cout << "[DEBUG] Using line direction = (" << dir.transpose() << "), norm = " << dir.norm() << std::endl;
     */
    const Vector dir_unit = dir / seg_len;
    auto line_result = ::Euclid::Geometry::intersect(Line<T, N>(seg.start, dir_unit), surface, tol);

    /*
    // [DEBUG] print line_result
    std::cout << "[DEBUG] Segment–Surface: line_result.intersects=" << line_result.intersects
              << ", point count=" << line_result.points.size() << std::endl;
     */

    if (line_result.intersects) {
        for (const auto& pt : line_result.points) {
            // [DEBUG] print candidate point
            //std::cout << "[DEBUG] Candidate point: (" << pt.coords.transpose() << ")" << std::endl;
            // Project to segment parameter t in [0,1] using normalized direction
            const T t = dir_unit.dot(pt.coords - p0) / seg_len;
            // [DEBUG] print computed t
            //std::cout << "[DEBUG] Computed t=" << t << ", seg_len=" << seg_len << std::endl;
            if (t >= -tol.evaluateEpsilon(seg_len) && t <= T(1) + tol.evaluateEpsilon(seg_len)) {
                // Clamp and verify proximity to the clamped point (paranoia check)
                const Vector seg_pt = p0 + std::clamp(t, T(0), T(1)) * dir;
                if ((pt.coords - seg_pt).norm() <= tol.evaluateEpsilon(seg_len) * T(10)) {
                    // [DEBUG] print accepted point
                    //std::cout << "[DEBUG] Accepted point: (" << pt.coords.transpose() << ")" << std::endl;
                    push_unique_point<T, N>(intersection_points, pt, seg_len, tol);
                }
            }
        }
    }

    // [DEBUG] after line–surface pass
    //std::cout << "[DEBUG] Line–Surface pass complete, accepted count=" << intersection_points.size() << std::endl;

    // Include segment endpoints that lie on the surface (handles grazing/endpoint cases)
    auto r0 = intersect(seg.start, surface, tol);
    if (r0.intersects) {
        for (const auto& pt : r0.points) {
            push_unique_point<T, N>(intersection_points, pt, seg_len, tol);
        }
    }
    auto r1 = intersect(seg.end, surface, tol);
    if (r1.intersects) {
        for (const auto& pt : r1.points) {
            push_unique_point<T, N>(intersection_points, pt, seg_len, tol);
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
