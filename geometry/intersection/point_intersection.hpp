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
// Point–Curve Intersection (closest point search)
// ======================================================
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Point<T, N>& p, const Curve<T, N>& c,
                             const Tolerance& tol = Tolerance()) {
    // Placeholder: requires iterative minimization
    // Idea: minimize ||c(t) - p|| over t ∈ [t0, t1]
    // For now, return an empty placeholder.
    return PointIntersectionResult<T, N>(false, "Point–Curve intersection not yet implemented");
}

// ======================================================
// Point–Surface Intersection (projection / closest point search)
// ======================================================
template <typename T, int N>
PointIntersectionResult<T, N> intersect(const Point<T, N>& p, const Surface<T, N>& s,
                             const Tolerance& tol = Tolerance()) {
    // Placeholder: requires 2D parametric minimization
    // Idea: minimize ||S(u, v) - p|| over domain
    // For now, return an empty placeholder.
    return PointIntersectionResult<T, N>(false, "Point–Surface intersection not yet implemented");
}

} // namespace Euclid::Geometry
