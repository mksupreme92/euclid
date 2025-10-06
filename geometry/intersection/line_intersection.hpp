#pragma once

#include "../geometry.hpp"
#include "../tolerance.hpp"
#include "intersection_result.hpp"
#include <optional>
#include <vector>
#include <string>

namespace Euclid::Geometry {

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

    // Generic ND implementation: check if lines are parallel or coincident
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
    Eigen::Matrix<T, 2, 1> b;
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

// ======================================================
// Line–Plane Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& line, const Plane<T, N>& plane,
                                      const Tolerance& tol = Tolerance()) {
    LineIntersectionResult<T, N> result;
    if constexpr (N == 3) {
        auto n = plane.normal;
        auto diff = line.point1().coords - plane.base.coords;
        T denom = n.dot(line.direction());
        T scale = line.direction().norm() + plane.normal.norm();
        if (std::abs(denom) < tol.evaluateEpsilon(scale)) {
            result.description = "Line is parallel to plane";
            return result;
        }
        T t = -n.dot(diff) / denom;
        Point<T, N> pt(line.point1().coords + t * line.direction());
        result.intersects = true;
        result.points = {pt};
        result.description = "Line intersects plane at a point";
        return result;
    }
    result.description = "ND Line–Plane intersection not implemented";
    return result;
}

// Reverse overload: Plane–Line
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Plane<T, N>& plane, const Line<T, N>& line,
                                      const Tolerance& tol = Tolerance()) {
    return intersect(line, plane, tol);
}

} // namespace Euclid::Geometry
