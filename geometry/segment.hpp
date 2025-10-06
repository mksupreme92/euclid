#pragma once
#include "point.hpp"
#include <Eigen/Dense>

namespace Euclid::Geometry {

template <typename T = double, int N = Eigen::Dynamic>
struct Segment {
    Point<T, N> start;
    Point<T, N> end;

    // Default constructor
    Segment() = default;

    // Constructor from two points
    Segment(const Point<T, N>& p1, const Point<T, N>& p2)
        : start(p1), end(p2) {}

    // Length of the segment
    T length() const {
        return (end.coords - start.coords).norm();
    }

    // Midpoint of the segment
    Point<T, N> midpoint() const {
        return Point<T, N>((start.coords + end.coords) * T(0.5));
    }

    // Vector from start to end
    Eigen::Matrix<T, N, 1> vector() const {
        return end.coords - start.coords;
    }

    // Check if a point lies on the segment (within tolerance)
    bool contains(const Point<T, N>& pt, const Tolerance& tol = Tolerance()) const {
        Eigen::Matrix<T, N, 1> v = vector();
        Eigen::Matrix<T, N, 1> w = pt.coords - start.coords;
        T lenSquared = v.squaredNorm();
        if (lenSquared == T(0)) return (pt.coords - start.coords).norm() <= tol.evaluateEpsilon(N);
        T t = v.dot(w) / lenSquared;
        if (t < T(0) || t > T(1)) return false;
        Point<T, N> projection(start.coords + t * v);
        return (projection.coords - pt.coords).norm() <= tol.evaluateEpsilon(N);
    }

    // Apply a linear transform (Eigen matrix)
    void applyTransform(const Eigen::Matrix<T, N, N>& mat) {
        start.coords = mat * start.coords;
        end.coords = mat * end.coords;
    }
};

} // namespace Euclid::Geometry
