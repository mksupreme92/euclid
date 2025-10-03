#pragma once
#include <Eigen/Dense>
#include <array>
#include <optional>
#include <cmath>
#include <vector>
#include "point.hpp"
#include "tolerance.hpp"

namespace Euclid::Geometry{

template<typename T, int N>
class Plane {
public:
    using VectorType = Eigen::Matrix<T, N, 1>;
    
    Point<T, N> base;
    VectorType normal;
    std::vector<Point<T, N>> definingPoints;
    
    Plane() = default;
    
    Plane(const Point<T, N>& base_, const VectorType& normal_)
    : base(base_), normal(normal_) {}
    
    static std::optional<Plane<T, N>> fromPoints(const std::array<Point<T, N>, N>& pts) {
        static_assert(N >= 2, "Plane requires at least 2 points");
        
        // For N=2, normal is perpendicular to the vector p1 - p0
        if constexpr (N == 2) {
            VectorType v = pts[1].coords - pts[0].coords;
            // For 2D, the normal vector perpendicular to (x, y) is (-y, x)
            VectorType n;
            n(0) = -v(1);
            n(1) = v(0);
            // Normalize normal vector
            T norm_val = n.norm();
            if (norm_val == T(0)) {
                return std::nullopt;
            }
            n /= norm_val;
            Plane<T, N> plane(pts[0], n);
            plane.definingPoints.assign(pts.begin(), pts.end());
            return plane;
        } else {
            // For N >= 3, compute the normal vector as the nullspace of the matrix formed by vectors (p_i - p_0)
            Eigen::Matrix<T, N, N - 1> mat;
            for (int i = 1; i < N; ++i) {
                mat.col(i - 1) = pts[i].coords - pts[0].coords;
            }
            
            // Compute QR decomposition
            Eigen::HouseholderQR<Eigen::Matrix<T, N, N - 1>> qr(mat);
            Eigen::Matrix<T, N, N> Q = qr.householderQ();
            
            // The normal vector is the last column of Q (orthogonal complement)
            VectorType n = Q.col(N - 1);
            
            // Normalize normal vector
            T norm_val = n.norm();
            if (norm_val == T(0)) {
                return std::nullopt;
            }
            n /= norm_val;
            
            Plane<T, N> plane(pts[0], n);
            plane.definingPoints.assign(pts.begin(), pts.end());
            return plane;
        }
    }
    
    static std::optional<Plane<T, N>> fromPointAndNormal(const Point<T, N>& base, const VectorType& normalVec) {
        T norm_val = normalVec.norm();
        if (norm_val == T(0)) {
            return std::nullopt;
        }
        VectorType n = normalVec / norm_val;
        return Plane<T, N>(base, n);
    }
    
    bool contains(const Point<T, N>& pt, const Euclid::Tolerance& tol = Euclid::Tolerance()) const {
        VectorType diff = pt.coords - base.coords;
        T dot_product = normal.dot(diff);
        double scale = std::max(diff.norm(), T(1));  // or any appropriate scale for this context
        return Euclid::equalWithinTolerance(dot_product, T(0), tol, scale);
    }
    
    T signedDistance(const Point<T, N>& pt) const {
        VectorType diff = pt.coords - base.coords;
        return normal.dot(diff);
    }
    
    template<typename Transform>
    Plane<T, N> applyTransform(const Transform& transform) const {
        if (!definingPoints.empty()) {
            std::array<Point<T, N>, N> transformedPoints;
            for (size_t i = 0; i < N; ++i) {
                transformedPoints[i] = transform.apply(definingPoints[i]);
            }
            return fromPoints(transformedPoints).value();
        } else {
            Point<T, N> newBase = transform.apply(base);
            VectorType newNormal = transform.getLinear().template cast<T>() * normal;
            newNormal.normalize();
            return Plane<T, N>(newBase, newNormal);
        }
    }
};
};

