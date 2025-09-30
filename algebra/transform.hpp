#pragma once
#include <Eigen/Dense>
#include <stdexcept>
#include <cmath>

namespace euclid::algebra {

// ---------------------------
// Linear (matrix only)
// ---------------------------
template <typename T, int N = Eigen::Dynamic>
struct Linear {
    Eigen::Matrix<T,N,N> A;

    Linear() : A(Eigen::Matrix<T,N,N>::Identity()) {}
    explicit Linear(const Eigen::Matrix<T,N,N>& mat) : A(mat) {}

    // Apply to vector
    Eigen::Matrix<T,N,1> apply(const Eigen::Matrix<T,N,1>& v) const {
        return A * v;
    }
};

// ---------------------------
// Affine (linear + translation)
// ---------------------------
template <typename T, int N = Eigen::Dynamic>
struct Affine : Linear<T,N> {
    Eigen::Matrix<T,N,1> b;

    Affine() : Linear<T,N>(), b(Eigen::Matrix<T,N,1>::Zero()) {}
    Affine(const Eigen::Matrix<T,N,N>& mat, const Eigen::Matrix<T,N,1>& t)
        : Linear<T,N>(mat), b(t) {}

    Eigen::Matrix<T,N,N> getLinear() const { return this->A; }
    Eigen::Matrix<T,N,1> getTranslation() const { return b; }

    // Apply to point
    Eigen::Matrix<T,N,1> apply(const Eigen::Matrix<T,N,1>& p) const {
        return this->A * p + b;
    }

    // Apply to Point<T,N>
    template <typename PointT>
    PointT apply(const PointT& p) const {
        return PointT(this->A * p.coords + b);
    }

    // Apply to vector (ignores translation)
    Eigen::Matrix<T,N,1> applyLinear(const Eigen::Matrix<T,N,1>& v) const {
        return Linear<T,N>::apply(v);
    }

    // Composition
    Affine operator*(const Affine& other) const {
        return { this->A * other.A, this->A * other.b + b };
    }

    // Inverse
    Affine inverse() const {
        Eigen::Matrix<T,N,N> A_inv = this->A.inverse();
        return { A_inv, -(A_inv * b) };
    }

    // Convenience constructors
    static Affine translate(const Eigen::Matrix<T,N,1>& t) {
        return Affine(Eigen::Matrix<T,N,N>::Identity(), t);
    }

    static Affine scale(T s) {
        return Affine(Eigen::Matrix<T,N,N>::Identity() * s, Eigen::Matrix<T,N,1>::Zero());
    }

    static Affine scale(const Eigen::Matrix<T,N,1>& s) {
        Eigen::Matrix<T,N,N> S = Eigen::Matrix<T,N,N>::Zero();
        for (typename Eigen::Matrix<T,N,1>::Index i = 0; i < s.size(); ++i) {
            S(i,i) = s(i);
        }
        return Affine(S, Eigen::Matrix<T,N,1>::Zero());
    }

};

// ---------------------------
// Orthogonal (metric-preserving)
// ---------------------------
template <typename T, int N = Eigen::Dynamic>
struct Orthogonal : Linear<T,N> {
    Orthogonal(const Eigen::Matrix<T,N,N>& mat) : Linear<T,N>(mat) {
        // optionally enforce AᵀA = I
    }

    // determinant check for orientation
    T det() const { return this->A.determinant(); }

    template <typename PointT>
    PointT apply(const PointT& p) const {
        return PointT(this->A * p.coords);
    }
};

// ---------------------------
// Special Orthogonal (rotation only)
// ---------------------------
template <typename T, int N = Eigen::Dynamic>
struct SpecialOrthogonal : Orthogonal<T,N> {
    SpecialOrthogonal(const Eigen::Matrix<T,N,N>& mat) : Orthogonal<T,N>(mat) {
        if (std::abs(mat.determinant() - T(1)) > 1e-6)
            throw std::runtime_error("Not a valid SO(N) matrix");
    }

    template <typename PointT>
    PointT apply(const PointT& p) const {
        return PointT(this->A * p.coords);
    }
};

// ---------------------------
// Projective (N+1 x N+1 homogeneous)
// ---------------------------
template <typename T, int N = Eigen::Dynamic>
struct Projective {
    Eigen::Matrix<T,N+1,N+1> P;

    Projective() : P(Eigen::Matrix<T,N+1,N+1>::Identity()) {}

    Eigen::Matrix<T,N,1> apply(const Eigen::Matrix<T,N,1>& p) const {
        Eigen::Matrix<T,N+1,1> hp;
        hp << p, 1;
        Eigen::Matrix<T,N+1,1> result = P * hp;
        return result.template head<N>() / result(N);
    }
};

} // namespace euclid::algebra
