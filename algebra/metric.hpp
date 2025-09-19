#pragma once
#include <Eigen/Dense>
#include <cmath>

namespace euclid::algebra {

// --- Inner product ⟨v, w⟩ for Euclidean metric ---
template <typename T, int N>
T innerProduct(const Eigen::Matrix<T, N, 1>& v, const Eigen::Matrix<T, N, 1>& w) {
    return v.dot(w);
}

// --- Euclidean distance between vectors v and w ---
template <typename T, int N>
T distance(const Eigen::Matrix<T, N, 1>& v, const Eigen::Matrix<T, N, 1>& w) {
    return (v - w).norm();
}

} // namespace euclid::algebra
