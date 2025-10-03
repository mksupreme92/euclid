#pragma once

#include <Eigen/Dense>
#include <cmath>
#include "tolerance.hpp"

namespace Euclid::Geometry {

template <typename T, int N>
struct Point {
    Eigen::Matrix<T, N, 1> coords;

    // === Constructors ===
    Point() : coords(Eigen::Matrix<T, N, 1>::Zero()) {
    }
    explicit Point(const Eigen::Matrix<T, N, 1>& v) : coords(v) {
    }

    // From initializer list
    Point(std::initializer_list<T> list) {
        int i = 0;
        for (auto v : list) coords[i++] = v;
    }

    // === Access ===
    T& operator[](size_t i) { return coords[i]; }
    const T& operator[](size_t i) const { return coords[i]; }

    // === Equality ===
    bool operator==(const Point& other) const {
        return isEqual(other);
    }
    bool operator!=(const Point& other) const { return !(*this == other); }

    bool isEqual(const Point& other, const Euclid::Tolerance& tol = Euclid::Tolerance()) const {
        T scale = std::max(coordinateMagnitude(), other.coordinateMagnitude());
        for (int i = 0; i < coords.size(); ++i) {
            if (!Euclid::equalWithinTolerance(coords[i], other.coords[i], tol, scale)) {
                return false;
            }
        }
        return true;
    }

    // === Arithmetic ===
    Point operator+(const Eigen::Matrix<T, N, 1>& v) const {
        return Point(coords + v);
    }
    Point operator-(const Eigen::Matrix<T, N, 1>& v) const {
        return Point(coords - v);
    }
    Eigen::Matrix<T, N, 1> operator-(const Point& other) const {
        return coords - other.coords; // yields a vector
    }

    // === Measurements ===
    T distanceTo(const Point& other) const {
        return (coords - other.coords).norm();
    }

    Point midpoint(const Point& other) const {
        return Point((coords + other.coords) / T(2));
    }

    // === Additional methods ===
    T coordinateMagnitude() const {
        return coords.norm();
    }

};

// Template deduction guide for dynamic points (Eigen::Dynamic)
template <typename T, typename... Ts>
Point(T, Ts...) -> Point<T, Eigen::Dynamic>;

} // namespace Euclid::Geometry
