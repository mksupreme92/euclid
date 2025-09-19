#pragma once

#include "config.hpp"
#include <Eigen/Dense>
#include <cmath>
#include "../algebra/metric.hpp"

namespace euclid::geometry {

template <typename T, int N>
struct Point {
    Eigen::Matrix<T, N, 1> coords;

    // === Constructors ===
    Point() : coords(Eigen::Matrix<T, N, 1>::Zero()) {
        if constexpr (N == Eigen::Dynamic) {
            assert(coords.size() == Euclid::getSpaceDimension() && "Point dimension mismatch");
        } else {
            assert(N == Euclid::getSpaceDimension() && "Point dimension mismatch");
        }
    }
    explicit Point(const Eigen::Matrix<T, N, 1>& v) : coords(v) {
        if constexpr (N == Eigen::Dynamic) {
            assert(coords.size() == Euclid::getSpaceDimension() && "Point dimension mismatch");
        } else {
            assert(N == Euclid::getSpaceDimension() && "Point dimension mismatch");
        }
    }

    // From initializer list
    Point(std::initializer_list<T> list) {
        assert(list.size() == Euclid::getSpaceDimension() && "Initializer list size mismatch for Point");
        int i = 0;
        for (auto v : list) coords[i++] = v;
        if constexpr (N == Eigen::Dynamic) {
            assert(coords.size() == Euclid::getSpaceDimension() && "Point dimension mismatch");
        } else {
            assert(N == Euclid::getSpaceDimension() && "Point dimension mismatch");
        }
    }

    // === Access ===
    T& operator[](size_t i) { return coords[i]; }
    const T& operator[](size_t i) const { return coords[i]; }

    // === Equality ===
    bool operator==(const Point& other) const {
        return coords.isApprox(other.coords);
    }
    bool operator!=(const Point& other) const { return !(*this == other); }

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
};

// Template deduction guide for dynamic points (Eigen::Dynamic)
template <typename T, typename... Ts>
Point(T, Ts...) -> Point<T, Eigen::Dynamic>;

} // namespace euclid::geometry
