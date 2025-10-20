#pragma once

#include <vector>
#include <string>

namespace Euclid::Geometry {

template <typename T, int N>
struct Point;

template <typename T, int N>
struct ParamHit {
    T t_line{};
    T u_curve{};
    T v_surface{};
    Point<T, N> p{};
    bool tangential = false;
};

template <typename T, int N>
struct IntersectionResult {
    bool intersects = false;
    std::vector<Point<T, N>> points;
    std::vector<ParamHit<T, N>> hits;
    std::string description;

    IntersectionResult() = default;

    IntersectionResult(bool inter, std::string desc = {})
        : intersects(inter), description(std::move(desc)) {}

    void addHit(const ParamHit<T, N>& h) {
        hits.push_back(h);
        points.push_back(h.p);
        intersects = true;
    }
};

} // namespace Euclid::Geometry
