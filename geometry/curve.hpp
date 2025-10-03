

#pragma once
#include <functional>
#include <utility>
#include "point.hpp"
#include "algebra/transform.hpp"

namespace Euclid::Geometry {

template <typename Scalar, int Dim = Eigen::Dynamic>
class Curve {
public:
    using PointType = Point<Scalar, Dim>;

    // General parametric curve constructor
    Curve(std::function<PointType(Scalar)> func, Scalar t0, Scalar t1)
        : curveFunc_(func), domain_{t0, t1} {}

    // Evaluate the curve at parameter t
    PointType evaluate(Scalar t) const {
        return curveFunc_(t);
    }

    // Curve parameter domain
    std::pair<Scalar, Scalar> domain() const { return domain_; }

    // Apply a transform to the curve
    template <typename Transform>
    Curve applyTransform(const Transform& T) const {
        auto newFunc = [T, this](Scalar t) -> PointType {
            return T.apply(this->evaluate(t));
        };
        return Curve(newFunc, domain_.first, domain_.second);
    }

    // Linear curve factory: p(t) = p0 + t*(p1 - p0)
    static Curve linearCurve(const PointType& p0, const PointType& p1) {
        auto func = [p0, p1](Scalar t) -> PointType {
            return p0 + (p1 - p0) * t;
        };
        return Curve(func, 0, 1);
    }

private:
    std::function<PointType(Scalar)> curveFunc_;
    std::pair<Scalar, Scalar> domain_;
};

} // namespace Euclid::Geometry
