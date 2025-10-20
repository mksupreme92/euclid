#pragma once
#include <functional>
#include <utility>
#include "point.hpp"
#include "tolerance.hpp"
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

    // Evaluate the derivative of the curve at parameter t (finite difference)
    Eigen::Matrix<Scalar, Dim, 1> evaluateDerivative(Scalar t, const Tolerance& tol) const {
        Scalar h = tol.evaluateEpsilon((curveFunc_(t)).coords.norm());
        return (curveFunc_(t + h).coords - curveFunc_(t - h).coords) / (2.0 * h);
    }

    // Estimate integral of curvature magnitude across the domain using adaptive sampling
    Scalar evaluateIntegral(const Tolerance& tol, int baseDivs = 64, int maxDepth = 10) const {
        Scalar t0 = domain_.first;
        Scalar t1 = domain_.second;
        std::vector<Scalar> samples;
        samples.reserve(baseDivs + 1);
        for (int i = 0; i <= baseDivs; ++i)
            samples.push_back(t0 + (t1 - t0) * (Scalar(i) / baseDivs));

        auto curvatureAt = [&](Scalar t) {
            const Scalar h = tol.evaluateEpsilon(curveFunc_(t).coords.norm());
            auto p0 = curveFunc_(t).coords;
            auto p1 = curveFunc_(std::min(t1, t + h)).coords;
            auto p_1 = curveFunc_(std::max(t0, t - h)).coords;
            auto first = (p1 - p_1) / (2 * h);
            auto second = (p1 - 2 * p0 + p_1) / (h * h);
            return second.norm() / std::pow(1.0 + first.squaredNorm(), 1.5);
        };

        // Adaptive subdivision loop
        for (int depth = 0; depth < maxDepth; ++depth) {
            std::vector<Scalar> newSamples;
            for (size_t i = 0; i + 1 < samples.size(); ++i) {
                Scalar a = samples[i];
                Scalar b = samples[i + 1];
                Scalar mid = (a + b) * 0.5;
                Scalar kA = curvatureAt(a);
                Scalar kB = curvatureAt(b);
                Scalar kM = curvatureAt(mid);
                Scalar avgK = (kA + kB + 4 * kM) / 6;
                if (avgK * std::abs(b - a) > tol.paramTol)
                    newSamples.push_back(mid);
            }
            if (newSamples.empty()) break;
            samples.insert(samples.end(), newSamples.begin(), newSamples.end());
            std::sort(samples.begin(), samples.end());
            samples.erase(std::unique(samples.begin(), samples.end(), [&](Scalar x, Scalar y){ return std::abs(x - y) < tol.paramTol; }), samples.end());
        }

        // Trapezoidal integration of curvature magnitude
        Scalar integral = 0;
        for (size_t i = 0; i + 1 < samples.size(); ++i) {
            Scalar a = samples[i];
            Scalar b = samples[i + 1];
            Scalar kA = curvatureAt(a);
            Scalar kB = curvatureAt(b);
            integral += (kA + kB) * 0.5 * (b - a);
        }

        return integral;
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
