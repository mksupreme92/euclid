#pragma once

#include <Eigen/Core>
#include "point.hpp"
#include "tolerance.hpp"

namespace Euclid::Geometry {

template <typename Scalar, int Dim = Eigen::Dynamic>
class Line {
public:
    using PointType = Point<Scalar, Dim>;
    using VectorType = Eigen::Matrix<Scalar, Dim, 1>;

    // Construct from two points
    Line(const PointType& p1, const PointType& p2)
        : m_point1(p1), m_point2(p2)
    {
        m_direction = (p2.coords - p1.coords);
        m_direction.normalize();
    }

    // Construct from a point and a direction vector
    Line(const PointType& p, const VectorType& dir, bool normalize = true)
        : m_point1(p), m_direction(dir)
    {
        if (normalize)
            m_direction.normalize();
        // create a dummy second point along the line
        m_point2 = PointType(p.coords + m_direction);
    }

    // Accessors
    const PointType& point1() const { return m_point1; }
    const PointType& point2() const { return m_point2; }
    const VectorType& direction() const { return m_direction; }

    // Apply a general transform to the line
    template <typename Transform>
    Line applyTransform(const Transform& T) const {
        // Evaluate any Eigen expressions to avoid lazy evaluation issues
        PointType newPoint = T.apply(m_point1);
        VectorType newDir  = T.applyLinear(m_direction).eval(); // force evaluation
        return Line(newPoint, newDir, true);         // normalize direction
    }

    // Measure angle (radians) between this line and another line
    Scalar measureAngle(const Line& other) const {
        const auto& d1 = this->direction();
        const auto& d2 = other.direction();
        Scalar cosTheta = d1.dot(d2) / (d1.norm() * d2.norm());
        return std::acos(std::clamp(cosTheta, Scalar(-1), Scalar(1)));
    }

    template <typename Scalar2 = Scalar, int Dim2 = Dim>
    struct LineIntersection {
        bool intersects;
        Point<Scalar2, Dim2> pointOnThis;
        Point<Scalar2, Dim2> pointOnOther;
        Scalar2 distance;
    };

    // Intersection method for ND lines (Dim >= 2)
    LineIntersection<Scalar, Dim> intersect(const Line<Scalar, Dim>& other) const {
        const auto& p1 = m_point1.coords;
        const auto& d1 = m_direction;
        const auto& p2 = other.m_point1.coords;
        const auto& d2 = other.m_direction;

        // Use dynamic-size Eigen matrix to avoid thin/full SVD crash
        Eigen::Matrix<Scalar, Eigen::Dynamic, 2> A(Dim, 2);
        A.col(0) = d1;
        A.col(1) = -d2;
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> b = p2 - p1;

        // Solve least-squares for t and s
        Eigen::Matrix<Scalar, 2, 1> ts = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);

        // Candidate points along each line
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> point_on_this  = p1 + ts(0) * d1;
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> point_on_other = p2 + ts(1) * d2;

        Scalar dist = (point_on_this - point_on_other).norm();
        Tolerance tolObj;
        Scalar tol = tolObj.evaluateEpsilon(static_cast<double>(Dim));
        bool does_intersect = dist <= tol;

        return {does_intersect,
                Point<Scalar, Dim>(point_on_this),
                Point<Scalar, Dim>(point_on_other),
                dist};
    }

private:
    PointType m_point1;
    PointType m_point2;
    VectorType m_direction;
};

} // namespace euclid::geometry
