#pragma once

#include <Eigen/Core>
#include "point.hpp"

namespace euclid::geometry {

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
    };

    // Intersection method for 2D lines
    template <int D = Dim>
    typename std::enable_if<(D == 2), LineIntersection<Scalar, 2>>::type
    intersect(const Line<Scalar, 2>& other) const {
        const auto& p = m_point1.coords;
        const auto& r = m_direction;
        const auto& q = other.m_point1.coords;
        const auto& s = other.m_direction;

        Scalar r_cross_s = r.x() * s.y() - r.y() * s.x();
        if (std::abs(r_cross_s) < Scalar(1e-12)) {
            // Lines are parallel or coincident
            return {false, Point<Scalar, 2>(Eigen::Matrix<Scalar, 2, 1>::Zero())};
        }

        Eigen::Matrix<Scalar, 2, 1> qp = q - p;
        Scalar t = (qp.x() * s.y() - qp.y() * s.x()) / r_cross_s;

        Point<Scalar, 2> intersection_point(p + t * r);
        return {true, intersection_point};
    }

    /*
    // TODO: ND intersection is disabled for now
    // Intersection method for ND lines (Dim > 2)
    template <int D = Dim>
    typename std::enable_if<(D > 2 || D == Eigen::Dynamic), LineIntersection<Scalar, Dim>>::type
    intersect(const Line<Scalar, Dim>& other) const {
        const auto& p1 = m_point1.coords;
        const auto& d1 = m_direction;
        const auto& p2 = other.m_point1.coords;
        const auto& d2 = other.m_direction;

        // Use dynamic-size Eigen matrix to avoid thin/full SVD crash
        Eigen::Matrix<Scalar, Eigen::Dynamic, 2> A(D, 2);
        A.col(0) = d1;
        A.col(1) = -d2;
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> b = p2 - p1;

        // Solve least-squares for t and s
        Eigen::Matrix<Scalar, 2, 1> ts = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);

        // Candidate points along each line
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> point_on_this  = p1 + ts(0) * d1;
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> point_on_other = p2 + ts(1) * d2;

        // Residual check for true intersection
        Scalar tol = Scalar(1e-12) * static_cast<Scalar>(D);
        if ((point_on_this - point_on_other).norm() > tol) {
            // Lines do not intersect
            return {false, Point<Scalar, Dim>()};
        }

        // Lines intersect
        return {true, Point<Scalar, Dim>(point_on_this)};
    }
    */

private:
    PointType m_point1;
    PointType m_point2;
    VectorType m_direction;
};

} // namespace euclid::geometry
