#pragma once

#include <vector>
#include <stdexcept>
#include <Eigen/Dense>
#include "point.hpp"

namespace Euclid::Geometry {

template<typename T, int N>
class Face {
public:
    using PointType = Point<T, N>;            // Mirror Plane API
    using VectorType = Eigen::Matrix<T, N, 1>;
    using Points = std::vector<PointType>;

    PointType base;
    VectorType normal;
    Points vertices;

    Face() = default;

    Face(const PointType& base_, const VectorType& normal_, const Points& verts)
        : base{base_}, normal{normal_}, vertices{verts} {}

    // Construct from a list of points (vertices)
    static Face fromPoints(const Points& pts) {
        if (pts.size() < 2)
            throw std::invalid_argument("Face requires at least 2 points");

        int dim = N;
        for (const auto& p : pts)
            if (int(p.coords.size()) != dim)
                throw std::invalid_argument("All points must have same dimension");

        VectorType n = computeNormal(pts, dim);
        n.normalize();
        PointType basePt{pts[0]};
        return Face{basePt, n, pts};
    }

    // Construct from base + normal + optional vertices
    static Face fromPointAndNormal(const PointType& base, const VectorType& normalVec, const Points& verts = {}) {
        T norm_val = normalVec.norm();
        if (norm_val == T(0))
            throw std::invalid_argument("Normal vector cannot be zero");
        VectorType n = normalVec / norm_val;
        return Face{base, n, verts};
    }

    // Check if a point is inside the face
    bool contains(const PointType& p, T tol = T(1e-8)) const {
        if (vertices.size() < 2)
            return false;

        if (std::abs(normal.dot(p.coords - base.coords)) > tol)
            return false;

        // Simple polygon/convex hull check in N dimensions
        VectorType centroid = VectorType::Zero(N);
        for (const auto& v : vertices)
            centroid += v.coords;
        centroid /= T(vertices.size());

        for (size_t i = 0; i < vertices.size(); ++i) {
            VectorType edge = vertices[(i + 1) % vertices.size()].coords - vertices[i].coords;
            VectorType edgeNormal = anyPerpendicular(edge);
            edgeNormal.normalize();

            T signCentroid = edgeNormal.dot(centroid - vertices[i].coords);
            T signPoint = edgeNormal.dot(p.coords - vertices[i].coords);

            if (signCentroid * signPoint < -tol)
                return false;
        }

        return true;
    }

    template<typename Transform>
    Face<T, N> applyTransform(const Transform& transform) const {
        if (!vertices.empty()) {
            Points transformedVerts;
            transformedVerts.reserve(vertices.size());
            for (const auto& vertex : vertices) {
                transformedVerts.push_back(transform.apply(vertex));
            }
            return Face::fromPoints(transformedVerts);
        } else {
            PointType transformedBase = transform.apply(base);
            VectorType transformedNormal = transform.applyLinear(normal).normalized();
            return Face::fromPointAndNormal(transformedBase, transformedNormal);
        }
    }

private:
    // Compute normal given vertices
    static VectorType computeNormal(const Points& pts, int dim) {
        if (dim == 2) {
            VectorType v = pts[1].coords - pts[0].coords;
            VectorType n;
            n(0) = -v(1);
            n(1) = v(0);
            return n;
        }
        if (int(pts.size()) < dim)
            throw std::invalid_argument("Not enough points to define hyperplane");

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A(dim, dim - 1);
        for (int j = 1; j < dim; ++j)
            A.col(j - 1) = pts[j].coords - pts[0].coords;

        Eigen::HouseholderQR<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> qr(A);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Q = qr.householderQ();
        return Q.col(dim - 1);
    }

    // Return a vector perpendicular to v
    static VectorType anyPerpendicular(const VectorType& v) {
        int dim = int(v.size());
        VectorType perp = VectorType::Zero(dim);

        int idx = 0;
        T minAbs = std::abs(v(0));
        for (int i = 1; i < dim; ++i) {
            T absVal = std::abs(v(i));
            if (absVal < minAbs) {
                minAbs = absVal;
                idx = i;
            }
        }
        perp(idx) = T(1);

        T vnorm2 = v.squaredNorm();
        if (vnorm2 > T(1e-16))
            perp -= (perp.dot(v) / vnorm2) * v;

        T perpNorm = perp.norm();
        if (perpNorm > T(0))
            perp /= perpNorm;
        else
            perp.setZero();

        return perp;
    }
};

} // namespace Euclid::Geometry
