#pragma once

#include "../geometry.hpp"
#include "../tolerance.hpp"
#include "intersection_result.hpp"
#include <optional>
#include <vector>
#include <string>

namespace Euclid::Geometry {

// ======================================================
// Line Intersection Result
// ======================================================
template <typename T, int N>
struct LineIntersectionResult : public IntersectionResult {
    std::vector<Point<T, N>> points; // intersection points (if any)
    std::vector<Line<T, N>> lines;   // intersection lines (if any)

    LineIntersectionResult() = default;

    LineIntersectionResult(bool intersects, const std::vector<Point<T, N>>& pts = {},
                           const std::vector<Line<T, N>>& lns = {}, std::string desc = "")
        : IntersectionResult{intersects, desc}, points(pts), lines(lns) {}
};

// ======================================================
// Line–Line Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& l1, const Line<T, N>& l2,
                                      const Tolerance& tol = Tolerance()) {
    LineIntersectionResult<T, N> result;

    T scale = l1.direction().norm() + l2.direction().norm();
    if constexpr (N == 2) {
        // 2D intersection using determinant
        T dx1 = l1.direction()[0], dy1 = l1.direction()[1];
        T dx2 = l2.direction()[0], dy2 = l2.direction()[1];
        T det = dx1 * dy2 - dy1 * dx2;

        if (std::abs(det) < tol.evaluateEpsilon(scale)) {
            // Lines are parallel — check if they are coincident
            auto diff = l2.point1().coords - l1.point1().coords;
            T eps = tol.evaluateEpsilon(scale);

            bool coincident = false;
            if (diff.norm() < eps) {
                coincident = true;
            } else {
                T dot = diff.normalized().dot(l1.direction().normalized());
                coincident = std::abs(std::abs(dot) - 1.0) < eps;
            }

            if (coincident) {
                result.intersects = true;
                result.lines = {l1};
                result.description = "Lines are coincident";
            } else {
                result.description = "Lines are parallel";
            }
            return result;
        }
        auto diff = l2.point1().coords - l1.point1().coords;
        T t1 = (diff[0] * dy2 - diff[1] * dx2) / det;
        Point<T, N> pt(l1.point1().coords + t1 * l1.direction());
        result.intersects = true;
        result.points = {pt};
        result.description = "Lines intersect at a point";
        return result;
    }

    if constexpr (N == 3) {
        // 3D intersection using closest point method
        auto p1 = l1.point1().coords;
        auto d1 = l1.direction();
        auto p2 = l2.point1().coords;
        auto d2 = l2.direction();
        auto r = p1 - p2;
        T a = d1.dot(d1), b = d1.dot(d2), c = d2.dot(d2);
        T d = d1.dot(r), e = d2.dot(r);
        T denom = a * c - b * b;

        if (std::abs(denom) < tol.evaluateEpsilon(scale)) {
            // Lines are parallel — check if they are coincident
            auto diff = l2.point1().coords - l1.point1().coords;
            T eps = tol.evaluateEpsilon(scale);

            bool coincident = false;
            if (diff.norm() < eps) {
                coincident = true;
            } else {
                T dot = diff.normalized().dot(l1.direction().normalized());
                coincident = std::abs(std::abs(dot) - 1.0) < eps;
            }

            if (coincident) {
                result.intersects = true;
                result.lines = {l1};
                result.description = "Lines are coincident";
            } else {
                result.description = "Lines are parallel";
            }
            return result;
        }
        T t1 = (b * e - c * d) / denom;
        T t2 = (a * e - b * d) / denom;
        Point<T, N> pt1(p1 + t1 * d1);
        Point<T, N> pt2(p2 + t2 * d2);
        T dist = (pt1 - pt2).norm();
        if (dist > tol.evaluateEpsilon(dist)) {
            result.description = "Lines are skew (do not intersect)";
            return result;
        }
        result.intersects = true;
        result.points = {pt1};
        result.description = "Lines intersect at a point";
        return result;
    }

    // Generic ND implementation: check if lines are parallel or coincident
    auto p1 = l1.point1().coords;
    auto d1 = l1.direction().normalized();
    auto p2 = l2.point1().coords;
    auto d2 = l2.direction().normalized();
    auto diff = p2 - p1;

    T dot_dir = d1.dot(d2);
    T eps = tol.evaluateEpsilon(scale);

    if (std::abs(std::abs(dot_dir) - 1) < eps) {
        // Lines are parallel, check if coincident
        T dist = diff.norm();
        if (dist < eps) {
            // Lines are coincident
            result.intersects = true;
            result.lines = {l1};
            result.description = "Lines are coincident";
        } else {
            // Check if diff is along the direction line (collinear)
            T proj = diff.dot(d1);
            Point<T, N> projected_point = Point<T, N>(p1 + proj * d1);
            if ((projected_point.coords - p2).norm() < eps) {
                result.intersects = true;
                result.lines = {l1};
                result.description = "Lines are coincident";
            } else {
                result.description = "Lines are parallel";
            }
        }
        return result;
    }

    // For ND, solve for intersection parameters if possible
    // Using least squares to solve for parameters t1 and t2 where:
    // p1 + t1*d1 = p2 + t2*d2

    // Construct system: t1*d1 - t2*d2 = p2 - p1
    // Rewrite as A * [t1; t2] = b, where A = [d1, -d2], b = diff

    Eigen::Matrix<T, N, 2> A;
    for (int i = 0; i < N; ++i) {
        A(i, 0) = d1[i];
        A(i, 1) = -d2[i];
    }
    Eigen::Matrix<T, 2, 1> b;
    for (int i = 0; i < N; ++i) {
        b(i, 0) = diff[i];
    }

    // Solve least squares
    Eigen::Matrix<T, 2, 1> t = A.colPivHouseholderQr().solve(b);

    Point<T, N> pt1 = Point<T, N>(p1 + t(0, 0) * d1);
    Point<T, N> pt2 = Point<T, N>(p2 + t(1, 0) * d2);
    T dist = (pt1 - pt2).norm();

    if (dist > tol.evaluateEpsilon(dist)) {
        result.description = "Lines do not intersect";
        return result;
    }

    result.intersects = true;
    result.points = {pt1};
    result.description = "Lines intersect at a point";
    return result;
}

// ======================================================
// Line–Curve Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& line, const Curve<T, N>& c,
                                      const Tolerance& tol = Tolerance(), int max_samples = -1) {
    LineIntersectionResult<T, N> result;

    // domain of curve
    auto domain = c.domain();
    T t0 = domain.first.first;
    T t1 = domain.first.second;
    if (t1 <= t0) return result;

    // adaptive sample count based on curve extent and tolerance
    T extent = std::abs(t1 - t0);
    T eps = tol.evaluateEpsilon(extent);
    int samples = max_samples > 0 ? max_samples : std::max(20, static_cast<int>(std::min(200.0, std::max(20.0, (double)(extent / eps) * 5.0))));

    // ensure the three canonical samples exist: t0, mid, t1
    std::vector<T> sampled_ts;
    sampled_ts.reserve(samples);
    sampled_ts.push_back(t0);
    sampled_ts.push_back((t0 + t1) / 2);
    sampled_ts.push_back(t1);

    for (int i = 0; i < samples; ++i) {
        T uu = t0 + (t1 - t0) * static_cast<T>(i) / static_cast<T>(samples - 1);
        // avoid duplicating the three canonical samples
        if (std::abs(uu - t0) < eps || std::abs(uu - t1) < eps || std::abs(uu - (t0 + t1) / 2) < eps) continue;
        sampled_ts.push_back(uu);
    }

    // Closest projection onto line helper
    const auto& P0 = line.point1().coords;
    const auto dir = line.direction();

    T bestDist = std::numeric_limits<T>::infinity();
    T best_u = t0;
    T best_t = T(0);
    Point<T, N> best_pt = c.evaluate(t0);

    for (size_t i = 0; i < sampled_ts.size(); ++i) {
        T u = sampled_ts[i];
        Point<T, N> cpt = c.evaluate(u);
        // projector parameter t along line
        T proj_t = dir.dot(cpt.coords - P0) / dir.squaredNorm();
        Point<T, N> lpt(P0 + proj_t * dir);
        T dist = (cpt - lpt).norm();
        if (dist < bestDist) {
            bestDist = dist;
            best_u = u;
            best_t = proj_t;
            best_pt = cpt;
        }
    }

    // Newton-like refinement solving for (u, t) that minimize ||C(u) - L(t)||^2.
    // Use numerical derivative for C'(u).
    int max_iters = 20;
    T alpha = T(0.2);
    T prevDist = bestDist;
    T u = best_u;
    T t = best_t;

    for (int iter = 0; iter < max_iters; ++iter) {
        Point<T, N> Cu = c.evaluate(u);
        // numerical derivative C'(u)
        T h = std::max((T)1e-7, (t1 - t0) * (T)1e-6);
        T up = std::min(t1, u + h);
        T um = std::max(t0, u - h);
        Eigen::Matrix<T, N, 1> deriv = (c.evaluate(up).coords - c.evaluate(um).coords) / (up - um);

        // residual and jacobian
        Eigen::Matrix<T, N, 1> R;
        for (int i = 0; i < N; ++i) R(i, 0) = Cu.coords[i] - (P0[i] + t * dir[i]);

        // J = [C'(u), -dir]
        Eigen::Matrix<T, N, 2> J;
        for (int i = 0; i < N; ++i) {
            J(i, 0) = deriv[i];
            J(i, 1) = -dir[i];
        }

        Eigen::Matrix<T, 2, 2> JTJ = J.transpose() * J;
        Eigen::Matrix<T, 2, 1> JTR = J.transpose() * R;

        // solve JTJ * delta = -JTR with damping fallback
        Eigen::Matrix<T, 2, 1> delta;
        Eigen::FullPivLU<Eigen::Matrix<T, 2, 2>> lu(JTJ);
        if (lu.isInvertible()) {
            delta = lu.solve(-JTR);
        } else {
            // damped least squares fallback
            Eigen::Matrix<T, 2, 2> JTJ_damped = JTJ;
            JTJ_damped(0, 0) += eps;
            JTJ_damped(1, 1) += eps;
            delta = JTJ_damped.colPivHouseholderQr().solve(-JTR);
        }

        T du = delta(0, 0);
        T dt = delta(1, 0);

        u += alpha * du;
        t += alpha * dt;

        // clamp
        if (u < t0) u = t0;
        if (u > t1) u = t1;

        // recompute distance
        Point<T, N> Cnew = c.evaluate(u);
        Point<T, N> Lnew(P0 + t * dir);
        T dist = (Cnew - Lnew).norm();

        if (std::abs(dist - prevDist) < eps * 0.5) break;
        prevDist = dist;
    }

    // final check
    Point<T, N> Cfinal = c.evaluate(u);
    Point<T, N> Lfinal(P0 + t * dir);
    T finalDist = (Cfinal - Lfinal).norm();
    if (finalDist <= tol.evaluateEpsilon(std::max((T)1.0, finalDist))) {
        result.intersects = true;
        result.points = {Point<T, N>( (Cfinal.coords + Lfinal.coords) * (T)0.5 )};
        result.description = "Line intersects curve (within tolerance)";
    } else {
        result.description = "Line does not intersect curve";
    }

    return result;
}


// ======================================================
// Line–Surface Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& line, const Surface<T, N>& s,
                                      const Tolerance& tol = Tolerance(), int max_grid = -1) {
    LineIntersectionResult<T, N> result;
    // Only 3D surfaces supported
    if constexpr (N != 3) {
        result.description = "Line–Surface only implemented for 3D";
        return result;
    }

    // domain
    auto domain = s.surfaceDomain; // ((u0,u1),(v0,v1))
    T u0 = domain.first.first;
    T u1 = domain.first.second;
    T v0 = domain.second.first;
    T v1 = domain.second.second;
    if (u1 <= u0 || v1 <= v0) return result;

    T u_extent = std::abs(u1 - u0);
    T v_extent = std::abs(v1 - v0);
    T extent = std::max(u_extent, v_extent);
    T eps = tol.evaluateEpsilon(extent);

    int M = max_grid > 0 ? max_grid : std::min(200, std::max(20, static_cast<int>((int)((u_extent + v_extent) / std::max(eps, (T)1e-9)) * 5)));

    // sample key points: corners and center
    std::vector<std::pair<T, T>> samples;
    samples.reserve(M * M);
    samples.push_back({u0, v0});
    samples.push_back({u0, v1});
    samples.push_back({u1, v0});
    samples.push_back({u1, v1});
    samples.push_back({(u0 + u1) / 2, (v0 + v1) / 2});

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            T uu = u0 + (u1 - u0) * static_cast<T>(i) / static_cast<T>(M - 1);
            T vv = v0 + (v1 - v0) * static_cast<T>(j) / static_cast<T>(M - 1);
            // skip duplicates
            if ((std::abs(uu - u0) < eps && std::abs(vv - v0) < eps) ||
                (std::abs(uu - u0) < eps && std::abs(vv - v1) < eps) ||
                (std::abs(uu - u1) < eps && std::abs(vv - v0) < eps) ||
                (std::abs(uu - u1) < eps && std::abs(vv - v1) < eps) ||
                (std::abs(uu - (u0+u1)/2) < eps && std::abs(vv - (v0+v1)/2) < eps)) continue;
            samples.push_back({uu, vv});
        }
    }

    const auto& P0 = line.point1().coords;
    const auto dir = line.direction();

    T bestDist = std::numeric_limits<T>::infinity();
    T best_u = u0;
    T best_v = v0;
    Point<T, 3> best_spt = s.evaluate(u0, v0);

    // coarse sampling
    for (size_t i = 0; i < samples.size(); ++i) {
        T uu = samples[i].first;
        T vv = samples[i].second;
        Point<T, 3> spt = s.evaluate(uu, vv);
        // project onto line
        T proj_t = dir.dot(spt.coords - P0) / dir.squaredNorm();
        Point<T, 3> lpt(P0 + proj_t * dir);
        T dist = (spt - lpt).norm();
        if (dist < bestDist) {
            bestDist = dist;
            best_u = uu;
            best_v = vv;
            best_spt = spt;
        }
    }

    // refine with Newton for (u,v,t) solving S(u,v) - L(t) = 0
    int max_iters = 30;
    T alpha = T(1.0);
    T prevDist = bestDist;
    T u = best_u;
    T v = best_v;
    T t = dir.dot(best_spt.coords - P0) / dir.squaredNorm();

    for (int iter = 0; iter < max_iters; ++iter) {
        Point<T, 3> Spt = s.evaluate(u, v);

        // numerical partials
        T uh = std::max((T)1e-7, (u1 - u0) * (T)1e-6);
        T vh = std::max((T)1e-7, (v1 - v0) * (T)1e-6);
        Point<T, 3> Sup = s.evaluate(std::min(u1, u + uh), v);
        Point<T, 3> Sum = s.evaluate(std::max(u0, u - uh), v);
        Point<T, 3> Svp = s.evaluate(u, std::min(v1, v + vh));
        Point<T, 3> Svm = s.evaluate(u, std::max(v0, v - vh));

        Eigen::Matrix<T, N, 1> Su = (Sup.coords - Sum.coords) / (Sup.coords[0] - Sum.coords[0] + (T)1e-12); // fallback to component-wise safe
        Eigen::Matrix<T, N, 1> Sv = (Svp.coords - Svm.coords) / (Svp.coords[1] - Svm.coords[1] + (T)1e-12);
        // The above is a simple numerical approximation; component-wise safe division prevents zero denom.

        // residual
        Eigen::Matrix<T, 3, 1> R;
        for (int k = 0; k < 3; ++k) R(k, 0) = Spt.coords[k] - (P0[k] + t * dir[k]);

        // Jacobian J = [Su Sv -dir]
        Eigen::Matrix<T, 3, 3> J;
        for (int k = 0; k < 3; ++k) {
            J(k, 0) = Su[k];
            J(k, 1) = Sv[k];
            J(k, 2) = -dir[k];
        }

        // solve J * delta = -R
        Eigen::Matrix<T, 3, 1> delta;
        Eigen::FullPivLU<Eigen::Matrix<T, 3, 3>> lu(J);
        if (lu.isInvertible()) {
            delta = lu.solve(-R);
        } else {
            // damped / least squares fallback
            Eigen::Matrix<T, 3, 3> JTJ = J.transpose() * J;
            Eigen::Matrix<T, 3, 1> JTR = J.transpose() * R;
            for (int d = 0; d < 3; ++d) JTJ(d, d) += eps;
            delta = JTJ.colPivHouseholderQr().solve(-JTR);
        }

        T du = delta(0, 0);
        T dv = delta(1, 0);
        T dt = delta(2, 0);

        u += alpha * du;
        v += alpha * dv;
        t += alpha * dt;

        // clamp
        if (u < u0) u = u0; if (u > u1) u = u1;
        if (v < v0) v = v0; if (v > v1) v = v1;

        Point<T, 3> Snew = s.evaluate(u, v);
        Point<T, 3> Lnew(P0 + t * dir);
        T dist = (Snew - Lnew).norm();
        if (std::abs(dist - prevDist) < eps * 0.5) break;
        prevDist = dist;
    }

    Point<T, 3> Sfinal = s.evaluate(u, v);
    Point<T, 3> Lfinal(P0 + t * dir);
    T finalDist = (Sfinal - Lfinal).norm();
    if (finalDist <= tol.evaluateEpsilon(std::max((T)1.0, finalDist))) {
        result.intersects = true;
        result.points = {Point<T, 3>((Sfinal.coords + Lfinal.coords) * (T)0.5)};
        result.description = "Line intersects surface (within tolerance)";
    } else {
        result.description = "Line does not intersect surface";
    }

    return result;
}

// ======================================================
// Line–Plane Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& line, const Plane<T, N>& plane,
                                      const Tolerance& tol = Tolerance()) {
    LineIntersectionResult<T, N> result;
    if constexpr (N == 3) {
        auto n = plane.normal;
        auto diff = line.point1().coords - plane.base.coords;
        T denom = n.dot(line.direction());
        T scale = line.direction().norm() + plane.normal.norm();
        if (std::abs(denom) < tol.evaluateEpsilon(scale)) {
            result.description = "Line is parallel to plane";
            return result;
        }
        T t = -n.dot(diff) / denom;
        Point<T, N> pt(line.point1().coords + t * line.direction());
        result.intersects = true;
        result.points = {pt};
        result.description = "Line intersects plane at a point";
        return result;
    }
    result.description = "ND Line–Plane intersection not implemented";
    return result;
}

// Reverse overload: Plane–Line
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Plane<T, N>& plane, const Line<T, N>& line,
                                      const Tolerance& tol = Tolerance()) {
    return intersect(line, plane, tol);
}

// ======================================================
// Line–Segment Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& line, const Segment<T, N>& seg,
                                      const Tolerance& tol = Tolerance()) {
    LineIntersectionResult<T, N> result;

    const auto& P0 = line.point1().coords;
    const auto dir = line.direction();
    const auto& A = seg.start.coords;
    const auto& B = seg.end.coords;
    const auto segDir = B - A;

    // Solve for t (line) and u (segment) such that P0 + t*dir = A + u*segDir
    Eigen::Matrix<T, N, 2> A_mat;
    for (int i = 0; i < N; ++i) {
        A_mat(i, 0) = dir[i];
        A_mat(i, 1) = -segDir[i];
    }
    Eigen::Matrix<T, N, 1> b = A - P0;

    Eigen::Matrix<T, 2, 1> x = A_mat.colPivHouseholderQr().solve(b);
    T t = x(0, 0);
    T u = x(1, 0);

    // Check if intersection lies within segment
    if (u < -tol.evaluateEpsilon(1.0) || u > 1.0 + tol.evaluateEpsilon(1.0)) {
        result.description = "Line and segment do not intersect (outside segment bounds)";
        return result;
    }

    Point<T, N> P_int(P0 + t * dir);
    Point<T, N> S_int(A + u * segDir);

    if ((P_int - S_int).norm() <= tol.evaluateEpsilon((P_int - S_int).norm())) {
        result.intersects = true;
        result.points = {P_int};
        result.description = "Line intersects segment at a point";
    } else {
        result.description = "Line and segment are skew (no intersection)";
    }

    return result;
}

// ======================================================
// Line–Face Intersection
// ======================================================
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Line<T, N>& line, const Face<T, N>& face,
                                      const Tolerance& tol = Tolerance()) {
    LineIntersectionResult<T, N> result;

    // First, intersect with the face’s supporting plane
    Plane<T, N> facePlane(face.base, face.normal);
    auto planeRes = intersect(line, facePlane, tol.evaluateEpsilon(1.0));
    if (!planeRes.intersects || planeRes.points.empty()) {
        result.description = "Line does not intersect face plane";
        return result;
    }

    Point<T, N> candidate = planeRes.points[0];

    // Check if intersection point lies within the face polygon
    if (face.contains(candidate, tol.evaluateEpsilon(1.0))) {
        result.intersects = true;
        result.points = {candidate};
        result.description = "Line intersects face at a point";
    } else {
        result.description = "Line intersects face plane but outside bounds";
    }

    return result;
}


// Reverse overload: Curve–Line
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Curve<T, N>& c, const Line<T, N>& line,
                                       const Tolerance& tol = Tolerance(), int max_samples = -1) {
    return intersect(line, c, tol, max_samples);
}

// Reverse overload: Surface–Line
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Surface<T, N>& s, const Line<T, N>& line,
                                       const Tolerance& tol = Tolerance(), int max_grid = -1) {
    return intersect(line, s, tol, max_grid);
}

// Reverse overload: Segment–Line
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Segment<T, N>& seg, const Line<T, N>& line,
                                       const Tolerance& tol = Tolerance()) {
    return intersect(line, seg, tol);
}

// Reverse overload: Face–Line
template <typename T, int N>
LineIntersectionResult<T, N> intersect(const Face<T, N>& face, const Line<T, N>& line,
                                       const Tolerance& tol = Tolerance()) {
    return intersect(line, face, tol);
}

} // namespace Euclid::Geometry
