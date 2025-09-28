#pragma once
#include <array>
#include <functional>
#include <optional>
#include <cmath>
#include <numeric>
#include <iostream>
#include <Eigen/Dense>
#include "point.hpp"
#include "face.hpp"

namespace euclid::geometry {

template<typename T, int N>
class Surface {
public:
    using PointType = Point<T, N>;

    // Parametric function: S(u,v) -> ℝⁿ
    std::function<PointType(T, T)> surfaceFunc;

    // Domain: ((umin, vmin), (umax, vmax))
    std::pair<std::pair<T,T>, std::pair<T,T>> surfaceDomain;

    Surface(const std::function<PointType(T,T)>& f,
            const std::pair<std::pair<T,T>, std::pair<T,T>>& domain)
        : surfaceFunc(f), surfaceDomain(domain) {}

    // Evaluate at (u,v)
    PointType evaluate(T u, T v) const {
        return surfaceFunc(u,v);
    }

    // Apply a transformation to the surface
    template<typename Transform>
    Surface<T,N> applyTransform(const Transform& transform) const {
        auto newFunc = [transform, f = surfaceFunc](T uVal, T vVal) -> PointType {
            return transform.apply(f(uVal, vVal));
        };
        return Surface<T,N>(newFunc, surfaceDomain);
    }

    // Sweep a cross-section along a spine curve to create a surface
    template<typename Curve, typename CrossSectionFunc>
    static Surface<T, 3> sweepSurface(
        const Curve& spineCurve,
        const CrossSectionFunc& crossSectionFunc,
        std::pair<std::pair<T, T>, std::pair<T, T>> domain = {})
    {
        T umin = spineCurve.domain().first;
        T umax = spineCurve.domain().second;

        if (domain.first.first == T() && domain.first.second == T() &&
            domain.second.first == T() && domain.second.second == T()) {
            domain = {{umin, T(0)}, {umax, T(1)}};
        }

        // --- Rotation-Minimizing Frame (RMF) ---
        int samples = 500;
        std::vector<T> uSamples(samples);
        std::vector<Eigen::Matrix<T,3,1>> spinePoints(samples);
        for(int i=0;i<samples;++i){
            T uVal = umin + (umax-umin)*T(i)/(samples-1);
            uSamples[i] = uVal;
            spinePoints[i] = spineCurve.evaluate(uVal).coords;
        }

        std::vector<Eigen::Matrix<T,3,1>> tangents(samples);
        std::vector<Eigen::Matrix<T,3,1>> normals(samples);
        std::vector<Eigen::Matrix<T,3,1>> binormals(samples);

        // Compute tangents
        for(int i=0;i<samples;++i){
            if(i<samples-1)
                tangents[i] = (spinePoints[i+1]-spinePoints[i]).normalized();
            else
                tangents[i] = (spinePoints[i]-spinePoints[i-1]).normalized();
        }

        // Initialize first normal
        Eigen::Matrix<T,3,1> arbitrary = Eigen::Matrix<T,3,1>::UnitY();
        if(std::abs(tangents[0].dot(arbitrary)) > T(0.9)) arbitrary = Eigen::Matrix<T,3,1>::UnitX();
        normals[0] = (arbitrary - tangents[0]*(tangents[0].dot(arbitrary))).normalized();
        binormals[0] = tangents[0].cross(normals[0]).normalized();

        // RMF propagation
        for(int i=1;i<samples;++i){
            Eigen::Matrix<T,3,1> vPrev = tangents[i-1];
            Eigen::Matrix<T,3,1> vCurr = tangents[i];
            Eigen::Matrix<T,3,1> axis = vPrev.cross(vCurr);
            T axisNorm = axis.norm();
            if(axisNorm < 1e-12){
                normals[i] = normals[i-1];
                binormals[i] = binormals[i-1];
                continue;
            }
            axis /= axisNorm;
            T angle = std::acos(std::clamp(vPrev.dot(vCurr), T(-1), T(1)));
            auto rotate = [&](const Eigen::Matrix<T,3,1>& vec){
                Eigen::Matrix<T,3,1> rotated = vec*std::cos(angle) + axis.cross(vec)*std::sin(angle) + axis*(axis.dot(vec))*(1-std::cos(angle));
                return rotated.normalized();
            };
            normals[i] = rotate(normals[i-1]);
            binormals[i] = rotate(binormals[i-1]);
        }

        // Frame interpolation function
        auto frameAt = [=](T uVal){
            T tNorm = (uVal-umin)/(umax-umin);
            int idx = std::min(int(tNorm*(samples-1)), samples-1);
            return std::make_tuple(tangents[idx], normals[idx], binormals[idx]);
        };

        auto surfFunc = [spineCurve, crossSectionFunc, frameAt](T uVal, T vVal){
            auto [tangent, normal, binormal] = frameAt(uVal);
            Eigen::Matrix<T,3,1> center = spineCurve.evaluate(uVal).coords;

            Eigen::Matrix<T,3,1> cs;
            if constexpr (std::is_invocable_v<decltype(crossSectionFunc), T>) {
                cs = crossSectionFunc(vVal).coords; // old 1-arg API
            } else {
                cs = crossSectionFunc(vVal, Point<T,3>(tangent), Point<T,3>(normal), Point<T,3>(binormal)).coords;
            }

            Eigen::Matrix<T,3,1> offset = normal*cs(0) + binormal*cs(1);
            return Point<T,3>(center + offset);
        };

        return Surface<T,3>(surfFunc, domain);
    }
};

// Mesh representation
template<typename T, int N>
class SurfaceMesh {
public:
    using PointType = Point<T,N>;

    std::vector<PointType> vertices;
    std::vector<std::vector<size_t>> faces;

    SurfaceMesh() = default;

    bool validate() const {
        size_t numVertices = vertices.size();
        for (const auto& f : faces) {
            if (f.size() < 3) return false;
            for (size_t idx : f) {
                if (idx >= numVertices) return false;
            }
        }
        return true;
    }

    // Apply a transformation to the mesh
    template<typename Transform>
    SurfaceMesh<T,N> applyTransform(const Transform& transform) const {
        SurfaceMesh<T,N> meshT;
        meshT.vertices.reserve(vertices.size());
        for (auto& v : vertices)
            meshT.vertices.push_back(transform.apply(v));
        meshT.faces = faces;
        return meshT;
    }

    // OBJ export only supports 3D
    std::optional<std::string> exportOBJ() const {
        if constexpr (N != 3) {
            return std::nullopt;
        } else {
            std::string s;
            for (auto& v : vertices) {
                s += "v ";
                for (int i=0;i<3;++i) s += std::to_string(v.coords(i)) + " ";
                s += "\n";
            }
            for (auto& f : faces) {
                s += "f ";
                for (auto idx : f) s += std::to_string(idx+1) + " ";
                s += "\n";
            }
            return s;
        }
    }

    // Approximate area for 3D surfaces
    T area() const {
        if constexpr (N != 3) return T(0);
        T total = 0;
        for (auto& f : faces) {
            if (f.size() != 3) continue;
            auto& a = vertices[f[0]].coords;
            auto& b = vertices[f[1]].coords;
            auto& c = vertices[f[2]].coords;
            auto ab = b - a;
            auto ac = c - a;
            total += (ab.cross(ac)).norm() / 2;
        }
        return total;
    }
};

// Generate a grid-based mesh
template<typename T, int N>
SurfaceMesh<T,N> generateSurfaceMesh(const Surface<T,N>& surf, int uSteps, int vSteps) {
    SurfaceMesh<T,N> mesh;
    auto& domain = surf.surfaceDomain;
    T umin = domain.first.first;
    T vmin = domain.first.second;
    T umax = domain.second.first;
    T vmax = domain.second.second;
    T stepU = (umax-umin)/T(uSteps-1);
    T stepV = (vmax-vmin)/T(vSteps-1);

    std::vector<T> us(uSteps);
    std::vector<T> vs(vSteps);
    for(int i=0;i<uSteps;++i) us[i] = umin + stepU*i;
    for(int j=0;j<vSteps;++j) vs[j] = vmin + stepV*j;

    // Vertices
    for(int j=0;j<vSteps;++j){
        for(int i=0;i<uSteps;++i){
            mesh.vertices.push_back(surf.evaluate(us[i], vs[j]));
        }
    }

    // Faces (triangles from quads)
    for(int j=0;j<vSteps-1;++j){
        for(int i=0;i<uSteps-1;++i){
            size_t idx = i + j*uSteps;
            size_t idxRight = idx+1;
            size_t idxDown = idx+uSteps;
            size_t idxDiag = idx+uSteps+1;
            mesh.faces.push_back({idx, idxRight, idxDiag});
            mesh.faces.push_back({idx, idxDiag, idxDown});
        }
    }
    return mesh;
}

} // namespace euclid::geometry
