#pragma once
#include <array>
#include <functional>
#include <optional>
#include <cmath>
#include <numeric>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <Eigen/Dense>
#include "point.hpp"

namespace Euclid::Geometry {

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
            Euclid::Tolerance tol;
            if(Euclid::equalWithinTolerance(axisNorm, T(0), tol, axisNorm)){
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

            Point<T,2> cs = crossSectionFunc(vVal);

            Eigen::Matrix<T,3,1> offset = normal*cs.coords(0) + binormal*cs.coords(1);
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

    struct TopologyCheckResult {
        bool isClosed;
        std::vector<std::vector<size_t>> boundarySubfaces;
    };

    TopologyCheckResult hasClosedTopology() const {
        TopologyCheckResult result{true, {}};
        if constexpr (N < 2) {
            // Not defined for dimension less than 2 (no faces)
            result.isClosed = false;
            return result;
        }

        // Determine subface size: (N-2)-faces have (N-1) vertices
        const int subfaceSize = N - 1;

        // Map from canonical subface (sorted indices) to count
        std::map<std::vector<size_t>, int> subfaceCount;

        // Helper lambda to check if face is degenerate (repeated vertices)
        auto isDegenerateFace = [&](const std::vector<size_t>& face) -> bool {
            std::set<size_t> uniqueVerts(face.begin(), face.end());
            if (uniqueVerts.size() < face.size()) return true;
            if constexpr (N == 3) {
                // For 3D triangles, also check zero area
                if (face.size() == 3) {
                    const auto& a = vertices[face[0]].coords;
                    const auto& b = vertices[face[1]].coords;
                    const auto& c = vertices[face[2]].coords;
                    auto ab = b - a;
                    auto ac = c - a;
                    auto cross = ab.cross(ac);
                    if (cross.norm() < 1e-12) return true;
                }
            }
            return false;
        };

        for (const auto& face : faces) {
            if (face.size() < subfaceSize) {
                // Skip degenerate or too small faces
                continue;
            }
            if (isDegenerateFace(face)) continue;

            if constexpr (N == 3) {
                // For 3D triangular mesh, (N-2)-faces are edges (pairs of vertices)
                // Each face has 3 edges
                for (size_t i = 0; i < 3; ++i) {
                    std::vector<size_t> edge = {face[i], face[(i+1)%3]};
                    std::sort(edge.begin(), edge.end());
                    subfaceCount[edge]++;
                }
            } else {
                // For general ND, generate all (N-2)-faces by removing one vertex from the face
                // Each (N-2)-face has subfaceSize = N-1 vertices
                for (size_t i = 0; i < face.size(); ++i) {
                    std::vector<size_t> subface;
                    subface.reserve(subfaceSize);
                    for (size_t j = 0; j < face.size(); ++j) {
                        if (j != i) {
                            subface.push_back(face[j]);
                        }
                    }
                    std::sort(subface.begin(), subface.end());
                    subfaceCount[subface]++;
                }
            }
        }

        // Identify boundary subfaces (those with count != 2)
        for (const auto& [subface, count] : subfaceCount) {
            if (count != 2) {
                result.isClosed = false;
                result.boundarySubfaces.push_back(subface);
            }
        }

        return result;
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
            mesh.faces.push_back({idx, idxDiag, idxRight});
            mesh.faces.push_back({idx, idxDown, idxDiag});
        }
    }

    // (Planar capping and boundary triangulation removed)

    return mesh;
}

template<typename T>
SurfaceMesh<T,3> generatePeriodicWrappedMesh(
    const Surface<T,3>& surf,
    int uSteps,
    int vSteps,
    bool uPeriodic = false,
    bool vPeriodic = false
) {
    SurfaceMesh<T,3> mesh;
    auto& domain = surf.surfaceDomain;
    T umin = domain.first.first;
    T vmin = domain.first.second;
    T umax = domain.second.first;
    T vmax = domain.second.second;

    T stepU = (umax - umin) / T(uSteps - 1);
    T stepV = (vmax - vmin) / T(vSteps - 1);

    // Determine effective number of vertices in u and v directions
    int vertCountU = uPeriodic ? uSteps - 1 : uSteps;
    int vertCountV = vPeriodic ? vSteps - 1 : vSteps;

    // Generate vertices with seam vertices shared if periodic
    for (int j = 0; j < vertCountV; ++j) {
        T v = vmin + stepV * j;
        for (int i = 0; i < vertCountU; ++i) {
            T u = umin + stepU * i;
            mesh.vertices.push_back(surf.evaluate(u, v));
        }
    }

    // Helper to get vertex index with wrapping
    auto vertexIndex = [&](int i, int j) -> size_t {
        if (uPeriodic) {
            i = i % vertCountU;
            if (i < 0) i += vertCountU;
        }
        if (vPeriodic) {
            j = j % vertCountV;
            if (j < 0) j += vertCountV;
        }
        return i + j * vertCountU;
    };

    // Generate faces
    int faceCountU = uPeriodic ? vertCountU : vertCountU - 1;
    int faceCountV = vPeriodic ? vertCountV : vertCountV - 1;

    for (int j = 0; j < faceCountV; ++j) {
        for (int i = 0; i < faceCountU; ++i) {
            size_t idx00 = vertexIndex(i, j);
            size_t idx10 = vertexIndex(i + 1, j);
            size_t idx01 = vertexIndex(i, j + 1);
            size_t idx11 = vertexIndex(i + 1, j + 1);

            mesh.faces.push_back({idx00, idx10, idx11});
            mesh.faces.push_back({idx00, idx11, idx01});
        }
    }

    return mesh;
}

} // namespace Euclid::Geometry
