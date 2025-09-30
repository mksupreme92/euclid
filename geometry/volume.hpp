#pragma once
#include <vector>
#include <array>
#include <functional>
#include <optional>
#include <Eigen/Dense>
#include "point.hpp"
#include <fstream>
#include <string>

namespace euclid::geometry {

// ---------------------------
// Parametric Volume
// ---------------------------
template<typename T, int N = Eigen::Dynamic>
class Volume {
public:
    using PointType = Point<T, N>;

    std::function<PointType(T,T,T)> volumeFunc;
    std::pair<std::tuple<T,T,T>, std::tuple<T,T,T>> volumeDomain;

    Volume(const std::function<PointType(T,T,T)>& f,
           const std::pair<std::tuple<T,T,T>, std::tuple<T,T,T>>& domain)
        : volumeFunc(f), volumeDomain(domain) {}

    // ---------------------------
    // Overload: Construct volume directly from a lambda (no RMF/spine/cross-section)
    // ---------------------------
    template<typename Lambda>
    static Volume<T,3> sweepVolume(
        const Lambda& lambda,
        std::pair<std::tuple<T,T,T>, std::tuple<T,T,T>> domain = { {T(0),T(0),T(0)}, {T(1),T(1),T(1)} }
    ) {
        // Accepts lambda of form Point<T,3>(T u, T v, T w)
        return Volume<T,3>(lambda, domain);
    }

    // ---------------------------
    // Sweep a 2D cross-section along a spine curve to create a tube
    // ---------------------------
    template<typename Curve, typename CrossSectionFunc>
    static Volume<T,3> sweepVolume(
        const Curve& spineCurve,
        const CrossSectionFunc& crossSectionFunc,
        std::pair<std::pair<T,T>, std::pair<T,T>> domain = {})
    {
        // Set default domain if not provided
        if(domain == std::pair<std::pair<T,T>, std::pair<T,T>>{}) {
            domain = {{T(0),T(1)}, {T(0),T(1)}};
        }
        auto [vRange, wRange] = domain;
        auto [vmin, vmax] = vRange;
        auto [wmin, wmax] = wRange;

        // Extract spine curve domain (assumes .domain() returns std::pair<T,T>)
        auto [curveUmin, curveUmax] = spineCurve.domain();
        T umin = curveUmin;
        T umax = curveUmax;

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

        auto frameAt = [=](T uVal){
            T tNorm = (uVal-umin)/(umax-umin);
            int idx = std::min(int(tNorm*(samples-1)), samples-1);
            return std::make_tuple(tangents[idx], normals[idx], binormals[idx]);
        };

        // Parametric volume function
        auto volFunc = [spineCurve, crossSectionFunc, frameAt, vmin, vmax, wmin, wmax](T u, T v, T w){
            auto [tangent, normal, binormal] = frameAt(u);
            Eigen::Matrix<T,3,1> center = spineCurve.evaluate(u).coords;

            // crossSectionFunc returns Point<T,2> with (v,w) in [vmin,vmax]x[wmin,wmax]
            // Map v,w from [vmin,vmax] and [wmin,wmax] to [0,1] for crossSectionFunc input
            T vNorm = (v - vmin) / (vmax - vmin);
            T wNorm = (w - wmin) / (wmax - wmin);
            Point<T,2> cs = crossSectionFunc(vNorm, wNorm);
            // Map 2D cross-section into 3D using RMF
            Eigen::Matrix<T,3,1> offset = normal*cs.coords(0) + binormal*cs.coords(1);
            return Point<T,3>(center + offset);
        };

        return Volume<T,3>(volFunc, {std::make_tuple(umin,vmin,wmin), std::make_tuple(umax,vmax,wmax)});
    }

    // Evaluate parametric volume at (u,v,w)
    PointType evaluate(T u, T v, T w) const {
        return volumeFunc(u,v,w);
    }

    template<typename Transform>
    Volume<T,N> applyTransform(const Transform& Xform) const {
        auto newFunc = [Xform, f = volumeFunc](T u, T v, T w) -> PointType {
            return Xform.apply(f(u,v,w));
        };
        return Volume<T,N>(newFunc, volumeDomain);
    }
};

// ---------------------------
// VolumeMesh and mesh generation remain unchanged
// ---------------------------
template<typename T, int N = Eigen::Dynamic>
class VolumeMesh {
public:
    using PointType = Point<T,N>;

    std::vector<PointType> vertices;
    std::vector<std::array<size_t,4>> tets;

    VolumeMesh() = default;

    // Validate tetrahedra
    bool validate() const {
        size_t n = vertices.size();
        for (const auto& tet : tets) {
            for (size_t i=0;i<4;++i) {
                if (tet[i] >= n) return false;
                for (size_t j=0;j<i;++j) if (tet[i]==tet[j]) return false;
            }
        }
        return true;
    }

    // Apply transform to mesh
    template<typename Transform>
    VolumeMesh<T,N> applyTransform(const Transform& Xform) const {
        VolumeMesh<T,N> meshT;
        meshT.vertices.reserve(vertices.size());
        for (const auto& v : vertices)
            meshT.vertices.push_back(Xform.apply(v));
        meshT.tets = tets;
        return meshT;
    }

    // Compute volume of a single tetrahedron (only works for N>=3)
    T tetVolume(const std::array<size_t,4>& tet) const {
        if constexpr (N != 3) return T(0);
        const auto& a = vertices[tet[0]].coords;
        const auto& b = vertices[tet[1]].coords;
        const auto& c = vertices[tet[2]].coords;
        const auto& d = vertices[tet[3]].coords;
        Eigen::Matrix<T,3,3> mat;
        mat.col(0) = b - a;
        mat.col(1) = c - a;
        mat.col(2) = d - a;
        return std::abs(mat.determinant())/T(6);
    }

    // Total volume
    T totalVolume() const {
        if constexpr (N != 3) return T(0);
        T vol = T(0);
        for (const auto& tet : tets)
            vol += tetVolume(tet);
        return vol;
    }

    // Export mesh to file in Gmsh .msh 2.2 ASCII format
    void exportMesh(const std::string& filename) const {
        std::ofstream ofs(filename);
        if(!ofs) return;

        // Write MeshFormat section
        ofs << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

        // Write Nodes section
        ofs << "$Nodes\n" << vertices.size() << "\n";
        for(size_t i = 0; i < vertices.size(); ++i) {
            const auto& v = vertices[i];
            ofs << (i+1) << " " << v.coords(0) << " " << v.coords(1) << " " << v.coords(2) << "\n";
        }
        ofs << "$EndNodes\n";

        // Write Elements section
        ofs << "$Elements\n" << tets.size() << "\n";
        for(size_t i = 0; i < tets.size(); ++i) {
            const auto& tet = tets[i];
            // elementType = 4 (TET4), numTags = 0
            ofs << (i+1) << " 4 0 " << (tet[0]+1) << " " << (tet[1]+1) << " " << (tet[2]+1) << " " << (tet[3]+1) << "\n";
        }
        ofs << "$EndElements\n";
    }
};

// ---------------------------
// Structured tetrahedral mesh from parametric volume
// ---------------------------
template<typename T, int N = Eigen::Dynamic>
VolumeMesh<T,N> generateStructuredVolumeMesh(const Volume<T,N>& vol,
                                             int uSteps, int vSteps, int wSteps,
                                             [[maybe_unused]] const std::optional<std::string>& exportToMeshFile = std::nullopt) {
    VolumeMesh<T,N> mesh;

    auto [umin_vmin_wmin, umax_vmax_wmax] = vol.volumeDomain;
    auto [umin,vmin,wmin] = umin_vmin_wmin;
    auto [umax,vmax,wmax] = umax_vmax_wmax;

    // Step sizes
    T stepU = (umax-umin)/T(uSteps-1);
    T stepV = (vmax-vmin)/T(vSteps-1);
    T stepW = (wmax-wmin)/T(wSteps-1);

    // Generate vertices
    mesh.vertices.reserve(uSteps*vSteps*wSteps);
    for(int k=0;k<wSteps;++k){
        T w = wmin + k*stepW;
        for(int j=0;j<vSteps;++j){
            T v = vmin + j*stepV;
            for(int i=0;i<uSteps;++i){
                T u = umin + i*stepU;
                mesh.vertices.push_back(vol.evaluate(u,v,w));
            }
        }
    }

    // Helper: linear index
    auto idx = [uSteps,vSteps](int i,int j,int k){ return i + j*uSteps + k*uSteps*vSteps; };

    // Generate tets: split each cube into 5 tets
    for(int k=0;k<wSteps-1;++k){
        for(int j=0;j<vSteps-1;++j){
            for(int i=0;i<uSteps-1;++i){
                size_t v000 = idx(i,j,k);
                size_t v100 = idx(i+1,j,k);
                size_t v010 = idx(i,j+1,k);
                size_t v110 = idx(i+1,j+1,k);
                size_t v001 = idx(i,j,k+1);
                size_t v101 = idx(i+1,j,k+1);
                size_t v011 = idx(i,j+1,k+1);
                size_t v111 = idx(i+1,j+1,k+1);

                mesh.tets.push_back({v000,v100,v010,v001});
                mesh.tets.push_back({v100,v110,v010,v111});
                mesh.tets.push_back({v100,v010,v001,v111});
                mesh.tets.push_back({v010,v001,v011,v111});
                mesh.tets.push_back({v100,v001,v101,v111});
            }
        }
    }


    return mesh;
}

} // namespace euclid::geometry
