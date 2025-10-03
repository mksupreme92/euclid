#pragma once

#include "geometry/curve.hpp"
#include "geometry/volume.hpp"
#include "tests/test_utilities.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>

using namespace Euclid::Geometry;
using namespace Euclid::Algebra;
using namespace Euclid::Tests;

inline void testVolumeMeshing() {
    std::cout << "\n∰ Testing Volume Meshing\n";
    
    using Point3 = Point<float,3>;
    using Volume3 = Volume<float,3>;
    using Curve3 = Curve<float,3>;
    
    // --- Unit cube ---
    auto cubeFunc = [](float u, float v, float w) -> Point3 {
        return Point3{u,v,w};
    };
    Volume3 cubeVolume(cubeFunc, {{0,0,0},{1,1,1}});
    auto cubeMesh = generateStructuredVolumeMesh(cubeVolume, 2,2,2);
    
    printTest("CubeMesh: correct number of vertices", cubeMesh.vertices.size() == 8);
    printTest("CubeMesh: correct number of tets", cubeMesh.tets.size() == 5);
    printTest("CubeMesh: all tets valid", cubeMesh.validate());
    printTest("CubeMesh: total volume == 1.0", std::abs(cubeMesh.totalVolume() - 1.0f) < 1e-6f);
    cubeMesh.exportMesh("output/cube_volume.msh");
    
    // --- Parabolic volume ---
    auto paraFunc = [](float u, float v, float w) -> Point3 {
        return Point3{u, v, u*v*w};
    };
    Volume3 paraVolume(paraFunc, {{0,0,0},{1,1,1}});
    auto paraMesh = generateStructuredVolumeMesh(paraVolume, 3,3,3);
    
    printTest("ParaMesh: correct number of vertices", paraMesh.vertices.size() == 27);
    printTest("ParaMesh: correct number of tets", paraMesh.tets.size() == 5*8);
    printTest("ParaMesh: all tets valid", paraMesh.validate());
    printTest("ParaMesh: total volume > 0", paraMesh.totalVolume() > 0.0f);
    paraMesh.exportMesh("output/parabolic_volume.msh");
    
    // --- Swept volume (linear tube) ---
    float r = 0.1f;
    // Use new sweepVolume API: lambda returns absolute 3D coordinates directly
    auto sweptVol = Volume3::sweepVolume([r](float u, float v, float w) -> Point3 {
        float z = u;
        float angle = 2.0f * M_PI * w;
        float radius = r * v;
        float x = radius * std::cos(angle);
        float y = radius * std::sin(angle);
        return Point3{x, y, z};
    });
    auto sweptMesh = generateStructuredVolumeMesh(sweptVol, 20, 10, 16); // uSteps, vSteps, wSteps
    
    printTest("SweepMesh: has vertices", sweptMesh.vertices.size() > 0);
    printTest("SweepMesh: has tets", sweptMesh.tets.size() > 0);
    printTest("SweepMesh: all tets valid", sweptMesh.validate());
    
    sweptMesh.exportMesh("output/swept_volume.msh");
    
    // --- Trefoil volume ---
    auto trefoilCurveFunc = [](float t) -> Point3 {
        float x = std::sin(t) + 2.0f * std::sin(2.0f * t);
        float y = std::cos(t) - 2.0f * std::cos(2.0f * t);
        float z = -std::sin(3.0f * t);
        return Point3{x, y, z};
    };
    Curve3 trefoilCurve(trefoilCurveFunc, 0.0f, 2.0f * M_PI);
    
    
    float trefoilRadius = 0.5f;
    int uSteps = 120;
    int vSteps = 10;
    int wSteps = 16;
    
    auto trefoilVol = Volume3::sweepVolume(trefoilCurve, [trefoilRadius](float v, float w) -> Point<float,2> {
        float angle = 2.0f * M_PI * w;
        float radius = trefoilRadius * v;
        return Point<float,2>{radius * std::cos(angle), radius * std::sin(angle)};
    });
    auto trefoilMesh = generateStructuredVolumeMesh(trefoilVol, uSteps, vSteps, wSteps);
    
    printTest("TrefoilMesh: has vertices", trefoilMesh.vertices.size() > 0);
    printTest("TrefoilMesh: has tets", trefoilMesh.tets.size() > 0);
    printTest("TrefoilMesh: all tets valid", trefoilMesh.validate());
    
    trefoilMesh.exportMesh("output/trefoil_volume.msh");
    
}

inline void testVolume() {
    std::cout << "\n∰ Testing Volume Primitive\n";
    testVolumeMeshing();
}
