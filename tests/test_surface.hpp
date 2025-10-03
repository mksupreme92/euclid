#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "geometry/surface.hpp"
#include "geometry/point.hpp"
#include "test_utilities.hpp"

using namespace Euclid::Geometry;
using namespace Euclid::Algebra;
using namespace Euclid::Tests;

inline void testSurfaceMeshing() {
    std::cout << "\n∯ Testing Surface Meshing\n";
    
    using Point3 = Point<float,3>;
    using Surface3 = Surface<float,3>;
    using SurfaceMesh3 = SurfaceMesh<float,3>;
    
    // --- Unit square surface ---
    auto squareFunc = [](float u, float v) -> Point3 {
        return Point3{u, v, 0.0f};
    };
    Surface3 squareSurface(squareFunc, {{0.0f,0.0f},{1.0f,1.0f}});
    SurfaceMesh3 squareMesh = generateSurfaceMesh(squareSurface, 2, 2);
    
    printTest("SquareMesh: correct number of vertices", squareMesh.vertices.size() == static_cast<size_t>(4));
    printTest("SquareMesh: correct number of faces", squareMesh.faces.size() == static_cast<size_t>(2));
    printTest("SquareMesh: valid faces", squareMesh.validate());
    printTest("SquareMesh: area == 1.0", std::abs(squareMesh.area() - 1.0f) < 1e-6f);

    std::filesystem::create_directories("output");
    std::string squarePath = "output/square.obj";
    std::ofstream squareFile(squarePath);
    if (squareFile) {
        if (auto objText = squareMesh.exportOBJ()) {
            squareFile << *objText;
            std::cout << "✅ Exported square mesh to " << squarePath << "\n";
        } else {
            std::cerr << "❌ Failed to export square mesh: vertices not 3D\n";
        }
        squareFile.close();
    } else {
        std::cerr << "❌ Failed to open file for writing: " << squarePath << "\n";
    }
    
    // --- Paraboloid surface ---
    auto paraFunc = [](float u, float v) -> Point3 {
        return Point3{u, v, u*u + v*v};
    };
    Surface3 paraSurface(paraFunc, {{0.0f,0.0f},{1.0f,1.0f}});
    SurfaceMesh3 paraMesh = generateSurfaceMesh(paraSurface, 5, 5);
    
    printTest("ParaboloidMesh: correct number of vertices", paraMesh.vertices.size() == static_cast<size_t>(25));
    printTest("ParaboloidMesh: correct number of faces", paraMesh.faces.size() == static_cast<size_t>(32));
    printTest("ParaboloidMesh: valid faces", paraMesh.validate());
    printTest("ParaboloidMesh: area > 0", paraMesh.area() > 0.0f);

    std::string paraPath = "output/paraboloid.obj";
    std::ofstream paraFile(paraPath);
    if (paraFile) {
        if (auto objText = paraMesh.exportOBJ()) {
            paraFile << *objText;
            std::cout << "✅ Exported paraboloid mesh to " << paraPath << "\n";
        } else {
            std::cerr << "❌ Failed to export paraboloid mesh: vertices not 3D\n";
        }
        paraFile.close();
    } else {
        std::cerr << "❌ Failed to open file for writing: " << paraPath << "\n";
    }


    // --- Sweep surface along a trefoil curve ---
    using Curve3 = Curve<float,3>;
    // Trefoil knot parametric equations (scaled)
    auto trefoilFunc = [](float t) -> Point3 {
        float a = 0.5f; // scale
        float x = a * (std::sin(t) + 2.0f * std::sin(2.0f * t));
        float y = a * (std::cos(t) - 2.0f * std::cos(2.0f * t));
        float z = a * (-std::sin(3.0f * t));
        return Point3{x, y, z};
    };
    Curve3 trefoilCurve(trefoilFunc, 0.0f, 2.0f * 3.14159265f);

    int trefoilSteps = 100;
    int trefoilCrossSectionSteps = 12;
    float trefoilRadius = 0.07f;

    // Cross-section function: circle in cross-section plane, returns 2D coordinates
    auto trefoilCrossSection = [trefoilRadius](float v) -> Point<float,2> {
        float angle = v * 2.0f * 3.14159265f;
        float cx = trefoilRadius * std::cos(angle);
        float cy = trefoilRadius * std::sin(angle);
        return Point<float,2>{cx, cy};
    };

    // Sweep along trefoil
    Surface3 trefoilSweepSurface = Surface3::sweepSurface(
        trefoilCurve,
        trefoilCrossSection,
        {{0.0f,0.0f},{2.0f * 3.14159265f, 1.0f}}
    );
    SurfaceMesh3 trefoilSweepMesh = generateSurfaceMesh(
        trefoilSweepSurface, trefoilSteps, trefoilCrossSectionSteps
    );

    printTest("TrefoilSweepMesh: correct number of vertices",
        trefoilSweepMesh.vertices.size() == static_cast<size_t>(trefoilSteps * trefoilCrossSectionSteps));
    printTest("TrefoilSweepMesh: correct number of faces",
        trefoilSweepMesh.faces.size() == static_cast<size_t>((trefoilSteps-1)*(trefoilCrossSectionSteps-1)*2));
    printTest("TrefoilSweepMesh: valid faces", trefoilSweepMesh.validate());
    printTest("TrefoilSweepMesh: area > 0", trefoilSweepMesh.area() > 0.0f);

    std::string trefoilSweepPath = "output/trefoil_sweep.obj";
    std::ofstream trefoilSweepFile(trefoilSweepPath);
    if (trefoilSweepFile) {
        if (auto objText = trefoilSweepMesh.exportOBJ()) {
            trefoilSweepFile << *objText;
            std::cout << "✅ Exported trefoil sweep mesh to " << trefoilSweepPath << "\n";
        } else {
            std::cerr << "❌ Failed to export trefoil sweep mesh: vertices not 3D\n";
        }
        trefoilSweepFile.close();
    } else {
        std::cerr << "❌ Failed to open file for writing: " << trefoilSweepPath << "\n";
    }

    // --- Sweep surface along a curve ---
    using Curve3 = Curve<float,3>;

    Point3 p0{0.0f, 0.0f, 0.0f};
    Point3 p1{1.0f, 0.0f, 0.0f};
    Curve3 linearCurve = Curve3::linearCurve(p0, p1);

    int crossSectionSteps = 10;
    float radius = 0.1f;

    auto sweepFunc = [linearCurve, radius](float u, float v) -> Point3 {
        Point3 center = linearCurve.evaluate(u);
        float angle = v * 2.0f * 3.14159265f;
        float y = radius * std::cos(angle);
        float z = radius * std::sin(angle);
        return Point3{center.coords(0), center.coords(1)+y, center.coords(2)+z};
    };

    Surface3 sweepSurface(sweepFunc, {{0.0f,0.0f},{1.0f,1.0f}});
    SurfaceMesh3 sweepMesh = generateSurfaceMesh(sweepSurface, 20, crossSectionSteps);

    printTest("SweepMesh: correct number of vertices", sweepMesh.vertices.size() == static_cast<size_t>(20*crossSectionSteps));
    printTest("SweepMesh: correct number of faces", sweepMesh.faces.size() == static_cast<size_t>((20-1)*(crossSectionSteps-1)*2));
    printTest("SweepMesh: valid faces", sweepMesh.validate());
    printTest("SweepMesh: area > 0", sweepMesh.area() > 0.0f);

    std::string sweepPath = "output/sweep.obj";
    std::ofstream sweepFile(sweepPath);
    if (sweepFile) {
        if (auto objText = sweepMesh.exportOBJ()) {
            sweepFile << *objText;
            std::cout << "✅ Exported sweep mesh to " << sweepPath << "\n";
        } else {
            std::cerr << "❌ Failed to export sweep mesh: vertices not 3D\n";
        }
        sweepFile.close();
    } else {
        std::cerr << "❌ Failed to open file for writing: " << sweepPath << "\n";
    }


    // --- Saddle surface ---
    auto saddleFunc = [](float u, float v) -> Point3 {
        return Point3{u, v, u*u - v*v};
    };
    Surface3 saddleSurface(saddleFunc, {{-1.0f,-1.0f},{1.0f,1.0f}});
    SurfaceMesh3 saddleMesh = generateSurfaceMesh(saddleSurface, 20, 20);

    printTest("SaddleMesh: correct number of vertices", saddleMesh.vertices.size() == static_cast<size_t>(400));
    printTest("SaddleMesh: correct number of faces", saddleMesh.faces.size() == static_cast<size_t>(722)); // 19x19 quads * 2 triangles
    printTest("SaddleMesh: valid faces", saddleMesh.validate());
    printTest("SaddleMesh: area > 0", saddleMesh.area() > 0.0f);

    std::string saddlePath = "output/saddle.obj";
    std::ofstream saddleFile(saddlePath);
    if (saddleFile) {
        if (auto objText = saddleMesh.exportOBJ()) {
            saddleFile << *objText;
            std::cout << "✅ Exported saddle mesh to " << saddlePath << "\n";
        } else {
            std::cerr << "❌ Failed to export saddle mesh: vertices not 3D\n";
        }
        saddleFile.close();
    } else {
        std::cerr << "❌ Failed to open file for writing: " << saddlePath << "\n";
    }


}

inline void testSurface() {

    std::cout << "\n∯ Testing Surface Primitive\n";

    // --- 1. Define a simple 3D parametric surface (plane z = 0) ---
    using Point3 = Point<float,3>;
    using Surface3 = Surface<float,3>;

    auto planeFunc = [](float u, float v) -> Point3 {
        return Point3{u, v, 0.0f};
    };

    Surface3 planeSurface{planeFunc, {{0.0f,0.0f},{1.0f,1.0f}}};

    auto mesh = generateSurfaceMesh(planeSurface, 2, 2);

    // --- Validate vertices ---
    printTest("SurfaceMesh: correct number of vertices", mesh.vertices.size() == 4);

    // --- Validate faces ---
    printTest("SurfaceMesh: correct number of faces", mesh.faces.size() == 2);

    // --- Optional: Check area of unit plane ---
    float area = mesh.area();
    printTest("SurfaceMesh: area == 1.0", std::abs(area - 1.0f) < 1e-6f);

    // --- Check closed topology for plane mesh ---
    printTest("SurfaceMesh: hasClosedTopology() == false", !mesh.hasClosedTopology().isClosed);

    // --- 2. Define a curved surface (paraboloid z = u^2 + v^2) ---
    auto paraboloidFunc = [](float u, float v) -> Point3 {
        return Point3{u, v, u*u + v*v};
    };

    Surface3 paraboloidSurface{paraboloidFunc, {{-1.0f,-1.0f},{1.0f,1.0f}}};

    // --- Generate mesh with 3x3 vertices ---
    auto paraboloidMesh = generateSurfaceMesh(paraboloidSurface, 3, 3);

    // --- Validate vertices and faces counts ---
    printTest("ParaboloidMesh: correct number of vertices", paraboloidMesh.vertices.size() == 9);
    printTest("ParaboloidMesh: correct number of faces", paraboloidMesh.faces.size() == 8);

    // --- Optional: Check approximate area (should be > 0) ---
    float paraboloidArea = paraboloidMesh.area();
    printTest("ParaboloidMesh: area > 0", paraboloidArea > 0.0f);

    // --- Check closed topology for paraboloid mesh ---
    printTest("ParaboloidMesh: hasClosedTopology() == false", !paraboloidMesh.hasClosedTopology().isClosed);

    // --- Create a small manually created open mesh (missing one face) ---
    {
        SurfaceMesh<float,3> openMesh;
        // 4 vertices forming a square
        openMesh.vertices = {
            Point3{0.0f, 0.0f, 0.0f},
            Point3{1.0f, 0.0f, 0.0f},
            Point3{1.0f, 1.0f, 0.0f},
            Point3{0.0f, 1.0f, 0.0f}
        };
        // Only one face (triangle), missing the other to form closed square
        openMesh.faces = {
            {0, 1, 2}
        };
        printTest("OpenMesh: hasClosedTopology() == false", !openMesh.hasClosedTopology().isClosed);
    }

    // --- 3. Test affine transformation (scale + translate) ---
    {
        
        Eigen::Matrix3f linear;
        linear << 2.0f, 0.0f, 0.0f,
                  0.0f, 3.0f, 0.0f,
                  0.0f, 0.0f, 1.0f;
        Eigen::Vector3f translation;
        translation << 1.0f, -1.0f, 0.5f;
        Affine<float,3> transform(linear, translation);

        // --- Apply transform to planeSurface ---
        Surface3 transformedPlane{
            [f=planeSurface.surfaceFunc, &transform](float u, float v) -> Point3 {
                return transform.apply(f(u,v));
            },
            planeSurface.surfaceDomain
        };

        auto transformedMesh = generateSurfaceMesh(transformedPlane, 2, 2);

        // --- Validate transformed vertices manually ---
        bool verticesCorrect = true;
        auto origVerts = mesh.vertices;
        for(size_t i=0;i<origVerts.size();++i){
            Point3 expected = transform.apply(origVerts[i]);
            for(int j=0;j<3;++j){
                if(std::abs(transformedMesh.vertices[i].coords(j) - expected.coords(j)) > 1e-6f){
                    verticesCorrect = false;
                    break;
                }
            }
        }
        printTest("Transformed SurfaceMesh: vertices match expected", verticesCorrect);
    }

    // --- 4. Test rotation transformation (SpecialOrthogonal) ---
    {
        // 90 degree rotation about z-axis
        float theta = static_cast<float>(M_PI) / 2.0f;
        Eigen::Matrix3f rotMat;
        rotMat << std::cos(theta), -std::sin(theta), 0.0f,
                  std::sin(theta),  std::cos(theta), 0.0f,
                  0.0f,             0.0f,            1.0f;
        SpecialOrthogonal<float,3> rotation(rotMat);

        // Plane
        Surface3 rotatedPlane{
            [f=planeSurface.surfaceFunc, &rotation](float u, float v) -> Point3 {
                return rotation.apply(f(u,v));
            },
            planeSurface.surfaceDomain
        };
        auto rotatedMesh = generateSurfaceMesh(rotatedPlane, 2, 2);
        bool rotVertsOK = true;
        auto origVerts = mesh.vertices;
        for(size_t i=0;i<origVerts.size();++i){
            Point3 expected = rotation.apply(origVerts[i]);
            for(int j=0;j<3;++j){
                if(std::abs(rotatedMesh.vertices[i].coords(j) - expected.coords(j)) > 1e-6f){
                    rotVertsOK = false;
                    break;
                }
            }
        }
        printTest("Rotated Plane SurfaceMesh: vertices match expected", rotVertsOK);

        // Paraboloid
        Surface3 rotatedParaboloid{
            [f=paraboloidSurface.surfaceFunc, &rotation](float u, float v) -> Point3 {
                return rotation.apply(f(u,v));
            },
            paraboloidSurface.surfaceDomain
        };
        auto rotatedParaboloidMesh = generateSurfaceMesh(rotatedParaboloid, 3, 3);
        bool rotParaOK = true;
        auto origParaVerts = paraboloidMesh.vertices;
        for(size_t i=0;i<origParaVerts.size();++i){
            Point3 expected = rotation.apply(origParaVerts[i]);
            for(int j=0;j<3;++j){
                if(std::abs(rotatedParaboloidMesh.vertices[i].coords(j) - expected.coords(j)) > 1e-6f){
                    rotParaOK = false;
                    break;
                }
            }
        }
        printTest("Rotated Paraboloid SurfaceMesh: vertices match expected", rotParaOK);
    }

    // --- 5. Test scaling transformation (Affine::scale) ---
    {
        Eigen::Vector3f scaleVec;
        scaleVec << 2.0f, 0.5f, 3.0f;
        auto scaleXform = Affine<float,3>::scale(scaleVec);

        // Plane
        Surface3 scaledPlane{
            [f=planeSurface.surfaceFunc, &scaleXform](float u, float v) -> Point3 {
                return scaleXform.apply(f(u,v));
            },
            planeSurface.surfaceDomain
        };
        auto scaledMesh = generateSurfaceMesh(scaledPlane, 2, 2);
        bool scaleVertsOK = true;
        auto origVerts = mesh.vertices;
        for(size_t i=0;i<origVerts.size();++i){
            Point3 expected = scaleXform.apply(origVerts[i]);
            for(int j=0;j<3;++j){
                if(std::abs(scaledMesh.vertices[i].coords(j) - expected.coords(j)) > 1e-6f){
                    scaleVertsOK = false;
                    break;
                }
            }
        }
        printTest("Scaled Plane SurfaceMesh: vertices match expected", scaleVertsOK);

        // Paraboloid
        Surface3 scaledParaboloid{
            [f=paraboloidSurface.surfaceFunc, &scaleXform](float u, float v) -> Point3 {
                return scaleXform.apply(f(u,v));
            },
            paraboloidSurface.surfaceDomain
        };
        auto scaledParaboloidMesh = generateSurfaceMesh(scaledParaboloid, 3, 3);
        bool scaleParaOK = true;
        auto origParaVerts = paraboloidMesh.vertices;
        for(size_t i=0;i<origParaVerts.size();++i){
            Point3 expected = scaleXform.apply(origParaVerts[i]);
            for(int j=0;j<3;++j){
                if(std::abs(scaledParaboloidMesh.vertices[i].coords(j) - expected.coords(j)) > 1e-6f){
                    scaleParaOK = false;
                    break;
                }
            }
        }
        printTest("Scaled Paraboloid SurfaceMesh: vertices match expected", scaleParaOK);
    }

    // --- New: Define and test a simple closed mesh (cube) ---
    {
        SurfaceMesh<float,3> cubeMesh;
        // 8 vertices of a unit cube
        cubeMesh.vertices = {
            Point3{0.0f, 0.0f, 0.0f}, // 0
            Point3{1.0f, 0.0f, 0.0f}, // 1
            Point3{1.0f, 1.0f, 0.0f}, // 2
            Point3{0.0f, 1.0f, 0.0f}, // 3
            Point3{0.0f, 0.0f, 1.0f}, // 4
            Point3{1.0f, 0.0f, 1.0f}, // 5
            Point3{1.0f, 1.0f, 1.0f}, // 6
            Point3{0.0f, 1.0f, 1.0f}  // 7
        };
        // 12 triangular faces (two per cube face)
        cubeMesh.faces = {
            {0,1,2}, {0,2,3}, // bottom
            {4,5,6}, {4,6,7}, // top
            {0,1,5}, {0,5,4}, // front
            {1,2,6}, {1,6,5}, // right
            {2,3,7}, {2,7,6}, // back
            {3,0,4}, {3,4,7}  // left
        };
        printTest("CubeMesh: hasClosedTopology() == true", cubeMesh.hasClosedTopology().isClosed);

        // Export cube mesh to OBJ file
        std::string cubePath = "output/cubeMesh.obj";
        auto cubeObjStr = cubeMesh.exportOBJ();
        if (cubeObjStr) {
            std::ofstream ofs(cubePath);
            if (ofs) {
                ofs << *cubeObjStr;
                ofs.close();
                std::cout << "Exported cube mesh to " << cubePath << std::endl;
            } else {
                std::cerr << "Failed to open file for writing: " << cubePath << std::endl;
            }
        } else {
            std::cerr << "Failed to export cube mesh OBJ string." << std::endl;
        }
    }

    // --- New: Define a closed parametric surface (torus) ---
    {
        auto torusFunc = [](float u, float v) -> Point3 {
            // u,v in [0,1]
            // Major radius R, minor radius r
            float R = 1.0f;
            float r = 0.3f;
            float twoPi = 2.0f * static_cast<float>(M_PI);
            float theta = u * twoPi; // angle around major circle
            float phi = v * twoPi;   // angle around minor circle
            float x = (R + r * std::cos(phi)) * std::cos(theta);
            float y = (R + r * std::cos(phi)) * std::sin(theta);
            float z = r * std::sin(phi);
            return Point3{x, y, z};
        };
        Surface3 torusSurface{torusFunc, {{0.0f, 0.0f}, {1.0f, 1.0f}}};

        // Generate mesh with 30x20 vertices for decent resolution, with periodic wrapping
        auto torusMesh = generatePeriodicWrappedMesh(torusSurface, 30, 20, true, true);

        // Validate closed topology for torus mesh
        printTest("TorusMesh: hasClosedTopology() == true", torusMesh.hasClosedTopology().isClosed);

        // Export torus mesh to OBJ file
        std::string torusPath = "output/torusMesh.obj";
        auto torusObjStr = torusMesh.exportOBJ();
        if (torusObjStr) {
            std::ofstream ofs(torusPath);
            if (ofs) {
                ofs << *torusObjStr;
                ofs.close();
                std::cout << "Exported torus mesh to " << torusPath << std::endl;
            } else {
                std::cerr << "Failed to open file for writing: " << torusPath << std::endl;
            }
        } else {
            std::cerr << "Failed to export torus mesh OBJ string." << std::endl;
        }
    }

    // --- New: Define and test a swept tube surface ---
    {
        using Curve3 = Curve<float,3>;

        // Cross-section function: circle in XY plane
        auto circleCrossSection = [](float v) -> Point<float,2> {
            float angle = 2.0f * static_cast<float>(M_PI) * v;
            return Point<float,2>{std::cos(angle), std::sin(angle)};
        };

        // Spine along Z-axis
        Curve3 pathCurve([](float s) -> Point3 {
            return Point3{0.0f, 0.0f, 2.0f * s};
        }, 0.0f, 1.0f);

        // Sweep the circle along the path
        auto sweptSurface = Surface3::sweepSurface(pathCurve, circleCrossSection);

        // Generate mesh with 30 segments along circle, 20 along path, closed topology true
        auto sweptTubeMesh = generateSurfaceMesh(sweptSurface, 30, 20);

        // Print closed topology status
        printTest("SweptTubeMesh: hasClosedTopology() == false", !sweptTubeMesh.hasClosedTopology().isClosed);

        // Export swept tube mesh to OBJ file
        std::string sweptTubePath = "output/sweptTube.obj";
        auto sweptTubeObjStr = sweptTubeMesh.exportOBJ();
        if (sweptTubeObjStr) {
            std::ofstream ofs(sweptTubePath);
            if (ofs) {
                ofs << *sweptTubeObjStr;
                ofs.close();
                std::cout << "Exported swept tube mesh to " << sweptTubePath << std::endl;
            } else {
                std::cerr << "Failed to open file for writing: " << sweptTubePath << std::endl;
            }
        } else {
            std::cerr << "Failed to export swept tube mesh OBJ string." << std::endl;
        }
    }
    
    testSurfaceMeshing();

}

