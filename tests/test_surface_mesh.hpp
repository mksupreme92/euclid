#include <iostream>
#include <fstream>
#include <filesystem>
#include "geometry/surface.hpp"

namespace euclid::geometry {

inline void testSurfaceMeshing() {
    std::cout << "\n∯ Testing Surface Meshing\n";
    
    using Point3 = euclid::geometry::Point<float,3>;
    using Surface3 = euclid::geometry::Surface<float,3>;
    using SurfaceMesh3 = euclid::geometry::SurfaceMesh<float,3>;
    
    // --- Unit square surface ---
    auto squareFunc = [](float u, float v) -> Point3 {
        return Point3{u, v, 0.0f};
    };
    Surface3 squareSurface(squareFunc, {{0.0f,0.0f},{1.0f,1.0f}});
    SurfaceMesh3 squareMesh = euclid::geometry::generateSurfaceMesh(squareSurface, 2, 2);
    
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
    SurfaceMesh3 paraMesh = euclid::geometry::generateSurfaceMesh(paraSurface, 5, 5);
    
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
    using Curve3 = euclid::geometry::Curve<float,3>;
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
    SurfaceMesh3 trefoilSweepMesh = euclid::geometry::generateSurfaceMesh(
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
    using Curve3 = euclid::geometry::Curve<float,3>;

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
    SurfaceMesh3 sweepMesh = euclid::geometry::generateSurfaceMesh(sweepSurface, 20, crossSectionSteps);

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
    SurfaceMesh3 saddleMesh = euclid::geometry::generateSurfaceMesh(saddleSurface, 20, 20);

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

}
