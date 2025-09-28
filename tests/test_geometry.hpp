#pragma once

#include <iostream>
#include "test_point.hpp"
#include "test_line.hpp"
#include "test_curve.hpp"
#include "test_plane.hpp"
#include "test_face.hpp"
#include "test_surface.hpp"
#include "test_surface_mesh.hpp"
#include "test_utilities.hpp"
#include "config.hpp"

using namespace euclid::geometry;
using namespace euclid::tests;


// Call in testGeometry
inline void testGeometry() {
    /*
    std::cout << "\n📐 Testing Geometry Logic\n";
    testPoint();
    testPointRotationAboutPivotCompact();
    testLine();
    testLineMeasureAngle();
    testCurve();
    testPlane();
    testFace();
    testSurface();
    */
    testSurfaceMeshing();
    
}
