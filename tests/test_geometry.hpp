#pragma once

#include <iostream>
#include "test_utilities.hpp"
#include "test_point.hpp"
#include "test_line.hpp"
#include "test_segment.hpp"
#include "test_curve.hpp"
#include "test_plane.hpp"
#include "test_face.hpp"
//#include "test_surface.hpp"
#include "test_volume.hpp"



inline void testGeometry() {
    std::cout << "\n📐 Testing Geometry Logic\n";
    //testPoint();
    //testLine();
    //testSegment();
    testCurve();
    //testPlane();
    //testFace();
    //testSurface();
    //testVolume();
    
    

}
