#pragma once
#include "test_geometry.hpp"


//using namespace Euclid::Tests;

inline void runAllTests() {
    std::cout << "\n🧪 Running Euclid Unit Tests...\n";
    
    testGeometry();
    
    
    std::cout << "--- End of Euclid Unit Tests ---\n";
}
