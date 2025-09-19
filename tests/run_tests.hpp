#pragma once
#include "test_algebra.hpp"
#include "test_transform.hpp"
#include "test_geometry.hpp"


namespace euclid::tests {

inline void runAllTests() {
    std::cout << "\n🧪 Running Euclid Unit Tests...\n";
    
    //testAlgebra();
    testTransform();
    testGeometry();
    
    
    std::cout << "--- End of Euclid Unit Tests ---\n";
}

} // namespace euclid::tests
