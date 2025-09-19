#pragma once

#include <iostream>
#include "algebra/transform.hpp"
#include "test_utilities.hpp"
#include <cmath>
#include <Eigen/Dense>

using namespace euclid::algebra;
using namespace euclid::tests;

inline void testTransform() {
    using Vec2 = Eigen::Matrix<float, 2, 1>;
    using Vec3 = Eigen::Matrix<float, 3, 1>;
    using Affine2 = Affine<float, 2>;
    using Affine3 = Affine<float, 3>;

    std::cout << "\n🔄 Testing Transform Logic\n\n";

    // Vector Translation test (2D)
    Vec2 v2{1.0f, 2.0f};
    Vec2 t2{3.0f, 5.0f};
    Affine2 trans2 = Affine2::translate(t2);
    Vec2 v2t = trans2.apply(v2);
    printTest("Vector 2D Translation", approxEqualVector(v2t, Vec2{4.0f, 7.0f}));

    // Vector Scaling test (2D)
    Vec2 s2{2.0f, 3.0f};
    Affine2 scale2 = Affine2::scale(s2);
    Vec2 v2s = scale2.apply(v2);
    printTest("Vector 2D Scaling", approxEqualVector(v2s, Vec2{2.0f, 6.0f}));

    // 2D Shear test using raw Eigen matrix
    Eigen::Matrix<float, 2, 2> shearMat;
    shearMat << 1.0f, 1.0f,
                0.0f, 1.0f;
    Vec2 shearTrans{0.0f, 0.0f};
    Affine2 shear2(shearMat, shearTrans);
    Vec2 v2sh = shear2.apply(v2);
    printTest("Vector 2D Shear", approxEqualVector(v2sh, Vec2{3.0f, 2.0f}));

    // 3D rotation around Z axis using raw Eigen matrix
    Vec3 v3{1.0f, 0.0f, 0.0f};
    float theta = float(M_PI) / 2.0f; // 90 degrees
    Eigen::Matrix<float, 3, 3> rotzMat;
    rotzMat << std::cos(theta), -std::sin(theta), 0.0f,
                std::sin(theta),  std::cos(theta), 0.0f,
                0.0f,            0.0f,           1.0f;
    Vec3 rotzTrans{0.0f, 0.0f, 0.0f};
    Affine3 rotz(rotzMat, rotzTrans);
    Vec3 v3r = rotz.apply(v3);
    printTest("Vector 3D Rotation Z", approxEqualVector(v3r, Vec3{0.0f, 1.0f, 0.0f}, 1e-5f));

    // 2D Composite translation + scaling test using raw Eigen matrix and Affine constructor
    Eigen::Matrix<float, 2, 2> scaleMat;
    scaleMat << s2[0], 0.0f,
                0.0f,  s2[1];
    Vec2 zeroVec2{0.0f, 0.0f};
    Affine2 scaleAffine(scaleMat, zeroVec2);
    Affine2 transAffine(Eigen::Matrix<float, 2, 2>::Identity(), t2);
    Eigen::Matrix<float, 2, 2> compositeLinear = transAffine.getLinear() * scaleAffine.getLinear();
    Vec2 compositeTranslation = transAffine.apply(scaleAffine.getTranslation());
    Affine2 composite2(compositeLinear, compositeTranslation);
    Vec2 v2comp = composite2.apply(v2);
    printTest("Vector 2D Composite Translation + Scaling", approxEqualVector(v2comp, Vec2{5.0f, 11.0f}));

    // Normal transform using applyLinear with rotz defined above
    Vec3 normal{0.0f, 0.0f, 1.0f};
    Vec3 transformedNormal = rotz.applyLinear(normal);
    printTest("Normal Transform applyLinear", approxEqualVector(transformedNormal, Vec3{0.0f, 0.0f, 1.0f}, 1e-5f));

    // Inverse transform test for the 2D composite
    Affine2 invComposite = composite2.inverse();
    Vec2 v2inv = invComposite.apply(v2comp);
    printTest("Inverse Transform 2D Composite", approxEqualVector(v2inv, v2, 1e-5f));
}
