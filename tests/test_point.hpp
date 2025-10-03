#include <iostream>
#include "geometry/point.hpp"
#include "algebra/transform.hpp"
#include "tests/test_utilities.hpp"

using namespace Euclid::Geometry;
using namespace Euclid::Algebra;
using namespace Euclid::Tests;

inline void testPoint() {
    std::cout << "\n📍 Testing Point Distance\n\n";

    using Point2 = Point<float, 2>;

    Point2 p1{1.0f, 2.0f};
    Point2 p2{4.0f, 6.0f};

    // Distance test
    float dist = p1.distanceTo(p2);

    Euclid::Tolerance tol;
    printTest("Distance between points",
              std::abs(dist - 5.0f) <= tol.evaluateEpsilon(dist));
    
    // Midpoint test
    Point2 mid = p1.midpoint(p2);
    printTest("Midpoint", mid == Point2{2.5f, 4.0f});

    // Vector difference test
    // Note - does not use tolerance model, uses Eigen's isApprox
    using Vec2 = Eigen::Matrix<float, 2, 1>;
    Vec2 diff = p2 - p1;
    printTest("Vector difference", diff.isApprox(Vec2{3.0f, 4.0f}));

    // --- Point Transform Tests ---
    Eigen::Matrix<float,2,2> linearPart;
    linearPart << 2.0f, 0.0f,
                  0.0f, 3.0f;
    
    Eigen::Matrix<float,2,1> translation;
    translation << 1.0f, -1.0f;

    Affine<float,2> transform(linearPart, translation);

    Point2 p_transformed = transform.apply(p1);
    Point2 expected{3.0f, 5.0f};
    printTest("Point transform (scale + translate)", p_transformed == expected);

    // --- Additional Point Transform Tests ---

    // 1. Pure translation
    Eigen::Matrix<float,2,2> identityMat = Eigen::Matrix<float,2,2>::Identity();
    Eigen::Matrix<float,2,1> translation2;
    translation2 << -2.0f, 3.0f;
    Affine<float,2> transformTranslate(identityMat, translation2);
    Point2 p_translated = transformTranslate.apply(p1);
    Point2 expectedTranslate{-1.0f, 5.0f};
    printTest("Point transform (pure translation)", p_translated == expectedTranslate);

    // 2. Pure scaling
    Eigen::Matrix<float,2,2> scaleOnly;
    scaleOnly << 3.0f, 0.0f,
                    0.0f, 0.5f;
    Eigen::Matrix<float,2,1> zeroTranslation = Eigen::Matrix<float,2,1>::Zero();
    Affine<float,2> transformScale(scaleOnly, zeroTranslation);
    Point2 p_scaled = transformScale.apply(p1);
    Point2 expectedScale{3.0f, 1.0f};
    printTest("Point transform (pure scaling)", p_scaled == expectedScale);

    // 3. Scale + translation combination
    Eigen::Matrix<float,2,2> scaleMat2;
    scaleMat2 << -1.0f, 0.0f,
                0.0f, 2.0f;
    Eigen::Matrix<float,2,1> translation3;
    translation3 << 5.0f, -3.0f;
    Affine<float,2> transformScaleTranslate(scaleMat2, translation3);
    Point2 p_transformed2 = transformScaleTranslate.apply(p1);
    Point2 expectedScaleTranslate{4.0f, 1.0f};
    printTest("Point transform (scale + translation combo)", p_transformed2 == expectedScaleTranslate);

    // --- Rotate point about another point ---

    Point2 pivot{2.0f, 1.0f};
    float angle = M_PI / 2.0f; // 90 degrees
    Eigen::Matrix2f rotMat;
    rotMat << std::cos(angle), -std::sin(angle),
                 std::sin(angle),  std::cos(angle);
    SpecialOrthogonal<float,2> rotation(rotMat);

    Eigen::Matrix2f linear = rotation.A;
    Eigen::Matrix<float,2,1> translationPivot = pivot.coords - linear * pivot.coords;
    Affine<float,2> rotateAroundPivot(linear, translationPivot);

    Point2 p_rotated = rotateAroundPivot.apply(p1);
    Point2 expectedRotated{
        static_cast<float>(pivot[0] - (p1[1] - pivot[1])),
        static_cast<float>(pivot[1] + (p1[0] - pivot[0]))
    };
    printTest("Point rotation about another point", p_rotated == expectedRotated);

    // --- Compact rotation about pivot test ---

    Point2 p{1.0f, 2.0f};
    Point2 pivotCompact{2.0f, 1.0f};
    float angleCompact = M_PI / 2.0f;
    Eigen::Matrix2f rotMatCompact;
    rotMatCompact << std::cos(angleCompact), -std::sin(angleCompact),
                     std::sin(angleCompact),  std::cos(angleCompact);
    SpecialOrthogonal<float,2> rotationCompact(rotMatCompact);
    Eigen::Matrix2f linearCompact = rotationCompact.A;
    Eigen::Matrix<float,2,1> translationCompact = pivotCompact.coords - linearCompact * pivotCompact.coords;
    Affine<float,2> rotateAroundPivotCompact(linearCompact, translationCompact);
    Point2 p_rotatedCompact = rotateAroundPivotCompact.apply(p);
    Point2 expectedCompact{
        static_cast<float>(pivotCompact[0] - (p[1] - pivotCompact[1])),
        static_cast<float>(pivotCompact[1] + (p[0] - pivotCompact[0]))
    };
    printTest("Point rotation about pivot compact", p_rotatedCompact == expectedCompact);
}
