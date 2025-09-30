#include <iostream>
#include "geometry/point.hpp"

inline void testPoint() {
    std::cout << "\n📍 Testing Point Primitive\n\n";

    // Save the old space dimension and set the new one

    using Point2 = Point<float, 2>;
    using Vec2 = Eigen::Matrix<float, 2, 1>;

    Point2 p1{1.0f, 2.0f};
    Point2 p2{4.0f, 6.0f};

    // Distance test
    float dist = p1.distanceTo(p2);
    printTest("Distance between points", std::abs(dist - 5.0f) < 1e-5f);

    // Midpoint test
    Point2 mid = p1.midpoint(p2);
    printTest("Midpoint", mid == Point2{2.5f, 4.0f});

    // Vector difference test
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
    printTest("Point transform (scale + translate)", approxEqual(p_transformed, expected));

    // --- Additional Point Transform Tests ---

    // 1. Pure translation
    Eigen::Matrix<float,2,2> identityMat = Eigen::Matrix<float,2,2>::Identity();
    Eigen::Matrix<float,2,1> translation2;
    translation2 << -2.0f, 3.0f;
    Affine<float,2> transformTranslate(identityMat, translation2);
    Point2 p_translated = transformTranslate.apply(p1);
    Point2 expectedTranslate{-1.0f, 5.0f};
    printTest("Point transform (pure translation)", approxEqual(p_translated, expectedTranslate));

    // 2. Pure scaling
    Eigen::Matrix<float,2,2> scaleOnly;
    scaleOnly << 3.0f, 0.0f,
                 0.0f, 0.5f;
    Eigen::Matrix<float,2,1> zeroTranslation = Eigen::Matrix<float,2,1>::Zero();
    Affine<float,2> transformScale(scaleOnly, zeroTranslation);
    Point2 p_scaled = transformScale.apply(p1);
    Point2 expectedScale{3.0f, 1.0f};
    printTest("Point transform (pure scaling)", approxEqual(p_scaled, expectedScale));

    // 3. Scale + translation combination
    Eigen::Matrix<float,2,2> scaleMat2;
    scaleMat2 << -1.0f, 0.0f,
                  0.0f, 2.0f;
    Eigen::Matrix<float,2,1> translation3;
    translation3 << 5.0f, -3.0f;
    Affine<float,2> transformScaleTranslate(scaleMat2, translation3);
    Point2 p_transformed2 = transformScaleTranslate.apply(p1);
    Point2 expectedScaleTranslate{4.0f, 1.0f};
    printTest("Point transform (scale + translation combo)", approxEqual(p_transformed2, expectedScaleTranslate));

    // --- Rotate point about another point ---

    Point2 pivot{2.0f, 1.0f};
    float angle = M_PI / 2.0f; // 90 degrees
    Eigen::Matrix2f rotMat;
    rotMat << std::cos(angle), -std::sin(angle),
              std::sin(angle),  std::cos(angle);
    euclid::algebra::SpecialOrthogonal<float,2> rotation(rotMat);

    Eigen::Matrix2f linear = rotation.A;
    Eigen::Matrix<float,2,1> translationPivot = pivot.coords - linear * pivot.coords;
    euclid::algebra::Affine<float,2> rotateAroundPivot(linear, translationPivot);

    Point2 p_rotated = rotateAroundPivot.apply(p1);
    Point2 expectedRotated{
        static_cast<float>(pivot[0] - (p1[1] - pivot[1])),
        static_cast<float>(pivot[1] + (p1[0] - pivot[0]))
    };
    printTest("Point rotation about another point", approxEqual(p_rotated, expectedRotated, 1e-6f));


    // Reset to old space dimension to avoid side effects
}

inline void testPointRotationAboutPivotCompact() {

    using Point2 = Point<float,2>;

    Point2 p{1.0f, 2.0f};
    Point2 pivot{2.0f, 1.0f};

    // 90 degrees CCW rotation
    float angle = M_PI / 2.0f;
    Eigen::Matrix2f rotMat;
    rotMat << std::cos(angle), -std::sin(angle),
              std::sin(angle),  std::cos(angle);

    SpecialOrthogonal<float,2> rotation(rotMat);

    // One-liner: rotate about pivot
    Eigen::Matrix2f linear = rotation.A;
    Eigen::Matrix<float,2,1> translationPivot = pivot.coords - linear * pivot.coords;
    Affine<float,2> rotateAroundPivot(linear, translationPivot);

    Point2 p_rotated = rotateAroundPivot.apply(p);

    // Expected manually computed
    Point2 expected{
        static_cast<float>(pivot[0] - (p[1] - pivot[1])),
        static_cast<float>(pivot[1] + (p[0] - pivot[0]))
    };

    if (approxEqual(p_rotated, expected, 1e-6f)) {
        std::cout << "✅ Compact rotation about pivot passed.\n";
    } else {
        std::cout << "❌ Compact rotation about pivot FAILED.\n";
        std::cout << "Got: [" << p_rotated.coords.transpose() << "] ";
        std::cout << "Expected: [" << expected.coords.transpose() << "]\n";
    }

}
