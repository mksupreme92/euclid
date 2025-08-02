-- 🜁 Euclid Geometry Kernel
-- Copyright (c) 2025 Joseph Shaw Petersen 🦅
-- Licensed under the BSD 3-Clause License. See LICENSE for details.

module Test.TestAlgebra (testAlgebra) where

import Algebra.Vector
import Algebra.Matrix
import Algebra.Metric
import Algebra.Transform
  ( Translatable(..)
  , Scalable(..)
  , Rotatable(..)
  , Skewable(..)
  , Reflectable(..)
  , Projectable(..)
  , rotationMatrix )
  
import Test.TestUtilities

import Control.Monad (when)

import Debug.Trace (trace, traceShow)


-- Euclidean metric for 3D (and general) tests
euclidean3D :: Metric Double
euclidean3D = Metric (identityMatrix 3)

testAlgebra :: IO ()
testAlgebra = do
  putStrLn "\n🧮 Testing Core Linear Algebra Modules\n"
  vectorTests
  matrixTests
  metricTests

vectorTests :: IO ()
vectorTests = do
  putStrLn "\n↑ Testing Vector Module\n"
  -- Explicitly type all vectors as Vector Double for clarity and type safety
  let v1, v2, v3 :: Vector Double
      v1 = vectorFromList [1.0, 2.0, 3.0]
      v2 = vectorFromList [4.0, 5.0, 6.0]
      v3 = vectorFromList [1.0, 2.0]

  -- vectorAdd
  printTest "vectorAdd valid" $
    vectorAdd v1 v2 == Just (vectorFromList [5.0, 7.0, 9.0])
  printTest "vectorAdd mismatched" $
    vectorAdd v1 v3 == Nothing

  -- vectorSub
  printTest "vectorSub valid" $
    vectorSub v2 v1 == Just (vectorFromList [3.0, 3.0, 3.0])
  printTest "vectorSub mismatched" $
    vectorSub v1 v3 == Nothing

  -- scalarMul
  printTest "scalarMul 2 * v1" $
    scalarMul 2 v1 == vectorFromList [2.0, 4.0, 6.0]

  -- dotProduct
  printTest "dotProduct valid" $
    dotProduct v1 v2 == Just 32.0
  printTest "dotProduct mismatched" $
    dotProduct v1 v3 == Nothing

  -- norm
  printTest "norm of v1" $
    approxEqual (norm v1) (sqrt 14)

  -- normalize
  let nv1 = normalize v1
  printTest "normalize v1 has norm 1" $
    approxEqual (norm nv1) 1.0

  -- translate (generic, type-consistent)
  let t1, t2 :: Vector Double
      t1 = vectorFromList [1.0, 1.0, 1.0]
      t2 = vectorFromList [1.0, 2.0, 3.0]
  printTest "translate valid" $
    case translate euclidean3D t1 t2 of
      Just v -> v == vectorFromList [2.0, 3.0, 4.0]
      Nothing -> False

  let t3, t4 :: Vector Double
      t3 = vectorFromList [1.0, 1.0]
      t4 = vectorFromList [1.0, 2.0, 3.0]
  printTest "translate mismatched" $
    case translate euclidean3D t3 t4 of
      Just _ -> False
      Nothing -> True

  -- scale (generic, type-consistent)
  let s1, s2 :: Vector Double
      s1 = vectorFromList [2.0, 0.5, -1.0]
      s2 = vectorFromList [1.0, 2.0, 3.0]
      identity2D = identityMatrix 3
  printTest "scale valid" $
    case scale (Metric identity2D) s1 s2 of
      Just v -> v == vectorFromList [2.0, 1.0, -3.0]
      Nothing -> False

  let s3, s4 :: Vector Double
      s3 = vectorFromList [2.0]
      s4 = vectorFromList [1.0, 2.0]
  printTest "scale mismatched" $
    case scale (Metric identity2D) s3 s4 of
      Just _ -> False
      Nothing -> True

  -- rotate (using explicit Matrix Double and Vector Double)
  -- 2D rotation by 90 degrees
  let angle2D = pi / 2 :: Double
      rotationMat2D :: Matrix Double
      rotationMat2D = matrixFromList [[cos angle2D, - sin angle2D],
                                      [sin angle2D,   cos angle2D]]
      v2d :: Vector Double
      v2d = vectorFromList [1.0, 0.0]
      expectedRotated2D = [0.0, 1.0]
  printTest "rotateVector (2D matrix rotation 90 deg)" $
    case rotate (Metric (identityMatrix 2)) rotationMat2D v2d of
      Just rotated -> approxEqualList (vectorToList rotated) expectedRotated2D
      Nothing -> False

  -- 3D rotation about Z axis by 90 degrees
  let angle3D = pi / 2 :: Double
      rotationMat3D :: Matrix Double
      rotationMat3D = matrixFromList
        [[cos angle3D, -sin angle3D, 0],
         [sin angle3D,  cos angle3D, 0],
         [0,            0,           1]]
      v3d :: Vector Double
      v3d = vectorFromList [1.0, 0.0, 0.0]
      expectedRotated3D = [0.0, 1.0, 0.0]
  printTest "rotateVector (3D Z-axis rotation 90 deg)" $
    case rotate (Metric (identityMatrix 3)) rotationMat3D v3d of
      Just rotated -> approxEqualList (vectorToList rotated) expectedRotated3D
      Nothing -> False

  -- 4D rotation in XY plane by 90 degrees
  let angle4D = pi / 2 :: Double
      rotationMat4D :: Matrix Double
      rotationMat4D = matrixFromList
        [[cos angle4D, -sin angle4D, 0, 0],
         [sin angle4D,  cos angle4D, 0, 0],
         [0,            0,           1, 0],
         [0,            0,           0, 1]]
      v4d :: Vector Double
      v4d = vectorFromList [1.0, 0.0, 0.0, 0.0]
      expectedRotated4D = [0.0, 1.0, 0.0, 0.0]
  printTest "rotateVector (4D XY-plane rotation 90 deg)" $
    case rotate (Metric (identityMatrix 4)) (matrixFromList rotationMat4D) v4d of
      Just rotated -> approxEqualList (vectorToList rotated) expectedRotated4D
      Nothing -> False

  -- 5D rotation in first two dims by 90 degrees
  let angle5D = pi / 2 :: Double
      rotationMat5D :: [[Double]]
      rotationMat5D =
        [[cos angle5D, -sin angle5D, 0, 0, 0],
         [sin angle5D,  cos angle5D, 0, 0, 0],
         [0,            0,           1, 0, 0],
         [0,            0,           0, 1, 0],
         [0,            0,           0, 0, 1]]
      v5d :: Vector Double
      v5d = vectorFromList [1.0, 0.0, 0.0, 0.0, 0.0]
      expectedRotated5D = [0.0, 1.0, 0.0, 0.0, 0.0]
  printTest "rotateVector (5D plane rotation 90 deg)" $
    case rotate (Metric (identityMatrix 5)) (matrixFromList rotationMat5D) v5d of
      Just rotated -> approxEqualList (vectorToList rotated) expectedRotated5D
      Nothing -> False
      
  -- scale then rotate composition
  let metric = Metric (identityMatrix 2)
      angle = pi / 2 :: Double
      sVec :: Vector Double
      sVec = vectorFromList [2.0, 1.0]
      xVec :: Vector Double
      xVec = vectorFromList [1.0, 0.0]
      rotMat = matrixFromList [[cos angle, -sin angle],
                               [sin angle,  cos angle]]
      expectedComposed :: Vector Double
      expectedComposed = vectorFromList [0.0, 2.0]

  printTest "scale then rotate composition" $
    case scale metric sVec xVec >>= rotate metric rotMat of
      Just result -> approxEqualVector result expectedComposed
      Nothing -> False

  -- Skew test (orthogonal: identity matrix)
  let metric = Metric [[1, 0], [0, 1]]
      shearMatrix = [[1.0, 0.0], [0.0, 1.0]]
      vSkew = Vector [2.0, 3.0]
      expectedSkewed = Vector [2.0, 3.0]
  printTest "Skew" $
    skew metric shearMatrix vSkew == Just expectedSkewed

  -- Skew should fail with non-metric-preserving matrix
  let invalidShearMatrix = [[2.0, 0.0], [0.0, 1.0]]  -- scales x-dimension
  printTest "Skew (non-preserving should fail)" $
    skew metric invalidShearMatrix vSkew == Nothing

  -- Reflect test (reflect over Y axis: x -> -x)
  let reflectMatrix = [[1, 0], [0, -1]]
      vReflect = Vector [2, 3]
      expectedReflected = Vector [2, -3]
  printTest "Reflect" $
    reflect metric reflectMatrix vReflect == Just expectedReflected

  -- Project test (orthogonal: identity matrix)
  let projectMatrix = [[1.0, 0.0], [0.0, 1.0]]
      vProject = Vector [2.0, 3.0]
      expectedProjected = Vector [2.0, 3.0]
  printTest "Project" $
    project metric projectMatrix vProject == Just expectedProjected

  -- Project should fail with non-metric-preserving matrix
  let invalidProjectMatrix = [[2.0, 0.0], [0.0, 2.0]]  -- scales both axes
  printTest "Project (non-preserving should fail)" $
    project metric invalidProjectMatrix vProject == Nothing

  -- Non-Euclidean metric transformation tests
  let metricNE = Metric [[1, 0], [0, 2]]
      vNE = Vector [2.0, 3.0]
      offsetNE = Vector [1.0, 1.0]
      scaleVecNE = Vector [2.0, -1.0]
      -- rotNE is a G-orthogonal for metricNE = diag(1,2)
      rotNE = matrixFromList [[1, 0], [0, -1]] -- G-orthogonal for metricNE = diag(1,2)
      -- shearNE is the identity matrix (metric-preserving placeholder)
      shearNE = [[1.0, 0.0],
                 [0.0, 1.0]]
      reflectNE = [[-1, 0], [0, 1]]
      -- projectNE is the identity (metric-compatible projection)
      projectNE = [[1.0, 0.0],
                   [0.0, 1.0]]

  printTest "translate (non-Euclidean)" $
    translate metricNE offsetNE vNE == Just (Vector [3.0, 4.0])

  -- Scaling with non-Euclidean metric that does NOT preserve the metric should fail
  printTest "scale (non-Euclidean, non-preserving)" $
      scale metricNE scaleVecNE vNE == Nothing
  -- Scaling with non-Euclidean metric that does preserve the metric should also fail
  let uniformScaleVecNE = Vector [2.0, 2.0]  -- uniform scaling preserves metric
  printTest "scale (non-Euclidean, preserving uniform)" $
      scale metricNE uniformScaleVecNE vNE == Nothing

  printTest "rotate (non-Euclidean)" $
    (let ok = isMetricPreserving metricNE rotNE
     in trace ("isMetricPreserving = " ++ show ok) $
        case rotate metricNE rotNE vNE of
          Just result ->
            let expected = Vector [2.0, -3.0]
            in traceShow result $ approxEqualVector result expected
          Nothing -> trace "rotate returned Nothing" False)

  -- shearNE and projectNE are identity, so output equals input
  printTest "skew (non-Euclidean)" $
    skew metricNE shearNE vNE == Just (Vector [2.0, 3.0])

  printTest "reflect (non-Euclidean)" $
    reflect metricNE reflectNE vNE == Just (Vector [-2.0, 3.0])

  printTest "project (non-Euclidean)" $
    project metricNE projectNE vNE == Just (Vector [2.0, 3.0])

matrixTests :: IO ()
matrixTests = do
  putStrLn "\n🔢 Testing Matrix Module\n"

  let m1 = matrixFromList [[1, 2], [3, 4]]
      m2 = matrixFromList [[5, 6], [7, 8]]
      id2 = identityMatrix 2
      v = [1, 1]

  -- identityMatrix
  printTest "identityMatrix 2x2" $
    id2 == [[1,0],[0,1]]

  -- matrixVectorProduct
  printTest "matrixVectorProduct valid" $
    matrixVectorProduct m1 v == Just [3, 7]

  printTest "matrixVectorProduct dim mismatch" $
    matrixVectorProduct m1 [1,1,1] == Nothing

  -- transpose
  printTest "transpose 2x2" $
    transpose m1 == [[1,3],[2,4]]

  -- matrixMultiply
  printTest "matrixMultiply valid" $
    matrixMultiply m1 m2 == Just [[19,22],[43,50]]

  printTest "matrixMultiply dim mismatch" $
    matrixMultiply m1 [[1]] == Nothing

metricTests :: IO ()
metricTests = do
  putStrLn "\n📐 Testing Metric Module\n"
  -- Metric tests
  let gMatrix = matrixFromList [[1.0, 0.0], [0.0, 1.0]]
      metric = Metric gMatrix
      u = vectorFromList [1.0, 2.0]
      w = vectorFromList [3.0, 4.0]

  -- innerProduct
  printTest "innerProduct (Euclidean)" $
    innerProduct metric u w == Just (1*3 + 2*4)  -- 11.0

  -- distance
  printTest "distance (Euclidean)" $
    case distance metric u w of
      Just d -> approxEqual d (sqrt ((3-1)^2 + (4-2)^2))  -- sqrt(8)
      Nothing -> False

  let gMatrixNonEuclid = matrixFromList [[2.0, 0.0], [0.0, 0.5]]
      metricNE = Metric gMatrixNonEuclid
      uNE = vectorFromList [1.0, 2.0]
      wNE = vectorFromList [3.0, 4.0]

  -- innerProduct with non-Euclidean metric
  printTest "innerProduct (non-Euclidean)" $
    innerProduct metricNE uNE wNE == Just (2*1*3 + 0.5*2*4)  -- 6 + 4 = 10

  -- distance with non-Euclidean metric
  printTest "distance (non-Euclidean)" $
    case distance metricNE uNE wNE of
      Just d -> approxEqual d (sqrt (2*(3-1)^2 + 0.5*(4-2)^2))  -- sqrt(8 + 2) = sqrt(10)
      Nothing -> False

  -- isEuclidean
  printTest "isEuclidean (true)" $
    isEuclidean metric == True

  printTest "isEuclidean (false)" $
    isEuclidean metricNE == False
