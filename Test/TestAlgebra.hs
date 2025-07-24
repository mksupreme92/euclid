module Test.TestAlgebra (testAlgebra) where

import Algebra.Vector
import Algebra.Matrix
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

-- Helper function for comparing vectors
approxEqual :: (Floating a, Ord a) => a -> a -> Bool
approxEqual a b = abs (a - b) < 1e-9

approxEqualList :: [Double] -> [Double] -> Bool
approxEqualList xs ys = and $ zipWith approxEqual xs ys

approxEqualVector :: Vector Double -> Vector Double -> Bool
approxEqualVector (Vector v1) (Vector v2) = approxEqualList v1 v2


testAlgebra :: IO ()
testAlgebra = do
  putStrLn "\n🧮 Testing core linear algebra logic\n"
  vectorTests
  matrixTests

vectorTests :: IO ()
vectorTests = do
  putStrLn "\n↑ Testing vector logic\n"
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
    case translate t1 t2 of
      Just v -> v == vectorFromList [2.0, 3.0, 4.0]
      Nothing -> False

  let t3, t4 :: Vector Double
      t3 = vectorFromList [1.0, 1.0]
      t4 = vectorFromList [1.0, 2.0, 3.0]
  printTest "translate mismatched" $
    case translate t3 t4 of
      Just _ -> False
      Nothing -> True

  -- scale (generic, type-consistent)
  let s1, s2 :: Vector Double
      s1 = vectorFromList [2.0, 0.5, -1.0]
      s2 = vectorFromList [1.0, 2.0, 3.0]
  printTest "scale valid" $
    case scale s1 s2 of
      Just v -> v == vectorFromList [2.0, 1.0, -3.0]
      Nothing -> False

  let s3, s4 :: Vector Double
      s3 = vectorFromList [2.0]
      s4 = vectorFromList [1.0, 2.0]
  printTest "scale mismatched" $
    case scale s3 s4 of
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
    case rotate rotationMat2D v2d of
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
    case rotate rotationMat3D v3d of
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
    case rotate rotationMat4D v4d of
      Just rotated -> approxEqualList (vectorToList rotated) expectedRotated4D
      Nothing -> False

  -- 5D rotation in first two dims by 90 degrees
  let angle5D = pi / 2 :: Double
      rotationMat5D :: Matrix Double
      rotationMat5D = matrixFromList
        [[cos angle5D, -sin angle5D, 0, 0, 0],
         [sin angle5D,  cos angle5D, 0, 0, 0],
         [0,            0,           1, 0, 0],
         [0,            0,           0, 1, 0],
         [0,            0,           0, 0, 1]]
      v5d :: Vector Double
      v5d = vectorFromList [1.0, 0.0, 0.0, 0.0, 0.0]
      expectedRotated5D = [0.0, 1.0, 0.0, 0.0, 0.0]
  printTest "rotateVector (5D plane rotation 90 deg)" $
    case rotate rotationMat5D v5d of
      Just rotated -> approxEqualList (vectorToList rotated) expectedRotated5D
      Nothing -> False
      
  -- scale then rotate composition
  let angle = pi / 2 :: Double
      sVec :: Vector Double
      sVec = vectorFromList [2.0, 1.0]
      xVec :: Vector Double
      xVec = vectorFromList [1.0, 0.0]
      rotMat = matrixFromList [[cos angle, -sin angle],
                               [sin angle,  cos angle]]
      expectedComposed :: Vector Double
      expectedComposed = vectorFromList [0.0, 2.0]

  printTest "scale then rotate composition" $
    case scale sVec xVec >>= rotate rotMat of
      Just result -> approxEqualVector result expectedComposed
      Nothing -> False

  -- Skew test (2D shear: x += y * 1)
  let shearMatrix = matrixFromList [[1, 1], [0, 1]]
  let vSkewed = skew shearMatrix v1
  printTest "skewVector (2D shear x += y)" (vSkewed == Just (vectorFromList [4.0, 2.0]))

  -- Reflect test (reflect over Y axis: x -> -x)
  let reflectMatrix = matrixFromList [[-1, 0], [0, 1]]
  let vReflected = reflect reflectMatrix v1
  printTest "reflectVector (reflect over Y axis)" (vReflected == Just (vectorFromList [-2.0, 2.0]))

  -- Project test (project onto X axis)
  let projectMatrix = matrixFromList [[1, 0], [0, 0]]
  let vProjected = project projectMatrix v1
  printTest "projectVector (onto X axis)" (vProjected == Just (vectorFromList [2.0, 0.0]))

matrixTests :: IO ()
matrixTests = do
  putStrLn "\n🔢 Testing matrix logic\n"

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
