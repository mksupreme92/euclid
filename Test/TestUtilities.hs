-- 🜁 Euclid Geometry Kernel
-- Copyright (c) 2025 Joseph Shaw Petersen 🦅
-- Licensed under the BSD 3-Clause License. See LICENSE for details.

module Test.TestUtilities (
    printTest
    , approxEqual
    , approxEqualList
    , approxEqualVector
    , vectorFromVertex
) where

import Algebra.Vector
import Geometry.Vertex (Vertex(..))

printTest :: String -> Bool -> IO ()
printTest name result = putStrLn $ (if result then "✅ " else "❌ ") ++ name


-- Helper function for comparing vectors
approxEqual :: (Floating a, Ord a) => a -> a -> Bool
approxEqual a b = abs (a - b) < 1e-9

approxEqualList :: [Double] -> [Double] -> Bool
approxEqualList xs ys = and $ zipWith approxEqual xs ys

approxEqualVector :: Vector Double -> Vector Double -> Bool
approxEqualVector (Vector v1) (Vector v2) = approxEqualList v1 v2


-- Convert a Vertex to its underlying Vector
vectorFromVertex :: Vertex a -> Vector a
vectorFromVertex (Vertex v) = v
