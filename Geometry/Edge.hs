-- 🜁 Euclid Geometry Kernel
-- Copyright (c) 2025 Joseph Shaw Petersen 🦅
-- Licensed under the BSD 3-Clause License. See LICENSE for details.

module Geometry.Edge (
    Edge(..),
    lineFromVertexAndVector,
    rayFromVertexAndVector,
    segmentFromVertices,
    defineEdge
) where

import Algebra.Vector
import Algebra.Matrix
import Algebra.Metric
import Geometry.Vertex


-- | An edge can be an infinite line, a ray, or a segment.
data Edge a
  = InfiniteLine (Vertex a) (Vector a)
  | Ray         (Vertex a) (Vector a)
  | Segment     (Vertex a) (Vertex a)
  deriving (Show, Eq)

lineFromVertexAndVector :: (Eq a, Num a) => Vertex a -> Vector a -> Maybe (Edge a)
lineFromVertexAndVector p dir
  | isValidVector dir = Just (InfiniteLine p dir)
  | otherwise         = Nothing

-- | Unsafe constructor: ray from origin vertex and direction
rayFromVertexAndVector :: (Eq a, Num a) => Vertex a -> Vector a -> Maybe (Edge a)
rayFromVertexAndVector p dir
  | isValidVector dir = Just (Ray p dir)
  | otherwise         = Nothing

-- | Unsafe constructor: segment from two vertices
segmentFromVertices :: Num a => Vertex a -> Vertex a -> Maybe (Edge a)
segmentFromVertices v1@(Vertex p1) v2@(Vertex p2) = do
  _ <- vectorSub p2 p1  -- just validate vector subtraction is possible
  return (Segment v1 v2)

-- | Safe constructor with space validation for all edge types
defineEdge :: Num a => Metric a -> Edge a -> Maybe (Edge a)
defineEdge (Metric m) edge =
  let (rows, cols) = matrixDimensions m
      check v = length (vectorToList v) == rows && cols == rows
  in case edge of
       InfiniteLine (Vertex p) dir
         | check p && check dir -> Just edge
         | otherwise            -> Nothing
       Ray (Vertex p) dir
         | check p && check dir -> Just edge
         | otherwise            -> Nothing
       Segment (Vertex p1) (Vertex p2)
         | check p1 && check p2 -> Just edge
         | otherwise            -> Nothing

-- | Utility to check vector is not zero-length (optional for stricter validation)
isValidVector :: (Eq a, Num a) => Vector a -> Bool
isValidVector (Vector xs) = any (/= 0) xs
