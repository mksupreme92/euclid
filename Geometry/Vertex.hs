-- 🜁 Euclid Geometry Kernel
-- Copyright (c) 2025 Joseph Shaw Petersen 🦅
-- Licensed under the BSD 3-Clause License. See LICENSE for details.

module Geometry.Vertex (
    Vertex(..),
    vertexFromList,
    vertexToVector,
    vertexDimension,
    isValidVertexForMetric,
    assertVertexInMetric,
    distanceBetween,
    defineVertex
) where

import Algebra.Vector
import Algebra.Metric (Metric(..), distance)
import qualified Algebra.Matrix as M

newtype Vertex a = Vertex (Vector a)
  deriving (Eq, Show)

-- Constructor from list
vertexFromList :: [a] -> Vertex a
vertexFromList = Vertex . vectorFromList

-- Unwrap vector from vertex
vertexToVector :: Vertex a -> Vector a
vertexToVector (Vertex v) = v

-- Dimension of a vertex
vertexDimension :: Vertex a -> Int
vertexDimension (Vertex v) = length (vectorToList v)

isValidVertexForMetric :: Metric a -> Vertex a -> Bool
isValidVertexForMetric (Metric m) (Vertex v) =
  let (rows, cols) = M.matrixDimensions m
  in rows == length (vectorToList v) && cols == length (vectorToList v)

-- | Assert that a vertex is valid under the given metric
assertVertexInMetric :: Metric a -> Vertex a -> Maybe (Vertex a)
assertVertexInMetric metric v
  | isValidVertexForMetric metric v = Just v
  | otherwise                       = Nothing


-- | DEPRECATED: This function will be replaced by generalized measurement logic in the Metric module.
distanceBetween :: Floating a => Metric a -> Vertex a -> Vertex a -> Maybe a
distanceBetween metric v1 v2
  | isValidVertexForMetric metric v1 && isValidVertexForMetric metric v2 =
      let Vertex vec1 = v1
          Vertex vec2 = v2
      in distance metric vec1 vec2
  | otherwise = Nothing

-- Validated constructor from list
defineVertex :: Metric a -> [a] -> Maybe (Vertex a)
defineVertex metric xs =
  let v = vertexFromList xs
  in if isValidVertexForMetric metric v then Just v else Nothing
