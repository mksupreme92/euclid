-- Geometry.Vertex
-- Euclid Geometry Kernel — Vertex Definition
-- A vertex is a point in space; it wraps a vector and supports transforms.

module Geometry.Vertex (
    Vertex(..),
    vertexFromList,
    vertexToVector,
    vertexDimension,
    isValidVertex,
    assertVertexInSpace,
    distanceBetween,
    defineVertex
) where

import Algebra.Vector
import Algebra.Space

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

-- Validity check against a space
isValidVertex :: Space a -> Vertex a -> Bool
isValidVertex space (Vertex v) = length (vectorToList v) == dim space

-- Validate a vertex is in the space
assertVertexInSpace :: Space a -> Vertex a -> Maybe (Vertex a)
assertVertexInSpace space v = if isValidVertex space v then Just v else Nothing



-- Distance between vertices in a space
distanceBetween :: Floating a => Space a -> Vertex a -> Vertex a -> Maybe a
distanceBetween space (Vertex v1) (Vertex v2) = distanceIn space v1 v2

-- Validated constructor from list
defineVertex :: Space a -> [a] -> Maybe (Vertex a)
defineVertex space xs =
  let v = vertexFromList xs
  in if isValidVertex space v then Just v else Nothing
