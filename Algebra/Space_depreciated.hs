module Algebra.Space
  ( Space(..)
  , normAt
  , distanceIn
  , isOrthogonal
  ) where

import Algebra.Vector
import Algebra.Metric

-- | A space is defined by its dimension and a metric
data Space a = Space
  { dim    :: Int
  , metric :: Metric a
  }

-- | Computes the norm of a vector in the space
normAt :: Floating a => Space a -> Vector a -> Maybe a
normAt (Space _ m) v = do
  ip <- innerProduct m v v
  return (sqrt ip)

-- | Computes the distance between two vectors in the space
distanceIn :: Floating a => Space a -> Vector a -> Vector a -> Maybe a
distanceIn (Space _ m) v w = distance m v w

-- | Determines if two vectors are orthogonal in the space
isOrthogonal :: (Eq a, Num a) => Space a -> Vector a -> Vector a -> Bool
isOrthogonal (Space _ m) v w =
  case innerProduct m v w of
    Just 0 -> True
    _      -> False
