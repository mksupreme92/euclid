module Algebra.Metric
  ( Metric(..)
  , innerProduct
  , distance
  ) where

import Algebra.Vector
import Algebra.Matrix

-- | A metric is just a constant matrix representing the inner product
newtype Metric a = Metric (Matrix a)

-- | Computes inner product ⟨v, w⟩ using the metric
innerProduct :: Num a => Metric a -> Vector a -> Vector a -> Maybe a
innerProduct (Metric g) v w = do
  gvList <- matrixVectorProduct g (vectorToList v)
  let gv = vectorFromList gvList
  dotProduct gv w

-- | Computes distance using the inner product (norm of v - w)
distance :: Floating a => Metric a -> Vector a -> Vector a -> Maybe a
distance metric v w = do
  delta <- vectorSub v w
  ip <- innerProduct metric delta delta
  return (sqrt ip)
