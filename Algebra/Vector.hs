module Algebra.Vector (
    Vector(..),
    vectorFromList,
    vectorToList,
    vectorAdd,
    vectorSub,
    scalarMul,
    dotProduct,
    norm,
    normalize
) where

import Data.List (zipWith)

newtype Vector a = Vector [a]
  deriving (Eq, Show)

-- Conversion helpers
vectorFromList :: [a] -> Vector a
vectorFromList = Vector

vectorToList :: Vector a -> [a]
vectorToList (Vector xs) = xs

-- Basic operations
vectorAdd :: Num a => Vector a -> Vector a -> Maybe (Vector a)
vectorAdd (Vector u) (Vector v)
  | length u == length v = Just $ Vector (zipWith (+) u v)
  | otherwise            = Nothing

vectorSub :: Num a => Vector a -> Vector a -> Maybe (Vector a)
vectorSub (Vector u) (Vector v)
  | length u == length v = Just $ Vector (zipWith (-) u v)
  | otherwise            = Nothing

scalarMul :: Num a => a -> Vector a -> Vector a
scalarMul k (Vector xs) = Vector (map (k *) xs)

dotProduct :: Num a => Vector a -> Vector a -> Maybe a
dotProduct (Vector u) (Vector v)
  | length u == length v = Just $ sum $ zipWith (*) u v
  | otherwise            = Nothing

norm :: Floating a => Vector a -> a
norm v = sqrt $ maybe 0 id (dotProduct v v)

normalize :: Floating a => Vector a -> Vector a
normalize v = scalarMul (1 / norm v) v
