module Algebra.Transform (
    Translatable(..),
    Scalable(..),
    Rotatable(..),
    Skewable(..),
    Reflectable(..),
    Projectable(..),
    rotationMatrix
) where

import Algebra.Vector
import Algebra.Matrix

-- | Objects that support translation by a vector
class Num a => Translatable t a where
  translate :: Vector a -> t a -> Maybe (t a)

-- | Objects that support scaling by a vector (componentwise)
class Num a => Scalable s a where
  scale :: Vector a -> s a -> Maybe (s a)

-- | Objects that support rotation via an orthonormal matrix
class Num a => Rotatable r a where
  rotate :: Matrix a -> r a -> Maybe (r a)

-- | Objects that support skewing via a matrix
class Skewable t where
  skew :: Num a => Matrix a -> t a -> Maybe (t a)

-- | Objects that support reflection via a matrix
class Reflectable t where
  reflect :: Num a => Matrix a -> t a -> Maybe (t a)

-- | Objects that support projection via a matrix
class Projectable t where
  project :: Num a => Matrix a -> t a -> Maybe (t a)

-- | Vector implements all three transforms
instance Num a => Translatable Vector a where
  translate (Vector t) (Vector v)
    | length t == length v = Just $ Vector (zipWith (+) t v)
    | otherwise            = Nothing

instance Num a => Scalable Vector a where
  scale (Vector s) (Vector v)
    | length s == length v = Just $ Vector (zipWith (*) s v)
    | otherwise            = Nothing

instance Num a => Rotatable Vector a where
  rotate m (Vector v) = Vector <$> matrixVectorProduct m v

instance Skewable Vector where
  skew m (Vector v) = Vector <$> matrixVectorProduct m v

instance Reflectable Vector where
  reflect m (Vector v) = Vector <$> matrixVectorProduct m v

instance Projectable Vector where
  project m (Vector v) = Vector <$> matrixVectorProduct m v

rotationMatrix :: [[a]] -> Matrix a
rotationMatrix = matrixFromList
