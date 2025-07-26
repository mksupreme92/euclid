module Algebra.Transform (
    Translatable(..),
    Scalable(..),
    Rotatable(..),
    Skewable(..),
    Reflectable(..),
    Projectable(..),
    RotatableAboutVertex(..),
    rotationMatrix
) where

import Algebra.Vector
import Algebra.Matrix
import Geometry.Vertex
import Geometry.Edge

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

instance Num a => Translatable Vector a where
  translate (Vector t) (Vector v)
    | length t == length v = Just $ Vector (zipWith (+) t v)
    | otherwise            = Nothing

instance Num a => Translatable Vertex a where
  translate offset (Vertex v) = fmap Vertex (translate offset v)

instance Num a => Translatable Edge a where
  translate offset (InfiniteLine p d) = InfiniteLine <$> translate offset p <*> translate offset d
  translate offset (Ray p d)          = Ray <$> translate offset p <*> translate offset d
  translate offset (Segment p1 p2)    = Segment <$> translate offset p1 <*> translate offset p2

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


-- | Objects that support rotation about a reference vertex
class Num a => RotatableAboutVertex g a where
  rotateAbout :: Vertex a -> Matrix a -> g a -> Maybe (g a)
-- | Instance: RotatableAboutVertex for Vertex
instance Num a => RotatableAboutVertex Vertex a where
  rotateAbout (Vertex center) rot (Vertex v) = do
    delta <- vectorSub v center
    rotatedDeltaList <- matrixVectorProduct rot (vectorToList delta)
    let rotatedDelta = vectorFromList rotatedDeltaList
    result <- vectorAdd rotatedDelta center
    return (Vertex result)

instance Num a => RotatableAboutVertex Edge a where
  rotateAbout c m (InfiniteLine p d) = InfiniteLine <$> rotateAbout c m p <*> rotate m d
  rotateAbout c m (Ray p d)          = Ray <$> rotateAbout c m p <*> rotate m d
  rotateAbout c m (Segment p1 p2)    = Segment <$> rotateAbout c m p1 <*> rotateAbout c m p2
