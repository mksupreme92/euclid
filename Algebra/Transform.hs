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
import Algebra.Metric
import Geometry.Vertex
import Geometry.Edge

-- | Objects that support translation by a vector
class Num a => Translatable t a where
  translate :: Metric a -> Vector a -> t a -> Maybe (t a)

-- | Objects that support scaling by a vector (componentwise)
class Num a => Scalable s a where
  scale :: Metric a -> Vector a -> s a -> Maybe (s a)

-- | Objects that support rotation via an orthonormal matrix
class Num a => Rotatable r a where
  rotate :: Metric a -> Matrix a -> r a -> Maybe (r a)

-- | Objects that support skewing via a matrix
class Num a => Skewable t a where
  skew :: Metric a -> Matrix a -> t a -> Maybe (t a)

-- | Objects that support reflection via a matrix
class Num a => Reflectable t a where
  reflect :: Metric a -> Matrix a -> t a -> Maybe (t a)

-- | Objects that support projection via a matrix
class Num a => Projectable t a where
  project :: Metric a -> Matrix a -> t a -> Maybe (t a)

instance (Num a, Eq a) => Translatable Vector a where
  translate metric (Vector t) (Vector v)
    | isEuclidean metric && length t == length v = Just $ Vector (zipWith (+) t v)
    | isEuclidean metric = Nothing
    | otherwise = Nothing

instance (Num a, Eq a) => Translatable Vertex a where
  translate metric offset (Vertex v) = fmap Vertex (translate metric offset v)

instance (Num a, Eq a) => Translatable Edge a where
  translate metric offset (InfiniteLine p d) = InfiniteLine <$> translate metric offset p <*> Just d
  translate metric offset (Ray p d)          = Ray <$> translate metric offset p <*> Just d
  translate metric offset (Segment p1 p2)    = Segment <$> translate metric offset p1 <*> translate metric offset p2

instance (Num a, Eq a) => Scalable Vector a where
  scale metric (Vector s) (Vector v)
    | isEuclidean metric && length s == length v = Just $ Vector (zipWith (*) s v)
    | isEuclidean metric = Nothing
    | otherwise = Nothing


instance (Num a, Eq a) => Rotatable Vector a where
  rotate metric m (Vector v)
    | isEuclidean metric = Vector <$> matrixVectorProduct m v
    | otherwise = Nothing


instance Num a => Skewable Vector a where
  skew _ m (Vector v) = Vector <$> matrixVectorProduct m v


instance Num a => Reflectable Vector a where
  reflect _ m (Vector v) = Vector <$> matrixVectorProduct m v


instance Num a => Projectable Vector a where
  project _ m (Vector v) = Vector <$> matrixVectorProduct m v


rotationMatrix :: [[a]] -> Matrix a
rotationMatrix = matrixFromList


-- | Objects that support rotation about a reference vertex
class Num a => RotatableAboutVertex g a where
  rotateAbout :: Metric a -> Vertex a -> Matrix a -> g a -> Maybe (g a)
-- | Instance: RotatableAboutVertex for Vertex
instance (Num a, Eq a) => RotatableAboutVertex Vertex a where
  rotateAbout metric (Vertex center) rot (Vertex v)
    | isEuclidean metric = do
        delta <- vectorSub v center
        rotatedDeltaList <- matrixVectorProduct rot (vectorToList delta)
        let rotatedDelta = vectorFromList rotatedDeltaList
        result <- vectorAdd rotatedDelta center
        return (Vertex result)
    | otherwise = Nothing

instance (Num a, Eq a) => RotatableAboutVertex Edge a where
  rotateAbout metric c m (InfiniteLine p d) = InfiniteLine <$> rotateAbout metric c m p <*> rotate metric m d
  rotateAbout metric c m (Ray p d)          = Ray <$> rotateAbout metric c m p <*> rotate metric m d
  rotateAbout metric c m (Segment p1 p2)    = Segment <$> rotateAbout metric c m p1 <*> rotateAbout metric c m p2
