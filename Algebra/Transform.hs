-- 🜁 Euclid Geometry Kernel
-- Copyright (c) 2025 Joseph Shaw Petersen 🦅
-- Licensed under the BSD 3-Clause License. See LICENSE for details.


-- NOTE (2025-07-28): Metric-specific transformation routing is partially implemented.
--   Functions like `rotate`, `skew`, `reflect`, and `project` now validate that the
--   transformation matrix preserves the metric using `isMetricPreserving`.
--   However, some geometric types still return `Nothing` when metric ≠ Euclidean
--   (e.g. `rotateAbout` for Vertex/Edge). These will need future handling.


module Algebra.Transform (
    Translatable(..),
    Scalable(..),
    Rotatable(..),
    RotatableAboutVertex(..),
    Skewable(..),
    Reflectable(..),
    Projectable(..),
    rotationMatrix
) where

import Algebra.Vector
import Algebra.Matrix
import Algebra.Metric
import Geometry.Vertex
import Geometry.Edge


-- | Instances that support translation by a vector ------------------
class Num a => Translatable t a where
  translate :: Metric a -> Vector a -> t a -> Maybe (t a)
  
instance (Num a, Eq a) => Translatable Vector a where
  translate metric (Vector t) (Vector v)
    | isEuclidean metric && length t == length v = Just $ Vector (zipWith (+) t v)
    | isEuclidean metric = Nothing
    | otherwise = Just $ Vector (zipWith (+) t v) -- Metric-aware translation stub (linear shift still valid in flat space)

instance (Num a, Eq a) => Translatable Vertex a where
  translate metric offset (Vertex v) = fmap Vertex (translate metric offset v)

instance (Num a, Eq a) => Translatable Edge a where
  translate metric offset (InfiniteLine p d) = InfiniteLine <$> translate metric offset p <*> Just d
  translate metric offset (Ray p d)          = Ray <$> translate metric offset p <*> Just d
  translate metric offset (Segment p1 p2)    = Segment <$> translate metric offset p1 <*> translate metric offset p2
  
-- To do:
    -- Translatable Curve
    -- Translatable Plane
    -- Translatable Face
    -- Translatable Surface
    -- Translatable Volume

---------------------------------------------------------------------

-- | Instances that support scaling by a vector (componentwise) -----
class Num a => Scalable s a where
  scale :: Metric a -> Vector a -> s a -> Maybe (s a)

instance (Num a, Eq a) => Scalable Vector a where
  scale metric (Vector s) (Vector v)
    | isEuclidean metric && length s == length v =
        Just $ Vector (zipWith (*) s v)
    | isEuclidean metric = Nothing
    | otherwise =
        let scaleMatrix = diag s
        in if isMetricPreserving metric scaleMatrix
              then Vector <$> matrixVectorProduct scaleMatrix v
              else Nothing
-- To do:
    -- Scalable Segment (instance of Edge)
    -- Scalable Curve
    -- Scalable Plane
    -- Scalable Face
    -- Scalable Surface
    -- Scalable Volume

---------------------------------------------------------------------

-- | Instances that support rotation via an orthonormal matrix ------
class Num a => Rotatable r a where
  rotate :: Metric a -> Matrix a -> r a -> Maybe (r a)
  
instance (Num a, Eq a) => Rotatable Vector a where
  rotate metric m (Vector v)
    | isMetricPreserving metric m = Vector <$> matrixVectorProduct m v
    | otherwise = Nothing
    
---------------------------------------------------------------------

-- | Instances that support rotation about a reference vertex -------
class Num a => RotatableAboutVertex g a where
  rotateAbout :: Metric a -> Vertex a -> Matrix a -> g a -> Maybe (g a)

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
  
-- To do:
    -- RotatableAboutVertex Curve
    -- RotatableAboutVertex Plane
    -- RotatableAboutVertex Face
    -- RotatableAboutVertex Surface
    -- RotatableAboutVertex Volume

  
---------------------------------------------------------------------

-- | Instances that support skewing via a matrix --------------------
class Num a => Skewable t a where
  skew :: Metric a -> Matrix a -> t a -> Maybe (t a)
  
instance (Num a, Eq a) => Skewable Vector a where
  skew metric m (Vector v)
    | isMetricPreserving metric m = Vector <$> matrixVectorProduct m v
    | otherwise = Nothing
    
-- To do:
    -- Skweable Face
    -- revist whether eforcing metric preservation is appropriate for skews
    
---------------------------------------------------------------------

-- | Instances that support reflection via a matrix -----------------
class Num a => Reflectable t a where
  reflect :: Metric a -> Matrix a -> t a -> Maybe (t a)

instance (Num a, Eq a) => Reflectable Vector a where
  reflect metric m (Vector v)
    | isMetricPreserving metric m = Vector <$> matrixVectorProduct m v
    | otherwise = Nothing

-- To do:
    -- ReflectableAboutPlane class for Geometry instances with defined locations in metric space (Vertex, Edge, Plane, Face, Surface, Volume)

---------------------------------------------------------------------

-- | Instances that support projection via a matrix -----------------
class Num a => Projectable t a where
  project :: Metric a -> Matrix a -> t a -> Maybe (t a)

instance (Num a, Eq a) => Projectable Vector a where
  project metric m (Vector v)
    | isMetricPreserving metric m = Vector <$> matrixVectorProduct m v
    | otherwise = Nothing

-- To do:
    -- ProjectAbleOntoPlane class for Geometry instances defined locations in metric space (Vertex, Edge, Plane, Face, Surface, Volume)

---------------------------------------------------------------------


rotationMatrix :: [[a]] -> Matrix a
rotationMatrix = matrixFromList
