module Test.TestGeometry (testGeometry) where

import Algebra.Vector
import Algebra.Metric
import Algebra.Matrix
import Algebra.Space
import Algebra.Transform
import Geometry.Vertex
import Geometry.Edge
import Test.TestUtilities
import Data.Maybe (fromJust)

testGeometry :: IO ()
testGeometry = do
    vertexTests
    edgeTests


vertexTests :: IO ()
vertexTests = do
  putStrLn "\n🧭 Testing Geometry.Vertex Module\n"

  let v1 = vertexFromList [1.0, 2.0]
      v2 = vertexFromList [4.0, 6.0]
      offset = vectorFromList [3.0, 4.0]
      space2D = Space 2 (Metric $ matrixFromList [[1,0],[0,1]])

  -- Dimension
  printTest "vertexDimension = 2" $
    vertexDimension v1 == 2

  -- isValidVertex
  printTest "isValidVertex (space2D)" $
    isValidVertex space2D v1

  -- assertVertexInSpace
  printTest "assertVertexInSpace valid" $
    assertVertexInSpace space2D v1 == Just v1

  -- translateVertex
  printTest "translateVertex +offset" $
    translateVertex space2D v1 offset == Just v2

  -- distanceBetween
  printTest "distanceBetween" $
    case distanceBetween space2D v1 v2 of
      Just d -> approxEqual d 5.0
      Nothing -> False

  -- Rotate vertex about center (2D)
  let v2D      :: Vertex Double
      v2D      = vertexFromList [1.0, 0.0]
      center2D = vertexFromList [0.0, 0.0]
      rot90_2D = [[0, -1], [1, 0]] :: Matrix Double

  printTest "rotateAbout 2D (90° around origin)" $
    rotateAbout center2D rot90_2D v2D == Just (vertexFromList [0.0, 1.0])

  -- Rotate vertex about center (3D)
  let v3D = vertexFromList [1.0, 0.0, 0.0]
      center3D = vertexFromList [0.0, 0.0, 0.0]
      rotZ_3D = matrixFromList
        [[0.0, -1.0, 0.0],
         [1.0,  0.0, 0.0],
         [0.0,  0.0, 1.0]]

  printTest "rotateAbout 3D (Z-axis 90°)" $
    rotateAbout center3D rotZ_3D v3D == Just (vertexFromList [0.0, 1.0, 0.0])

  -- Rotate vertex about center (5D)
  let v5D = vertexFromList [1.0, 0.0, 0.0, 0.0, 0.0]
      center5D = vertexFromList [0.0, 0.0, 0.0, 0.0, 0.0]
      rotXY_5D = [[0, -1, 0, 0, 0],
                  [1,  0, 0, 0, 0],
                  [0,  0, 1, 0, 0],
                  [0,  0, 0, 1, 0],
                  [0,  0, 0, 0, 1]] :: Matrix Double

  printTest "rotateAbout 5D (XY-plane 90°)" $
    rotateAbout center5D rotXY_5D v5D == Just (vertexFromList [0.0, 1.0, 0.0, 0.0, 0.0])


edgeTests :: IO ()
edgeTests = do
  -- 🧵 Edge Tests
  putStrLn "\n📏 Testing Geometry.Edge Module\n"

  let p1     = vertexFromList [1.0, 2.0] :: Vertex Double
      p2     = vertexFromList [4.0, 6.0] :: Vertex Double
      dir    = vectorFromList [3.0, 4.0] :: Vector Double
      offset = vectorFromList [3.0, 4.0] :: Vector Double
      v2D      = vertexFromList [0.0, 1.0] :: Vertex Double
      center2D = vertexFromList [0.0, 0.0] :: Vertex Double

  -- Edge Constructors
  printTest "lineFromVertexAndVector valid" $
    lineFromVertexAndVector p1 dir == Just (InfiniteLine p1 dir)

  printTest "rayFromVertexAndVector valid" $
    rayFromVertexAndVector p1 dir == Just (Ray p1 dir)

  printTest "segmentFromVertices valid" $
    segmentFromVertices p1 p2 == Just (Segment p1 p2)

  -- Edge Transformations
  let rot90_2D = matrixFromList [[0, -1], [1, 0]] :: Matrix Double

  printTest "translate Edge (InfiniteLine)" $
    fmap (\(InfiniteLine v d) -> (v, d)) (translate offset =<< lineFromVertexAndVector p1 dir)
      == Just (fromJust (translate offset p1), dir)

  -- printTest "scale Edge (Ray)" $
  --   fmap (\(Ray v d) -> (v, d)) (scale (vectorFromList [2.0, 2.0]) =<< rayFromVertexAndVector p1 dir)
  --     == Just (p1, vectorFromList [6.0, 8.0])

  do
    let seg = Segment v2D center2D
        Just rotated = rotateAbout center2D rot90_2D seg
    printTest "rotateAbout Edge (Segment 90°)" $
      rotated == Segment
        (fromJust $ rotateAbout center2D rot90_2D v2D)
        (fromJust $ rotateAbout center2D rot90_2D center2D)

  printTest "rotateAbout Edge (InfiniteLine)" $
    fmap (\(InfiniteLine v d) -> (v, d)) (rotateAbout center2D rot90_2D =<< lineFromVertexAndVector v2D (vectorFromList [1.0, 0.0]))
      == Just (vertexFromList [0.0, 1.0], vectorFromList [0.0, 1.0])
