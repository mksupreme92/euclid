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



  -- distanceBetween
  printTest "distanceBetween" $
    case distanceBetween space2D v1 v2 of
      Just d -> approxEqual d 5.0
      Nothing -> False

  -- translateVertex
  printTest "translate Vertex in 2D" $
    fromJust (translate offset v1) == v2

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

  -- Updated test: translate InfiniteLine, check base point and direction vector
  printTest "translate Edge (InfiniteLine)" $
    let pT = Vertex (Vector [1.0, 2.0])
        dT = Vector [1.0, 0.0]
        offsetT = Vector [0.0, 1.0]
        expectedP = Vector [1.0, 3.0]
        expectedD = dT
        result = translate offsetT (InfiniteLine pT dT)
    in case result of
        Just (InfiniteLine p' d') ->
          approxEqualVector (vectorFromVertex p') expectedP &&
          approxEqualVector d' expectedD

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

  -- Test rotateAbout for InfiniteLine: apply 90 degree rotation about origin
  printTest "rotateAbout Edge (InfiniteLine 90°)" $
    case rotateAbout centerR rot90_2D (InfiniteLine pR dR) of
      Just (InfiniteLine p' d') ->
        approxEqualVector (vectorFromVertex p') expectedP &&
        approxEqualVector d' expectedD
      _ -> False
    where
      pR = vertexFromList [1.0, 0.0]
      dR = vectorFromList [1.0, 0.0]
      centerR = vertexFromList [0.0, 0.0]
      expectedP = vectorFromList [0.0, 1.0]
      expectedD = vectorFromList [0.0, 1.0]
