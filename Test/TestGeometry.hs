-- 🜁 Euclid Geometry Kernel
-- Copyright (c) 2025 Joseph Shaw Petersen 🦅
-- Licensed under the BSD 3-Clause License. See LICENSE for details.

module Test.TestGeometry (testGeometry) where

import Algebra.Vector
import Algebra.Metric
import Algebra.Matrix
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
      euclideanMetric = Metric (identityMatrix 2)
      metric2D = Metric (identityMatrix 2)

  -- Dimension
  printTest "vertexDimension = 2" $
    vertexDimension v1 == 2

  -- isValidVertex
  printTest "isValidVertex (metric2D)" $
    isValidVertexForMetric euclideanMetric v1

  printTest "assertVertexInMetric" $
    assertVertexInMetric euclideanMetric v1 == Just v1


  printTest "defineVertex" $
    defineVertex euclideanMetric [1, 2] == Just v1

  printTest "defineVertex invalid" $
    defineVertex euclideanMetric [1, 2, 3] == Nothing

  -- translateVertex
  printTest "translate Vertex in 2D" $
    fromJust (translate euclideanMetric offset v1) == v2

  -- Rotate vertex about center (2D)
  let v2D      :: Vertex Double
      v2D      = vertexFromList [1.0, 0.0]
      center2D = vertexFromList [0.0, 0.0]
      rot90_2D = [[0, -1], [1, 0]] :: Matrix Double

  printTest "rotateAbout 2D (90° around origin)" $
    rotateAbout metric2D center2D rot90_2D v2D == Just (vertexFromList [0.0, 1.0])

  -- Rotate vertex about center (3D)
  let v3D = vertexFromList [1.0, 0.0, 0.0]
      center3D = vertexFromList [0.0, 0.0, 0.0]
      rotZ_3D = matrixFromList
        [[0.0, -1.0, 0.0],
         [1.0,  0.0, 0.0],
         [0.0,  0.0, 1.0]]

  printTest "rotateAbout 3D (Z-axis 90°)" $
    rotateAbout metric2D center3D rotZ_3D v3D == Just (vertexFromList [0.0, 1.0, 0.0])

  -- Rotate vertex about center (5D)
  let v5D = vertexFromList [1.0, 0.0, 0.0, 0.0, 0.0]
      center5D = vertexFromList [0.0, 0.0, 0.0, 0.0, 0.0]
      rotXY_5D = [[0, -1, 0, 0, 0],
                  [1,  0, 0, 0, 0],
                  [0,  0, 1, 0, 0],
                  [0,  0, 0, 1, 0],
                  [0,  0, 0, 0, 1]] :: Matrix Double
      metric5D = Metric (identityMatrix 5)

  printTest "rotateAbout 5D (XY-plane 90°)" $
    rotateAbout metric5D center5D rotXY_5D v5D == Just (vertexFromList [0.0, 1.0, 0.0, 0.0, 0.0])


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
      euclideanMetric = Metric (identityMatrix 2)
      metric2D = Metric (identityMatrix 2)

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
        result = translate euclideanMetric offsetT (InfiniteLine pT dT)
    in case result of
        Just (InfiniteLine p' d') ->
          approxEqualVector (vectorFromVertex p') expectedP &&
          approxEqualVector d' expectedD

  -- printTest "scale Edge (Ray)" $
  --   fmap (\(Ray v d) -> (v, d)) (scale euclideanMetric (vectorFromList [2.0, 2.0]) =<< rayFromVertexAndVector p1 dir)
  --     == Just (p1, vectorFromList [6.0, 8.0])

  do
    let seg = Segment v2D center2D
        Just rotated = rotateAbout metric2D center2D rot90_2D seg
    printTest "rotateAbout Edge (Segment 90°)" $
      rotated == Segment
        (fromJust $ rotateAbout metric2D center2D rot90_2D v2D)
        (fromJust $ rotateAbout metric2D center2D rot90_2D center2D)

  -- Test rotateAbout for InfiniteLine: apply 90 degree rotation about origin
  printTest "rotateAbout Edge (InfiniteLine 90°)" $
    case rotateAbout metric2D centerR rot90_2D (InfiniteLine pR dR) of
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
