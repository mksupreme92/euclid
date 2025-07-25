module Test.TestGeometry (testGeometry) where

import Algebra.Vector
import Algebra.Metric
import Algebra.Matrix
import Algebra.Space
import Geometry.Vertex
import Test.TestUtilities

testGeometry :: IO ()
testGeometry = do
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
