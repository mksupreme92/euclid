module Test.TestAlgebra (testAlgebra) where

import Algebra.Vector
import Test.TestUtilities

import Control.Monad (when)

-- Helper function for comparing vectors
approxEqual :: (Floating a, Ord a) => a -> a -> Bool
approxEqual a b = abs (a - b) < 1e-9


testAlgebra :: IO ()
testAlgebra = do
  putStrLn "\n🧮 Testing core linear algebra logic\n"
  vectorTests

vectorTests :: IO ()
vectorTests = do
  putStrLn "\n↑ Testing vector logic\n"
  let v1 = vectorFromList [1.0, 2.0, 3.0]
      v2 = vectorFromList [4.0, 5.0, 6.0]
      v3 = vectorFromList [1.0, 2.0]
  
  -- vectorAdd
  printTest "vectorAdd valid" $
    vectorAdd v1 v2 == Just (vectorFromList [5.0, 7.0, 9.0])
  printTest "vectorAdd mismatched" $
    vectorAdd v1 v3 == Nothing

  -- vectorSub
  printTest "vectorSub valid" $
    vectorSub v2 v1 == Just (vectorFromList [3.0, 3.0, 3.0])
  printTest "vectorSub mismatched" $
    vectorSub v1 v3 == Nothing

  -- scalarMul
  printTest "scalarMul 2 * v1" $
    scalarMul 2 v1 == vectorFromList [2.0, 4.0, 6.0]

  -- dotProduct
  printTest "dotProduct valid" $
    dotProduct v1 v2 == Just 32.0
  printTest "dotProduct mismatched" $
    dotProduct v1 v3 == Nothing

  -- norm
  printTest "norm of v1" $
    approxEqual (norm v1) (sqrt 14)

  -- normalize
  let nv1 = normalize v1
  printTest "normalize v1 has norm 1" $
    approxEqual (norm nv1) 1.0
    
