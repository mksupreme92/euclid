module Test.RunTests where

import Test.TestAlgebra
import Test.TestGeometry

main :: IO ()
main = do
  putStrLn "🧪 Running Euclid Unit Tests..."
  testAlgebra
  testGeometry
