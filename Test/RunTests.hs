module Test.RunTests where

import Test.TestAlgebra

main :: IO ()
main = do
  putStrLn "🧪 Running Euclid Unit Tests..."
  testAlgebra
