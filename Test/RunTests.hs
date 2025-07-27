-- 🜁 Euclid Geometry Kernel
-- Copyright (c) 2025 Joseph Shaw Petersen 🦅
-- Licensed under the BSD 3-Clause License. See LICENSE for details.

module Test.RunTests where

import Test.TestAlgebra
import Test.TestGeometry

main :: IO ()
main = do
  putStrLn "🧪 Running Euclid Unit Tests..."
  testAlgebra
  testGeometry
