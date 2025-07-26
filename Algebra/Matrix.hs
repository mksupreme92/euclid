module Algebra.Matrix
  ( Matrix
  , matrixFromList
  , identityMatrix
  , matrixVectorProduct
  , transpose
  , matrixMultiply
  , matrixDimensions
  ) where

type Matrix a = [[a]]

matrixFromList :: [[a]] -> Matrix a
matrixFromList = id

identityMatrix :: Num a => Int -> Matrix a
identityMatrix n =
  [ [ fromIntegral (fromEnum (i == j)) | j <- [0 .. n - 1] ] | i <- [0 .. n - 1] ]

-- | Safe matrix-vector product: m · v
-- Returns Nothing if dimensions mismatch
matrixVectorProduct :: Num a => Matrix a -> [a] -> Maybe [a]
matrixVectorProduct m v =
  if all (\row -> length row == length v) m
     then Just $ map (sum . zipWith (*) v) m
     else Nothing

transpose :: Matrix a -> Matrix a
transpose ([]:_) = []
transpose x = map head x : transpose (map tail x)

matrixMultiply :: Num a => Matrix a -> Matrix a -> Maybe (Matrix a)
matrixMultiply a b =
  let bT = transpose b
  in if all (\row -> length row == length bT) a
       then Just [ [ sum $ zipWith (*) row col | col <- bT ] | row <- a ]
       else Nothing

matrixDimensions :: Matrix a -> (Int, Int)
matrixDimensions m = (length m, if null m then 0 else length (head m))
