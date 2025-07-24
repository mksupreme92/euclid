module Test.TestUtilities (
    printTest,
) where

printTest :: String -> Bool -> IO ()
printTest name result = putStrLn $ (if result then "✅ " else "❌ ") ++ name
