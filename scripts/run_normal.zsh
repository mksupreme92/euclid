#!/usr/bin/env zsh
set -e

echo "🧩 Building Euclid (normal mode)..."
rm -rf build bin
mkdir -p bin

clang++ -std=c++20 -Wall -Wextra \
  -I. -I./dependencies/eigen-3.4.0 \
  main.cpp -o bin/euclid

echo "🚀 Running Euclid..."
./bin/euclid
