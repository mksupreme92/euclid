#!/usr/bin/env zsh
set -e

echo "🧩 Building Euclid (timed mode)..."
rm -rf build bin timing.log summary.log
mkdir -p bin

clang++ -std=c++20 -Wall -Wextra \
  -I. -I./dependencies/eigen-3.4.0 \
  -DEUCLID_DEBUG_TIMING_CURVE -DEUCLID_DEBUG_TIMING_INTERSECT \
  main.cpp -o bin/euclid

echo "🚀 Running Euclid with timing instrumentation..."
./bin/euclid 2> timing.log | tee output.log

echo
echo "📊 Generating timing summary..."
{
  echo "──── Timing Summary ───────────────────────────────────────────────"
  grep -E "\[CURVE_TIMING\]|\[INTERSECT_TIMING\]" timing.log | \
  awk '{sum[$2]+=$4; count[$2]++} \
  END{printf "%-30s %12s %12s %8s\n", "Label", "Avg (µs)", "Total (ms)", "Count"; \
  for(k in sum) printf "%-30s %12.3f %12.3f %8d\n", k, sum[k]/count[k], sum[k]/1000, count[k]}' | sort
  echo
  echo "──── Cache Summary ────────────────────────────────────────────────"
  grep "INTERSECT_CACHE" timing.log | \
  awk '{
    for (i=1; i<=NF; i++) {
      if ($i=="hits")    hits+=$(i+1);
      if ($i=="misses")  misses+=$(i+1);
    }
  }
  END{
    total=hits+misses;
    if (total>0)
      printf "Cache hits: %d\nCache misses: %d\nHit ratio: %.2f%%\n",
        hits, misses, 100*hits/total;
    else
      print "No cache statistics found.";
  }'
} | tee summary.log

echo "✅ Timing results saved to summary.log"
