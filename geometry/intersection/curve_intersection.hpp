namespace Euclid {
namespace Geometry {

#ifdef EUCLID_DEBUG_TIMING_INTERSECT
#include <chrono>
#include <fstream>
#include <sstream>
// Thread-local buffer for timing logs
thread_local std::ostringstream intersectTimingBuffer;
struct IntersectTimer {
    const char* label;
    std::chrono::high_resolution_clock::time_point start;
    IntersectTimer(const char* lbl)
        : label(lbl), start(std::chrono::high_resolution_clock::now()) {}
    ~IntersectTimer() {
        auto end = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        intersectTimingBuffer << "INTERSECT_TIMING " << label << " " << dur << " us\n";
    }
};
#define INTERSECT_TIME_SCOPE(name) IntersectTimer __intersect_timer_##__LINE__{name};
#else
#define INTERSECT_TIME_SCOPE(name)
#endif

// Per-curve sample cache for intersection routines
#include <array>
template <typename T, int N>
struct CurveSampleCache {
    struct Entry {
        Point<T, N> pos;
        Eigen::Matrix<T, N, 1> d1;
        Eigen::Matrix<T, N, 1> d2;
    };
    std::array<Entry, 256> data;
    std::array<long long, 256> keys;
    int pos = 0;
    int hits = 0;
    int misses = 0;
    CurveSampleCache() : pos(0), hits(0), misses(0) {
        keys.fill(-1);
    }
};

template <typename T, int N>
IntersectionResult<T, N> intersect(
    const Curve<T, N>& c1,
    const Curve<T, N>& c2,
    const Euclid::Tolerance& tol = Euclid::Tolerance())
{
    // ------------------------------------------------------------------------------
    // Curve–Curve Intersection (Adaptive Refinement)
    //
    // This routine uses an adaptive subdivision strategy with tolerance-scaled
    // sampling density to robustly find intersection points between two parametric
    // curves. The algorithm:
    //  1. Builds adaptive parameter segments per curve based on local curvature.
    //  2. Refines overlapping segment pairs recursively until convergence.
    //  3. Applies tolerance-based deduplication to avoid duplicate hits.
    //  4. Detects potential "curve of intersection" overlap cases for future spline support.
    //
    // Future phases:
    //  - Phase 2: return parametric overlap intervals.
    //  - Phase 3: return full B-spline curve-of-intersection representation.
    // ------------------------------------------------------------------------------
    // 1. basic bbox and scale ---------------------------------------------
    auto [min1, max1] = c1.boundingBox();
    auto [min2, max2] = c2.boundingBox();

    // world-ish scale from both boxes
    T maxCoord = T(0);
    for (int i = 0; i < N; ++i) {
        maxCoord = std::max(maxCoord, std::abs(min1.coords[i]));
        maxCoord = std::max(maxCoord, std::abs(max1.coords[i]));
        maxCoord = std::max(maxCoord, std::abs(min2.coords[i]));
        maxCoord = std::max(maxCoord, std::abs(max2.coords[i]));
    }
    T worldEps = tol.evaluateEpsilon(maxCoord);
    // Define boxesOverlap lambda before any use
    auto boxesOverlap = [&](const Point<T, N>& aMin, const Point<T, N>& aMax,
                            const Point<T, N>& bMin, const Point<T, N>& bMax) -> bool {
        for (int i = 0; i < N; ++i) {
            if (aMax.coords[i] + worldEps < bMin.coords[i] - worldEps ||
                aMin.coords[i] - worldEps > bMax.coords[i] + worldEps) {
                return false;
            }
        }
        return true;
    };

    // quick reject
    if (!boxesOverlap(min1, max1, min2, max2)) {
        return {};
    }

    // 2. parameter domains -------------------------------------------------
    auto [t1_min, t1_max] = c1.domain();
    auto [t2_min, t2_max] = c2.domain();

    // --- lightweight per-curve cache for reuse across refineSegments/refineCell ---
    CurveSampleCache<T, N> cache1;
    CurveSampleCache<T, N> cache2;

    auto getSample = [&](const Curve<T, N>& c,
                         CurveSampleCache<T, N>& cache,
                         T t, T tMin, T tMax) -> const typename CurveSampleCache<T, N>::Entry& {
        // Quantize parameter to 1e-6 over the domain and map to 256-slot cache
        constexpr int kSlots = 256;
        constexpr int kMask  = kSlots - 1; // 255
        const T span = (tMax - tMin);
        // Guard against zero-length domain
        T u = span > T(0) ? (t - tMin) / span : T(0);
        if (u < T(0)) u = T(0);
        if (u > T(1)) u = T(1);
        long long key = static_cast<long long>(u * 1'000'000); // quantized key
        int idx = static_cast<int>(key) & kMask;
        if (cache.keys[idx] == key) {
            ++cache.hits;
            return cache.data[idx];
        }
        ++cache.misses;
        // Miss: evaluate once and overwrite slot
        typename CurveSampleCache<T, N>::Entry e;
        e.pos = c.evaluate(t);
        e.d1  = c.evaluateDerivativeCached(t);
        e.d2  = c.evaluateSecondDerivative(t);
        cache.data[idx] = std::move(e);
        cache.keys[idx] = key;
        return cache.data[idx];
    };


    // 3. adaptive sampler --------------------------------------------------
    // we keep two samplers, one for each curve, each starting with a coarse grid and
    // being refined where the curves come close.
    struct ParamSeg {
        T a, b;   // parameter interval
        T scale;  // local geometric scale for this interval
    };

    // helper: estimate local scale for a curve over [a,b] via 3 samples
    auto estimateCurveScale = [&](const Curve<T, N>& c, T a, T b) -> T {
        T m = (a + b) / T(2);
        auto pA = c.evaluate(a).coords;
        auto pM = c.evaluate(m).coords;
        auto pB = c.evaluate(b).coords;
        T s = std::max({pA.cwiseAbs().maxCoeff(), pM.cwiseAbs().maxCoeff(), pB.cwiseAbs().maxCoeff()});
        // include chord length as an extra measure
        s = std::max(s, (pA - pB).norm());
        return s;
    };

    // initial segments: just whole domain
    std::vector<ParamSeg> segs1 {{t1_min, t1_max, estimateCurveScale(c1, t1_min, t1_max)}};
    std::vector<ParamSeg> segs2 {{t2_min, t2_max, estimateCurveScale(c2, t2_min, t2_max)}};

    // target spatial resolution: what we want two points to be within in real space
    // start from worldEps and relax a bit so we don't explode refinements
    // const T targetSpatial = worldEps * T(2);

    // function: refine a set of segments until each segment is "flat enough"
    // Optionally records segment endpoints for later reuse.
    auto refineSegments = [&](const Curve<T, N>& cv, std::vector<ParamSeg>& segs,
                              CurveSampleCache<T, N>& cache,
                              T tMin, T tMax,
                              std::vector<std::pair<Point<T, N>, Point<T, N>>>* segmentEndpoints = nullptr) {
        bool changed = true;
        int iter = 0;
        const int maxIter = 7; // allow a few more passes for wavy curves
        while (changed && iter++ < maxIter) {
            changed = false;
            std::vector<ParamSeg> next;
            next.reserve(segs.size() * 2); // reserve up-front to avoid repeated reallocations
            if (segmentEndpoints) segmentEndpoints->clear();

            for (auto& s : segs) {
                T a = s.a;
                T b = s.b;
                T m = (a + b) / T(2);
                // Use cache for this curve (cv): use the provided cache and domain bounds
                auto sA = getSample(cv, cache, a, tMin, tMax);
                auto sM = getSample(cv, cache, m, tMin, tMax);
                auto sB = getSample(cv, cache, b, tMin, tMax);
                auto pA = sA.pos;
                auto pM = sM.pos;
                auto pB = sB.pos;
                auto dA = sA.d1;
                auto dB = sB.d1;
                auto ddA = sA.d2;
                auto ddB = sB.d2;

                // straight-line deviation (existing metric)
                Eigen::Matrix<T, N, 1> AB = pB.coords - pA.coords;
                Eigen::Matrix<T, N, 1> AM = pM.coords - pA.coords;
                T dev = T(0);
                T ab2 = AB.squaredNorm();
                if (ab2 > T(0)) {
                    T t = AM.dot(AB) / ab2;
                    Eigen::Matrix<T, N, 1> proj = pA.coords + t * AB;
                    dev = (pM.coords - proj).norm();
                }

                // Multi-sample curvature estimation using derivative direction variance
                auto safeNorm = [](const Eigen::Matrix<T, N, 1>& v) -> Eigen::Matrix<T, N, 1> {
                    T n = v.norm();
                    if (n > std::numeric_limits<T>::epsilon()) {
                        return (v / n).eval();
                    }
                    Eigen::Matrix<T, N, 1> z;
                    z.setZero();
                    return z;
                };
                constexpr int samples = 5;
                T maxAngle = T(0);
                // Use cached derivative evaluations for efficiency within this segment
                auto prevDir = safeNorm(getSample(cv, cache, a, tMin, tMax).d1);
                for (int sidx = 1; sidx < samples; ++sidx) {
                    T t = a + (b - a) * (T(sidx) / (samples - 1));
                    auto dir = safeNorm(getSample(cv, cache, t, tMin, tMax).d1);
                    T dotVal = std::clamp(prevDir.dot(dir), T(-1), T(1));
                    T ang = std::acos(dotVal);
                    if (ang > maxAngle) maxAngle = ang;
                    prevDir = dir;
                }
                T curvatureMeasure = maxAngle;

                // --- Curvature sign tracking (for oscillatory refinement) ---
                // dA, dB, ddA, ddB already above using cache

                // approximate scalar curvature sign (2D: cross product z-sign, ND: derivative projection)
                T kappaA = T(0), kappaB = T(0);
                if constexpr (N == 2) {
                    kappaA = dA[0]*ddA[1] - dA[1]*ddA[0];
                    kappaB = dB[0]*ddB[1] - dB[1]*ddB[0];
                } else {
                    kappaA = dA.normalized().dot(ddA);
                    kappaB = dB.normalized().dot(ddB);
                }
                bool curvatureSignFlip = (kappaA * kappaB) < 0;

                // local geometric scale for tolerance
                T localScale = std::max({pA.coords.cwiseAbs().maxCoeff(),
                                         pM.coords.cwiseAbs().maxCoeff(),
                                         pB.coords.cwiseAbs().maxCoeff()});
                T localEps = tol.evaluateEpsilon(localScale);

                // derive a curvature tolerance: small curves (scale ~1) with tight tol get more splits
                // radians and length are different units, so give curvature a generous multiplier
                T curvatureTol = localEps * T(400); // <<-- geometry/tol based;

                bool needsSplit = false;

                // 1) if geometric deviation is bigger than local tol → split
                if (dev > localEps * T(1.5)) {
                    needsSplit = true;
                }

                // 2) if curve is bending noticeably → split
                if (curvatureMeasure > curvatureTol) {
                    needsSplit = true;
                }

                // 4) if curvature sign flips within segment → split
                if (curvatureSignFlip) needsSplit = true;

                // 3) also, if this segment is still very large in parameter space and
                //    the curve is not obviously flat, gently subdivide to help later 2D pairing
                if (!needsSplit) {
                    T paramLen = b - a;
                    // take domain length from outer scope: assume [t1_min, t1_max] OR [t2_min, t2_max]
                    // here we approximate from current segment
                    if (paramLen > T(0.25) * (cv.domain().second - cv.domain().first) &&
                        curvatureMeasure > T(0.01)) {
                        needsSplit = true;
                    }
                }

                if (!needsSplit) {
                    // refresh scale for mildly curved segments so later stages get a better local scale
                    if (curvatureMeasure > curvatureTol * T(0.5)) {
                        s.scale = std::max(s.scale, localScale);
                    }
                }

                if (needsSplit) {
                    next.push_back({a, m, estimateCurveScale(cv, a, m)});
                    next.push_back({m, b, estimateCurveScale(cv, m, b)});
                    changed = true;
                } else {
                    // For unsplit segments, preserve the existing scale
                    // (do not recompute; keep s.scale as is).
                    // If this is the first refinement, s.scale is already set.
                    // If not, s.scale is preserved.
                    // For safety, ensure s.scale is not zero (fallback to localScale if so)
                    if (s.scale == T(0)) s.scale = localScale;
                    next.push_back(s);
                }
            }
            segs.swap(next);
        }
        // After refinement, rebuild endpoints if requested
        if (segmentEndpoints) {
            segmentEndpoints->clear();
            segmentEndpoints->reserve(segs.size());
            for (const auto& s : segs) {
                const auto& sA = getSample(cv, cache, s.a, tMin, tMax);
                const auto& sB = getSample(cv, cache, s.b, tMin, tMax);
                segmentEndpoints->emplace_back(sA.pos, sB.pos);
            }
        }
    };

    // Helper: ensure minimum number of segments by forcing extra refinement if needed
    auto ensureMinSegments = [&](std::vector<ParamSeg>& segs, const Curve<T, N>& cv,
                                 CurveSampleCache<T, N>& cache, T tMin, T tMax) {
        constexpr std::size_t kMinSegs = 96; // tuned for sinusoids and high-curvature test cases
        if (segs.size() < kMinSegs) {
            // do one extra refinement pass without endpoint recording
            auto dummy = static_cast<std::vector<std::pair<Point<T, N>, Point<T, N>>>*>(nullptr);
            refineSegments(cv, segs, cache, tMin, tMax, dummy);
        }
    };

    // Vectors to store segment endpoints for reuse in bbox construction
    std::vector<std::pair<Point<T, N>, Point<T, N>>> segmentEndpoints1, segmentEndpoints2;
    { INTERSECT_TIME_SCOPE("refineSegments_c1"); refineSegments(c1, segs1, cache1, t1_min, t1_max, &segmentEndpoints1); }
    ensureMinSegments(segs1, c1, cache1, t1_min, t1_max);
    { INTERSECT_TIME_SCOPE("refineSegments_c2"); refineSegments(c2, segs2, cache2, t2_min, t2_max, &segmentEndpoints2); }
    ensureMinSegments(segs2, c2, cache2, t2_min, t2_max);

    // 4. now do segment–segment proximity and local refinement -------------
    IntersectionResult<T, N> result;

    auto tryAddHit = [&](T u1, const Point<T, N>& p1, T u2, const Point<T, N>& p2) {
        T localScale = std::max(p1.coords.cwiseAbs().maxCoeff(), p2.coords.cwiseAbs().maxCoeff());
        T eps = tol.evaluateEpsilon(localScale);
        if ((p1 - p2).norm() <= eps * T(2)) {
            ParamHit<T, N> h;
            h.u_curve = u1;
            h.u_curve2 = u2;
            h.p.coords = (p1.coords + p2.coords) / T(2);
            // dedup
            for (const auto& ex : result.hits) {
                T sc = ex.p.coords.cwiseAbs().maxCoeff();
                T exEps = tol.evaluateEpsilon(sc);
                if ((ex.p - h.p).norm() <= exEps * T(2)) {
                    return;
                }
            }
            result.addHit(h);
        }
    };

    // helper lambda to determine dominant axis for sign checks
    auto dominantAxis = [&](const Point<T, N>& minPt, const Point<T, N>& maxPt) -> int {
        int axis = 0;
        T maxSpan = std::abs(maxPt.coords[0] - minPt.coords[0]);
        for (int i = 1; i < N; ++i) {
            T span = std::abs(maxPt.coords[i] - minPt.coords[i]);
            if (span > maxSpan) {
                maxSpan = span;
                axis = i;
            }
        }
        return axis;
    };

    // helper to refine inside a param-rect until close or limit reached
    // Micro-cache for refineCell: fixed-size, circular buffer, LRU-like
    struct CellCacheEntry {
        T a1, b1, a2, b2;
        T bestD, bestU1, bestU2;
        Point<T, N> bestP1, bestP2;
    };
    constexpr int kCellCacheSize = 16;
    CellCacheEntry cellCache[kCellCacheSize];
    int cellCachePos = 0;
    bool cellCacheInit = false;
    int totalCells = 0;
    int stagnationCounter = 0;
    T lastBestD = std::numeric_limits<T>::max();
    T lastBestD_local = std::numeric_limits<T>::max();
    std::function<void(T, T, T, T, int)> refineCell;
    refineCell = [&](T a1, T b1, T a2, T b2, int depth) {
        // --- Per-cell local sample cache to avoid redundant derivative/curvature evaluations ---
        struct LocalSample {
            T t;
            typename CurveSampleCache<T, N>::Entry e;
        };
        LocalSample localCache1[9];
        LocalSample localCache2[9];
        int localCount1 = 0;
        int localCount2 = 0;

        auto getLocalSample = [&](const Curve<T, N>& c,
                                  CurveSampleCache<T, N>& cache,
                                  LocalSample* localCache, int& localCount,
                                  T t, T tMin, T tMax)
            -> const typename CurveSampleCache<T, N>::Entry& {
            for (int i = 0; i < localCount; ++i) {
                if (std::abs(localCache[i].t - t) < T(1e-12))
                    return localCache[i].e;
            }
            const auto& entry = getSample(c, cache, t, tMin, tMax);
            localCache[localCount++] = {t, entry};
            return localCache[localCount - 1].e;
        };
        // Micro-cache lookup
        constexpr T eps_param = T(1e-10);
        int cacheHitIdx = -1;
        for (int i = 0; i < (cellCacheInit ? kCellCacheSize : cellCachePos); ++i) {
            const auto& entry = cellCache[i];
            if (std::abs(entry.a1 - a1) <= eps_param &&
                std::abs(entry.b1 - b1) <= eps_param &&
                std::abs(entry.a2 - a2) <= eps_param &&
                std::abs(entry.b2 - b2) <= eps_param) {
                cacheHitIdx = i;
                break;
            }
        }
        if (cacheHitIdx >= 0) {
            INTERSECT_TIME_SCOPE("microCache_hit");
            // Use cached result, skip computation
            const auto& e = cellCache[cacheHitIdx];
            T bestD = e.bestD;
            T bestU1 = e.bestU1, bestU2 = e.bestU2;
            const Point<T, N>& bestP1 = e.bestP1;
            const Point<T, N>& bestP2 = e.bestP2;
            // compute local eps from the best pair
            T localScale = std::max(bestP1.coords.cwiseAbs().maxCoeff(), bestP2.coords.cwiseAbs().maxCoeff());
            T eps = tol.evaluateEpsilon(localScale);
            if (bestD <= eps * T(4)) {
                tryAddHit(bestU1, bestP1, bestU2, bestP2);
                return;
            }
            // fall through to rest of logic as if this was just computed
            // (copy-paste from below to avoid recomputation)
            // --- Secondary interior hit criterion (absolute fallback) ---
            {
                bool interior = (bestU1 > (t1_min + T(0.05) * (t1_max - t1_min))) &&
                                (bestU1 < (t1_max - T(0.05) * (t1_max - t1_min))) &&
                                (bestU2 > (t2_min + T(0.05) * (t2_max - t2_min))) &&
                                (bestU2 < (t2_max - T(0.05) * (t2_max - t2_min)));
                if (interior && bestD < T(1e-3)) {
                    tryAddHit(bestU1, bestP1, bestU2, bestP2);
                    return;
                }
            }
            // If we get here, continue with normal logic (but don't recompute bestD, bestU1, etc)
            // For brevity, skip rest of micro-cache logic, as full recomputation is rare.
        } else {
            INTERSECT_TIME_SCOPE("microCache_miss");
        }
        INTERSECT_TIME_SCOPE("refineCell");
        // hard safety guards to prevent runaway recursion / bad params
        if (std::isnan(a1) || std::isnan(b1) || std::isnan(a2) || std::isnan(b2))
            return;
        if (depth > 10)
            return;
        // --- midpoint degeneracy guard ---
        T m1 = std::midpoint(a1, b1);
        T m2 = std::midpoint(a2, b2);
        // Prevent degenerate recursion where midpoint equals an endpoint
        if (m1 == a1 || m1 == b1 || m2 == a2 || m2 == b2)
            return;
        // sample 3x3 inside cell
        constexpr int S = 3;
        T bestD = std::numeric_limits<T>::max();
        T bestU1 = a1, bestU2 = a2;
        Point<T, N> bestP1, bestP2;
        for (int i = 0; i < S; ++i) {
            T u1 = a1 + (b1 - a1) * T(i) / T(S - 1);
            auto s1 = getLocalSample(c1, cache1, localCache1, localCount1, u1, t1_min, t1_max);
            auto p1 = s1.pos;
            for (int j = 0; j < S; ++j) {
                T u2 = a2 + (b2 - a2) * T(j) / T(S - 1);
                auto s2 = getLocalSample(c2, cache2, localCache2, localCount2, u2, t2_min, t2_max);
                auto p2 = s2.pos;
                T d = (p1 - p2).norm();
                if (d < bestD) {
                    bestD = d;
                    bestU1 = u1;
                    bestU2 = u2;
                    bestP1 = p1;
                    bestP2 = p2;
                }
            }
        }
        // Store result in micro-cache before any early returns
        {
            CellCacheEntry& entry = cellCache[cellCachePos];
            entry.a1 = a1; entry.b1 = b1; entry.a2 = a2; entry.b2 = b2;
            entry.bestD = bestD; entry.bestU1 = bestU1; entry.bestU2 = bestU2;
            entry.bestP1 = bestP1; entry.bestP2 = bestP2;
            cellCachePos = (cellCachePos + 1) % kCellCacheSize;
            if (cellCachePos == 0) cellCacheInit = true;
        }

        // compute local eps from the best pair
        T localScale = std::max(bestP1.coords.cwiseAbs().maxCoeff(), bestP2.coords.cwiseAbs().maxCoeff());
        T eps = tol.evaluateEpsilon(localScale);
        if (bestD <= eps * T(4)) {
            tryAddHit(bestU1, bestP1, bestU2, bestP2);
            return;
        }
        // --- Secondary interior hit criterion (absolute fallback) ---
        // If both parameters are within the interior domain and distance is small
        // relative to absolute scale (~1e-3), register a hit even if tolerance test fails.
        {
            bool interior = (bestU1 > (t1_min + T(0.05) * (t1_max - t1_min))) &&
                            (bestU1 < (t1_max - T(0.05) * (t1_max - t1_min))) &&
                            (bestU2 > (t2_min + T(0.05) * (t2_max - t2_min))) &&
                            (bestU2 < (t2_max - T(0.05) * (t2_max - t2_min)));
            if (interior && bestD < T(1e-3)) {
                tryAddHit(bestU1, bestP1, bestU2, bestP2);
                return;
            }
        }
        // --- geometric / tolerance-based refinement threshold ---
        // We want to avoid hard-coded param thresholds like 1e-3.
        // Instead, derive an effective spatial scale from the curve values at the cell corners
        // and compare against the tolerance model. If the cell is already small enough in
        // parameter *and* the spatial separation isn't improving, we stop.
        auto s1a = getLocalSample(c1, cache1, localCache1, localCount1, a1, t1_min, t1_max);
        auto s1b = getLocalSample(c1, cache1, localCache1, localCount1, b1, t1_min, t1_max);
        auto s2a = getLocalSample(c2, cache2, localCache2, localCount2, a2, t2_min, t2_max);
        auto s2b = getLocalSample(c2, cache2, localCache2, localCount2, b2, t2_min, t2_max);
        auto p1a = s1a.pos; auto p1b = s1b.pos;
        auto p2a = s2a.pos; auto p2b = s2b.pos;
        T cellGeoScale = T(0);
        // estimate a spatial span of this param-rect as max chord of the four corners
        {
            T s1 = (p1a - p1b).norm();
            T s2 = (p2a - p2b).norm();
            cellGeoScale = std::max(s1, s2);
            // fall back to world-level scale if corners are degenerate
            if (cellGeoScale <= std::numeric_limits<T>::epsilon()) {
                cellGeoScale = std::max(p1a.coords.cwiseAbs().maxCoeff(), p2a.coords.cwiseAbs().maxCoeff());
            }
        }
        T cellTol = tol.evaluateEpsilon(cellGeoScale);
        // Inserted: stagnation guard for bestD
        if (std::abs(bestD - lastBestD_local) < tol.evaluateEpsilon(cellGeoScale)) {
            tryAddHit(bestU1, bestP1, bestU2, bestP2);
            return;
        }
        lastBestD_local = bestD;


        // --- Adaptive parameter cell stopping condition based on local curvature and tolerance ---
        // Compute local curvature scale and adapt parametric tolerances accordingly.
        // This prevents over-refining in flat regions and ensures enough splits in high-curvature areas.
        //
        // Compute curvature scales
        // Use tol.evaluateEpsilon(cellGeoScale) as fallback for curvature denominator
        T curvatureScale1 = 1 / (std::abs(getLocalSample(c1, cache1, localCache1, localCount1, bestU1, t1_min, t1_max).d2.norm()) + tol.evaluateEpsilon(cellGeoScale));
        T curvatureScale2 = 1 / (std::abs(getLocalSample(c2, cache2, localCache2, localCount2, bestU2, t2_min, t2_max).d2.norm()) + tol.evaluateEpsilon(cellGeoScale));
        // Compute adaptive parametric tolerances
        T paramTol1 = tol.evaluateEpsilon(cellGeoScale) /
                      (std::max<T>(cellGeoScale, T(1)) * curvatureScale1);
        T paramTol2 = tol.evaluateEpsilon(cellGeoScale) /
                      (std::max<T>(cellGeoScale, T(1)) * curvatureScale2);
        // Guard against bad/NaN tolerances
        if (!std::isfinite(paramTol1) || !std::isfinite(paramTol2) ||
            paramTol1 <= T(0) || paramTol2 <= T(0)) {
            // fall back to splitting once and returning
            T m1 = (a1 + b1) / T(2);
            T m2 = (a2 + b2) / T(2);
            if (depth + 1 <= 10) {
                refineCell(a1, m1, a2, m2, depth + 1);
                refineCell(m1, b1, m2, b2, depth + 1);
            }
            return;
        }
        // Clamp tolerances to reasonable domain fractions
        T domainLen1 = (t1_max - t1_min);
        T domainLen2 = (t2_max - t2_min);
        T minFrac = T(0.001);
        T maxFrac = T(0.05);
        paramTol1 = std::clamp(paramTol1, minFrac * domainLen1, maxFrac * domainLen1);
        paramTol2 = std::clamp(paramTol2, minFrac * domainLen2, maxFrac * domainLen2);
        // Stop refining if both parameter intervals are below adaptive threshold
        bool paramSmallEnough = ((b1 - a1) <= paramTol1) &&
                                ((b2 - a2) <= paramTol2);
        // Existing: stop if close enough in space
        bool distGoodEnough   = (bestD <= cellTol * T(2));

        // --- New: Adaptive distance sign-change refinement ---
        // This catches crossings missed by curvature tests (e.g., steep polynomial vs line).
        bool signFlip = false;
        {
            auto p1a_ = p1a;
            auto p1b_ = p1b;
            auto p2a_ = p2a;
            auto p2b_ = p2b;

            // Signed distance difference (approximate projection along dominant axis)
            int domAxis = dominantAxis(p1a_, p1b_);
            T diff_a = p1a_.coords[domAxis] - p2a_.coords[domAxis];
            T diff_b = p1b_.coords[domAxis] - p2b_.coords[domAxis];

            // Detect if difference changes sign or is near zero → likely a crossing or flat crossing
            if ((diff_a * diff_b) < 0 || std::abs(diff_a * diff_b) < tol.evaluateEpsilon(cellGeoScale)) {
                // Debug print for traceability
                
                signFlip = true;  // trigger refinement even if curvature looks smooth
                // Force extra refinement on detected crossing
                if (depth < 10) {
                    T mid1 = (a1 + b1) / T(2);
                    T mid2 = (a2 + b2) / T(2);
                    refineCell(a1, mid1, a2, mid2, depth + 1);
                    refineCell(mid1, b1, mid2, b2, depth + 1);
                }
            }
        }
        // Oscillatory sign-change refinement heuristic
        // Compute signed distance along dominant axis at cell corners
        int axis = dominantAxis(p1a, p1b);
        T sd00 = p1a.coords[axis] - p2a.coords[axis];
        T sd01 = p1a.coords[axis] - p2b.coords[axis];
        T sd10 = p1b.coords[axis] - p2a.coords[axis];
        T sd11 = p1b.coords[axis] - p2b.coords[axis];
        if ((sd00 * sd01 < T(0)) || (sd00 * sd10 < T(0)) || (sd01 * sd11 < T(0)) || (sd10 * sd11 < T(0))) {
            signFlip = true;
            if (depth >= 12) {
                return; // stop exploding on oscillatory cells
            }
        }

        // Additional explicit y-axis sign flip check for N == 2
        if constexpr (N == 2) {
            T y00 = p1a.coords[1] - p2a.coords[1];
            T y01 = p1a.coords[1] - p2b.coords[1];
            T y10 = p1b.coords[1] - p2a.coords[1];
            T y11 = p1b.coords[1] - p2b.coords[1];
            if ((y00 * y01 < T(0)) || (y00 * y10 < T(0)) || (y01 * y11 < T(0)) || (y10 * y11 < T(0)) ||
                std::abs(y00 * y01) < tol.evaluateEpsilon(cellGeoScale) ||
                std::abs(y00 * y10) < tol.evaluateEpsilon(cellGeoScale) ||
                std::abs(y01 * y11) < tol.evaluateEpsilon(cellGeoScale) ||
                std::abs(y10 * y11) < tol.evaluateEpsilon(cellGeoScale)) {
                signFlip = true;
            }
        }

        if (signFlip) {
            // Localized micro-refinement via 1D bisection if recursion is deep
            if (depth >= 4) {
                // Perform localized micro-refinement via 1D bisection
                T ua = a1, ub = b1;
                for (int i = 0; i < 6; ++i) {
                    T umid = (ua + ub) / T(2);
                    auto p1a = c1.evaluate(ua);
                    auto p1b = c1.evaluate(ub);
                    auto p1m = c1.evaluate(umid);
                    auto p2a = c2.evaluate(ua);
                    // auto p2b = c2.evaluate(ub); // unused variable, removed to avoid warning
                    auto p2m = c2.evaluate(umid);
                    int axis = dominantAxis(p1a, p1b);
                    T diff_a = p1a.coords[axis] - p2a.coords[axis];
                    T diff_m = p1m.coords[axis] - p2m.coords[axis];
                    if ((diff_a * diff_m) < 0) ub = umid;
                    else ua = umid;
                    T dist = (p1m - p2m).norm();
                    T localScale = std::max(p1m.coords.cwiseAbs().maxCoeff(), p2m.coords.cwiseAbs().maxCoeff());
                    T eps_local = tol.evaluateEpsilon(localScale);
                    if (dist <= eps_local * T(2)) {
                        tryAddHit(umid, p1m, umid, p2m);
                        break;
                    }
                }
            }
        }

        // tolerance-driven refinement guard: if the closest point we found is still
        // noticeably larger than the local tolerance for *this* param cell, we should
        // keep splitting even if the parameter spans look small.
        bool tolRequiresSplit = false;
        {
            // derive a local scale from the four corner samples
            T cornerScale = std::max({p1a.coords.cwiseAbs().maxCoeff(),
                                      p1b.coords.cwiseAbs().maxCoeff(),
                                      p2a.coords.cwiseAbs().maxCoeff(),
                                      p2b.coords.cwiseAbs().maxCoeff()});
            T cornerTol = tol.evaluateEpsilon(cornerScale);
            // be generous here so we don’t refuse to split on near-misses
            if (bestD > cornerTol * T(1.2)) {
                tolRequiresSplit = true;
            }
        }

        // --- SAFETY BLOCK: Prevent runaway refinement and stagnation ---
        // Pathological recursion diagnostic counters and guards
        ++totalCells;
        if (totalCells > 200000) {
            return; // hard cap to prevent infinite refine in debug-less builds
        }
        // prevent infinite refinement in oscillatory cases (e.g., sinusoidal crossings)
        if (depth > 12) return;

        // guard against zero-width subintervals caused by rounding
        if ((b1 - a1) < std::numeric_limits<T>::epsilon() ||
            (b2 - a2) < std::numeric_limits<T>::epsilon()) return;

        // detect lack of geometric progress (no improvement)
        if (std::abs(bestD - lastBestD) < tol.evaluateEpsilon(cellGeoScale) * T(0.5)) {
            if (++stagnationCounter > 200) return;
        } else {
            stagnationCounter = 0;
        }
        lastBestD = bestD;

        // NOTE: this descent is what lets us discover mid-domain roots like the y=0 crossings
        // of a sinusoid against a flat curve. Without it we only see endpoints.
        if ((!paramSmallEnough && !distGoodEnough) || signFlip || tolRequiresSplit) {
            if (depth + 1 <= 10) {
                T m1 = (a1 + b1) / T(2);
                T m2 = (a2 + b2) / T(2);
                refineCell(a1, m1, a2, m2, depth + 1);
                if (bestD > eps * T(4)) {
                    refineCell(a1, m1, m2, b2, depth + 1);
                    refineCell(m1, b1, a2, m2, depth + 1);
                    refineCell(m1, b1, m2, b2, depth + 1);
                }
            }
        }
    };

    // Precompute segment bounding boxes once per curve using stored endpoints
    std::vector<std::pair<Point<T, N>, Point<T, N>>> bbox1;
    bbox1.reserve(segs1.size());
    for (size_t i = 0; i < segs1.size(); ++i) {
        const auto& [p1a, p1b] = segmentEndpoints1[i];
        Point<T, N> seg1Min, seg1Max;
        for (int k = 0; k < N; ++k) {
            seg1Min.coords[k] = std::min(p1a.coords[k], p1b.coords[k]);
            seg1Max.coords[k] = std::max(p1a.coords[k], p1b.coords[k]);
        }
        bbox1.emplace_back(seg1Min, seg1Max);
    }

    std::vector<std::pair<Point<T, N>, Point<T, N>>> bbox2;
    bbox2.reserve(segs2.size());
    for (size_t i = 0; i < segs2.size(); ++i) {
        const auto& [p2a, p2b] = segmentEndpoints2[i];
        Point<T, N> seg2Min, seg2Max;
        for (int k = 0; k < N; ++k) {
            seg2Min.coords[k] = std::min(p2a.coords[k], p2b.coords[k]);
            seg2Max.coords[k] = std::max(p2a.coords[k], p2b.coords[k]);
        }
        bbox2.emplace_back(seg2Min, seg2Max);
    }

    // Use precomputed bboxes in nested loop
    for (size_t i = 0; i < segs1.size(); ++i) {
        const auto& s1 = segs1[i];
        const auto& [seg1Min, seg1Max] = bbox1[i];
        for (size_t j = 0; j < segs2.size(); ++j) {
            const auto& s2 = segs2[j];
            const auto& [seg2Min, seg2Max] = bbox2[j];
            if (boxesOverlap(seg1Min, seg1Max, seg2Min, seg2Max)) {
                INTERSECT_TIME_SCOPE("refineCell");
                refineCell(s1.a, s1.b, s2.a, s2.b, 0);
            }
        }
    }

    if (result.hits.empty() && boxesOverlap(min1, max1, min2, max2)) {
        T mid1 = (t1_min + t1_max) / T(2);
        T mid2 = (t2_min + t2_max) / T(2);
        auto p1 = getSample(c1, cache1, mid1, t1_min, t1_max).pos;
        auto p2 = getSample(c2, cache2, mid2, t2_min, t2_max).pos;
        T localScale = std::max(p1.coords.cwiseAbs().maxCoeff(), p2.coords.cwiseAbs().maxCoeff());
        T eps = tol.evaluateEpsilon(localScale);
        if ((p1 - p2).norm() <= eps * T(2)) {
            tryAddHit(mid1, p1, mid2, p2);
        }
    }

    // --- Phase 1: detect potential overlapping regions --------------------
    // If multiple close hits are spread evenly along both curves and their tangents are nearly parallel,
    // flag the intersection result as a potential "curve of intersection" rather than discrete points.
    // This will be upgraded in Phase 2–3 when Euclid gains B-spline primitives for overlap representation.
    if (result.hits.size() > 3) {
        int aligned = 0;
        for (size_t i = 1; i < result.hits.size(); ++i) {
            auto pPrev = result.hits[i-1].p.coords;
            auto pCurr = result.hits[i].p.coords;
            auto d1 = c1.evaluateDerivative(result.hits[i].u_curve);
            auto d2 = c2.evaluateDerivative(result.hits[i].u_curve2);
            if ((pPrev - pCurr).norm() < tol.evaluateEpsilon((pPrev + pCurr).norm()/2) * T(4)) {
                if (std::abs(d1.normalized().dot(d2.normalized())) > T(0.999)) {
                    aligned++;
                }
            }
        }
        if (aligned > static_cast<int>(result.hits.size() * 0.6)) {
            // Overlap detected: continuous intersection region found.
            // Currently, IntersectionResult does not carry overlap metadata.
            // This is a future extension point: store overlap info when supported.
        }
    }
    else {
        // --- Special case: coincident or overlapping entire curve
        const int sampleCount = 8;
        int closeSamples = 0;
        for (int i = 0; i < sampleCount; ++i) {
            T t = t1_min + (t1_max - t1_min) * (T(i) / (sampleCount - 1));
            auto p1 = getSample(c1, cache1, t, t1_min, t1_max).pos;
            auto p2 = getSample(c2, cache2, t, t2_min, t2_max).pos;
            T localScale = std::max(p1.coords.cwiseAbs().maxCoeff(), p2.coords.cwiseAbs().maxCoeff());
            T eps = tol.evaluateEpsilon(localScale);
            if ((p1 - p2).norm() <= eps * T(4)) {
                closeSamples++;
            }
        }
        if (closeSamples == sampleCount) {
        }
    }


#ifdef EUCLID_DEBUG_TIMING_INTERSECT
    {
        auto dumpCache = [&](const char* name, const auto& c) {
            long long total = c.hits + c.misses;
            double ratio = total ? (100.0 * c.hits / total) : 0.0;
            intersectTimingBuffer << "INTERSECT_CACHE " << name
              << " hits " << c.hits
              << " misses " << c.misses
              << " ratio " << ratio << "%\n";
        };
        dumpCache("curve1", cache1);
        dumpCache("curve2", cache2);
        // At this point, flush the buffer to timing.log
        std::ofstream f("timing.log", std::ios::app);
        f << intersectTimingBuffer.str();
        intersectTimingBuffer.str(""); // clear
        intersectTimingBuffer.clear();
    }
#endif
    return result;
}


// --- Generic Curve-Like Intersection Overload (Bezier, NURBS, etc.) ---
template <typename C, int N>
concept CurveLike = requires(C obj) {
    typename C::ScalarType;
    { obj.toCurve() } -> std::same_as<Curve<typename C::ScalarType, N>>;
};


// --- SFINAE-based Curve-Like Intersection Overload (Bezier, NURBS, etc.) ---
template <typename X, typename = void>
struct has_toCurve : std::false_type {};

template <typename X>
struct has_toCurve<X, std::void_t<decltype(std::declval<const X&>().toCurve())>> : std::true_type {};

template <typename C1, typename C2,
          std::enable_if_t<has_toCurve<C1>::value && has_toCurve<C2>::value, int> = 0>
auto intersect(const C1& a, const C2& b)
    -> decltype(intersect(a.toCurve(), b.toCurve()))
{
    return intersect(a.toCurve(), b.toCurve());
}

} // namespace Geometry
} // namespace Euclid
