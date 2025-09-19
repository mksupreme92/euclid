#pragma once
#include <cassert>

namespace Euclid {
    inline size_t SpaceDimension = 3;  // default

    inline void setSpaceDimension(size_t dim) {
        assert(dim > 0);  // still catch negative/zero
        SpaceDimension = dim;
    }

    inline size_t getSpaceDimension() { return SpaceDimension; }
}
