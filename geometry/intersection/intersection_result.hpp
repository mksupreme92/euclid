#pragma once

// ======================================================
// Intersection Base Struct
// ======================================================
struct IntersectionResult {
    bool intersects = false;
    std::string description;         // optional textual description

    IntersectionResult() = default;
    IntersectionResult(bool i, std::string desc = "")
        : intersects(i), description(std::move(desc)) {}
};
