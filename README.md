# 🜁 Project Euclid

**Euclid** is a modern, modular geometry kernel built from first principles for high-performance computational design. It is designed to be the mathematical foundation for a next-generation CAD ecosystem.

## 🏗️ TODO

## ⚙️ Core Linear Algebra Infrastructure

- [ ] `Metric`
    - [x] Core infrastrucure (inner products, etc.)
    - [ ] Measurement operations for all geometry primitives (Point, Line, Plane, Face, Surface, Volume)

## 📐 Geometry Primitives

- [x] `Point`
- [ ] `Line`
    - need to revisit logic for intersection tests
    - need to implement line segments
- [x] `Curve` - Parametric Linear & Nonlinear Curves
- [ ] `Plane`
- [ ] `Face`
    - [ ] Validate convexity, winding, etc. upon construction
- [ ] `Surface` 
- [ ] `Volume` ((ℝ³ → ℝⁿ))

## 🧮 Operations on Geometry Primitives
- [x] `Transformations` - Affine, Linear, Orthagonal, Special Orthagonal, Projective
- [ ] `Boolean Operations` — Boolean operations (union, subtract, intersect)
- [ ] `Collisions` — Collision and intersection logic


## 🔢 Numerical Methods

- [ ] `Sampling` - sampling/surface and volume division for meshing
- [ ] `Integration`
- [ ] `Differentiation`

## ✒️ Elementa Programing Language API via LLVM

- [ ] `Elementa` —  layer for input to kernel 
