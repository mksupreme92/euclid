# 🜁 Project Euclid

**Euclid** is a modern, modular geometry kernel built from first principles for high-performance computational design. It is designed to be the mathematical foundation for a next-generation CAD ecosystem—lean, precise, and free from the legacy bloat of traditional systems like SolidWorks or OpenCascade.

Written in **Haskell** for mathematical rigor and extensibility with meshing and rendering optimized in **Rust** for performace and memory safety, Euclid cleanly separates the geometry kernel from rendering and UI concerns. It supports arbitrary-dimensional spaces, Euclidean & non-Euclidean metrics, and type-safe operations on vectors, matrices, points/vertices, curves, surfaces, and volumes. Transformation logic is typeclass-driven and adaptable to metric context.

**Architecture**:
- ✒️ **Elementa**: A Ruby/OCaml-inspired DSL for expressive, readable geometry construction.
- 🧠 **Haskell Algebra and Geometry Core**: Functional backbone with strong typing, purity, and unit-tested primitives.
- ⚙️ **Rust Rendering and Meshing**: High-performance meshing and rendering.

Euclid is for builders, researchers, and rebels—those who want to *design reality* using clean code and solid math, not overpriced GUI clickfests.

---

## ✅ Status: Version 2 Refactor in Progress

---

## 🏗️ TODO

## ⚙️ Core Linear Algebra Infrastructure

- [x] `vector.hs` — Vector math (add, sub, norm, dot, etc.)
- [x] `matrix.hs` — Matrix operations
- [ ] `metric.hs`
    - [x] Core infrastrucure (inner products, etc.)
    - [ ] Measurement operations for vectors and all geometry primitives

## 📐 Geometry Primitives

- [x] `vertex.hs`
- [x] `edge.hs`
- [ ] `curve.hs` — Parametric curves, splines
    - [ ] Linear Curves
        - test for 2D–5D, negatives, irrationals, and extrapolation
    - [ ] Parametric Nonlinear Curves
- [ ] `plane.hs`
- [ ] `face.hs` — Polygon face with convexity, winding, etc.
- [ ] `surface.hs` — 
    - [ ] parametric surface
    - [ ] meshing & area integration in Rust via FFI
- [ ] `volume.hs` 
    — [ ] Parametric volume logic (ℝ³ → ℝⁿ)
    - [ ] tetrahedral meshing & volume integration in Rust via FFI

## 🧮 Algebraic Operations

- [ ] `transform.hs` — Rotations, translations, projections, scaling, etc. (note: metric based operations are WIP)
- [ ] `boolean.hs`— Boolean operations (union, subtract, intersect)
- [ ] `collision.hs` (optional) — Collision and intersection logic


## 🔢 Numerical Methods

- [ ] 'sampling.hs' - sampling/surface and volume division for meshing
- [x] `integrate.hs` — Runge-Kutta 4 (RK4) and general-purpose ODE integration
- [ ] `numerics/differentiate.hs` — Numerical differentiation (finite difference, etc.)
- [ ] `numerics/jacobian.hs` — Jacobian matrix and determinant computations for volume integration

## ✒️ Elementa Progaming Language

- [ ] `elementa.hs` — Elementa DSL language layer
    - [ ] Grammar
    - [ ] Syntax
- [ ] `parser.hs` — Parser/interpreter for Elementa DSL

## 🎬 Rendering and Visualization

- [ ] Rendering engine in Rust
- [ ] Real time animation visualization
