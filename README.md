# 🜁 Project Euclid

**Euclid** is a modern, modular geometry kernel built from first principles for high-performance computational design. It is designed to be the mathematical foundation for a next-generation CAD ecosystem—lean, precise, and free from the legacy bloat of traditional systems like SolidWorks or OpenCascade.  

![Trefoil Knot Mesh](trefoil.png "Trefoil Knot Mesh")  

Written in **Haskell** for mathematical rigor and extensibility with meshing and rendering optimized in **Rust** for performance and memory safety, Euclid cleanly separates the geometry kernel from rendering and UI concerns. It supports arbitrary-dimensional spaces, Euclidean & non-Euclidean metrics, and type-safe operations on vectors, matrices, points/vertices, curves, surfaces, and volumes. Transformation logic is typeclass-driven and adaptable to metric context.

**Architecture**:
- ✒️ **Elementa**: A Ruby/OCaml-inspired DSL for expressive, readable geometry construction.
- 🧠 **Algebra and Geometry Core written in [Haskell](https://www.haskell.org)**: Functional backbone with strong typing, purity, and unit-tested primitives.
- ⚙️ **Rendering and Meshing in [Rust](https://www.rust-lang.org)**: High-performance meshing and rendering.

Euclid is for builders, researchers, and rebels—those who want to *design reality* using clean code and solid math, not overpriced GUI clickfests.

---

📦 **Current Phase**: V1 Refactor in Progress  
🛠️ Rebuilding core infrastructure for better composability, transformation logic, and Rust interop.  

- 📘 See [euclid-legacy](https://github.com/mksupreme92/euclid-legacy) repo for V0 prototype with most functionality completed: 
  - Geometry including points, lines, parametric curves, planes/faces, parametric surfaces, & parametric volumes in arbitrary dimensions  
  - Surface and Volume meshing  
---

## 🏗️ TODO

## ⚙️ Core Linear Algebra Infrastructure

- [x] `Algebra/Vector.hs` — Vector math (add, sub, norm, dot, etc.)
- [x] `Algebra/Matrix.hs` — Matrix operations
- [ ] `Algebra/Metric.hs`
    - [x] Core infrastrucure (inner products, etc.)
    - [ ] Measurement operations for vectors and all geometry primitives (Vertex, Edge, Plane, Face, Surface, Volume)

## 📐 Geometry Primitives

- [ ] `Geometry/Vertex.hs`
    - [x] Defined typeclass and constructor function
    - [ ] Defined instances in Transform.hs:
        - [x] Translateable
        - [x] RotatableAboutVertex
        - [ ] ReflectableAboutPlane
        - [ ] ProjectAbleOntoPlane
    - [ ] Defined instances in Metric.hs for measurement operations
    
- [x] `Geometry/Edge.hs`
    - [x] Defined typeclass and constructor function
        - infinite lines, rays, and segments
    - [ ] Defined instances in Transform.hs:
        - [x] Translateable
        - [ ] Scaleable (for Segment)
        - [x] RotatableAboutVertex
        - [ ] ReflectableAboutPlane
        - [ ] ProjectAbleOntoPlane
    - [ ] Defined instances in Metric.hs for measurement operations
    
- [ ] `Geometry/Curve.hs`
    - [ ] Defined typeclass and constructor functions
        - [ ] Linear Curves
            - test for 2D–5D, negatives, irrationals, and extrapolation
        - [ ] Parametric Nonlinear Curves
    - [ ] Defined instances in Transform.hs:
        - [ ] Translateable
        - [ ] Scaleable 
        - [ ] RotatableAboutVertex
        - [ ] ReflectableAboutPlane
        - [ ] ProjectAbleOntoPlane
    - [ ] Defined instances in Metric.hs for measurement operations
    
- [ ] `Geometry/Plane.hs`
    - [ ] Defined typeclass and constructor functions
    - [ ] Defined instances in Transform.hs:
        - [ ] Translateable
        - [ ] Scaleable 
        - [ ] RotatableAboutVertex
        - [ ] ReflectableAboutPlane
        - [ ] ProjectAbleOntoPlane
    - [ ] Defined instances in Metric.hs for measurement operations

- [ ] `Geometry/Face.hs` — Polygon face with 
    - [ ] Defined typeclass and constructor function
    - [ ] Validate convexity, winding, etc. upon construction
    - [ ] Defined instances in Transform.hs:
        - [ ] Translateable
        - [ ] Scaleable
        - [ ] Skewable
        - [ ] RotatableAboutVertex
        - [ ] ReflectableAboutPlane
        - [ ] ProjectAbleOntoPlane
    - [ ] Defined instances in Metric.hs for measurement operations

- [ ] `Geometry/Surface.hs` 
    - [ ] Defined typeclass and constructor function for parametric surfaces
    - [ ] Defined instances in Transform.hs:
        - [ ] Translateable
        - [ ] Scaleable
        - [ ] RotatableAboutVertex
        - [ ] ReflectableAboutPlane
        - [ ] ProjectAbleOntoPlane
    - [ ] Defined instances in Metric.hs for measurement operations
    - [ ] Meshing & area integration in Rust via FFI
    
- [ ] `Geometry/Volume.hs`
    - [ ] Parametric volume logic (ℝ³ → ℝⁿ)
    - [ ] Defined instances in Transform.hs:
        - [ ] Translateable
        - [ ] Scaleable
        - [ ] RotatableAboutVertex
        - [ ] ReflectableAboutPlane
        - [ ] ProjectAbleOntoPlane
    - [ ] Defined instances in Metric.hs for measurement operations
    - [ ] Tetrahedral meshing & volume integration in Rust via FFI

## 🧮 Algebraic Operations

- [ ] `Algebra/Transform.hs` — Rotations, translations, projections, scaling, etc.
    - [ ] ReflectableAboutPlane class for Geometry instances with defined locations in metric space (Vertex, Edge, Plane, Face, Surface, Volume)
    - [ ] ProjectAbleOntoPlane class for Geometry instances defined locations in metric space (Vertex, Edge, Plane, Face, Surface, Volume)
- [ ] `Algebra/Boolean.hs` — Boolean operations (union, subtract, intersect)
- [ ] `Algebra/Collision.hs` — Collision and intersection logic


## 🔢 Numerical Methods

- [ ] `Numerics/Sampling.hs` - sampling/surface and volume division for meshing
- [ ] `Numerics/Integrate.hs` — Runge-Kutta 4 integration
- [ ] `Numerics/Differentiate.hs` — Numerical differentiation (finite difference, etc.)
- [ ] `Numerics/Jacobian.hs` — Jacobian matrix and determinant computations for volume integration

## ✒️ Elementa Programing Language

- [ ] `Elementa.hs` — Elementa language layer for input to kernel
    - [ ] Grammar
    - [ ] Syntax
- [ ] `Parser.hs` — Parser/interpreter for Elementa DSL

## 🎬 Rendering and Visualization

- [ ] Rendering engine in Rust
- [ ] Real time animation visualization
