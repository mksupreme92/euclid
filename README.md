# 🜁 Project Euclid

**Euclid** is a modern, modular geometry kernel built from first principles for high-performance computational design. It is designed to be the mathematical foundation for a next-generation CAD and multiphysics simulation ecosystem.

## 🏗️ TODO

## ⚙️ Core Infrastructure

- [x] `Tolearance Model`

## 📐 Geometry Primitives

- [x] `Point`
- [x] `Line`
- [ ] `Curve`
    - [x] Parametric Curve
    - [ ] NURBS Curve
- [x] `Plane`
- [x] `Face`
- [ ] `Surface`
    - [x] Parametric Surface
    - [ ] NURBS Surface
- [x] `Volume` ((ℝ³ → ℝⁿ))


## 🧮 Operations on Geometry Primitives
- [x] `Transformations` - Affine, Linear, Orthagonal, Special Orthagonal, Projective

### Intersections (Required for NURBS Support)

#### Point Intersections
- [x] `Point - Point`
- [x] `Point - Line` → returns closest point
- [ ] `Point - Segment` → returns closest point
- [x] `Point - Plane` → returns projected point
- [ ] `Point - Face` → returns projected point
- [ ] `Point - Curve` → returns closest point on curve
- [ ] `Point - Surface` → returns projected point on surface

#### Line Intersections
- [x] `Line - Line` → returns point (or line if coincident)
- [ ] `Line - Segment` → returns point (or segment if colinear)
- [x] `Line - Plane` → returns point or line
- [ ] `Line - Face` → returns point or segment
- [ ] `Line - Curve` → returns points of intersection
- [ ] `Line - Surface` → returns point(s) of intersection

#### Segment Intersections
- [ ] `Segment - Segment` → returns point (or segment if colinear)
- [ ] `Segment - Plane` → returns point or segment
- [ ] `Segment - Face` → returns point or segment
- [ ] `Segment - Curve` → returns points of intersection
- [ ] `Segment - Surface` → returns point(s) of intersection

#### Curve Intersections
- [ ] `Curve - Curve` → returns points or overlapping segments
- [ ] `Curve - Surface` → returns curve(s) of intersection

#### Plane Intersections
- [ ] `Plane - Plane` → returns line (or plane if coincident)
- [ ] `Plane - Face` → returns line (or face if coplanar)
- [ ] `Plane - Surface` → returns curve of intersection

#### Face Intersections
- [ ] `Face - Face` → returns line (or face if coplanar)
- [ ] `Face - Surface` → returns curve(s) of intersection

#### Surface Intersection
- [ ] `Surface - Surface` → returns curve(s) of intersection
