# multicool

A library for numerical root-finding on systems of multivariate polynomial equations.
If you have a system of equations that looks like the following, then this library is for you:

$$
\begin{align}
x^2 + (y + 0.25)^2 &= 0 \\
(x - 1)^2 + y^2 &= 0 \\
\end{align}
$$

where $x \in [x_{min}, x_{max}], y \in [y_{min}, y_{max}]$

Requires:

- At least as many equations as unknowns
- Initial bounding region on the unknowns

# Example

```rust
use multicool::{Monomial, MultivarPoly, MultivarPolySystem}

//       x^2 + (y + 0.25)^2 = 1
// x^2 + y^2 + 0.5y - 0.875 = 0
let circle1 = MultivarPoly::new().add_monomials([
    Monomial::new(1.0, [2, 0]),
    Monomial::new(1.0, [0, 2]),
    Monomial::new(0.5, [0, 1]),
    Monomial::new(-0.875, [0, 0]),
]);

// (x - 1)^2 + y^2 = 1
//  x^2 - 2x + y^2 = 0
let circle2 = MultivarPoly::new().add_monomials([
    Monomial::new(1.0, [2, 0]),
    Monomial::new(-2.0, [1, 0]),
    Monomial::new(1.0, [0, 2]),
]);

let system =
    MultivarPolySystem::from_polys(
        [circle1, circle2],
        [(-1.0, 2.0), (-1.0, 1.0)]
    ).unwrap();

const EPS: f64 = 1e-9;
let roots: Vec<[f64; 2]> = system.roots(EPS).unwrap();
```

# Motivation

This library is heavily motivated by computational geometry problems like finding the intersections of curves and surfaces. In particular CAD and CAM applications.

# References

The algorithms in this library were developed from the following publications.

- [Shape Interrogation for Computer Aided Design and Manufacturing](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/)
- [Solving multivariate polynomial systems using hyperplane arithmetic and linear programming](https://www.sciencedirect.com/science/article/abs/pii/S0010448513001620)

# License

Dual-licensed to be compatible with the Rust project.

Licensed under the Apache License, Version 2.0
http://www.apache.org/licenses/LICENSE-2.0 or the MIT license
http://opensource.org/licenses/MIT, at your
option. This file may not be copied, modified, or distributed
except according to those terms.
