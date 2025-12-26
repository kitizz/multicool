use polycool::Poly;
use smallvec::SmallVec;

use crate::{BezierSurface, MulticoolError, bounding_box::BoundingBox};

pub struct MultivarPolySystem<const NUM_VARS: usize> {
    surfaces: Vec<BezierSurface<NUM_VARS>>,
}

impl<const NUM_VARS: usize> MultivarPolySystem<NUM_VARS> {
    pub fn new(surfaces: Vec<BezierSurface<NUM_VARS>>) -> Self {
        const {
            assert!(NUM_VARS >= 1, "Number of variables must be at least 1");
        }
        Self { surfaces }
    }

    pub fn from_polys(
        polys: impl IntoIterator<Item = crate::MultivarPoly<NUM_VARS>>,
        domain: [(f64, f64); NUM_VARS],
    ) -> Result<Self, MulticoolError> {
        let surfaces = polys
            .into_iter()
            .map(|poly| poly.bezier_surface(domain))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self { surfaces })
    }

    /// Find the roots of the system using the Bezier Projected Polyhedra method.
    ///
    /// This method projects the Bezier surfaces onto each (x_i, x_{NUM_VARS}) plane
    /// and finds the intersection regions where x_{NUM_VARS} = 0.
    pub fn roots_pp(&self, tol: f64) -> Result<Vec<[f64; NUM_VARS]>, MulticoolError> {
        let refiner = refinement::Refiner::ProjectedPolyhedra(
            refinement::ProjectedPolyhedraRefiner::new(self.surfaces.len()),
        );
        self.roots(tol, refiner)
    }

    /// Find the roots of the system using the Bezier Linear Programming method.
    ///
    /// This method projects hyperplanes that bound each surface onto the x_{NUM_VARS} = 0
    /// hyperplane. Linear programming is used to find the axis-aligned bounding region
    /// of those projections.
    /// Two hyperplanes are calculated for each surface. They are parallel to the tangent
    /// of the surface at the center of the valid region, and are offset such that they
    /// bound all control points of the surface.
    pub fn roots_lp(&self, tol: f64) -> Result<Vec<[f64; NUM_VARS]>, MulticoolError> {
        log::info!("Num surfaces: {}", self.surfaces.len());
        let refiner = refinement::Refiner::HyperplanesLp(refinement::HyperplanesLpRefiner::new(
            self.surfaces.len(),
        ));
        self.roots(tol, refiner)
    }

    pub fn roots(
        &self,
        tol: f64,
        mut refiner: refinement::Refiner<NUM_VARS>,
    ) -> Result<Vec<[f64; NUM_VARS]>, MulticoolError> {
        // Check that we have at least as many equations as variables
        if self.surfaces.len() < NUM_VARS {
            return Err(MulticoolError::UnderdefinedSystem {
                num_vars: NUM_VARS,
                num_equations: self.surfaces.len(),
            });
        }

        // First, compute the valid regions for each variable by intersecting the convex hulls
        // of each surface's control points with the hyperplane at x_{NUM_VARS} = 0
        let mut region_queue = vec![self.surfaces[0].domain];
        let mut root_regions = smallvec::smallvec![];

        // Re-use surfaces vector to avoid reallocations
        let mut surfaces = self.surfaces.clone();

        let mut iteration = 0;
        while let Some(mut valid_regions) = region_queue.pop() {
            for s in 0..surfaces.len() {
                surfaces[s] = self.surfaces[s].clone().subsection(valid_regions);
            }

            // All surfaces should have the same domain
            log::info!("It {iteration}: valid region: {:#?}", valid_regions);
            iteration += 1;

            // Try to reduce the size of valid_regions.
            refiner.refine_regions(&surfaces, &mut valid_regions)?;

            log::info!("Refined region: {:?}", valid_regions);

            // If any valid_regions are invalid, return no roots
            if !valid_regions.is_valid() {
                continue;
            }

            // If the largest valid region is sufficiently small, stop
            let (largest_d, largest_size) = valid_regions.largest_side();
            if largest_size < tol {
                root_regions.push(valid_regions);
                continue;
            }

            // Subdivide each surface along the largest dimension at the midpoint
            let (lower_box, upper_box) = valid_regions.subdivide(largest_d);

            log::info!(
                "Subdividing into: lower {:?}, upper {:?}",
                lower_box,
                upper_box
            );
            region_queue.push(lower_box);
            region_queue.push(upper_box);
        }

        log::info!("Initial root regions: {:#?}", root_regions);

        merge_overlapping_regions(&mut root_regions);

        let mut roots = Vec::with_capacity(root_regions.len());
        for region in root_regions {
            let mut root = [0.0; NUM_VARS];
            for d in 0..NUM_VARS {
                root[d] = 0.5 * (region[d].0 + region[d].1);
            }
            roots.push(root);
        }
        Ok(roots)
    }
}

impl MultivarPolySystem<2> {
    pub fn from_beziers<const A_MAX: usize, const B_MAX: usize>(
        a: [[f64; 2]; A_MAX],
        b: [[f64; 2]; B_MAX],
        domain: [(f64, f64); 2],
    ) -> Self {
        let (a_poly_x, a_poly_y) = Self::bezier_poly(a);
        let (b_poly_x, b_poly_y) = Self::bezier_poly(b);

        let eq1 = Self::multivar_eq(a_poly_x, b_poly_x);
        let eq2 = Self::multivar_eq(a_poly_y, b_poly_y);

        Self::from_polys([eq1, eq2], domain).unwrap()
    }

    fn bezier_poly<const N: usize>(points: [[f64; 2]; N]) -> (Poly<N>, Poly<N>) {
        let xs = points.map(|p| p[0]);
        let ys = points.map(|p| p[1]);
        let poly_x = Self::bezier_poly_1d(xs);
        let poly_y = Self::bezier_poly_1d(ys);
        (poly_x, poly_y)
    }

    fn bezier_poly_1d<const N: usize>(v: [f64; N]) -> Poly<N> {
        // p_i = sum_{k=0 to i} (b_k * binomial(n, i) * binomial(i, k) * (-1)^(i-k))
        let mut coeffs = [0.0; N];
        let n = (N - 1) as u8;
        for i in 0..N {
            let binom_n_i = crate::binomial_coefficient(n, i as u8);
            for k in 0..=i {
                let binom_i_k = crate::binomial_coefficient(i as u8, k as u8);
                let sign = if (i - k) % 2 == 0 { 1.0 } else { -1.0 };
                coeffs[i] += v[k] * (binom_n_i as f64) * (binom_i_k as f64) * sign;
            }
        }
        Poly::new(coeffs)
    }

    fn multivar_eq<const N: usize, const M: usize>(
        lhs: Poly<N>,
        rhs: Poly<M>,
    ) -> crate::MultivarPoly<2> {
        let mut monomials = Vec::with_capacity(N + M);
        for (i, &coeff_lhs) in lhs.coeffs().iter().enumerate() {
            let mon = crate::Monomial::new(coeff_lhs, [i as u8, 0]);
            monomials.push(mon);
        }
        for (j, &coeff_rhs) in rhs.coeffs().iter().enumerate() {
            let mon = crate::Monomial::new(-coeff_rhs, [0, j as u8]);
            monomials.push(mon);
        }
        crate::MultivarPoly::<2>::new().add_monomials(monomials)
    }
}

/// Merge overlapping regions in-place by intersecting them.
///
/// There is no guarantee on the order of the resulting regions.
fn merge_overlapping_regions<const NUM_VARS: usize>(
    regions: &mut SmallVec<[BoundingBox<NUM_VARS>; 4]>,
) {
    if regions.len() <= 1 {
        return;
    }
    let mut i = 0;
    while i < regions.len() - 1 {
        let mut j = i + 1;
        while j < regions.len() {
            if regions[i].overlaps(&regions[j]) {
                // Merge regions i and j by intersecting them
                for d in 0..NUM_VARS {
                    regions[i][d].0 = regions[i][d].0.max(regions[j][d].0);
                    regions[i][d].1 = regions[i][d].1.min(regions[j][d].1);
                }
                regions.swap_remove(j);
            } else {
                j += 1;
            }
        }
        i += 1;
    }
}

mod refinement;

#[cfg(test)]
mod tests {
    // use pretty_assertions as pa;

    use crate::{Monomial, MultivarPoly};

    use super::*;

    #[test]
    fn system_roots_cross() {
        crate::test_utils::init_test_logger();

        // x - y = 0
        let line1 = MultivarPoly::new()
            .add_monomials([Monomial::new(1.0, [1, 0]), Monomial::new(-1.0, [0, 1])]);
        // x + y = 0
        let line2 = MultivarPoly::new()
            .add_monomials([Monomial::new(1.0, [1, 0]), Monomial::new(1.0, [0, 1])]);

        let system =
            MultivarPolySystem::from_polys([line1, line2], [(-1.0, 1.0), (-1.0, 1.0)]).unwrap();

        const EPS: f64 = 1e-6;
        let roots = system.roots_pp(EPS).unwrap();
        approx::assert_abs_diff_eq!(
            ApproxRoots(roots),
            ApproxRoots(vec![[0.0, 0.0]]),
            epsilon = EPS
        );
    }

    #[test]
    fn system_roots_circles() {
        crate::test_utils::init_test_logger();

        // x^2 + y^2 = 1
        let circle1 = MultivarPoly::new().add_monomials([
            Monomial::new(1.0, [2, 0]),
            Monomial::new(1.0, [0, 2]),
            Monomial::new(-1.0, [0, 0]),
        ]);

        // (x - 1)^2 + y^2 = 1
        let circle2 = MultivarPoly::new().add_monomials([
            Monomial::new(1.0, [2, 0]),
            Monomial::new(-2.0, [1, 0]),
            Monomial::new(1.0, [0, 2]),
        ]);

        let system =
            MultivarPolySystem::from_polys([circle1, circle2], [(-1.0, 2.0), (-1.0, 1.0)]).unwrap();

        const EPS: f64 = 1e-6;
        let roots = system.roots_pp(EPS).unwrap();
        let root_y = 3f64.sqrt() / 2.0;
        approx::assert_abs_diff_eq!(
            ApproxRoots(roots),
            ApproxRoots(vec![[0.5, root_y], [0.5, -root_y]]),
            epsilon = EPS
        );
    }

    #[test]
    fn system_roots_beziers() {
        crate::test_utils::init_test_logger();

        let system = MultivarPolySystem::from_beziers(
            [[0.0, 0.0], [1.0, 2.0], [2.0, -1.0], [3.0, 1.0]],
            [[0.0, 1.0], [1.0, -1.0], [2.0, 2.0], [3.0, 0.0]],
            [(0.0, 1.0), (0.0, 1.0)],
        );

        const EPS: f64 = 1e-9;
        let roots = system.roots_pp(EPS).unwrap();
        approx::assert_abs_diff_eq!(
            ApproxRoots(roots),
            ApproxRoots(vec![
                [0.11270166537925831, 0.11270166537925831],
                [0.5, 0.5],
                [0.8872983346207417, 0.8872983346207417],
            ]),
            epsilon = EPS
        );
    }

    #[test]
    fn system_roots_lp_cross() {
        crate::test_utils::init_test_logger();

        // x - y = 0
        let line1 = MultivarPoly::new()
            .add_monomials([Monomial::new(1.0, [1, 0]), Monomial::new(-1.0, [0, 1])]);
        // x + y = 0
        let line2 = MultivarPoly::new()
            .add_monomials([Monomial::new(1.0, [1, 0]), Monomial::new(1.0, [0, 1])]);

        let system =
            MultivarPolySystem::from_polys([line1, line2], [(-1.0, 1.0), (-1.0, 1.0)]).unwrap();

        const EPS: f64 = 1e-6;
        let roots = system.roots_lp(EPS).unwrap();
        approx::assert_abs_diff_eq!(
            ApproxRoots(roots),
            ApproxRoots(vec![[0.0, 0.0]]),
            epsilon = EPS
        );
    }

    #[test]
    fn system_roots_lp_circles() {
        crate::test_utils::init_test_logger();

        // x^2 + y^2 = 1
        let circle1 = MultivarPoly::new().add_monomials([
            Monomial::new(1.0, [2, 0]),
            Monomial::new(1.0, [0, 2]),
            Monomial::new(-1.0, [0, 0]),
        ]);

        // (x - 1)^2 + y^2 = 1
        let circle2 = MultivarPoly::new().add_monomials([
            Monomial::new(1.0, [2, 0]),
            Monomial::new(-2.0, [1, 0]),
            Monomial::new(1.0, [0, 2]),
        ]);

        let system =
            MultivarPolySystem::from_polys([circle1, circle2], [(-1.0, 2.0), (-1.0, 1.0)]).unwrap();

        const EPS: f64 = 1e-6;
        let roots = system.roots_lp(EPS).unwrap();
        let root_y = 3f64.sqrt() / 2.0;
        approx::assert_abs_diff_eq!(
            ApproxRoots(roots),
            ApproxRoots(vec![[0.5, root_y], [0.5, -root_y]]),
            epsilon = EPS
        );
    }

    #[test]
    fn system_roots_lp_beziers() {
        crate::test_utils::init_test_logger();

        let system = MultivarPolySystem::from_beziers(
            [[0.0, 0.0], [1.0, 2.0], [2.0, -1.0], [3.0, 1.0]],
            [[0.0, 1.0], [1.0, -1.0], [2.0, 2.0], [3.0, 0.0]],
            [(0.0, 1.0), (0.0, 1.0)],
        );

        const EPS: f64 = 1e-9;
        let roots = system.roots_lp(EPS).unwrap();
        approx::assert_abs_diff_eq!(
            ApproxRoots(roots),
            ApproxRoots(vec![
                [0.11270166537925831, 0.11270166537925831],
                [0.5, 0.5],
                [0.8872983346207417, 0.8872983346207417],
            ]),
            epsilon = EPS
        );
    }

    #[test]
    fn system_roots_lp_cylinders_3d() {
        crate::test_utils::init_test_logger();

        const DIMS: usize = 3;

        let mut polys = Vec::new();
        for d in 0..DIMS {
            polys.push(crate::test_utils::cylinder_poly(d));
        }

        let system =
            MultivarPolySystem::from_polys(polys, [(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)]).unwrap();

        let expected_roots = {
            let root_val = (1f64 / (DIMS - 1) as f64).sqrt();
            let pm = [-root_val, root_val];
            let mut roots = Vec::new();
            for (x, y, z) in itertools::iproduct!(pm, pm, pm) {
                roots.push([x, y, z]);
            }
            roots
        };

        const EPS: f64 = 1e-9;
        let roots = system.roots_lp(EPS).unwrap();
        approx::assert_abs_diff_eq!(
            ApproxRoots(roots),
            ApproxRoots(expected_roots),
            epsilon = EPS
        );
    }

    #[derive(Clone)]
    struct ApproxRoots<const V: usize>(Vec<[f64; V]>);

    impl<const V: usize> core::fmt::Debug for ApproxRoots<V> {
        fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
            self.0.fmt(f)
        }
    }

    impl<const V: usize> approx::AbsDiffEq for ApproxRoots<V> {
        type Epsilon = f64;

        fn default_epsilon() -> Self::Epsilon {
            1e-6
        }

        fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
            let a = &self.0;
            let b = &other.0;
            if a.len() != b.len() {
                return false;
            }

            let mut distances = Vec::with_capacity(a.len() * b.len());

            for ind_a in 0..a.len() {
                for ind_b in 0..b.len() {
                    let mut dist = 0.0;
                    for d in 0..V {
                        let delta = a[ind_a][d] - b[ind_b][d];
                        dist += delta * delta;
                    }
                    distances.push((dist.sqrt(), ind_a, ind_b));
                }
            }
            distances.sort_by(|(dist_a, _, _), (dist_b, _, _)| dist_a.partial_cmp(dist_b).unwrap());

            // Greedily match closest pairs
            let mut used_a = vec![false; a.len()];
            let mut used_b = vec![false; b.len()];
            let mut matches = 0;
            for (dist, ind_a, ind_b) in distances {
                if dist > epsilon || matches >= a.len() {
                    break;
                }
                if !used_a[ind_a] && !used_b[ind_b] {
                    used_a[ind_a] = true;
                    used_b[ind_b] = true;
                    matches += 1;
                }
            }
            matches == a.len()
        }
    }

    impl<const V: usize> PartialEq for ApproxRoots<V> {
        fn eq(&self, other: &Self) -> bool {
            self.0 == other.0
        }
    }
}
