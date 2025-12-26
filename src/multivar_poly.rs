use smallvec::{SmallVec, smallvec};

use crate::{BezierSurface, Monomial, MulticoolError, binomial_product};

#[derive(Debug, Clone)]
pub struct MultivarPoly<const NUM_VAR: usize> {
    pub(crate) terms: SmallVec<[Monomial<NUM_VAR>; 8]>,
}

impl<const NUM_VAR: usize> MultivarPoly<NUM_VAR> {
    pub fn new() -> Self {
        Self {
            terms: SmallVec::new(),
        }
    }

    pub fn from_monomials(monomials: impl IntoIterator<Item = Monomial<NUM_VAR>>) -> Self {
        Self {
            terms: monomials.into_iter().collect(),
        }
    }

    fn terms(&self) -> impl Iterator<Item = &Monomial<NUM_VAR>> {
        self.terms.iter()
    }

    /// Evaluate the multivariate polynomial at the given variable values
    pub fn eval(&self, vars: &[f64; NUM_VAR]) -> f64 {
        let mut result = 0.0;
        for monomial in &self.terms {
            result += monomial.eval(vars);
        }
        result
    }

    pub fn sub_polys<const N: usize>(
        &self,
        polys: &[[f64; N]; NUM_VAR],
    ) -> Result<Self, MulticoolError> {
        let mut output = MultivarPoly::<NUM_VAR>::new();

        for monomial in &self.terms {
            let sub_poly = monomial.sub_polys::<N>(polys)?;
            output = output.merge(sub_poly);
        }

        Ok(output)
    }

    /// Calculate the Bernstein control points for a given multivariate polynomial
    ///
    /// Note that the original polynomial cannot be reconstructed from these control points alone.
    ///
    /// # Parameters
    /// - `PDIM`: The number of dimensions for the Bernstein basis. Requires `PDIM <= NUM_VAR + 1`.
    ///    If `PDIM < NUM_VAR + 1`, the extra variables are truncated as if they didn't exist.
    ///
    pub fn bezier_surface(
        &self,
        bbox: [(f64, f64); NUM_VAR],
    ) -> Result<BezierSurface<NUM_VAR>, MulticoolError> {
        // Perform affine parameter transformation from bbox to unit box [0, 1]^NUM_VAR
        // by substituting x_i = (t_i - min_i) / (max_i - min_i) = 1/Δi * t_i + (-min_i/Δi)
        let affines: [[f64; 2]; NUM_VAR] = std::array::from_fn(|i| {
            let (min_i, max_i) = bbox[i];
            [min_i, max_i - min_i]
        });
        let unit_poly = self.sub_polys(&affines)?;

        // Ref: https://adrianbowyer.com/Publications/berchtold2000.pdf
        // 1. Calculate max exp/degree in each variable across all monomials
        let max_degs = {
            let mut max_degs = [0u8; NUM_VAR];
            for monomial in unit_poly.terms() {
                // By only taking the first PDIM - 1 variables, we truncate any extra variables.
                for (i, &exp) in monomial.exp.iter().enumerate() {
                    max_degs[i] = max_degs[i].max(exp);
                }
            }
            max_degs
        };

        let grid_size = std::array::from_fn(|i| (max_degs[i] + 1).max(2));
        let mut surface = BezierSurface::zeros(grid_size.map(|s| s as usize), bbox);

        // 2. Sum over Eq (3) in the paper
        for control_inds in crate::gridex_excl(grid_size) {
            let bern_coeff = unit_poly
                .terms()
                .filter(|m| m.exp_all_le(&control_inds))
                .map(|m| {
                    let num = binomial_product(control_inds, m.exp) as f64;
                    let den = binomial_product(max_degs, m.exp) as f64;
                    m.coeff * num / den
                })
                .sum();

            surface.set_coeff(control_inds.map(usize::from), bern_coeff);
        }
        Ok(surface)
    }

    pub fn merge(mut self, other: Self) -> Self {
        let mut merged = SmallVec::with_capacity(self.terms.len() + other.terms.len());

        let mut a_iter = self.terms.into_iter().peekable();
        let mut b_iter = other.terms.into_iter().peekable();

        // The inputs should already be sorted. Keep them sorted while merging.
        while let (Some(a), Some(b)) = (a_iter.peek(), b_iter.peek()) {
            if a.exp < b.exp {
                merged.push(a_iter.next().unwrap());
            } else if a.exp > b.exp {
                merged.push(b_iter.next().unwrap());
            } else {
                // Equal exponents, combine coefficients
                let combined_coeff = a.coeff + b.coeff;
                if combined_coeff != 0.0 {
                    merged.push(Monomial::new(combined_coeff, a.exp));
                }
                a_iter.next();
                b_iter.next();
            }
        }
        for a in a_iter {
            merged.push(a);
        }
        for b in b_iter {
            merged.push(b);
        }
        self.terms = merged;
        self
    }

    pub fn add_monomials(self, monomials: impl IntoIterator<Item = Monomial<NUM_VAR>>) -> Self {
        let mut other = smallvec![];
        for m in monomials {
            other.push(m);
        }
        other.sort();
        self.merge(MultivarPoly { terms: other })
    }
}

#[cfg(test)]
impl<const NUM_VAR: usize> approx::AbsDiffEq for MultivarPoly<NUM_VAR> {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let a_terms = self.terms().collect::<Vec<_>>();
        let b_terms = other.terms().collect::<Vec<_>>();
        if a_terms.len() != b_terms.len() {
            return false;
        }

        // Try to find a full unordered match between a_terms and b_terms
        // Assumes no duplicate exponents within each polynomial
        let a_matched = a_terms
            .iter()
            .all(|a| b_terms.iter().any(|b| a.abs_diff_eq(b, epsilon)));
        let b_matched = b_terms
            .iter()
            .all(|b| a_terms.iter().any(|a| a.abs_diff_eq(b, epsilon)));
        a_matched && b_matched
    }
}

impl<const NUM_VAR: usize> PartialEq for MultivarPoly<NUM_VAR> {
    fn eq(&self, other: &Self) -> bool {
        let a_terms = self.terms().collect::<Vec<_>>();
        let b_terms = other.terms().collect::<Vec<_>>();
        if a_terms.len() != b_terms.len() {
            return false;
        }
        use std::iter::zip;
        zip(a_terms, b_terms).all(|(a, b)| a == b)
    }
}

impl<const NUM_VAR: usize> std::ops::Add for MultivarPoly<NUM_VAR> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.merge(rhs)
    }
}

#[cfg(test)]
mod tests {
    use crate::{Monomial, test_utils::unit_box};

    use super::*;
    use approx::assert_abs_diff_eq;
    use assertables::assert_ok;
    use pretty_assertions as pa;

    #[test]
    fn bezier_surface_1d_linear() {
        let poly = MultivarPoly::from_monomials([Monomial {
            coeff: 1.0,
            exp: [1],
        }]);

        let controls = assert_ok!(poly.bezier_surface(unit_box()));

        let expected = BezierSurface::new(&[0.0, 1.0], [2], unit_box());
        assert_abs_diff_eq!(controls, expected);
    }

    #[test]
    fn bezier_surface_1d_linear_offset() {
        let poly = MultivarPoly::from_monomials([
            //
            Monomial::new(2.0, [0]),
            Monomial::new(1.0, [1]),
        ]);

        let controls = assert_ok!(poly.bezier_surface(unit_box()));

        // Points interpolate a plane.
        let expected = BezierSurface::new(&[2.0, 3.0], [2], unit_box());
        assert_abs_diff_eq!(controls, expected);
    }

    #[test]
    fn bezier_surface_1d_linear_offset_affine() {
        let poly = MultivarPoly::from_monomials([
            //
            Monomial::new(2.0, [0]),
            Monomial::new(1.0, [1]),
        ]);

        let controls = poly.bezier_surface([(1.0, 3.0)]);

        // Points interpolate a plane.
        let expected = BezierSurface::new(&[3.0, 5.0], [2], [(1.0, 3.0)]);
        assert_abs_diff_eq!(controls.unwrap(), expected);
    }

    #[test]
    fn bezier_surface_2d_linear() {
        let poly = MultivarPoly::from_monomials([
            Monomial::new(1.0, [1, 0]),
            Monomial::new(2.0, [0, 1]),
            Monomial::new(3.0, [0, 0]),
        ]);

        let controls = assert_ok!(poly.bezier_surface(unit_box()));

        // Points interpolate a plane.
        let expected = BezierSurface::new(&[3.0, 5.0, 4.0, 6.0], [2, 2], unit_box());
        assert_abs_diff_eq!(controls, expected);
    }

    #[test]
    fn bezier_surface_1d_constant() {
        let poly = MultivarPoly::from_monomials([Monomial::new(2.0, [0])]);

        let controls = assert_ok!(poly.bezier_surface(unit_box()));

        // Points interpolate a plane.
        let expected = BezierSurface::new(&[2.0, 2.0], [2], unit_box());
        assert_abs_diff_eq!(controls, expected);
    }

    #[test]
    fn bezier_surface_2d_quadratic() {
        // Unit circle represented by implicit polynomial function:
        //   f(x, y) = -1 + x^2 + y^2
        let circle = MultivarPoly::new().add_monomials([
            Monomial::new(-1.0, [0, 0]),
            Monomial::new(1.0, [2, 0]),
            Monomial::new(1.0, [0, 2]),
        ]);

        let controls = circle.bezier_surface([(-1.0, 1.0), (-1.0, 1.0)]);

        let expected = BezierSurface::new(
            &[
                1.0, -1.0, 1.0, //
                -1.0, -3.0, -1.0, //
                1.0, -1.0, 1.0,
            ],
            [3, 3],
            [(-1.0, 1.0), (-1.0, 1.0)],
        );

        assert_abs_diff_eq!(controls.unwrap(), expected);
    }

    #[test]
    fn eval_circle() {
        // Unit circle represented by implicit polynomial function:
        //   f(x, y) = -1 + x^2 + y^2
        let circle = MultivarPoly::new().add_monomials([
            Monomial {
                coeff: -1.0,
                exp: [0, 0],
            },
            Monomial {
                coeff: 1.0,
                exp: [2, 0],
            },
            Monomial {
                coeff: 1.0,
                exp: [0, 2],
            },
        ]);

        pa::assert_eq!(circle.eval(&[0.0, 0.0]), -1.0);
        pa::assert_eq!(circle.eval(&[-1.0, 0.0]), 0.0);
        pa::assert_eq!(circle.eval(&[1.0, 0.0]), 0.0);
        pa::assert_eq!(circle.eval(&[0.0, -1.0]), 0.0);
        pa::assert_eq!(circle.eval(&[0.0, 1.0]), 0.0);
        pa::assert_eq!(circle.eval(&[1.0, 1.0]), 1.0);
    }

    #[test]
    fn sub_polys_circle_affine() {
        // Unit circle represented by implicit polynomial function:
        //   f(x, y) = -1 + x^2 + y^2
        let circle = MultivarPoly::new().add_monomials([
            Monomial::new(-1.0, [0, 0]),
            Monomial::new(1.0, [2, 0]),
            Monomial::new(1.0, [0, 2]),
        ]);

        // Affine transformation to bigger radius and x offset.
        let affine = [
            [-3.0, 2.0], // x
            [0.0, 2.0],  // y
        ];

        let result_multivar = assert_ok!(circle.sub_polys(&affine));

        let expected = MultivarPoly::from_monomials([
            Monomial::new(8.0, [0, 0]),
            Monomial::new(-12.0, [1, 0]),
            Monomial::new(4.0, [2, 0]),
            Monomial::new(4.0, [0, 2]),
        ]);

        approx::assert_abs_diff_eq!(result_multivar, expected);
    }

    #[test]
    fn sub_polys_linear() {
        crate::test_utils::init_test_logger();

        let line = MultivarPoly::new()
            .add_monomials([Monomial::new(1.0, [1, 0]), Monomial::new(-1.0, [0, 1])]);

        let affine = [
            [-1.0, 2.0], // x
            [-1.0, 2.0], // y
        ];

        let result_multivar = assert_ok!(line.sub_polys(&affine));

        let expected =
            MultivarPoly::from_monomials([Monomial::new(2.0, [1, 0]), Monomial::new(-2.0, [0, 1])]);

        approx::assert_abs_diff_eq!(result_multivar, expected);
    }
}
