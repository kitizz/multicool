use arrayvec::ArrayVec;
use fmtastic::{Subscript, Superscript};
use smallvec::SmallVec;

use crate::{MulticoolError, MultivarPoly, gridex_excl};

///
/// V: Number of variables
#[derive(Clone)]
pub struct Monomial<const V: usize = 4> {
    /// Coefficient of the monomial
    ///
    /// For example, in 3.5 * x^2 * y^1 * z^0, the coefficient is 3.5
    pub(crate) coeff: f64,

    /// Exponents for each variable
    ///
    /// For example, if variables are ordered (x, y, z), and the monomial is 3 * x^2 * y^1 * z^0,
    /// then exp would be [2, 1, 0]
    pub(crate) exp: [u8; V],
}

impl<const V: usize> core::fmt::Debug for Monomial<V> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        if self.exp.iter().any(|&e| e != 0) {
            if self.coeff != 1.0 {
                write!(f, "{} * ", self.coeff)?;
            }
        } else {
            write!(f, "{}", self.coeff)?;
            return Ok(());
        }

        for i in 0..V {
            if self.exp[i] == 0 {
                continue;
            }
            write!(f, "x{}{} ", Subscript(i), Superscript(self.exp[i]))?;
        }
        Ok(())
    }
}

// Type alias for a vector of polynomial coefficients
// Use SmallVec to reduce heap allocations for small polynomials
type PolyVec = SmallVec<[f64; 16]>;

impl<const V: usize> Monomial<V> {
    pub fn new(coeff: f64, exp: [u8; V]) -> Self {
        Self { coeff, exp }
    }

    pub fn eval(&self, vars: &[f64; V]) -> f64 {
        let mut result = self.coeff;
        for i in 0..V {
            result *= vars[i].powi(self.exp[i] as i32);
        }
        result
    }

    /// Check if all exponents in self are less than or equal to those in other
    /// See Section 2.1 in [1]
    ///
    /// [1] https://adrianbowyer.com/Publications/berchtold2000.pdf
    pub fn exp_all_le(&self, other_exp: &[u8; V]) -> bool {
        for i in 0..V {
            if self.exp[i] > other_exp[i] {
                return false;
            }
        }
        true
    }

    /// Substitute each variable in this monomial with the given polynomial.
    ///
    /// Note that each polynomial's variable is separate from the others,
    /// and separate from the monomial's variables.
    ///
    /// Example:
    /// For monomial f(x, y) = 3 * x^2 * y and vars [
    ///     x = 1 + 2u,
    ///     y = 3 + 6v^2,
    /// ], the code looks like:
    /// ```rust
    /// use multicool::{Monomial, MultivarPoly};
    ///
    /// let mon = Monomial::new(3.0, [2, 1]);
    /// let polys = [
    ///     [1.0, 2.0, 0.0],
    ///     [3.0, 0.0, 6.0],
    /// ];
    /// let result_multivar = mon.sub_polys(&polys).unwrap();
    ///
    /// let expected = MultivarPoly::stack_from_monomials([
    ///     Monomial::new(9.0, [0, 0]),
    ///     Monomial::new(18.0, [0, 2]),
    ///     Monomial::new(36.0, [1, 0]),
    ///     Monomial::new(72.0, [1, 2]),
    ///     Monomial::new(36.0, [2, 0]),
    ///     Monomial::new(72.0, [2, 2]),
    /// ]);
    /// assert_eq!(result_multivar, expected);
    /// ```
    /// Results in the multivariate polynomial:
    /// g(u, v) = 9 + 18(v^2) + 36(u) + 72(u)(v^2) + 36(u^2) + 72(u^2)(v^2)
    ///
    /// # Parameters
    /// - `N`: Maximum number of coefficients in the input polynomials
    ///
    pub fn sub_polys<const N: usize>(
        &self,
        rhs_polys: &[[f64; N]; V],
    ) -> Result<MultivarPoly<V>, MulticoolError> {
        let mut var_polys = ArrayVec::<PolyVec, V>::new();
        let mut num_coeffs = [0u8; V];

        for i in 0..V {
            let var_pow = poly_pow(&rhs_polys[i], self.exp[i])?;
            num_coeffs[i] = (var_pow.len()) as u8;
            var_polys.push(var_pow);
        }

        // Expand the product of all var_polys into its monomial terms.
        let mut monomials = smallvec::smallvec![];
        for exp in gridex_excl(num_coeffs) {
            let mut sub_coeff = 1.0;
            for i in 0..V {
                sub_coeff *= var_polys[i][exp[i] as usize];
            }
            if sub_coeff == 0.0 {
                // Don't waste space on zero coefficient terms
                continue;
            }
            let term = Monomial::new(self.coeff * sub_coeff, exp);
            monomials.push(term);
        }
        monomials.sort();
        Ok(MultivarPoly { terms: monomials })
    }
}

/// Find the size of this array if trailing zeros are ignored.
#[allow(dead_code)] // For future use
fn size<const N: usize>(coeffs: &[f64; N]) -> usize {
    for idx in (N - 1)..=0 {
        if coeffs[idx] != 0.0 {
            return idx + 1;
        }
    }
    N
}

fn poly_pow<const N: usize>(coeffs: &[f64; N], exp: u8) -> Result<PolyVec, MulticoolError> {
    if exp == 0 {
        let mut result = smallvec::smallvec![];
        result.push(1.0);
        return Ok(result);
    }

    let mut result = smallvec::smallvec![];
    for i in 0..N {
        result.push(coeffs[i]);
    }

    for _ in 1..exp {
        result = convolve(result, coeffs)?;
    }
    Ok(result)
}

fn convolve<const N: usize>(a: PolyVec, b: &[f64; N]) -> Result<PolyVec, MulticoolError> {
    let mut result = smallvec::smallvec![0.0; a.len() + b.len() - 1];

    for i in 0..a.len() {
        for j in 0..b.len() {
            result[i + j] += a[i] * b[j];
        }
    }
    Ok(result)
}

use std::cmp::Ordering;

impl<const V: usize> Ord for Monomial<V> {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.exp.cmp(&other.exp) {
            Ordering::Equal => self
                .coeff
                .partial_cmp(&other.coeff)
                .unwrap_or(Ordering::Equal),
            ord => ord,
        }
    }
}

impl<const V: usize> PartialOrd for Monomial<V> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.exp.cmp(&other.exp))
    }
}

impl<const V: usize> PartialEq for Monomial<V> {
    fn eq(&self, other: &Self) -> bool {
        self.exp == other.exp && self.coeff == other.coeff
    }
}

impl<const V: usize> Eq for Monomial<V> {}

// Implement approximate equality for testing purposes
#[cfg(test)]
impl<const V: usize> approx::AbsDiffEq for Monomial<V> {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.exp == other.exp && (self.coeff - other.coeff).abs() < epsilon
    }
}

#[cfg(test)]
mod tests {
    use assertables::assert_ok;

    use super::*;

    #[test]
    fn build_monomial() {
        let monomial = Monomial {
            coeff: 3.0,
            exp: [2, 1, 0, 0],
        };
        assert_eq!(monomial.coeff, 3.0);
        assert_eq!(monomial.exp, [2, 1, 0, 0]);
    }

    #[test]
    fn exp_all_le_true() {
        let other_exp = [1, 2, 3];

        for i in 0..=1 {
            for j in 0..=2 {
                for k in 0..=3 {
                    let monomial = Monomial {
                        coeff: 2.0,
                        exp: [i, j, k],
                    };
                    assert!(monomial.exp_all_le(&other_exp));
                }
            }
        }
    }

    #[test]
    fn exp_all_le_false() {
        let other_exp = [1, 2, 3];

        for i in 0..=1 {
            for j in 0..=2 {
                for k in 4..8 {
                    let monomial = Monomial {
                        coeff: 2.0,
                        exp: [i, j, k],
                    };
                    assert!(!monomial.exp_all_le(&other_exp));
                }
            }
        }

        for i in 0..=1 {
            for j in 3..5 {
                for k in 0..=3 {
                    let monomial = Monomial {
                        coeff: 2.0,
                        exp: [i, j, k],
                    };
                    assert!(!monomial.exp_all_le(&other_exp));
                }
            }
        }

        for i in 2..4 {
            for j in 0..=2 {
                for k in 0..=3 {
                    let monomial = Monomial {
                        coeff: 2.0,
                        exp: [i, j, k],
                    };
                    assert!(!monomial.exp_all_le(&other_exp));
                }
            }
        }
    }

    #[test]
    fn sub_polys_identity() {
        let monomial = Monomial {
            coeff: 2.0,
            exp: [2, 3, 1],
        };
        let vars = [
            [0.0, 1.0], // x
            [0.0, 1.0], // y
            [0.0, 1.0], // z
        ];
        let result_multivar = assert_ok!(monomial.sub_polys(&vars));

        let expected = MultivarPoly::from_monomials([Monomial {
            coeff: 2.0,
            exp: [2, 3, 1],
        }]);
        assert_eq!(result_multivar, expected);
    }

    #[test]
    fn sub_polys_constant() {
        let monomial = Monomial {
            coeff: 5.0,
            exp: [2],
        };
        let vars = [
            [3.0], // x
        ];
        let result_multivar = assert_ok!(monomial.sub_polys(&vars));

        let expected = MultivarPoly::from_monomials([Monomial {
            coeff: 45.0,
            exp: [0],
        }]);
        assert_eq!(result_multivar, expected);
    }
}
