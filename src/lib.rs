// TODO:
// - [ ] Benchmark methods for common CAD/CAM problems
// - [ ] Add methods for building a system from explicit and implicit hyper-surfaces
//   - [ ] Eg. p(t) = ... and q(x,y,z) = 0  =>  q(p(t)) = 0
//   - [ ] Benchmark overhead

mod bezier_surface;
mod binomial;
mod bounding_box;
mod gridex;
mod monomial;
mod multivar_poly;
mod multivar_poly_system;
mod vector_promoted;

#[cfg(test)]
mod test_utils;

pub use bezier_surface::*;
pub use binomial::*;
pub use gridex::*;
pub use monomial::Monomial;
pub use multivar_poly::*;
pub use multivar_poly_system::*;

use snafu::prelude::*;

#[derive(Debug, Snafu)]
pub enum MulticoolError {
    #[snafu(display("Exceeded maximum number of terms allowed: {} > {}", actual, max))]
    ExceededMaxTerms { max: usize, actual: usize },

    #[snafu(display(
        "Underdefined system: {} variables but only {} equations",
        num_vars,
        num_equations
    ))]
    UnderdefinedSystem {
        num_vars: usize,
        num_equations: usize,
    },

    #[snafu(display("Algorithm error (bug in library): {}", message))]
    AlgorithmError { message: String },
}
