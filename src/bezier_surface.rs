use std::{cell::RefCell, ops::IndexMut};

use smallvec::{SmallVec, smallvec};

use crate::{bounding_box::BoundingBox, gridex_excl, vector_promoted::VectorPromoted};

///
/// # Parameters
/// - `NDIM`: Number of dimensions (eg. 2 for a 2D manifold in 3d space)
#[derive(Clone, Debug)]
pub struct BezierSurface<const NDIM: usize> {
    // The coefficients/values of the control points in the Bezier hypersurface hypergrid.
    pub(crate) coeffs: SmallVec<[f64; 32]>,

    // The number of values along each dimension of the Bezier hypergrid.
    pub(crate) grid_size: [usize; NDIM],

    // The parameter space domain for each dimension.
    pub(crate) domain: BoundingBox<NDIM>,

    // The stride lengths for each dimension in the flattened coeffs array.
    // This is used to convert between flat indexes and grid indexes (gridexes).
    strides: [usize; NDIM],

    // Temporary storage for surface evaluation to avoid repeated allocations.
    eval_coeffs: RefCell<Vec<f64>>,
}

impl<const NDIM: usize> BezierSurface<NDIM> {
    /// Create a new Bezier surface with the given control points and grid size.
    ///
    /// # Arguments
    /// - `coeffs`: Flattened array of control point coefficients/values. Organized such that
    ///    each dimension has decreasing stride length. Ending with 1. For example, for a 2D
    ///    surface with grid size [3, 4], the grid indexes would be:
    ///        `[[0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1], ... [2, 3]]`
    /// - `grid_size`: Number of control points along each dimension (== degree + 1)
    /// - `domain`: Parameter space domain for each dimension. Such that the result of
    ///    `self.controls()` will have coordinates in the range defined by `domain`.
    ///
    /// # Panics
    /// Panics if the product of the product of the grid size does not equal coeffs.len().
    pub fn new(
        coeffs_arr: &[f64],
        grid_size: [usize; NDIM],
        domain: impl Into<BoundingBox<NDIM>>,
    ) -> Self {
        assert!(
            grid_size.iter().product::<usize>() == coeffs_arr.len(),
            "Grid size mismatch with number of control points"
        );
        assert!(
            grid_size.iter().all(|&s| s >= 2),
            "Grid size must be at least 2 in each dimension"
        );

        let num_coeffs = grid_size.iter().product();
        let mut coeffs = SmallVec::with_capacity(num_coeffs);
        for &c in coeffs_arr.iter() {
            coeffs.push(c);
        }
        Self {
            eval_coeffs: RefCell::new(vec![0.0; num_coeffs]),
            coeffs,
            grid_size,
            domain: domain.into(),
            strides: Self::calc_strides(grid_size),
        }
    }

    pub fn zeros(grid_size: [usize; NDIM], domain: impl Into<BoundingBox<NDIM>>) -> Self {
        assert!(
            grid_size.iter().all(|&s| s >= 2),
            "Grid size must be at least 2 in each dimension"
        );

        let num_coeffs = grid_size.iter().product();
        Self {
            coeffs: smallvec![0.0; num_coeffs],
            eval_coeffs: RefCell::new(vec![0.0; num_coeffs]),
            grid_size,
            domain: domain.into(),
            strides: Self::calc_strides(grid_size),
        }
    }

    pub fn set_coeff(&mut self, coord: [usize; NDIM], coeff: f64) {
        let flat_index = self.flat_index(coord);
        self.coeffs[flat_index] = coeff;
    }

    pub fn eval(&self, param_coord: [f64; NDIM]) -> f64 {
        // Re-using eval_coeffs leads to significant speed-ups by avoiding allocations.
        let mut eval_coeffs = self.eval_coeffs.borrow_mut();
        eval_coeffs.clone_from_slice(&self.coeffs);

        let mut sub_grid_size = self.grid_size;

        for d in 0..NDIM {
            let (min, max) = self.domain[d];
            let t = (param_coord[d] - min) / (max - min);
            let stride = self.strides[d];
            let size = self.grid_size[d];

            // With each loop, we are reducing the dimension by 1.
            // For example, a 2D manifold evaluates to a 1D curve in the first pass.
            // And then we run de Casteljau again on that curve.
            // Setting the current d grid_size to one, skips over all the intermediate points
            // calculated in the previous iteration.
            // Side note: we want these intermediate points in subsection(), but not in eval().
            sub_grid_size[d] = 1;

            for start_gridex in gridex_excl(sub_grid_size) {
                let start = self.flat_index(start_gridex);
                de_casteljau_upper_nd(&mut (*eval_coeffs), start, stride, size, t);
            }
        }

        eval_coeffs[0]
    }

    pub fn eval_grad(&self, param_coord: [f64; NDIM]) -> (f64, [f64; NDIM]) {
        // Re-using eval_coeffs leads to significant speed-ups by avoiding allocations.
        let mut eval_coeffs = self.eval_coeffs.borrow_mut();
        eval_coeffs.clone_from_slice(&self.coeffs);

        let mut sub_grid_size = self.grid_size;
        let mut corner_gridex = [0usize; NDIM];
        // let mut ts = [0.0f64; NDIM];

        for d in 0..NDIM {
            let (min, max) = self.domain[d];
            let t = if max == min {
                0.5
            } else {
                (param_coord[d] - min) / (max - min)
            };

            let stride = self.strides[d];
            let size = self.grid_size[d];

            sub_grid_size[d] = 1;

            // Use a threshold of 0.95 to bias branch prediction.
            if t < 0.95 {
                for mut start_gridex in gridex_excl(sub_grid_size) {
                    for i in 0..d {
                        if corner_gridex[i] > 0 {
                            start_gridex[i] = corner_gridex[i] - start_gridex[i];
                        }
                    }
                    let start = self.flat_index(start_gridex);
                    de_casteljau_upper_nd(&mut (*eval_coeffs), start, stride, size, t);
                }
            } else {
                // Use the lower de Casteljau algorithm to avoid numerical instability near t=1.
                for mut start_gridex in gridex_excl(sub_grid_size) {
                    for i in 0..d {
                        if corner_gridex[i] > 0 {
                            start_gridex[i] = corner_gridex[i] - start_gridex[i];
                        }
                    }
                    let start = self.flat_index(start_gridex);
                    de_casteljau_lower_nd(&mut (*eval_coeffs), start, stride, size, t);
                }
                corner_gridex[d] = size - 1;
            }

            // Ensures we have the correct intermediate values for derivative calculation
            sub_grid_size[d] = self.grid_size[d].min(2);
        }

        let corner_index = self.flat_index(corner_gridex);
        let coeff = eval_coeffs[corner_index];
        log::debug!("  coeffs: {:?}", eval_coeffs);
        log::debug!("  corner_gridex: {:?}", corner_gridex);
        log::debug!("  corner_index: {}", self.flat_index(corner_gridex));

        let mut derivs = [0.0f64; NDIM];
        for d in 0..NDIM {
            let (min, max) = self.domain[d];
            let length = max - min;
            if length == 0.0 {
                // This can be arbitrarily non-zero. Since length only collapses to zero
                // when that dimension is entirely intersecting x_NDIM = 0.
                // A gradient of 1.0 helps improve the condition of hyperplane intersections.
                derivs[d] = 1.0;
                continue;
            }

            let t_length = param_coord[d] - min;
            let scale_length = (self.grid_size[d] - 1) as f64;
            log::debug!("  d={} t={} scale={}", d, t_length, scale_length);
            log::debug!("  length={}", length);

            let mut next_gridex = corner_gridex;
            let (delta, w_length) = if corner_gridex[d] == 0 {
                // Forward difference
                next_gridex[d] += 1;
                let next_ind = self.flat_index(next_gridex);
                let delta = eval_coeffs[next_ind] - coeff;
                (delta, length - t_length)
            } else {
                // Backward difference
                next_gridex[d] -= 1;
                let next_ind = self.flat_index(next_gridex);
                let delta = coeff - eval_coeffs[next_ind];
                (delta, t_length)
            };
            debug_assert!(w_length != 0.0);
            derivs[d] = delta * scale_length / w_length;
        }
        (coeff, derivs)
    }

    pub fn numeric_gradient(&self, param_coord: [f64; NDIM], epsilon: f64) -> [f64; NDIM] {
        let mut gradient = [0.0f64; NDIM];
        for d in 0..NDIM {
            let mut param_plus = param_coord;
            param_plus[d] += epsilon;
            let f_plus = self.eval(param_plus);

            let mut param_minus = param_coord;
            param_minus[d] -= epsilon;
            let f_minus = self.eval(param_minus);

            gradient[d] = (f_plus - f_minus) / (2.0 * epsilon);
        }
        gradient
    }

    pub fn num_controls(&self) -> usize {
        self.grid_size.iter().product()
    }

    pub fn controls(&self) -> impl Iterator<Item = ControlPoint<NDIM>> {
        (0..self.num_controls()).map(move |i| {
            // Build coordinate and control point
            let index = self.unflat_index(i);
            let mut coord = [0.0; NDIM];
            for d in 0..NDIM {
                let t = index[d] as f64 / (self.grid_size[d] - 1) as f64;
                coord[d] = (1.0 - t) * self.domain[d].0 + t * self.domain[d].1;
            }

            ControlPoint::new(coord, self.coeffs[i], index)
        })
    }

    /// Apply de Casteljau subdivision to the Bezier surface to extract the subsection
    /// defined by the given bounding box in parameter space.
    pub fn subsection(mut self, bbox: impl Into<BoundingBox<NDIM>>) -> Self {
        let bbox = bbox.into();

        for d in 0..NDIM {
            let (min, max) = self.domain[d];
            let (t0, t1) = bbox[d];
            let t0 = (t0 - min) / (max - min);
            let t1 = (t1 - min) / (max - min);

            let stride = self.strides[d];
            let size = self.grid_size[d];

            // Setting the current dimension to zero for the sub-iteration
            // of the grid "plane" at gridex[d] == 0.
            let mut sub_grid_size = self.grid_size;
            sub_grid_size[d] = 0;

            for start_gridex in gridex_excl(sub_grid_size) {
                // There might be a faster way to calculate the start indexes,
                // But this seems easy to understand and reason about.
                let start = self.flat_index(start_gridex);
                de_casteljau_lower_nd(&mut self.coeffs, start, stride, size, t1);
                // Don't forget to re-normalize t0
                de_casteljau_upper_nd(&mut self.coeffs, start, stride, size, t0 / t1);
            }
        }
        self.domain = bbox;
        self
    }

    fn flat_index(&self, coord: [usize; NDIM]) -> usize {
        let mut index = 0;
        for i in 0..NDIM {
            index += coord[i] * self.strides[i];
        }
        index
    }

    fn unflat_index(&self, index: usize) -> [usize; NDIM] {
        let mut coord = [0usize; NDIM];
        let mut remainder = index;
        for d in (0..NDIM).rev() {
            coord[d] = remainder % self.grid_size[d];
            remainder /= self.grid_size[d];
        }
        coord
    }

    fn calc_strides(grid_size: [usize; NDIM]) -> [usize; NDIM] {
        let mut strides = [1usize; NDIM];
        for i in (0..NDIM - 1).rev() {
            strides[i] = strides[i + 1] * grid_size[i + 1];
        }
        strides
    }
}

/// Represents a Bezier control point in (D+1)-dimensional space
///
/// The first D components are the position, and the last component is the surface "height".
/// For example, a a surface in 3D space would have D=2, with (x,y) as the first two components
/// and z as the last component.
#[derive(Debug, Clone)]
pub struct ControlPoint<const D: usize> {
    /// The coordinate of this control point in "parameter" space, and its coefficient/value.
    pub point: VectorPromoted<D>,

    /// The integer index of this control point in the Bezier grid
    pub index: [usize; D],
}

impl<const D: usize> ControlPoint<D> {
    pub fn new(coord: [f64; D], coeff: f64, index: [usize; D]) -> Self {
        Self {
            point: VectorPromoted::new(coord, coeff),
            index,
        }
    }

    pub fn axis_value(&self, axis: usize) -> f64 {
        if axis < D {
            self.point.base()[axis]
        } else if axis == D {
            self.point.extra()
        } else {
            panic!("Axis index out of bounds for ControlPoint")
        }
    }

    pub fn coeff(&self) -> f64 {
        self.point.extra()
    }
}

fn de_casteljau_lower_nd(
    control_points: &mut impl IndexMut<usize, Output = f64>,
    start: usize,
    stride: usize,
    size: usize,
    t: f64,
) {
    let s = 1.0 - t;
    for n in 1..size {
        for i in (n..size).rev() {
            let j = start + i * stride;
            let j_prev = j - stride;
            control_points[j] = s * control_points[j_prev] + t * control_points[j];
        }
    }
}

fn de_casteljau_upper_nd(
    control_points: &mut impl IndexMut<usize, Output = f64>,
    start: usize,
    stride: usize,
    size: usize,
    t: f64,
) {
    let s = 1.0 - t;
    for k in 1..size {
        for i in 0..(size - k) {
            let j = start + i * stride;
            let j_next = j + stride;
            control_points[j] = s * control_points[j] + t * control_points[j_next];
        }
    }
}

// Equality for testing purposes
// Consider using approx crate instead (https://docs.rs/approx/latest/approx/)
#[cfg(test)]
impl<const NDIM: usize> approx::AbsDiffEq for BezierSurface<NDIM> {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-10
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        if self.grid_size != other.grid_size {
            return false;
        }

        if !self.domain.abs_diff_eq(&other.domain, epsilon) {
            return false;
        }

        let num_controls: usize = self.grid_size.iter().product();
        (0..num_controls).all(|i| (self.coeffs[i] - other.coeffs[i]).abs() <= epsilon)
    }
}

#[cfg(test)]
impl<const NDIM: usize> PartialEq for BezierSurface<NDIM> {
    fn eq(&self, other: &Self) -> bool {
        approx::AbsDiffEq::abs_diff_eq(self, other, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use assertables::assert_lt;

    use crate::{
        binomial_coefficient,
        test_utils::{linspace, unit_box},
    };

    use super::*;

    #[test]
    fn index_mapping() {
        const NDIM: usize = 3;
        const GRID_SIZE: [usize; NDIM] = [4, 3, 2];
        const MAX_CONTROLS: usize = 4 * 3 * 2;
        let bezier = BezierSurface::zeros(GRID_SIZE, unit_box());

        let mut coords = HashSet::new();
        for i in 0..MAX_CONTROLS {
            let coord = bezier.unflat_index(i);
            assert_lt!(coord[0], 4);
            assert_lt!(coord[1], 3);
            assert_lt!(coord[2], 2);
            coords.insert(coord);

            let flat_index = bezier.flat_index(coord);
            assert_eq!(i, flat_index);
        }

        // Ensure all coordinates were visited
        assert_eq!(coords.len(), MAX_CONTROLS);
    }

    #[test]
    fn surface_eval_quadratic() {
        crate::test_utils::init_test_logger();

        let grid_size = [3, 3];
        // Control points for the surface z = x^2 + y^2 - 1
        let coeffs = [
            1.0, -1.0, 1.0, //
            -1.0, -3.0, -1.0, //
            1.0, -1.0, 1.0,
        ];
        let bezier = BezierSurface::new(&coeffs, grid_size, [(-1.0, 1.0), (-1.0, 1.0)]);

        for x in linspace(-1.0, 1.0, 5) {
            for y in linspace(-1.0, 1.0, 5) {
                let z = x * x + y * y - 1.0;
                let eval_point = bezier.eval([x, y]);
                assert_lt!((eval_point - z).abs(), 1e-10);
            }
        }
    }

    #[test]
    fn surface_eval_grad_quadratic_1d() {
        crate::test_utils::init_test_logger();

        let grid_size = [3];
        // Control points for the surface f(x) = x^2 - 1
        let coeffs = [0.0, -2.0, 0.0];
        let bezier = BezierSurface::new(&coeffs, grid_size, [(-1.0, 1.0)]);

        for x in linspace(-1.0, 1.0, 5) {
            let f = x * x - 1.0;
            let df = [2.0 * x];
            let (eval_f, eval_df) = bezier.eval_grad([x]);
            assert_lt!((eval_f - f).abs(), 1e-10);
            assert_lt!(
                distance(eval_df, df),
                1e-10,
                "eval_df: {:?}, df: {:?}",
                eval_df,
                df
            );
        }
    }

    #[test]
    fn surface_eval_grad_cubic_1d() {
        crate::test_utils::init_test_logger();

        let grid_size = [4];
        // Control points for the surface f(x) = ?
        let coeffs = [2.0, -2.0, 1.0, -10.0];
        let bezier = BezierSurface::new(&coeffs, grid_size, [(-2.0, 1.0)]);

        for x in linspace(-2.0, 1.0, 20) {
            let f = bezier.eval([x]);
            let df = {
                // Numerical derivative, not very accurate.
                const EPS: f64 = 1e-8;
                let x1 = x - EPS;
                let x2 = x + EPS;
                [(bezier.eval([x2]) - bezier.eval([x1])) / (2.0 * EPS)]
            };

            let (eval_f, eval_df) = bezier.eval_grad([x]);
            assert_lt!((eval_f - f).abs(), 1e-10);
            assert_lt!(
                distance(eval_df, df),
                1e-6,
                "eval_df: {:?}, df: {:?}",
                eval_df,
                df
            );
        }
    }

    #[test]
    fn surface_eval_grad_2d() {
        crate::test_utils::init_test_logger();

        let grid_size = [3, 4];
        // Control points for the surface f(x, y) = ?
        let coeffs = [
            2.0, -2.0, 1.0, -10.0, //
            5.0, 0.0, -6.0, 4.0, //
            -1.0, 3.0, -2.0, 1.0, //
        ];
        let bezier = BezierSurface::new(&coeffs, grid_size, [(-2.0, 1.0), (-1.5, 1.0)]);
        // let bezier = BezierSurface::new(&coeffs, grid_size, [(0.0, 1.0), (0.0, 1.0)]);

        // for x in linspace(0.0, 1.0, 20) {
        //     for y in linspace(0.0, 1.0, 20) {
        for x in linspace(-2.0, 1.0, 20) {
            for y in linspace(-1.5, 1.0, 20) {
                log::info!("Testing point ({}, {})", x, y);
                let f = bezier.eval([x, y]);
                let df = bezier.numeric_gradient([x, y], 1e-8);

                let (eval_f, eval_df) = bezier.eval_grad([x, y]);
                assert_lt!((eval_f - f).abs(), 1e-10, "eval_f: {}, f: {}", eval_f, f);
                assert_lt!(
                    distance(eval_df, df),
                    1e-6,
                    "eval_df: {:?}, df: {:?}",
                    eval_df,
                    df
                );
            }
        }
    }

    #[test]
    fn surface_subsection_quadratic() {
        let grid_size = [3, 3];
        // Control points for the surface z = x^2 + y^2 - 1
        let coeffs = [
            1.0, -1.0, 1.0, //
            -1.0, -3.0, -1.0, //
            1.0, -1.0, 1.0,
        ];
        let bezier = BezierSurface::new(&coeffs, grid_size, [(-1.0, 1.0), (-1.0, 1.0)]);

        // Set control points to their flat index for easy verification
        let sub_domain = [(0.25, 0.75), (0.0, 0.5)];
        let sub_bezier = bezier.subsection(sub_domain);

        for d in 0..2 {
            assert_lt!(distance(sub_bezier.domain[d], sub_domain[d]), 1e-10);
        }

        for x in linspace(0.25, 0.75, 5) {
            for y in linspace(0.0, 0.5, 5) {
                let z = x * x + y * y - 1.0;

                let z_eval = sub_bezier.eval([x, y]);
                assert_lt!((z_eval - z).abs(), 1e-10);
            }
        }
    }

    //
    // Helper func testing
    //

    #[test]
    fn de_casteljau_lower_linear() {
        let control_points = [[0.0], [1.0], [2.0]];
        let t0 = 0.5;
        let new_points = de_casteljau_lower(control_points, t0);

        for i in 0..20 {
            let t_new = i as f64 / 20.0;
            let t_old = t_new * t0;
            let pt_new = eval_bezier(&new_points, t_new);
            let pt_old = eval_bezier(&control_points, t_old);
            assert_lt!(distance(pt_new, pt_old), 1e-10);
        }
    }

    #[test]
    fn de_casteljau_lower_quadratic() {
        let control_points = [[-1.0, 1.0], [0.0, -1.0], [1.0, 1.0]];
        let t0 = 0.45;
        let new_points = de_casteljau_lower(control_points, t0);

        for i in 0..20 {
            let t_new = i as f64 / 20.0;
            let t_old = t_new * t0;
            let pt_new = eval_bezier(&new_points, t_new);
            let pt_old = eval_bezier(&control_points, t_old);
            assert_lt!(distance(pt_new, pt_old), 1e-10);
        }
    }

    #[test]
    fn de_casteljau_upper_linear() {
        let control_points = [[0.0], [1.0], [2.0]];
        let t0 = 0.5;
        let new_points = de_casteljau_upper(control_points, t0);

        for i in 0..20 {
            let t_new = i as f64 / 20.0;
            let t_old = t0 + t_new * (1.0 - t0);
            let pt_new = eval_bezier(&new_points, t_new);
            let pt_old = eval_bezier(&control_points, t_old);
            assert_lt!(distance(pt_new, pt_old), 1e-10);
        }
    }

    #[test]
    fn de_casteljau_upper_quadratic() {
        let control_points = [[-1.0, 1.0], [0.0, -1.0], [1.0, 1.0]];
        let t0 = 0.45;
        let new_points = de_casteljau_upper(control_points, t0);

        for i in 0..20 {
            let t_new = i as f64 / 20.0;
            let t_old = t0 + t_new * (1.0 - t0);
            let pt_new = eval_bezier(&new_points, t_new);
            let pt_old = eval_bezier(&control_points, t_old);
            assert_lt!(distance(pt_new, pt_old), 1e-10);
        }
    }

    #[test]
    fn de_casteljau_subsection_linear() {
        let control_points = [[0.0], [1.0], [2.0]];
        let t0 = 0.5;
        let t1 = 0.8;
        let new_points = de_casteljau_subsection(control_points, t0, t1);

        dbg!(&new_points);

        for i in 0..20 {
            let t_new = i as f64 / 20.0;
            let t_old = t0 + t_new * (t1 - t0);
            let pt_new = eval_bezier(&new_points, t_new);
            let pt_old = eval_bezier(&control_points, t_old);
            assert_lt!(distance(pt_new, pt_old), 1e-10);
        }
    }

    #[test]
    fn de_casteljau_subsection_quadratic() {
        let control_points = [[-1.0, 1.0], [0.0, -1.0], [1.0, 1.0]];
        let t0 = 0.45;
        let t1 = 0.75;
        let new_points = de_casteljau_subsection(control_points, t0, t1);

        for i in 0..20 {
            let t_new = i as f64 / 20.0;
            let t_old = t0 + t_new * (t1 - t0);
            let pt_new = eval_bezier(&new_points, t_new);
            let pt_old = eval_bezier(&control_points, t_old);
            assert_lt!(distance(pt_new, pt_old), 1e-10);
        }
    }

    #[test]
    fn bezier_eval_linear() {
        let control_points = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]];

        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let result = eval_bezier(&control_points, t);
            assert_lt!(distance(result, [3.0 * t, 3.0 * t]), 1e-10);
        }
    }

    #[test]
    fn bezier_eval_quadratic() {
        let control_points = [[-1.0, 1.0], [0.0, -1.0], [1.0, 1.0]];

        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let result = eval_bezier(&control_points, t);

            let x = 2.0 * t - 1.0;
            assert_lt!(distance(result, [x, x * x]), 1e-10);
        }
    }

    fn distance<const D: usize>(a: impl Into<[f64; D]>, b: impl Into<[f64; D]>) -> f64 {
        let a = a.into();
        let b = b.into();
        let mut sum = 0.0;
        for i in 0..D {
            let diff = a[i] - b[i];
            sum += diff * diff;
        }
        sum.sqrt()
    }

    //
    // Play and testing
    //
    //
    // 1D De Casteljau algorithms
    //

    fn de_casteljau_subsection<const D: usize, const N: usize>(
        mut control_points: [[f64; D]; N],
        t0: f64,
        t1: f64,
    ) -> [[f64; D]; N] {
        assert!(t0 < t1, "t0 must be less than t1");
        assert!(t0 >= 0.0 && t1 <= 1.0, "t0 and t1 must be in [0, 1]");

        // NOTE: There's probably a way to do this that keeps good data locality/cache coherence.
        // I thiiiink I might just require the smooshing of the two inner loops next to each other.
        // For now, just do two separate passes.
        control_points = de_casteljau_lower(control_points, t1);
        // Don't forget to re-normalize t0
        de_casteljau_upper(control_points, t0 / t1)
    }

    fn _de_casteljau_split<const D: usize, const N: usize>(
        control_points: [[f64; D]; N],
        t: f64,
    ) -> ([[f64; D]; N], [[f64; D]; N]) {
        (
            de_casteljau_lower(control_points, t),
            de_casteljau_upper(control_points, t),
        )
    }

    fn de_casteljau_lower<const D: usize, const N: usize>(
        mut control_points: [[f64; D]; N],
        t: f64,
    ) -> [[f64; D]; N] {
        let s = 1.0 - t;
        for n in 1..N {
            for i in (n..N).rev() {
                for d in 0..D {
                    control_points[i][d] = s * control_points[i - 1][d] + t * control_points[i][d];
                }
            }
        }
        control_points
    }

    fn de_casteljau_upper<const D: usize, const N: usize>(
        mut control_points: [[f64; D]; N],
        t: f64,
    ) -> [[f64; D]; N] {
        let s = 1.0 - t;
        for n in 1..N {
            for i in 0..(N - n) {
                for d in 0..D {
                    control_points[i][d] = s * control_points[i][d] + t * control_points[i + 1][d];
                }
            }
        }
        control_points
    }

    fn eval_bezier<const D: usize, const N: usize>(
        control_points: &[[f64; D]; N],
        t: f64,
    ) -> [f64; D] {
        let s = 1.0 - t;
        let mut result = [0.0; D];
        for i in 0..N {
            let coeff = binomial_coefficient((N - 1) as u8, i as u8) as f64
                * s.powi((N - 1 - i) as i32)
                * t.powi(i as i32);
            for d in 0..D {
                result[d] += coeff * control_points[i][d];
            }
        }
        result
    }
}
