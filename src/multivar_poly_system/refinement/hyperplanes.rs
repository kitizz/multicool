use crate::{
    BezierSurface, MulticoolError, bounding_box::BoundingBox, vector_promoted::VectorPromoted,
};

pub struct HyperplanesRefiner<const D: usize> {
    bounding_hyperplanes: Vec<BoundingHyperplanes<D>>,
}

impl<const D: usize> HyperplanesRefiner<D> {
    pub fn new(num_surfaces: usize) -> Self {
        Self {
            bounding_hyperplanes: vec![BoundingHyperplanes::default(); num_surfaces],
        }
    }

    pub fn refine_regions(
        &mut self,
        surfaces: &[BezierSurface<D>],
        valid_regions: &mut BoundingBox<D>,
    ) -> Result<(), MulticoolError> {
        // Calculate bounding hyperplanes for each surface
        self.bounding_hyperplanes
            .resize(surfaces.len(), BoundingHyperplanes::default());

        for (i, surface) in surfaces.iter().enumerate() {
            // The following gives us the normal of the gradient at the region center.
            // We use this normal as the normal for the parallel bounding hyperplanes.
            let (_, gradient) = surface.eval_grad(valid_regions.center());
            let normal = VectorPromoted::new(gradient, -1.0).normalized();

            // Find the control points with the min and max dot product with the gradient
            let mut bmin = f64::INFINITY;
            let mut bmax = f64::NEG_INFINITY;
            for control in surface.controls() {
                // Project each control point onto the hyperplanes normal.
                let projected = normal.dot(&control.point);
                bmin = bmin.min(projected);
                bmax = bmax.max(projected);
            }
            let bounding_hyperplane = &mut self.bounding_hyperplanes[i];
            bounding_hyperplane.normal = normal;
            bounding_hyperplane.bmin = bmin;
            bounding_hyperplane.bmax = bmax;
        }

        let hyperplane_bounds = intersection_region(&self.bounding_hyperplanes);

        match hyperplane_bounds {
            None => {}
            Some(hyperplane_bounds) => {
                for d in 0..D {
                    valid_regions[d].0 = valid_regions[d].0.max(hyperplane_bounds[d].0);
                    valid_regions[d].1 = valid_regions[d].1.min(hyperplane_bounds[d].1);
                }
            }
        };

        Ok(())
    }
}

/// BoundingHyperplanes represent two parallel hyperplanes that bound some space.
///
/// The space is defined by:
///   { v | bmin <= normal • v <= bmax }
/// where (v, normal) ∈ R^D  and  bmin, bmax ∈ R
#[derive(Debug, Clone, Default)]
pub struct BoundingHyperplanes<const D: usize> {
    normal: VectorPromoted<D>,
    bmin: f64,
    bmax: f64,
}

/// Compute the bounding box that contains the intersection of BoundingHyperplanes
fn intersection_region<const D: usize>(
    planes: &[BoundingHyperplanes<D>],
) -> Option<BoundingBox<D>> {
    match D {
        1 => intersection_region_1(planes),
        2 => intersection_region_2(planes),
        3 => intersection_region_3(planes),
        4 => intersection_region_4(planes),
        _ => None,
    }
}

macro_rules! impl_intersection_region {
($fn_name: ident, $D:expr) => {
    fn $fn_name<const D: usize>(planes: &[BoundingHyperplanes<D>]) -> Option<BoundingBox<D>> {
        use nalgebra as na;
        // We want to find the intersection point of the mid-hyperplanes.
        // Then the bounding region is computed around this central point.
        debug_assert!(planes.len() == $D);

        let lhs = na::SMatrix::<f64, $D, $D>::from_fn(|r, c| planes[r].normal.base()[c]);
        let mut rhs = na::SVector::<f64, $D>::from_fn(|r, _c| {
            0.5 * (planes[r].bmin + planes[r].bmax)
        });

        let Some(center_point) = lhs.lu().solve(&rhs) else {
            return None;
        };

        // Compute the bounding box around the center point
        let mut bounds_extents = [0.0; $D];
        for d in 0..$D {
            let half_plane_disp = 0.5 * (planes[d].bmax - planes[d].bmin);
            rhs.fill(0.0);
            rhs[d] = half_plane_disp;
            let has_solution = lhs.lu().solve_mut(&mut rhs);
            if !has_solution {
                panic!(
                    "Bounding hyperplanes are not linearly independent. Coincident curves?\n{:#?}\n{:?}", lhs, rhs
                );
            }
            for dim in 0..$D {
                bounds_extents[dim] += rhs[dim].abs();
            }
        }

        let mut bbox = [(0.0, 0.0); D];
        for d in 0..$D {
            bbox[d] = (
                center_point[d] - bounds_extents[d],
                center_point[d] + bounds_extents[d],
            );
        }
        Some(BoundingBox(bbox))
    }
};
}

impl_intersection_region!(intersection_region_3, 3);
impl_intersection_region!(intersection_region_4, 4);

/// Specialized version for 1D case.
fn intersection_region_1<const D: usize>(
    planes: &[BoundingHyperplanes<D>],
) -> Option<BoundingBox<D>> {
    debug_assert!(planes.len() == 1);

    Some(BoundingBox(std::array::from_fn(|d| {
        let inv_z = 1.0 / planes[d].normal[0];
        let x_min = planes[d].bmin * inv_z;
        let x_max = planes[d].bmax * inv_z;
        (x_min.min(x_max), x_min.max(x_max))
    })))
}

fn intersection_region_2<const D: usize>(
    planes: &[BoundingHyperplanes<D>],
) -> Option<BoundingBox<D>> {
    debug_assert!(
        planes.len() == 2,
        "Expected 2 bounding hyperplanes for 2D case, got {}",
        planes.len()
    );

    let mut lhs = Mat2::zeros();
    lhs.data[0] = planes[0].normal[0];
    lhs.data[1] = planes[0].normal[1];
    lhs.data[2] = planes[1].normal[0];
    lhs.data[3] = planes[1].normal[1];

    let mut rhs = [
        0.5 * (planes[0].bmin + planes[0].bmax),
        0.5 * (planes[1].bmin + planes[1].bmax),
    ];
    let Some(center_point) = lhs.solve(rhs) else {
        return None;
    };

    // Compute the bounding box around the center point
    let mut bounds_extents = [0.0; 2];
    for d in 0..2 {
        let half_plane_disp = 0.5 * (planes[d].bmax - planes[d].bmin);
        rhs.fill(0.0);
        rhs[d] = half_plane_disp;
        let Some(solution) = lhs.solve(rhs) else {
            panic!(
                "Bounding hyperplanes are not linearly independent. Coincident curves?\n{:#?}\n{:?}",
                lhs, rhs
            );
        };
        for dim in 0..2 {
            bounds_extents[dim] += solution[dim].abs();
        }
    }

    let mut bbox = [(0.0, 0.0); D];
    for d in 0..2 {
        bbox[d] = (
            center_point[d] - bounds_extents[d],
            center_point[d] + bounds_extents[d],
        );
    }
    Some(BoundingBox(bbox))
}

#[derive(Debug, Clone)]
struct Mat2 {
    data: [f64; 4],
}

impl Mat2 {
    pub fn zeros() -> Self {
        Self { data: [0.0; 4] }
    }

    pub fn solve(&self, rhs: [f64; 2]) -> Option<[f64; 2]> {
        let det = self.data[0] * self.data[3] - self.data[1] * self.data[2];
        if det.abs() < 1e-12 {
            return None;
        }
        let inv_det = 1.0 / det;
        let x0 = inv_det * (self.data[3] * rhs[0] - self.data[1] * rhs[1]);
        let x1 = inv_det * (-self.data[2] * rhs[0] + self.data[0] * rhs[1]);
        Some([x0, x1])
    }
}
