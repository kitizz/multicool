//
// Projected Polyhedra refinement
//

use smallvec::SmallVec;

use crate::{BezierSurface, bounding_box::BoundingBox};

pub struct ProjectedPolyhedraRefiner {
    hulls: Vec<HullIntersector>,
}

impl ProjectedPolyhedraRefiner {
    pub fn new(num_surfaces: usize) -> Self {
        Self {
            hulls: vec![HullIntersector::new(); num_surfaces],
        }
    }

    pub fn refine_regions<const D: usize>(
        &mut self,
        surfaces: &[BezierSurface<D>],
        valid_regions: &mut BoundingBox<D>,
    ) -> Result<(), crate::MulticoolError> {
        debug_assert!(self.hulls.len() == D);

        for (s_i, surface) in surfaces.iter().enumerate() {
            for d in 0..D {
                self.hulls[d].reset();
            }
            for control_point in surface.controls() {
                for d in 0..D {
                    self.hulls[d].add_point(control_point.axis_value(d), control_point.coeff());
                }
            }
            log::info!("Surface {} hull:", s_i);
            for d in 0..D {
                let (min_x, max_x) = self.hulls[d].intersect_zero_n2();

                log::info!("  {}: [{}, {}]", d, min_x, max_x);
                valid_regions[d].0 = valid_regions[d].0.max(min_x);
                valid_regions[d].1 = valid_regions[d].1.min(max_x);
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct HullIntersector {
    above: SmallVec<[[f64; 2]; 16]>,
    below: SmallVec<[[f64; 2]; 16]>,
    on_min: f64,
    on_max: f64,
}

impl HullIntersector {
    pub fn new() -> Self {
        Self {
            above: SmallVec::new(),
            below: SmallVec::new(),
            on_min: f64::INFINITY,
            on_max: -f64::INFINITY,
        }
    }

    pub fn reset(&mut self) {
        self.above.clear();
        self.below.clear();
        self.on_min = f64::INFINITY;
        self.on_max = -f64::INFINITY;
    }

    pub fn add_point(&mut self, x: f64, y: f64) {
        debug_assert!(x.is_finite() && y.is_finite(), "Invalid point ({x}, {y})");
        if y > 0.0 {
            self.above.push([x, y]);
        } else if y < 0.0 {
            self.below.push([x, y]);
        } else {
            self.on_min = self.on_min.min(x);
            self.on_max = self.on_max.max(x);
        }
    }

    /// Find the intersection region of the convex hull with the "x"-axis
    pub fn intersect_zero_n2(&mut self) -> (f64, f64) {
        let mut min_x = self.on_min;
        let mut max_x = self.on_max;

        for &[x1, z1] in &self.above {
            for &[x2, z2] in &self.below {
                let t = z1 / (z1 - z2);
                let x = x1 + t * (x2 - x1);
                min_x = min_x.min(x);
                max_x = max_x.max(x);
            }
        }
        (min_x, max_x)
    }
}
