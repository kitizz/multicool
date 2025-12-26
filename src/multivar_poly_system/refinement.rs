use crate::{BezierSurface, MulticoolError, bounding_box::BoundingBox};

mod hyperplanes_lp;
mod projected_polyhedra;

pub use hyperplanes_lp::HyperplanesLpRefiner;
pub use projected_polyhedra::ProjectedPolyhedraRefiner;

pub enum Refiner<const D: usize> {
    ProjectedPolyhedra(ProjectedPolyhedraRefiner),
    HyperplanesLp(HyperplanesLpRefiner<D>),
}

impl<const D: usize> Refiner<D> {
    pub fn refine_regions(
        &mut self,
        surfaces: &[BezierSurface<D>],
        valid_regions: &mut BoundingBox<D>,
    ) -> Result<(), MulticoolError> {
        match self {
            Refiner::HyperplanesLp(r) => r.refine_regions(surfaces, valid_regions),
            Refiner::ProjectedPolyhedra(r) => r.refine_regions(surfaces, valid_regions),
        }
    }
}
