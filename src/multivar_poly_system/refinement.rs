use crate::{BezierSurface, MulticoolError, bounding_box::BoundingBox};

mod hyperplanes;
mod projected_polyhedra;

pub use hyperplanes::HyperplanesRefiner;
pub use projected_polyhedra::ProjectedPolyhedraRefiner;

pub enum Refiner<const D: usize> {
    ProjectedPolyhedra(ProjectedPolyhedraRefiner),
    Hyperplanes(HyperplanesRefiner<D>),
}

impl<const D: usize> Refiner<D> {
    pub fn refine_regions(
        &mut self,
        surfaces: &[BezierSurface<D>],
        valid_regions: &mut BoundingBox<D>,
    ) -> Result<(), MulticoolError> {
        match self {
            Refiner::Hyperplanes(r) => r.refine_regions(surfaces, valid_regions),
            Refiner::ProjectedPolyhedra(r) => r.refine_regions(surfaces, valid_regions),
        }
    }
}
