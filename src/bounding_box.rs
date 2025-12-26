#[derive(Clone, Copy, Debug)]
pub struct BoundingBox<const N: usize>(pub [(f64, f64); N]);

impl<const N: usize> BoundingBox<N> {
    pub fn at(&self, dim: usize) -> (f64, f64) {
        self.0[dim]
    }

    pub fn center(&self) -> [f64; N] {
        std::array::from_fn(|i| {
            let (min, max) = self.0[i];
            0.5 * (min + max)
        })
    }

    pub fn is_valid(&self) -> bool {
        for i in 0..N {
            let (min, max) = self.0[i];
            if min > max {
                return false;
            }
        }
        true
    }

    pub fn side_length(&self, dim: usize) -> f64 {
        let (min, max) = self.0[dim];
        max - min
    }

    pub fn largest_side(&self) -> (usize, f64) {
        let mut largest_dim = 0;
        let mut largest_length = 0.0;
        for i in 0..N {
            let length = self.side_length(i);
            if length > largest_length {
                largest_length = length;
                largest_dim = i;
            }
        }
        (largest_dim, largest_length)
    }

    pub fn subdivide(&self, dim: usize) -> (Self, Self) {
        let (min, max) = self.0[dim];
        let mid = 0.5 * (min + max);

        let mut box1 = *self;
        let mut box2 = *self;

        box1.0[dim] = (min, mid);
        box2.0[dim] = (mid, max);

        (box1, box2)
    }

    pub fn overlaps(&self, other: &Self) -> bool {
        for i in 0..N {
            let (a_min, a_max) = self.0[i];
            let (b_min, b_max) = other.0[i];
            if a_max < b_min || b_max < a_min {
                return false;
            }
        }
        true
    }
}

impl<const N: usize> From<[(f64, f64); N]> for BoundingBox<N> {
    fn from(bounds: [(f64, f64); N]) -> Self {
        BoundingBox(bounds)
    }
}

impl<const N: usize> std::ops::Index<usize> for BoundingBox<N> {
    type Output = (f64, f64);

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<const N: usize> std::ops::IndexMut<usize> for BoundingBox<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<const N: usize> PartialEq for BoundingBox<N> {
    fn eq(&self, other: &Self) -> bool {
        self.0
            .iter()
            .zip(other.0.iter())
            .all(|(&(a_min, a_max), &(b_min, b_max))| a_min == b_min && a_max == b_max)
    }
}

#[cfg(test)]
impl<const N: usize> approx::AbsDiffEq for BoundingBox<N> {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0
            .iter()
            .zip(other.0.iter())
            .all(|(&(a_min, a_max), &(b_min, b_max))| {
                (a_min - b_min).abs() <= epsilon && (a_max - b_max).abs() <= epsilon
            })
    }
}
