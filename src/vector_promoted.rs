use std::ops::{Add, Mul};

/// A vector that been promoted from size N to size N+1 by adding an extra element.
///
/// This is currently needed since (stable) Rust does not yet support const generics arithmetic.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VectorPromoted<const N: usize, T: Copy = f64> {
    base: [T; N],
    extra: T,
}

impl<const N: usize, T: Copy> VectorPromoted<N, T> {
    /// Create a new promoted vector from the base vector and the extra element.
    pub fn new(base: [T; N], extra: T) -> Self {
        Self { base, extra }
    }

    /// Get the base vector (first N elements).
    pub fn base(&self) -> &[T; N] {
        &self.base
    }

    /// Get the extra element (N+1-th element).
    pub fn extra(&self) -> T {
        self.extra
    }

    pub fn dot(&self, other: &Self) -> T
    where
        T: Add<Output = T> + Mul<Output = T> + Default,
    {
        let mut sum = T::default();
        for i in 0..N {
            sum = sum + self.base[i] * other.base[i];
        }
        sum + self.extra * other.extra
    }
}

impl<const N: usize> VectorPromoted<N, f64> {
    /// Normalize the vector in-place.
    pub fn normalized(mut self) -> Self {
        let inv_norm = 1.0 / (self.dot(&self)).sqrt();
        for i in 0..N {
            self.base[i] *= inv_norm;
        }
        self.extra *= inv_norm;
        self
    }
}

impl<const N: usize, T> std::ops::Sub<&VectorPromoted<N, T>> for VectorPromoted<N, T>
where
    T: std::ops::Sub<Output = T> + Copy,
{
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        let mut base = self.base;
        for i in 0..N {
            base[i] = base[i] - rhs.base[i];
        }
        let extra = self.extra - rhs.extra;
        Self { base, extra }
    }
}

impl<const N: usize, T: Copy> std::ops::Index<usize> for VectorPromoted<N, T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        if index < N {
            &self.base[index]
        } else if index == N {
            &self.extra
        } else {
            panic!("Index out of bounds for VectorPromoted");
        }
    }
}

impl<const N: usize, T: Copy> Default for VectorPromoted<N, T>
where
    T: Default,
{
    fn default() -> Self {
        Self {
            base: [T::default(); N],
            extra: T::default(),
        }
    }
}
