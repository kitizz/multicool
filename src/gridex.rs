use num_traits::{Num, NumAssignOps};

/// Iterates over all multi-indices [i_1, i_2, ..., i_NUM_VAR] such that
/// If boundary == INCLUSIVE:
///     0 <= i_k <= max_degs[k] for each variable k.
/// If boundary == EXCLUSIVE:
///     0 <= i_k < max_degs[k] for each variable k.
pub fn gridex<const NUM_VAR: usize, T: Num + NumAssignOps + Copy + PartialOrd>(
    grid_size: [T; NUM_VAR],
    boundary: Boundary,
) -> impl Iterator<Item = [T; NUM_VAR]> {
    let mut indices = [T::zero(); NUM_VAR];
    let mut done = false;

    let grid_size = match boundary {
        Boundary::INCLUSIVE => grid_size,
        Boundary::EXCLUSIVE => {
            let mut gs = grid_size;
            for i in 0..NUM_VAR {
                if gs[i] > T::zero() {
                    gs[i] -= T::one();
                }
            }
            gs
        }
    };

    std::iter::from_fn(move || {
        if done {
            return None;
        }

        let current = indices;

        // Increment indices
        for i in (0..NUM_VAR).rev() {
            if indices[i] < grid_size[i] {
                indices[i] += T::one();
                break;
            } else {
                indices[i] = T::zero();
                if i == 0 {
                    done = true;
                }
            }
        }

        Some(current)
    })
}

pub fn gridex_incl<const NUM_VAR: usize, T: Num + NumAssignOps + Copy + PartialOrd>(
    grid_size: [T; NUM_VAR],
) -> impl Iterator<Item = [T; NUM_VAR]> {
    gridex::<NUM_VAR, T>(grid_size, Boundary::INCLUSIVE)
}

pub fn gridex_excl<const NUM_VAR: usize, T: Num + NumAssignOps + Copy + PartialOrd>(
    grid_size: [T; NUM_VAR],
) -> impl Iterator<Item = [T; NUM_VAR]> {
    gridex::<NUM_VAR, T>(grid_size, Boundary::EXCLUSIVE)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    INCLUSIVE,
    EXCLUSIVE,
}

#[cfg(test)]
mod tests {
    use pretty_assertions as pa;

    use super::*;

    #[test]
    fn multi_index_iterator_general() {
        let max_degs: [u8; 2] = [2, 1];
        let indexes = gridex_incl(max_degs).collect::<Vec<_>>();

        pa::assert_eq!(
            indexes,
            vec![[0, 0], [0, 1], [1, 0], [1, 1], [2, 0], [2, 1],]
        );
    }

    #[test]
    fn multi_index_iterator_zero() {
        let max_degs: [u8; 2] = [0, 0];
        let indexes = gridex_incl(max_degs).collect::<Vec<_>>();

        pa::assert_eq!(indexes, vec![[0, 0]]);
    }
}
