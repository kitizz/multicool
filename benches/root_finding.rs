use criterion::{Criterion, criterion_group, criterion_main};
use multicool::{Monomial, MultivarPoly, MultivarPolySystem};
use std::hint::black_box;

fn root_find_circles(c: &mut Criterion) {
    // x^2 + (y + 0.25)^2 = 1
    let circle1 = MultivarPoly::new().add_monomials([
        Monomial::new(1.0, [2, 0]),
        Monomial::new(1.0, [0, 2]),
        Monomial::new(0.5, [0, 1]),
        Monomial::new(-0.875, [0, 0]),
    ]);

    // (x - 1)^2 + y^2 = 1
    let circle2 = MultivarPoly::new().add_monomials([
        Monomial::new(1.0, [2, 0]),
        Monomial::new(-2.0, [1, 0]),
        Monomial::new(1.0, [0, 2]),
    ]);

    c.bench_function("root_find_circles", |b| {
        b.iter(|| {
            let circle1 = black_box(circle1.clone());
            let circle2 = black_box(circle2.clone());
            let system =
                MultivarPolySystem::from_polys([circle1, circle2], [(-1.0, 2.0), (-1.0, 1.0)])
                    .unwrap();

            const EPS: f64 = 1e-9;
            let roots = system.roots_pp(EPS).unwrap();
            black_box(roots)
        })
    });
}

fn root_find_beziers(c: &mut Criterion) {
    c.bench_function("root_find_beziers", |b| {
        b.iter(|| {
            let system = MultivarPolySystem::from_beziers(
                [[0.0, 0.0], [1.1, 2.0], [2.0, -1.0], [3.0, 1.0]],
                [[0.0, 1.0], [1.0, -1.0], [2.0, 2.0], [3.0, 0.0]],
                [(0.0, 1.0), (0.0, 1.0)],
            );

            const EPS: f64 = 1e-9;
            let roots = system.roots_pp(EPS).unwrap();
            black_box(roots)
        })
    });
}

fn root_find_cylinders(c: &mut Criterion) {
    const DIMS: usize = 3;

    let mut polys = Vec::new();
    for d in 0..DIMS {
        polys.push(cylinder_poly(d));
    }

    c.bench_function("root_find_cylinders", |b| {
        b.iter(|| {
            let system = MultivarPolySystem::from_polys(
                polys.clone(),
                [(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)],
            )
            .unwrap();

            const EPS: f64 = 1e-9;
            let roots = system.roots_pp(EPS).unwrap();
            black_box(roots)
        })
    });
}

fn root_lp_find_circles(c: &mut Criterion) {
    // x^2 + (y + 0.25)^2 = 1
    let circle1 = MultivarPoly::new().add_monomials([
        Monomial::new(1.0, [2, 0]),
        Monomial::new(1.0, [0, 2]),
        Monomial::new(0.5, [0, 1]),
        Monomial::new(-0.875, [0, 0]),
    ]);

    // (x - 1)^2 + y^2 = 1
    let circle2 = MultivarPoly::new().add_monomials([
        Monomial::new(1.0, [2, 0]),
        Monomial::new(-2.0, [1, 0]),
        Monomial::new(1.0, [0, 2]),
    ]);

    c.bench_function("root_lp_find_circles", |b| {
        b.iter(|| {
            let circle1 = black_box(circle1.clone());
            let circle2 = black_box(circle2.clone());
            let system =
                MultivarPolySystem::from_polys([circle1, circle2], [(-1.0, 2.0), (-1.0, 1.0)])
                    .unwrap();

            const EPS: f64 = 1e-9;
            let roots = system.roots_hp(EPS).unwrap();
            black_box(roots)
        })
    });
}

fn root_lp_find_beziers(c: &mut Criterion) {
    c.bench_function("root_lp_find_beziers", |b| {
        b.iter(|| {
            let system = MultivarPolySystem::from_beziers(
                [[0.0, 0.0], [1.1, 2.0], [2.0, -1.0], [3.0, 1.0]],
                [[0.0, 1.0], [1.0, -1.0], [2.0, 2.0], [3.0, 0.0]],
                [(0.0, 1.0), (0.0, 1.0)],
            );

            const EPS: f64 = 1e-9;
            let roots = system.roots_hp(EPS).unwrap();
            black_box(roots)
        })
    });
}

fn root_lp_find_cylinders(c: &mut Criterion) {
    const DIMS: usize = 3;

    let mut polys = Vec::new();
    for d in 0..DIMS {
        polys.push(cylinder_poly(d));
    }

    c.bench_function("root_lp_find_cylinders", |b| {
        b.iter(|| {
            let system = MultivarPolySystem::from_polys(
                polys.clone(),
                [(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)],
            )
            .unwrap();

            const EPS: f64 = 1e-9;
            let roots = system.roots_hp(EPS).unwrap();
            black_box(roots)
        })
    });
}

/// Create a cylinder polynomial aligned along the specified axis.
///
/// The polynomial equation exists in D dimensions and represents a unit cylinder
/// centered at the origin, extending infinitely along the specified axis.
fn cylinder_poly<const D: usize>(axis: usize) -> crate::MultivarPoly<D> {
    use crate::{Monomial, MultivarPoly};
    // use smallvec::{SmallVec, smallvec};

    let mut terms = vec![];
    for d in 0..D {
        if d == axis {
            continue;
        }
        let mut exp = [0u8; D];
        exp[d] = 2;
        terms.push(Monomial::new(1.0, exp));
    }
    terms.push(Monomial::new(-1.0, [0u8; D]));

    MultivarPoly::from_monomials(terms)
}

criterion_group!(
    benches,
    root_find_circles,
    root_find_beziers,
    root_find_cylinders,
    root_lp_find_circles,
    root_lp_find_beziers,
    root_lp_find_cylinders
);
criterion_main!(benches);
