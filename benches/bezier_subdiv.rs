use criterion::{Criterion, criterion_group, criterion_main};
use multicool::{BezierSurface, Monomial, MultivarPoly};
use polycool::Poly;
use std::hint::black_box;

fn bezier_subsection_quadratic(c: &mut Criterion) {
    let surface = BezierSurface::new(
        &[
            1.0, -1.0, 1.0, //
            -1.0, -3.0, -1.0, //
            1.0, -1.0, 1.0,
        ],
        [3, 3],
        [(-1.0, 1.0), (-1.0, 1.0)],
    );
    c.bench_function("bezier_subsection_quadratic", |b| {
        b.iter(|| {
            let surf = black_box(surface.clone());
            let sub = surf.subsection([(0.65, 0.9), (-0.8, 0.3)]);
            black_box(sub)
        })
    });
}

fn bezier_surface(c: &mut Criterion) {
    let circle = MultivarPoly::new().add_monomials([
        Monomial::new(-1.0, [0, 0]),
        Monomial::new(1.0, [2, 0]),
        Monomial::new(1.0, [0, 2]),
    ]);
    let bezier = [
        Poly::new([-3.0, 1.0, 2.5, -10.0]), // x
        Poly::new([4.0, 2.0, -5.0, -3.0]),  // y
    ];
    let combined = circle.sub_polys(black_box(&bezier)).unwrap();

    c.bench_function("bezier_surface", |b| {
        b.iter(|| {
            //
            let surf = combined
                .bezier_surface(black_box([(2.0, 5.0), (-2.0, 2.0)]))
                .unwrap();
            black_box(surf)
        })
    });
}

fn bezier_subsection_larger(c: &mut Criterion) {
    let circle = MultivarPoly::new().add_monomials([
        Monomial::new(-1.0, [0, 0]),
        Monomial::new(1.0, [2, 0]),
        Monomial::new(1.0, [0, 2]),
    ]);
    let bezier = [
        Poly::new([-3.0, 1.0, 2.5, -10.0]), // x
        Poly::new([4.0, 2.0, -5.0, -3.0]),  // y
    ];
    let combined = circle.sub_polys(black_box(&bezier)).unwrap();
    let surface = combined.bezier_surface([(2.0, 5.0), (-2.0, 2.0)]).unwrap();

    c.bench_function("bezier_subsection_larger", |b| {
        b.iter(|| {
            let surf = black_box(surface.clone());
            let sub = surf.subsection([(-1.0, 1.0), (-1.0, 1.0)]);
            black_box(sub)
        })
    });
}

criterion_group!(
    benches,
    bezier_subsection_quadratic,
    bezier_surface,
    bezier_subsection_larger
);
criterion_main!(benches);
