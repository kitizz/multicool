use criterion::{Criterion, criterion_group, criterion_main};
use multicool::{Monomial, MultivarPoly};
use std::hint::black_box;

fn x2_sub_bezier(c: &mut Criterion) {
    let x2 = Monomial::new(1.0, [2, 0]);
    let bezier = [
        [-3.0, 1.0, 2.5, -10.0], // x
        [4.0, 2.0, -5.0, -3.0],  // y
    ];
    c.bench_function("x2_sub_bezier", |b| {
        b.iter(|| x2.sub_polys(black_box(&bezier)).unwrap())
    });
}

fn circle_sub_bezier(c: &mut Criterion) {
    let circle = MultivarPoly::new().add_monomials([
        Monomial::new(-1.0, [0, 0]),
        Monomial::new(1.0, [2, 0]),
        Monomial::new(1.0, [0, 2]),
    ]);
    let bezier = [
        [-3.0, 1.0, 2.5, -10.0], // x
        [4.0, 2.0, -5.0, -3.0],  // y
    ];
    c.bench_function("circle_sub_bezier", |b| {
        b.iter(|| circle.sub_polys(black_box(&bezier)).unwrap())
    });
}

fn scale_to_unit_box(c: &mut Criterion) {
    let circle = MultivarPoly::new().add_monomials([
        Monomial::new(-1.0, [0, 0]),
        Monomial::new(1.0, [2, 0]),
        Monomial::new(1.0, [0, 2]),
    ]);
    let bezier = [
        [-3.0, 1.0, 2.5, -10.0], // x
        [4.0, 2.0, -5.0, -3.0],  // y
    ];
    let mixed = circle.sub_polys(&bezier).unwrap(); // 14 terms
    let source_box = [[-3.0, 4.0], [5.0, 5.1]];

    // dbg!(&mixed);

    c.bench_function("scale_to_unit_box", |b| {
        b.iter(|| {
            // Scale the mixed polynomial to the unit box
            let affine_coeffs: [[f64; 2]; 2] = core::array::from_fn(|i| {
                let [min, max] = source_box[i];
                let inv_delta = 1.0 / (max - min);
                [-min * inv_delta, inv_delta]
            });
            black_box(mixed.sub_polys(black_box(&affine_coeffs)).unwrap())
        })
    });
}

criterion_group!(benches, x2_sub_bezier, circle_sub_bezier, scale_to_unit_box);
criterion_main!(benches);
