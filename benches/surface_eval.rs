use criterion::{Criterion, criterion_group, criterion_main};
use multicool::BezierSurface;
use rand::{Rng as _, SeedableRng as _};
use std::hint::black_box;

fn surface_eval_deg2(c: &mut Criterion) {
    let surface = BezierSurface::new(
        &[
            1.0, -1.0, 1.0, //
            -1.0, -3.0, -1.0, //
            1.0, -1.0, 1.0,
        ],
        [3, 3],
        [(-1.0, 1.0), (-1.0, 1.0)],
    );
    c.bench_function("surface_eval_deg2", |b| {
        // Generate 1000 random (u,v) pairs in domain.
        let mut rng = rand::rng();
        let samples: Vec<(f64, f64)> = (0..1000)
            .map(|_| (rng.random_range(-1.0..1.0), rng.random_range(-1.0..1.0)))
            .collect();

        b.iter(|| {
            for &(u, v) in &samples {
                let val = black_box(surface.eval([u, v]));
                black_box(val);
            }
        })
    });
}

fn surface_eval_deg2_grad(c: &mut Criterion) {
    let surface = BezierSurface::new(
        &[
            1.0, -1.0, 1.0, //
            -1.0, -3.0, -1.0, //
            1.0, -1.0, 1.0,
        ],
        [3, 3],
        [(-1.0, 1.0), (-1.0, 1.0)],
    );
    c.bench_function("surface_eval_deg2_grad", |b| {
        // Generate 1000 random (u,v) pairs in domain.
        let mut rng = rand::rng();
        let samples: Vec<(f64, f64)> = (0..1000)
            .map(|_| (rng.random_range(-1.0..1.0), rng.random_range(-1.0..1.0)))
            .collect();

        b.iter(|| {
            for &(u, v) in &samples {
                let val = black_box(surface.eval_grad([u, v]));
                black_box(val);
            }
        })
    });
}

fn make_surface<const NDIM: usize>(degree: usize) -> BezierSurface<NDIM> {
    let seed = [1u8; 32];
    let mut rng = rand::rngs::SmallRng::from_seed(seed);

    let coeffs = (0..degree.pow(NDIM as u32))
        .map(|_| rng.random_range(-5.0..5.0))
        .collect::<Vec<_>>();

    let grid_size = [degree; NDIM];
    let domain = [(-1.0, 1.0); NDIM];
    BezierSurface::new(&coeffs, grid_size, domain)
}

fn surface_eval_deg7(c: &mut Criterion) {
    let surface = make_surface::<2>(7);

    c.bench_function("surface_eval_deg7", |b| {
        // Generate 1000 random (u,v) pairs in domain.
        let mut rng = rand::rng();
        let samples: Vec<(f64, f64)> = (0..1000)
            .map(|_| (rng.random_range(-1.0..1.0), rng.random_range(-1.0..1.0)))
            .collect();

        b.iter(|| {
            for &(u, v) in &samples {
                let val = black_box(surface.eval([u, v]));
                black_box(val);
            }
        })
    });
}

fn surface_eval_deg7_grad(c: &mut Criterion) {
    let surface = make_surface::<2>(7);

    c.bench_function("surface_eval_deg7_grad", |b| {
        // Generate 1000 random (u,v) pairs in domain.
        let mut rng = rand::rng();
        let samples: Vec<(f64, f64)> = (0..1000)
            .map(|_| (rng.random_range(-1.0..1.0), rng.random_range(-1.0..1.0)))
            .collect();

        b.iter(|| {
            for &(u, v) in &samples {
                let val = black_box(surface.eval_grad([u, v]));
                black_box(val);
            }
        })
    });
}

fn surface_eval_deg7_dim3(c: &mut Criterion) {
    let surface = make_surface::<3>(7);

    c.bench_function("surface_eval_deg7_dim3", |b| {
        // Generate 1000 random (u,v) pairs in domain.
        let mut rng = rand::rng();
        let samples: Vec<(f64, f64, f64)> = (0..1000)
            .map(|_| {
                (
                    rng.random_range(-1.0..1.0),
                    rng.random_range(-1.0..1.0),
                    rng.random_range(-1.0..1.0),
                )
            })
            .collect();

        b.iter(|| {
            for &(u, v, w) in &samples {
                let val = black_box(surface.eval([u, v, w]));
                black_box(val);
            }
        })
    });
}

fn surface_eval_deg7_dim3_grad(c: &mut Criterion) {
    let surface = make_surface::<3>(7);

    c.bench_function("surface_eval_deg7_dim3_grad", |b| {
        // Generate 1000 random (u,v) pairs in domain.
        let mut rng = rand::rng();
        let samples: Vec<(f64, f64, f64)> = (0..1000)
            .map(|_| {
                (
                    rng.random_range(-1.0..1.0),
                    rng.random_range(-1.0..1.0),
                    rng.random_range(-1.0..1.0),
                )
            })
            .collect();

        b.iter(|| {
            for &(u, v, w) in &samples {
                let val = black_box(surface.eval_grad([u, v, w]));
                black_box(val);
            }
        })
    });
}

criterion_group!(
    benches,
    surface_eval_deg2,
    surface_eval_deg2_grad,
    surface_eval_deg7,
    surface_eval_deg7_grad,
    surface_eval_deg7_dim3,
    surface_eval_deg7_dim3_grad,
);
criterion_main!(benches);
