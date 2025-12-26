pub fn unit_box<const N: usize>() -> [(f64, f64); N] {
    std::array::from_fn(|_| (0.0, 1.0))
}

pub fn linspace(start: f64, end: f64, num: usize) -> impl Iterator<Item = f64> {
    let step = if num > 1 {
        (end - start) / (num - 1) as f64
    } else {
        0.0
    };
    (0..num).map(move |i| start + i as f64 * step)
}

#[allow(dead_code)]
pub fn init_test_logger() {
    use std::io::Write as _;
    let _ = env_logger::builder()
        .is_test(true)
        .format(|buf, record| {
            // Ansi256 ref: https://hexdocs.pm/color_palette/ansi_color_codes.html
            let bg = anstyle::Ansi256Color(240);
            let level_style = buf
                .default_level_style(record.level())
                .bg_color(Some(bg.into()));
            let grey = anstyle::Ansi256Color(255).on(bg);

            let filepath = match record.file() {
                Some(f) => {
                    // Get just the file name, not the full path.
                    let path = std::path::Path::new(f);
                    match path.file_name() {
                        Some(name) => name.to_string_lossy(),
                        None => "unknown".into(),
                    }
                }
                None => "unknown".into(),
            };
            writeln!(
                buf,
                "{grey}[{grey:#}{level_style}{}{level_style:#}{grey} {}:{}]{grey:#}   {}",
                record.level(),
                filepath,
                record.line().unwrap_or(0),
                record.args()
            )
        })
        .try_init();
}

/// Create a cylinder polynomial aligned along the specified axis.
///
/// The polynomial equation exists in D dimensions and represents a unit cylinder
/// centered at the origin, extending infinitely along the specified axis.
pub fn cylinder_poly<const D: usize>(axis: usize) -> crate::MultivarPoly<D> {
    use crate::{Monomial, MultivarPoly};
    // use smallvec::{SmallVec, smallvec};

    let mut terms = vec![];
    for d in 0..D {
        if d == axis {
            continue;
        }
        let mut exp = [0u8; D];
        exp[d] = 2;
        terms.push(Monomial { coeff: 1.0, exp });
    }
    terms.push(Monomial {
        coeff: -1.0,
        exp: [0u8; D],
    });

    MultivarPoly::from_monomials(terms)
}
