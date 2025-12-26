pub fn binomial_product<const N: usize>(ns: [u8; N], ks: [u8; N]) -> u64 {
    let mut product = 1;
    for i in 0..ns.len() {
        product *= binomial_coefficient(ns[i], ks[i]);
    }
    product
}

pub fn binomial_coefficient(n: u8, mut k: u8) -> u64 {
    if k > n {
        return 0;
    }
    k = k.min(n - k);
    if k == 0 {
        return 1;
    }
    let mut res = 1u64;
    let n = n as u64;
    let k = k as u64;
    for i in 0..k {
        res = res * (n - i) / (i + 1);
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binomial_coefficient_0k_test() {
        for n in 0..20 {
            assert_eq!(binomial_coefficient(n, 0), 1);
        }
    }

    #[test]
    fn binomial_coefficient_symmetric_test() {
        for n in 0..20 {
            for k in 0..=n / 2 {
                assert_eq!(binomial_coefficient(n, k), binomial_coefficient(n, n - k));
            }
        }
    }

    #[test]
    fn binomial_coefficient_1k_test() {
        for n in 1..20 {
            assert_eq!(binomial_coefficient(n, 1), n as u64);
        }
    }

    #[test]
    fn binomial_coefficient_general_test() {
        assert_eq!(binomial_coefficient(5, 2), 10);
        assert_eq!(binomial_coefficient(10, 3), 120);
        assert_eq!(binomial_coefficient(6, 4), 15);
        assert_eq!(binomial_coefficient(7, 5), 21);
    }

    #[test]
    fn binomial_product_test() {
        let ns = [5, 6, 7];
        let ks = [2, 3, 4];
        let expected =
            binomial_coefficient(5, 2) * binomial_coefficient(6, 3) * binomial_coefficient(7, 4);
        assert_eq!(binomial_product(ns, ks), expected);
    }
}
