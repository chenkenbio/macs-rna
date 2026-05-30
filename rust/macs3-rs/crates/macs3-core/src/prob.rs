//! Statistical helpers ported from `MACS3/Signal/Prob.py` (MACS3 3.0.4).
//!
//! Types follow the Cython source exactly: observation counts are `u32`
//! (`uint32_t`), trial/coordinate counts are `i64` (`int64_t`), and all
//! probability arithmetic is `f64` (`float64_t`). The only place values are
//! truncated to lower precision is the explicit 5-digit `round()` MACS applies
//! to the log10 large-lambda tails — replicated here with Python's
//! round-half-to-even semantics.

/// `LSTEP` from Prob.py (number of natural-log units per scaling step).
pub const LSTEP: i32 = 200;
/// `EXPTHRES = exp(LSTEP)`.
pub const EXPTHRES: f64 = 7.225973768125749e86; // exp(200)
/// `EXPSTEP = exp(-LSTEP)`.
pub const EXPSTEP: f64 = 1.3838965267367376e-87; // exp(-200)
/// `bigx` threshold used in the chi-square routines.
pub const BIGX: f64 = 20.0;

/// Natural log of 10, matching libc `M_LN10`.
const M_LN10: f64 = std::f64::consts::LN_10;

/// Python `round(x, ndigits)` — round-half-to-even at `ndigits` decimals.
///
/// CPython's `round` for floats uses a correctly-rounded decimal algorithm.
/// MACS calls it only with `ndigits = 5` on already-finite values, so we
/// reproduce it via a scale / round-half-even / unscale that matches CPython
/// for these magnitudes. We implement banker's rounding on the scaled value.
fn py_round5(x: f64) -> f64 {
    py_round(x, 5)
}

/// Generic Python-style `round(x, ndigits)` (round-half-to-even).
fn py_round(x: f64, ndigits: i32) -> f64 {
    if !x.is_finite() {
        return x;
    }
    // Use the same approach as CPython's float.__round__: format with the
    // requested number of digits (which performs correct round-half-even),
    // then parse back. This avoids the double-rounding error of scale/round.
    let s = format!("{:.*}", ndigits.max(0) as usize, x);
    s.parse::<f64>().unwrap_or(x)
}

/// `logspace_add(logx, logy)` = `log(exp(logx) + exp(logy))`.
///
/// Ported from `Prob.py::logspace_add`.
pub fn logspace_add(logx: f64, logy: f64) -> f64 {
    if logx > logy {
        logx + (logy - logx).exp().ln_1p()
    } else {
        logy + (logx - logy).exp().ln_1p()
    }
}

/// `factorial(n)` returning an f64, matching `Prob.py::factorial`.
///
/// Computes the product `2 * 3 * ... * n` in f64 (overflows to `inf` beyond
/// ~170, exactly like the Python loop).
pub fn factorial(n: u32) -> f64 {
    // n is unsigned, so the `n < 0` branch in the source is dead.
    let mut fact: f64 = 1.0;
    let mut i: u64 = 2;
    while i <= n as u64 {
        fact *= i as f64;
        i += 1;
    }
    fact
}

/// `poisson_cdf(n, lam, lower, log10)` dispatcher (`Prob.py::poisson_cdf`).
///
/// # Panics
/// Panics if `lam <= 0.0` (the source asserts `lam > 0`).
pub fn poisson_cdf(n: u32, lam: f64, lower: bool, log10: bool) -> f64 {
    assert!(lam > 0.0, "Lambda must > 0, however we got {lam}");

    if log10 {
        if lower {
            return log10_poisson_cdf_p_large_lambda(n, lam);
        } else {
            return log10_poisson_cdf_q_large_lambda(n, lam);
        }
    }

    if lower {
        if lam > 700.0 {
            poisson_cdf_large_lambda(n, lam)
        } else {
            poisson_cdf_small(n, lam)
        }
    } else if lam > 700.0 {
        poisson_cdf_q_large_lambda(n, lam)
    } else {
        poisson_cdf_q(n, lam)
    }
}

/// `__poisson_cdf` — lower tail, small lambda.
fn poisson_cdf_small(k: u32, a: f64) -> f64 {
    let mut nextcdf = (-a).exp();
    let mut cdf = nextcdf;
    for i in 1..=k {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / i as f64;
        cdf += nextcdf;
    }
    if cdf > 1.0 {
        1.0
    } else {
        cdf
    }
}

/// `__poisson_cdf_large_lambda` — lower tail, lambda > 700.
fn poisson_cdf_large_lambda(k: u32, a: f64) -> f64 {
    debug_assert!(a > 700.0);
    let mut num_parts: i32 = (a / LSTEP as f64) as i32;
    let mut lastexp = (-(a % LSTEP as f64)).exp();
    let mut nextcdf = EXPSTEP;
    let mut cdf: f64 = 0.0;

    num_parts -= 1;

    for i in 1..=k {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / i as f64;
        cdf += nextcdf;
        if nextcdf > EXPTHRES || cdf > EXPTHRES {
            if num_parts >= 1 {
                cdf *= EXPSTEP;
                nextcdf *= EXPSTEP;
                num_parts -= 1;
            } else {
                cdf *= lastexp;
                lastexp = 1.0;
            }
        }
    }

    for _ in 0..num_parts {
        cdf *= EXPSTEP;
    }
    cdf *= lastexp;
    cdf
}

/// `__poisson_cdf_Q` — upper tail, small lambda.
fn poisson_cdf_q(k: u32, a: f64) -> f64 {
    let mut nextcdf = (-a).exp();
    let mut cdf: f64 = 0.0;

    for i in 1..=k {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / i as f64;
    }
    let mut i = k + 1;
    while nextcdf > 0.0 {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / i as f64;
        cdf += nextcdf;
        i += 1;
    }
    cdf
}

/// `__poisson_cdf_Q_large_lambda` — upper tail, lambda > 700.
///
/// # Panics
/// Panics ("Unexpected error") in the same situation the Python source raises.
fn poisson_cdf_q_large_lambda(k: u32, a: f64) -> f64 {
    debug_assert!(a > 700.0);
    let mut num_parts: i64 = (a / LSTEP as f64) as i64;
    let mut lastexp = (-(a % LSTEP as f64)).exp();
    let mut nextcdf = EXPSTEP;
    let mut cdf: f64 = 0.0;

    num_parts -= 1;

    for i in 1..=k {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / i as f64;
        if nextcdf > EXPTHRES {
            if num_parts >= 1 {
                nextcdf *= EXPSTEP;
                num_parts -= 1;
            } else {
                panic!("Unexpected error");
            }
        }
    }

    let mut i = k + 1;
    while nextcdf > 0.0 {
        let lastcdf = nextcdf;
        nextcdf = lastcdf * a / i as f64;
        cdf += nextcdf;
        i += 1;
        if nextcdf > EXPTHRES || cdf > EXPTHRES {
            if num_parts >= 1 {
                cdf *= EXPSTEP;
                nextcdf *= EXPSTEP;
                num_parts -= 1;
            } else {
                cdf *= lastexp;
                lastexp = 1.0;
            }
        }
    }

    for _ in 0..num_parts {
        cdf *= EXPSTEP;
    }
    cdf *= lastexp;
    cdf
}

/// `log10_poisson_cdf_P_large_lambda` — lower tail in log10 space.
fn log10_poisson_cdf_p_large_lambda(k: u32, lbd: f64) -> f64 {
    let ln_lbd = lbd.ln();
    let mut m: i64 = k as i64;
    let mut sum_ln_m: f64 = 0.0;
    let mut i: i64 = 1;
    while i <= m {
        sum_ln_m += (i as f64).ln();
        i += 1;
    }
    let mut logx = m as f64 * ln_lbd - sum_ln_m;
    let mut residue = logx;

    while m > 1 {
        m -= 1;
        let logy = logx - ln_lbd + (m as f64).ln();
        let pre_residue = residue;
        residue = logspace_add(pre_residue, logy);
        if (pre_residue - residue).abs() < 1e-10 {
            break;
        }
        logx = logy;
    }

    py_round5((residue - lbd) / M_LN10)
}

/// `log10_poisson_cdf_Q_large_lambda` — upper tail in log10 space.
fn log10_poisson_cdf_q_large_lambda(k: u32, lbd: f64) -> f64 {
    let ln_lbd = lbd.ln();
    let mut m: i64 = k as i64 + 1;
    let mut sum_ln_m: f64 = 0.0;
    let mut i: i64 = 1;
    while i <= m {
        sum_ln_m += (i as f64).ln();
        i += 1;
    }
    let mut logx = m as f64 * ln_lbd - sum_ln_m;
    let mut residue = logx;

    loop {
        m += 1;
        let logy = logx + ln_lbd - (m as f64).ln();
        let pre_residue = residue;
        residue = logspace_add(pre_residue, logy);
        if (pre_residue - residue).abs() < 1e-5 {
            break;
        }
        logx = logy;
    }

    py_round5((residue - lbd) / M_LN10)
}

/// `poisson_cdf_inv(cdf, lam, maximum=1000)` (`Prob.py::poisson_cdf_inv`).
///
/// # Panics
/// Panics if `lam >= 740` or `cdf` is outside `[0, 1]`.
pub fn poisson_cdf_inv(cdf: f64, lam: f64, maximum: i32) -> i32 {
    assert!(lam < 740.0);
    if cdf < 0.0 || cdf > 1.0 {
        panic!("CDF must >= 0 and <= 1");
    } else if cdf == 0.0 {
        return 0;
    }

    let mut newval = (-lam).exp();
    let mut sum2 = newval;
    for i in 1..=maximum {
        let sumold = sum2;
        let lastval = newval;
        newval = lastval * lam / i as f64;
        sum2 += newval;
        if sumold <= cdf && cdf <= sum2 {
            return i;
        }
    }
    maximum
}

/// `poisson_cdf_Q_inv(cdf, lam, maximum=1000)` (identical body to
/// `poisson_cdf_inv` in the source).
///
/// # Panics
/// Panics if `lam >= 740` or `cdf` is outside `[0, 1]`.
pub fn poisson_cdf_q_inv(cdf: f64, lam: f64, maximum: i32) -> i32 {
    assert!(lam < 740.0);
    if cdf < 0.0 || cdf > 1.0 {
        panic!("CDF must >= 0 and <= 1");
    } else if cdf == 0.0 {
        return 0;
    }

    let mut newval = (-lam).exp();
    let mut sum2 = newval;
    for i in 1..=maximum {
        let sumold = sum2;
        let lastval = newval;
        newval = lastval * lam / i as f64;
        sum2 += newval;
        if sumold <= cdf && cdf <= sum2 {
            return i;
        }
    }
    maximum
}

/// `chisq_pvalue_e(x, df)` — even-df upper-tail chi-square CDF.
///
/// Ported from `Prob.py::chisq_pvalue_e`. Assumes `df` is an even number > 1.
pub fn chisq_pvalue_e(x: f64, df: u32) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    let a = 0.5 * x;
    let y = ex20(-a);
    let mut s = y;
    if df > 2 {
        let xx = 0.5 * (df as f64 - 1.0);
        let mut z = 1.0;
        if a > BIGX {
            let mut e = 0.0;
            let c = a.ln();
            while z <= xx {
                e = z.ln() + e;
                s += ex20(c * z - a - e);
                z += 1.0;
            }
            s
        } else {
            let mut e = 1.0;
            let mut c = 0.0;
            while z <= xx {
                e *= a / z;
                c += e;
                z += 1.0;
            }
            c * y + s
        }
    } else {
        s
    }
}

/// `chisq_logp_e(x, df, log10)` — even-df upper-tail chi-square in log space.
///
/// Ported from `Prob.py::chisq_logp_e`. Returns `-logp` (or `-log10p`).
pub fn chisq_logp_e(x: f64, df: u32, log10: bool) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let a = 0.5 * x;
    let y = (-a).exp();
    let mut s = -a;
    if df > 2 {
        let xx = 0.5 * (df as f64 - 1.0);
        let mut z = 1.0;
        if a > BIGX {
            let mut e = 0.0;
            let c = a.ln();
            while z <= xx {
                e += z.ln();
                s = logspace_add(s, c * z - a - e);
                z += 1.0;
            }
        } else {
            let mut e = 1.0;
            let mut c = 0.0;
            while z <= xx {
                e *= a / z;
                c += e;
                z += 1.0;
            }
            s = (y + c * y).ln();
        }
    }
    if log10 {
        -s / (10.0_f64).ln()
    } else {
        -s
    }
}

/// `ex20(x)` — `exp(x)` but 0 below -20 (`Prob.py::ex20`).
fn ex20(x: f64) -> f64 {
    if x < -20.0 {
        0.0
    } else {
        x.exp()
    }
}

/// `binomial_pdf(x, a, b)` (`Prob.py::binomial_pdf`), by H. Gene Shin.
pub fn binomial_pdf(x: i64, a: i64, b: f64) -> f64 {
    if a < 1 {
        return 0.0;
    } else if x < 0 || a < x {
        return 0.0;
    } else if b == 0.0 {
        return if x == 0 { 1.0 } else { 0.0 };
    } else if b == 1.0 {
        return if x == a { 1.0 } else { 0.0 };
    }

    let (p, mn, mx): (f64, i64, i64);
    if x > a - x {
        p = 1.0 - b;
        mn = a - x;
        mx = x;
    } else {
        p = b;
        mn = x;
        mx = a - x;
    }
    let mut pdf: f64 = 1.0;
    let mut t: i64 = 0;
    for q in 1..=mn {
        pdf *= (a - q + 1) as f64 * p / (mn - q + 1) as f64;
        if pdf < 1e-100 {
            while pdf < 1e-3 {
                pdf /= 1.0 - p;
                t -= 1;
            }
        }
        if pdf > 1e100 {
            while pdf > 1e3 && t < mx {
                pdf *= 1.0 - p;
                t += 1;
            }
        }
    }
    for _ in 0..(mx - t) {
        pdf *= 1.0 - p;
    }
    pdf
}

/// `binomial_cdf(x, a, b, lower)` (`Prob.py::binomial_cdf`).
pub fn binomial_cdf(x: i64, a: i64, b: f64, lower: bool) -> f64 {
    if lower {
        binomial_cdf_f(x, a, b)
    } else {
        binomial_cdf_r(x, a, b)
    }
}

/// `binomial_sf(x, a, b, lower)` (`Prob.py::binomial_sf`) — survival function.
pub fn binomial_sf(x: i64, a: i64, b: f64, lower: bool) -> f64 {
    if lower {
        1.0 - binomial_cdf_f(x, a, b)
    } else {
        1.0 - binomial_cdf_r(x, a, b)
    }
}

/// `_binomial_cdf_r` — upper tail.
fn binomial_cdf_r(x: i64, a: i64, b: f64) -> f64 {
    let argmax: i64 = (a as f64 * b) as i64;

    if x < 0 {
        return 1.0;
    } else if a < x {
        return 0.0;
    } else if b == 0.0 {
        return 0.0;
    } else if b == 1.0 {
        return 1.0;
    }

    if x < argmax {
        let seedpdf = binomial_pdf(argmax, a, b);
        let mut pdf = seedpdf;
        let mut cdf = pdf;
        let mut i = argmax - 1;
        while i > x {
            pdf /= (a - i) as f64 * b / (1.0 - b) / (i + 1) as f64;
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
            i -= 1;
        }

        pdf = seedpdf;
        let mut i = argmax;
        loop {
            pdf *= (a - i) as f64 * b / (1.0 - b) / (i + 1) as f64;
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
            i += 1;
        }
        cdf.min(1.0)
    } else {
        let mut pdf = binomial_pdf(x + 1, a, b);
        let mut cdf = pdf;
        let mut i = x + 1;
        loop {
            pdf *= (a - i) as f64 * b / (1.0 - b) / (i + 1) as f64;
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
            i += 1;
        }
        cdf.min(1.0)
    }
}

/// `_binomial_cdf_f` — lower tail.
fn binomial_cdf_f(x: i64, a: i64, b: f64) -> f64 {
    let argmax: i64 = (a as f64 * b) as i64;

    if x < 0 {
        return 0.0;
    } else if a < x {
        return 1.0;
    } else if b == 0.0 {
        return 1.0;
    } else if b == 1.0 {
        return 0.0;
    }

    if x > argmax {
        let seedpdf = binomial_pdf(argmax, a, b);
        let mut pdf = seedpdf;
        let mut cdf = pdf;
        let mut i = argmax - 1;
        while i >= 0 {
            pdf /= (a - i) as f64 * b / (1.0 - b) / (i + 1) as f64;
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
            i -= 1;
        }

        pdf = seedpdf;
        let mut i = argmax;
        while i < x {
            pdf *= (a - i) as f64 * b / (1.0 - b) / (i + 1) as f64;
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
            i += 1;
        }
        cdf.min(1.0)
    } else {
        let mut pdf = binomial_pdf(x, a, b);
        let mut cdf = pdf;
        let mut i = x - 1;
        while i >= 0 {
            pdf /= (a - i) as f64 * b / (1.0 - b) / (i + 1) as f64;
            if pdf == 0.0 {
                break;
            }
            cdf += pdf;
            i -= 1;
        }
        cdf.min(1.0)
    }
}

/// `binomial_cdf_inv(cdf, a, b)` — invert the binomial CDF (lower tail only).
///
/// # Panics
/// Panics if `cdf` is outside `[0, 1]`.
pub fn binomial_cdf_inv(cdf: f64, a: i64, b: f64) -> i64 {
    if cdf < 0.0 || cdf > 1.0 {
        panic!("CDF must >= 0 or <= 1");
    }
    let mut cdf2: f64 = 0.0;
    for x in 0..=a {
        let pdf = binomial_pdf(x, a, b);
        cdf2 += pdf;
        if cdf < cdf2 {
            return x;
        }
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(got: f64, want: f64, ctx: &str) {
        if want == 0.0 {
            assert!(got.abs() < 1e-300, "{ctx}: got {got}, want {want}");
        } else {
            let rel = ((got - want) / want).abs();
            assert!(
                rel < 1e-12 || (got - want).abs() < 1e-300,
                "{ctx}: got {got}, want {want}, rel {rel}"
            );
        }
    }

    #[test]
    fn factorial_known() {
        assert_eq!(factorial(0), 1.0);
        assert_eq!(factorial(1), 1.0);
        assert_eq!(factorial(10), 3628800.0);
    }

    /// Ported assertions from `references/.../test/test_Prob.py`.
    #[test]
    fn poisson_cdf_test_prob_py() {
        // n1 = (80, 100)
        let a = (poisson_cdf(80, 100.0, false, false) * 1e5).round() / 1e5;
        let b = (poisson_cdf(80, 100.0, true, false) * 1e5).round() / 1e5;
        assert_eq!(a, (0.9773508_f64 * 1e5).round() / 1e5);
        assert_eq!(b, (0.02264918_f64 * 1e5).round() / 1e5);
    }

    #[test]
    fn binomial_cdf_inv_test_prob_py() {
        assert_eq!(binomial_cdf_inv(0.1, 1000, 0.01), 6);
        assert_eq!(binomial_cdf_inv(0.01, 1000, 0.01), 3);
    }

    /// Compare against the MACS3-generated golden fixture for all functions.
    /// See `tests/fixtures/gen_prob_fixture.py`.
    #[test]
    fn matches_macs3_golden() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/prob_golden.tsv"
        );
        let text = std::fs::read_to_string(path).expect("read prob_golden.tsv");
        let mut n = 0usize;
        for line in text.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let func = f[0];
            match func {
                "poisson_cdf" => {
                    let nn: u32 = f[1].parse().unwrap();
                    let lam: f64 = f[2].parse().unwrap();
                    let lower = f[3] == "1";
                    let log10 = f[4] == "1";
                    let want = f64::from_bits(u64::from_str_radix(f[5], 16).unwrap());
                    let got = poisson_cdf(nn, lam, lower, log10);
                    // log10 large-lambda path is rounded to 5 digits in MACS,
                    // so equality is exact; others within tight relative tol.
                    if log10 {
                        assert_eq!(got, want, "poisson_cdf {nn} {lam} log10 path");
                    } else {
                        assert_close(got, want, &format!("poisson_cdf {nn} {lam} {lower}"));
                    }
                }
                "factorial" => {
                    let nn: u32 = f[1].parse().unwrap();
                    let want = f64::from_bits(u64::from_str_radix(f[5], 16).unwrap());
                    assert_eq!(factorial(nn), want, "factorial {nn}");
                }
                "chisq_logp_e" => {
                    let x: f64 = f[1].parse().unwrap();
                    let df: u32 = f[2].parse().unwrap();
                    let log10 = f[4] == "1";
                    let want = f64::from_bits(u64::from_str_radix(f[5], 16).unwrap());
                    assert_close(chisq_logp_e(x, df, log10), want, &format!("chisq_logp_e {x} {df}"));
                }
                "binomial_cdf" => {
                    let x: i64 = f[1].parse().unwrap();
                    let a: i64 = f[2].parse().unwrap();
                    let b: f64 = f[3].parse().unwrap();
                    let lower = f[4] == "1";
                    let want = f64::from_bits(u64::from_str_radix(f[5], 16).unwrap());
                    assert_close(binomial_cdf(x, a, b, lower), want, &format!("binomial_cdf {x} {a} {b} {lower}"));
                }
                "binomial_pdf" => {
                    let x: i64 = f[1].parse().unwrap();
                    let a: i64 = f[2].parse().unwrap();
                    let b: f64 = f[3].parse().unwrap();
                    let want = f64::from_bits(u64::from_str_radix(f[5], 16).unwrap());
                    assert_close(binomial_pdf(x, a, b), want, &format!("binomial_pdf {x} {a} {b}"));
                }
                "binomial_cdf_inv" => {
                    let cdf: f64 = f[1].parse().unwrap();
                    let a: i64 = f[2].parse().unwrap();
                    let b: f64 = f[3].parse().unwrap();
                    let want: i64 = f[5].parse().unwrap();
                    assert_eq!(binomial_cdf_inv(cdf, a, b), want, "binomial_cdf_inv {cdf} {a} {b}");
                }
                other => panic!("unknown func {other}"),
            }
            n += 1;
        }
        assert!(n > 100, "expected many fixture rows, got {n}");
    }
}
