//! C `printf`-compatible float formatters for f32 values widened to f64.
//!
//! MACS3 writes pileup / score values with C `printf` family conversions:
//!   - `%.5f` (see `cPosValCalculation.c::write_pv_array_to_bedGraph` and
//!     `IO/BedGraphIO.py::write_bedGraph`),
//!   - `%.6g` (see `IO/PeakIO.py::write_to_narrowPeak` / `write_to_broadPeak`
//!     for the `fc`/`pscore`/`qscore` columns).
//!
//! Because all internal values are stored as `f32` (C `float`), the C runtime
//! widens them to `double` (f64) before formatting. We replicate that exactly:
//! widen `f32 -> f64`, then format with correctly-rounded fixed / general
//! conversions. Rust's standard `{:.5}` / `{:.*}` formatter performs the same
//! IEEE round-half-to-even decimal rounding that glibc / CPython use, so we
//! build `%g` on top of `{:.*e}` and `{:.*}` rather than pulling in a
//! shortest-round-trip formatter (`ryu` is explicitly forbidden — it would not
//! match C output for values like `1.0`).
//!
//! The `String`-returning functions are the frozen public contract. The
//! `write_*` helpers append into a `String` for hot-loop reuse.

/// Format `v` (widened to f64) like C `printf("%.5f", (double)v)`.
///
/// Fixed notation, exactly 5 fractional digits, round-half-to-even. `-0.0`
/// prints as `-0.00000` (matching C / CPython).
///
/// # Examples
/// ```
/// use macs3_core::fmt::format_f5;
/// assert_eq!(format_f5(1.0_f32), "1.00000");
/// assert_eq!(format_f5(0.5_f32), "0.50000");
/// ```
pub fn format_f5(v: f32) -> String {
    let mut s = String::with_capacity(16);
    write_f5(&mut s, v);
    s
}

/// Append the `%.5f` rendering of `v` to `out` (no allocation of a new String).
pub fn write_f5(out: &mut String, v: f32) {
    use std::fmt::Write;
    let d: f64 = v as f64; // exact widening, as C does for varargs `float`.
    // Rust's `{:.5}` matches C `%.5f` for finite values including the sign of
    // negative zero. Non-finite handling below matches glibc lowercase forms.
    if d.is_nan() {
        out.push_str("nan");
        return;
    }
    if d.is_infinite() {
        out.push_str(if d < 0.0 { "-inf" } else { "inf" });
        return;
    }
    let _ = write!(out, "{:.5}", d);
}

/// Format `v` (widened to f64) like C `printf("%.6g", (double)v)`.
///
/// `%g` with precision `P` (here 6 significant digits, with `P` treated as 1
/// when 0): let `X` be the decimal exponent that `%e` would use. If
/// `X >= -4 && X < P`, use `%f`-style with precision `P-1-X`; otherwise use
/// `%e`-style with precision `P-1`. Then strip trailing zeros in the fraction
/// and a trailing decimal point (the C standard `#` flag is *not* set here).
///
/// # Examples
/// ```
/// use macs3_core::fmt::format_g6;
/// assert_eq!(format_g6(1.0_f32), "1");
/// assert_eq!(format_g6(0.5_f32), "0.5");
/// ```
pub fn format_g6(v: f32) -> String {
    let mut s = String::with_capacity(16);
    write_g6(&mut s, v);
    s
}

/// Append the `%.6g` rendering of `v` to `out`.
pub fn write_g6(out: &mut String, v: f32) {
    write_g(out, v as f64, 6);
}

/// Core `%g`-style formatter (precision = number of significant digits).
fn write_g(out: &mut String, d: f64, precision: usize) {
    use std::fmt::Write;

    if d.is_nan() {
        out.push_str("nan");
        return;
    }
    if d.is_infinite() {
        out.push_str(if d < 0.0 { "-inf" } else { "inf" });
        return;
    }

    // C: precision 0 is treated as 1 for %g.
    let p = if precision == 0 { 1 } else { precision };

    // Determine the exponent X that %e would produce (with p-1 fractional
    // digits, which is what controls the choice of style after rounding).
    // We obtain it from Rust's scientific formatting, which rounds the same
    // way C does, so the post-rounding exponent (e.g. 9.9999996e-1 -> 1e0)
    // is consistent with C's choice.
    let mut e_buf = String::with_capacity(24);
    let _ = write!(e_buf, "{:.*e}", p - 1, d);
    // e_buf looks like "1.23457e2" or "-9.99999e-1" or "0e0".
    let exp_pos = e_buf
        .as_bytes()
        .iter()
        .position(|&c| c == b'e')
        .expect("scientific format always contains 'e'");
    let x: i32 = e_buf[exp_pos + 1..]
        .parse()
        .expect("exponent is a valid integer");

    if x >= -4 && x < p as i32 {
        // Use %f style with precision p - 1 - X.
        let frac = (p as i32 - 1 - x).max(0) as usize;
        let mut f_buf = String::with_capacity(32);
        let _ = write!(f_buf, "{:.*}", frac, d);
        strip_g_trailing(out, &f_buf);
    } else {
        // Use %e style with precision p - 1, then strip trailing zeros in the
        // mantissa, and re-emit the exponent in C's `e[+-]NN` form (>= 2 digits).
        emit_e_style(out, &e_buf, exp_pos, x);
    }
}

/// Strip trailing zeros (and a dangling decimal point) from a `%f`-style
/// string, appending the cleaned result to `out`. Used for the `%g` no-`#`
/// case.
fn strip_g_trailing(out: &mut String, f_buf: &str) {
    if let Some(dot) = f_buf.find('.') {
        let bytes = f_buf.as_bytes();
        let mut end = f_buf.len();
        while end > dot + 1 && bytes[end - 1] == b'0' {
            end -= 1;
        }
        if end == dot + 1 {
            // all fractional digits stripped -> drop the '.' too
            end = dot;
        }
        out.push_str(&f_buf[..end]);
    } else {
        out.push_str(f_buf);
    }
}

/// Re-render a Rust scientific string into C `%e`-style with trailing-zero
/// stripping in the mantissa and a `e[+-]NN` exponent (minimum two digits).
fn emit_e_style(out: &mut String, e_buf: &str, exp_pos: usize, x: i32) {
    use std::fmt::Write;
    // mantissa is everything before 'e'
    let mantissa = &e_buf[..exp_pos];
    // strip trailing zeros / dot in the mantissa
    let mut m = String::with_capacity(mantissa.len());
    strip_g_trailing(&mut m, mantissa);
    out.push_str(&m);
    // C exponent: 'e', sign, at least two digits.
    let sign = if x < 0 { '-' } else { '+' };
    let mag = x.unsigned_abs();
    let _ = write!(out, "e{}{:02}", sign, mag);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spot_checks_f5() {
        assert_eq!(format_f5(0.0_f32), "0.00000");
        assert_eq!(format_f5(-0.0_f32), "-0.00000");
        assert_eq!(format_f5(1.0_f32), "1.00000");
        assert_eq!(format_f5(0.5_f32), "0.50000");
    }

    #[test]
    fn spot_checks_g6() {
        assert_eq!(format_g6(0.0_f32), "0");
        assert_eq!(format_g6(1.0_f32), "1");
        assert_eq!(format_g6(0.5_f32), "0.5");
        assert_eq!(format_g6(-1.0_f32), "-1");
    }

    /// Compare against the Python-generated golden fixture. Each row carries
    /// the f32 bit pattern (so we read the exact same value), the C `%.5f`
    /// string and the C `%.6g` string. See
    /// `tests/fixtures/gen_fmt_fixture.py` for generation.
    #[test]
    fn matches_python_golden() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/fmt_golden.tsv"
        );
        let text = std::fs::read_to_string(path).expect("read fmt_golden.tsv");
        let mut n = 0usize;
        for line in text.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let mut it = line.split('\t');
            let bits_hex = it.next().expect("bits");
            let want5 = it.next().expect("%.5f");
            let want6 = it.next().expect("%.6g");
            let bits = u32::from_str_radix(bits_hex, 16).expect("hex u32");
            let v = f32::from_bits(bits);
            let got5 = format_f5(v);
            let got6 = format_g6(v);
            assert_eq!(
                got5, want5,
                "format_f5 mismatch for bits {bits_hex} (v={v:?})"
            );
            assert_eq!(
                got6, want6,
                "format_g6 mismatch for bits {bits_hex} (v={v:?})"
            );
            n += 1;
        }
        assert!(n > 100, "expected many fixture rows, got {n}");
    }
}
