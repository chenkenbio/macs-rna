//! Pileup primitives, ported from `MACS3/Signal/cPosValCalculation.c` and
//! `MACS3/Signal/Pileup.py`.
//!
//! Phase 1A owns the algorithm bodies here; the foundation phase fixes the
//! `PosVal` layout and the public function signatures that the pileup, callpeak
//! and bdgcmp engines depend on.

/// Position/value pair, mirroring `struct PosVal { int pos; float value; }`
/// in `cPosValCalculation.c`. A run `[prev_pos, pos)` carries `val`; the first
/// run starts at coordinate 0.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PosVal {
    /// Right-open end coordinate of the interval (0-based).
    pub pos: i32,
    /// Pileup / score value for the interval (C `float`).
    pub val: f32,
}

impl PosVal {
    /// Construct a new `PosVal`.
    pub fn new(pos: i32, val: f32) -> Self {
        PosVal { pos, val }
    }
}

/// Pile up single-end tags, mirroring
/// `cPosValCalculation.c::single_end_pileup`.
///
/// `plus_tags` / `minus_tags` are the 5'-end positions on each strand.
/// `five_shift` / `three_shift` extend each tag to a fragment; coordinates are
/// clamped to `[leftmost_coord, rightmost_coord]`. `scale_factor` multiplies
/// the pileup and `baseline_value` is the per-base floor.
///
/// Returns the coalesced `PosVal` array in bedGraph order.
pub fn single_end_pileup(
    plus_tags: &[i32],
    minus_tags: &[i32],
    five_shift: i32,
    three_shift: i32,
    leftmost_coord: i32,
    rightmost_coord: i32,
    scale_factor: f32,
    baseline_value: f32,
) -> Vec<PosVal> {
    let l = plus_tags.len() + minus_tags.len();
    let mut start_poss: Vec<i32> = Vec::with_capacity(l);
    let mut end_poss: Vec<i32> = Vec::with_capacity(l);

    // Plus strand: fragment spans [tag - five_shift, tag + three_shift).
    for &t in plus_tags {
        start_poss.push(t - five_shift);
        end_poss.push(t + three_shift);
    }
    // Minus strand: the 5' end is the rightmost coordinate, so the shifts swap.
    for &t in minus_tags {
        start_poss.push(t - three_shift);
        end_poss.push(t + five_shift);
    }

    start_poss.sort_unstable();
    end_poss.sort_unstable();

    fix_coordinates(&mut start_poss, leftmost_coord, rightmost_coord);
    fix_coordinates(&mut end_poss, leftmost_coord, rightmost_coord);

    quick_pileup(&start_poss, &end_poss, scale_factor, baseline_value)
}

/// Clamp a sorted coordinate array into `[leftmost, rightmost]`, mirroring
/// `cPosValCalculation.c::fix_coordinates`. Because the input is sorted, the
/// out-of-range values cluster at the two ends, so the scans stop early.
fn fix_coordinates(poss: &mut [i32], leftmost: i32, rightmost: i32) {
    for v in poss.iter_mut() {
        if *v < leftmost {
            *v = leftmost;
        } else {
            break;
        }
    }
    for v in poss.iter_mut().rev() {
        if *v > rightmost {
            *v = rightmost;
        } else {
            break;
        }
    }
}

/// Pile up from pre-sorted, pre-clamped start/end position arrays, mirroring
/// `cPosValCalculation.c::quick_pileup`.
///
/// Both arrays must have equal length `l`. A continuous run `[pre_p, p)`
/// carries `max(pileup * scale_factor, baseline_value)`; the leading run
/// `[0, pre_p)` (only emitted when `pre_p != 0`) carries `max(0, baseline)`.
pub fn quick_pileup(
    start_poss: &[i32],
    end_poss: &[i32],
    scale_factor: f32,
    baseline_value: f32,
) -> Vec<PosVal> {
    let l = start_poss.len();
    debug_assert_eq!(l, end_poss.len());
    // The C code dereferences start_poss[0]/end_poss[0] unconditionally; MACS
    // never pileups an empty chromosome. We return an empty array to be safe,
    // matching `se_all_in_one_pileup`'s `shape[0] == 0` early return.
    if l == 0 {
        return Vec::new();
    }

    let mut out: Vec<PosVal> = Vec::with_capacity(2 * l);

    let mut pileup: i32 = 0;
    let mut pre_p: i32 = start_poss[0].min(end_poss[0]);

    if pre_p != 0 {
        // first chunk of zeros: max(0, baseline)
        out.push(PosVal::new(pre_p, 0.0_f32.max(baseline_value)));
    }

    let mut i_s: usize = 0;
    let mut i_e: usize = 0;

    while i_s < l && i_e < l {
        let s = start_poss[i_s];
        let e = end_poss[i_e];
        if s < e {
            let p = s;
            if p != pre_p {
                out.push(PosVal::new(
                    p,
                    (pileup as f32 * scale_factor).max(baseline_value),
                ));
                pre_p = p;
            }
            pileup += 1;
            i_s += 1;
        } else if s > e {
            let p = e;
            if p != pre_p {
                out.push(PosVal::new(
                    p,
                    (pileup as f32 * scale_factor).max(baseline_value),
                ));
                pre_p = p;
            }
            pileup -= 1;
            i_e += 1;
        } else {
            i_s += 1;
            i_e += 1;
        }
    }

    // add the rest of end positions
    while i_e < l {
        let p = end_poss[i_e];
        if p != pre_p {
            out.push(PosVal::new(
                p,
                (pileup as f32 * scale_factor).max(baseline_value),
            ));
            pre_p = p;
        }
        pileup -= 1;
        i_e += 1;
    }

    out
}

/// Element-wise maximum of two `PosVal` arrays, mirroring
/// `cPosValCalculation.c::max_over_two_pv_array`.
///
/// Both inputs are bedGraph-ordered `PosVal` arrays. The merge walks them by
/// `pos`, emitting at each breakpoint the per-position maximum of the two
/// values. Note the C loop ends when *either* array is exhausted (it does not
/// drain a longer tail), which we replicate.
pub fn max_over_two_pv_array(pva1: &[PosVal], pva2: &[PosVal]) -> Vec<PosVal> {
    let l1 = pva1.len();
    let l2 = pva2.len();
    let mut out: Vec<PosVal> = Vec::with_capacity(l1 + l2);

    let mut i1: usize = 0;
    let mut i2: usize = 0;

    while i1 < l1 && i2 < l2 {
        let a = pva1[i1];
        let b = pva2[i2];
        let v = a.val.max(b.val);
        if a.pos < b.pos {
            out.push(PosVal::new(a.pos, v));
            i1 += 1;
        } else if a.pos > b.pos {
            out.push(PosVal::new(b.pos, v));
            i2 += 1;
        } else {
            out.push(PosVal::new(a.pos, v));
            i1 += 1;
            i2 += 1;
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reconstruct `(start, end, value)` triples from a bedGraph-style `PosVal`
    /// array the way `test_Pileup.py` does (pre starts at 0).
    fn to_triples(pv: &[PosVal]) -> Vec<(i32, i32, f32)> {
        let mut pre = 0;
        let mut out = Vec::new();
        for e in pv {
            out.push((pre, e.pos, e.val));
            pre = e.pos;
        }
        out
    }

    /// Golden values from `test_Pileup.py::Test_SE_Pileup::test_pileup_1`.
    #[test]
    fn se_pileup_golden_1() {
        let plus = [0, 1, 3];
        let minus = [8, 9, 10];
        // five_shift=0, three_shift=5, scale=0.5, baseline=0, rlength=100
        let pv = single_end_pileup(&plus, &minus, 0, 5, 0, 100, 0.5, 0.0);
        let expect = vec![
            (0, 1, 0.5),
            (1, 3, 1.0),
            (3, 4, 2.0),
            (4, 6, 2.5),
            (6, 8, 2.0),
            (8, 9, 1.0),
            (9, 10, 0.5),
        ];
        assert_eq!(to_triples(&pv), expect);
    }

    /// Golden values from `test_Pileup.py::Test_SE_Pileup::test_pileup_2`.
    #[test]
    fn se_pileup_golden_2() {
        let plus = [0, 1, 3];
        let minus = [8, 9, 10];
        // five_shift=0, three_shift=10, scale=2, baseline=8, rlength=100
        let pv = single_end_pileup(&plus, &minus, 0, 10, 0, 100, 2.0, 8.0);
        let expect = vec![
            (0, 1, 8.0),
            (1, 3, 10.0),
            (3, 8, 12.0),
            (8, 9, 10.0),
            (9, 10, 8.0),
            (10, 11, 8.0),
            (11, 13, 8.0),
        ];
        assert_eq!(to_triples(&pv), expect);
    }

    /// Golden values from `test_Pileup.py::Test_Quick_Pileup::test_pileup_1`.
    #[test]
    fn quick_pileup_golden() {
        let starts = [0, 1, 3, 3, 4, 5];
        let ends = [5, 6, 8, 8, 9, 10];
        let pv = quick_pileup(&starts, &ends, 0.5, 0.0);
        let expect = vec![
            (0, 1, 0.5),
            (1, 3, 1.0),
            (3, 4, 2.0),
            (4, 6, 2.5),
            (6, 8, 2.0),
            (8, 9, 1.0),
            (9, 10, 0.5),
        ];
        assert_eq!(to_triples(&pv), expect);
    }

    /// Golden values from `test_Pileup.py::Test_Over_Two_PV_Array::test_max`.
    /// `max_over_two_pv_array` corresponds to `over_two_pv_array(func="max")`.
    #[test]
    fn max_over_two_golden() {
        let pv1 = [
            PosVal::new(2, 1.0),
            PosVal::new(5, 2.0),
            PosVal::new(7, 3.0),
            PosVal::new(8, 4.0),
            PosVal::new(9, 3.0),
            PosVal::new(12, 2.0),
        ];
        let pv2 = [
            PosVal::new(1, 5.0),
            PosVal::new(4, 3.0),
            PosVal::new(6, 2.0),
            PosVal::new(8, 1.0),
            PosVal::new(10, 0.0),
            PosVal::new(11, 3.0),
        ];
        let pv = max_over_two_pv_array(&pv1, &pv2);
        let expect = vec![
            (0, 1, 5.0),
            (1, 2, 3.0),
            (2, 4, 3.0),
            (4, 5, 2.0),
            (5, 6, 3.0),
            (6, 7, 3.0),
            (7, 8, 4.0),
            (8, 9, 3.0),
            (9, 10, 2.0),
            (10, 11, 3.0),
        ];
        assert_eq!(to_triples(&pv), expect);
    }
}
