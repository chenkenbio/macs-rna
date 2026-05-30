//! Two-condition differential score track, ported from
//! `MACS3/Signal/ScoreTrack.py::TwoConditionScores` (used by `bdgdiff`).
//!
//! Per chromosome MACS stores four parallel arrays: `pos` (`int32`), and three
//! `float32` log-likelihood-ratio scores: `t1 vs c1`, `t2 vs c2`, `t1 vs t2`.
//! We hold them in [`ChromTwoCond`], keyed by raw chromosome bytes in a
//! `BTreeMap` for bytewise-sorted iteration.
//!
//! Phase 0 owns the struct, constructor and accessors. Phase 1C fills `build`
//! (aligning the four input `bedGraphTrackI`s and computing the three logLR
//! columns), `call_peaks`, and `write_bedGraph` (which uses the `>=1e-6` break
//! predicate from `bedgraph_io`).

use std::collections::BTreeMap;

use crate::bedgraph::BedGraphTrack;
use crate::peak_io::PeakIO;

/// `LOG10_E = log10(e)` as a `cython.float` (`ScoreTrack.py::LOG10_E`).
const LOG10_E: f32 = 0.43429448190325176_f32;

/// Asymmetric `log10` likelihood ratio (`ScoreTrack.py::logLR_asym`).
/// The `log`-bearing expression evaluates in `double` then narrows to f32.
fn log_lr_asym(x: f32, y: f32) -> f32 {
    let xd = x as f64;
    let yd = y as f64;
    let log10_e = LOG10_E as f64;
    let s: f64 = if x > y {
        (xd * (xd.ln() - yd.ln()) + yd - xd) * log10_e
    } else if x < y {
        (xd * (-xd.ln() + yd.ln()) - yd + xd) * log10_e
    } else {
        0.0
    };
    s as f32
}

/// Symmetric `log10` likelihood ratio (`ScoreTrack.py::logLR_sym`).
fn log_lr_sym(x: f32, y: f32) -> f32 {
    let xd = x as f64;
    let yd = y as f64;
    let log10_e = LOG10_E as f64;
    let s: f64 = if x > y {
        (xd * (xd.ln() - yd.ln()) + yd - xd) * log10_e
    } else if y > x {
        (yd * (xd.ln() - yd.ln()) + yd - xd) * log10_e
    } else {
        0.0
    };
    s as f32
}

/// One chromosome's parallel differential-score arrays
/// (`TwoConditionScores.data[chrom]`).
#[derive(Debug, Clone, Default)]
pub struct ChromTwoCond {
    /// Right-open interval end coordinates (C `int32`).
    pub pos: Vec<i32>,
    /// logLR of condition-1 treatment vs control (`t1 vs c1`).
    pub llr_1: Vec<f32>,
    /// logLR of condition-2 treatment vs control (`t2 vs c2`).
    pub llr_2: Vec<f32>,
    /// symmetric logLR of condition-1 vs condition-2 treatment (`t1 vs t2`).
    pub llr_3: Vec<f32>,
}

/// Differential score track (`ScoreTrack.py::TwoConditionScores`).
///
/// Holds the four source bedGraph tracks plus the computed scores. The source
/// tracks are owned so `build` can iterate them; Phase 1C decides whether to
/// keep or drop them after building.
///
/// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
#[derive(Debug, Clone)]
pub struct TwoConditionScores {
    /// Per-chromosome computed scores, bytewise-sorted by chromosome name.
    data: BTreeMap<Vec<u8>, ChromTwoCond>,
    /// Scaling factor applied to condition-1 tracks (`cond1_factor`).
    pub cond1_factor: f32,
    /// Scaling factor applied to condition-2 tracks (`cond2_factor`).
    pub cond2_factor: f32,
    /// Pseudocount added before taking logs (`pseudocount`).
    pub pseudocount: f32,
    /// Cutoff most recently used by `call_peaks` (`cutoff`).
    pub cutoff: f32,
    /// Condition-1 treatment bedGraph (`t1bdg`).
    pub t1bdg: BedGraphTrack,
    /// Condition-1 control bedGraph (`c1bdg`).
    pub c1bdg: BedGraphTrack,
    /// Condition-2 treatment bedGraph (`t2bdg`).
    pub t2bdg: BedGraphTrack,
    /// Condition-2 control bedGraph (`c2bdg`).
    pub c2bdg: BedGraphTrack,
}

impl TwoConditionScores {
    /// Create a differential track from the four input bedGraphs
    /// (`TwoConditionScores.__init__`).
    pub fn new(
        t1bdg: BedGraphTrack,
        c1bdg: BedGraphTrack,
        t2bdg: BedGraphTrack,
        c2bdg: BedGraphTrack,
        cond1_factor: f32,
        cond2_factor: f32,
        pseudocount: f32,
    ) -> Self {
        TwoConditionScores {
            data: BTreeMap::new(),
            cond1_factor,
            cond2_factor,
            pseudocount,
            cutoff: 0.0,
            t1bdg,
            c1bdg,
            t2bdg,
            c2bdg,
        }
    }

    /// Update the pseudocount (`set_pseudocount`).
    pub fn set_pseudocount(&mut self, pseudocount: f32) {
        self.pseudocount = pseudocount;
    }

    /// Immutable access to a chromosome's score arrays (`get_data_by_chr`).
    pub fn get_data_by_chr(&self, chromosome: &[u8]) -> Option<&ChromTwoCond> {
        self.data.get(chromosome)
    }

    /// All chromosome names, bytewise-sorted (`get_chr_names`).
    pub fn get_chr_names(&self) -> Vec<&[u8]> {
        self.data.keys().map(|k| k.as_slice()).collect()
    }

    /// Append a precomputed interval (used by [`build`](Self::build) once it is
    /// implemented). Exposed so Phase 1C can populate `data` from this file.
    pub fn push_interval(&mut self, chrom: &[u8], endpos: i32, llr1: f32, llr2: f32, llr3: f32) {
        let c = self.data.entry(chrom.to_vec()).or_default();
        c.pos.push(endpos);
        c.llr_1.push(llr1);
        c.llr_2.push(llr2);
        c.llr_3.push(llr3);
    }

    /// Compute one interval's three logLR columns and append it
    /// (`ScoreTrack.py::TwoConditionScores.add`). `t1`/`c1`/`t2`/`c2` are the
    /// raw pileup values; the pseudocount is added and each pair scaled by its
    /// condition factor before the log-likelihood ratios are taken.
    fn add(&mut self, chrom: &[u8], endpos: i32, t1: f32, c1: f32, t2: f32, c2: f32) {
        let pc = self.pseudocount;
        let f1 = self.cond1_factor;
        let f2 = self.cond2_factor;
        let llr1 = log_lr_asym((t1 + pc) * f1, (c1 + pc) * f1);
        let llr2 = log_lr_asym((t2 + pc) * f2, (c2 + pc) * f2);
        let llr3 = log_lr_sym((t1 + pc) * f1, (t2 + pc) * f2);
        self.push_interval(chrom, endpos, llr1, llr2, llr3);
    }

    // ------------------------------------------------------------------
    // Build / IO — Phase 1C.
    // ------------------------------------------------------------------

    /// Align the four bedGraphs and compute the three logLR columns
    /// (`ScoreTrack.py::TwoConditionScores.build` + `build_chromosome`).
    ///
    /// For each chromosome common to all four tracks, a 4-way two-pointer merge
    /// emits one interval per boundary. **Subtle MACS detail:** the source adds
    /// the interval with `endpos = pre_p` (the *previous* boundary) carrying the
    /// *current* values, then sets `pre_p = minp` and advances. The first add
    /// therefore uses `endpos = 0`, and the final boundary is **never emitted**
    /// because advancing the last pointer raises `StopIteration` before the next
    /// `add`. We replicate this exactly (off-by-one is intentional).
    pub fn build(&mut self) {
        let common = self.common_chrs();
        for chrom in common {
            let (c1tp, c1tv, c1cp, c1cv, c2tp, c2tv, c2cp, c2cv) = {
                let d1t = self.t1bdg.get_data_by_chr(&chrom).expect("t1 chrom");
                let d1c = self.c1bdg.get_data_by_chr(&chrom).expect("c1 chrom");
                let d2t = self.t2bdg.get_data_by_chr(&chrom).expect("t2 chrom");
                let d2c = self.c2bdg.get_data_by_chr(&chrom).expect("c2 chrom");
                (
                    d1t.pos.clone(), d1t.val.clone(),
                    d1c.pos.clone(), d1c.val.clone(),
                    d2t.pos.clone(), d2t.val.clone(),
                    d2c.pos.clone(), d2c.val.clone(),
                )
            };

            if c1tp.is_empty() || c1cp.is_empty() || c2tp.is_empty() || c2cp.is_empty() {
                continue;
            }

            // Indices into the four position/value lists.
            let (mut i1t, mut i1c, mut i2t, mut i2c) = (0usize, 0usize, 0usize, 0usize);
            let mut pre_p: i32 = 0;
            loop {
                let p1t = c1tp[i1t];
                let p1c = c1cp[i1c];
                let p2t = c2tp[i2t];
                let p2c = c2cp[i2c];
                let v1t = c1tv[i1t];
                let v1c = c1cv[i1c];
                let v2t = c2tv[i2t];
                let v2c = c2cv[i2c];

                let minp = p1t.min(p1c).min(p2t).min(p2c);
                // add at the PREVIOUS boundary with the CURRENT values.
                self.add(&chrom, pre_p, v1t, v1c, v2t, v2c);
                pre_p = minp;

                let mut done = false;
                if p1t == minp {
                    i1t += 1;
                    if i1t >= c1tp.len() { done = true; }
                }
                if p1c == minp {
                    i1c += 1;
                    if i1c >= c1cp.len() { done = true; }
                }
                if p2t == minp {
                    i2t += 1;
                    if i2t >= c2tp.len() { done = true; }
                }
                if p2c == minp {
                    i2c += 1;
                    if i2c >= c2cp.len() { done = true; }
                }
                if done {
                    break;
                }
            }
        }
    }

    /// Chromosome names common to all four input tracks, bytewise-sorted.
    fn common_chrs(&self) -> Vec<Vec<u8>> {
        use std::collections::BTreeSet;
        let s1: BTreeSet<&[u8]> = self.t1bdg.get_chr_names().into_iter().collect();
        let s2: BTreeSet<&[u8]> = self.c1bdg.get_chr_names().into_iter().collect();
        let s3: BTreeSet<&[u8]> = self.t2bdg.get_chr_names().into_iter().collect();
        let s4: BTreeSet<&[u8]> = self.c2bdg.get_chr_names().into_iter().collect();
        s1.iter()
            .filter(|c| s2.contains(*c) && s3.contains(*c) && s4.contains(*c))
            .map(|c| c.to_vec())
            .collect()
    }

    /// Call differential peaks above `cutoff` (`ScoreTrack.py::call_peaks`).
    ///
    /// Returns three [`PeakIO`]s: condition-1-unique (`t1 vs t2 >= cutoff` and
    /// `t1 vs c1 >= cutoff`), condition-2-unique (`t1 vs t2 <= -cutoff` and
    /// `t2 vs c2 >= cutoff`), and common (`|t1 vs t2| <= cutoff` with both
    /// conditions significant). Each peak's score is the length-weighted mean
    /// LLR for its category.
    pub fn call_peaks(
        &mut self,
        cutoff: f32,
        min_length: i32,
        max_gap: i32,
    ) -> (PeakIO, PeakIO, PeakIO) {
        let mut cat1_peaks = PeakIO::default();
        let mut cat2_peaks = PeakIO::default();
        let mut cat3_peaks = PeakIO::default();
        self.cutoff = cutoff;

        let chroms: Vec<Vec<u8>> = self.data.keys().cloned().collect();
        for chrom in &chroms {
            let c = self.data.get(chrom.as_slice()).expect("chrom present");
            let pos = &c.pos;
            let t1c1 = &c.llr_1;
            let t2c2 = &c.llr_2;
            let t1t2 = &c.llr_3;
            let ln = pos.len();

            // category index lists, mirroring the numpy boolean masks.
            let mut cat1: Vec<usize> = Vec::new();
            let mut cat2: Vec<usize> = Vec::new();
            let mut cat3: Vec<usize> = Vec::new();
            for i in 0..ln {
                let cond1_over = t1t2[i] >= cutoff;
                let cond2_over = t1t2[i] <= -cutoff;
                let cond_equal = t1t2[i] >= -cutoff && t1t2[i] <= cutoff;
                let cond1_sig = t1c1[i] >= cutoff;
                let cond2_sig = t2c2[i] >= cutoff;
                if cond1_sig && cond1_over {
                    cat1.push(i);
                }
                if cond2_over && cond2_sig {
                    cat2.push(i);
                }
                if cond1_sig && cond2_sig && cond_equal {
                    cat3.push(i);
                }
            }

            // cat1 scores: t1_vs_t2; cat2: -t1_vs_t2; cat3: |t1_vs_t2|.
            Self::add_a_peak(&mut cat1_peaks, chrom, &cat1, pos, t1t2, max_gap, min_length, |v| v);
            Self::add_a_peak(&mut cat2_peaks, chrom, &cat2, pos, t1t2, max_gap, min_length, |v| -v);
            Self::add_a_peak(&mut cat3_peaks, chrom, &cat3, pos, t1t2, max_gap, min_length, |v| v.abs());
        }

        (cat1_peaks, cat2_peaks, cat3_peaks)
    }

    /// Merge nearby significant intervals into peaks and add them with a
    /// length-weighted mean LLR (`ScoreTrack.py::__add_a_peak` +
    /// `mean_from_peakcontent`). `score_fn` maps the raw `t1_vs_t2` value to the
    /// category-specific score (identity, negation, or absolute value).
    #[allow(clippy::too_many_arguments)]
    fn add_a_peak<F: Fn(f32) -> f32>(
        peaks: &mut PeakIO,
        chrom: &[u8],
        indices: &[usize],
        pos: &[i32],
        t1t2: &[f32],
        max_gap: i32,
        min_length: i32,
        score_fn: F,
    ) {
        if indices.is_empty() {
            return;
        }
        // peak_content: (start, end, score)
        let mut peak_content: Vec<(i32, i32, f32)> = Vec::new();

        // start = pos[idx-1]; MACS clamps the very first index to 0 when it is 0
        // (otherwise pos[-1] would wrap to the last element).
        let first_idx = indices[0];
        let first_start = if first_idx == 0 { 0 } else { pos[first_idx - 1] };
        peak_content.push((first_start, pos[first_idx], score_fn(t1t2[first_idx])));
        for &idx in indices.iter().skip(1) {
            let startpos = pos[idx - 1];
            if startpos - peak_content.last().unwrap().1 <= max_gap {
                peak_content.push((startpos, pos[idx], score_fn(t1t2[idx])));
            } else {
                Self::close_diff_peak(peaks, chrom, &peak_content, min_length);
                peak_content.clear();
                peak_content.push((startpos, pos[idx], score_fn(t1t2[idx])));
            }
        }
        if !peak_content.is_empty() {
            Self::close_diff_peak(peaks, chrom, &peak_content, min_length);
        }
    }

    /// Emit one differential peak if it is at least `min_length` long, with the
    /// length-weighted mean LLR (`mean_from_peakcontent`).
    fn close_diff_peak(
        peaks: &mut PeakIO,
        chrom: &[u8],
        peak_content: &[(i32, i32, f32)],
        min_length: i32,
    ) {
        let first = peak_content[0];
        let last = peak_content[peak_content.len() - 1];
        if last.1 - first.0 < min_length {
            return;
        }
        let mean_llr = Self::mean_from_peakcontent(peak_content);
        peaks.add(chrom, first.0, last.1, -1, mean_llr, 0.0, 0.0, 0.0, 0.0, b"");
    }

    /// Length-weighted mean LLR (`ScoreTrack.py::mean_from_peakcontent`).
    ///
    /// MACS accumulates `sum_v += tmp_v * (tmp_e - tmp_s)` where `tmp_v` is the
    /// score **cast to `cython.long`** (truncated toward zero), then returns
    /// `cast(float, sum_v / ln)`. Cython's default `cdivision=False` makes the
    /// `long / long` a **Python true division** (a float), so we divide as f64
    /// (matching CPython's `int / int` -> `double`) and narrow to f32.
    fn mean_from_peakcontent(peak_content: &[(i32, i32, f32)]) -> f32 {
        let mut sum_v: i64 = 0;
        let mut ln: i64 = 0;
        for &(s, e, v) in peak_content {
            let tmp_v = v as i64; // cython.long cast truncates toward zero
            sum_v += tmp_v * (e - s) as i64;
            ln += (e - s) as i64;
        }
        (sum_v as f64 / ln as f64) as f32
    }

    /// Write column `column` (1=`t1 vs c1`, 2=`t2 vs c2`, 3=`t1 vs t2`) as
    /// bedGraph using the `>=1e-6` break predicate
    /// (`ScoreTrack.py::TwoConditionScores.write_bedGraph`).
    ///
    /// MACS comments out the track line for this writer, so `name`/`description`
    /// are accepted (for API parity) but unused.
    pub fn write_bedgraph<W: std::io::Write>(
        &self,
        fhd: &mut W,
        _name: &str,
        _description: &str,
        column: u8,
    ) -> std::io::Result<()> {
        assert!((1..=3).contains(&column), "column should be between 1, 2 or 3.");
        for (chrom, c) in &self.data {
            let value: &[f32] = match column {
                1 => &c.llr_1,
                2 => &c.llr_2,
                _ => &c.llr_3,
            };
            crate::bedgraph_io::write_bedgraph_predicate(
                chrom,
                &c.pos,
                value,
                &crate::bedgraph_io::Ge1e6,
                fhd,
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn push_and_accessors() {
        let mut t = TwoConditionScores::new(
            BedGraphTrack::new(0.0),
            BedGraphTrack::new(0.0),
            BedGraphTrack::new(0.0),
            BedGraphTrack::new(0.0),
            1.0,
            1.0,
            0.01,
        );
        t.push_interval(b"chr1", 100, 1.0, 2.0, 0.5);
        let c = t.get_data_by_chr(b"chr1").unwrap();
        assert_eq!(c.pos, vec![100]);
        assert_eq!(c.llr_1, vec![1.0_f32]);
        assert_eq!(c.llr_2, vec![2.0_f32]);
        assert_eq!(c.llr_3, vec![0.5_f32]);
        assert_eq!(t.get_chr_names(), vec![&b"chr1"[..]]);
    }

    /// Build a TwoConditionScores from the `test_ScoreTrack.py` regions and
    /// cross-check against MACS v3.0.4 ground truth (computed via the macs3v304
    /// env): `pos=[0,70,75,80,90]` and the three logLR columns; the final 150
    /// boundary is intentionally dropped (MACS off-by-one in build_chromosome).
    ///
    /// The logLR columns are compared to within `1e-4` absolute tolerance. A few
    /// values differ by exactly 1 ULP from compiled MACS because MACS calls the
    /// C `log` directly (`libc.math.log`), which can differ by 1 ULP from the
    /// libm `ln` Rust/CPython use; the diff (~5e-7) is far inside the harness's
    /// 1e-5 bedGraph tolerance.
    #[test]
    fn golden_build_and_call_peaks() {
        let regions1 = [(0, 70, 0.00_f32, 0.01_f32), (70, 80, 7.00, 0.5), (80, 150, 0.00, 0.02)];
        let regions2 = [(0, 75, 20.0_f32, 4.00_f32), (75, 90, 35.0, 6.00), (90, 150, 10.0, 15.00)];
        let mut t1 = BedGraphTrack::new(0.0);
        let mut c1 = BedGraphTrack::new(0.0);
        let mut t2 = BedGraphTrack::new(0.0);
        let mut c2 = BedGraphTrack::new(0.0);
        for &(s, e, tv, cv) in &regions1 {
            t1.add_loc(b"chrY", s, e, tv);
            c1.add_loc(b"chrY", s, e, cv);
        }
        for &(s, e, tv, cv) in &regions2 {
            t2.add_loc(b"chrY", s, e, tv);
            c2.add_loc(b"chrY", s, e, cv);
        }
        let mut tc = TwoConditionScores::new(t1, c1, t2, c2, 1.0, 1.0, 0.01);
        tc.build();

        let d = tc.get_data_by_chr(b"chrY").unwrap();
        assert_eq!(d.pos, vec![0, 70, 75, 80, 90]);
        // Expected logLR values from MACS v3.0.4 (macs3v304 env), as f32.
        let approx = |got: &[f32], want: &[f32]| {
            assert_eq!(got.len(), want.len());
            for (i, (&g, &w)) in got.iter().zip(want).enumerate() {
                assert!(
                    (g as f64 - w as f64).abs() <= 1e-4,
                    "col elem {i}: got {g} want {w}"
                );
            }
        };
        approx(&d.llr_1, &[-0.00133, 5.1555, 5.1555, -0.00391, -0.00391]);
        approx(&d.llr_2, &[7.02032, 7.02032, 14.19923, 14.19923, -0.41025]);
        approx(&d.llr_3, &[-57.37207, -3.46931, -12.29333, -108.88185, -25.6914]);

        // call_peaks: only cond2-unique (cat2) is non-empty:
        //   chrY 0 80, score 10.125
        let (cat1, cat2, cat3) = tc.call_peaks(3.0, 10, 10);
        assert!(cat1.is_empty());
        assert!(cat3.is_empty());
        let p2 = cat2.peaks_by_chr(b"chrY").expect("cat2 chrY");
        assert_eq!(p2.len(), 1);
        assert_eq!((p2[0].start, p2[0].end), (0, 80));
        assert_eq!(p2[0].score, 10.125_f32);
    }
}
