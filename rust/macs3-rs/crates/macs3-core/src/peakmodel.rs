//! Shifting-size model, ported from `MACS3/Signal/PeakModel.py::PeakModel`.
//!
//! `predictd` (and `callpeak` without `--nomodel`) builds a cross-correlation
//! model from paired plus/minus strand tag peaks to estimate the fragment size
//! `d`. Phase 0 fixes the field layout and the `build` entry point; Phase 1E
//! fills the cross-correlation algorithm (`__find_paired_peaks`,
//! `__naive_find_peaks`, `__find_pair_center`, `__paired_peak_model`).
//!
//! # Float policy
//! MACS pileups (`naive_quick_pileup`'s `ret_v`) are `float32`; the per-strand
//! profile lines (`plus_line`/`minus_line`) are `int32`. The cross-correlation
//! pipeline (`mean`/`std`/`np.correlate`/`smooth`) runs in numpy's default
//! `float64`, and the lag axis `xcorr = np.linspace(...)` is `float64`. The
//! winning lag is `d = int(xcorr[argmax])` — a C truncation toward zero of a
//! `float64`. We therefore run the whole correlation pipeline in `f64` and only
//! widen-then-narrow into the public `Vec<f32>` fields at the very end, so the
//! argmax (and hence `d`) is computed exactly as MACS does. No fast-math.

use crate::track_fw::FwTrack;

/// Raised when fewer than 100 paired plus/minus peaks are found
/// (`PeakModel.py::NotEnoughPairsException`).
#[derive(Debug, Clone)]
pub struct NotEnoughPairs {
    /// Diagnostic message.
    pub value: String,
}

impl std::fmt::Display for NotEnoughPairs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl std::error::Error for NotEnoughPairs {}

/// Tunable inputs for the model, gathered from the CLI options
/// (`PeakModel.__init__` reads these off `opt`).
#[derive(Debug, Clone, Copy)]
pub struct PeakModelOptions {
    /// Effective genome size (`opt.gsize`).
    pub gsize: f64,
    /// Upper mfold bound (`opt.umfold`).
    pub umfold: i32,
    /// Lower mfold bound (`opt.lmfold`).
    pub lmfold: i32,
    /// Band width (`opt.bw`).
    pub bw: i32,
    /// Minimum acceptable predicted fragment size (`opt.d_min`).
    pub d_min: i32,
}

/// Cross-correlation shifting model (`PeakModel.py::PeakModel`).
///
/// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
#[derive(Debug, Clone, Default)]
pub struct PeakModel {
    /// Effective genome size (`gz`).
    pub gz: f64,
    /// Maximum pairs used while building (`max_pairnum`).
    pub max_pairnum: i32,
    /// Upper mfold bound (`umfold`).
    pub umfold: i32,
    /// Lower mfold bound (`lmfold`).
    pub lmfold: i32,
    /// Band width (`bw`).
    pub bw: i32,
    /// Minimum acceptable predicted fragment size (`d_min`).
    pub d_min: i32,
    /// Each tag is expanded to this many bp centred on its cut point
    /// (`tag_expansion_size`, fixed at 10 in MACS).
    pub tag_expansion_size: i32,

    /// Lower tag-count bound per candidate peak (`min_tags`).
    pub min_tags: i32,
    /// Upper tag-count bound per candidate peak (`max_tags`).
    pub max_tags: i32,
    /// Candidate peak window size, `2 * bw` (`peaksize`).
    pub peaksize: i32,

    /// Estimated fragment size (`d`).
    pub d: i32,
    /// Scan window used downstream, derived from `d` (`scan_window`).
    pub scan_window: i32,
    /// Alternative `d` candidates from the correlation (`alternative_d`).
    pub alternative_d: Vec<i32>,

    /// Averaged plus-strand profile (`plus_line`).
    pub plus_line: Vec<f32>,
    /// Averaged minus-strand profile (`minus_line`).
    pub minus_line: Vec<f32>,
    /// Shifted combined profile (`shifted_line`).
    pub shifted_line: Vec<f32>,
    /// Correlation lag axis (`xcorr`).
    pub xcorr: Vec<f32>,
    /// Correlation value axis (`ycorr`).
    pub ycorr: Vec<f32>,

    /// Human-readable summary (`__str__` summary).
    pub summary: String,
}

impl PeakModel {
    /// Create a model from CLI options (`PeakModel.__init__`). `max_pairnum`
    /// defaults to 500 in MACS; the treatment track is supplied to
    /// [`build`](Self::build) by Phase 1E rather than stored here, to keep this
    /// foundation type free of a track-type dependency.
    pub fn new(opt: PeakModelOptions, max_pairnum: i32) -> Self {
        PeakModel {
            gz: opt.gsize,
            max_pairnum,
            umfold: opt.umfold,
            lmfold: opt.lmfold,
            bw: opt.bw,
            d_min: opt.d_min,
            tag_expansion_size: 10,
            ..Default::default()
        }
    }

    /// Build the cross-correlation model and set `d` / `scan_window`
    /// (`PeakModel.build`, `PeakModel.py:98`).
    ///
    /// `treatment` is the finalized [`FwTrack`] of 5' cut positions (the Cython
    /// `treatment` member). Steps mirror the source exactly:
    /// 1. compute `peaksize`, `min_tags`, `max_tags`;
    /// 2. find paired +/- strand peaks (`__find_paired_peaks`);
    /// 3. require >= 100 paired peaks (else [`NotEnoughPairs`]);
    /// 4. estimate `d` via cross-correlation (`__paired_peak_model`).
    ///
    /// Note: MACS 3.0.4 stores `max_pairnum` but never uses it inside `build`,
    /// so we do not cap the pair count (bit-exact with the source).
    pub fn build(&mut self, treatment: &FwTrack) -> Result<(), NotEnoughPairs> {
        // PeakModel.py:112-120
        self.peaksize = 2 * self.bw;
        // min/max unique hits on a single strand, decided by lmfold/umfold.
        // Python: int(round(float(total) * fold * peaksize / gz / 2)).
        // float(total)*fold*peaksize is f64; the chained `/` are f64; Python's
        // round() is banker's rounding (round-half-to-even).
        self.min_tags = py_round(
            treatment.total as f64 * self.lmfold as f64 * self.peaksize as f64 / self.gz / 2.0,
        ) as i32;
        self.max_tags = py_round(
            treatment.total as f64 * self.umfold as f64 * self.peaksize as f64 / self.gz / 2.0,
        ) as i32;

        // PeakModel.py:124-135
        let paired_peakpos = self.find_paired_peaks(treatment);
        let num_paired_peakpos: usize = paired_peakpos.iter().map(|(_, v)| v.len()).sum();

        if num_paired_peakpos < 100 {
            return Err(NotEnoughPairs {
                value: "No enough pairs to build model".to_string(),
            });
        }

        // PeakModel.py:138
        self.paired_peak_model(treatment, &paired_peakpos);
        Ok(())
    }

    /// `PeakModel.py::__str__` (`PeakModel.py:140`). Filled into `self.summary`
    /// by the CLI when needed.
    pub fn str_summary(&self) -> String {
        format!(
            "\nSummary of Peak Model:\n  Baseline: {}\n  Upperline: {}\n  Fragment size: {}\n  Scan window size: {}\n",
            self.min_tags, self.max_tags, self.d, self.scan_window
        )
    }

    // ------------------------------------------------------------------
    // __find_paired_peaks  (PeakModel.py:152)
    // ------------------------------------------------------------------

    /// Call paired +/- peaks per chromosome, returning, for each chromosome
    /// (in bytewise-sorted order, matching `chrs.sort()`), the list of paired
    /// peak centers (`__find_paired_peaks`).
    fn find_paired_peaks(&self, treatment: &FwTrack) -> Vec<(Vec<u8>, Vec<i32>)> {
        // `get_chr_names()` is already bytewise-sorted; the Python code does an
        // explicit `chrs.sort()` over the same keys.
        let mut out: Vec<(Vec<u8>, Vec<i32>)> = Vec::new();
        for chrom in treatment.get_chr_names() {
            let loc = match treatment.get_locations_by_chr(chrom) {
                Some(l) => l,
                None => continue,
            };
            let plus_peaksinfo = self.naive_find_peaks(&loc.plus);
            let minus_peaksinfo = self.naive_find_peaks(&loc.minus);
            // PeakModel.py:187-192 -- skip chrom if either strand has no peaks.
            if plus_peaksinfo.is_empty() || minus_peaksinfo.is_empty() {
                continue;
            }
            let centers = self.find_pair_center(&plus_peaksinfo, &minus_peaksinfo);
            out.push((chrom.to_vec(), centers));
        }
        out
    }

    // ------------------------------------------------------------------
    // __naive_find_peaks  (PeakModel.py:195)
    // ------------------------------------------------------------------

    /// Naively call peaks on a single strand's sorted tag list, returning
    /// `(summit_pos, height)` pairs (`__naive_find_peaks`).
    fn naive_find_peaks(&self, taglist: &[i32]) -> Vec<(i32, f32)> {
        // PeakModel.py:210 -- fewer than 2 tags -> no peaks.
        if taglist.len() < 2 {
            return Vec::new();
        }
        // PeakModel.py:214 -- extension = int(self.peaksize/2). peaksize is an
        // int; Python `int(int/2)` truncates toward zero (peaksize >= 0).
        let extension = self.peaksize / 2;
        let pv = naive_quick_pileup(taglist, extension);
        // PeakModel.py:215 -- naive_call_peaks(pv, min_tags, max_tags) with the
        // function's own defaults max_gap=50, min_length=200.
        naive_call_peaks(&pv, self.min_tags as f32, self.max_tags as f32, 50, 200)
    }

    // ------------------------------------------------------------------
    // __find_pair_center  (PeakModel.py:400)
    // ------------------------------------------------------------------

    /// Pair up plus/minus peak summits whose centers overlap within `peaksize`
    /// and whose tag counts are comparable, returning the midpoints
    /// (`__find_pair_center`).
    fn find_pair_center(
        &self,
        pluspeaks: &[(i32, f32)],
        minuspeaks: &[(i32, f32)],
    ) -> Vec<i32> {
        let mut ip: usize = 0; // index for plus peaks
        let mut im: usize = 0; // index for minus peaks
        let mut im_prev: usize = 0; // index for minus peaks in previous plus peak
        let mut pair_centers: Vec<i32> = Vec::new();
        let ip_max = pluspeaks.len();
        let im_max = minuspeaks.len();
        let mut flag_find_overlap = false;
        let peaksize = self.peaksize;
        while ip < ip_max && im < im_max {
            let (pp, pn) = pluspeaks[ip];
            let (mp, mn) = minuspeaks[im];
            if pp - peaksize > mp {
                // move minus
                im += 1;
            } else if pp + peaksize < mp {
                // move plus
                ip += 1;
                im = im_prev; // search minus peaks from previous index
                flag_find_overlap = false;
            } else {
                // overlap!
                if !flag_find_overlap {
                    flag_find_overlap = true;
                    im_prev = im; // only the first index is recorded
                }
                // number of tags in plus & minus peak regions are comparable.
                // Cython operands are C float; replicate in f32.
                let ratio = pn / mn;
                if ratio < 2.0 && ratio > 0.5 && pp < mp {
                    // Python: (pp+mp)//2 -- floor division of two ints.
                    pair_centers.push(floor_div2(pp + mp));
                }
                im += 1;
            }
        }
        pair_centers
    }

    // ------------------------------------------------------------------
    // __paired_peak_model  (PeakModel.py:221)
    // ------------------------------------------------------------------

    /// Build the per-strand profile lines from the paired centers and estimate
    /// `d` via cross-correlation (`__paired_peak_model`).
    fn paired_peak_model(
        &mut self,
        treatment: &FwTrack,
        paired_peakpos: &[(Vec<u8>, Vec<i32>)],
    ) {
        let window_size = (1 + 2 * self.peaksize + self.tag_expansion_size) as usize;
        // Accumulator arrays of difference deltas (start += / end -=), int32.
        let mut plus_start = vec![0i32; window_size];
        let mut plus_end = vec![0i32; window_size];
        let mut minus_start = vec![0i32; window_size];
        let mut minus_end = vec![0i32; window_size];

        for (chrom, centers) in paired_peakpos {
            let loc = match treatment.get_locations_by_chr(chrom) {
                Some(l) => l,
                None => continue,
            };
            self.model_add_line(centers, &loc.plus, &mut plus_start, &mut plus_end);
            self.model_add_line(centers, &loc.minus, &mut minus_start, &mut minus_end);
        }

        // __count: prefix-sum start+end into the int line.
        let plus_line_i = count_line(&plus_start, &plus_end);
        let minus_line_i = count_line(&minus_start, &minus_end);

        // --- Cross-correlation (numpy float64 throughout) ---
        let n = window_size;
        let plus_f: Vec<f64> = plus_line_i.iter().map(|&v| v as f64).collect();
        let minus_f: Vec<f64> = minus_line_i.iter().map(|&v| v as f64).collect();

        // normalize first: (line - mean)/(std*len)   [population std, ddof=0]
        let plus_data = normalize(&plus_f);
        let minus_data = normalize(&minus_f);

        // ycorr = np.correlate(minus_data, plus_data, "full")[ws-peaksize : ws+peaksize]
        // np.correlate(a, v, "full") == np.convolve(a, v[::-1], "full").
        let full = correlate_full(&minus_data, &plus_data);
        let lo = n - self.peaksize as usize;
        let hi = n + self.peaksize as usize;
        let ycorr_raw: Vec<f64> = full[lo..hi].to_vec();

        // xcorr = np.linspace(len//2*-1, len//2, num=len)  -- float64 axis.
        let len_y = ycorr_raw.len();
        let xcorr = linspace(
            -((len_y / 2) as f64),
            (len_y / 2) as f64,
            len_y,
        );

        // smooth correlation (flat window, default window_len=11).
        let ycorr = smooth_flat(&ycorr_raw, 11);

        // local maxima:
        //   i_l_max = np.r_[False, y[1:]>y[:-1]] & np.r_[y[:-1]>y[1:], False]
        let m = ycorr.len();
        let mut idx: Vec<usize> = Vec::new();
        for i in 0..m {
            let left = if i == 0 {
                false
            } else {
                ycorr[i] > ycorr[i - 1]
            };
            let right = if i + 1 >= m {
                false
            } else {
                ycorr[i] > ycorr[i + 1]
            };
            if left && right {
                idx.push(i);
            }
        }

        // i_l_max = i_l_max[xcorr[i_l_max] > d_min]
        idx.retain(|&i| xcorr[i] > self.d_min as f64);

        // i_l_max = i_l_max[argsort(ycorr[i_l_max])[::-1]]
        // numpy argsort is a stable ascending sort by default; [::-1] reverses
        // the whole permutation (so equal keys end up in reverse-of-original
        // order). Replicate stable-ascending-then-reverse.
        idx.sort_by(|&a, &b| {
            ycorr[a]
                .partial_cmp(&ycorr[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        idx.reverse();

        // alternative_d = sorted(int(xcorr[i]) for i in i_l_max)
        let mut alt: Vec<i32> = idx.iter().map(|&i| trunc_i32(xcorr[i])).collect();
        alt.sort_unstable();
        self.alternative_d = alt;

        // assert len(alternative_d) > 0 -- mirror MACS's assertion. If no peak
        // survives, leave d at 0; callers (predictd) treat empty alt as failure.
        // PeakModel.py:308 -- d = xcorr[i_l_max[0]] (cython.int truncation).
        if let Some(&best) = idx.first() {
            self.d = trunc_i32(xcorr[best]);
        }

        // store profiles / correlation as f32 (public field type).
        self.plus_line = plus_line_i.iter().map(|&v| v as f32).collect();
        self.minus_line = minus_line_i.iter().map(|&v| v as f32).collect();
        self.ycorr = ycorr.iter().map(|&v| v as f32).collect();
        self.xcorr = xcorr.iter().map(|&v| v as f32).collect();

        // scan_window = max(d, tag_expansion_size)*2
        self.scan_window = self.d.max(self.tag_expansion_size) * 2;
    }

    // ------------------------------------------------------------------
    // __model_add_line  (PeakModel.py:319)
    // ------------------------------------------------------------------

    /// Project each tag in `pos2` that falls inside a center's window onto the
    /// difference arrays `start`/`end` (`__model_add_line`).
    fn model_add_line(
        &self,
        pos1: &[i32],
        pos2: &[i32],
        start: &mut [i32],
        end: &mut [i32],
    ) {
        let mut i1: usize = 0;
        let mut i2: usize = 0;
        let mut i2_prev: usize = 0;
        let i1_max = pos1.len();
        let i2_max = pos2.len();
        let mut flag_find_overlap = false;
        let max_index = (start.len() - 1) as i32;
        // half window: peaksize + tag_expansion_size // 2 (floor div).
        let psize_adjusted1 = self.peaksize + self.tag_expansion_size / 2;
        let half_exp = self.tag_expansion_size / 2; // int(tes/2), tes>=0

        while i1 < i1_max && i2 < i2_max {
            let p1 = pos1[i1];
            let p2 = pos2[i2];
            if p1 - psize_adjusted1 > p2 {
                i2 += 1;
            } else if p1 + psize_adjusted1 < p2 {
                i1 += 1;
                i2 = i2_prev;
                flag_find_overlap = false;
            } else {
                // overlap!
                if !flag_find_overlap {
                    flag_find_overlap = true;
                    i2_prev = i2;
                }
                // s = max(int(p2 - tes/2 - p1 + psize_adjusted1), 0)
                // tes/2 is integer here (tag_expansion_size is int); int(...)
                // of an int is identity, so this is exact integer arithmetic.
                let s = (p2 - half_exp - p1 + psize_adjusted1).max(0);
                start[s as usize] += 1;
                let e = (p2 + half_exp - p1 + psize_adjusted1).min(max_index);
                end[e as usize] -= 1;
                i2 += 1;
            }
        }
    }
}

// ----------------------------------------------------------------------
// Free helpers replicating MACS3.Signal.Pileup functions used by PeakModel.
// ----------------------------------------------------------------------

/// `Pileup.py::naive_quick_pileup` (`Pileup.py:643`). Symmetric extension of
/// each (sorted) position by `extension`; returns coalesced `(pos, value)`
/// pairs. `value` is `float32` in MACS (`ret_v` dtype `f4`).
fn naive_quick_pileup(sorted_poss: &[i32], extension: i32) -> (Vec<i32>, Vec<f32>) {
    let l = sorted_poss.len();
    // MACS raises if l == 0, but callers guarantee l >= 2.
    let mut start_poss: Vec<i32> = Vec::with_capacity(l);
    for &p in sorted_poss {
        let s = p - extension;
        start_poss.push(if s < 0 { 0 } else { s });
    }
    let end_poss: Vec<i32> = sorted_poss.iter().map(|&p| p + extension).collect();

    let mut ret_p: Vec<i32> = Vec::with_capacity(2 * l);
    let mut ret_v: Vec<f32> = Vec::with_capacity(2 * l);

    let mut pileup: i32 = 0;
    let mut i_s: usize = 0;
    let mut i_e: usize = 0;

    let mut pre_p = start_poss[0].min(end_poss[0]);
    if pre_p != 0 {
        ret_p.push(pre_p);
        ret_v.push(0.0);
    }

    while i_s < l && i_e < l {
        let sp = start_poss[i_s];
        let ep = end_poss[i_e];
        if sp < ep {
            let p = sp;
            if p != pre_p {
                ret_p.push(p);
                ret_v.push(pileup as f32);
                pre_p = p;
            }
            pileup += 1;
            i_s += 1;
        } else if sp > ep {
            let p = ep;
            if p != pre_p {
                ret_p.push(p);
                ret_v.push(pileup as f32);
                pre_p = p;
            }
            pileup -= 1;
            i_e += 1;
        } else {
            i_s += 1;
            i_e += 1;
        }
    }

    // remaining end positions
    while i_e < l {
        let p = end_poss[i_e];
        if p != pre_p {
            ret_p.push(p);
            ret_v.push(pileup as f32);
            pre_p = p;
        }
        pileup -= 1;
        i_e += 1;
    }

    (ret_p, ret_v)
}

/// `Pileup.py::naive_call_peaks` (`Pileup.py:882`). Returns `(summit, height)`
/// pairs. `min_v`/`max_v` compared in `f64` (Cython reads `v` as `double`).
fn naive_call_peaks(
    pv: &(Vec<i32>, Vec<f32>),
    min_v: f32,
    max_v: f32,
    max_gap: i32,
    min_length: i32,
) -> Vec<(i32, f32)> {
    let (ps, vs) = pv;
    let n = ps.len();
    let mut ret_peaks: Vec<(i32, f32)> = Vec::new();
    // peak_content: (pre_p, p, v)
    let mut peak_content: Vec<(i32, i32, f32)> = Vec::new();

    let min_v64 = min_v as f64;
    let mut pre_p: i32 = 0;
    let mut x: usize = 0;

    // find the first region above min_v
    while x < n {
        let p = ps[x];
        let v = vs[x];
        x += 1;
        if (v as f64) > min_v64 {
            peak_content.push((pre_p, p, v));
            pre_p = p;
            break;
        } else {
            pre_p = p;
        }
    }

    // continue scanning the rest
    for i in x..n {
        let p = ps[i];
        let v = vs[i];
        if (v as f64) <= min_v64 {
            pre_p = p;
            continue;
        }
        // gap = pre_p - peak_content[-1].1
        let last_end = peak_content.last().unwrap().1;
        if pre_p - last_end <= max_gap {
            peak_content.push((pre_p, p, v));
        } else {
            let first_start = peak_content[0].0;
            let last_end2 = peak_content.last().unwrap().1;
            if last_end2 - first_start >= min_length {
                close_peak(&peak_content, &mut ret_peaks, max_v);
            }
            peak_content = vec![(pre_p, p, v)];
        }
        pre_p = p;
    }

    // save the last peak
    if !peak_content.is_empty() {
        let first_start = peak_content[0].0;
        let last_end = peak_content.last().unwrap().1;
        if last_end - first_start >= min_length {
            close_peak(&peak_content, &mut ret_peaks, max_v);
        }
    }

    ret_peaks
}

/// `Pileup.py::__close_peak` (`Pileup.py:966`). Finds the summit (midpoint of
/// the highest-value sub-region; ties resolved at the middle index) and pushes
/// `(summit, height)` if `height < max_v`.
fn close_peak(
    peak_content: &[(i32, i32, f32)],
    peaks: &mut Vec<(i32, f32)>,
    max_v: f32,
) {
    let mut tsummit: Vec<i32> = Vec::new();
    let mut summit_value: f32 = 0.0;
    for &(tstart, tend, tvalue) in peak_content {
        // `not summit_value` is true while summit_value == 0.0
        if summit_value == 0.0 || summit_value < tvalue {
            tsummit = vec![midpoint(tstart, tend)];
            summit_value = tvalue;
        } else if summit_value == tvalue {
            tsummit.push(midpoint(tstart, tend));
        }
    }
    // summit = tsummit[int((len(tsummit)+1)/2)-1]
    let idx = ((tsummit.len() + 1) / 2) - 1;
    let summit = tsummit[idx];
    if summit_value < max_v {
        peaks.push((summit, summit_value));
    }
}

/// `int((tend+tstart)/2)` -- Python true division then truncation toward zero.
#[inline]
fn midpoint(tstart: i32, tend: i32) -> i32 {
    trunc_i32((tend as f64 + tstart as f64) / 2.0)
}

/// `__count` (`PeakModel.py:384`): prefix-sum of `start[i] + end[i]`.
fn count_line(start: &[i32], end: &[i32]) -> Vec<i32> {
    let mut line = vec![0i32; start.len()];
    let mut pileup: i64 = 0;
    for i in 0..start.len() {
        pileup += start[i] as i64 + end[i] as i64;
        line[i] = pileup as i32;
    }
    line
}

/// `(line - line.mean())/(line.std()*len(line))` with numpy population std
/// (`ddof=0`). All in `float64`.
fn normalize(line: &[f64]) -> Vec<f64> {
    let n = line.len();
    if n == 0 {
        return Vec::new();
    }
    let mean = line.iter().sum::<f64>() / n as f64;
    let var = line.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / n as f64;
    let std = var.sqrt();
    let denom = std * n as f64;
    line.iter().map(|&v| (v - mean) / denom).collect()
}

/// `np.correlate(a, v, mode="full")` for real input, which equals
/// `np.convolve(a, v[::-1], mode="full")`. Output length `len(a)+len(v)-1`.
fn correlate_full(a: &[f64], v: &[f64]) -> Vec<f64> {
    let la = a.len();
    let lv = v.len();
    if la == 0 || lv == 0 {
        return Vec::new();
    }
    // vrev = v reversed
    let out_len = la + lv - 1;
    let mut out = vec![0.0f64; out_len];
    // correlate full: out[k] = sum_j a[k - (lv-1) + j] * v[j], j in [0,lv)
    // equivalently convolve(a, vrev): out[k] = sum_j a[j]*vrev[k-j]
    for k in 0..out_len {
        let mut acc = 0.0f64;
        // a index = j, vrev index = k - j; vrev[m] = v[lv-1-m]
        let j_lo = k.saturating_sub(lv - 1);
        let j_hi = (k + 1).min(la);
        for j in j_lo..j_hi {
            let m = k - j; // index into vrev
            acc += a[j] * v[lv - 1 - m];
        }
        out[k] = acc;
    }
    out
}

/// `np.linspace(start, stop, num)` with `endpoint=True` (numpy default). For
/// `num > 1`, `step = (stop-start)/(num-1)` and `out[num-1]` is forced to
/// `stop` exactly.
fn linspace(start: f64, stop: f64, num: usize) -> Vec<f64> {
    if num == 0 {
        return Vec::new();
    }
    if num == 1 {
        return vec![start];
    }
    let div = (num - 1) as f64;
    let step = (stop - start) / div;
    let mut out = Vec::with_capacity(num);
    for i in 0..num {
        out.push(start + step * i as f64);
    }
    // numpy sets the last sample exactly to `stop`.
    out[num - 1] = stop;
    out
}

/// `PeakModel.py::smooth(x, window_len=11, window="flat")` (`PeakModel.py:453`).
///
/// Reflected padding `s = [x[wl-1:0:-1], x, x[-1:-wl:-1]]`, flat (boxcar)
/// window `w = ones(wl)/wl`, `y = convolve(w, s, "valid")`, return
/// `y[wl//2 : len(y)-wl//2]`. Output length equals input length.
fn smooth_flat(x: &[f64], window_len: usize) -> Vec<f64> {
    // MACS short-circuits if window_len < 3.
    if window_len < 3 {
        return x.to_vec();
    }
    let wl = window_len;
    // s = np.r_[x[wl-1:0:-1], x, x[-1:-wl:-1]]
    let n = x.len();
    let mut s: Vec<f64> = Vec::with_capacity(2 * (wl - 1) + n);
    // x[wl-1 : 0 : -1]  -> indices wl-1, wl-2, ..., 1  (length wl-1)
    for i in (1..wl).rev() {
        s.push(x[i]);
    }
    // x
    s.extend_from_slice(x);
    // x[-1 : -wl : -1] -> indices n-1, n-2, ..., n-wl+1 (length wl-1)
    for k in 1..wl {
        s.push(x[n - k]);
    }

    // w = ones(wl)/wl
    let w = 1.0f64 / wl as f64;

    // y = np.convolve(w, s, mode="valid"); valid length = len(s) - wl + 1.
    // Since w is constant, each valid output is the mean of a wl-wide window.
    let s_len = s.len();
    let valid_len = s_len - wl + 1;
    let mut y = vec![0.0f64; valid_len];
    // running form would change rounding vs np.convolve's direct sum; numpy
    // computes a direct dot product per output, so do the same.
    for i in 0..valid_len {
        let mut acc = 0.0f64;
        for j in 0..wl {
            acc += w * s[i + j];
        }
        y[i] = acc;
    }

    // return y[wl//2 : len(y)-wl//2]
    let lo = wl / 2;
    let hi = y.len() - wl / 2;
    y[lo..hi].to_vec()
}

/// Truncate a `float64` toward zero into `i32` (C `(int)` cast / Cython
/// `cython.int`).
#[inline]
fn trunc_i32(x: f64) -> i32 {
    x.trunc() as i32
}

/// Floor division by 2 of an `i32` (`(a) // 2` in Python).
#[inline]
fn floor_div2(a: i32) -> i32 {
    a.div_euclid(2)
}

/// Python `round()` (round-half-to-even / banker's rounding) for a non-negative
/// `f64`, returning the nearest even integer on ties.
#[inline]
fn py_round(x: f64) -> i64 {
    // f64::round_ties_even is exactly Python 3's round() semantics for floats.
    x.round_ties_even() as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // Goldens ported from references/MACS/test/test_PeakModel.py
    // ------------------------------------------------------------------

    /// `test_build_raises_when_not_enough_pairs`: empty track -> NotEnoughPairs.
    #[test]
    fn build_raises_when_not_enough_pairs() {
        let opt = PeakModelOptions {
            gsize: 1_000_000.0,
            umfold: 30,
            lmfold: 10,
            bw: 100,
            d_min: 20,
        };
        let mut model = PeakModel::new(opt, 500);
        let mut fw = FwTrack::new(0, String::new(), 100_000);
        // DummyTreatment(total=1000) with a single empty chr1.
        fw.add_loc(b"chr1", 1, 0); // give it a chrom, then clear via finalize
        // Make it genuinely empty of usable pairs: 1 plus tag only -> no peaks.
        fw.finalize();
        fw.total = 1000; // match DummyTreatment.total used for thresholds
        let res = model.build(&fw);
        assert!(res.is_err());
    }

    /// `test_build_computes_thresholds_before_pairing`: peaksize / min_tags /
    /// max_tags are computed before the (failing) pairing step.
    #[test]
    fn build_computes_thresholds_before_pairing() {
        let opt = PeakModelOptions {
            gsize: 500_000.0,
            umfold: 15,
            lmfold: 5,
            bw: 120,
            d_min: 20,
        };
        let mut model = PeakModel::new(opt, 500);
        let mut fw = FwTrack::new(0, String::new(), 100_000);
        fw.add_loc(b"chr1", 1, 0);
        fw.finalize();
        fw.total = 300; // DummyTreatment(total=300)
        let _ = model.build(&fw);

        let peaksize = 2 * opt.bw;
        assert_eq!(model.peaksize, peaksize);
        let expected_min =
            (300.0_f64 * opt.lmfold as f64 * peaksize as f64 / opt.gsize / 2.0).round_ties_even()
                as i32;
        let expected_max =
            (300.0_f64 * opt.umfold as f64 * peaksize as f64 / opt.gsize / 2.0).round_ties_even()
                as i32;
        assert_eq!(model.min_tags, expected_min);
        assert_eq!(model.max_tags, expected_max);
    }

    /// `test_str_representation_reflects_summary_fields`.
    #[test]
    fn str_representation_reflects_summary_fields() {
        let opt = PeakModelOptions {
            gsize: 1_000_000.0,
            umfold: 30,
            lmfold: 10,
            bw: 100,
            d_min: 20,
        };
        let mut model = PeakModel::new(opt, 500);
        model.min_tags = 3;
        model.max_tags = 9;
        model.d = 147;
        model.scan_window = 300;
        let s = model.str_summary();
        assert!(s.contains("Baseline: 3"));
        assert!(s.contains("Upperline: 9"));
        assert!(s.contains("Fragment size: 147"));
        assert!(s.contains("Scan window size: 300"));
    }

    // ------------------------------------------------------------------
    // End-to-end golden against the real MACS3 PeakModel on a synthetic
    // FWTrack (seeded reproducible cluster layout, TRUE_D = 150). Computed
    // with macs3v304 env: d=149, min_tags=3, scan_window=298,
    // alternative_d=[149], plus_line len=411, plus_line sum=80000.
    // See report for the generation script.
    // ------------------------------------------------------------------

    /// Deterministic LCG matching the synthetic generator's tag jitter so the
    /// Rust track is byte-identical to the Python one fed to MACS3.
    fn build_synthetic_track() -> FwTrack {
        // We can't reuse numpy's RNG, so the golden was regenerated with this
        // exact deterministic jitter pattern (see report). Layout: 200 clusters,
        // spacing 5000, 40 tags/peak, jitter in [-8, 8], TRUE_D = 150.
        const N_CLUSTERS: i32 = 200;
        const SPACING: i32 = 5000;
        const TAGS: i32 = 40;
        const TRUE_D: i32 = 150;
        // Simple deterministic jitter: cycle through [-8..=8].
        let jitter: Vec<i32> = (-8..=8).collect();
        let mut fw = FwTrack::new(50, String::new(), 100_000);
        let mut ji = 0usize;
        for k in 0..N_CLUSTERS {
            let base = 10_000 + k * SPACING;
            for _ in 0..TAGS {
                let p = base + jitter[ji % jitter.len()];
                ji += 1;
                fw.add_loc(b"chr1", p, 0);
            }
            for _ in 0..TAGS {
                let m = base + TRUE_D + jitter[ji % jitter.len()];
                ji += 1;
                fw.add_loc(b"chr1", m, 1);
            }
        }
        fw.finalize();
        fw
    }

    /// End-to-end build on the deterministic synthetic track. All goldens were
    /// produced by running THIS EXACT track (same jitter cycle, same add order)
    /// through MACS3's `PeakModel.build` in the macs3v304 conda env:
    ///   d=149, min_tags=3, scan_window=298, alternative_d=[149],
    ///   plus_line len=411, plus_line sum=80000, minus_line sum=80000.
    #[test]
    fn build_recovers_fragment_size_synthetic() {
        let opt = PeakModelOptions {
            gsize: 1_000_000.0,
            umfold: 50,
            lmfold: 2,
            bw: 100,
            d_min: 10,
        };
        let mut model = PeakModel::new(opt, 1000);
        let fw = build_synthetic_track();
        model
            .build(&fw)
            .expect("synthetic track should build a model");

        assert_eq!(model.peaksize, 200);
        assert_eq!(model.min_tags, 3);
        assert_eq!(model.max_tags, 80);
        // plus_line length = 1 + 2*peaksize + tag_expansion_size = 411
        assert_eq!(model.plus_line.len(), 411);
        assert_eq!(model.plus_line.iter().sum::<f32>() as i64, 80_000);
        assert_eq!(model.minus_line.iter().sum::<f32>() as i64, 80_000);
        // EXACT golden d from MACS3.
        assert_eq!(model.d, 149);
        assert_eq!(model.scan_window, 298);
        assert_eq!(model.alternative_d, vec![149]);
    }
}
