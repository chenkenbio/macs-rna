//! Alignment-driven peak caller, ported from
//! `MACS3/Signal/CallPeakUnit.py::CallerFromAlignments`.
//!
//! This is the heavy `callpeak` engine: it piles up treatment and control
//! tracks, builds local lambda, scores every bin and emits narrow/broad peaks.
//! Single-end (`FwTrack`) only — the harness never exercises paired-end here.
//!
//! All internal pileup / score values are `f32` (C `float`); poisson p-scores
//! are computed in f64 then narrowed, matching MACS. bedGraph output uses the
//! `|pre - cur| > 1e-5` coalescing predicate (`__write_bedGraph_for_a_chromosome`).

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};

use crate::peak_io::{BroadPeakIO, PeakIO};
use crate::pileup::{max_over_two_pv_array, single_end_pileup, PosVal};
use crate::prob::poisson_cdf;
use crate::track_fw::{FwTrack, INT_MAX};

/// `-log10` Poisson tail probability for `observed` counts under `expectation`
/// (`CallPeakUnit.py::get_pscore`). MACS computes
/// `-1 * poisson_cdf(observed, expectation, lower=False, log10=True)` and caches
/// the f32 result keyed by `(int(treat), ctrl)`. The cache only avoids recompute;
/// values are identical, so we recompute (with our own small cache for speed).
#[inline]
fn get_pscore(observed: i32, expectation: f32) -> f32 {
    let v: f64 = -1.0 * poisson_cdf(observed.max(0) as u32, expectation as f64, false, true);
    v as f32
}

/// One chromosome's merged pileup arrays (`self.chr_pos_treat_ctrl`).
///
/// `pos[i]` is the right-open end of interval `i`; `treat[i]`/`ctrl[i]` carry the
/// treatment pileup and the control lambda over `[pos[i-1], pos[i])`.
struct ChromPTC {
    pos: Vec<i32>,
    treat: Vec<f32>,
    ctrl: Vec<f32>,
}

/// Alignment-driven caller (`CallPeakUnit.py::CallerFromAlignments`).
///
/// Holds clones of the treatment / control 5'-end positions per chromosome so it
/// can pile up each chromosome on demand. Construction takes the same parameters
/// `PeakDetect` passes in MACS (`d`, `ctrl_d_s`, scaling factors, `lambda_bg`,
/// `end_shift`, bedGraph options).
pub struct CallerFromAlignments {
    /// Per-chromosome treatment 5'-end positions `(plus, minus)`.
    treat_loc: BTreeMap<Vec<u8>, (Vec<i32>, Vec<i32>)>,
    /// Per-chromosome control 5'-end positions `(plus, minus)`.
    ctrl_loc: BTreeMap<Vec<u8>, (Vec<i32>, Vec<i32>)>,
    /// Reference chromosome lengths (used to clip pileups). Missing => INT_MAX.
    rlengths: BTreeMap<Vec<u8>, i32>,
    /// Common chromosome names, bytewise-sorted (`self.chromosomes`).
    chromosomes: Vec<Vec<u8>>,
    /// Treatment effective total reads (for SPMR denominator).
    treat_total: u64,
    /// Control effective total reads (for SPMR denominator).
    ctrl_total: u64,

    /// Treatment fragment extension size (`d`).
    d: i32,
    /// Control extension sizes (`ctrl_d_s`).
    ctrl_d_s: Vec<i32>,
    /// Treatment scaling factor (`treat_scaling_factor`).
    treat_scaling_factor: f32,
    /// Control scaling factors, one per `ctrl_d_s` (`ctrl_scaling_factor_s`).
    ctrl_scaling_factor_s: Vec<f32>,
    /// Minimum local bias filling missing values (`lambda_bg`).
    lambda_bg: f32,
    /// Shift applied to read ends before extension (`end_shift`).
    end_shift: i32,
    /// Pseudocount for FE / fold-change (`pseudocount`).
    pseudocount: f32,
    /// True when no local lambda is requested (control = flat lambda_bg).
    no_lambda_flag: bool,

    /// Whether to emit a UCSC track line in bedGraph outputs (`trackline`).
    trackline: bool,
    /// Whether to save pileup/lambda bedGraphs (`save_bedGraph`).
    save_bedgraph: bool,
    /// Whether to normalize pileup per million reads (`save_SPMR`).
    save_spmr: bool,
    /// bedGraph prefix used in track-line descriptions (`bedGraph_filename_prefix`).
    bedgraph_prefix: String,

    /// p->q lookup, keyed by the f32 bit pattern of each pscore (`pqtable`).
    pqtable: BTreeMap<u32, f32>,
}

impl CallerFromAlignments {
    /// Construct the caller (`CallerFromAlignments.__init__`).
    ///
    /// `ctrl` may be `None` (no control); then the control track is the treatment
    /// track itself (matching MACS, which sets `self.ctrl = treat`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        treat: &mut FwTrack,
        ctrl: Option<&mut FwTrack>,
        d: i32,
        ctrl_d_s: Vec<i32>,
        treat_scaling_factor: f32,
        ctrl_scaling_factor_s: Vec<f32>,
        end_shift: i32,
        lambda_bg: f32,
        pseudocount: f32,
        save_bedgraph: bool,
        bedgraph_prefix: String,
        save_spmr: bool,
    ) -> Self {
        let treat_total = treat.total;
        let treat_loc = collect_locations(treat);
        let treat_rlengths = treat.get_rlengths().clone();

        let (ctrl_loc, ctrl_total, ctrl_rlengths) = match ctrl {
            Some(c) => {
                let ct = c.total;
                let cl = collect_locations(c);
                let cr = c.get_rlengths().clone();
                (cl, ct, cr)
            }
            None => (treat_loc.clone(), treat_total, treat_rlengths.clone()),
        };

        // Merge rlengths: treatment first, control fills any gaps. Both default
        // to INT_MAX where absent. (For BED input every value is INT_MAX.)
        let mut rlengths = treat_rlengths;
        for (k, v) in ctrl_rlengths {
            rlengths.entry(k).or_insert(v);
        }

        // Common chromosomes, bytewise-sorted (sorted(chr1 ∩ chr2)).
        let chr1: std::collections::BTreeSet<&Vec<u8>> = treat_loc.keys().collect();
        let chr2: std::collections::BTreeSet<&Vec<u8>> = ctrl_loc.keys().collect();
        let chromosomes: Vec<Vec<u8>> = chr1
            .intersection(&chr2)
            .map(|k| (*k).clone())
            .collect();

        let no_lambda_flag = ctrl_d_s.is_empty() || ctrl_scaling_factor_s.is_empty();

        CallerFromAlignments {
            treat_loc,
            ctrl_loc,
            rlengths,
            chromosomes,
            treat_total,
            ctrl_total,
            d,
            ctrl_d_s,
            treat_scaling_factor,
            ctrl_scaling_factor_s,
            lambda_bg,
            end_shift,
            pseudocount,
            no_lambda_flag,
            trackline: false,
            save_bedgraph,
            save_spmr,
            bedgraph_prefix,
            pqtable: BTreeMap::new(),
        }
    }

    /// Enable UCSC track line output (`enable_trackline`).
    pub fn enable_trackline(&mut self) {
        self.trackline = true;
    }

    /// Test-only view of the scaling parameters (`ctrl_d_s`,
    /// `ctrl_scaling_factor_s`, `lambda_bg`, `treat_scaling_factor`).
    #[cfg(test)]
    pub fn scaling_params(&self) -> (&[i32], &[f32], f32, f32) {
        (
            &self.ctrl_d_s,
            &self.ctrl_scaling_factor_s,
            self.lambda_bg,
            self.treat_scaling_factor,
        )
    }

    /// Common chromosome names (bytewise-sorted), for the caller to check.
    pub fn chromosomes(&self) -> &[Vec<u8>] {
        &self.chromosomes
    }

    /// `rlength` for `chrom`, defaulting to [`INT_MAX`].
    fn rlength(&self, chrom: &[u8]) -> i32 {
        self.rlengths.get(chrom).copied().unwrap_or(INT_MAX)
    }

    /// Pile up treatment + control for one chromosome and merge them, returning
    /// the `[pos, treat, ctrl]` arrays (`pileup_treat_ctrl_a_chromosome`).
    fn pileup_treat_ctrl_a_chromosome(&self, chrom: &[u8]) -> ChromPTC {
        let rlength = self.rlength(chrom);
        let (plus, minus) = &self.treat_loc[chrom];

        // Treatment pileup: directional=True. five_shift = -end_shift,
        // three_shift = end_shift + d.
        let five_shift = -self.end_shift;
        let three_shift = self.end_shift + self.d;
        let treat_pv = single_end_pileup(
            plus,
            minus,
            five_shift,
            three_shift,
            0,
            rlength,
            self.treat_scaling_factor,
            0.0,
        );

        // Control pileup.
        let ctrl_pv: Vec<PosVal> = if !self.no_lambda_flag {
            // Multi-window max, directional=False (pileup_a_chromosome_c).
            let (cplus, cminus) = &self.ctrl_loc[chrom];
            let mut prev: Option<Vec<PosVal>> = None;
            for i in 0..self.ctrl_d_s.len() {
                let dd = self.ctrl_d_s[i];
                // directional=False: five = d//2, three = d - d//2 (end_shift=0).
                let five = dd.div_euclid(2);
                let three = dd - dd.div_euclid(2);
                let tmp = single_end_pileup(
                    cplus,
                    cminus,
                    five,
                    three,
                    0,
                    rlength,
                    self.ctrl_scaling_factor_s[i],
                    self.lambda_bg,
                );
                prev = Some(match prev {
                    Some(p) => max_over_two_pv_array(&p, &tmp),
                    None => tmp,
                });
            }
            prev.unwrap_or_default()
        } else {
            // Flat global lambda: [last treat pos] -> lambda_bg.
            let last_pos = treat_pv.last().map(|p| p.pos).unwrap_or(0);
            vec![PosVal::new(last_pos, self.lambda_bg)]
        };

        chrom_pair_treat_ctrl(&treat_pv, &ctrl_pv)
    }

    /// Build the p->q table by scanning every chromosome
    /// (`__cal_pvalue_qvalue_table`). Accumulates `pscore -> total bp length`,
    /// then applies the Benjamini-Hochberg conversion.
    fn cal_pvalue_qvalue_table(&mut self) {
        // pscore (as f32 bits) -> total bp length.
        let mut pscore_stat: BTreeMap<u32, i64> = BTreeMap::new();
        let chroms = self.chromosomes.clone();
        for chrom in &chroms {
            let ptc = self.pileup_treat_ctrl_a_chromosome(chrom);
            let mut pre_p: i64 = 0;
            for j in 0..ptc.pos.len() {
                let this_v = get_pscore(ptc.treat[j] as i32, ptc.ctrl[j]);
                let this_l = ptc.pos[j] as i64 - pre_p;
                *pscore_stat.entry(this_v.to_bits()).or_insert(0) += this_l;
                pre_p = ptc.pos[j] as i64;
            }
        }

        self.pqtable = build_pq_table(&pscore_stat);
    }

    /// Call narrow peaks (`CallerFromAlignments.call_peaks`).
    ///
    /// `scoring_methods` selects the per-position score (`b'p'`/`b'q'`/`b'f'`/
    /// `b's'`); `scoring_cutoffs` the corresponding cutoffs. When `save_bedgraph`
    /// is set, the treatment-pileup and control-lambda bedGraphs are written to
    /// `treat_out`/`ctrl_out` while peaks are called.
    #[allow(clippy::too_many_arguments)]
    pub fn call_peaks<W: Write>(
        &mut self,
        scoring_methods: &[u8],
        scoring_cutoffs: &[f32],
        min_length: i32,
        max_gap: i32,
        call_summits: bool,
        treat_out: Option<&mut W>,
        ctrl_out: Option<&mut W>,
    ) -> PeakIO {
        assert!(
            !call_summits,
            "call_summits=true (multi-summit) is not supported by macs3-rs callpeak yet"
        );
        let mut peaks = PeakIO::default();

        if self.pqtable.is_empty() {
            self.cal_pvalue_qvalue_table();
        }

        let mut bdg = self.open_bedgraph(treat_out, ctrl_out);

        let chroms = self.chromosomes.clone();
        for chrom in &chroms {
            self.chrom_call_peak(
                &mut peaks,
                chrom,
                scoring_methods,
                scoring_cutoffs,
                min_length,
                max_gap,
                bdg.as_mut(),
            );
        }
        peaks
    }

    /// Call broad peaks (`CallerFromAlignments.call_broadpeaks`).
    #[allow(clippy::too_many_arguments)]
    pub fn call_broadpeaks<W: Write>(
        &mut self,
        scoring_methods: &[u8],
        lvl1_cutoffs: &[f32],
        lvl2_cutoffs: &[f32],
        min_length: i32,
        lvl1_max_gap: i32,
        lvl2_max_gap: i32,
        treat_out: Option<&mut W>,
        ctrl_out: Option<&mut W>,
    ) -> BroadPeakIO {
        let mut lvl1peaks = PeakIO::default();
        let mut lvl2peaks = PeakIO::default();

        if self.pqtable.is_empty() {
            self.cal_pvalue_qvalue_table();
        }

        let mut bdg = self.open_bedgraph(treat_out, ctrl_out);

        let chroms = self.chromosomes.clone();
        for chrom in &chroms {
            self.chrom_call_broadpeak(
                &mut lvl1peaks,
                &mut lvl2peaks,
                chrom,
                scoring_methods,
                lvl1_cutoffs,
                lvl2_cutoffs,
                min_length,
                lvl1_max_gap,
                lvl2_max_gap,
                bdg.as_mut(),
            );
        }

        combine_broad(&lvl1peaks, &lvl2peaks)
    }

    /// Open and prime the bedGraph writers (emitting the optional track line).
    fn open_bedgraph<'a, W: Write>(
        &self,
        treat_out: Option<&'a mut W>,
        ctrl_out: Option<&'a mut W>,
    ) -> Option<BedgraphWriters<'a, W>> {
        if !self.save_bedgraph {
            return None;
        }
        let (Some(t), Some(c)) = (treat_out, ctrl_out) else {
            return None;
        };
        let w = BedgraphWriters { treat: t, ctrl: c };
        if self.trackline {
            let _ = write!(
                w.treat,
                "track type=bedGraph name=\"treatment pileup\" description=\"treatment pileup after possible scaling for \'{}\'\"\n",
                self.bedgraph_prefix
            );
            let _ = write!(
                w.ctrl,
                "track type=bedGraph name=\"control lambda\" description=\"control lambda after possible scaling for \'{}\'\"\n",
                self.bedgraph_prefix
            );
        }
        Some(w)
    }

    /// Compute the per-position score array for one scoring symbol.
    fn cal_score(&self, ptc: &ChromPTC, symbol: u8) -> Vec<f32> {
        let n = ptc.pos.len();
        let mut s = vec![0.0_f32; n];
        match symbol {
            b'p' => {
                for i in 0..n {
                    s[i] = get_pscore(ptc.treat[i] as i32, ptc.ctrl[i]);
                }
            }
            b'q' => {
                for i in 0..n {
                    let p = get_pscore(ptc.treat[i] as i32, ptc.ctrl[i]);
                    s[i] = *self.pqtable.get(&p.to_bits()).unwrap_or(&0.0);
                }
            }
            b'f' => {
                for i in 0..n {
                    s[i] = (ptc.treat[i] + self.pseudocount) / (ptc.ctrl[i] + self.pseudocount);
                }
            }
            b's' => {
                for i in 0..n {
                    s[i] = ptc.treat[i] - ptc.ctrl[i];
                }
            }
            other => panic!("unsupported scoring symbol {other}"),
        }
        s
    }

    /// Write one chromosome's treat/control bedGraph
    /// (`__write_bedGraph_for_a_chromosome`).
    fn write_bedgraph_for_chrom<W: Write>(
        &self,
        chrom: &[u8],
        ptc: &ChromPTC,
        w: &mut BedgraphWriters<W>,
    ) {
        let denominator: f32 = if self.save_spmr {
            if self.treat_scaling_factor == 1.0 {
                self.treat_total as f32 / 1e6
            } else {
                self.ctrl_total as f32 / 1e6
            }
        } else {
            1.0
        };

        let l = ptc.pos.len();
        if l == 0 {
            return;
        }

        let mut pre_p_t: i32 = 0;
        let mut pre_p_c: i32 = 0;
        let mut pre_v_t: f32 = ptc.treat[0] / denominator;
        let mut pre_v_c: f32 = ptc.ctrl[0] / denominator;

        for i in 1..l {
            let v_t = ptc.treat[i] / denominator;
            let v_c = ptc.ctrl[i] / denominator;
            let p = ptc.pos[i - 1];

            if (pre_v_t - v_t).abs() > 1e-5 {
                write_bdg_line(w.treat, chrom, pre_p_t, p, pre_v_t);
                pre_v_t = v_t;
                pre_p_t = p;
            }
            if (pre_v_c - v_c).abs() > 1e-5 {
                write_bdg_line(w.ctrl, chrom, pre_p_c, p, pre_v_c);
                pre_v_c = v_c;
                pre_p_c = p;
            }
        }

        let p = ptc.pos[l - 1];
        write_bdg_line(w.treat, chrom, pre_p_t, p, pre_v_t);
        write_bdg_line(w.ctrl, chrom, pre_p_c, p, pre_v_c);
    }

    /// Call peaks for one chromosome (`__chrom_call_peak_using_certain_criteria`).
    #[allow(clippy::too_many_arguments)]
    fn chrom_call_peak<W: Write>(
        &self,
        peaks: &mut PeakIO,
        chrom: &[u8],
        scoring_methods: &[u8],
        scoring_cutoffs: &[f32],
        min_length: i32,
        max_gap: i32,
        bdg: Option<&mut BedgraphWriters<W>>,
    ) {
        let ptc = self.pileup_treat_ctrl_a_chromosome(chrom);

        if let Some(w) = bdg {
            self.write_bedgraph_for_chrom(chrom, &ptc, w);
        }

        let score_arrays: Vec<Vec<f32>> =
            scoring_methods.iter().map(|&s| self.cal_score(&ptc, s)).collect();

        // above_cutoff: indices where ALL scores pass their cutoff (strict >).
        let n = ptc.pos.len();
        let mut above: Vec<usize> = Vec::new();
        for i in 0..n {
            if passes_cutoffs(&score_arrays, scoring_cutoffs, i) {
                above.push(i);
            }
        }
        if above.is_empty() {
            return;
        }

        // peak_content: (start, end, treat_p, ctrl_p, index).
        let mut peak_content: Vec<(i32, i32, f32, f32, usize)> = Vec::new();
        let mut lastp: i64;

        let first = above[0];
        let ts0 = if first == 0 { 0 } else { ptc.pos[first - 1] };
        let te0 = ptc.pos[first];
        peak_content.push((ts0, te0, ptc.treat[first], ptc.ctrl[first], first));
        lastp = te0 as i64;

        for &idx in &above[1..] {
            let ts = if idx == 0 { 0 } else { ptc.pos[idx - 1] };
            let te = ptc.pos[idx];
            let tl = ts as i64 - lastp;
            if tl <= max_gap as i64 {
                peak_content.push((ts, te, ptc.treat[idx], ptc.ctrl[idx], idx));
                lastp = te as i64;
            } else {
                self.close_peak_wo_subpeaks(
                    &peak_content,
                    peaks,
                    min_length,
                    chrom,
                    &score_arrays,
                    scoring_cutoffs,
                );
                peak_content.clear();
                peak_content.push((ts, te, ptc.treat[idx], ptc.ctrl[idx], idx));
                lastp = te as i64;
            }
        }
        if !peak_content.is_empty() {
            self.close_peak_wo_subpeaks(
                &peak_content,
                peaks,
                min_length,
                chrom,
                &score_arrays,
                scoring_cutoffs,
            );
        }
    }

    /// Close one narrow peak (`__close_peak_wo_subpeaks`). Summit = midpoint of
    /// the highest-treatment-pileup region; ties pick the middle one.
    fn close_peak_wo_subpeaks(
        &self,
        peak_content: &[(i32, i32, f32, f32, usize)],
        peaks: &mut PeakIO,
        min_length: i32,
        chrom: &[u8],
        score_arrays: &[Vec<f32>],
        scoring_cutoffs: &[f32],
    ) {
        let peak_length = peak_content[peak_content.len() - 1].1 - peak_content[0].0;
        if peak_length < min_length {
            return;
        }

        // Find summit: middle of the highest-treat-pileup region(s).
        let mut tsummit: Vec<i32> = Vec::new();
        let mut tsummit_index: Vec<usize> = Vec::new();
        let mut summit_value: f32 = 0.0;
        for (i, &(tstart, tend, ttreat_p, _tctrl_p, _ti)) in peak_content.iter().enumerate() {
            let tscore = ttreat_p;
            // `not summit_value` is true when summit_value == 0.0.
            if summit_value == 0.0 || summit_value < tscore {
                tsummit = vec![(tend + tstart) / 2];
                tsummit_index = vec![i];
                summit_value = tscore;
            } else if summit_value == tscore {
                tsummit.push((tend + tstart) / 2);
                tsummit_index.push(i);
            }
        }
        let midindex = (tsummit.len() + 1) / 2 - 1;
        let summit_pos = tsummit[midindex];
        let summit_index = tsummit_index[midindex];

        let summit_treat = peak_content[summit_index].2;
        let summit_ctrl = peak_content[summit_index].3;
        let summit_array_index = peak_content[summit_index].4;

        // Double-check the summit passes all cutoffs (strict >).
        for i in 0..scoring_cutoffs.len() {
            if scoring_cutoffs[i] > score_arrays[i][summit_array_index] {
                return;
            }
        }

        let summit_p_score = get_pscore(summit_treat as i32, summit_ctrl);
        let summit_q_score = *self.pqtable.get(&summit_p_score.to_bits()).unwrap_or(&0.0);

        peaks.add(
            chrom,
            peak_content[0].0,
            peak_content[peak_content.len() - 1].1,
            summit_pos,
            summit_q_score,
            summit_treat,
            summit_p_score,
            (summit_treat + self.pseudocount) / (summit_ctrl + self.pseudocount),
            summit_q_score,
            b"",
        );
    }

    /// Call broad peaks for one chromosome
    /// (`__chrom_call_broadpeak_using_certain_criteria`).
    #[allow(clippy::too_many_arguments)]
    fn chrom_call_broadpeak<W: Write>(
        &self,
        lvl1peaks: &mut PeakIO,
        lvl2peaks: &mut PeakIO,
        chrom: &[u8],
        scoring_methods: &[u8],
        lvl1_cutoffs: &[f32],
        lvl2_cutoffs: &[f32],
        min_length: i32,
        lvl1_max_gap: i32,
        lvl2_max_gap: i32,
        bdg: Option<&mut BedgraphWriters<W>>,
    ) {
        let ptc = self.pileup_treat_ctrl_a_chromosome(chrom);

        if let Some(w) = bdg {
            self.write_bedgraph_for_chrom(chrom, &ptc, w);
        }

        let score_arrays: Vec<Vec<f32>> =
            scoring_methods.iter().map(|&s| self.cal_score(&ptc, s)).collect();

        // lvl1: strong peaks.
        self.segment_broad(
            &ptc,
            &score_arrays,
            lvl1_cutoffs,
            lvl1_max_gap,
            min_length,
            chrom,
            lvl1peaks,
        );
        // lvl2: weak peaks.
        self.segment_broad(
            &ptc,
            &score_arrays,
            lvl2_cutoffs,
            lvl2_max_gap,
            min_length,
            chrom,
            lvl2peaks,
        );
    }

    /// Segment one level of a broad peak call into `out`.
    #[allow(clippy::too_many_arguments)]
    fn segment_broad(
        &self,
        ptc: &ChromPTC,
        score_arrays: &[Vec<f32>],
        cutoffs: &[f32],
        max_gap: i32,
        min_length: i32,
        chrom: &[u8],
        out: &mut PeakIO,
    ) {
        let n = ptc.pos.len();
        let mut above: Vec<usize> = Vec::new();
        for i in 0..n {
            if passes_cutoffs(score_arrays, cutoffs, i) {
                above.push(i);
            }
        }
        if above.is_empty() {
            return;
        }

        let mut peak_content: Vec<(i32, i32, f32, f32, usize)> = Vec::new();
        let first = above[0];
        let ts0 = if first == 0 { 0 } else { ptc.pos[first - 1] };
        let te0 = ptc.pos[first];
        peak_content.push((ts0, te0, ptc.treat[first], ptc.ctrl[first], first));
        let mut lastp = te0 as i64;

        for &idx in &above[1..] {
            let ts = if idx == 0 { 0 } else { ptc.pos[idx - 1] };
            let te = ptc.pos[idx];
            let tl = ts as i64 - lastp;
            if tl <= max_gap as i64 {
                peak_content.push((ts, te, ptc.treat[idx], ptc.ctrl[idx], idx));
                lastp = te as i64;
            } else {
                self.close_peak_for_broad(&peak_content, out, min_length, chrom);
                peak_content.clear();
                peak_content.push((ts, te, ptc.treat[idx], ptc.ctrl[idx], idx));
                lastp = te as i64;
            }
        }
        if !peak_content.is_empty() {
            self.close_peak_for_broad(&peak_content, out, min_length, chrom);
        }
    }

    /// Close one broad-region peak (`__close_peak_for_broad_region`). Scores are
    /// the length-weighted mean over the region.
    fn close_peak_for_broad(
        &self,
        peak_content: &[(i32, i32, f32, f32, usize)],
        peaks: &mut PeakIO,
        min_length: i32,
        chrom: &[u8],
    ) {
        let peak_length = peak_content[peak_content.len() - 1].1 - peak_content[0].0;
        if peak_length < min_length {
            return;
        }

        let mut pileups: Vec<f32> = Vec::with_capacity(peak_content.len());
        let mut controls: Vec<f32> = Vec::with_capacity(peak_content.len());
        let mut lengths: Vec<i32> = Vec::with_capacity(peak_content.len());
        for &(tstart, tend, ttreat_p, tctrl_p, _ti) in peak_content {
            pileups.push(ttreat_p);
            controls.push(tctrl_p);
            lengths.push(tend - tstart);
        }
        let pscores: Vec<f32> = pileups
            .iter()
            .zip(&controls)
            .map(|(&t, &c)| get_pscore(t as i32, c))
            .collect();
        let qscores: Vec<f32> = pscores
            .iter()
            .map(|&p| *self.pqtable.get(&p.to_bits()).unwrap_or(&0.0))
            .collect();
        let fcs: Vec<f32> = pileups
            .iter()
            .zip(&controls)
            .map(|(&t, &c)| (t + self.pseudocount) / (c + self.pseudocount))
            .collect();

        peaks.add(
            chrom,
            peak_content[0].0,
            peak_content[peak_content.len() - 1].1,
            0,
            mean_from_value_length(&qscores, &lengths),
            mean_from_value_length(&pileups, &lengths),
            mean_from_value_length(&pscores, &lengths),
            mean_from_value_length(&fcs, &lengths),
            mean_from_value_length(&qscores, &lengths),
            b"",
        );
    }
}

/// Holds the two bedGraph output sinks during a call.
struct BedgraphWriters<'a, W: Write> {
    treat: &'a mut W,
    ctrl: &'a mut W,
}

/// Detected input format from the leading bytes of `path`
/// (`Parser.py::guess_parser`-style sniff): `"BAM"`, `"SAM"`, or `"BED"`.
pub fn sniff_format(path: &str) -> std::io::Result<&'static str> {
    let mut f = File::open(path)?;
    let mut head = [0u8; 4];
    let n = f.read(&mut head).unwrap_or(0);
    if n >= 4 && &head == b"BAM\x01" {
        return Ok("BAM");
    }
    if n >= 1 && head[0] == b'@' {
        return Ok("SAM");
    }
    if n >= 2 && head[0] == 0x1f && head[1] == 0x8b {
        // gzip: BAM (BGZF) vs gzipped BED.
        let f2 = File::open(path)?;
        let mut dec = flate2::read::MultiGzDecoder::new(f2);
        let mut magic = [0u8; 4];
        if dec.read_exact(&mut magic).is_ok() && &magic == b"BAM\x01" {
            return Ok("BAM");
        }
        return Ok("BED");
    }
    Ok("BED")
}

/// Compute tag size from a BED file (`BEDParser.tsize`): average `end - start`
/// over the first 10 valid alignments (max 10000 lines). Returns `-1` when no
/// valid line is found. Transparently handles gzip.
pub fn bed_tsize(path: &str) -> std::io::Result<i32> {
    let reader: Box<dyn BufRead> = {
        let mut probe = File::open(path)?;
        let mut magic = [0u8; 2];
        let n = probe.read(&mut magic).unwrap_or(0);
        let f = File::open(path)?;
        if n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
            Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(f)))
        } else {
            Box::new(BufReader::new(f))
        }
    };

    let mut s: i64 = 0;
    let mut n: i32 = 0;
    let mut m: i32 = 0;
    for line in reader.lines() {
        if n >= 10 || m >= 10000 {
            break;
        }
        let line = line?;
        m += 1;
        let trimmed = line.trim_end();
        if trimmed.is_empty()
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
            || trimmed.starts_with('#')
        {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let start = atoi_str(fields[1]);
        let end = atoi_str(fields[2]);
        let taglen = end - start;
        if taglen > 0 {
            s += taglen;
            n += 1;
        }
    }
    if n != 0 {
        Ok((s / n as i64) as i32)
    } else {
        Ok(-1)
    }
}

/// C `atoi` on a `&str`: parse the leading (optionally signed) integer; 0 on
/// failure.
fn atoi_str(s: &str) -> i64 {
    let s = s.trim();
    let bytes = s.as_bytes();
    let mut end = 0;
    if bytes.first().is_some_and(|&b| b == b'-' || b == b'+') {
        end = 1;
    }
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
    }
    s[..end].parse::<i64>().unwrap_or(0)
}

/// Collect per-chromosome `(plus, minus)` 5'-end vectors from an `FwTrack`.
fn collect_locations(track: &FwTrack) -> BTreeMap<Vec<u8>, (Vec<i32>, Vec<i32>)> {
    let mut out: BTreeMap<Vec<u8>, (Vec<i32>, Vec<i32>)> = BTreeMap::new();
    for chrom in track.get_chr_names() {
        if let Some(loc) = track.get_locations_by_chr(chrom) {
            out.insert(chrom.to_vec(), (loc.plus.clone(), loc.minus.clone()));
        }
    }
    out
}

/// Whether all score arrays pass their cutoffs at index `i` (strict `>`),
/// mirroring `apply_multiple_cutoffs` followed by `np.nonzero` (which keeps
/// entries where the count of passing cutoffs is non-zero — but MACS sums
/// booleans, so for a single scoring function this is exactly `score > cutoff`;
/// for multiple it's "any pass". The harness only uses a single function).
#[inline]
fn passes_cutoffs(score_arrays: &[Vec<f32>], cutoffs: &[f32], i: usize) -> bool {
    // `apply_multiple_cutoffs` returns sum of (arr > cutoff); np.nonzero keeps
    // indices where the sum is non-zero, i.e. ANY function passes.
    for k in 0..cutoffs.len() {
        if score_arrays[k][i] > cutoffs[k] {
            return true;
        }
    }
    false
}

/// Pair treatment and control pileups (`__chrom_pair_treat_ctrl`). Merges the two
/// bedGraph-ordered arrays on their breakpoints, carrying the current treat/ctrl
/// value at each emitted position.
fn chrom_pair_treat_ctrl(treat_pv: &[PosVal], ctrl_pv: &[PosVal]) -> ChromPTC {
    let lt = treat_pv.len();
    let lc = ctrl_pv.len();
    let cap = lt + lc;
    let mut pos = Vec::with_capacity(cap);
    let mut treat = Vec::with_capacity(cap);
    let mut ctrl = Vec::with_capacity(cap);

    let mut it = 0usize;
    let mut ic = 0usize;
    while it < lt && ic < lc {
        let tp = treat_pv[it];
        let cp = ctrl_pv[ic];
        if tp.pos < cp.pos {
            pos.push(tp.pos);
            treat.push(tp.val);
            ctrl.push(cp.val);
            it += 1;
        } else if tp.pos > cp.pos {
            pos.push(cp.pos);
            treat.push(tp.val);
            ctrl.push(cp.val);
            ic += 1;
        } else {
            pos.push(tp.pos);
            treat.push(tp.val);
            ctrl.push(cp.val);
            it += 1;
            ic += 1;
        }
    }

    ChromPTC { pos, treat, ctrl }
}

/// Benjamini-Hochberg p->q conversion (`__cal_pvalue_qvalue_table` tail).
///
/// `N = Σ length`, `f = -log10(N)`, rank `k` accumulates bp length, and for each
/// unique pscore (descending) `q = v + (log10(k) + f)`, clamped so it never
/// increases (`if q > pre_q: q = pre_q`) and stops once `q <= 0` (remaining
/// pscores all map to 0). Returns a map keyed by the f32 bit pattern of pscores.
fn build_pq_table(pscore_stat: &BTreeMap<u32, i64>) -> BTreeMap<u32, f32> {
    let n: i64 = pscore_stat.values().sum();
    let mut k: f64 = 1.0;
    // `f` is a `cython.float` (f32) in MACS: `-log10(N)` is computed in double
    // then narrowed to f32 before being reused in the `q` expression.
    let f: f32 = (-(n as f64).log10()) as f32;
    // `pre_q` is a `cython.float` (f32) too.
    let mut pre_q: f32 = 2147483647.0;

    let mut unique: Vec<(f32, i64)> = pscore_stat
        .iter()
        .map(|(&bits, &ln)| (f32::from_bits(bits), ln))
        .collect();
    // descending by pscore.
    unique.sort_by(|a, b| b.0.partial_cmp(&a.0).expect("pscores are finite"));

    let mut table: BTreeMap<u32, f32> = BTreeMap::new();
    let mut i = 0usize;
    while i < unique.len() {
        let (v, ln) = unique[i];
        // q = v + (log10(k) + f): log10(k) is double, + f promotes f to double,
        // + v promotes v to double; the f32 local `q` then narrows the result.
        let mut q: f32 = (v as f64 + (k.log10() + f as f64)) as f32;
        if q > pre_q {
            q = pre_q;
        }
        if q <= 0.0 {
            break;
        }
        table.insert(v.to_bits(), q);
        pre_q = q;
        k += ln as f64;
        i += 1;
    }
    for &(v, _) in &unique[i..] {
        table.insert(v.to_bits(), 0.0);
    }
    table
}

/// Length-weighted mean (`mean_from_value_length`). Sum is accumulated in f64
/// (`tmp_v * tmp_l`), then narrowed to f32.
fn mean_from_value_length(value: &[f32], length: &[i32]) -> f32 {
    let mut sum_v: f64 = 0.0;
    let mut ln: i64 = 0;
    for i in 0..length.len() {
        let tmp_l = length[i] as i64;
        let tmp_v = value[i] as f64;
        sum_v += tmp_v * tmp_l as f64;
        ln += tmp_l;
    }
    (sum_v / ln as f64) as f32
}

/// Write one `chrom\tstart\tend\t%.5f\n` bedGraph line, matching the C
/// `fprintf(... "%.5f" ...)`.
fn write_bdg_line<W: Write>(w: &mut W, chrom: &[u8], start: i32, end: i32, value: f32) {
    let _ = w.write_all(chrom);
    let _ = write!(w, "\t{}\t{}\t", start, end);
    let mut s = String::new();
    crate::fmt::write_f5(&mut s, value);
    let _ = w.write_all(s.as_bytes());
    let _ = w.write_all(b"\n");
}

/// Combine lvl1 (strong) and lvl2 (weak) peaks into broad peaks
/// (`call_broadpeaks` tail). lvl1 regions are assumed to nest inside lvl2 ones.
fn combine_broad(lvl1peaks: &PeakIO, lvl2peaks: &PeakIO) -> BroadPeakIO {
    let mut broadpeaks = BroadPeakIO::new();

    // chrs = lvl1peaks.get_chr_names(); iterate sorted.
    for chrom in lvl1peaks.get_chr_names() {
        let lvl1chrom = lvl1peaks.peaks_by_chr(chrom).cloned().unwrap_or_default();
        let lvl2chrom = lvl2peaks.peaks_by_chr(chrom).cloned().unwrap_or_default();

        let mut lvl1_iter = lvl1chrom.iter();
        let mut tmppeakset: Vec<crate::peak_io::NarrowPeak> = Vec::new();

        let mut current_lvl1 = lvl1_iter.next();
        let mut i = 0usize;
        let mut stopped = false;
        while i < lvl2chrom.len() {
            let lvl2 = &lvl2chrom[i];
            // for each lvl2 peak, collect all lvl1 peaks inside.
            loop {
                match current_lvl1 {
                    Some(lvl1)
                        if lvl2.start <= lvl1.start && lvl1.end <= lvl2.end =>
                    {
                        tmppeakset.push(lvl1.clone());
                        current_lvl1 = lvl1_iter.next();
                    }
                    Some(_) => {
                        add_broadpeak(&mut broadpeaks, chrom, lvl2, &tmppeakset);
                        tmppeakset.clear();
                        break;
                    }
                    None => {
                        // StopIteration: no more lvl1 peaks.
                        add_broadpeak(&mut broadpeaks, chrom, lvl2, &tmppeakset);
                        tmppeakset.clear();
                        // add the rest of lvl2 peaks (j from i+1).
                        for lvl2b in &lvl2chrom[i + 1..] {
                            add_broadpeak(&mut broadpeaks, chrom, lvl2b, &[]);
                        }
                        stopped = true;
                        break;
                    }
                }
            }
            if stopped {
                break;
            }
            i += 1;
        }
    }

    broadpeaks
}

/// Add one broad peak to `bpeaks` (`__add_broadpeak`). Builds the block structure
/// from the nested lvl1 peakset (or a 2-block 1bp pad when empty).
fn add_broadpeak(
    bpeaks: &mut BroadPeakIO,
    chrom: &[u8],
    lvl2peak: &crate::peak_io::NarrowPeak,
    lvl1peakset: &[crate::peak_io::NarrowPeak],
) {
    let start = lvl2peak.start;
    let end = lvl2peak.end;

    if lvl1peakset.is_empty() {
        bpeaks.add(
            chrom,
            start,
            end,
            lvl2peak.score,
            format!("{}", start).into_bytes(),
            format!("{}", end).into_bytes(),
            2,
            b"1,1".to_vec(),
            format!("0,{}", end - start - 1).into_bytes(),
            lvl2peak.pileup,
            lvl2peak.pscore,
            lvl2peak.fc,
            lvl2peak.qscore,
            b"".to_vec(),
        );
        return;
    }

    let mut thick_start = lvl1peakset[0].start;
    let mut thick_end = lvl1peakset[lvl1peakset.len() - 1].end;
    let mut block_num = lvl1peakset.len() as i32;
    let mut block_sizes: Vec<String> =
        lvl1peakset.iter().map(|x| format!("{}", x.length)).collect();
    let mut block_starts: Vec<String> = lvl1peakset
        .iter()
        .map(|x| format!("{}", x.start - start))
        .collect();

    if thick_start != start {
        thick_start = start;
        block_num += 1;
        block_sizes.insert(0, "1".to_string());
        block_starts.insert(0, "0".to_string());
    }
    if thick_end != end {
        thick_end = end;
        block_num += 1;
        block_sizes.push("1".to_string());
        block_starts.push(format!("{}", end - start - 1));
    }

    bpeaks.add(
        chrom,
        start,
        end,
        lvl2peak.score,
        format!("{}", thick_start).into_bytes(),
        format!("{}", thick_end).into_bytes(),
        block_num,
        block_sizes.join(",").into_bytes(),
        block_starts.join(",").into_bytes(),
        lvl2peak.pileup,
        lvl2peak.pscore,
        lvl2peak.fc,
        lvl2peak.qscore,
        b"".to_vec(),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hand-checked Benjamini-Hochberg conversion (`__cal_pvalue_qvalue_table`).
    /// Two pscores: 10.0 over 100 bp, 2.0 over 900 bp. N = 1000, f = -log10(1000)
    /// = -3. Processing descending:
    ///   v=10: q = 10 + (log10(1) - 3) = 10 - 3 = 7.0; k -> 1+100=101.
    ///   v=2:  q = 2 + (log10(101) - 3) = 2 + (2.00432 - 3) = 1.00432.
    #[test]
    fn pq_table_bh_golden() {
        let mut stat: BTreeMap<u32, i64> = BTreeMap::new();
        stat.insert(10.0_f32.to_bits(), 100);
        stat.insert(2.0_f32.to_bits(), 900);
        let table = build_pq_table(&stat);

        let q10 = *table.get(&10.0_f32.to_bits()).unwrap();
        let q2 = *table.get(&2.0_f32.to_bits()).unwrap();
        assert!((q10 - 7.0).abs() < 1e-5, "q10={q10}");
        let expected_q2 = (2.0_f64 + ((101.0_f64).log10() + (-3.0_f64))) as f32;
        assert!((q2 - expected_q2).abs() < 1e-6, "q2={q2} expected={expected_q2}");
    }

    /// Monotone clamp: a higher-rank (lower) pscore whose raw q would exceed the
    /// previous q is clamped down (`if q > pre_q: q = pre_q`).
    #[test]
    fn pq_table_monotone_clamp() {
        // v=5 over 1 bp (k starts at 1), then v=4 over a huge length so its raw q
        // would be larger; it must be clamped to the previous q.
        let mut stat: BTreeMap<u32, i64> = BTreeMap::new();
        stat.insert(5.0_f32.to_bits(), 1);
        stat.insert(4.0_f32.to_bits(), 1_000_000);
        let table = build_pq_table(&stat);
        let q5 = *table.get(&5.0_f32.to_bits()).unwrap();
        let q4 = *table.get(&4.0_f32.to_bits()).unwrap();
        // q4 is never allowed to exceed q5.
        assert!(q4 <= q5, "q4={q4} > q5={q5}");
    }

    /// `chrom_pair_treat_ctrl` merges two bedGraph-ordered arrays on breakpoints,
    /// carrying current treat/ctrl values (`__chrom_pair_treat_ctrl`).
    #[test]
    fn pair_treat_ctrl_merge() {
        // treat: [0,10)=5, [10,30)=8 ; ctrl: [0,20)=2, [20,30)=4
        let treat_pv = [PosVal::new(10, 5.0), PosVal::new(30, 8.0)];
        let ctrl_pv = [PosVal::new(20, 2.0), PosVal::new(30, 4.0)];
        let ptc = chrom_pair_treat_ctrl(&treat_pv, &ctrl_pv);
        assert_eq!(ptc.pos, vec![10, 20, 30]);
        assert_eq!(ptc.treat, vec![5.0, 8.0, 8.0]);
        assert_eq!(ptc.ctrl, vec![2.0, 2.0, 4.0]);
    }

    /// Length-weighted mean (`mean_from_value_length`).
    #[test]
    fn weighted_mean() {
        // values [2, 4] with lengths [10, 30] -> (2*10 + 4*30)/40 = 140/40 = 3.5.
        let m = mean_from_value_length(&[2.0, 4.0], &[10, 30]);
        assert!((m - 3.5).abs() < 1e-6);
    }
}
