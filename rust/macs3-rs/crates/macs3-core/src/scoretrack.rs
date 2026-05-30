//! Per-interval treatment/control score track, ported from
//! `MACS3/Signal/ScoreTrack.py::ScoreTrackII`.
//!
//! Per chromosome MACS stores four parallel arrays: `pos` (`int32` right-open
//! interval ends), `treat` (`float32` treatment pileup), `ctrl` (`float32`
//! control pileup) and `score` (`float32`, computed by `change_score_method`).
//! We hold them in parallel `Vec`s inside [`ChromScore`], keyed by raw
//! chromosome bytes in a `BTreeMap` so iteration is bytewise-sorted (matching
//! `sorted(self.data.keys())`).
//!
//! Phase 0 owns the struct, constructor, `add_chromosome`/`add`/`finalize`,
//! accessors and `get_chr_names`. Phase 1C / Phase 2 fill the scoring methods
//! (`ppois`/`qpois`/`FE`/`logFE`/`logLR`/`slogLR`/`subtract`/`max`),
//! `make_pq_table`, `call_peaks`, `call_broadpeaks`, `cutoff_analysis` and
//! `write_bedgraph` (which uses the `>1e-5` predicate from `bedgraph_io`).

use std::collections::BTreeMap;

use crate::bedgraph::BedGraphTrack;
use crate::peak_io::{BroadPeakIO, PeakIO};
use crate::prob::poisson_cdf;

/// `LOG10_E = log10(e)` as a `cython.float` constant (`ScoreTrack.py::LOG10_E`).
/// MACS stores this as a C `float`; in the likelihood expressions it is promoted
/// back to `double` (C arithmetic with the `double`-returning `log`), so we keep
/// the f32-rounded value but compute in f64 like the C runtime does.
const LOG10_E: f32 = 0.43429448190325176_f32;

/// `-log10` Poisson tail probability for `observed` counts under `expectation`
/// (`ScoreTrack.py::get_pscore`). MACS computes
/// `-1 * poisson_cdf(observed, expectation, lower=False, log10=True)` and caches
/// the f32 result keyed by `(observed, expectation)`. The cache only avoids
/// recomputation; values are identical with or without it, so we recompute.
///
/// `observed` is a C `int` in MACS; we take it as `i32` and pass it to
/// `poisson_cdf` as a `u32` (the count is non-negative in every call site).
fn get_pscore(observed: i32, expectation: f64) -> f32 {
    let score: f64 = -1.0 * poisson_cdf(observed as u32, expectation, false, true);
    score as f32
}

/// Asymmetric `log10` likelihood ratio between `x` and `y`
/// (`ScoreTrack.py::logLR_asym`). The C expression promotes the `cython.float`
/// operands to `double` (because `log` returns `double`), evaluates in `double`,
/// then narrows the result to the `cython.float` local `s`. We mirror that:
/// widen `x`/`y` to f64, compute, then `as f32`.
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

/// Symmetric `log10` likelihood ratio between `x` and `y`
/// (`ScoreTrack.py::logLR_sym`). Same `double` promotion as [`log_lr_asym`].
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

/// One chromosome's parallel score arrays (`ScoreTrackII.data[chrom]`).
#[derive(Debug, Clone, Default)]
pub struct ChromScore {
    /// Right-open interval end coordinates (C `int32`).
    pub pos: Vec<i32>,
    /// Treatment pileup per interval (C `float32`).
    pub treat: Vec<f32>,
    /// Control pileup per interval (C `float32`).
    pub ctrl: Vec<f32>,
    /// Computed score per interval (C `float32`).
    pub score: Vec<f32>,
}

/// Treatment/control score track (`ScoreTrack.py::ScoreTrackII`).
///
/// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
#[derive(Debug, Clone)]
pub struct ScoreTrack2 {
    /// Per-chromosome arrays, bytewise-sorted by chromosome name.
    data: BTreeMap<Vec<u8>, ChromScore>,
    /// Effective treatment depth (millions of filtered reads), `treat_edm`.
    pub treat_edm: f32,
    /// Effective control depth (millions of filtered reads), `ctrl_edm`.
    pub ctrl_edm: f32,
    /// Current scoring method byte (`scoring_method`; `b'N'` = unset).
    pub scoring_method: u8,
    /// Current normalization method byte (`normalization_method`; `b'N'`).
    pub normalization_method: u8,
    /// Pseudocount used by logLR / FE / logFE (`pseudocount`).
    pub pseudocount: f32,
    /// Whether a UCSC track line should be emitted (`trackline`).
    pub trackline: bool,
    /// Cutoff most recently used by `call_peaks` (`cutoff`).
    pub cutoff: f32,
    /// `-log10(pvalue)` -> total covered base pairs, for the p->q table
    /// (`pvalue_stat`). Keyed by the f32 bit pattern to match the float keys.
    pub pvalue_stat: BTreeMap<u32, i64>,
}

impl ScoreTrack2 {
    /// Create an empty track (`ScoreTrackII.__init__`).
    pub fn new(treat_depth: f32, ctrl_depth: f32, pseudocount: f32) -> Self {
        ScoreTrack2 {
            data: BTreeMap::new(),
            treat_edm: treat_depth,
            ctrl_edm: ctrl_depth,
            scoring_method: b'N',
            normalization_method: b'N',
            pseudocount,
            trackline: false,
            cutoff: 0.0,
            pvalue_stat: BTreeMap::new(),
        }
    }

    /// Update the pseudocount (`set_pseudocount`).
    pub fn set_pseudocount(&mut self, pseudocount: f32) {
        self.pseudocount = pseudocount;
    }

    /// Enable UCSC track line output (`enable_trackline`).
    pub fn enable_trackline(&mut self) {
        self.trackline = true;
    }

    /// Ensure storage exists for `chrom` (`add_chromosome`). Unlike the Cython
    /// version we grow the vectors lazily in `add`, so `chrom_max_len` is only
    /// used as a capacity hint.
    pub fn add_chromosome(&mut self, chrom: &[u8], chrom_max_len: usize) {
        self.data.entry(chrom.to_vec()).or_insert_with(|| ChromScore {
            pos: Vec::with_capacity(chrom_max_len),
            treat: Vec::with_capacity(chrom_max_len),
            ctrl: Vec::with_capacity(chrom_max_len),
            score: Vec::with_capacity(chrom_max_len),
        });
    }

    /// Append a `(endpos, chip, control)` interval to `chromosome` (`add`).
    /// The score slot is initialised to `0.0` and filled by the chosen scoring
    /// method later.
    pub fn add(&mut self, chromosome: &[u8], endpos: i32, chip: f32, control: f32) {
        let c = self
            .data
            .get_mut(chromosome)
            .expect("add_chromosome must precede add");
        c.pos.push(endpos);
        c.treat.push(chip);
        c.ctrl.push(control);
        c.score.push(0.0);
    }

    /// No-op trim hook (`finalize`). Our vectors are already exactly sized
    /// because `add` pushes; kept for parity with the MACS call sequence.
    pub fn finalize(&mut self) {}

    /// Immutable access to a chromosome's arrays (`get_data_by_chr`).
    pub fn get_data_by_chr(&self, chromosome: &[u8]) -> Option<&ChromScore> {
        self.data.get(chromosome)
    }

    /// Mutable access to a chromosome's arrays (used by Phase 1C scoring).
    pub fn get_data_by_chr_mut(&mut self, chromosome: &[u8]) -> Option<&mut ChromScore> {
        self.data.get_mut(chromosome)
    }

    /// All chromosome names, bytewise-sorted (`get_chr_names`).
    pub fn get_chr_names(&self) -> Vec<&[u8]> {
        self.data.keys().map(|k| k.as_slice()).collect()
    }

    /// Total number of intervals across all chromosomes (`total`).
    pub fn total(&self) -> i64 {
        self.data.values().map(|c| c.pos.len() as i64).sum()
    }

    /// Build a [`ScoreTrack2`] by merge-walking two bedGraph tracks, porting
    /// `BedGraph.py::bedGraphTrackI.make_ScoreTrackII_for_macs`.
    ///
    /// For each chromosome common to `treat` and `ctrl`, the two `(pos, val)`
    /// interval lists are merged on their boundaries: at each step the smaller
    /// of the two current end positions becomes the next interval end, carrying
    /// the *current* treatment value and the *current* control value. `depth1`
    /// and `depth2` seed `treat_edm` / `ctrl_edm` (bdgcmp passes
    /// `1/scaling_factor` for both). `finalize` is a no-op here (vectors are
    /// already exactly sized).
    pub fn from_two_bedgraphs(
        treat: &BedGraphTrack,
        ctrl: &BedGraphTrack,
        depth1: f32,
        depth2: f32,
    ) -> Self {
        let mut ret = ScoreTrack2::new(depth1, depth2, 1.0);

        // common chromosomes, bytewise-sorted (both accessors return sorted).
        let chr1: std::collections::BTreeSet<&[u8]> = treat.get_chr_names().into_iter().collect();
        let chr2: std::collections::BTreeSet<&[u8]> = ctrl.get_chr_names().into_iter().collect();
        for &chrom in chr1.intersection(&chr2) {
            let d1 = treat.get_data_by_chr(chrom).expect("chrom present in treat");
            let d2 = ctrl.get_data_by_chr(chrom).expect("chrom present in ctrl");
            let p1s = &d1.pos;
            let v1s = &d1.val;
            let p2s = &d2.pos;
            let v2s = &d2.val;

            let chrom_max_len = p1s.len() + p2s.len();
            ret.add_chromosome(chrom, chrom_max_len);

            // Two-pointer merge over interval boundaries; both lists are
            // non-empty in MACS (a track has at least one interval per chrom).
            if p1s.is_empty() || p2s.is_empty() {
                continue;
            }
            let mut i1 = 0usize;
            let mut i2 = 0usize;
            loop {
                let p1 = p1s[i1];
                let v1 = v1s[i1];
                let p2 = p2s[i2];
                let v2 = v2s[i2];
                if p1 < p2 {
                    ret.add(chrom, p1, v1, v2);
                    i1 += 1;
                    if i1 >= p1s.len() {
                        break;
                    }
                } else if p2 < p1 {
                    ret.add(chrom, p2, v1, v2);
                    i2 += 1;
                    if i2 >= p2s.len() {
                        break;
                    }
                } else {
                    // p1 == p2
                    ret.add(chrom, p1, v1, v2);
                    i1 += 1;
                    i2 += 1;
                    if i1 >= p1s.len() || i2 >= p2s.len() {
                        break;
                    }
                }
            }
        }
        ret.finalize();
        ret
    }

    // ------------------------------------------------------------------
    // Stubbed scoring / IO — Phase 1C / Phase 2 fill these in this file.
    // ------------------------------------------------------------------

    /// Dispatch a scoring method by its MACS byte (`change_score_method`):
    /// `p q l s f F d M m`.
    pub fn change_score_method(&mut self, scoring_method: u8) {
        match scoring_method {
            b'p' => self.compute_pvalue(),
            b'q' => {
                if self.scoring_method != b'p' {
                    self.compute_pvalue();
                }
                self.compute_qvalue();
            }
            b'l' => self.compute_likelihood(),
            b's' => self.compute_sym_likelihood(),
            b'f' => self.compute_logfe(),
            b'F' => self.compute_foldenrichment(),
            b'd' => self.compute_subtraction(),
            b'm' => self.compute_spmr(),
            b'M' => self.compute_max(),
            other => panic!("change_score_method: unsupported method byte {other}"),
        }
    }

    /// `compute_pvalue` (`ScoreTrack.py::compute_pvalue`): per-interval
    /// `-log10 p` Poisson score, accumulating `pvalue_stat[score] += length`.
    fn compute_pvalue(&mut self) {
        let pseudocount = self.pseudocount;
        // Split-borrow disjoint fields so the per-chrom loop and pvalue_stat
        // updates don't alias `self`.
        let data = &mut self.data;
        let stat = &mut self.pvalue_stat;
        for chrom in data.values_mut() {
            let mut prev_pos: i64 = 0;
            for i in 0..chrom.pos.len() {
                // observed = cast(int, p[i] + pseudocount); expectation =
                // c[i] + pseudocount (passed as f64 per the plan).
                let observed = (chrom.treat[i] + pseudocount) as i32;
                let expectation = (chrom.ctrl[i] + pseudocount) as f64;
                let v = get_pscore(observed, expectation);
                chrom.score[i] = v;
                let tmp_l = chrom.pos[i] as i64 - prev_pos;
                *stat.entry(v.to_bits()).or_insert(0) += tmp_l;
                prev_pos = chrom.pos[i] as i64;
            }
        }
        self.scoring_method = b'p';
    }

    /// `compute_qvalue` (`ScoreTrack.py::compute_qvalue`): convert each pscore to
    /// its qscore via the p->q table. Requires p-values computed first.
    fn compute_qvalue(&mut self) {
        debug_assert_eq!(self.scoring_method, b'p');
        let pqtable = self.make_pq_table();
        for chrom in self.data.values_mut() {
            for i in 0..chrom.score.len() {
                let key = chrom.score[i].to_bits();
                chrom.score[i] = *pqtable.get(&key).expect("pscore present in pq table");
            }
        }
        self.scoring_method = b'q';
    }

    /// `compute_likelihood` (`ScoreTrack.py::compute_likelihood`): asymmetric
    /// logLR of treatment vs control with the pseudocount added.
    fn compute_likelihood(&mut self) {
        let pseudocount = self.pseudocount;
        for chrom in self.data.values_mut() {
            for i in 0..chrom.score.len() {
                chrom.score[i] = log_lr_asym(chrom.treat[i] + pseudocount, chrom.ctrl[i] + pseudocount);
            }
        }
        self.scoring_method = b'l';
    }

    /// `compute_sym_likelihood` (`ScoreTrack.py::compute_sym_likelihood`).
    fn compute_sym_likelihood(&mut self) {
        let pseudocount = self.pseudocount;
        for chrom in self.data.values_mut() {
            for i in 0..chrom.score.len() {
                chrom.score[i] = log_lr_sym(chrom.treat[i] + pseudocount, chrom.ctrl[i] + pseudocount);
            }
        }
        self.scoring_method = b's';
    }

    /// `compute_logFE` (`ScoreTrack.py::compute_logFE`): `log10((p+pc)/(c+pc))`.
    /// MACS computes `log10(x/y)` where `x`,`y` are `cython.float`; the division
    /// is f32 and `log10` returns `double`, narrowed to the f32 score slot.
    fn compute_logfe(&mut self) {
        let pseudocount = self.pseudocount;
        for chrom in self.data.values_mut() {
            for i in 0..chrom.score.len() {
                let x = chrom.treat[i] + pseudocount;
                let y = chrom.ctrl[i] + pseudocount;
                chrom.score[i] = (x / y).log10();
            }
        }
        self.scoring_method = b'f';
    }

    /// `compute_foldenrichment` (`ScoreTrack.py::compute_foldenrichment`):
    /// linear `(p+pc)/(c+pc)`, all f32.
    fn compute_foldenrichment(&mut self) {
        let pseudocount = self.pseudocount;
        for chrom in self.data.values_mut() {
            for i in 0..chrom.score.len() {
                chrom.score[i] = (chrom.treat[i] + pseudocount) / (chrom.ctrl[i] + pseudocount);
            }
        }
        self.scoring_method = b'F';
    }

    /// `compute_subtraction` (`ScoreTrack.py::compute_subtraction`): `p - c`.
    fn compute_subtraction(&mut self) {
        for chrom in self.data.values_mut() {
            for i in 0..chrom.score.len() {
                chrom.score[i] = chrom.treat[i] - chrom.ctrl[i];
            }
        }
        self.scoring_method = b'd';
    }

    /// `compute_SPMR` (`ScoreTrack.py::compute_SPMR`): treatment pileup divided
    /// by the relevant depth scale (per the current normalization method).
    fn compute_spmr(&mut self) {
        let scale: f32 = match self.normalization_method {
            b'T' | b'N' => self.treat_edm,
            b'C' => self.ctrl_edm,
            b'M' => 1.0,
            _ => self.treat_edm,
        };
        for chrom in self.data.values_mut() {
            for i in 0..chrom.score.len() {
                chrom.score[i] = chrom.treat[i] / scale;
            }
        }
        self.scoring_method = b'm';
    }

    /// `compute_max` (`ScoreTrack.py::compute_max`): element-wise max of
    /// treatment and control pileups.
    fn compute_max(&mut self) {
        for chrom in self.data.values_mut() {
            for i in 0..chrom.score.len() {
                let p = chrom.treat[i];
                let c = chrom.ctrl[i];
                // Python's max(p, c) returns the first arg on ties (p == c).
                chrom.score[i] = if c > p { c } else { p };
            }
        }
        self.scoring_method = b'M';
    }

    /// Scale treatment and control pileups in place (`ScoreTrack.py::normalize`).
    fn normalize(&mut self, treat_scale: f32, control_scale: f32) {
        for chrom in self.data.values_mut() {
            for i in 0..chrom.pos.len() {
                chrom.treat[i] *= treat_scale;
                chrom.ctrl[i] *= control_scale;
            }
        }
    }

    /// Change normalization method (`change_normalization_method`). Dispatches
    /// the same transition table MACS uses; bdgcmp only ever calls this with
    /// `'M'` (to override SPMR by multiplying by the scaling factor).
    pub fn change_normalization_method(&mut self, normalization_method: u8) {
        let cur = self.normalization_method;
        match normalization_method {
            b'T' => {
                match cur {
                    b'T' => {}
                    b'C' => self.normalize(self.treat_edm / self.ctrl_edm, self.treat_edm / self.ctrl_edm),
                    b'M' => self.normalize(self.treat_edm, self.treat_edm),
                    b'N' => self.normalize(1.0, self.treat_edm / self.ctrl_edm),
                    _ => unimplemented!(),
                }
                self.normalization_method = b'T';
            }
            b'C' => {
                match cur {
                    b'T' => self.normalize(self.ctrl_edm / self.treat_edm, self.ctrl_edm / self.treat_edm),
                    b'C' => {}
                    b'M' => self.normalize(self.ctrl_edm, self.ctrl_edm),
                    b'N' => self.normalize(self.ctrl_edm / self.treat_edm, 1.0),
                    _ => unimplemented!(),
                }
                self.normalization_method = b'C';
            }
            b'M' => {
                match cur {
                    b'T' => self.normalize(1.0 / self.treat_edm, 1.0 / self.treat_edm),
                    b'C' => self.normalize(1.0 / self.ctrl_edm, 1.0 / self.ctrl_edm),
                    b'M' => {}
                    b'N' => self.normalize(1.0 / self.treat_edm, 1.0 / self.ctrl_edm),
                    _ => unimplemented!(),
                }
                self.normalization_method = b'M';
            }
            b'N' => {
                match cur {
                    b'T' => self.normalize(self.treat_edm, self.treat_edm),
                    b'C' => self.normalize(self.ctrl_edm, self.ctrl_edm),
                    b'M' => self.normalize(self.treat_edm, self.ctrl_edm),
                    b'N' => {}
                    _ => unimplemented!(),
                }
                self.normalization_method = b'N';
            }
            other => panic!("change_normalization_method: unsupported byte {other}"),
        }
    }

    /// Build the `-log10 p` -> `-log10 q` lookup table (`make_pq_table`),
    /// porting `ScoreTrack.py::make_pq_table` (the callpeak BH path).
    ///
    /// `N = Σ length`, `f = -log10(N)`, rank `k` accumulates bp length, and for
    /// each unique pscore (descending) `q = v + (log10(k) + f)`, clamped so it
    /// never increases (`if q > pre_q: q = pre_q`) and stops once `q <= 0`
    /// (remaining pscores all map to 0). Returns a map keyed by the f32 bit
    /// pattern of each pscore.
    pub fn make_pq_table(&self) -> BTreeMap<u32, f32> {
        debug_assert_eq!(self.scoring_method, b'p');

        // N = sum of all lengths.
        let n: i64 = self.pvalue_stat.values().sum();
        let mut k: f64 = 1.0; // rank, accumulated bp length
        let f: f64 = -(n as f64).log10();
        let mut pre_q: f64 = 2147483647.0; // previous q-value

        // unique pscores sorted descending. pvalue_stat is keyed by f32 bits,
        // so decode each key to an f32 and sort by value (descending).
        let mut unique: Vec<(f32, i64)> = self
            .pvalue_stat
            .iter()
            .map(|(&bits, &ln)| (f32::from_bits(bits), ln))
            .collect();
        unique.sort_by(|a, b| b.0.partial_cmp(&a.0).expect("pscores are finite"));

        let mut table: BTreeMap<u32, f32> = BTreeMap::new();
        let mut i = 0usize;
        while i < unique.len() {
            let (v, ln) = unique[i];
            let q: f64 = v as f64 + (k.log10() + f);
            let mut q = q;
            if q > pre_q {
                q = pre_q;
            }
            if q <= 0.0 {
                // q = 0; break — remaining entries handled below.
                break;
            }
            table.insert(v.to_bits(), q as f32);
            pre_q = q;
            k += ln as f64;
            i += 1;
        }
        // bottom-rank pscores all get qscore 0 (from index i onward).
        for &(v, _) in &unique[i..] {
            table.insert(v.to_bits(), 0.0);
        }
        table
    }

    /// Call narrow peaks above `cutoff` (`ScoreTrack.py::call_peaks`).
    ///
    /// Scans each chromosome for runs of intervals whose score is `>= cutoff`
    /// (note: `>=`, plan #7), merges runs separated by `<= max_gap`, and reports
    /// peaks at least `min_length` long. `call_summits = true` (signal-processed
    /// multi-summit calling) is Phase 2; only the single-summit `__close_peak`
    /// path is implemented here.
    pub fn call_peaks(
        &mut self,
        cutoff: f32,
        min_length: i32,
        max_gap: i32,
        call_summits: bool,
    ) -> PeakIO {
        assert!(
            !call_summits,
            "ScoreTrack2::call_peaks call_summits=true is Phase 2 (needs SignalProcessing)"
        );
        let mut peaks = PeakIO::default();
        self.cutoff = cutoff;

        let chroms: Vec<Vec<u8>> = self.data.keys().cloned().collect();
        for chrom in &chroms {
            // Collect indices where score >= cutoff, with their start/end/score
            // and sample-pileup, replicating the numpy gather in the source.
            let c = self.data.get(chrom.as_slice()).expect("chrom present");
            let ln = c.pos.len();
            // peak_content: (start, end, value, sample_value, index)
            let mut peak_content: Vec<(i32, i32, f32, f32, usize)> = Vec::new();
            let mut started = false;

            for idx in 0..ln {
                if c.score[idx] < cutoff {
                    continue;
                }
                let endpos = c.pos[idx];
                // start position = pos[idx-1]; for idx==0 the source clamps it
                // to 0 (above_cutoff_startpos[0] = 0 when above_cutoff[0]==0).
                let startpos = if idx == 0 { 0 } else { c.pos[idx - 1] };
                let value = c.score[idx];
                let sv = c.treat[idx];

                if !started {
                    peak_content.push((startpos, endpos, value, sv, idx));
                    started = true;
                } else if startpos - peak_content.last().unwrap().1 <= max_gap {
                    peak_content.push((startpos, endpos, value, sv, idx));
                } else {
                    self.close_peak(&peak_content, &mut peaks, min_length, chrom);
                    peak_content.clear();
                    peak_content.push((startpos, endpos, value, sv, idx));
                }
            }
            if !peak_content.is_empty() {
                self.close_peak(&peak_content, &mut peaks, min_length, chrom);
            }
        }
        peaks
    }

    /// Close one peak region (`ScoreTrack.py::__close_peak`). The summit is the
    /// midpoint of the highest-`sample_value` region; ties pick the middle of
    /// the tied midpoints (`midindex = (n+1)/2 - 1`). Adds to `peaks` only when
    /// the region is at least `min_length` long.
    fn close_peak(
        &self,
        peak_content: &[(i32, i32, f32, f32, usize)],
        peaks: &mut PeakIO,
        min_length: i32,
        chrom: &[u8],
    ) {
        let first = peak_content[0];
        let last = peak_content[peak_content.len() - 1];
        let peak_length = last.1 - first.0;
        if peak_length < min_length {
            return;
        }

        // Find the summit: middle of the highest-sample-value region(s).
        // MACS: `if not summit_value or summit_value < tsummitvalue:` — note
        // that a `summit_value` of exactly 0.0 is "falsy" and forces a reset, so
        // we replicate `summit_value == 0.0 || summit_value < tsummitvalue`.
        let mut tsummit: Vec<i32> = Vec::new();
        let mut tsummit_index: Vec<usize> = Vec::new();
        let mut summit_value: f32 = 0.0;
        for &(tstart, tend, _tvalue, tsummitvalue, tindex) in peak_content {
            if summit_value == 0.0 || summit_value < tsummitvalue {
                tsummit = vec![(tend + tstart) / 2];
                tsummit_index = vec![tindex];
                summit_value = tsummitvalue;
            } else if summit_value == tsummitvalue {
                tsummit.push((tend + tstart) / 2);
                tsummit_index.push(tindex);
            }
        }

        let midindex = ((tsummit.len() + 1) / 2) - 1;
        let summit_pos = tsummit[midindex];
        let summit_index = tsummit_index[midindex];

        let c = self.data.get(chrom).expect("chrom present");
        let qscore = if self.scoring_method == b'q' {
            c.score[summit_index]
        } else {
            -1.0
        };
        let pileup = c.treat[summit_index];
        let ctrl = c.ctrl[summit_index];
        let pscore = get_pscore(pileup as i32, ctrl as f64);
        let fold_change = (pileup + self.pseudocount) / (ctrl + self.pseudocount);

        peaks.add(
            chrom,
            first.0,
            last.1,
            summit_pos,
            c.score[summit_index],
            pileup,
            pscore,
            fold_change,
            qscore,
            b"",
        );
    }

    /// Call broad peaks (`call_broadpeaks`).
    pub fn call_broadpeaks(
        &mut self,
        _lvl1_cutoff: f32,
        _lvl2_cutoff: f32,
        _min_length: i32,
        _lvl1_max_gap: i32,
        _lvl2_max_gap: i32,
    ) -> BroadPeakIO {
        todo!("Phase 2: call_broadpeaks")
    }

    /// Cutoff-analysis report over a range of score thresholds
    /// (`cutoff_analysis`).
    pub fn cutoff_analysis(
        &self,
        _max_gap: i32,
        _min_length: i32,
        _steps: i32,
        _min_score: f32,
        _max_score: f32,
    ) -> String {
        todo!("Phase 2: cutoff_analysis")
    }

    /// Write column `column` (1=chip, 2=control, 3=score) as bedGraph using the
    /// `>1e-5` break predicate (`ScoreTrack.py::write_bedGraph`).
    ///
    /// Emits the optional UCSC track line when `self.trackline` is set, then for
    /// each chromosome (bytewise-sorted) walks the selected column with
    /// [`crate::bedgraph_io::Gt1e5`] coalescing via
    /// [`crate::bedgraph_io::write_bedgraph_predicate`].
    pub fn write_bedgraph<W: std::io::Write>(
        &self,
        fhd: &mut W,
        name: &str,
        description: &str,
        column: u8,
    ) -> std::io::Result<()> {
        assert!((1..=3).contains(&column), "column should be between 1, 2 or 3.");

        if self.trackline {
            // MACS writes: track type=bedGraph name="%s" description="%s"\n
            // with the raw (un-escaped) name/description here.
            write!(
                fhd,
                "track type=bedGraph name=\"{name}\" description=\"{description}\"\n"
            )?;
        }

        for (chrom, c) in &self.data {
            let value: &[f32] = match column {
                1 => &c.treat,
                2 => &c.ctrl,
                _ => &c.score,
            };
            crate::bedgraph_io::write_bedgraph_predicate(
                chrom,
                &c.pos,
                value,
                &crate::bedgraph_io::Gt1e5,
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
    fn add_and_accessors() {
        let mut s = ScoreTrack2::new(1.0, 1.0, 1.0);
        s.add_chromosome(b"chr1", 4);
        s.add(b"chr1", 100, 5.0, 2.0);
        s.add(b"chr1", 200, 7.0, 3.0);
        let c = s.get_data_by_chr(b"chr1").unwrap();
        assert_eq!(c.pos, vec![100, 200]);
        assert_eq!(c.treat, vec![5.0_f32, 7.0]);
        assert_eq!(c.ctrl, vec![2.0_f32, 3.0]);
        assert_eq!(c.score, vec![0.0_f32, 0.0]);
        assert_eq!(s.total(), 2);
    }

    #[test]
    fn chr_names_sorted() {
        let mut s = ScoreTrack2::new(1.0, 1.0, 1.0);
        s.add_chromosome(b"chr2", 1);
        s.add_chromosome(b"chr1", 1);
        assert_eq!(s.get_chr_names(), vec![&b"chr1"[..], &b"chr2"[..]]);
    }

    /// Test regions from `test_ScoreTrack.py::Test_ScoreTrackII.setUp`.
    /// `(endpos, chip, control)`.
    fn build_track() -> ScoreTrack2 {
        let regions = [(10, 100.0, 10.0), (60, 10.0, 10.0), (110, 15.0, 20.0), (160, 5.0, 20.0), (210, 20.0, 5.0)];
        let mut s = ScoreTrack2::new(10.0, 5.0, 1.0);
        s.add_chromosome(b"chrY", 5);
        for (endpos, chip, ctrl) in regions {
            s.add(b"chrY", endpos, chip, ctrl);
        }
        s
    }

    fn rounded2(v: &[f32]) -> Vec<f32> {
        v.iter().map(|x| (x * 100.0).round() / 100.0).collect()
    }

    /// Golden `compute_scores` (p/q/l/f/d/m) from `test_ScoreTrack.py`.
    #[test]
    fn golden_compute_scores() {
        let mut s = build_track();
        s.set_pseudocount(1.0);

        s.change_score_method(b'p');
        assert_eq!(
            rounded2(&s.get_data_by_chr(b"chrY").unwrap().score),
            vec![60.49, 0.38, 0.08, 0.0, 6.41]
        );

        s.change_score_method(b'q');
        assert_eq!(
            rounded2(&s.get_data_by_chr(b"chrY").unwrap().score),
            vec![58.17, 0.0, 0.0, 0.0, 5.13]
        );

        s.change_score_method(b'l');
        assert_eq!(
            rounded2(&s.get_data_by_chr(b"chrY").unwrap().score),
            vec![58.17, 0.0, -0.28, -3.25, 4.91]
        );

        s.change_score_method(b'f');
        assert_eq!(
            rounded2(&s.get_data_by_chr(b"chrY").unwrap().score),
            vec![0.96, 0.00, -0.12, -0.54, 0.54]
        );

        s.change_score_method(b'd');
        assert_eq!(
            rounded2(&s.get_data_by_chr(b"chrY").unwrap().score),
            vec![90.00, 0.0, -5.00, -15.00, 15.00]
        );

        s.change_score_method(b'm');
        assert_eq!(
            rounded2(&s.get_data_by_chr(b"chrY").unwrap().score),
            vec![10.00, 1.00, 1.50, 0.50, 2.00]
        );
    }

    /// Golden `normalize` (T/C/M/N) from `test_ScoreTrack.py`.
    #[test]
    fn golden_normalize() {
        // Build without pseudocount-1 scoring; treat_edm=10, ctrl_edm=5.
        let mut s = build_track();
        // norm_T: treat & ctrl scaled to treatment depth (10/5 etc.). The
        // golden arrays are [pos, treat, ctrl, score] transposed.
        s.change_normalization_method(b'T');
        let c = s.get_data_by_chr(b"chrY").unwrap();
        assert_eq!(c.treat, vec![100.0, 10.0, 15.0, 5.0, 20.0]);
        assert_eq!(c.ctrl, vec![20.0, 20.0, 40.0, 40.0, 10.0]);

        s.change_normalization_method(b'C');
        let c = s.get_data_by_chr(b"chrY").unwrap();
        assert_eq!(c.treat, vec![50.0, 5.0, 7.5, 2.5, 10.0]);
        assert_eq!(c.ctrl, vec![10.0, 10.0, 20.0, 20.0, 5.0]);

        s.change_normalization_method(b'M');
        let c = s.get_data_by_chr(b"chrY").unwrap();
        assert_eq!(c.treat, vec![10.0, 1.0, 1.5, 0.5, 2.0]);
        assert_eq!(c.ctrl, vec![2.0, 2.0, 4.0, 4.0, 1.0]);

        s.change_normalization_method(b'N');
        let c = s.get_data_by_chr(b"chrY").unwrap();
        assert_eq!(c.treat, vec![100.0, 10.0, 15.0, 5.0, 20.0]);
        assert_eq!(c.ctrl, vec![10.0, 10.0, 20.0, 20.0, 5.0]);
    }

    /// Golden `write_bedGraph` columns 1/2/3 from `test_ScoreTrack.py`.
    #[test]
    fn golden_write_bedgraph() {
        let mut s = build_track();
        s.change_score_method(b'p');

        let mut out1 = Vec::new();
        s.write_bedgraph(&mut out1, "NAME", "DESC", 1).unwrap();
        assert_eq!(
            String::from_utf8(out1).unwrap(),
            "chrY\t0\t10\t100.00000\n\
             chrY\t10\t60\t10.00000\n\
             chrY\t60\t110\t15.00000\n\
             chrY\t110\t160\t5.00000\n\
             chrY\t160\t210\t20.00000\n"
        );

        let mut out2 = Vec::new();
        s.write_bedgraph(&mut out2, "NAME", "DESC", 2).unwrap();
        assert_eq!(
            String::from_utf8(out2).unwrap(),
            "chrY\t0\t60\t10.00000\n\
             chrY\t60\t160\t20.00000\n\
             chrY\t160\t210\t5.00000\n"
        );

        let mut out3 = Vec::new();
        s.write_bedgraph(&mut out3, "NAME", "DESC", 3).unwrap();
        assert_eq!(
            String::from_utf8(out3).unwrap(),
            "chrY\t0\t10\t60.48912\n\
             chrY\t10\t60\t0.37599\n\
             chrY\t60\t110\t0.07723\n\
             chrY\t110\t160\t0.00006\n\
             chrY\t160\t210\t6.40804\n"
        );
    }

    /// Golden `call_peaks` (peak boundaries, summit, score) from
    /// `test_ScoreTrack.py::test_callpeak`.
    #[test]
    fn golden_call_peaks() {
        use crate::fmt::format_g6;
        let mut s = build_track();
        s.change_score_method(b'p');
        let p = s.call_peaks(0.10, 10, 10, false);

        let peaks = p.peaks_by_chr(b"chrY").expect("chrY peaks");
        assert_eq!(peaks.len(), 2);
        // peak1: chrY 0 60, summit 5, score %.6g = 60.4891
        assert_eq!((peaks[0].start, peaks[0].end, peaks[0].summit), (0, 60, 5));
        assert_eq!(format_g6(peaks[0].score), "60.4891");
        // peak2: chrY 160 210, summit 185, score %.6g = 6.40804
        assert_eq!((peaks[1].start, peaks[1].end, peaks[1].summit), (160, 210, 185));
        assert_eq!(format_g6(peaks[1].score), "6.40804");
    }

    /// Merge-walk constructor from two bedGraphs (`make_ScoreTrackII_for_macs`).
    #[test]
    fn from_two_bedgraphs_merge_walk() {
        use crate::bedgraph::BedGraphTrack;
        let mut t = BedGraphTrack::new(0.0);
        let mut c = BedGraphTrack::new(0.0);
        // treat: [0,10)=5, [10,30)=8 ; ctrl: [0,20)=2, [20,30)=4
        t.add_loc(b"chr1", 0, 10, 5.0);
        t.add_loc(b"chr1", 10, 30, 8.0);
        c.add_loc(b"chr1", 0, 20, 2.0);
        c.add_loc(b"chr1", 20, 30, 4.0);
        let s = ScoreTrack2::from_two_bedgraphs(&t, &c, 1.0, 1.0);
        let d = s.get_data_by_chr(b"chr1").unwrap();
        // boundaries: 10 (treat), 20 (ctrl), 30 (both)
        assert_eq!(d.pos, vec![10, 20, 30]);
        assert_eq!(d.treat, vec![5.0, 8.0, 8.0]);
        assert_eq!(d.ctrl, vec![2.0, 2.0, 4.0]);
    }
}
