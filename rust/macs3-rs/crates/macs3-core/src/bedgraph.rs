//! Sparse bedGraph signal track, ported from
//! `MACS3/Signal/BedGraph.py::bedGraphTrackI`.
//!
//! Data layout (per chromosome): two parallel vectors `pos: Vec<i32>` (right-
//! open transition end points) and `val: Vec<f32>`. Interval `i` spans
//! `[prev_end, pos[i])` where `prev_end` is `pos[i-1]` (and `0` for `i == 0`).
//! `baseline_value` fills uncovered bases when new segments are appended.
//!
//! Phase 1B owns the heavy algorithm methods (`apply_func`/`overlie`/`p2q`/
//! `call_peaks`/`call_broadpeaks`/...). The struct, `add_loc`, `add_loc_wo_merge`
//! and the accessors below are the **frozen** public API.

use std::collections::BTreeMap;

use rayon::prelude::*;

/// One chromosome's transition points and values (parallel vectors).
#[derive(Debug, Clone, Default)]
pub struct ChromData {
    /// Right-open end coordinates of each interval (0-based, ascending).
    pub pos: Vec<i32>,
    /// Value for each interval (C `float`).
    pub val: Vec<f32>,
}

/// Sparse bedGraph track, the `bedGraphTrackI` analogue.
///
/// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
#[derive(Debug, Clone)]
pub struct BedGraphTrack {
    // Chromosomes are stored in a BTreeMap keyed by raw bytes so iteration is
    // bytewise-sorted, matching MACS's `sorted(self.__data.keys())`.
    data: BTreeMap<Vec<u8>, ChromData>,
    /// Largest value seen (initialised tiny so the first `add_loc` updates it).
    pub maxvalue: f32,
    /// Smallest value seen (initialised large so the first `add_loc` updates it).
    pub minvalue: f32,
    /// Baseline value inserted for the gap before the first interval.
    pub baseline_value: f32,
}

impl BedGraphTrack {
    /// Create an empty track with the given baseline value
    /// (`bedGraphTrackI.__init__`). MACS seeds `maxvalue`/`minvalue` with
    /// `-10000000`/`10000000` so the first inserted value overrides them.
    ///
    /// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
    pub fn new(baseline_value: f32) -> Self {
        BedGraphTrack {
            data: BTreeMap::new(),
            maxvalue: -10_000_000.0,
            minvalue: 10_000_000.0,
            baseline_value,
        }
    }

    /// Append interval `[startpos, endpos)` with `value` on `chromosome`,
    /// merging with the preceding interval when their values are exactly equal
    /// (`bedGraphTrackI.add_loc`).
    ///
    /// Caller guarantees regions are added sorted and non-overlapping. Matches
    /// the source exactly: `endpos <= 0` returns; `startpos < 0` is clamped to
    /// 0; a brand-new chromosome whose first interval does not start at 0 gets
    /// a leading baseline block.
    ///
    /// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
    pub fn add_loc(&mut self, chromosome: &[u8], mut startpos: i32, endpos: i32, value: f32) {
        if endpos <= 0 {
            return;
        }
        if startpos < 0 {
            startpos = 0;
        }

        match self.data.get_mut(chromosome) {
            None => {
                let mut c = ChromData::default();
                if startpos != 0 {
                    c.pos.push(startpos);
                    c.val.push(self.baseline_value);
                }
                c.pos.push(endpos);
                c.val.push(value);
                self.data.insert(chromosome.to_vec(), c);
            }
            Some(c) => {
                let pre_v = *c.val.last().expect("chromosome data is never empty");
                if pre_v == value {
                    *c.pos.last_mut().unwrap() = endpos;
                } else {
                    c.pos.push(endpos);
                    c.val.push(value);
                }
            }
        }

        if value > self.maxvalue {
            self.maxvalue = value;
        }
        if value < self.minvalue {
            self.minvalue = value;
        }
    }

    /// Append interval `[startpos, endpos)` without merging identical
    /// neighbours (`bedGraphTrackI.add_loc_wo_merge`). Values below the
    /// baseline are raised to the baseline before insertion.
    ///
    /// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
    pub fn add_loc_wo_merge(
        &mut self,
        chromosome: &[u8],
        mut startpos: i32,
        endpos: i32,
        mut value: f32,
    ) {
        if endpos <= 0 {
            return;
        }
        if startpos < 0 {
            startpos = 0;
        }
        if value < self.baseline_value {
            value = self.baseline_value;
        }

        let c = self.data.entry(chromosome.to_vec()).or_insert_with(|| {
            let mut c = ChromData::default();
            if startpos != 0 {
                c.pos.push(startpos);
                c.val.push(self.baseline_value);
            }
            c
        });
        c.pos.push(endpos);
        c.val.push(value);

        if value > self.maxvalue {
            self.maxvalue = value;
        }
        if value < self.minvalue {
            self.minvalue = value;
        }
    }

    /// Replace a chromosome's data with explicit position/value vectors
    /// (`bedGraphTrackI.add_chrom_data`). The two vectors must have equal
    /// length.
    ///
    /// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
    pub fn add_chrom_data(&mut self, chromosome: &[u8], pos: Vec<i32>, val: Vec<f32>) {
        debug_assert_eq!(pos.len(), val.len());
        if let Some(&maxv) = val.iter().reduce(|a, b| if b > a { b } else { a }) {
            if maxv > self.maxvalue {
                self.maxvalue = maxv;
            }
        }
        if let Some(&minv) = val.iter().reduce(|a, b| if b < a { b } else { a }) {
            if minv < self.minvalue {
                self.minvalue = minv;
            }
        }
        self.data.insert(chromosome.to_vec(), ChromData { pos, val });
    }

    /// Return all chromosome names, **bytewise-sorted**
    /// (`bedGraphTrackI.get_chr_names`).
    ///
    /// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
    pub fn get_chr_names(&self) -> Vec<&[u8]> {
        self.data.keys().map(|k| k.as_slice()).collect()
    }

    /// Return the data for `chromosome`, or `None` if absent
    /// (`bedGraphTrackI.get_data_by_chr`).
    ///
    /// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
    pub fn get_data_by_chr(&self, chromosome: &[u8]) -> Option<&ChromData> {
        self.data.get(chromosome)
    }

    /// Mutable access to a chromosome's data (used by Phase 1B algorithms).
    ///
    /// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
    pub fn get_data_by_chr_mut(&mut self, chromosome: &[u8]) -> Option<&mut ChromData> {
        self.data.get_mut(chromosome)
    }

    /// Reset the baseline value, clamp sub-baseline regions and re-merge
    /// (`bedGraphTrackI.reset_baseline`).
    ///
    /// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
    pub fn reset_baseline(&mut self, baseline_value: f32) {
        self.baseline_value = baseline_value;
        // filter_score(cutoff=baseline_value) + merge_regions() — stubbed for
        // Phase 1B. read_bedgraph relies on baseline 0 (the default), where no
        // clamping is needed, so this is safe to defer.
        // todo!("Phase 1B: filter_score + merge_regions for non-zero baseline")
    }

    /// Total number of intervals across all chromosomes
    /// (`bedGraphTrackI.total`).
    ///
    /// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
    pub fn total(&self) -> usize {
        self.data.values().map(|c| c.pos.len()).sum()
    }

    /// Whether the track holds no chromosomes.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    // ------------------------------------------------------------------
    // Algorithm methods — Phase 1B.
    // ------------------------------------------------------------------

    /// Coalesce adjacent segments that share the same value
    /// (`bedGraphTrackI.merge_regions`). Two adjacent intervals with exactly
    /// equal values (`==` on `f32`) are merged into one.
    fn merge_regions(&mut self) {
        for c in self.data.values_mut() {
            if c.pos.is_empty() {
                continue;
            }
            let mut new_pos: Vec<i32> = vec![c.pos[0]];
            let mut new_val: Vec<f32> = vec![c.val[0]];
            let mut pre_v = new_val[0];
            for i in 1..c.pos.len() {
                let pos = c.pos[i];
                let value = c.val[i];
                if value == pre_v {
                    *new_pos.last_mut().unwrap() = pos;
                } else {
                    new_pos.push(pos);
                    new_val.push(value);
                    pre_v = value;
                }
            }
            c.pos = new_pos;
            c.val = new_val;
        }
    }

    /// Apply `func` to every value in the track (`bedGraphTrackI.apply_func`).
    ///
    /// MACS does **not** merge adjacent regions that become equal after the
    /// function is applied; it only rewrites values in place. `maxvalue` /
    /// `minvalue` are updated by applying `func` to the previous extremes
    /// (matching the source's `self.maxvalue = func(self.maxvalue)`).
    ///
    /// Each chromosome's value vector is rewritten in parallel; the extreme
    /// updates are serial and applied to the stored extremes (not a reduction).
    pub fn apply_func(&mut self, func: impl Fn(f32) -> f32 + Sync) {
        self.data.par_iter_mut().for_each(|(_, c)| {
            for v in c.val.iter_mut() {
                *v = func(*v);
            }
        });
        self.maxvalue = func(self.maxvalue);
        self.minvalue = func(self.minvalue);
    }

    /// Convert p-value scores to q-value scores (`bedGraphTrackI.p2q`).
    ///
    /// GOTCHA (plan #2, VERIFIED): in `BedGraph.py::p2q` the `# k += l` line is
    /// COMMENTED OUT, so the rank `k` stays `1` for every unique value. The
    /// q-score for unique value `v` is therefore
    ///   `q = max(0, min(pre_q, v + log10(1) + (-log10(N))))`
    ///       `= max(0, min(pre_q, v - log10(N)))`
    /// where `N = Σ interval_length` over the whole track, walking values in
    /// strictly DECREASING order with `pre_q` clamping the result monotonic
    /// (`pre_q` starts at `2147483647`). The lookup table maps each pvalue to
    /// its qvalue; afterwards adjacent equal regions are merged.
    pub fn p2q(&mut self) {
        // 1. frequencies of each p-score (total length per unique value).
        //    Iterate chromosomes in bytewise-sorted order like the source.
        let mut pvalue_stat: BTreeMap<u32, i64> = BTreeMap::new();
        let mut order: Vec<f32> = Vec::new(); // unique values, insertion order
        let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
        for c in self.data.values() {
            let mut pre_p: i64 = 0;
            for i in 0..c.pos.len() {
                let this_p = c.pos[i] as i64;
                let this_v = c.val[i];
                let this_l = this_p - pre_p;
                let key = this_v.to_bits();
                if seen.insert(key) {
                    order.push(this_v);
                }
                *pvalue_stat.entry(key).or_insert(0) += this_l;
                pre_p = this_p;
            }
        }

        let n: i64 = pvalue_stat.values().sum();
        let f: f64 = -(n as f64).log10();

        // 2. q-score for each p-score, walking unique values in DECREASING
        //    order with monotonic clamping. k stays 1 -> log10(k) = 0.
        let mut uniq = order.clone();
        uniq.sort_by(|a, b| b.partial_cmp(a).unwrap()); // descending
        let mut pqtable: std::collections::HashMap<u32, f32> = std::collections::HashMap::new();
        let mut pre_q: f32 = 2147483647.0;
        for &v in &uniq {
            // q = v + (log10(k) + f); log10(1) == 0
            let mut q = v as f64 + f;
            // q = max(0, min(pre_q, q))
            q = q.min(pre_q as f64).max(0.0);
            let q32 = q as f32;
            pqtable.insert(v.to_bits(), q32);
            pre_q = q32;
        }

        // 3. convert pscores to qscores in place.
        for c in self.data.values_mut() {
            for v in c.val.iter_mut() {
                *v = pqtable[&v.to_bits()];
            }
        }

        self.merge_regions();
    }

    /// Call narrow peaks where signal stays `>=` `cutoff`
    /// (`bedGraphTrackI.call_peaks`). GOTCHA (plan #7): the comparison is `>=`
    /// (not strict `>`). Segments above the cutoff are merged when separated by
    /// gaps no larger than `max_gap`; peaks shorter than `min_length` are
    /// dropped. The `call_summits` flag is accepted for API parity but unused
    /// (MACS always computes a single summit here).
    ///
    /// **Phase C (parallel):** each chromosome's segmentation runs independently
    /// and produces a `Vec<NarrowPeak>` in ascending position order. **Phase D
    /// (serial):** results are collected into `PeakIO` in bytewise-sorted
    /// chromosome order, preserving peak-number correctness.
    pub fn call_peaks(
        &self,
        cutoff: f32,
        min_length: i32,
        max_gap: i32,
        _call_summits: bool,
    ) -> crate::peak_io::PeakIO {
        let chroms = self.get_chr_names();

        // Phase C: parallel per-chrom segmentation.
        let per_chrom: Vec<Vec<crate::peak_io::NarrowPeak>> = chroms
            .par_iter()
            .map(|&chrom| self.chrom_call_peaks(chrom, cutoff, min_length, max_gap))
            .collect();

        // Phase D: serial collect into PeakIO in bytewise-sorted chromosome order.
        let mut peaks = crate::peak_io::PeakIO::default();
        for (chrom, chrom_peaks) in chroms.iter().zip(per_chrom.into_iter()) {
            for peak in chrom_peaks {
                peaks.add_peak_content(chrom, peak);
            }
        }
        peaks
    }

    /// Per-chromosome narrow-peak segmentation (Phase C helper).
    ///
    /// Returns peaks in ascending position order (insertion order == coordinate
    /// order here, matching the serial algorithm). No shared mutable state —
    /// safe to call in parallel across chromosomes.
    fn chrom_call_peaks(
        &self,
        chrom: &[u8],
        cutoff: f32,
        min_length: i32,
        max_gap: i32,
    ) -> Vec<crate::peak_io::NarrowPeak> {
        let c = self.data.get(chrom).unwrap();
        let ps = &c.pos;
        let vs = &c.val;
        let mut result: Vec<crate::peak_io::NarrowPeak> = Vec::new();
        // peak_content: list of (start, end, value)
        let mut peak_content: Vec<(i32, i32, f32)> = Vec::new();
        let mut pre_p: i32 = 0;
        // Find the first region above cutoff.
        let mut x: usize = 0;
        let mut started = false;
        while x < ps.len() {
            let p = ps[x];
            let v = vs[x];
            x += 1;
            if v >= cutoff {
                peak_content.push((pre_p, p, v));
                pre_p = p;
                started = true;
                break;
            } else {
                pre_p = p;
            }
        }
        if !started {
            return result;
        }
        // Continue scanning the remaining regions.
        for i in x..ps.len() {
            let p = ps[i];
            let v = vs[i];
            if v < cutoff {
                pre_p = p;
                continue;
            }
            // For points above cutoff: gap = pre_p - last segment end.
            if pre_p - peak_content.last().unwrap().1 <= max_gap {
                peak_content.push((pre_p, p, v));
            } else {
                Self::close_peak(&peak_content, &mut result, min_length, chrom);
                peak_content = vec![(pre_p, p, v)];
            }
            pre_p = p;
        }
        if !peak_content.is_empty() {
            Self::close_peak(&peak_content, &mut result, min_length, chrom);
        }
        result
    }

    /// Convert buffered segments into a peak entry if `peak_length >=
    /// min_length` (`bedGraphTrackI.__close_peak`).
    ///
    /// Summit = midpoint `(tstart+tend)/2` (C int truncation) of the
    /// maximum-value segment; on ties the chosen index is `(len+1)/2 - 1`
    /// (C int division). The peak score column is the summit value.
    ///
    /// Pushes into the caller's local `Vec<NarrowPeak>` so this function is
    /// free of shared-state side effects (Phase C safe).
    fn close_peak(
        peak_content: &[(i32, i32, f32)],
        peaks: &mut Vec<crate::peak_io::NarrowPeak>,
        min_length: i32,
        chrom: &[u8],
    ) {
        let peak_length = peak_content.last().unwrap().1 - peak_content[0].0;
        if peak_length < min_length {
            return;
        }
        let mut tsummit: Vec<i32> = Vec::new();
        let mut summit_value: f32 = 0.0;
        for &(tstart, tend, tvalue) in peak_content {
            // `if not summit_value or summit_value < tvalue`: the `not
            // summit_value` is true exactly when summit_value == 0.0 (Python
            // truthiness of a float).
            if summit_value == 0.0 || summit_value < tvalue {
                tsummit = vec![((tend + tstart) / 2)];
                summit_value = tvalue;
            } else if summit_value == tvalue {
                tsummit.push((tend + tstart) / 2);
            }
        }
        // tie index: (len+1)//2 - 1 (C int division).
        let idx = ((tsummit.len() as i32 + 1) / 2 - 1) as usize;
        let summit = tsummit[idx];
        peaks.push(crate::peak_io::NarrowPeak::new(
            chrom.to_vec(),
            peak_content[0].0,
            peak_content.last().unwrap().1,
            summit,
            summit_value,
            0.0,
            0.0,
            0.0,
            0.0,
            b"".to_vec(),
        ));
    }

    /// Call broad peaks from high- and low-stringency thresholds
    /// (`bedGraphTrackI.call_broadpeaks`).
    ///
    /// Runs `call_peaks` twice (level-1 with `lvl1_cutoff`/`lvl1_max_gap`,
    /// level-2 with `lvl2_cutoff`/`lvl2_max_gap`), then nests level-1 peaks
    /// inside the level-2 linking regions via `add_broadpeak`.
    pub fn call_broadpeaks(
        &self,
        lvl1_cutoff: f32,
        lvl2_cutoff: f32,
        min_length: i32,
        lvl1_max_gap: i32,
        lvl2_max_gap: i32,
    ) -> crate::peak_io::BroadPeakIO {
        assert!(
            lvl1_cutoff > lvl2_cutoff,
            "level 1 cutoff should be larger than level 2."
        );
        assert!(
            lvl1_max_gap < lvl2_max_gap,
            "level 2 maximum gap should be larger than level 1."
        );
        let lvl1_peaks = self.call_peaks(lvl1_cutoff, min_length, lvl1_max_gap, false);
        let lvl2_peaks = self.call_peaks(lvl2_cutoff, min_length, lvl2_max_gap, false);

        let mut broadpeaks = crate::peak_io::BroadPeakIO::new();
        // Iterate chromosomes from lvl1_peaks in bytewise-sorted order.
        for chrom in lvl1_peaks.get_chr_names() {
            let lvl1chrom: &[crate::peak_io::NarrowPeak] =
                lvl1_peaks.peaks_by_chr(chrom).map(|v| v.as_slice()).unwrap_or(&[]);
            let lvl2chrom: &[crate::peak_io::NarrowPeak] =
                lvl2_peaks.peaks_by_chr(chrom).map(|v| v.as_slice()).unwrap_or(&[]);

            let mut tmppeakset: Vec<&crate::peak_io::NarrowPeak> = Vec::new();
            let mut lvl1_iter = lvl1chrom.iter();
            // Mirror the Python try/except over the lvl1 iterator: when lvl1 is
            // exhausted (StopIteration) we flush the current lvl2 and the rest.
            let mut lvl1: Option<&crate::peak_io::NarrowPeak> = lvl1_iter.next();
            let mut exhausted = lvl1.is_none();
            let mut last_i: usize = 0;
            if !exhausted {
                for i in 0..lvl2chrom.len() {
                    last_i = i;
                    let lvl2 = &lvl2chrom[i];
                    loop {
                        let l1 = lvl1.unwrap();
                        if lvl2.start <= l1.start && l1.end <= lvl2.end {
                            tmppeakset.push(l1);
                            lvl1 = lvl1_iter.next();
                            if lvl1.is_none() {
                                exhausted = true;
                                break;
                            }
                        } else {
                            Self::add_broadpeak(&mut broadpeaks, chrom, lvl2, &tmppeakset);
                            tmppeakset.clear();
                            break;
                        }
                    }
                    if exhausted {
                        break;
                    }
                }
            }
            if exhausted && !lvl2chrom.is_empty() {
                // StopIteration branch: flush lvl2chrom[last_i] with current
                // tmppeakset, then the remaining lvl2 peaks with empty sets.
                let lvl2 = &lvl2chrom[last_i];
                Self::add_broadpeak(&mut broadpeaks, chrom, lvl2, &tmppeakset);
                tmppeakset.clear();
                for j in (last_i + 1)..lvl2chrom.len() {
                    Self::add_broadpeak(&mut broadpeaks, chrom, &lvl2chrom[j], &tmppeakset);
                }
            }
        }
        broadpeaks
    }

    /// Append a broad peak assembled from a core (level-1) peak set
    /// (`bedGraphTrackI.__add_broadpeak`).
    fn add_broadpeak(
        bpeaks: &mut crate::peak_io::BroadPeakIO,
        chrom: &[u8],
        lvl2peak: &crate::peak_io::NarrowPeak,
        lvl1peakset: &[&crate::peak_io::NarrowPeak],
    ) {
        let start = lvl2peak.start;
        let end = lvl2peak.end;

        if lvl1peakset.is_empty() {
            // No strong/lvl1 peaks inside: complement with 1bp start/end blocks.
            bpeaks.add(
                chrom,
                start,
                end,
                lvl2peak.score,
                format!("{start}").into_bytes(),
                format!("{end}").into_bytes(),
                2,
                b"1,1".to_vec(),
                format!("0,{}", end - start - 1).into_bytes(),
                lvl2peak.pileup,
                lvl2peak.pscore,
                lvl2peak.fc,
                lvl2peak.qscore,
                b"NA".to_vec(),
            );
            return;
        }

        let mut thick_start_i = lvl1peakset[0].start;
        let mut thick_end_i = lvl1peakset.last().unwrap().end;
        let mut block_num = lvl1peakset.len() as i32;
        let mut block_sizes: Vec<u8> = lvl1peakset
            .iter()
            .map(|x| x.length.to_string())
            .collect::<Vec<_>>()
            .join(",")
            .into_bytes();
        let mut block_starts: Vec<u8> = lvl1peakset
            .iter()
            .map(|x| (x.start - start).to_string())
            .collect::<Vec<_>>()
            .join(",")
            .into_bytes();

        if thick_start_i != start {
            // add 1bp left block
            thick_start_i = start;
            block_num += 1;
            let mut nb = b"1,".to_vec();
            nb.extend_from_slice(&block_sizes);
            block_sizes = nb;
            let mut ns = b"0,".to_vec();
            ns.extend_from_slice(&block_starts);
            block_starts = ns;
        }
        if thick_end_i != end {
            // add 1bp right block
            thick_end_i = end;
            block_num += 1;
            block_sizes.extend_from_slice(b",1");
            block_starts.extend_from_slice(b",");
            block_starts.extend_from_slice(format!("{}", end - start - 1).as_bytes());
        }

        bpeaks.add(
            chrom,
            start,
            end,
            lvl2peak.score,
            format!("{thick_start_i}").into_bytes(),
            format!("{thick_end_i}").into_bytes(),
            block_num,
            block_sizes,
            block_starts,
            lvl2peak.pileup,
            lvl2peak.pscore,
            lvl2peak.fc,
            lvl2peak.qscore,
            b"NA".to_vec(),
        );
    }

    /// Cutoff-analysis report over a range of score thresholds
    /// (`bedGraphTrackI.cutoff_analysis`). Not needed by the macs-rna paths
    /// exercised here; left unimplemented.
    pub fn cutoff_analysis(
        &self,
        _max_gap: i32,
        _min_length: i32,
        _steps: i32,
        _min_score: f32,
        _max_score: f32,
    ) -> String {
        todo!("Phase 1B: cutoff_analysis (unused by macs-rna paths)")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_loc_new_chrom_from_zero() {
        let mut t = BedGraphTrack::new(0.0);
        t.add_loc(b"chr1", 0, 100, 3.0);
        let c = t.get_data_by_chr(b"chr1").unwrap();
        assert_eq!(c.pos, vec![100]);
        assert_eq!(c.val, vec![3.0]);
        assert_eq!(t.maxvalue, 3.0);
        assert_eq!(t.minvalue, 3.0);
    }

    #[test]
    fn add_loc_new_chrom_with_leading_baseline() {
        let mut t = BedGraphTrack::new(0.0);
        t.add_loc(b"chr1", 50, 100, 3.0);
        let c = t.get_data_by_chr(b"chr1").unwrap();
        assert_eq!(c.pos, vec![50, 100]);
        assert_eq!(c.val, vec![0.0, 3.0]);
    }

    #[test]
    fn add_loc_merges_equal_values() {
        let mut t = BedGraphTrack::new(0.0);
        t.add_loc(b"chr1", 0, 100, 3.0);
        t.add_loc(b"chr1", 100, 200, 3.0); // same value -> extend
        t.add_loc(b"chr1", 200, 300, 5.0); // different -> new interval
        let c = t.get_data_by_chr(b"chr1").unwrap();
        assert_eq!(c.pos, vec![200, 300]);
        assert_eq!(c.val, vec![3.0, 5.0]);
    }

    #[test]
    fn add_loc_endpos_nonpositive_returns() {
        let mut t = BedGraphTrack::new(0.0);
        t.add_loc(b"chr1", 0, 0, 3.0);
        assert!(t.get_data_by_chr(b"chr1").is_none());
    }

    #[test]
    fn chr_names_bytewise_sorted() {
        let mut t = BedGraphTrack::new(0.0);
        t.add_loc(b"chr2", 0, 10, 1.0);
        t.add_loc(b"chr10", 0, 10, 1.0);
        t.add_loc(b"chr1", 0, 10, 1.0);
        // bytewise: "chr1" < "chr10" < "chr2"
        let names: Vec<&[u8]> = t.get_chr_names();
        assert_eq!(names, vec![&b"chr1"[..], &b"chr10"[..], &b"chr2"[..]]);
    }

    #[test]
    fn apply_func_multiply_updates_extremes() {
        let mut t = BedGraphTrack::new(0.0);
        t.add_loc(b"chr1", 0, 10, 2.0);
        t.add_loc(b"chr1", 10, 20, 4.0);
        t.apply_func(|x| x * 1.5);
        let c = t.get_data_by_chr(b"chr1").unwrap();
        assert_eq!(c.val, vec![3.0, 6.0]);
        // maxvalue was 4.0, minvalue 2.0 -> scaled.
        assert_eq!(t.maxvalue, 6.0);
        assert_eq!(t.minvalue, 3.0);
    }

    /// Ground-truthed against MACS3 v3.0.4: regions 5,3,1,5 with total length
    /// 50 -> q = max(0, min(pre_q, v - log10(50))) walked descending.
    #[test]
    fn p2q_matches_macs3() {
        let mut t = BedGraphTrack::new(0.0);
        for (s, e, v) in [(0, 10, 5.0), (10, 20, 3.0), (20, 30, 1.0), (30, 50, 5.0)] {
            t.add_loc(b"chr1", s, e, v);
        }
        t.p2q();
        let c = t.get_data_by_chr(b"chr1").unwrap();
        // MACS3: pos [10,20,30,50] val [3.30103, 1.30103, 0.0, 3.30103]
        assert_eq!(c.pos, vec![10, 20, 30, 50]);
        let want = [3.30103_f32, 1.30103, 0.0, 3.30103];
        for (g, w) in c.val.iter().zip(want.iter()) {
            assert!((g - w).abs() < 1e-4, "p2q got {g} want {w}");
        }
    }

    /// Ground-truthed against MACS3: cutoff 2, min_length 50, max_gap 30.
    /// add_loc merges equal neighbours, so the two 5.0 and two 6.0 blocks
    /// coalesce -> two peaks [100,300] (score 5, summit 200) and [400,700]
    /// (score 6, summit 550).
    #[test]
    fn call_peaks_matches_macs3() {
        let mut t = BedGraphTrack::new(0.0);
        for (s, e, v) in [
            (0, 100, 1.0),
            (100, 200, 5.0),
            (200, 300, 5.0),
            (300, 400, 1.0),
            (400, 500, 6.0),
            (500, 700, 6.0),
        ] {
            t.add_loc(b"chr1", s, e, v);
        }
        let peaks = t.call_peaks(2.0, 50, 30, false);
        let v = peaks.peaks_by_chr(b"chr1").unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!((v[0].start, v[0].end, v[0].summit), (100, 300, 200));
        assert_eq!(v[0].score, 5.0);
        assert_eq!((v[1].start, v[1].end, v[1].summit), (400, 700, 550));
        assert_eq!(v[1].score, 6.0);
    }

    #[test]
    fn call_peaks_uses_ge_not_gt() {
        // A block exactly equal to the cutoff must be included (>=).
        let mut t = BedGraphTrack::new(0.0);
        t.add_loc(b"chr1", 0, 300, 2.0); // exactly cutoff
        let peaks = t.call_peaks(2.0, 50, 30, false);
        assert_eq!(peaks.total, 1);
    }

    #[test]
    fn call_peaks_min_length_rejects_short() {
        let mut t = BedGraphTrack::new(0.0);
        t.add_loc(b"chr1", 0, 40, 5.0); // length 40 < 50
        let peaks = t.call_peaks(2.0, 50, 30, false);
        assert!(peaks.is_empty());
    }

    /// Ground-truthed against MACS3 gappedPeak output: a single broad region
    /// [100,1400] with two lvl1 cores [300,600] and [1000,1400], yielding
    /// blockNum 3 (one synthetic 1bp left block), sizes 1,300,400, starts
    /// 0,200,900, integer score int(10*3.0)=30.
    #[test]
    fn call_broadpeaks_matches_macs3() {
        let mut t = BedGraphTrack::new(0.0);
        for (s, e, v) in [
            (0, 100, 0.5),
            (100, 300, 1.5),
            (300, 600, 3.0),
            (600, 800, 1.2),
            (800, 1000, 0.2),
            (1000, 1400, 2.5),
        ] {
            t.add_loc(b"chr1", s, e, v);
        }
        let bp = t.call_broadpeaks(2.0, 1.0, 50, 30, 400);
        let v = bp.peaks_by_chr(b"chr1").unwrap();
        assert_eq!(v.len(), 1);
        let p = &v[0];
        assert_eq!((p.start, p.end), (100, 1400));
        assert_eq!(p.score, 3.0);
        assert_eq!(p.block_num, 3);
        assert_eq!(p.block_sizes, b"1,300,400".to_vec());
        assert_eq!(p.block_starts, b"0,200,900".to_vec());
        assert_eq!(p.thick_start, b"100".to_vec());
        assert_eq!(p.thick_end, b"1400".to_vec());
    }
}
