//! Peak record containers, ported from `MACS3/IO/PeakIO.py`.
//!
//! Mirrors `PeakContent` / `PeakIO` (narrow peaks) and
//! `BroadPeakContent` / `BroadPeakIO` (broad / gapped peaks) field-for-field.
//! The structs, constructors, accessors and `add` methods are the **frozen**
//! public API (Phase 0). The writers (`write_to_narrowPeak`,
//! `write_to_broadPeak`, `write_to_gappedPeak`, `write_to_xls`, `write_to_bed`)
//! are stubbed with `todo!()` and filled by Phase 1B.
//!
//! Storage: peaks are kept in a `BTreeMap<Vec<u8>, Vec<...>>` keyed by raw
//! chromosome bytes, so iteration is bytewise-sorted to match MACS's
//! `sorted(self.peaks.keys())`. Within a chromosome the per-peak order is the
//! insertion order (MACS appends to a Python list); callers that need
//! coordinate order call `sort()`.

use std::collections::BTreeMap;
use std::io::Write;

use crate::fmt::write_g6;

/// Alphabetical label for a zero-based subpeak index
/// (`PeakIO.py::subpeak_letters`): `a, b, ..., z, aa, ab, ...`.
fn subpeak_letters(i: i32) -> String {
    if i < 26 {
        ((97 + i) as u8 as char).to_string()
    } else {
        let mut s = subpeak_letters(i / 26);
        s.push((97 + (i % 26)) as u8 as char);
        s
    }
}

/// Interpolate a `%s`-style bytes prefix: replace the first `%s` in
/// `name_prefix` with `name` (Python `name_prefix % name`). If there is no
/// `%s`, the prefix is returned unchanged (mirroring the `try/except` fallback).
fn interp_prefix(name_prefix: &[u8], name: &[u8]) -> Vec<u8> {
    if let Some(pos) = find_subslice(name_prefix, b"%s") {
        let mut out = Vec::with_capacity(name_prefix.len() + name.len());
        out.extend_from_slice(&name_prefix[..pos]);
        out.extend_from_slice(name);
        out.extend_from_slice(&name_prefix[pos + 2..]);
        out
    } else {
        name_prefix.to_vec()
    }
}

/// Same `%s` interpolation used for the description field.
fn interp_desc(description: &[u8], name: &[u8]) -> Vec<u8> {
    interp_prefix(description, name)
}

/// Find the first occurrence of `needle` in `haystack`.
fn find_subslice(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || needle.len() > haystack.len() {
        return None;
    }
    haystack
        .windows(needle.len())
        .position(|w| w == needle)
}

/// `int(10*score)` truncated toward zero, matching Python `int(10*x)` on a C
/// `float` widened to `double`.
fn int10(score: f32) -> i64 {
    (10.0_f64 * score as f64).trunc() as i64
}

/// A single narrow peak and its summit statistics
/// (`PeakIO.py::PeakContent`). All score fields are C `float` (`f32`).
#[derive(Debug, Clone, PartialEq)]
pub struct NarrowPeak {
    /// Chromosome name (raw bytes).
    pub chrom: Vec<u8>,
    /// 0-based inclusive start coordinate.
    pub start: i32,
    /// 0-based exclusive end coordinate.
    pub end: i32,
    /// Peak length in bp (`end - start`).
    pub length: i32,
    /// 0-based summit position.
    pub summit: i32,
    /// Peak score reported by MACS3.
    pub score: f32,
    /// Tag pileup at the summit.
    pub pileup: f32,
    /// `-log10(pvalue)` at the summit.
    pub pscore: f32,
    /// Fold enrichment at the summit.
    pub fc: f32,
    /// `-log10(qvalue)` at the summit.
    pub qscore: f32,
    /// Optional peak identifier.
    pub name: Vec<u8>,
}

impl NarrowPeak {
    /// Construct a peak record, deriving `length = end - start`
    /// (`PeakContent.__init__`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        chrom: Vec<u8>,
        start: i32,
        end: i32,
        summit: i32,
        peak_score: f32,
        pileup: f32,
        pscore: f32,
        fold_change: f32,
        qscore: f32,
        name: Vec<u8>,
    ) -> Self {
        NarrowPeak {
            chrom,
            start,
            end,
            length: end - start,
            summit,
            score: peak_score,
            pileup,
            pscore,
            fc: fold_change,
            qscore,
            name,
        }
    }
}

/// In-memory collection of narrow peaks (`PeakIO.py::PeakIO`).
///
/// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
#[derive(Debug, Clone)]
pub struct PeakIO {
    /// Per-chromosome peak lists, bytewise-sorted by chromosome name.
    peaks: BTreeMap<Vec<u8>, Vec<NarrowPeak>>,
    /// Whether peaks have been coordinate-sorted (`CO_sorted`).
    pub co_sorted: bool,
    /// Total number of peaks across all chromosomes.
    pub total: i64,
    /// Collection name used in output (`name`).
    pub name: Vec<u8>,
}

impl Default for PeakIO {
    fn default() -> Self {
        PeakIO::new(b"MACS3".to_vec())
    }
}

impl PeakIO {
    /// Create an empty collection (`PeakIO.__init__`, default name `b"MACS3"`).
    pub fn new(name: Vec<u8>) -> Self {
        PeakIO {
            peaks: BTreeMap::new(),
            co_sorted: false,
            total: 0,
            name,
        }
    }

    /// Append a peak described by raw coordinates and scores (`PeakIO.add`).
    #[allow(clippy::too_many_arguments)]
    pub fn add(
        &mut self,
        chromosome: &[u8],
        start: i32,
        end: i32,
        summit: i32,
        peak_score: f32,
        pileup: f32,
        pscore: f32,
        fold_change: f32,
        qscore: f32,
        name: &[u8],
    ) {
        self.peaks
            .entry(chromosome.to_vec())
            .or_default()
            .push(NarrowPeak::new(
                chromosome.to_vec(),
                start,
                end,
                summit,
                peak_score,
                pileup,
                pscore,
                fold_change,
                qscore,
                name.to_vec(),
            ));
        self.total += 1;
        self.co_sorted = false;
    }

    /// Append an existing `NarrowPeak` (`PeakIO.add_PeakContent`).
    pub fn add_peak_content(&mut self, chromosome: &[u8], peak: NarrowPeak) {
        self.peaks.entry(chromosome.to_vec()).or_default().push(peak);
        self.total += 1;
        self.co_sorted = false;
    }

    /// Return (creating if absent) the peak list for `chrom`
    /// (`PeakIO.get_data_from_chrom`).
    pub fn get_data_from_chrom(&mut self, chrom: &[u8]) -> &mut Vec<NarrowPeak> {
        self.peaks.entry(chrom.to_vec()).or_default()
    }

    /// Immutable access to a chromosome's peaks.
    pub fn peaks_by_chr(&self, chrom: &[u8]) -> Option<&Vec<NarrowPeak>> {
        self.peaks.get(chrom)
    }

    /// All chromosome names, **bytewise-sorted** (`PeakIO.get_chr_names`;
    /// the source returns a set but iterates it sorted everywhere).
    pub fn get_chr_names(&self) -> Vec<&[u8]> {
        self.peaks.keys().map(|k| k.as_slice()).collect()
    }

    /// Sort peaks on each chromosome by ascending start position
    /// (`PeakIO.sort`).
    pub fn sort(&mut self) {
        if self.co_sorted {
            return;
        }
        for v in self.peaks.values_mut() {
            v.sort_by_key(|p| p.start);
        }
        self.co_sorted = true;
    }

    /// Iterate chromosome / peak-list pairs in bytewise-sorted order.
    pub fn iter(&self) -> impl Iterator<Item = (&[u8], &Vec<NarrowPeak>)> {
        self.peaks.iter().map(|(k, v)| (k.as_slice(), v))
    }

    /// Whether the collection holds no peaks.
    pub fn is_empty(&self) -> bool {
        self.total == 0
    }

    // ------------------------------------------------------------------
    // Writers — Phase 1B.
    // ------------------------------------------------------------------

    /// Select a peak's value for `score_column` (`PeakContent.__getitem__`).
    fn col(peak: &NarrowPeak, score_column: &str) -> f32 {
        match score_column {
            "score" => peak.score,
            "pileup" => peak.pileup,
            "pscore" => peak.pscore,
            "fc" => peak.fc,
            "qscore" => peak.qscore,
            other => panic!("unsupported score_column {other:?}"),
        }
    }

    /// Write ENCODE narrowPeak (BED6+4) (`PeakIO.write_to_narrowPeak`).
    ///
    /// Columns: chrom, start, end, name, `int(10*score_column)`, `.`, fc
    /// (`%.6g`), pscore (`%.6g`), qscore (`%.6g`), summit offset
    /// (`summit - start`, or `-1` when summit is `-1`). Peaks sharing the same
    /// `end` (call-summits) are grouped and lettered.
    pub fn write_to_narrow_peak<W: Write>(
        &self,
        fhd: &mut W,
        name_prefix: &[u8],
        name: &[u8],
        score_column: &str,
        trackline: bool,
    ) -> std::io::Result<()> {
        let peakprefix = interp_prefix(name_prefix, name);
        if trackline {
            let nm = String::from_utf8_lossy(name);
            write!(
                fhd,
                "track type=narrowPeak name=\"{nm}\" description=\"{nm}\" nextItemButton=on\n"
            )?;
        }
        let mut n_peak: i32 = 0;
        let mut g6 = String::with_capacity(16);
        for (chrom, peaks) in self.iter() {
            for group in group_by_end(peaks) {
                n_peak += 1;
                let multi = group.len() > 1;
                for (i, peak) in group.iter().enumerate() {
                    let peakname = if multi {
                        format!(
                            "{}{}{}",
                            String::from_utf8_lossy(&peakprefix),
                            n_peak,
                            subpeak_letters(i as i32)
                        )
                    } else {
                        format!("{}{}", String::from_utf8_lossy(&peakprefix), n_peak)
                    };
                    let s = if peak.summit == -1 {
                        -1
                    } else {
                        peak.summit - peak.start
                    };
                    fhd.write_all(chrom)?;
                    write!(
                        fhd,
                        "\t{}\t{}\t{}\t{}\t.\t",
                        peak.start,
                        peak.end,
                        peakname,
                        int10(Self::col(peak, score_column)),
                    )?;
                    g6.clear();
                    write_g6(&mut g6, peak.fc);
                    fhd.write_all(g6.as_bytes())?;
                    fhd.write_all(b"\t")?;
                    g6.clear();
                    write_g6(&mut g6, peak.pscore);
                    fhd.write_all(g6.as_bytes())?;
                    fhd.write_all(b"\t")?;
                    g6.clear();
                    write_g6(&mut g6, peak.qscore);
                    fhd.write_all(g6.as_bytes())?;
                    write!(fhd, "\t{s}\n")?;
                }
            }
        }
        Ok(())
    }

    /// Write the MACS3 `.xls` report (`PeakIO.write_to_xls`).
    ///
    /// Coordinates are 1-based (`start+1`, `summit+1`); `pileup` is rounded to
    /// 2 decimals before `%.6g`; score columns use `%.6g`.
    pub fn write_to_xls<W: Write>(
        &self,
        ofhd: &mut W,
        name_prefix: &[u8],
        name: &[u8],
    ) -> std::io::Result<()> {
        writeln!(
            ofhd,
            "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"
        )?;
        let peakprefix = interp_prefix(name_prefix, name);
        let mut n_peak: i32 = 0;
        let mut g6 = String::with_capacity(16);
        for (chrom, peaks) in self.iter() {
            for group in group_by_end(peaks) {
                n_peak += 1;
                let multi = group.len() > 1;
                for (i, peak) in group.iter().enumerate() {
                    let peakname = if multi {
                        format!(
                            "{}{}{}",
                            String::from_utf8_lossy(&peakprefix),
                            n_peak,
                            subpeak_letters(i as i32)
                        )
                    } else {
                        format!("{}{}", String::from_utf8_lossy(&peakprefix), n_peak)
                    };
                    ofhd.write_all(chrom)?;
                    write!(
                        ofhd,
                        "\t{}\t{}\t{}\t{}\t",
                        peak.start + 1,
                        peak.end,
                        peak.length,
                        peak.summit + 1,
                    )?;
                    g6.clear();
                    write_g6(&mut g6, round2(peak.pileup));
                    ofhd.write_all(g6.as_bytes())?;
                    ofhd.write_all(b"\t")?;
                    g6.clear();
                    write_g6(&mut g6, peak.pscore);
                    ofhd.write_all(g6.as_bytes())?;
                    ofhd.write_all(b"\t")?;
                    g6.clear();
                    write_g6(&mut g6, peak.fc);
                    ofhd.write_all(g6.as_bytes())?;
                    ofhd.write_all(b"\t")?;
                    g6.clear();
                    write_g6(&mut g6, peak.qscore);
                    ofhd.write_all(g6.as_bytes())?;
                    write!(ofhd, "\t{peakname}\n")?;
                }
            }
        }
        Ok(())
    }

    /// Write BED5 peaks (`PeakIO._to_bed`): chrom, start, end, name,
    /// `score_column` (`%.6g`). Optional UCSC trackline.
    pub fn write_to_bed<W: Write>(
        &self,
        fhd: &mut W,
        name_prefix: &[u8],
        name: &[u8],
        description: &[u8],
        score_column: &str,
        trackline: bool,
    ) -> std::io::Result<()> {
        let peakprefix = interp_prefix(name_prefix, name);
        let desc = interp_desc(description, name);
        if trackline {
            let nm = escape_quotes(name);
            let ds = escape_quotes(&desc);
            write!(
                fhd,
                "track name=\"{} (peaks)\" description=\"{}\" visibility=1\n",
                String::from_utf8_lossy(&nm),
                String::from_utf8_lossy(&ds),
            )?;
        }
        let mut n_peak: i32 = 0;
        let mut g6 = String::with_capacity(16);
        for (chrom, peaks) in self.iter() {
            for group in group_by_end(peaks) {
                n_peak += 1;
                let multi = group.len() > 1;
                for (i, peak) in group.iter().enumerate() {
                    fhd.write_all(chrom)?;
                    if multi {
                        write!(
                            fhd,
                            "\t{}\t{}\t{}{}{}\t",
                            peak.start,
                            peak.end,
                            String::from_utf8_lossy(&peakprefix),
                            n_peak,
                            subpeak_letters(i as i32),
                        )?;
                    } else {
                        write!(
                            fhd,
                            "\t{}\t{}\t{}{}\t",
                            peak.start,
                            peak.end,
                            String::from_utf8_lossy(&peakprefix),
                            n_peak,
                        )?;
                    }
                    g6.clear();
                    write_g6(&mut g6, Self::col(peak, score_column));
                    fhd.write_all(g6.as_bytes())?;
                    fhd.write_all(b"\n")?;
                }
            }
        }
        Ok(())
    }

    /// Write BED5 summits (`PeakIO._to_summits_bed`): each summit is the 1bp
    /// interval `[summit, summit+1)` with `score_column` (`%.6g`).
    pub fn write_to_summit_bed<W: Write>(
        &self,
        fhd: &mut W,
        name_prefix: &[u8],
        name: &[u8],
        description: &[u8],
        score_column: &str,
        trackline: bool,
    ) -> std::io::Result<()> {
        let peakprefix = interp_prefix(name_prefix, name);
        let desc = interp_desc(description, name);
        if trackline {
            let nm = escape_quotes(name);
            let ds = escape_quotes(&desc);
            write!(
                fhd,
                "track name=\"{} (summits)\" description=\"{}\" visibility=1\n",
                String::from_utf8_lossy(&nm),
                String::from_utf8_lossy(&ds),
            )?;
        }
        let mut n_peak: i32 = 0;
        let mut g6 = String::with_capacity(16);
        for (chrom, peaks) in self.iter() {
            for group in group_by_end(peaks) {
                n_peak += 1;
                let multi = group.len() > 1;
                for (i, peak) in group.iter().enumerate() {
                    let summit_p = peak.summit;
                    fhd.write_all(chrom)?;
                    if multi {
                        write!(
                            fhd,
                            "\t{}\t{}\t{}{}{}\t",
                            summit_p,
                            summit_p + 1,
                            String::from_utf8_lossy(&peakprefix),
                            n_peak,
                            subpeak_letters(i as i32),
                        )?;
                    } else {
                        write!(
                            fhd,
                            "\t{}\t{}\t{}{}\t",
                            summit_p,
                            summit_p + 1,
                            String::from_utf8_lossy(&peakprefix),
                            n_peak,
                        )?;
                    }
                    g6.clear();
                    write_g6(&mut g6, Self::col(peak, score_column));
                    fhd.write_all(g6.as_bytes())?;
                    fhd.write_all(b"\n")?;
                }
            }
        }
        Ok(())
    }

    /// Filter peaks by minimum `-log10(qvalue)` (`PeakIO.filter_qscore`).
    pub fn filter_qscore(&mut self, qscore_cut: f64) {
        self.total = 0;
        for v in self.peaks.values_mut() {
            v.retain(|p| p.qscore as f64 >= qscore_cut);
            self.total += v.len() as i64;
        }
        self.co_sorted = true;
        self.sort();
    }

    /// Filter peaks by fold-change range (`PeakIO.filter_fc`): keep `fc >=
    /// fc_low`, and additionally `fc < fc_up` when `fc_up > 0 && fc_up >
    /// fc_low`.
    pub fn filter_fc(&mut self, fc_low: f32, fc_up: f32) {
        self.total = 0;
        let upper = fc_up > 0.0 && fc_up > fc_low;
        for v in self.peaks.values_mut() {
            if upper {
                v.retain(|p| p.fc >= fc_low && p.fc < fc_up);
            } else {
                v.retain(|p| p.fc >= fc_low);
            }
            self.total += v.len() as i64;
        }
        self.co_sorted = true;
        self.sort();
    }
}

/// Group a chromosome's peaks into runs of consecutive entries that share the
/// same `end` coordinate, mirroring `itertools.groupby(peaks, key=end)`.
/// Without `--call-summits` every group has length 1.
fn group_by_end(peaks: &[NarrowPeak]) -> Vec<&[NarrowPeak]> {
    let mut groups: Vec<&[NarrowPeak]> = Vec::new();
    let mut i = 0;
    while i < peaks.len() {
        let end = peaks[i].end;
        let mut j = i + 1;
        while j < peaks.len() && peaks[j].end == end {
            j += 1;
        }
        groups.push(&peaks[i..j]);
        i = j;
    }
    groups
}

/// `round(x, 2)` — Python banker's rounding to 2 decimals on the f64 widening.
fn round2(x: f32) -> f32 {
    let d = x as f64;
    let scaled = d * 100.0;
    let r = scaled.round_ties_even();
    (r / 100.0) as f32
}

/// Backslash-escape `"` in a byte string (Python `.replace(b"\"", b"\\\"")`).
fn escape_quotes(s: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(s.len());
    for &b in s {
        if b == b'"' {
            out.push(b'\\');
        }
        out.push(b);
    }
    out
}

/// A single broad peak with block structure (`PeakIO.py::BroadPeakContent`).
#[derive(Debug, Clone, PartialEq)]
pub struct BroadPeak {
    /// 0-based inclusive start coordinate.
    pub start: i32,
    /// 0-based exclusive end coordinate.
    pub end: i32,
    /// Peak length in bp (`end - start`).
    pub length: i32,
    /// Average score across blocks.
    pub score: f32,
    /// Start of the high-enrichment segment, or `b"."` (kept as raw bytes
    /// because MACS stores it as `b"%d"`-formatted bytes / `b"."`).
    pub thick_start: Vec<u8>,
    /// End of the high-enrichment segment, or `b"."`.
    pub thick_end: Vec<u8>,
    /// Number of sub-blocks composing the region.
    pub block_num: i32,
    /// Comma-separated block sizes (raw bytes).
    pub block_sizes: Vec<u8>,
    /// Comma-separated block starts (raw bytes).
    pub block_starts: Vec<u8>,
    /// Median pileup within the region.
    pub pileup: f32,
    /// Median `-log10(pvalue)`.
    pub pscore: f32,
    /// Median fold change.
    pub fc: f32,
    /// Median `-log10(qvalue)`.
    pub qscore: f32,
    /// Optional region identifier.
    pub name: Vec<u8>,
}

impl BroadPeak {
    /// Construct a broad-peak record (`BroadPeakContent.__init__`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        start: i32,
        end: i32,
        score: f32,
        thick_start: Vec<u8>,
        thick_end: Vec<u8>,
        block_num: i32,
        block_sizes: Vec<u8>,
        block_starts: Vec<u8>,
        pileup: f32,
        pscore: f32,
        fold_change: f32,
        qscore: f32,
        name: Vec<u8>,
    ) -> Self {
        BroadPeak {
            start,
            end,
            length: end - start,
            score,
            thick_start,
            thick_end,
            block_num,
            block_sizes,
            block_starts,
            pileup,
            pscore,
            fc: fold_change,
            qscore,
            name,
        }
    }
}

/// Gapped peak alias.
///
/// MACS3 has no distinct `GappedPeakContent`; the gappedPeak format is written
/// from a [`BroadPeak`] (see `BroadPeakIO.write_to_gappedPeak`). This alias
/// exists so downstream code can name the concept explicitly while the storage
/// stays identical to `BroadPeak`.
pub type GappedPeak = BroadPeak;

/// In-memory collection of broad peaks (`PeakIO.py::BroadPeakIO`).
///
/// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
#[derive(Debug, Clone, Default)]
pub struct BroadPeakIO {
    /// Per-chromosome broad-peak lists, bytewise-sorted by chromosome name.
    peaks: BTreeMap<Vec<u8>, Vec<BroadPeak>>,
}

impl BroadPeakIO {
    /// Create an empty collection (`BroadPeakIO.__init__`).
    pub fn new() -> Self {
        BroadPeakIO {
            peaks: BTreeMap::new(),
        }
    }

    /// Append a broad-peak record (`BroadPeakIO.add`).
    #[allow(clippy::too_many_arguments)]
    pub fn add(
        &mut self,
        chromosome: &[u8],
        start: i32,
        end: i32,
        score: f32,
        thick_start: Vec<u8>,
        thick_end: Vec<u8>,
        block_num: i32,
        block_sizes: Vec<u8>,
        block_starts: Vec<u8>,
        pileup: f32,
        pscore: f32,
        fold_change: f32,
        qscore: f32,
        name: Vec<u8>,
    ) {
        self.peaks
            .entry(chromosome.to_vec())
            .or_default()
            .push(BroadPeak::new(
                start,
                end,
                score,
                thick_start,
                thick_end,
                block_num,
                block_sizes,
                block_starts,
                pileup,
                pscore,
                fold_change,
                qscore,
                name,
            ));
    }

    /// Immutable access to a chromosome's broad peaks.
    pub fn peaks_by_chr(&self, chrom: &[u8]) -> Option<&Vec<BroadPeak>> {
        self.peaks.get(chrom)
    }

    /// All chromosome names, bytewise-sorted.
    pub fn get_chr_names(&self) -> Vec<&[u8]> {
        self.peaks.keys().map(|k| k.as_slice()).collect()
    }

    /// Iterate chromosome / broad-peak-list pairs in bytewise-sorted order.
    pub fn iter(&self) -> impl Iterator<Item = (&[u8], &Vec<BroadPeak>)> {
        self.peaks.iter().map(|(k, v)| (k.as_slice(), v))
    }

    /// Total number of broad peaks (`BroadPeakIO.total`).
    pub fn total(&self) -> i64 {
        self.peaks.values().map(|v| v.len() as i64).sum()
    }

    /// Whether the collection holds no peaks.
    pub fn is_empty(&self) -> bool {
        self.peaks.values().all(|v| v.is_empty())
    }

    // ------------------------------------------------------------------
    // Writers — Phase 1B.
    // ------------------------------------------------------------------

    /// Select a broad peak's value for `score_column`
    /// (`BroadPeakContent.__getitem__`).
    fn col(peak: &BroadPeak, score_column: &str) -> f32 {
        match score_column {
            "score" => peak.score,
            "pileup" => peak.pileup,
            "pscore" => peak.pscore,
            "fc" => peak.fc,
            "qscore" => peak.qscore,
            other => panic!("unsupported score_column {other:?}"),
        }
    }

    /// Write gappedPeak (BED12+3) (`BroadPeakIO.write_to_gappedPeak`).
    ///
    /// Layout: chrom, start, end, `name_prefix+n`, `int(10*score_column)`, `.`,
    /// `0`, `0`, `0`, blockNum, blockSizes, blockStarts, fc, pscore, qscore
    /// (the last three `%.6g`). Peaks with `thickStart == b"."` are skipped
    /// (they would violate gappedPeak); a plain per-peak loop, no groupby.
    pub fn write_to_gapped_peak<W: Write>(
        &self,
        fhd: &mut W,
        name_prefix: &[u8],
        name: &[u8],
        description: &[u8],
        score_column: &str,
        trackline: bool,
    ) -> std::io::Result<()> {
        let peakprefix = interp_prefix(name_prefix, name);
        let desc = interp_desc(description, name);
        if trackline {
            write!(
                fhd,
                "track name=\"{}\" description=\"{}\" type=gappedPeak nextItemButton=on\n",
                String::from_utf8_lossy(name),
                String::from_utf8_lossy(&desc),
            )?;
        }
        let mut n_peak: i32 = 0;
        let mut g6 = String::with_capacity(16);
        for (chrom, peaks) in self.iter() {
            for peak in peaks {
                n_peak += 1;
                if peak.thick_start == b"." {
                    continue;
                }
                fhd.write_all(chrom)?;
                write!(
                    fhd,
                    "\t{}\t{}\t{}{}\t{}\t.\t0\t0\t0\t{}\t{}\t{}\t",
                    peak.start,
                    peak.end,
                    String::from_utf8_lossy(&peakprefix),
                    n_peak,
                    int10(Self::col(peak, score_column)),
                    peak.block_num,
                    String::from_utf8_lossy(&peak.block_sizes),
                    String::from_utf8_lossy(&peak.block_starts),
                )?;
                g6.clear();
                write_g6(&mut g6, peak.fc);
                fhd.write_all(g6.as_bytes())?;
                fhd.write_all(b"\t")?;
                g6.clear();
                write_g6(&mut g6, peak.pscore);
                fhd.write_all(g6.as_bytes())?;
                fhd.write_all(b"\t")?;
                g6.clear();
                write_g6(&mut g6, peak.qscore);
                fhd.write_all(g6.as_bytes())?;
                fhd.write_all(b"\n")?;
            }
        }
        Ok(())
    }

    /// Write ENCODE broadPeak (BED6+3) (`BroadPeakIO.write_to_broadPeak`).
    ///
    /// Groups by `end`; one line per group (first peak): chrom, start, end,
    /// name, `int(10*score_column)`, `.`, fc, pscore, qscore (`%.6g`).
    pub fn write_to_broad_peak<W: Write>(
        &self,
        fhd: &mut W,
        name_prefix: &[u8],
        name: &[u8],
        _description: &[u8],
        score_column: &str,
        trackline: bool,
    ) -> std::io::Result<()> {
        let peakprefix = interp_prefix(name_prefix, name);
        if trackline {
            let nm = String::from_utf8_lossy(name);
            write!(
                fhd,
                "track type=broadPeak name=\"{nm}\" description=\"{nm}\" nextItemButton=on\n"
            )?;
        }
        let mut n_peak: i32 = 0;
        let mut g6 = String::with_capacity(16);
        for (chrom, peaks) in self.iter() {
            for group in group_by_end_broad(peaks) {
                n_peak += 1;
                let peak = &group[0];
                fhd.write_all(chrom)?;
                write!(
                    fhd,
                    "\t{}\t{}\t{}{}\t{}\t.\t",
                    peak.start,
                    peak.end,
                    String::from_utf8_lossy(&peakprefix),
                    n_peak,
                    int10(Self::col(peak, score_column)),
                )?;
                g6.clear();
                write_g6(&mut g6, peak.fc);
                fhd.write_all(g6.as_bytes())?;
                fhd.write_all(b"\t")?;
                g6.clear();
                write_g6(&mut g6, peak.pscore);
                fhd.write_all(g6.as_bytes())?;
                fhd.write_all(b"\t")?;
                g6.clear();
                write_g6(&mut g6, peak.qscore);
                fhd.write_all(g6.as_bytes())?;
                fhd.write_all(b"\n")?;
            }
        }
        Ok(())
    }

    /// Write the MACS3 `.xls` report (`BroadPeakIO.write_to_xls`).
    pub fn write_to_xls<W: Write>(
        &self,
        ofhd: &mut W,
        name_prefix: &[u8],
        name: &[u8],
    ) -> std::io::Result<()> {
        writeln!(
            ofhd,
            "chr\tstart\tend\tlength\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"
        )?;
        let peakprefix = interp_prefix(name_prefix, name);
        let mut n_peak: i32 = 0;
        let mut g6 = String::with_capacity(16);
        for (chrom, peaks) in self.iter() {
            for group in group_by_end_broad(peaks) {
                n_peak += 1;
                let peak = &group[0];
                ofhd.write_all(chrom)?;
                write!(
                    ofhd,
                    "\t{}\t{}\t{}\t",
                    peak.start + 1,
                    peak.end,
                    peak.length,
                )?;
                g6.clear();
                write_g6(&mut g6, round2(peak.pileup));
                ofhd.write_all(g6.as_bytes())?;
                ofhd.write_all(b"\t")?;
                g6.clear();
                write_g6(&mut g6, peak.pscore);
                ofhd.write_all(g6.as_bytes())?;
                ofhd.write_all(b"\t")?;
                g6.clear();
                write_g6(&mut g6, peak.fc);
                ofhd.write_all(g6.as_bytes())?;
                ofhd.write_all(b"\t")?;
                g6.clear();
                write_g6(&mut g6, peak.qscore);
                ofhd.write_all(g6.as_bytes())?;
                write!(
                    ofhd,
                    "\t{}{}\n",
                    String::from_utf8_lossy(&peakprefix),
                    n_peak
                )?;
            }
        }
        Ok(())
    }
}

/// Group a chromosome's broad peaks into runs of consecutive entries sharing
/// the same `end` (`itertools.groupby(peaks, key=end)`).
fn group_by_end_broad(peaks: &[BroadPeak]) -> Vec<&[BroadPeak]> {
    let mut groups: Vec<&[BroadPeak]> = Vec::new();
    let mut i = 0;
    while i < peaks.len() {
        let end = peaks[i].end;
        let mut j = i + 1;
        while j < peaks.len() && peaks[j].end == end {
            j += 1;
        }
        groups.push(&peaks[i..j]);
        i = j;
    }
    groups
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_and_count_narrow() {
        let mut p = PeakIO::default();
        p.add(b"chr2", 10, 20, 15, 5.0, 3.0, 2.0, 1.5, 1.0, b"");
        p.add(b"chr1", 0, 10, 5, 4.0, 3.0, 2.0, 1.5, 1.0, b"");
        assert_eq!(p.total, 2);
        // bytewise-sorted chromosome order
        assert_eq!(p.get_chr_names(), vec![&b"chr1"[..], &b"chr2"[..]]);
        let np = &p.peaks_by_chr(b"chr1").unwrap()[0];
        assert_eq!(np.length, 10);
        assert_eq!(np.summit, 5);
    }

    #[test]
    fn sort_by_start() {
        let mut p = PeakIO::default();
        p.add(b"chr1", 100, 200, 150, 1.0, 0.0, 0.0, 0.0, 0.0, b"");
        p.add(b"chr1", 0, 50, 25, 1.0, 0.0, 0.0, 0.0, 0.0, b"");
        p.sort();
        let v = p.peaks_by_chr(b"chr1").unwrap();
        assert_eq!(v[0].start, 0);
        assert_eq!(v[1].start, 100);
    }

    #[test]
    fn broad_total_and_length() {
        let mut b = BroadPeakIO::new();
        b.add(
            b"chr1", 100, 500, 10.0, b".".to_vec(), b".".to_vec(), 1,
            b"400".to_vec(), b"0".to_vec(), 3.0, 2.0, 1.5, 1.0, b"NA".to_vec(),
        );
        assert_eq!(b.total(), 1);
        assert_eq!(b.peaks_by_chr(b"chr1").unwrap()[0].length, 400);
    }

    fn out(s: &[u8]) -> String {
        String::from_utf8(s.to_vec()).unwrap()
    }

    /// Ground-truthed against MACS3 `write_to_narrowPeak`.
    #[test]
    fn narrowpeak_matches_macs3() {
        let mut p = PeakIO::default();
        p.add(b"chr1", 100, 300, 200, 5.0, 12.3456, 4.5, 3.14159, 2.71828, b"");
        let mut buf: Vec<u8> = Vec::new();
        p.write_to_narrow_peak(&mut buf, b"test_narrowPeak", b"test", "score", true)
            .unwrap();
        let want = "track type=narrowPeak name=\"test\" description=\"test\" nextItemButton=on\n\
                    chr1\t100\t300\ttest_narrowPeak1\t50\t.\t3.14159\t4.5\t2.71828\t100\n";
        assert_eq!(out(&buf), want);
    }

    /// Ground-truthed against MACS3 `write_to_xls` (1-based, pileup round2).
    #[test]
    fn xls_matches_macs3() {
        let mut p = PeakIO::default();
        p.add(b"chr1", 100, 300, 200, 5.0, 12.3456, 4.5, 3.14159, 2.71828, b"");
        let mut buf: Vec<u8> = Vec::new();
        p.write_to_xls(&mut buf, b"%s_peak_", b"test").unwrap();
        let want = "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname\n\
                    chr1\t101\t300\t200\t201\t12.35\t4.5\t3.14159\t2.71828\ttest_peak_1\n";
        assert_eq!(out(&buf), want);
    }

    #[test]
    fn bed_and_summit_match_macs3() {
        let mut p = PeakIO::default();
        p.add(b"chr1", 100, 300, 200, 5.0, 12.3456, 4.5, 3.14159, 2.71828, b"");
        let mut bed: Vec<u8> = Vec::new();
        p.write_to_bed(&mut bed, b"%s_peak_", b"test", b"%s", "score", false)
            .unwrap();
        assert_eq!(out(&bed), "chr1\t100\t300\ttest_peak_1\t5\n");

        let mut sm: Vec<u8> = Vec::new();
        p.write_to_summit_bed(&mut sm, b"%s_peak_", b"test", b"%s", "score", false)
            .unwrap();
        assert_eq!(out(&sm), "chr1\t200\t201\ttest_peak_1\t5\n");
    }

    #[test]
    fn subpeak_letters_wraps() {
        // Faithful port of MACS3's recursion s(i//26)+chr(97+i%26): the
        // carry digit starts at 'b' (1), so index 26 -> "ba" (not "aa").
        assert_eq!(subpeak_letters(0), "a");
        assert_eq!(subpeak_letters(25), "z");
        assert_eq!(subpeak_letters(26), "ba");
        assert_eq!(subpeak_letters(27), "bb");
        assert_eq!(subpeak_letters(52), "ca");
    }
}
