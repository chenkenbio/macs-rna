//! Fixed-width fragment track, ported from
//! `MACS3/Signal/FixWidthTrack.py::FWTrack`.
//!
//! Per chromosome MACS stores two `int32` arrays of 5'-end cut positions, one
//! per strand (index `0` = plus, `1` = minus). We hold them in [`ChromLoc`],
//! keyed by raw chromosome bytes in a `BTreeMap` for bytewise-sorted iteration
//! (matching `set(sorted(self.locations.keys()))`).
//!
//! Phase 0 owns the struct, `add_loc`, `finalize`, `sort`, `total`, the
//! accessors and the trivial `filter_dup(maxnum<0)` short-circuit. The full
//! duplicate-filtering body and pileup helpers are filled by Phase 1A / Phase 2.

use std::collections::BTreeMap;

/// `INT32_MAX`, used as the default reference length for chromosomes whose true
/// length is unknown (`FixWidthTrack.py` imports `INT32_MAX as INT_MAX`).
pub const INT_MAX: i32 = i32::MAX;

/// Plus / minus strand 5'-end positions for one chromosome
/// (`FWTrack.locations[chrom]`).
#[derive(Debug, Clone, Default)]
pub struct ChromLoc {
    /// Plus-strand 5' cut positions (strand index `0`).
    pub plus: Vec<i32>,
    /// Minus-strand 5' cut positions (strand index `1`).
    pub minus: Vec<i32>,
}

/// Fixed-width fragment track (`FixWidthTrack.py::FWTrack`).
///
/// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
#[derive(Debug, Clone)]
pub struct FwTrack {
    /// Per-chromosome strand positions, bytewise-sorted by chromosome name.
    locations: BTreeMap<Vec<u8>, ChromLoc>,
    /// Chromosome names in first-seen (insertion) order. MACS's `FWTrack`
    /// stores `locations` in a plain Python `dict`, so iterating it yields
    /// chromosomes in the order they were first added. `pileup_and_write_se`
    /// emits chromosomes in `get_rlengths()` insertion order, which — when
    /// `set_rlengths` was never called (BED/SAM input) — equals this order.
    /// We track it separately because our `locations` map is a `BTreeMap`.
    insertion_order: Vec<Vec<u8>>,
    /// Reference chromosome lengths (`rlengths`); missing ones default to
    /// [`INT_MAX`].
    rlengths: BTreeMap<Vec<u8>, i32>,
    /// Whether `rlengths` was populated by an explicit `set_rlengths` call
    /// (BAM input carries header lengths). Mirrors MACS where BED/SAM never
    /// call `set_rlengths` while BAM does. Controls `rlengths_order`.
    rlengths_set: bool,
    /// Chromosome names in the order `get_rlengths()` would yield them once
    /// `set_rlengths` has run: sorted valid chromosomes first, then sorted
    /// missing ones — replicating MACS `set_rlengths`'s two `sorted(...)`
    /// loops over a Python dict.
    rlengths_order: Vec<Vec<u8>>,
    /// Whether per-strand arrays are sorted ascending (`is_sorted`).
    pub is_sorted: bool,
    /// Total number of retained 5' positions across both strands (`total`).
    pub total: u64,
    /// Fixed fragment width in bp (`fw`).
    pub fw: i32,
    /// `fw * total` after `finalize` (`length`).
    pub length: u64,
    /// Growth chunk size hint (`buffer_size`).
    pub buffer_size: i64,
    /// Free-text annotation (`annotation`).
    pub annotation: String,
}

impl FwTrack {
    /// Create an empty track (`FWTrack.__init__`).
    pub fn new(fw: i32, anno: String, buffer_size: i64) -> Self {
        FwTrack {
            locations: BTreeMap::new(),
            insertion_order: Vec::new(),
            rlengths: BTreeMap::new(),
            rlengths_set: false,
            rlengths_order: Vec::new(),
            is_sorted: false,
            total: 0,
            fw,
            length: 0,
            buffer_size,
            annotation: anno,
        }
    }

    /// Append a 5' cut position (`add_loc`). `strand` is `0` (plus) or `1`
    /// (minus).
    pub fn add_loc(&mut self, chromosome: &[u8], fiveendpos: i32, strand: i32) {
        if !self.locations.contains_key(chromosome) {
            self.insertion_order.push(chromosome.to_vec());
        }
        let c = self.locations.entry(chromosome.to_vec()).or_default();
        if strand == 0 {
            c.plus.push(fiveendpos);
        } else {
            c.minus.push(fiveendpos);
        }
    }

    /// Trim, sort per strand, and refresh `total`/`length` (`finalize`).
    pub fn finalize(&mut self) {
        self.total = 0;
        for c in self.locations.values_mut() {
            c.plus.sort_unstable();
            c.minus.sort_unstable();
            self.total += (c.plus.len() + c.minus.len()) as u64;
        }
        self.is_sorted = true;
        self.length = self.fw as u64 * self.total;
    }

    /// Sort per-strand coordinate arrays for every chromosome (`sort`).
    pub fn sort(&mut self) {
        for c in self.locations.values_mut() {
            c.plus.sort_unstable();
            c.minus.sort_unstable();
        }
        self.is_sorted = true;
    }

    /// Attach reference chromosome lengths (`set_rlengths`); chromosomes present
    /// in the track but absent from `rlengths` get [`INT_MAX`].
    pub fn set_rlengths(&mut self, rlengths: &BTreeMap<Vec<u8>, i32>) {
        // Replicate `FixWidthTrack.py::set_rlengths`: split the track's
        // chromosomes into those present in `rlengths` (valid) and those
        // absent (missed), filling each in sorted order. The resulting dict
        // insertion order is `sorted(valid) ++ sorted(missed)` — which we
        // record in `rlengths_order` so `pileup` can iterate it.
        self.rlengths_order.clear();
        let mut valid: Vec<&Vec<u8>> = Vec::new();
        let mut missed: Vec<&Vec<u8>> = Vec::new();
        for chrom in self.locations.keys() {
            if rlengths.contains_key(chrom) {
                valid.push(chrom);
            } else {
                missed.push(chrom);
            }
        }
        // `self.locations` is a BTreeMap, so `keys()` is already bytewise
        // sorted; `valid`/`missed` therefore retain sorted order.
        for chrom in valid {
            let v = rlengths.get(chrom).copied().unwrap_or(INT_MAX);
            self.rlengths.insert(chrom.clone(), v);
            self.rlengths_order.push(chrom.clone());
        }
        for chrom in missed {
            self.rlengths.insert(chrom.clone(), INT_MAX);
            self.rlengths_order.push(chrom.clone());
        }
        self.rlengths_set = true;
    }

    /// Return reference chromosome lengths, defaulting all to [`INT_MAX`] when
    /// none were set (`get_rlengths`).
    pub fn get_rlengths(&mut self) -> &BTreeMap<Vec<u8>, i32> {
        if self.rlengths.is_empty() {
            for k in self.locations.keys() {
                self.rlengths.insert(k.clone(), INT_MAX);
            }
        }
        &self.rlengths
    }

    /// Chromosome names in the order `get_rlengths()` / `pileup_and_write_se`
    /// would iterate them (`list(chrlengths.keys())`).
    ///
    /// * If `set_rlengths` ran (BAM input): `sorted(valid) ++ sorted(missed)`.
    /// * Otherwise (BED/SAM input): first-seen insertion order, matching the
    ///   plain Python `dict` built lazily inside `get_rlengths`.
    pub fn rlengths_iter_order(&self) -> Vec<&[u8]> {
        if self.rlengths_set {
            self.rlengths_order.iter().map(|k| k.as_slice()).collect()
        } else {
            self.insertion_order.iter().map(|k| k.as_slice()).collect()
        }
    }

    /// Strand arrays for a chromosome (`get_locations_by_chr`).
    pub fn get_locations_by_chr(&self, chromosome: &[u8]) -> Option<&ChromLoc> {
        self.locations.get(chromosome)
    }

    /// All chromosome names, bytewise-sorted (`get_chr_names`).
    pub fn get_chr_names(&self) -> Vec<&[u8]> {
        self.locations.keys().map(|k| k.as_slice()).collect()
    }

    // ------------------------------------------------------------------
    // Stubbed heavy paths — Phase 1A / Phase 2 fill these in this file.
    // ------------------------------------------------------------------

    /// Limit duplicate 5' positions to `maxnum` per coordinate per strand
    /// (`filter_dup`). The `maxnum < 0` short-circuit ("do nothing") is
    /// implemented here; the actual de-duplication scan is Phase 1A.
    pub fn filter_dup(&mut self, maxnum: i32) -> u64 {
        if maxnum < 0 {
            return self.total; // do nothing
        }
        if !self.is_sorted {
            self.sort();
        }

        self.total = 0;
        self.length = 0;

        for c in self.locations.values_mut() {
            filter_strand(&mut c.plus, maxnum, &mut self.total);
            filter_strand(&mut c.minus, maxnum, &mut self.total);
        }

        self.length = self.fw as u64 * self.total;
        self.total
    }
}

/// In-place duplicate filter for one strand's sorted position array, mirroring
/// the per-strand block of `FixWidthTrack.py::filter_dup`.
///
/// Replicates a quirk of the reference: when a strand has `<= 1` element it is
/// left untouched **and its count is not added to `total`** (MACS only does
/// `self.total += i_new` inside the `else` branch). Filtering keeps the first
/// position unconditionally, then keeps each subsequent position only while the
/// running run-length `n` at the current coordinate is `<= maxnum`.
fn filter_strand(arr: &mut Vec<i32>, maxnum: i32, total: &mut u64) {
    let size = arr.len();
    if size <= 1 {
        return; // do nothing; not counted into total (matches MACS)
    }
    let mut new_arr: Vec<i32> = Vec::with_capacity(size + 1);
    new_arr.push(arr[0]); // first item
    let mut n: i32 = 1;
    let mut current_loc = arr[0];
    for i_old in 1..size {
        let p = arr[i_old];
        if p == current_loc {
            n += 1;
        } else {
            current_loc = p;
            n = 1;
        }
        if n <= maxnum {
            new_arr.push(p);
        }
    }
    let i_new = new_arr.len() as u64;
    *total += i_new;
    *arr = new_arr;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_finalize_total() {
        let mut t = FwTrack::new(0, String::new(), 100000);
        t.add_loc(b"chr1", 30, 0);
        t.add_loc(b"chr1", 10, 0);
        t.add_loc(b"chr1", 20, 1);
        t.finalize();
        let c = t.get_locations_by_chr(b"chr1").unwrap();
        assert_eq!(c.plus, vec![10, 30]); // sorted
        assert_eq!(c.minus, vec![20]);
        assert_eq!(t.total, 3);
        assert!(t.is_sorted);
    }

    #[test]
    fn filter_dup_negative_is_noop() {
        let mut t = FwTrack::new(0, String::new(), 100000);
        t.add_loc(b"chr1", 10, 0);
        t.finalize();
        assert_eq!(t.filter_dup(-1), t.total);
    }

    #[test]
    fn chr_names_sorted() {
        let mut t = FwTrack::new(0, String::new(), 100000);
        t.add_loc(b"chr2", 1, 0);
        t.add_loc(b"chr1", 1, 0);
        assert_eq!(t.get_chr_names(), vec![&b"chr1"[..], &b"chr2"[..]]);
    }

    /// Golden values ported from `references/MACS/test/test_FixWidthTrack.py`
    /// (`Test_FWTrack::test_filter_dup`). Filtering is destructive and applied
    /// successively, so the same track is filtered with shrinking thresholds.
    #[test]
    fn filter_dup_golden() {
        let fw = 50;
        let input_regions: &[(&[u8], i32, i32)] = &[
            (b"chrY", 0, 0),
            (b"chrY", 90, 0),
            (b"chrY", 150, 0),
            (b"chrY", 70, 0),
            (b"chrY", 80, 0),
            (b"chrY", 85, 0),
            (b"chrY", 85, 0),
            (b"chrY", 85, 0),
            (b"chrY", 85, 0),
            (b"chrY", 90, 1),
            (b"chrY", 150, 1),
            (b"chrY", 70, 1),
            (b"chrY", 80, 1),
            (b"chrY", 80, 1),
            (b"chrY", 80, 1),
            (b"chrY", 85, 1),
            (b"chrY", 90, 1),
        ];
        let mut t = FwTrack::new(fw, String::new(), 100000);
        for &(c, p, s) in input_regions {
            t.add_loc(c, p, s);
        }
        t.finalize();
        assert_eq!(t.total, 17);
        assert_eq!(t.length, 17 * fw as u64);

        // filter out more than 3 tags -> one chrY:85:+ removed
        assert_eq!(t.filter_dup(3), 16);
        assert_eq!(t.total, 16);
        // filter out more than 2 tags -> one chrY:85:+ and one chrY:80:- removed
        assert_eq!(t.filter_dup(2), 14);
        // filter out more than 1 tag -> chrY:85:+, chrY:80:-, chrY:90:- removed
        assert_eq!(t.filter_dup(1), 11);
    }

    #[test]
    fn filter_dup_negative_returns_total() {
        let mut t = FwTrack::new(0, String::new(), 100000);
        t.add_loc(b"chr1", 10, 0);
        t.add_loc(b"chr1", 10, 0);
        t.finalize();
        assert_eq!(t.filter_dup(-1), 2);
    }

    /// `rlengths_iter_order` falls back to first-seen insertion order when
    /// `set_rlengths` was never called (BED/SAM path).
    #[test]
    fn insertion_order_without_set_rlengths() {
        let mut t = FwTrack::new(0, String::new(), 100000);
        t.add_loc(b"chr2", 1, 0);
        t.add_loc(b"chr1", 1, 0);
        t.add_loc(b"chr2", 2, 0);
        assert_eq!(t.rlengths_iter_order(), vec![&b"chr2"[..], &b"chr1"[..]]);
    }

    /// After `set_rlengths`, ordering is `sorted(valid) ++ sorted(missed)`.
    #[test]
    fn rlengths_order_after_set() {
        let mut t = FwTrack::new(0, String::new(), 100000);
        t.add_loc(b"chr2", 1, 0);
        t.add_loc(b"chr1", 1, 0);
        t.add_loc(b"chr3", 1, 0);
        let mut rl = BTreeMap::new();
        rl.insert(b"chr3".to_vec(), 1000);
        rl.insert(b"chr1".to_vec(), 2000);
        // chr2 is missing -> goes last
        t.set_rlengths(&rl);
        assert_eq!(
            t.rlengths_iter_order(),
            vec![&b"chr1"[..], &b"chr3"[..], &b"chr2"[..]]
        );
    }
}
