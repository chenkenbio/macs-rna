//! Paired-end fragment track skeleton, ported from
//! `MACS3/Signal/PairedEndTrack.py::PETrackI`.
//!
//! macs-rna does not use paired-end mode, so this is a minimal compile-only
//! skeleton: it fixes the public type and the `add_loc`/`finalize`/accessor
//! signatures so the rest of the workspace can name `PeTrack`, but the heavy
//! bodies (and pileup helpers) are stubbed. Phase 2 can flesh it out only if PE
//! support is ever required.
//!
//! Per chromosome MACS stores an `(start, end)` array of fragments; we mirror
//! that with parallel `Vec<i32>`s keyed by raw chromosome bytes.

use std::collections::BTreeMap;

/// One chromosome's `(start, end)` fragment intervals (`PETrackI.locations`).
#[derive(Debug, Clone, Default)]
pub struct ChromPe {
    /// Fragment start coordinates.
    pub starts: Vec<i32>,
    /// Fragment end coordinates.
    pub ends: Vec<i32>,
}

/// Paired-end fragment track (`PairedEndTrack.py::PETrackI`).
///
/// FROZEN API — Phase 1 depends on this; add methods, do not change signatures.
#[derive(Debug, Clone)]
pub struct PeTrack {
    /// Per-chromosome fragment intervals, bytewise-sorted by chromosome name.
    locations: BTreeMap<Vec<u8>, ChromPe>,
    /// Whether arrays are sorted (`is_sorted`).
    pub is_sorted: bool,
    /// Total number of fragments (`total`).
    pub total: u64,
    /// Sum of all fragment lengths (`length`).
    pub length: u64,
    /// Mean template length (`average_template_length`).
    pub average_template_length: f32,
    /// Growth chunk size hint (`buffer_size`).
    pub buffer_size: i64,
    /// Free-text annotation (`annotation`).
    pub annotation: String,
}

impl PeTrack {
    /// Create an empty PE track (`PETrackI.__init__`).
    pub fn new(anno: String, buffer_size: i64) -> Self {
        PeTrack {
            locations: BTreeMap::new(),
            is_sorted: false,
            total: 0,
            length: 0,
            average_template_length: 0.0,
            buffer_size,
            annotation: anno,
        }
    }

    /// Append a fragment `[start, end)` to `chromosome` (`add_loc`).
    pub fn add_loc(&mut self, chromosome: &[u8], start: i32, end: i32) {
        let c = self.locations.entry(chromosome.to_vec()).or_default();
        c.starts.push(start);
        c.ends.push(end);
        self.length += (end - start) as u64;
    }

    /// Fragment intervals for a chromosome (`get_locations_by_chr`).
    pub fn get_locations_by_chr(&self, chromosome: &[u8]) -> Option<&ChromPe> {
        self.locations.get(chromosome)
    }

    /// All chromosome names, bytewise-sorted (`get_chr_names`).
    pub fn get_chr_names(&self) -> Vec<&[u8]> {
        self.locations.keys().map(|k| k.as_slice()).collect()
    }

    /// Trim / sort / recompute totals (`finalize`). Stubbed — Phase 2 fills the
    /// PE sort-by-(start,end) and total/average-length recomputation.
    pub fn finalize(&mut self) {
        todo!("Phase 2 (PE): finalize")
    }

    /// Limit duplicate fragments to `maxnum` (`filter_dup`). Stubbed.
    pub fn filter_dup(&mut self, _maxnum: i32) -> u64 {
        todo!("Phase 2 (PE): filter_dup")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_loc_accumulates_length() {
        let mut t = PeTrack::new(String::new(), 100000);
        t.add_loc(b"chr1", 10, 30);
        t.add_loc(b"chr1", 0, 5);
        let c = t.get_locations_by_chr(b"chr1").unwrap();
        assert_eq!(c.starts, vec![10, 0]);
        assert_eq!(c.ends, vec![30, 5]);
        assert_eq!(t.length, 25);
    }
}
