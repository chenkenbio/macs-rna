//! Project-wide constants and defaults, ported from
//! `MACS3/Utilities/Constants.py` and `MACS3/Utilities/OptValidator.py`.
//!
//! Values are reproduced verbatim from MACS3 3.0.4 so downstream phases share
//! one authoritative source. Where MACS reads a value from argparse defaults,
//! the constant name mirrors the argparse `dest`.

/// MACS3 version string (`Constants.py::MACS_VERSION`).
pub const MACS_VERSION: &str = "3.0.4";

/// Maximum number of pairs used while building the shifting model
/// (`Constants.py::MAX_PAIRNUM`).
pub const MAX_PAIRNUM: i32 = 1000;

/// Upper clamp on local lambda (`Constants.py::MAX_LAMBDA`).
pub const MAX_LAMBDA: i32 = 100000;

/// Fold-enrichment step for the FE cutoff sweep (`Constants.py::FESTEP`).
pub const FESTEP: i32 = 20;

/// Growth chunk size for internal arrays (`Constants.py::BUFFER_SIZE`).
pub const BUFFER_SIZE: usize = 100000;

/// Read buffer size in bytes (`Constants.py::READ_BUFFER_SIZE`).
pub const READ_BUFFER_SIZE: usize = 10000000;

/// Default number of processes (`Constants.py::N_MP`).
pub const N_MP: i32 = 2;

/// Effective genome sizes, keyed by the shortcut accepted on the command line
/// (`Constants.py::EFFECTIVEGS`). Collected from deepTools.
///
/// * `hs` — GRCh38 human
/// * `mm` — GRCm38 mouse
/// * `ce` — WBcel235 C. elegans
/// * `dm` — dm6 fruitfly
pub const EFFECTIVE_GS: [(&str, u64); 4] = [
    ("hs", 2_913_022_398),
    ("mm", 2_652_783_500),
    ("ce", 100_286_401),
    ("dm", 142_573_017),
];

/// Resolve an effective genome-size shortcut (`hs`/`mm`/`ce`/`dm`) to its
/// numeric value as an `f64` (MACS stores `gsize` as a float). Returns `None`
/// for unknown shortcuts so the caller can fall back to `parse::<f64>()`,
/// mirroring `OptValidator.py`'s `efgsize[...]` / `float(...)` logic.
pub fn effective_gsize(shortcut: &str) -> Option<f64> {
    EFFECTIVE_GS
        .iter()
        .find(|(k, _)| *k == shortcut)
        .map(|(_, v)| *v as f64)
}

// ------------------------------------------------------------------
// Default option values (from bin/macs3 argparse + OptValidator.py).
// ------------------------------------------------------------------

/// Default small local region for dynamic lambda (`--slocal`).
pub const DEFAULT_SLOCAL: i32 = 1000;

/// Default large local region for dynamic lambda (`--llocal`).
pub const DEFAULT_LLOCAL: i32 = 10000;

/// Default extension size (`--extsize`).
pub const DEFAULT_EXTSIZE: i32 = 200;

/// Default band width for the shifting model (`--bw`).
pub const DEFAULT_BW: i32 = 300;

/// Default minimum fragment size (`--d-min`).
pub const DEFAULT_D_MIN: i32 = 20;

/// Default q-value cutoff (`-q`).
pub const DEFAULT_QVALUE: f64 = 0.05;

/// Default broad-region cutoff (`--broad-cutoff`).
pub const DEFAULT_BROAD_CUTOFF: f64 = 0.1;

/// Default fold-enrichment cutoff for post-filtering (`--fe-cutoff`).
pub const DEFAULT_FECUTOFF: f64 = 1.0;

/// Default mfold lower / upper bounds (`-m`).
pub const DEFAULT_MFOLD: (i32, i32) = (5, 50);

/// Default effective genome-size shortcut (`-g`).
pub const DEFAULT_GSIZE: &str = "hs";

/// Default `bdgpeakcall` cutoff (`-c`).
pub const DEFAULT_BDGPEAKCALL_CUTOFF: f64 = 5.0;

/// Default `bdgpeakcall` minimum length (`-l`).
pub const DEFAULT_BDGPEAKCALL_MINLEN: i32 = 200;

/// Default `bdgpeakcall` maximum gap (`-g`).
pub const DEFAULT_BDGPEAKCALL_MAXGAP: i32 = 30;

/// Default `bdgbroadcall` peak cutoff (`-c`).
pub const DEFAULT_BDGBROADCALL_CUTOFF_PEAK: f64 = 2.0;

/// Default `bdgbroadcall` link cutoff (`-C`).
pub const DEFAULT_BDGBROADCALL_CUTOFF_LINK: f64 = 1.0;

/// Default `bdgbroadcall` level-1 maximum gap (`-g`).
pub const DEFAULT_BDGBROADCALL_LVL1_MAXGAP: i32 = 30;

/// Default `bdgbroadcall` level-2 maximum gap (`-G`).
pub const DEFAULT_BDGBROADCALL_LVL2_MAXGAP: i32 = 800;

/// Default `bdgdiff` log10 likelihood-ratio cutoff (`-C`).
pub const DEFAULT_BDGDIFF_CUTOFF: f64 = 3.0;

/// Default `bdgdiff` minimum length (`-l`).
pub const DEFAULT_BDGDIFF_MINLEN: i32 = 200;

/// Default `bdgdiff` maximum gap (`-g`).
pub const DEFAULT_BDGDIFF_MAXGAP: i32 = 100;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gsize_lookup() {
        assert_eq!(effective_gsize("hs"), Some(2_913_022_398.0));
        assert_eq!(effective_gsize("dm"), Some(142_573_017.0));
        assert_eq!(effective_gsize("nope"), None);
    }
}
