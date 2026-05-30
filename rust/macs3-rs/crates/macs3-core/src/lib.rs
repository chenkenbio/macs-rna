//! macs3-core: pure-Rust reimplementation of MACS3 (3.0.4) signal processing.
//!
//! Module layout mirrors the MACS3 source it ports from (see the plan's
//! module-mapping table). The foundation phase fills in `fmt`, `prob`,
//! `bedgraph_io`, `constants`, and the core data structures; the fan-out
//! phases fill in the per-subcommand engines.
//!
//! Bit-exactness contract: all internal pileup/score values are `f32`; output
//! formatting replicates C `printf("%.5f")` / `%.6g`. Never use fast-math.

// --- Foundation (Phase 0): shared contracts, frozen after this phase ---
pub mod constants;
pub mod fmt;
pub mod prob;

// --- Core data structures + per-subcommand engines ---
// Ownership for Phase 1 (see ARCHITECTURE.md):
//   A: parser, pileup, track_fw, track_pe   B: bedgraph, peak_io
//   C: scoretrack, twocond                  E: peakmodel
//   Phase 2: caller, peakdetect
// bedgraph_io is frozen (read-only) after Phase 0.
pub mod bedgraph;
pub mod bedgraph_io;
pub mod caller;
pub mod parser;
pub mod peak_io;
pub mod peakdetect;
pub mod peakmodel;
pub mod pileup;
pub mod scoretrack;
pub mod track_fw;
pub mod track_pe;
pub mod twocond;
