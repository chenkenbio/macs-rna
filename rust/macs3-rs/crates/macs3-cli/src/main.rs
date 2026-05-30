//! macs3-rs CLI entry point.
//!
//! Mirrors MACS3's `bin/macs3` subcommand surface for the eight commands
//! macs-rna needs: `callpeak`, `predictd`, `pileup`, `bdgcmp`, `bdgopt`,
//! `bdgpeakcall`, `bdgbroadcall`, `bdgdiff`. Each subcommand wraps the per-
//! command `Args` struct from `cmd/<name>.rs` and dispatches to
//! `cmd::<name>::run(&args)`.
//!
//! This file and `cmd/mod.rs` are frozen (read-only) after Phase 0; the Phase 1
//! agents fill in the individual `cmd/<name>.rs::run` bodies.

use clap::{Parser, Subcommand};

pub mod cmd;

/// Top-level CLI (`bin/macs3` argparse root).
#[derive(Parser, Debug)]
#[command(
    name = "macs3-rs",
    version = macs3_core::constants::MACS_VERSION,
    about = "Fast Rust reimplementation of MACS3 hot paths (drop-in for macs3 3.0.4)"
)]
struct Cli {
    /// The MACS3 subcommand to run.
    #[command(subcommand)]
    command: Commands,
}

/// One variant per supported MACS3 subcommand.
#[derive(Subcommand, Debug)]
enum Commands {
    /// Call peaks from alignment results (`callpeak`).
    Callpeak(cmd::callpeak::Args),
    /// Predict fragment size `d` from alignment results (`predictd`).
    Predictd(cmd::predictd::Args),
    /// Pile up alignments into a bedGraph (`pileup`).
    Pileup(cmd::pileup::Args),
    /// Compare treatment vs control bedGraphs (`bdgcmp`).
    Bdgcmp(cmd::bdgcmp::Args),
    /// Modify the score column of a bedGraph (`bdgopt`).
    Bdgopt(cmd::bdgopt::Args),
    /// Call narrow peaks from a score bedGraph (`bdgpeakcall`).
    Bdgpeakcall(cmd::bdgpeakcall::Args),
    /// Call broad peaks from a score bedGraph (`bdgbroadcall`).
    Bdgbroadcall(cmd::bdgbroadcall::Args),
    /// Differential peak calling between two conditions (`bdgdiff`).
    Bdgdiff(cmd::bdgdiff::Args),
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Callpeak(a) => cmd::callpeak::run(a),
        Commands::Predictd(a) => cmd::predictd::run(a),
        Commands::Pileup(a) => cmd::pileup::run(a),
        Commands::Bdgcmp(a) => cmd::bdgcmp::run(a),
        Commands::Bdgopt(a) => cmd::bdgopt::run(a),
        Commands::Bdgpeakcall(a) => cmd::bdgpeakcall::run(a),
        Commands::Bdgbroadcall(a) => cmd::bdgbroadcall::run(a),
        Commands::Bdgdiff(a) => cmd::bdgdiff::run(a),
    }
}
