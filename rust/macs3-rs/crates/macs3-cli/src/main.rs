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
    /// `-T/--threads`: rayon worker threads for per-chromosome parallelism.
    /// `global = true` so it is accepted before OR after any subcommand. When
    /// unset, the default is resolved by [`resolve_threads`] (env-aware, and
    /// considerate on shared/cluster nodes — see that function).
    #[arg(short = 'T', long = "threads", global = true, value_name = "N")]
    threads: Option<usize>,

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

/// Resolve the number of rayon worker threads, with a server-friendly default.
///
/// Precedence (highest first):
/// 1. `--threads N` (explicit CLI override — wins over everything).
/// 2. `RAYON_NUM_THREADS` env (lets users/wrappers cap threads the rayon way).
/// 3. `SLURM_CPUS_ON_NODE` env (respect the Slurm allocation on cluster jobs).
/// 4. `min(8, available_parallelism)` — a polite default that never grabs a
///    whole 80-core shared node just because it can.
///
/// Any value `<= 0` from CLI/env is ignored in favour of the next source, and a
/// zero CLI value is bumped to 1 (rayon requires `>= 1`).
fn resolve_threads(cli_threads: Option<usize>) -> usize {
    if let Some(t) = cli_threads {
        return t.max(1);
    }
    for var in ["RAYON_NUM_THREADS", "SLURM_CPUS_ON_NODE"] {
        if let Ok(s) = std::env::var(var) {
            if let Ok(n) = s.trim().parse::<usize>() {
                if n >= 1 {
                    return n;
                }
            }
        }
    }
    let nproc = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    nproc.min(8)
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    // Initialise the rayon global pool ONCE, before any parallel work. Every
    // parallelized engine (callpeak + the bdg/pileup commands) uses the global
    // pool, so capping it here makes a single, predictable thread budget for the
    // whole process. `build_global` may only be called once; ignore the error if
    // a pool somehow already exists (it never does in this single entry point).
    let nthreads = resolve_threads(cli.threads);
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build_global();

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
